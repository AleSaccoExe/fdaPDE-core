#ifndef __SIMPLIFICATION_H__
#define __SIMPLIFICATION_H__

#include "mesh.h"
#include "StructuredGridSearch.h"
#include "connections.h"

namespace fdapde{
namespace core{
template<int M, int N>
class Simplification{
	using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
	using FacetType = std::array<int, Mesh<M, N>::n_vertices_per_facet>;
private:
	StructuredGridSearch<M, N> sgs_;
	Connections connections_;
	std::vector<Element<M, N>> elems_vec_;
	DMatrix<double> nodes_;
	DMatrix<int> elems_mat_;
	std::vector<FacetType> facets_;
	DMatrix<int> boundary_;

	DMatrix<double> data_;
	ConnectionsType data_to_elems_;
    ConnectionsType elem_to_data_;

    std::multimap<double, std::pair<unsigned, SVector<N> >> costs_map_;
    // mappa utile per l'update di collapse_costs dopo un collapse
    std::unordered_map<unsigned, double> facets_cost_;

    unsigned n_nodes_;

    // metodi privati

    // calcola il costo per ogni facet
    template<typename CostType_>
    void compute_costs(CostType_ cost_obj);
    // trova i possibili punti per il collapse di facet
    std::vector<SVector<N>> get_collapse_points(const FacetType & facet) const;
public:
	Simplification(const Mesh<M, N> & mesh);
	template<typename CostType_>
	void simplify(CostType_ cost_obj);
};

// ===============
// implementazione
// ===============

// costruttore 
template<int M, int N>
Simplification<M, N>::Simplification(const Mesh<M, N> & mesh):
	sgs_(mesh), connections_(mesh), elems_mat_(mesh.elements()), nodes_(mesh.nodes()),
	elem_to_data_(elems_mat_.rows()), n_nodes_(mesh.n_nodes()), boundary_(mesh.boundary())
{
	// viene costruito il vettore di elementi
	elems_vec_.reserve(mesh.n_elements());
	for(unsigned elem_id = 0; elem_id < mesh.n_elements(); ++elem_id)
		elems_vec_.push_back(mesh.element(elem_id));
	// vengono costruite le facce
	facets_.reserve(mesh.n_facets());
    for(unsigned facet_id = 0; facet_id < mesh.n_facets(); ++facet_id)
        facets_.push_back(mesh.facet(facet_id).node_ids());
    // vengono costruite le connessioni elementi - dati
    data_ = mesh.nodes();
    data_to_elems_.resize(data_.rows());
    for(unsigned i = 0; i < data_.rows(); ++i) // loop su tutti i dati
    {
        auto tmp_set = connections_.get_node_to_elems(i);
        data_to_elems_[i] = tmp_set;
        for(unsigned elem_id : tmp_set) // loop sugli elementi connessi al dato i
            elem_to_data_[elem_id].insert(i);
    }

}

template<int M, int N>
template<typename CostType_>
void Simplification<M, N>::compute_costs(CostType_ cost_obj)
{
	unsigned facet_id = 0;
    for(const auto & facet : facets_){
	    if(connections_.nodes_on_facet(facet).size()==2)
	    {
	    std::vector<SVector<N>> collapse_points = get_collapse_points(facet);
	    // per ogni punto in collapse_points si calcola il costo della contrazione
	    std::map<double, SVector<N>> tmp_costs_map;
	    auto elems_to_modify = connections_.elems_modified_in_collapse(facet);
	    auto elems_to_erase = connections_.elems_erased_in_collapse(facet);
	    for(auto collapse_point : collapse_points)
	    {
	        // viene creato un vettore degli elementi involved nel collapse di facet
	        std::vector<Element<M, N>> elems_tmp;
	        for(unsigned elem_id : elems_to_modify)
	            elems_tmp.push_back(elems_vec_[elem_id]);
	        for(unsigned elem_id : elems_to_erase)
	            elems_tmp.push_back(elems_vec_[elem_id]);
	        // calcolato il costo del collapse
	        tmp_costs_map[cost_obj(elems_tmp, collapse_point)] = collapse_point;
	    }
	    // ora si prende il costo minore e si controllano le intersezioni
	    for(auto it = tmp_costs_map.begin(); it != tmp_costs_map.end(); ++it)
	    {
	        // viene creato un vettore degli elementi modificati nel collapse di facet
	        std::vector<Element<M, N>> elems_tmp;
	        for(unsigned elem_to_modify : elems_to_modify)
	        {
	            auto elem_tmp = elems_vec_[elem_to_modify];
	            std::array<int, ct_nvertices(M)> node_ids;
	            std::array<SVector<N>, ct_nvertices(M)> coords; 
	            // il vettore degli id viene riempito
	            for(unsigned i = 0; i < ct_nvertices(M); ++i){
	                node_ids[i] = elems_mat_(elem_to_modify, i);
	                coords[i] = nodes_.row(node_ids[i]);
	            }
	            for(int i = 0; i < Mesh<M, N>::n_vertices; ++i)
	                for(int j = 1; j < Mesh<M, N>::local_dimension; ++j)
	                    if(facet[j] == node_ids[i]){
	                        node_ids[i] = facet[0];
	                        coords[i] = it->second;
	            }
	            // si costruisce ora un nuovo elemento con l'id del nodo modificato
	            elems_tmp.emplace_back(elem_tmp.ID(), node_ids, coords, elem_tmp.neighbors(), elem_tmp.is_on_boundary());
	        }
	        // ora passo a sgs il vettore di elementi elems_tmp, che è formato di elementi da modificare
	        // e il set elems_to_erase. La classe provvede quindi a capire le possibili intersezioni date dal collapse di facet
	        sgs_.update(elems_tmp, elems_to_erase);
	        bool valid_collapse = true;
	        for(auto & elem : elems_tmp){
	            auto elems_to_check = sgs_.get_neighbouring_elements(elem);
	            for(unsigned elem_id : elems_to_check)
	                if(elems_vec_[elem_id].intersection(elem))
	                {
	                    // std::cout<<"trovata intersezione\n";
	                    valid_collapse = false;
	                    break;
	                }  
	            if(!valid_collapse)
	                break;
	        }
	        // se non ci sono intersezioni si procede ad aggiungere le informazioni nelle strutture adeguate
	        if(valid_collapse)
	        {
	            costs_map_.insert({it->first, {facet_id, it->second} });
	            facets_cost_[facet_id] = it->first;
	            break;
	        }

	    }
	    std::vector<Element<M, N>> elems_tmp; // questo vettore viene ora riempito di elementi da aggiungere
	    // a questo punto si ripristina la classe sgs al suo stato iniziale
	    for(unsigned elem_id : elems_to_erase) // loop sugli elementi che erano stati eliminati
	        elems_tmp.push_back(elems_vec_[elem_id]);
	    sgs_.add_elements(elems_tmp);
	    elems_tmp.clear(); // ora viene riempito di elementi da modificare
	    for(unsigned elem_id : elems_to_modify)
	        elems_tmp.push_back(elems_vec_[elem_id]);
	    sgs_.update_f(elems_tmp);

	    } // if(connections_.nodes_on_facet(facet).size()==2)
	    facet_id++;
    }
}

template<int M, int N>
std::vector<SVector<N>> Simplification<M, N>::get_collapse_points(const FacetType & facet) const
{
	int boundary_nodes = 0;
	// viene contato il numero di nodi sul bordo
	for(unsigned node_id : facet) 
		boundary_nodes += boundary_(node_id);
	// se ci sono almeno due nodi sul bordo facet non può essere contratta
	if(boundary_nodes > 1)
		return {};
	std::vector<SVector<N>> collapse_points;
	// se non ci sono nodi sul bordo i possibili punti di contrazione sono:
	// i vertici di facet i il suo punto medio
	if(boundary_nodes == 0){
		for(unsigned id : facet)
			collapse_points.push_back(nodes_.row(id));
		SMatrix<N, M> J;
		for (std::size_t j = 0; j < M; ++j) { 
			SVector<N> col = nodes_.row(facet[j+1]) - nodes_.row(facet[0]);
			J.col(j) = col; 
		}
		SVector<M> barycentric_mid_point;
	    barycentric_mid_point.fill(1.0 / (M + 1));
		collapse_points.push_back(J * barycentric_mid_point + ref );
	}
	// c'è un solo vertice sul bordo: è l'unico punto di contrazione
	if(boundary_nodes == 1)
	{
		for(unsigned id : facet)
			if(boundary_(id) == 1)
				collapse_points.push_back(nodes_.row(id));
	}
	return collapse_points;

}


// ========================
// implementazione simplify
// ========================
template<int M, int N>
template<typename CostType_>
void Simplification<M, N>::simplify(CostType_ cost_obj)
{
	// vengono calcolati i costi per ogni facet valida
	compute_costs(cost_obj);

}

} // fdapde
} // core

#endif // __SIMPLIFICATION_H__