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
    void compute_costs(const std::set<unsigned> & facet_ids,CostType_ cost_obj);
    // trova i possibili punti per il collapse di facet
    std::vector<SVector<N>> get_collapse_points(const FacetType & facet) const;
    // elimina le facce in input
    void erase_facets(const std::unordered_set<unsigned> & facet_ids);
    // modifica le facce in input usando le informazioni in facet
    void modify_facets(const std::unordered_set<unsigned> & facet_ids, const FacetType & facet);
    std::vector<Element<M, N>> modify_elements(std::unordered_set<unsigned> & elem_ids, const FacetType & facet, SVector<N> new_coord) const;
public:
	Simplification(const Mesh<M, N> & mesh);
	template<typename CostType_>
	void simplify(unsigned n_nodes, CostType_ cost_obj);
	// costruisce la mesh semplificata
	Mesh<M, N> build_mesh() const;
	// temporaneo
	const DMatrix<double> & get_data() const {return data_;}
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
std::vector<Element<M, N>> Simplification<M, N>::modify_elements(std::unordered_set<unsigned> & elem_ids, const FacetType & facet, SVector<N> new_coord) const
{
	std::vector<Element<M, N>> elems;
	for(unsigned elem_to_modify : elem_ids)
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
            for(int j = 0; j < Mesh<M, N>::local_dimension; ++j)
                if(facet[j] == node_ids[i]){
                    node_ids[i] = facet[0];
                    coords[i] = new_coord;
        }
        // si costruisce ora un nuovo elemento con l'id del nodo modificato
        elems.emplace_back(elem_tmp.ID(), node_ids, coords, elem_tmp.neighbors(), elem_tmp.is_on_boundary());
    }
    return elems;

}

template<int M, int N>
template<typename CostType_>
void Simplification<M, N>::compute_costs(const std::set<unsigned> & facet_ids, CostType_ cost_obj)
{
	for(unsigned facet_id : facet_ids){
		const auto & facet = facets_[facet_id];
    	std::vector<SVector<N>> collapse_points = get_collapse_points(facet);
    	// viene aggiunto un costo se: sulla faccia ci sono 2 nodi
	    if(connections_.nodes_on_facet(facet).size()==2 && collapse_points.size() != 0)
	    {
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
	        std::vector<Element<M, N>> elems_tmp = modify_elements(elems_to_modify, facet, it->second);
	        // ora passo a sgs il vettore di elementi elems_tmp, che è formato di elementi da modificare
	        // e il set elems_to_erase. La classe provvede quindi a capire le possibili intersezioni date dal collapse di facet
	        sgs_.update(elems_tmp, elems_to_erase, false);
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
} // compute_costs

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
		// collapse_points.push_back(J * barycentric_mid_point + nodes_.row(facet[0]) );
	}
	// c'è un solo vertice sul bordo: è l'unico punto di contrazione
	else if(boundary_nodes == 1)
	{
		for(unsigned id : facet)
			if(boundary_(id) == 1)
				collapse_points.push_back(nodes_.row(id));
	}
	return collapse_points;
} // get_collapse_points

template<int M, int N>
void Simplification<M, N>::erase_facets(const std::unordered_set<unsigned> & facet_ids)
{
	// assert(facet_ids.size() == 2);
	for(unsigned facet_id : facet_ids){
		if(facets_cost_.find(facet_id) != facets_cost_.end()){
			double cost = facets_cost_.at(facet_id);
			auto it = costs_map_.find(cost);
			assert(it != costs_map_.end());
			while(it->second.first != facet_id)
				++it;
			costs_map_.erase(it);
			facets_cost_.erase(facet_id);
		}
	}
}

template<int M, int N>
void Simplification<M, N>::modify_facets(const std::unordered_set<unsigned> & facet_ids, const FacetType & facet)
{
	for(unsigned facet_id : facet_ids) // loop sulle facce da modificare
    {
        auto & facet_to_modify = facets_[facet_id];
        for(int i = 0; i < Mesh<M, N>::n_vertices_per_facet; ++i)
            for(int j = 0; j < Mesh<M, N>::n_vertices_per_facet; ++j)
                if(facet_to_modify[i] == facet[j]){  facet_to_modify[i] = facet[0]; }
    }

}


// ========================
// implementazione simplify
// ========================
template<int M, int N>
template<typename CostType_>
void Simplification<M, N>::simplify(unsigned n_nodes, CostType_ cost_obj)
{
	if(n_nodes >= n_nodes_)
	{
		std::cout<<"required nodes: "<<n_nodes<<", nodes in the mesh: "<<n_nodes_<<std::endl;
		return;
	}

	// vengono calcolati i costi per ogni facet valida
	std::set<unsigned> all_facets;
	for(unsigned i = 0; i < facets_.size(); ++i) {all_facets.insert(i);}
	compute_costs(all_facets, cost_obj);
	// viene semplificata la faccia con costo minore
	for(auto collapse_info : costs_map_)
	{
		auto facet = facets_[collapse_info.second.first]; // viene presa la faccia da contrarre
		auto elems_to_modify = connections_.elems_modified_in_collapse(facet);
        auto elems_to_erase = connections_.elems_erased_in_collapse(facet);
        // si modificano gli elementi
        // il vettore viene riempito con gli elementi modificati
		std::vector<Element<M, N>> elems_tmp = modify_elements(elems_to_modify, facet, collapse_info.second.second);
		// vengono aggiornate le informazioni della classe sgs_
		sgs_.update(elems_tmp, elems_to_erase, true);
		// vengono presi i dati da proiettare
        std::unordered_set<unsigned> data_ids;
        for(unsigned elem_id : elems_to_erase)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        for(unsigned elem_id : elems_to_modify)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        // vengono proiettati i dati sugli elementi già modificati
        auto new_data_positions = project(elems_tmp, data_, data_ids);
        // ora si può eseguire il collapse e dopo si cambiano le informazioni di data_to_elems e elem_to_data
        for(auto & elem : elems_tmp){ // loop sugli elementi modificati dalla contrazione
            // viene sostituito nel vettore l'elemento modificato
            elems_vec_[elem.ID()] = elem;
            // nella matrice degli elementi vengono modificati gli id dei nodi
            for(unsigned i = 0; i < ct_nvertices(M); ++i)
                elems_mat_(elem.ID(), i) = elem.node_ids()[i];
        }
       	// vengono cambiate le coordinate del nodo del collapse
       	nodes_.row(facet[0]) = collapse_info.second.second;
       	elems_tmp.clear(); // questo vettore viene ora riempito di elementi da eliminare
        for(unsigned elem_id : elems_to_erase) // loop sugli elementi da eliminare
            elems_tmp.push_back(elems_vec_[elem_id]);
        // vengono modificate le connessioni e prese le informazioni necessarie a modificare
        // le facce
        auto facets_pair = connections_.collapse_facet(facet, elems_tmp);
        assert(facets_pair.first.size()==2);
        // vengono eliminate le informazioni sulle facce che non esistono più
        erase_facets(facets_pair.first);
        // e anche tutte le altre facce con costi da ricalcolare
        auto facets_to_update = connections_.facets_to_update(facet[0]);
        erase_facets({facets_to_update.begin(), facets_to_update.end()});
        // vengono modificati gli id dei vertici delle facce
        modify_facets(facets_pair.second, facet);
        // vengono ora modificate le connessioni tra dati e elementi
        unsigned i = 0; // id del dato
        for(unsigned data_id : data_ids)
        {
            // prendo i vecchi id degli elementi connessi al dato data_id
            auto old_data_to_elems = data_to_elems_[data_id];
            // rimuovo le vecchie connessioni
            for(unsigned elem_id : old_data_to_elems) // loop sugli elementi
                elem_to_data_[elem_id].erase(data_id);
            auto data_info = new_data_positions[i];
            if(data_info.first == -2) // -2 vuol dire che il dato è sull'elemento
            {
                data_to_elems_[data_id] = {static_cast<unsigned>(data_info.second)};
                elem_to_data_[data_info.second].insert(data_id);
            }
            else if(data_info.first == -1) // -1 vuol dire che il dato è su un nodo
            {
                auto conn_elems = connections_.get_node_to_elems(data_info.second); // elementi connessi al dato
                data_to_elems_[data_id] = conn_elems;
                for(unsigned elem_id : conn_elems)
                    elem_to_data_[elem_id].insert(data_id);
            }
            else // ultimo caso: il dato è dentro ad un edge
            {
                auto conn_elems1 = connections_.get_node_to_elems(data_info.first);
                auto conn_elems2 = connections_.get_node_to_elems(data_info.second);
                // intersezione tra i due set
                for (auto it = conn_elems1.begin(); it != conn_elems1.end();) {
                    if (!conn_elems2.count(*it)) { it = conn_elems1.erase(it); }
                    else              { ++it; }
                }
                // vengono ora aggiornati le connessioni
                // assert(conn_elems1.size()==2); devono essere per forza 2?
                data_to_elems_[data_id] = conn_elems1;
                for(unsigned elem_id : conn_elems1) {elem_to_data_[elem_id].insert(data_id);}
            }
            ++i;
        }
        // vengono trovate le facce il cui costo deve essere calcolato di nuovo
        
        compute_costs(facets_to_update, cost_obj);

        // finito il collapse vengono tolti M - 1 nodi
        n_nodes_ = n_nodes_ - (M-1);
        if(n_nodes_ < n_nodes)
        	break;

	}

} // simplify

template<int M, int N>
Mesh<M, N> Simplification<M, N>::build_mesh() const
{
	// vengono presi gli elementi e i nodi attivi
	const auto & active_elems = connections_.get_active_elements();
    const auto & active_nodes = connections_.get_active_nodes();
    DMatrix<double> new_nodes(active_nodes.size(), N);
    DMatrix<int> new_elems(active_elems.size(), Mesh<M, N>::n_vertices);
    DMatrix<int> new_boundary(active_nodes.size(), 1);
    new_boundary.setZero();
    // viene creata la matrice dei nodi
    std::unordered_map<unsigned, unsigned> node_ids_map; // da utilizzare per mappare i vecchi id dei nodi ai nuovi id
    unsigned new_id = 0;
    for(unsigned old_id : active_nodes){
        for(unsigned j = 0; j < N; ++j)
            new_nodes(new_id, j) = nodes_(old_id, j); 
        node_ids_map[old_id] = new_id;
        new_id++;
    }
    // viene creata la matrice degli elementi
    new_id = 0;
    for(unsigned old_id : active_elems){
        for(unsigned j = 0; j < Mesh<M, N>::n_vertices; ++j)
            new_elems(new_id, j) = node_ids_map.at(elems_mat_(old_id, j));
        new_id++;
    }
    return Mesh<M, N>(new_nodes, new_elems, new_boundary);
}

} // fdapde
} // core

#endif // __SIMPLIFICATION_H__