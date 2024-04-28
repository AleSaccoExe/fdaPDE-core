#ifndef __SIMPLIFICATION_H__
#define __SIMPLIFICATION_H__

#include "mesh.h"
#include "StructuredGridSearch.h"
#include "connections.h"
#include <chrono>

namespace fdapde{
namespace core{
template<int M, int N>
class Simplification{
	using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
	using FacetType = std::array<int, Mesh<M, N>::n_vertices_per_facet>;
private:
	double cost_threshold_ = 1.0;

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

    bool does_intersect(unsigned facet_id, SVector<N> collapse_point);
    // calcola il costo per ogni facet
    template<std::size_t K, typename... Args>
    void compute_costs(const std::set<unsigned> & facet_ids, const std::array<double, K> & w, Args&&... cost_obj);
    // trova i possibili punti per il collapse di facet
    std::vector<SVector<N>> get_collapse_points(const FacetType & facet) const;
    // elimina le facce in input
    void erase_facets(const std::unordered_set<unsigned> & facet_ids);
    // modifica le facce in input usando le informazioni in facet
    void modify_facets(const std::unordered_set<unsigned> & facet_ids, const FacetType & facet);
    
    void update_boundary(const FacetType & facet);
    // metodo per sommare i costi con i pesi
    /*
    template<std::size_t K, typename CostType_, typename... CostTypes_>
    double sum_costs(const std::array<double, K>& w, 
    				 const std::vector<Element<M, N>> & elems_to_modify, 
    				 const std::vector<Element<M, N>> & elems_to_erase, 
    				 const std::vector<Element<M, N>> & elems_modified, 
    				 SVector<3> collapse_point, 
    				 const std::unordered_set<unsigned> & data_ids,
    				 CostType_&& cost_obj, CostTypes_&&... cost_objs) const
    {
	std::array<double, K-1> ww;
    std::copy(w.begin()+1, w.end(), ww.begin());
    return  sum_costs(ww, elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids , cost_objs...) 
    		+ w[0]*cost_obj(elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids);
    }
    template<std::size_t K, typename CostType_>
    double sum_costs(const std::array<double, K>& w, 
    				 const std::vector<Element<M, N>> & elems_to_modify, 
    				 const std::vector<Element<M, N>> & elems_to_erase, 
    				 const std::vector<Element<M, N>> & elems_modified, 
    				 SVector<3> collapse_point, 
    				 const std::unordered_set<unsigned> & data_ids,
    				 CostType_&& cost_obj) const
    { return w[0]*cost_obj(elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids); }
    */
     template<std::size_t K>
     double sum_costs(const std::array<double, K>& w, 
    				 const std::vector<Element<M, N>> & elems_to_modify, 
    				 const std::vector<Element<M, N>> & elems_to_erase, 
    				 const std::vector<Element<M, N>> & elems_modified, 
    				 SVector<N> collapse_point, 
    				 const std::unordered_set<unsigned> & data_ids) const
    {
    	return 0.0;
    }
    template<std::size_t K, typename CostType_, typename... CostTypes_>
    double sum_costs(const std::array<double, K>& w, 
    				 const std::vector<Element<M, N>> & elems_to_modify, 
    				 const std::vector<Element<M, N>> & elems_to_erase, 
    				 const std::vector<Element<M, N>> & elems_modified, 
    				 SVector<N> collapse_point, 
    				 const std::unordered_set<unsigned> & data_ids,
    				 CostType_&& cost_obj, CostTypes_&&... cost_objs) const
    {
    	double weighted_cost = cost_obj(elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids)*w[K-sizeof...(cost_objs)-1];
    	double rest = sum_costs(w, elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids , cost_objs...);
    	return  rest + weighted_cost;
    }

    // per gli oggetti costo vengono aggiornati i costi massimi
    template<typename... CostTypes_>
    void update_max_costs(const std::set<unsigned>& facet_ids, CostTypes_&&... cost_objs) const;
public:
	Simplification(const Mesh<M, N> & mesh, const DMatrix<double> & data = DMatrix<double>(0, N));
	template<std::size_t K, typename... Args>
	void simplify(unsigned n_nodes, std::array<double, K> w, Args&&... cost_objs);
	// costruisce la mesh semplificata
	Mesh<M, N> build_mesh() const;
	// temporaneo
	const DMatrix<double> & get_data() const {return data_;}

	void controlli_vari() const;
	std::vector<Element<M, N>> modify_elements(std::unordered_set<unsigned> & elem_ids, const FacetType & facet, SVector<N> new_coord) const;

	// getters
	std::set<unsigned> active_elems() const {return connections_.get_active_elements();}
	std::unordered_set<unsigned> data_to_elems(unsigned datum_id) const {return data_to_elems_[datum_id];}
	std::unordered_set<unsigned> elem_to_data(unsigned elem_id) const {return elem_to_data_[elem_id];}
};

// ===============
// implementazione
// ===============

// costruttore 
template<int M, int N>
Simplification<M, N>::Simplification(const Mesh<M, N> & mesh, const DMatrix<double> & data):
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
    if(data.rows() == 0){
    	std::cout<<"matrice dei dati vuota\n";
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
	else{
		data_ = data;
		data_to_elems_.resize(data_.rows());
		auto locations = mesh.locate(data);
		for(unsigned datum_id = 0; datum_id < data.rows(); ++datum_id)
		{
			unsigned elem_id = locations[datum_id];
			data_to_elems_[datum_id].insert(elem_id);
			elem_to_data_[elem_id].insert(datum_id);
		}
	}
}

template<int M, int N>
std::vector<Element<M, N>> Simplification<M, N>::modify_elements(std::unordered_set<unsigned> & elem_ids, const FacetType & facet, SVector<N> new_coord) const
{
	std::vector<Element<M, N>> elems;
	for(unsigned elem_to_modify : elem_ids)
    {
        const auto & elem_tmp = elems_vec_[elem_to_modify];
        std::array<int, ct_nvertices(M)> node_ids = elem_tmp.node_ids();
        std::array<SVector<N>, ct_nvertices(M)> coords = elem_tmp.coords(); 
        // il vettore degli id viene riempito
        /*for(unsigned i = 0; i < ct_nvertices(M); ++i){
            node_ids[i] = elems_mat_(elem_to_modify, i);
            coords[i] = nodes_.row(node_ids[i]);
        }*/
        for(int i = 0; i < Mesh<M, N>::n_vertices; ++i)
            for(int j = 0; j < Mesh<M, N>::local_dimension; ++j)
                if(facet[j] == node_ids[i]){
                    node_ids[i] = facet[0];
                    coords[i] = new_coord;
        		}
        // si costruisce ora un nuovo elemento con l'id del nodo modificato
        elems.emplace_back(elem_to_modify, node_ids, coords, elem_tmp.neighbors(), elem_tmp.is_on_boundary());
    }
    return elems;

}

// =============================
// implementazione compute_costs
// =============================
template<int M, int N>
template<std::size_t K, typename... Args>
void Simplification<M, N>::compute_costs(const std::set<unsigned> & facet_ids, const std::array<double, K>& w, Args&&... cost_objs)
{
	using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::duration_cast;

	bool max_to_update = false;
	for(unsigned facet_id : facet_ids){
		const auto & facet = facets_[facet_id];
    	std::vector<SVector<N>> collapse_points = get_collapse_points(facet);
    	// viene aggiunto un costo se: sulla faccia ci sono 2 nodi
	    if(connections_.nodes_on_facet(facet).size()==2 && collapse_points.size() != 0)
	    {
	    // per ogni punto in collapse_points si calcola il costo della contrazione
	    std::map<double, SVector<N>> tmp_costs_map;
	    auto elems_to_modify_ids = connections_.elems_modified_in_collapse(facet);
	    auto elems_to_erase_ids = connections_.elems_erased_in_collapse(facet);
	    // viene creato un vettore degli elementi involved nel collapse di facet
        std::vector<Element<M, N>> elems_to_modify;
        elems_to_modify.reserve(elems_to_modify_ids.size());
        for(unsigned elem_id : elems_to_modify_ids) {elems_to_modify.push_back(elems_vec_[elem_id]); }
        std::vector<Element<M, N>> elems_to_erase;
        for(unsigned elem_id : elems_to_erase_ids) {elems_to_erase.push_back(elems_vec_[elem_id]); }
        // vengono presi i dati da proiettare
        std::unordered_set<unsigned> data_ids;
        for(unsigned elem_id : elems_to_erase_ids)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        for(unsigned elem_id : elems_to_modify_ids)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        // vengono salvate le vecchie normali agli elementi
        std::vector<SVector<N>> old_normals;
        old_normals.reserve(elems_to_modify_ids.size());
        for(unsigned elem_id : elems_to_modify_ids) { old_normals.push_back(elems_vec_[elem_id].hyperplane().normal()); }

        auto start_for_collapse_points = Clock::now();
	    for(auto& collapse_point : collapse_points)
	    {
	    	std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify_ids, facet, collapse_point);
	    	// vengono calcolate le nuove normali agli elementi
	    	std::vector<SVector<N>> new_normals;
	    	new_normals.reserve(elems_modified.size());
	    	for(const auto & elem : elems_modified) { new_normals.push_back(elem.hyperplane().normal()); }
	    	// controllo sulla validità del collapse
	    	bool valid_collapse = true;
	    	for(unsigned i = 0; i < elems_modified.size() && valid_collapse; ++i)
	        	{ valid_collapse = valid_collapse && ( elems_modified[i].measure()>DOUBLE_TOLERANCE )
	        								  && ( new_normals[i].dot(old_normals[i]) > DOUBLE_TOLERANCE ); }
	     	if(valid_collapse)
	     	{
		        // calcolato il costo del collapse
		        double cost = sum_costs(w, elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids,
		        						cost_objs...);
		        // double cost = (... + cost_objs(elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids));
		        tmp_costs_map[cost] = collapse_point;
		    }
	    }
	    auto end_for_collapse_points = Clock::now();
	    // std::cout<<"tempo for su collapse_points: "<<duration_cast<duration<double>>(end_for_collapse_points - start_for_collapse_points).count()<<"\n";
	    // ora si prende il costo minore e si controllano le intersezioni
	    // prima dalla classe sgs_ vengono eliminate le informazioni su elementi modificati ed elementi eliminati
	    sgs_.erase_elements(elems_to_modify_ids);
	    sgs_.erase_elements(elems_to_erase_ids);
	    auto start_for_tmp_costs_map = Clock::now();
	    for(auto it = tmp_costs_map.begin(); it != tmp_costs_map.end(); ++it)
	    {
	    	bool valid_collapse = true;
	    	if constexpr(K > 1){
		    	if(it->first > cost_threshold_ || (... || cost_objs.check_update()))
		    	{
		    		// cost_threshold_ = it->first;
		    		cost_threshold_ = std::numeric_limits<double>::max();
		    		(cost_objs.set_threshold(std::numeric_limits<double>::max()), ...);
		    		std::cout<<"costi massimi da aggiornare: "<<it->first<<"\n";
		    		max_to_update = true;
		    		break;
		    	}
		    }
	        // viene creato un vettore degli elementi modificati nel collapse di facet
	        std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify_ids, facet, it->second);
	        // ora passo a sgs il vettore di elementi elems_modified, che è formato di elementi da modificare
	        // La classe provvede quindi a capire le possibili intersezioni date dal collapse di facet
	        // sgs_.update(elems_modified, elems_to_erase_ids, false);
	        int n_controlli = 0;
	        // std::cout<<"elementi modificati: "<<elems_modified.size()<<"\n";
	        /*for(const auto & elem : elems_modified){
	        	unsigned stop;
	            auto elems_to_check = sgs_.get_neighbouring_elements(elem);
	            n_controlli+=elems_to_check.size();
	            // std::cout<<"elementi da controllare per le intersezioni: "<<elems_to_check.size()<<"\n";	           
	            for(unsigned elem_id : elems_to_check)
	        	{ 
	        		auto start_intersezione = Clock::now();
	                if(elems_vec_[elem_id].intersection(elem))
	                {
	                    std::cout<<"trovata intersezione\n";
	                    valid_collapse = false;
	                    break;
	                    // std::cin>>stop;
	                } 
	            }
	            if(!valid_collapse) 
	            	break;
	    	}*/
	    	// std::cout<<"numero controlli: "<<n_controlli<<"\n";
			// std::cout<<"tempo intersezioni: "<<duration_cast<duration<double>>(end_intersezioni - start_intersezioni).count()<<"\n";
	        // se non ci sono intersezioni si procede ad aggiungere le informazioni nelle strutture adeguate
	        if(valid_collapse)
	        {
	            costs_map_.insert({it->first, {facet_id, it->second} });
	            facets_cost_[facet_id] = it->first;
	            break;
	        }

	    }
	    auto end_for_tmp_costs_map = Clock::now();
	    // std::cout<<"tempo for su tmp_costs_map: "<<duration_cast<duration<double>>(end_for_tmp_costs_map - start_for_tmp_costs_map).count()<<"\n";
	    // a questo punto si ripristina la classe sgs al suo stato iniziale
	    sgs_.add_elements(elems_to_modify);
	    // sgs_.update_f(elems_tmp);
	    sgs_.add_elements(elems_to_erase);
	    // se i costi massimi sono da aggiornare non si prosegue con il controllo sulle intersezioni
	    if constexpr(K > 1) { if(max_to_update) {break;} }
	    } // if(connections_.nodes_on_facet(facet).size()==2)
    }
    if constexpr(K > 1){
    if(max_to_update)
	    {
	    	std::set<unsigned> all_facets;
	    	for(const auto & pair : facets_cost_) {all_facets.insert(pair.first);}
	    	all_facets.insert(facet_ids.begin(), facet_ids.end());
	    	facets_cost_.clear();
	    	costs_map_.clear();
	    	// update_max_costs(all_facets, cost_objs...);
	    	update_max_costs(facet_ids, cost_objs...);
	    	compute_costs(all_facets, w, cost_objs...);
	    	std::cout<<"finito l'aggiornamento dei costi\n";
	    	cost_threshold_ = 1.0;
	    	(cost_objs.set_threshold(1.3), ...);

	    }
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
		for (std::size_t j = 0; j < M-1; ++j) {
			SVector<N> col = nodes_.row(facet[j+1]) - nodes_.row(facet[0]);
			J.col(j) = col; 
		}
		SVector<M> barycentric_mid_point;
	    barycentric_mid_point.fill(1.0 / (M + 1));
	    SVector<N> mid_point = J * barycentric_mid_point + static_cast<SVector<N>>(nodes_.row(facet[0]));
		// collapse_points.emplace_back(0.5*(nodes_.row(facet[0]) + nodes_.row(facet[1])));
		collapse_points.push_back(mid_point);
		if constexpr(M == 2 && N == 3)
		{
			auto elem_ids1 = connections_.elems_modified_in_collapse(facet);
			auto elem_ids2 = connections_.elems_erased_in_collapse(facet);
			SVector<10> Q;
			Q.setZero();
			for(unsigned elem_id : elem_ids1)
			{
				SVector<3> n = elems_vec_[elem_id].hyperplane().normal(); // normale all'elemento
				double d = -n.dot(elems_vec_[elem_id].coords()[0]);
				// Construct matrix K
				SVector<10> K;
				K <<  n[0]*n[0], n[0]*n[1], n[0]*n[2], n[0]*d, n[1]*n[1], n[1]*n[2], n[1]*d, n[2]*n[2], n[2]*d, d*d;
				Q += K;
			}
			for(unsigned elem_id : elem_ids2)
			{
				SVector<3> n = elems_vec_[elem_id].hyperplane().normal(); // normale all'elemento
				double d = -n.dot(elems_vec_[elem_id].coords()[0]);
				// Construct matrix K
				SVector<10> K;
				K <<  n[0]*n[0], n[0]*n[1], n[0]*n[2], n[0]*d, n[1]*n[1], n[1]*n[2], n[1]*d, n[2]*n[2], n[2]*d, d*d;
				Q += 2*K;
			}
			Eigen::Matrix3d A;
			A << Q[0], Q[1], Q[2], 
				Q[1], Q[4], Q[5],
				Q[2], Q[5], Q[7];
				
			// Right-hand side
			Eigen::Vector3d b(-Q[3], -Q[6], -Q[8]);
			
			//
			// Solve the linear system
			//
			// Exploit QR decomposition with column pivoting
			// This choice should be a good compromise bewteen
			// performance and accuracy
			
			Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
			
			//
			// Check if the solution exists and return
			//
			// As suggested on the Eigen wiki, to know if the
			// solution exists one may compute the relative
			// a posteriori error and check it is below an
			// user-defined tolerance
			
			auto err = (A*x - b).norm() / b.norm();
			if (err < DOUBLE_TOLERANCE)
				collapse_points.emplace_back(x(0), x(1), x(2));
		}
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
template<typename... CostTypes_>
void Simplification<M, N>::update_max_costs(const std::set<unsigned> & facet_ids, CostTypes_&&... cost_objs) const
{
	for(unsigned facet_id : facet_ids){
		const auto & facet = facets_[facet_id];
    	std::vector<SVector<N>> collapse_points = get_collapse_points(facet);
    	// viene aggiunto un costo se: sulla faccia ci sono 2 nodi
	    if(connections_.nodes_on_facet(facet).size()==2 && collapse_points.size() != 0)
	    {
	    // per ogni punto in collapse_points si calcola il costo della contrazione
	    auto elems_to_modify_ids = connections_.elems_modified_in_collapse(facet);
	    auto elems_to_erase_ids = connections_.elems_erased_in_collapse(facet);
	    // viene creato un vettore degli elementi involved nel collapse di facet
        std::vector<Element<M, N>> elems_to_modify;
        for(unsigned elem_id : elems_to_modify_ids) {elems_to_modify.push_back(elems_vec_[elem_id]); }
        std::vector<Element<M, N>> elems_to_erase;
        for(unsigned elem_id : elems_to_erase_ids) {elems_to_erase.push_back(elems_vec_[elem_id]); }
        // vengono presi i dati da proiettare
        std::unordered_set<unsigned> data_ids;
        for(unsigned elem_id : elems_to_erase_ids)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        for(unsigned elem_id : elems_to_modify_ids)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        std::vector<SVector<N>> old_normals; // da utilizzare per controllare le inversioni degli elementi
   		old_normals.reserve(elems_to_modify_ids.size());
        for(unsigned elem_id : elems_to_modify_ids) { old_normals.push_back(elems_vec_[elem_id].hyperplane().normal()); }
        bool go_update_max = false;
	    for(auto collapse_point : collapse_points)
	    {
	    	std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify_ids, facet, collapse_point);
	    	std::vector<SVector<N>> new_normals; // contiene le normali degli elementi modificati
	    	new_normals.reserve(elems_modified.size());
	    	for(const auto & elem : elems_modified) { new_normals.push_back(elem.hyperplane().normal()); }
	    	bool valid_collapse = true;
	        // controllo sull'area dei triangoli e sulle normali:
	        for(unsigned i = 0; i < elems_modified.size() && valid_collapse; ++i)
	        {
	        	valid_collapse = valid_collapse && ( elems_modified[i].measure()>DOUBLE_TOLERANCE )
	        									&& ( new_normals[i].dot(old_normals[i]) > DOUBLE_TOLERANCE );
	        }
	        // calcolato il costo del collapse
	        if(valid_collapse){
	        	go_update_max = true;
	        	(cost_objs.update_min(elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids),...);
	        }
	    }
	    if(go_update_max) { (cost_objs.update_max(), ...); }

	    } // if(connections_.nodes_on_facet(facet).size()==2)
    }
} // update_max_costs


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
			assert(it != costs_map_.end());
			assert(it->second.first == facet_id);
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

template<int M, int N>
void Simplification<M, N>::update_boundary(const FacetType & facet)
{
	int nodes_on_boundary = 0;
	for(unsigned node_id : facet) { nodes_on_boundary += boundary_(node_id); }
	if(nodes_on_boundary != 0)
		boundary_(facet[0]) = 1;
}

template<int M, int N>
bool Simplification<M, N>::does_intersect(unsigned facet_id, SVector<N> collapse_point)
{
	const auto & facet = facets_[facet_id];
	auto elems_to_modify_ids = connections_.elems_modified_in_collapse(facet);
	auto elems_to_erase_ids = connections_.elems_erased_in_collapse(facet);
	sgs_.erase_elements(elems_to_erase_ids);
	sgs_.erase_elements(elems_to_modify_ids);
	std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify_ids, facet, collapse_point);
	std::vector<Element<M, N>> elems_to_modify;
	for(unsigned elem_to_modify : elems_to_modify_ids) {elems_to_modify.push_back(elems_vec_[elem_to_modify]);}
	std::vector<Element<M, N>> elems_to_erase;
	for(unsigned elem_to_erase : elems_to_erase_ids) {elems_to_erase.push_back(elems_vec_[elem_to_erase]);}
	for(const auto & elem : elems_modified){
        auto elems_to_check = sgs_.get_neighbouring_elements(elem);
        // std::cout<<"elementi da controllare per le intersezioni: "<<elems_to_check.size()<<"\n";	           
        for(unsigned elem_id : elems_to_check)
    	{ 
            if(elems_vec_[elem_id].intersection(elem)){
            	sgs_.add_elements(elems_to_modify);
            	sgs_.add_elements(elems_to_erase);
            	return true; 
            }
        }
	}
	sgs_.add_elements(elems_to_modify);
	sgs_.add_elements(elems_to_erase);
	return false;
}


// ========================
// implementazione simplify
// ========================

template<int M, int N>
template<std::size_t K, typename... Args>
void Simplification<M, N>::simplify(unsigned n_nodes, std::array<double, K> w, Args&&... cost_objs)
{
	using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::duration_cast;

	static_assert(K == sizeof...(cost_objs));
	if(n_nodes >= n_nodes_)
	{
		std::cout<<"required nodes: "<<n_nodes<<", nodes in the mesh: "<<n_nodes_<<std::endl;
		return;
	}
	// setup dei costi
	(cost_objs.setup(this), ...);
	std::cout<<"vengono calcolati i massimi...\n";

	// vengono calcolati i costi per ogni facet valida
	std::set<unsigned> all_facets;
	for(unsigned i = 0; i < facets_.size(); ++i) {all_facets.insert(i);}
	(cost_objs.set_threshold(std::numeric_limits<double>::max()), ...);
	update_max_costs(all_facets, cost_objs...);
	std::cout<<"vengono calcolati i costi...\n";
	compute_costs(all_facets, w, cost_objs...);
	(cost_objs.set_threshold(1.3), ...);
	// viene semplificata la faccia con costo minore
	unsigned numero_semplificazione = 1;
	std::cout<<"inizia la semplificazione\n";
	std::cout<<"costs_map_.size()="<<costs_map_.size()<<"\n";
	while(n_nodes_ > n_nodes && costs_map_.size() > 0)
	{
		auto start_collapse = Clock::now();
		std::cout<<"inizio semplificazione "<<numero_semplificazione<<"\n";
		numero_semplificazione++;
		auto it = costs_map_.begin();
		// auto start_intersezioni = Clock::now();
		/*while(does_intersect(it->second.first, it->second.second)) {
			erase_facets({it->second.first});
			it = costs_map_.begin();
		}*/
		// auto end_intersezioni = Clock::now();
		// std::cout<<"tempo intersezioni: "<<duration_cast<duration<double>>(end_intersezioni - start_intersezioni).count()<<"\n";
		auto collapse_info = *it;
		std::cout<<"facet_id: "<<collapse_info.second.first<<"\n";
		auto facet = facets_[collapse_info.second.first]; // viene presa la faccia da contrarre
		auto elems_to_modify = connections_.elems_modified_in_collapse(facet);
        auto elems_to_erase = connections_.elems_erased_in_collapse(facet);
        // si modificano gli elementi
        // il vettore viene riempito con gli elementi modificati
		std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify, facet, collapse_info.second.second);
		// vengono aggiornate le informazioni della classe sgs_
		// sgs_.update(elems_modified, elems_to_erase, true);
		sgs_.erase_elements(elems_to_modify);
		sgs_.erase_elements(elems_to_erase);
		sgs_.add_elements(elems_modified);
		// vengono presi i dati da proiettare
        std::unordered_set<unsigned> data_ids;
        for(unsigned elem_id : elems_to_erase)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        for(unsigned elem_id : elems_to_modify)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        // vengono proiettati i dati sugli elementi già modificati
        auto start_project = Clock::now();
        auto new_elem_to_data = project(elems_modified, data_, data_ids);
        auto end_project = Clock::now();
        // std::cout<<"tempo project: "<<duration_cast<duration<double>>(end_project-start_project).count()<<"\n";
        // ora si può eseguire il collapse e dopo si cambiano le informazioni di data_to_elems e elem_to_data
        for(const auto & elem : elems_modified){ // loop sugli elementi modificati dalla contrazione
            // viene sostituito nel vettore l'elemento modificato
            elems_vec_[elem.ID()] = elem;
            // nella matrice degli elementi vengono modificati gli id dei nodi
            for(unsigned i = 0; i < ct_nvertices(M); ++i)
                elems_mat_(elem.ID(), i) = elem.node_ids()[i];
        }
       	// vengono cambiate le coordinate del nodo del collapse
       	nodes_.row(facet[0]) = collapse_info.second.second;
       	// elems_tmp.clear(); // questo vettore viene ora riempito di elementi da eliminare
       	std::vector<Element<M, N>> elems_erased;
        for(unsigned elem_id : elems_to_erase) // loop sugli elementi da eliminare
            elems_erased.push_back(elems_vec_[elem_id]);
        // vengono modificate le connessioni e prese le informazioni necessarie a modificare
        // le facce
        auto start_collapse_facet = Clock::now();
        auto facets_pair = connections_.collapse_facet(facet, elems_erased);
        auto end_collapse_facet = Clock::now();
        // std::cout<<"tempo collapse_facet: "<<duration_cast<duration<double>>(end_collapse_facet - start_collapse_facet).count()<<"\n";
        // assert(facets_pair.first.size()==2);
        // vengono eliminate le informazioni sulle facce che non esistono più
        facets_pair.first.insert(collapse_info.second.first);
        erase_facets(facets_pair.first);
        // e anche tutte le altre facce con costi da ricalcolare
        // auto start_facets_to_update = Clock::now();
        auto facets_to_update = connections_.facets_to_update(facet[0]);
        // auto end_facets_to_update = Clock::now();
        // std::cout<<"tempo facets_to_update: "<<duration_cast<duration<double>>(end_facets_to_update - start_facets_to_update).count()<<"\n";
		modify_facets(facets_pair.second, facet);
        erase_facets({facets_to_update.begin(), facets_to_update.end()});
        // vengono modificati gli id dei vertici delle facce
        // vengono ora modificate le connessioni tra dati e elementi
        std::map<unsigned, std::set<unsigned>> new_data_to_elem;
		for(unsigned i = 0; i < new_elem_to_data.size(); ++i)
		{
			auto data_on_elem = new_elem_to_data[i];
			for(unsigned datum_id : data_on_elem) {new_data_to_elem[datum_id].insert(i);}
		}
		for(unsigned i = 0; i < elems_modified.size(); ++i)
		{
			unsigned elem_id = elems_modified[i].ID();
			elem_to_data_[elem_id] = new_elem_to_data[i];
		}
		for(unsigned datum_id : data_ids)
		{
			auto elems_on_datum = data_to_elems_[datum_id];
			for(unsigned elem_id : elems_to_modify) {elems_on_datum.erase(elem_id);}
			for(unsigned elem_id : elems_to_erase)  {elems_on_datum.erase(elem_id);}
			for(unsigned i : new_data_to_elem.at(datum_id)) {elems_on_datum.insert(elems_modified[i].ID());}
			data_to_elems_[datum_id] = elems_on_datum;
		}
		for(unsigned elem_id : elems_to_erase) {elem_to_data_[elem_id].clear();}
		// aggiornata la classe sgs_ (se necessario)
		bool to_refresh = false;
		for(const auto & elem : elems_modified) { to_refresh = to_refresh || sgs_.check_elem_size(elem); }
        if(to_refresh){
        	std::cout<<"sgs da refreshare\n";
        	auto start_refresh = Clock::now();
        	sgs_.refresh(elems_vec_, connections_.get_active_elements(), connections_.get_active_nodes().size());
        	auto end_refresh = Clock::now();
        	std::cout<<"tempo sgs_.refresh : "<<duration_cast<duration<double>>(end_refresh - start_refresh).count()<<"\n";
        }
        // aggiornate le informazioni sul bordo
        update_boundary(facet);
        // aggiornate le funzioni costo:
        (cost_objs.update(elems_erased, elems_modified), ...);
        // aggiornato il costo delle facce
        std::cout<<"facets_to_update.size() = "<<facets_to_update.size()<<"\n";
        auto start_compute_costs = Clock::now();
        compute_costs(facets_to_update, w, cost_objs...);
        auto end_compute_costs = Clock::now();
    	std::cout << "tempo compute_costs: "<<duration_cast<duration<double>>(end_compute_costs - start_compute_costs).count()<<"\n";

        // finito il collapse vengono tolti M - 1 nodi
        n_nodes_ = n_nodes_ - (M-1);
        auto end_collapse = Clock::now();
        std::cout<<"tempo collapse: "<<duration_cast<duration<double>>(end_collapse - start_collapse).count()<<"\n";
	}

} // simplify

template<int M, int N>
Mesh<M, N> Simplification<M, N>::build_mesh() const
{
	controlli_vari();
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

template<int M, int N>
void Simplification<M, N>::controlli_vari() const
{
	auto elems_ids = connections_.get_active_elements();
	for(const auto & elems_on_datum : data_to_elems_)
	{
		for(unsigned elem_id : elems_on_datum)
			if(elems_ids.find(elem_id) == elems_ids.end())
				std::cout<<"a qualche dato sono connessi elementi eliminati\n";
	}
	std::unordered_set<unsigned> all_data;
	for(const auto & data_on_elem : elem_to_data_)
		all_data.insert(data_on_elem.begin(), data_on_elem.end());
	if(all_data.size() != data_.rows())
		std::cout<<"dati trovati nelle connessioni: "<<all_data.size()<<", dati totali: "<<data_.rows()<<"\n";
	std::cout<<"fine controlli\n";
}

} // fdapde
} // core

#endif // __SIMPLIFICATION_H__