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
	// private members:

	// used to check for self-intersections:
	bool check_intersections_ = false;
	// used to check for self-intersections:
	StructuredGridSearch<M, N> sgs_;
	// contains the connectivity and changes it when needed:
	Connections connections_;
	// contains all the elements in the mesh:
	std::vector<Element<M, N>> elems_vec_;
	// contains the coordinates of the mesh' nodes:
	DMatrix<double> nodes_;
	DMatrix<int> elems_mat_;
	// contains all the facets in the mesh (defined by the ids of the vertices):
	std::vector<FacetType> facets_;
	// specifies whether a node is on the boundary:
	DMatrix<int> boundary_;
	// contains the coordinates of all data points:
	DMatrix<double> data_;
	// contains for each datum the elements it lies on: 
	ConnectionsType data_to_elems_;
	// contains for each element the data points that lie on it:
    ConnectionsType elem_to_data_;
   	// the facets are ordered base on their collapsing cost:
    std::multimap<double, std::pair<unsigned, SVector<N> >> costs_map_;
    // stores for each facet the cost of its collapse:
    std::unordered_map<unsigned, double> facets_cost_;
    // stores the number of nodes in the current mesh:
    unsigned n_nodes_;

    // private methods:

    bool does_intersect(unsigned facet_id, SVector<N> collapse_point);
    // computes the collapsing cost for each facet in input :
    template<std::size_t K, typename... Args>
    void compute_costs(const std::set<unsigned> & facet_ids, const std::array<double, K> & w, Args&&... cost_obj);
    // finds all possible collapsing points for a given facet:
    std::vector<SVector<N>> get_collapse_points(const FacetType & facet) const;
    // deletes the facets in input:
    void erase_facets(const std::unordered_set<unsigned> & facet_ids);
    // modifies the facets in input:
    void modify_facets(const std::unordered_set<unsigned> & facet_ids, const FacetType & facet);
    // updates the boundary assuming the facet in input has to be collapsed:
    void update_boundary(const FacetType & facet);
    // sum all the normalized cost functions
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
    // sum all the normalized cost functions
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

    // computes the maxiumum values of the cost functions
    template<typename... CostTypes_>
    void update_max_costs(const std::set<unsigned>& facet_ids, CostTypes_&&... cost_objs) const;
public:

	// constructor
	Simplification(const Mesh<M, N> & mesh, const DMatrix<double> & data = DMatrix<double>(0, N));
	// set the check for self-intersections
	void set_check_intersections(bool check_intersections) {check_intersections_ = check_intersections;}
	// the simplify method
	template<std::size_t K, typename... Args>
	void simplify(unsigned n_nodes, std::array<double, K> w, Args&&... cost_objs);
	// simplify method to be used in the case of only one cost function:
	template<typename CostObj>
	void simplify(unsigned n_nodes, CostObj&& cost_obj)
	{
		std::array<double, 1> w = {0.5};
		cost_obj.max_ = 1.0;
		this->simplify(n_nodes, w, std::forward<CostObj>(cost_obj));
	}
	// build the simplified mesh
	Mesh<M, N> build_mesh() const;

	std::vector<Element<M, N>> modify_elements(std::unordered_set<unsigned> & elem_ids, const FacetType & facet, SVector<N> new_coord) const;

	// getters
	std::set<unsigned> active_elems() const {return connections_.get_active_elements();}
	std::unordered_set<unsigned> data_to_elems(unsigned datum_id) const {return data_to_elems_[datum_id];}
	std::unordered_set<unsigned> elem_to_data(unsigned elem_id) const {return elem_to_data_[elem_id];}
	const DMatrix<double> & get_data() const {return data_;}
};

// ======================
// IMPLEMENTATIVE DETAILS
// ======================

// costruttore 
template<int M, int N>
Simplification<M, N>::Simplification(const Mesh<M, N> & mesh, const DMatrix<double> & data):
	sgs_(mesh), connections_(mesh), elems_mat_(mesh.elements()), nodes_(mesh.nodes()),
	elem_to_data_(elems_mat_.rows()), n_nodes_(mesh.n_nodes()), boundary_(mesh.boundary())
{
	// The vector of elements is built:
	elems_vec_.reserve(mesh.n_elements());
	for(unsigned elem_id = 0; elem_id < mesh.n_elements(); ++elem_id)
		elems_vec_.push_back(mesh.element(elem_id));
	// All the facets are built:
	facets_.reserve(mesh.n_facets());
    for(unsigned facet_id = 0; facet_id < mesh.n_facets(); ++facet_id)
        facets_.push_back(mesh.facet(facet_id).node_ids());
    // The connections between data and elements are initialized:
    if(data.rows() == 0){ // If no data in input: assume the data to lie on the nodes
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
	else{ // if data positions are passed in input...
		data_ = data;
		data_to_elems_.resize(data_.rows());
		auto locations = mesh.locate(data);
		// ... find for each datum the element it lies on:
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
	// Fill the vector with the elements in input modified by the collapse of the input 
	// facet into the point given by new_coord
	std::vector<Element<M, N>> elems;
	for(unsigned elem_to_modify : elem_ids)
    {
        const auto & elem_tmp = elems_vec_[elem_to_modify];
        std::array<int, ct_nvertices(M)> node_ids = elem_tmp.node_ids();
        std::array<SVector<N>, ct_nvertices(M)> coords = elem_tmp.coords(); 
  
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

// ===============================
// implementation of compute_costs
// ===============================
template<int M, int N>
template<std::size_t K, typename... Args>
void Simplification<M, N>::compute_costs(const std::set<unsigned> & facet_ids, const std::array<double, K>& w, Args&&... cost_objs)
{
	bool max_to_update = false;
	// Compute the collapsing cost for all the facets in input:
	for(unsigned facet_id : facet_ids){
		const auto & facet = facets_[facet_id];
		// Get the list of possible collapsing points:
    	std::vector<SVector<N>> collapse_points = get_collapse_points(facet);
    	// The cost is computed if only two nodes insist on the current facet
	    if(connections_.nodes_on_facet(facet).size()==2 && collapse_points.size() != 0)
	    {
	    // Compute the collapsing cost for each point in collapse_points
	    std::map<double, SVector<N>> tmp_costs_map;
	    auto elems_to_modify_ids = connections_.elems_modified_in_collapse(facet);
	    auto elems_to_erase_ids = connections_.elems_erased_in_collapse(facet);
	    // A vector containing all the involved elements in the collpase is created:
        std::vector<Element<M, N>> elems_to_modify;
        elems_to_modify.reserve(elems_to_modify_ids.size());
        for(unsigned elem_id : elems_to_modify_ids) {elems_to_modify.push_back(elems_vec_[elem_id]); }
        std::vector<Element<M, N>> elems_to_erase;
        for(unsigned elem_id : elems_to_erase_ids) {elems_to_erase.push_back(elems_vec_[elem_id]); }
        // The data to project are selected:
        std::unordered_set<unsigned> data_ids;
        for(unsigned elem_id : elems_to_erase_ids)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        for(unsigned elem_id : elems_to_modify_ids)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        // The normal unit vectors for the element before the collapse are saved:
        std::vector<SVector<N>> old_normals;
        old_normals.reserve(elems_to_modify_ids.size());
        for(unsigned elem_id : elems_to_modify_ids) { old_normals.push_back(elems_vec_[elem_id].hyperplane().normal()); }
		// For each point in collapse_points check the validity of the collapse:
	    for(int i = 0; i < collapse_points.size(); ++i)
	    {
	    	const SVector<N>& collapse_point = collapse_points[i];
	    	// Create a vector of the elements obtained after the collapse:
	    	std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify_ids, facet, collapse_point);
	    	// Compute the normal unit vectors of these elements:
	    	std::vector<SVector<N>> new_normals;
	    	new_normals.reserve(elems_modified.size());
	    	for(const auto & elem : elems_modified) { 
	    		SVector<N> new_normal = elem.hyperplane().normal(); 
	    		new_normals.push_back(new_normal); 
	    	}
	    	// Check if inverted elements are created by the collapse:
	    	bool valid_collapse = true;
	    	for(unsigned i = 0; i < elems_modified.size() && valid_collapse; ++i)
	        	{ valid_collapse = valid_collapse && ( elems_modified[i].measure()>100.0*DOUBLE_TOLERANCE )
	        								  && ( new_normals[i].dot(old_normals[i]) > DOUBLE_TOLERANCE ); }
	        // If the collapse is valid, compute and save the collapsing cost:
	     	if(valid_collapse)
	     	{
		        double cost = sum_costs(w, elems_to_modify, elems_to_erase, elems_modified, collapse_point, data_ids,
		        						cost_objs...);
		        tmp_costs_map[cost] = collapse_point;
		    }
	    }
	    // Take the minimum cost and check self-intersections (if needed)
	    // prima dalla classe sgs_ vengono eliminate le informazioni su elementi modificati ed elementi eliminati
	    sgs_.erase_elements(elems_to_modify_ids);
	    sgs_.erase_elements(elems_to_erase_ids);
	    for(auto it = tmp_costs_map.begin(); it != tmp_costs_map.end(); ++it)
	    {
	    	bool valid_collapse = true;

	    	int n_controlli = 0;
	        std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify_ids, facet, it->second);
	        // check intersections only for surface meshes
	        if constexpr(M == 2 && N == 3)
	        	if(check_intersections_)
	        	{
	        		// loop over the modified elements:
	        		for(const auto & elem : elems_modified){
			        	unsigned stop;
			        	// Get the smallest set that might intersect the current element:
			            auto elems_to_check = sgs_.get_neighbouring_elements(elem);
			            // check if any intersection if detected:
			            for(auto elem_id = elems_to_check.begin(); elem_id != elems_to_check.end() && valid_collapse; ++elem_id)
			        	{ 
			        		valid_collapse = valid_collapse && !(elems_vec_[*elem_id].intersection(elem));
			            }
			            // If an intersection was detected, the collapse is not valid
			            if(!valid_collapse){break;}
			    	}
	        } // if(check_intersections_)
	        // If more then one cost function is used, check if the maximum values have to be computed again:
	    	if constexpr(K > 1){
		    	if(valid_collapse && (... || cost_objs.check_update()))
		    	{
		    		(cost_objs.set_threshold(std::numeric_limits<double>::max()), ...);
		    		max_to_update = true;
		    		break;
		    	}
		    }

	        // ora passo a sgs il vettore di elementi elems_modified, che è formato di elementi da modificare
	        // La classe provvede quindi a capire le possibili intersezioni date dal collapse di facet
	        // sgs_.update(elems_modified, elems_to_erase_ids, false);

	        // If the collapse is actually valid, add the necessary informations to costs_map_ and facets_cost_,
	        // and exit the loop.
	        if(valid_collapse)
	        {
	            costs_map_.insert({it->first, {facet_id, it->second} });
	            facets_cost_[facet_id] = it->first;
	            break;
	        }

	    }
	    sgs_.add_elements(elems_to_modify);
	    sgs_.add_elements(elems_to_erase);
	    // If the maximum costs have to be updated, the excecution of compute_costs is stopped:
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
	    	(cost_objs.set_threshold(1.8), ...);

	    }
	}	
} // compute_costs

template<int M, int N>
std::vector<SVector<N>> Simplification<M, N>::get_collapse_points(const FacetType & facet) const
{
	int boundary_nodes = 0;
	// Count the number of nodes on the boundary:
	for(unsigned node_id : facet) 
		boundary_nodes += boundary_(node_id);
	// If the facet has more then one vertex on the boundary, 
	// it can not be collapsed, therefore return an empty vector
	if(boundary_nodes > 1)
		return {};
	std::vector<SVector<N>> collapse_points;
	// If the facet has no nodes on the boundary, the possible collapse points are:
	// the facet's vertices and its middle point (or its center of gravity in the case of a 3D mesh)
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
	// The facet has only one vertex on the boundary: it is the only possible collapsing point
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
	// Implementation similar to compute_costs, but this method does not save 
	// any information of the collapsing costs.
	for(unsigned facet_id : facet_ids){
		const auto & facet = facets_[facet_id];
    	std::vector<SVector<N>> collapse_points = get_collapse_points(facet);
	    if(connections_.nodes_on_facet(facet).size()==2 && collapse_points.size() != 0)
	    {
	    // For each point in collapse_points, the cost of the contraction is computed
	    auto elems_to_modify_ids = connections_.elems_modified_in_collapse(facet);
	    auto elems_to_erase_ids = connections_.elems_erased_in_collapse(facet);
	    // A vector containing the elements involved in the collapse is creatd
        std::vector<Element<M, N>> elems_to_modify;
        for(unsigned elem_id : elems_to_modify_ids) {elems_to_modify.push_back(elems_vec_[elem_id]); }
        std::vector<Element<M, N>> elems_to_erase;
        for(unsigned elem_id : elems_to_erase_ids) {elems_to_erase.push_back(elems_vec_[elem_id]); }
        // Take the data to project:
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
	        // Check on the normals
	        for(unsigned i = 0; i < elems_modified.size() && valid_collapse; ++i)
	        {
	        	valid_collapse = valid_collapse && ( elems_modified[i].measure()>100.0*DOUBLE_TOLERANCE )
	        									&& ( new_normals[i].dot(old_normals[i]) > DOUBLE_TOLERANCE );
	        }
	        // Now, after the check is passed, the cost of the collapse can be computed
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
	// Loop over the facets in input:
	for(unsigned facet_id : facet_ids){
		if(facets_cost_.find(facet_id) != facets_cost_.end()){
			double cost = facets_cost_.at(facet_id);
			auto it = costs_map_.find(cost);
			// assert(it != costs_map_.end());
			while(it->second.first != facet_id)
				++it;
			// assert(it != costs_map_.end());
			// assert(it->second.first == facet_id);
			costs_map_.erase(it);
			facets_cost_.erase(facet_id);
		}
	}
}

template<int M, int N>
void Simplification<M, N>::modify_facets(const std::unordered_set<unsigned> & facet_ids, const FacetType & facet)
{
	for(unsigned facet_id : facet_ids) // Loop over all the facets to modify
    {
        auto & facet_to_modify = facets_[facet_id];
        // The current facet share a node with the facet in input.
        // Find the node and modify its id in facet[0] (the collapsing point)
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


// ==========================
// implementation of simplify
// ==========================

template<int M, int N>
template<std::size_t K, typename... Args>
void Simplification<M, N>::simplify(unsigned n_nodes, std::array<double, K> w, Args&&... cost_objs)
{
	static_assert(K == sizeof...(cost_objs));
	if(n_nodes >= n_nodes_)
	{
		#ifdef __TEST__
			std::cout<<"required nodes: "<<n_nodes<<", nodes in the mesh: "<<n_nodes_<<std::endl;
		#endif
		return;
	}
	// Setup of the costs functions:
	(cost_objs.setup(this), ...);
	#ifdef __TEST__
		std::cout<<"Computing maximum costs...\n";
	#endif

	// Compute the collapsing cost for all the facets:
	std::set<unsigned> all_facets;
	for(unsigned i = 0; i < facets_.size(); ++i) {all_facets.insert(i);}
	(cost_objs.set_threshold(std::numeric_limits<double>::max()), ...);
	if constexpr(sizeof...(cost_objs)>1) { update_max_costs(all_facets, cost_objs...); }
	#ifdef __TEST__
		std::cout<<"Computing collapsing costs...\n";
	#endif
	(cost_objs.set_threshold(std::numeric_limits<double>::max()), ...);
	compute_costs(all_facets, w, cost_objs...);
	(cost_objs.set_threshold(1.5), ...);
	unsigned numero_semplificazione = 1;
	#ifdef __TEST__
		std::cout<<"Starting simplification...\n";
		unsigned n0 = n_nodes_;
		unsigned nf = n_nodes;
		double m = -1./static_cast<double>(n0-nf);
		double q = +1.*n0/static_cast<double>(n0-nf);
		int barWidth = 70; // Larghezza della barra di progresso
	#endif
	// The simplification continues until the nodes in the current mesh is 
	// higher then the nodes required
	while(n_nodes_ > n_nodes && costs_map_.size() > 0)
	{
		#ifdef __TEST__
			double progress = m*n_nodes_+q;
			std::cout << "Simplification process        [";
			unsigned pos(barWidth * progress);
			for (unsigned i = 0; i < barWidth; ++i) 
			{
				if (i < pos) 
					std::cout << "=";
				else if (i == pos) 
					std::cout << ">";
				else 
					std::cout << " ";
			}
			std::cout << "] " << unsigned(progress * 100.0) << " %\r";
			std::cout.flush();
		#endif
		// Take the first element in costs_map_, i.e. the facets with lower collapsing cost
		auto it = costs_map_.begin();
		auto collapse_info = *it;
		// In collapse_info.second.first it is contained the id of the facet to simplify:
		auto facet = facets_[collapse_info.second.first]; 
		// Find the elements erased by the collapse of the current facet:
        auto elems_to_erase = connections_.elems_erased_in_collapse(facet);
        // Build a vector containing the modified elements:
        auto elems_to_modify = connections_.elems_modified_in_collapse(facet);
		std::vector<Element<M, N>> elems_modified = modify_elements(elems_to_modify, facet, collapse_info.second.second);
		// Erase and add the proper element in sgs_ in order to update it:
		// sgs_.update(elems_modified, elems_to_erase, true);
		sgs_.erase_elements(elems_to_modify);
		sgs_.erase_elements(elems_to_erase);
		sgs_.add_elements(elems_modified);
		// Take the data that lie in the involved elements:
        std::unordered_set<unsigned> data_ids;
        for(unsigned elem_id : elems_to_erase)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        for(unsigned elem_id : elems_to_modify)
            data_ids.insert(elem_to_data_[elem_id].begin(), elem_to_data_[elem_id].end());
        // Project the data:
        auto new_elem_to_data = project(elems_modified, data_, data_ids);
        // Now the actual collapse can be done:
        for(const auto & elem : elems_modified){ // Loop on the elements modified by the collapse
            // Replace the modified element in the vector elems_vec_
            elems_vec_[elem.ID()] = elem;
            // In the elemets matrix, modify the ids of the vertices:
            for(unsigned i = 0; i < ct_nvertices(M); ++i)
                elems_mat_(elem.ID(), i) = elem.node_ids()[i];
        }
       	// Modify the coordinates of the collapsing node:
       	nodes_.row(facet[0]) = collapse_info.second.second;
       	// Crete a vector containing the elements erased by the collapse:
       	std::vector<Element<M, N>> elems_erased;
        for(unsigned elem_id : elems_to_erase) // loop sugli elementi da eliminare
            elems_erased.push_back(elems_vec_[elem_id]);
        // Now the connectivity can be properly updated:
        auto facets_pair = connections_.collapse_facet(facet, elems_erased);
        // Erase the facets that no longer exist after the collapse:
        facets_pair.first.insert(collapse_info.second.first);
        erase_facets(facets_pair.first);
        // Find the facets whose cost has to be computed again:
        auto facets_to_update = connections_.facets_to_update(facet[0]);
        // Update the ids of the facets' vertices:
		modify_facets(facets_pair.second, facet);
		// Erase the costs of the facets whose cost has to be computed again:
        erase_facets({facets_to_update.begin(), facets_to_update.end()});
        // Update the connections between data and elemetns:
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
		// check if the member sgs_ has to be updated:
		bool to_refresh = false;
		for(const auto & elem : elems_modified) { to_refresh = to_refresh || sgs_.check_elem_size(elem); }
        if(to_refresh){
        	sgs_.refresh(elems_vec_, connections_.get_active_elements(), connections_.get_active_nodes().size());
        }
        // Some nodes might now be on the boundary:
        update_boundary(facet);
        // Update the cost functions:
        (cost_objs.update(elems_erased, elems_modified), ...);
        // Compute again the cost for some facets:
        compute_costs(facets_to_update, w, cost_objs...);

        // Now the collapse is done. Update the number of nodes in the current mesh:
        n_nodes_ = n_nodes_ - (M-1);
	}

} // simplify

template<int M, int N>
Mesh<M, N> Simplification<M, N>::build_mesh() const
{
	const auto & active_elems = connections_.get_active_elements();
    const auto & active_nodes = connections_.get_active_nodes();
    DMatrix<double> new_nodes(active_nodes.size(), N);
    DMatrix<int> new_elems(active_elems.size(), Mesh<M, N>::n_vertices);
    DMatrix<int> new_boundary(active_nodes.size(), 1);
    new_boundary.setZero();
    // creates node matrix
    std::unordered_map<unsigned, unsigned> node_ids_map; // used to map the old node ids to the new node ids
    unsigned new_id = 0;
    for(unsigned old_id : active_nodes){
        for(unsigned j = 0; j < N; ++j)
            new_nodes(new_id, j) = nodes_(old_id, j); 
        node_ids_map[old_id] = new_id;
        new_id++;
    }
    // create element matrix
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