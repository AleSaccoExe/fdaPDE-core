#ifndef __STRUCTURED_GRID_SEARCH_H__
#define __STRUCTURED_GRID_SEARCH_H__

#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <array>

#include "../utils/symbols.h"
#include "../utils/intersection.h"
#include "mesh_utils.h"


// This class is used by Simplification. In the case of a 2.5D mesh, the class is used to find 
// possible self intersections among the mesh. However, the implementation is general, and given for 
// any (compatible) values of M, N

namespace fdapde{
namespace core{


// M: local dimension, N: embedding dimension
template<int M, int N>
class StructuredGridSearch{
	private:
		using VectorType = SVector<N>;
		// The type of the element bounding box
		using BoundingBoxType = std::array<VectorType, 2>;
		// Nortg-Est point of the entire mesh
		VectorType global_NE;
		// South-West point of the entire mesh
		VectorType global_SW;
		// This  variable contains the size of the cells (dx, dy, and, eventually, dz)
		std::array<double, N> cell_size = {0.0};
		// Each element in the mesh is assign to an unsigned index indicating its position
		// inside the grid.
		// The key of the map is the index, the value (the std::set) containg the 
		// ids of the elements assigned to that index.
		std::unordered_map<unsigned, std::set<unsigned>> idx_map_;
		// The number of the cells in each directions
		std::array<unsigned, N> n_cells;
		// The vector contains the bounding box of each element in the mesh
		std::vector< BoundingBoxType > bounding_boxes_;
		// Used to check if the bounding boxes and the elements' indeces have to be computed again
		bool to_refresh_ = false;

		// private methods
		
	public:
		// Public methods:

		// constructor
		StructuredGridSearch(const Mesh<M, N> & mesh);

		// Check if the size of the input elements is big enough to set to_refresh_ to true
		bool check_elem_size(const Element<M, N> & elem) const;
		// Compute the position index of the point in input
		unsigned compute_index(const SVector<N> & x) const; 
		// Main method: given the element el, find all elements that might intersect el
		std::unordered_set<unsigned> get_neighbouring_elements(const Element<M, N> & el) const; 
		// Compute the position of the cell containing the input point
		std::array<unsigned, N> compute_coordinates(const SVector<N> & x) const;
		// Update the size of the cells using the size saved in cell_size
		void update_cell_size();
		// Erase from idx_map_ the elements in input
		void erase_elements(const std::unordered_set<unsigned> & el_ids);
		// Add in idx_map_ the elements in input, and update the bounding box in bounding_boxes_
		void add_elements(const std::vector<Element<M, N>> & elems);
		// Temporary method: it is not used by Simplification
		// Updates the information of the input elements 
		// The "f" stands for fast, since the method does not check if the structure has to be refreshed
		void update_f(const std::vector<Element<M, N>> & elements);
		// Temporary method: it is not used by Simplification
		// Erase some elements, update informations of other elements and check if the structure
		// has to be refreshed
		void update(const std::vector<Element<M, N>> & elems_to_update, 
			const std::unordered_set<unsigned> & elems_to_erase, bool refresh_check = false);
		// refresh the structure: compute again cell_size, n_cells, and the index for all elements
		void refresh(const std::vector<Element<M, N>> & elems, const std::set<unsigned> & active_elems, unsigned n_nodes); 


		// getters
		BoundingBoxType get_global_bounding_box() const {return {global_SW, global_NE};}
		const std::array<double, N> & get_cell_size() const {return cell_size;}
		const std::array<unsigned, N> & get_n_cells() const {return n_cells;}
		BoundingBoxType get_bounding_box(unsigned el_id) const;
		bool to_refresh() const {return to_refresh_;}



};

// ======================
// IMPELMENTATIVE DETAILS
// ======================


// constructor
template<int M, int N>
StructuredGridSearch<M, N>::StructuredGridSearch(const Mesh<M, N> & mesh)
{
	unsigned n_nodes = mesh.n_nodes();
	unsigned n_elements = mesh.n_elements();

	// global_NE, global_SW initialized
	global_SW.setZero();
	global_NE.setZero();
	global_SW = mesh.node(0);
	global_NE = mesh.node(0);
	
	// global_NE, global_SW are computed
	for(unsigned i = 0; i < n_nodes; ++i) // loop over all nodes
	{
		for(unsigned j = 0; j < N; ++j) // loop over all coordinates of node i
		{
			if(mesh.node(i)[j] < global_SW[j])
				global_SW[j] = mesh.node(i)[j];
			if(mesh.node(i)[j] > global_NE[j])
				global_NE[j] = mesh.node(i)[j];
		}

	}

	// dx, dy (and dz) are computed

	for(unsigned id_el = 0; id_el < n_elements; ++id_el) // loop over all elements
	{
		for(unsigned v1 = 0; v1 < ct_nvertices(M); ++v1) // loop over all edges (v1, v2)
		{
			for(unsigned v2 = v1 + 1; v2 < ct_nvertices(M); ++v2)
			{
				unsigned id_v1 = mesh.elements()(id_el, v1);
				unsigned id_v2 = mesh.elements()(id_el, v2);
				for(unsigned k = 0; k < N; ++k)
				{
					double coord_v1 =  mesh.node( id_v1 )[k];
					double coord_v2 =  mesh.node( id_v2 )[k];
					double diff = std::abs( coord_v1 - coord_v2 );
					if(diff > cell_size[k])
						cell_size[k] = diff;
				}
				
			}
		}	
	}
	/*std::cout<<"numero di celle"<<std::endl;
	for(unsigned i = 0; i < N; ++i){
		n_cells[i] = ( global_NE[i] - global_SW[i] )/cell_size[i];
		std::cout<<n_cells[i]<<", ";
	}
	std::cout<<std::endl;*/
	// cell_size is updated using the computed number of cells
	update_cell_size();
	// std::cout<<"grandezza celle:"<<std::endl<<cell_size[0]<<"  "<<cell_size[1]<<"  "<<cell_size[2]<<std::endl;

	// compute index for all elements
	for(unsigned id_el = 0; id_el < n_elements; ++id_el) // Loop over the elements
	{
		// Take the element's bounding box:
		auto bounding_box = mesh.element(id_el).bounding_box();
		// Add it to the vector of all bounding boxes:
		bounding_boxes_.push_back({bounding_box.first, bounding_box.second});
		// Take the middle point of the bounding box:
		SVector<N> middle_point =  ( bounding_box.first + bounding_box.second )*0.5;
		// Compute the index and add it to the map:
		idx_map_[compute_index(middle_point)].insert(id_el);
	}
}

template<int M, int N>
void StructuredGridSearch<M, N>::update_f(const std::vector<Element<M, N>> & elements)
{
	for (const auto & element : elements)
	{
		unsigned el_id = element.ID();
		//
		// Remove old bounding box
		//
		
		// Extract element index
		// if(boxes_map.find(el_id) != boxes_map.end())
		{
			// auto bounding_box = boxes_map.at(el_id); // get the bounding box of the element
			const auto& bounding_box = bounding_boxes_[el_id];
			// compute the middle point of the bb to get the old index of the element
			SVector<N> middle_point = 0.5*(bounding_box[0] + bounding_box[1]);
			unsigned idx = compute_index(middle_point); // old index of the element
			// the bounding box of the element is erased and then computed again
			// boxes_map.erase(el_id);
			// the element is erased
			// se è tutto giusto non c'e bisogno di un controllo
			// assert( idx_map.at(idx).find(el_id) != idx_map.at(idx).end() );
			idx_map_.at(idx).erase(el_id);
			// now the new bounding box is computed
			auto new_bb = element.bounding_box();
			SVector<N> new_middle_point = 0.5*(new_bb.first + new_bb.second);
			bounding_boxes_[el_id] = new_bb;
			unsigned new_idx = compute_index(new_middle_point); // compute the new index
			idx_map_[new_idx].insert(el_id); // insert it in the map

		}
		
	}
}

template<int M, int N>
typename StructuredGridSearch<M, N>::BoundingBoxType StructuredGridSearch<M, N>::get_bounding_box(unsigned el_id) const
{
	return bounding_boxes_[el_id];
}

template<int M, int N>
void StructuredGridSearch<M, N>::erase_elements(const std::unordered_set<unsigned> & el_ids)
{
	// Loop over the input element ids
	for(unsigned el_id : el_ids)
		{
			// Take the current element's bounding box:
			const auto& bounding_box = bounding_boxes_[el_id];
			SVector<N> middle_point = 0.5*( bounding_box[0] + bounding_box[1] );
			// Compute the index of the element:
			unsigned el_idx = compute_index(middle_point);
			assert(idx_map_[el_idx].find(el_id) != idx_map_[el_idx].end());
			// Erase from the map the id of the current element:
			idx_map_[el_idx].erase(el_id);
		}
}

template<int M, int N>
void StructuredGridSearch<M, N>::update_cell_size()
{
	for(unsigned i = 0; i < N; ++i)
		cell_size[i] = ( global_NE[i] - global_SW[i] )/n_cells[i];
}


template<int M, int N>
std::array<unsigned, N> StructuredGridSearch<M, N>::compute_coordinates(const VectorType & x) const
{
	std::array<unsigned, N> coordinates;
	for(unsigned j = 0; j < N; ++j)
		coordinates[j] = static_cast<unsigned>( floor( (x[j] - global_SW[j] )/cell_size[j] ) );
	return coordinates;
}


template<int M, int N>
unsigned StructuredGridSearch<M, N>::compute_index(const VectorType & x) const
{
	unsigned idx = 0;
	for(unsigned j = 0; j < N; ++j)
	{
		unsigned m = 1;
		for(unsigned i = 0; i < j; ++i)
			m*=n_cells[i];
		idx += m*static_cast<unsigned>( floor( (x[j] - global_SW[j] )/cell_size[j] ) );
	}
	return idx;
}


template<int M, int N>
bool StructuredGridSearch<M, N>::check_elem_size(const Element<M, N> & elem) const
{
	double toll = 1.5;
	// loop over all edges of the element
	for(unsigned i = 0; i < Element<M, N>::n_vertices; ++i)
		for(unsigned j = i+1; j < Element<M, N>::n_vertices; ++j)
		{
			SVector<N> v_i_coords = elem.coords()[i];
			SVector<N> v_j_coords = elem.coords()[j];
			for(unsigned k = 0; k < N; ++k) // loop over all coordinates
				if( std::abs(v_j_coords[k] - v_i_coords[k]) > toll*cell_size[k] )
					return true; // return true if an edge is too long
		}
	// The element is not too big. Return false
	return false;
}

template<int M, int N>
void StructuredGridSearch<M, N>::update(const std::vector<Element<M, N>> & elems_to_update, 
	const std::unordered_set<unsigned> & elems_to_erase, bool refresh_check)
{
	erase_elements({elems_to_erase.begin(), elems_to_erase.end()});
	update_f(elems_to_update);
	// check if the elements are getting too big
	if(refresh_check)
		for(const auto & elem : elems_to_update)
			to_refresh_ = to_refresh_ || check_elem_size(elem);

}

template<int M, int N>
void StructuredGridSearch<M, N>::add_elements(const std::vector<Element<M, N>> & elems)
{
	for(const auto & elem : elems)
	{
		const auto& bounding_box = elem.bounding_box(); 
		auto middle_point = 0.5*(bounding_box.first + bounding_box.second);
		unsigned idx = compute_index(middle_point); // compute the index of the element
		idx_map_[idx].insert(elem.ID());
		bounding_boxes_[elem.ID()] = {bounding_box.first, bounding_box.second};
	}
}

template<int M, int N>
void StructuredGridSearch<M, N>::refresh(const std::vector<Element<M, N>> & elems, 
	const std::set<unsigned> & active_elems, unsigned n_nodes)
{
	to_refresh_ = false;
		
	unsigned n_elements = active_elems.size();
	// global_NE, global_SW initialized
	global_SW.setZero();
	global_NE.setZero();
	unsigned tmp_id = *active_elems.begin();
	global_NE = elems[tmp_id].coords()[0];
	global_SW = elems[tmp_id].coords()[0];
	// global_NE, global_SW are computed
	std::unordered_set<unsigned> visited_nodes;
	for(unsigned i = 0; i < N; ++i) {cell_size[i] = 0.;}
	// Loop over (active) elements
	for(unsigned id_elem : active_elems) 
	{
		const auto & node_ids = elems[id_elem].node_ids(); // get the ids of the element's vertices
		for(unsigned i = 0; i < Element<M, N>::n_vertices; ++i) // loop over all nodes of the element
		{
			if(visited_nodes.find(node_ids[i]) == visited_nodes.end()) // if the node was not visited
			{
				// visited_nodes.insert(node_ids[i]);
				const auto & coords = elems[id_elem].coords()[i]; // get the coords of the i-th vertex of element id_elem
				for(unsigned k = 0; k < N; ++k) // loop over dimensions of the point
				{
					if(coords[k] < global_SW[k]) {global_SW[k] = coords[k];} // update global_SW if necessary
					if(coords[k] > global_NE[k]) {global_NE[k] = coords[k];} // update global_NE if necessary
				}
			}
		}
	}
	// dx, dy (and dz) are computed
	for(unsigned id_el : active_elems) // loop over all active elements
	{
		for(unsigned v1 = 0; v1 < ct_nvertices(M); ++v1) // loop over all edges (v1, v2)
		{
			for(unsigned v2 = v1 + 1; v2 < ct_nvertices(M); ++v2)
			{

				auto coords_v1 = elems[id_el].coords()[v1];
				auto coords_v2 = elems[id_el].coords()[v2];
				for(unsigned k = 0; k < N; ++k){ // loop over dimension of the points
					double diff = std::abs(coords_v1[k] - coords_v2[k]);
					if( diff > cell_size[k] ) // update cell_size if necessary
						cell_size[k] = diff;
				}
			}
		}
	}
	/*std::cout<<"numero di celle:"<<std::endl;
	for(unsigned k = 0; k < N; ++k) { 
		n_cells[k] = ( global_NE[k] - global_SW[k] )/cell_size[k]; 
		if(n_cells[k]==0) {n_cells[k]=1;}
		std::cout<<n_cells[k]<<", ";
	}
	std::cout<<std::endl;*/
	// cell_size is updated using the computed number of cells
	update_cell_size();
	// std::cout<<"grandezza celle:"<<std::endl<<cell_size[0]<<"  "<<cell_size[1]<<"  "<<cell_size[2]<<std::endl;
	// idx_map has to be rebuilt from scratch
	idx_map_.clear();
	// compute index for all elements
	// bounding_boxes_.clear();
	for(unsigned id_el : active_elems) 
	{
		auto bounding_box = elems[id_el].bounding_box();
		// bounding_boxes_.push_back({bounding_box.first, bounding_box.second});
		bounding_boxes_[id_el] = {bounding_box.first, bounding_box.second};
		SVector<N> middle_point =  ( bounding_box.first + bounding_box.second )*0.5;
		idx_map_[compute_index(middle_point)].insert(id_el);
	}
}

template<int M, int N>
std::unordered_set<unsigned> StructuredGridSearch<M, N>::get_neighbouring_elements(const Element<M, N> & el) const
{
	using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using namespace std::chrono;
    // This method is used only in the case of a 2.5D mesh
	if constexpr((M == 2 && N == 2) || M == 3) {
       return {};
	} else {
	// compute the middle point of the element's bounding box
	const auto& bounding_box = el.bounding_box();
	BoundingBoxType bb_element = {bounding_box.first, bounding_box.second};
	// compute the coordinates of the SW and NE points of the bb
	std::array<unsigned, N> coordinates_NE = compute_coordinates(bb_element[1]);
	std::array<unsigned, N> coordinates_SW = compute_coordinates(bb_element[0]);

		// Find intersecting boxes
		//
		// Scan all the cells intersecting the reference box
		// plus an extra layer. This should ensure that all
		// possible intersecting elements are taken into account
		// since an element cannot span more than one cell.
					
		std::unordered_set<unsigned> res;
		
		// Determine indices range
		unsigned i_start = (coordinates_SW[0] == 0) ? 0 : coordinates_SW[0]-1;
		unsigned i_stop = coordinates_NE[0]+1;
		unsigned j_start = (coordinates_SW[1] == 0) ? 0 : coordinates_SW[1]-1;
		unsigned j_stop = coordinates_NE[1]+1;
		unsigned k_start = (coordinates_SW[2] == 0) ? 0 : coordinates_SW[2]-1;
		unsigned k_stop = coordinates_NE[2]+1;


		for (unsigned i = i_start; i <= i_stop; ++i)
			for (unsigned j = j_start; j <= j_stop; ++j)
				for (unsigned k = k_start; k <= k_stop; ++k)
				{
					// Compute scalar index
					unsigned idx = i + j * n_cells[0] + 
						k * n_cells[0] * n_cells[1];
												
					// Extract all boxes having that index
					auto it = idx_map_.find(idx);
					if(it != idx_map_.end())
					{
						std::set<unsigned> range = it->second;
						for(unsigned id : range){
							// check if the bounding boxes intersect
							if(boxes_intersection(bounding_boxes_[id], bb_element)) 
								res.insert(id);
						}
					}
				}
		
		// Remove Id of reference element from the set,
		// then convert to a vector
		res.erase(el.ID());
		return res;
	}
}


}
}

#endif // __BOUNDING_BOXES_H__