#ifndef __BOUNDING_BOXES_H__
#define __BOUNDING_BOXES_H__

#include <cmath>
#include <unordered_set>

#include "../utils/symbols.h"
#include "../utils/intersection.h"
#include "mesh_utils.h"
#include <unordered_map>


// idea generale: in questa classe sono contenute informazioni aggiuntive utili per trovare possibili
// intersezioni tra elementi della mesh. In particolare in questa classe non sono contenuti elementi (triangoli/tetraedri)
// ma ogni elemento è rappresentato dal suo bounding box (quadrati nel caso N=2, parrallelepipedi nel caso N=3)

namespace fdapde{
namespace core{


// M: local dimension, N: embedding dimension
template<int M, int N>
class StructuredGridSearch{
	private:
		using VectorType = SVector<N>;
		using BoundingBoxType = std::pair<VectorType, VectorType>;
		VectorType global_NE;
		VectorType global_SW;
		// questa variabile rappresenta dx, dy (e dz)
		std::array<double, N> cell_size = {0.0};
		// ogni elemento deve avere il suo bounding box e quindi un indice che rappresenta la posizione 
		// del bounding box all'interno dela mesh.
		// unsigned è l'indice di una cella. set<unsigned> è l'insieme degli Id degli elementi associati
		// a quella cella.
		std::unordered_map<unsigned, std::set<unsigned>> idx_map;
		// il numero delle celle nelle N direzioni
		std::array<unsigned, N> n_cells;
		// questa mappa contiene per ogni elemento (riconosciuto tramite id) il suo bounding box
		// il bounding box è definito come una std::pair dei punti che lo definiscono
		std::vector< BoundingBoxType > boxes_map;
		// serve per capire quando è necessario ricalcolare i bounding_box e gli indici degli elementi
		bool to_refresh_ = false;

		// private methods
		
	public:

		// constructor
		StructuredGridSearch(const Mesh<M, N> & mesh);

		// calcola l'indice della cella in cui cade 
		bool check_elem_size(const Element<M, N> & elem) const;
		
		unsigned compute_index(const SVector<N> & x) const; 
		// DOMANDA: HA SENSO SOLO PER N=3?
		// Main method: given the element el, find all elements that might intersect el
		std::unordered_set<unsigned> get_neighbouring_elements(const Element<M, N> & el) const; 
		// vengono calcolate le coordinate, cioè la posizione della cella che contiene il punto
		std::array<unsigned, N> compute_coordinates(const SVector<N> & x) const;
		// aggiorna la dimensione delle celle in base al numero di celle
		void update_cell_size();
		// elimina gli elementi dati in input
		void erase_elements(const std::unordered_set<unsigned> & el_ids);
		// add input elements 
		void add_elements(const std::vector<Element<M, N>> & elems);
		// aggiorna le informazioni per gli elementi dati in input
		// la f sinifica che non viene fatto un controllo per il refresh della struttura
		void update_f(const std::vector<Element<M, N>> & elements);
		// erase some elements, update informations of other elements and check if the structure
		// has to be refreshed
		void update(const std::vector<Element<M, N>> & elems_to_update, 
			const std::unordered_set<unsigned> & elems_to_erase, bool refresh_check = false);
		// refresh the structure: compute again cell_size, n_cells, and the index for all elements
		void refresh(const std::vector<Element<M, N>> & elems, const std::set<unsigned> & active_elems, unsigned n_nodes); 


		// getters
		BoundingBoxType get_global_bounding_box() const {return std::make_pair(global_SW, global_NE);}
		const std::array<double, N> & get_cell_size() const {return cell_size;}
		const std::array<unsigned, N> & get_n_cells() const {return n_cells;}
		std::pair<SVector<N>, SVector<N>> get_bounding_box(unsigned el_id) const;
		bool to_refresh() const {return to_refresh_;}



};

// ===============
// implementazione
// ===============


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
				for(unsigned i = 0; i < N; ++i)
				{
					double coord_v1 =  mesh.node( id_v1 )[i];
					double coord_v2 =  mesh.node( id_v2 )[i];
					double diff = std::abs( coord_v1 - coord_v2 );
					if(diff > cell_size[i])
						cell_size[i] = diff;
				}
				
			}
		}
		
	}
	for(unsigned i = 0; i < N; ++i){
		if(( global_NE[i] - global_SW[i] )<cell_size[i])
			std::cout<<"CELL_SIZE > GLOBAL_NE - GLOBAL_sw ???\n";
		n_cells[i] = ( global_NE[i] - global_SW[i] )/cell_size[i];
	}
	// cell_size is updated using the computed number of cells
	std::cout<<"cell_size: "<<cell_size[0]<<", "<<cell_size[1]<<", "<<cell_size[2]<<"\n";
	update_cell_size();

	// compute index for all elements
	for(unsigned id_el = 0; id_el < n_elements; ++id_el) 
	{
		auto bounding_box = mesh.element(id_el).bounding_box();
		// boxes_map[id_el] = bounding_box;
		boxes_map.push_back(bounding_box);
		SVector<N> middle_point =  ( bounding_box.first + bounding_box.second )*0.5;
		idx_map[compute_index(middle_point)].insert(id_el);
	}
	std::cout<<"num_cells: "<<n_cells[0]<<", "<<n_cells[1]<<", "<<n_cells[2]<<"\n";
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
			auto bounding_box = boxes_map[el_id];
			// compute the middle point of the bb to get the old index of the element
			SVector<N> middle_point = 0.5*(bounding_box.first + bounding_box.second);
			unsigned idx = compute_index(middle_point); // old index of the element
			// the bounding box of the element is erased and then computed again
			// boxes_map.erase(el_id);
			// the element is erased
			// se è tutto giusto non c'e bisogno di un controllo
			assert( idx_map.at(idx).find(el_id) != idx_map.at(idx).end() );
			idx_map.at(idx).erase(el_id);
			// now the new bounding box is computed
			auto new_bb = element.bounding_box();
			SVector<N> new_middle_point = 0.5*(new_bb.first + new_bb.second);
			boxes_map[el_id] = new_bb;
			unsigned new_idx = compute_index(new_middle_point); // compute the new index
			idx_map[new_idx].insert(el_id); // insert it in the map

		}
		
	}
}

template<int M, int N>
typename StructuredGridSearch<M, N>::BoundingBoxType StructuredGridSearch<M, N>::get_bounding_box(unsigned el_id) const
{
	/*
	if(boxes_map.find(el_id) != boxes_map.end())
		return boxes_map.at(el_id);
	// se l'elemento non è presente viene ritornato il global bb giusto per far tornare qualcosa
	return std::make_pair(global_SW, global_NE);
	*/
	return boxes_map[el_id];
}

template<int M, int N>
void StructuredGridSearch<M, N>::erase_elements(const std::unordered_set<unsigned> & el_ids)
{
	for(unsigned el_id : el_ids)
		// if(boxes_map.find(el_id) != boxes_map.end())
		{
			auto bounding_box = boxes_map[el_id];
			SVector<N> middle_point = 0.5*( bounding_box.first + bounding_box.second );
			unsigned el_idx = compute_index(middle_point);
			// se tutto è giusto nell'indice el_idx dovrebbe trovarsi l'id dell'elemento corrente
			// da eliminare
			// assert(idx_map[el_idx].find(el_id) != idx_map[el_idx].end());
			idx_map[el_idx].erase(el_id);
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
	double toll = 1.3;
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
		auto bounding_box = elem.bounding_box();
		auto middle_point = 0.5*(bounding_box.first + bounding_box.second);
		unsigned idx = compute_index(middle_point); // compute the index of the element
		idx_map[idx].insert(elem.ID());
		boxes_map[elem.ID()] = bounding_box;
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
	/*
	for(unsigned j = 0; j < N; ++j)
	{
		global_SW[j] = elems[tmp_id].coords()[j];
		global_NE[j] = elems[tmp_id].coords()[j];
	}*/
	global_NE = elems[tmp_id].coords()[0];
	global_SW = elems[tmp_id].coords()[0];
	// global_NE, global_SW are computed
	std::unordered_set<unsigned> visited_nodes;
	for(unsigned id_elem : active_elems) // loop over active elements
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
	std::cout<<"NE: "<<global_NE<<"\n";
	std::cout<<"SW: "<<global_SW<<"\n";


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
	std::cout<<"cell_size: "<<cell_size[0]<<", "<<cell_size[1]<<", "<<cell_size[2]<<"\n";
	for(unsigned k = 0; k < N; ++k) { 
		if(( global_NE[k] - global_SW[k] )<cell_size[k])
			std::cout<<"CELL_SIZE > GLOBAL_NE - GLOBAL_sw ???\n";
		if(n_cells[k]==0) {n_cells[k]=1;}
		n_cells[k] = ( global_NE[k] - global_SW[k] )/cell_size[k]; }
	// cell_size is updated using the computed number of cells
	update_cell_size();
	std::cout<<"num_cells: "<<n_cells[0]<<", "<<n_cells[1]<<", "<<n_cells[2]<<"\n";
	// idx_map has to be rebuilt from scratch
	idx_map.clear();

	// compute index for all elements
	boxes_map.clear();
	for(unsigned id_el : active_elems) 
	{
		auto bounding_box = elems[id_el].bounding_box();
		// boxes_map[id_el] = bounding_box;
		boxes_map.push_back(bounding_box);
		SVector<N> middle_point =  ( bounding_box.first + bounding_box.second )*0.5;
		idx_map[compute_index(middle_point)].insert(id_el);
	}
}

// per ora implemento con N=3
template<int M, int N>
std::unordered_set<unsigned> StructuredGridSearch<M, N>::get_neighbouring_elements(const Element<M, N> & el) const
{
	if constexpr(M == 2 && N == 2) {
       return {};
	} else {
	// compute the middle point of the element's bounding box
	auto bounding_box = el.bounding_box();
	// PROBLEMA: QUALI TRA FIRST E SECOND SONO NE E SW??
	// compute the coordinates of the SW and NE points of the bb
	std::array<unsigned, N> coordinates_NE = compute_coordinates(bounding_box.second);
	std::array<unsigned, N> coordinates_SW = compute_coordinates(bounding_box.first);

		// Find intersecting boxes
		//
		// We scan all the cells intersecting the reference box
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
					auto it = idx_map.find(idx);
					if(it != idx_map.end())
					{
						std::set<unsigned> range = it->second;
						for(unsigned id : range)
							if(boxes_intersection(boxes_map[id], bounding_box))
								res.insert(id);
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