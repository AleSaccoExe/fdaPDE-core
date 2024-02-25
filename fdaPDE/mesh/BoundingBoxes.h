#ifndef __BOUNDING_BOXES_H__
#define __BOUNDING_BOXES_H__

#include <cmath>
#include <unordered_set>

#include "../utils/symbols.h"
#include "../utils/intersection.h"
#include "mesh_utils.h"
#include "mesh.h"



// idea generale: in questa classe sono contenute informazioni aggiuntive utili per trovare possibili
// intersezioni tra elementi della mesh. In particolare in questa classe non sono contenuti elementi (triangoli/tetraedri)
// ma ogni elemento è rappresentato dal suo bounding box (quadrati nel caso N=2, parrallelepipedi nel caso N=3)

namespace fdapde{
namespace core{


// M: local dimension, N: embedding dimension
template<int M, int N>
class BoundingBoxes{
	private:
		SVector<N> global_NE;
		SVector<N> global_SW;
		// questa variabile rappresenta dx, dy (e dz)
		std::array<double, N> deltas = {0.0};
		// ogni elemento deve avere il suo bounding box e quindi un indice che rappresenta la posizione 
		// del bounding box all'interno dela mesh.
		// unsigned è l'indice di una cella. set<unsigned> è l'insieme degli Id degli elementi associati
		// a quella cella.
		std::unordered_map<unsigned, std::set<unsigned>> idx_map;
		// il numero delle celle nelle N direzioni
		std::array<unsigned, N> n_cells;
		// questa mappa contiene per ogni elemento (riconosciuto tramite id) il suo bounding box
		// il bounding box è definito come una std::pair dei punti che lo definiscono
		std::unordered_map< unsigned, std::pair<SVector<N>, SVector<N>> > boxes;
	public:

		// costruttori
		BoundingBoxes(const Mesh<M, N> & mesh);


		// calcola l'indice della cella in cui cade x
		unsigned compute_index(const SVector<N> & x) const; 
		// trova i possibili elementi che potrebbero intersecare l'elemento el
		// DOMANDA: HA SENSO SOLO PER N=3? 
		std::unordered_set<unsigned> get_neighbouring_elements(const Element<M, N> & el) const; 
		// vengono calcolate le coordinate, cioè la posizione della cella che contiene il punto
		std::array<unsigned, N> compute_coordinates(const SVector<N> & x) const;

		// getters
		std::pair<SVector<N>, SVector<N>> get_global_bounding_box() const {return std::make_pair(global_SW, global_NE);}
		const std::array<double, N> & get_deltas() const {return deltas;}
		const std::array<unsigned, N> & get_n_cells() const {return n_cells;}



};

// ===============
// implementazione
// ===============


template<int M, int N>
std::array<unsigned, N> BoundingBoxes<M, N>::compute_coordinates(const SVector<N> & x) const
{
	std::array<unsigned, N> coordinates;
	for(unsigned j = 0; j < N; ++j)
		coordinates[j] =static_cast<unsigned>( floor( (x[j] - global_SW[j] )/deltas[j] ) );
	return coordinates;
}

// costruttore
template<int M, int N>
BoundingBoxes<M, N>::BoundingBoxes(const Mesh<M, N> & mesh)
{
	
	// da fare: vedere tutti gli elementi della mesh in modo da calcolare dx, dy (e dz)
	// oltre a global_NE, global_SW.
	// poi una volta calcolati per ogni elemento bisogna calcolare il suo indice.
	unsigned n_nodes = mesh.n_nodes();
	unsigned n_elements = mesh.n_elements();

	// adesso calcolo global_NE e global_SW
	// global_NE = mesh.node(0).row(0);
	// global_SW = mesh.node(0).row(0);
	global_SW.setZero();
	global_NE.setZero();
	

	for(unsigned i = 0; i < n_nodes; ++i)
	{
		for(unsigned j = 0; j < N; ++j)
		{
			if(mesh.node(i)[j] < global_SW[j])
				global_SW[j] = mesh.node(i)[j];
			if(mesh.node(i)[j] > global_NE[j])
				global_NE[j] = mesh.node(i)[j];
		}

	}

	// ora calcolo dx, dy (e dz)

	for(unsigned id_el = 0; id_el < n_elements; ++id_el)
	{
		for(unsigned v1 = 0; v1 < ct_nvertices(M); ++v1)
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
					if(diff > deltas[i])
						deltas[i] = diff;
				}
				
			}
		}
		
	}

	for(unsigned i = 0; i < N; ++i)
		n_cells[i] = ( global_NE[i] - global_SW[i] )/deltas[i];


	// ora calcolo l'indice di ogni elemento

	for(unsigned id_el = 0; id_el < n_elements; ++id_el)
	{
		auto bounding_box = mesh.element(id_el).bounding_box();
		boxes[id_el] = bounding_box;
		SVector<N> middle_point =  ( bounding_box.first + bounding_box.second )*0.5;
		idx_map[compute_index(middle_point)].insert(id_el);
	}
	

}

template<int M, int N>
unsigned BoundingBoxes<M, N>::compute_index(const SVector<N> & x) const
{
	unsigned idx = 0;
	for(unsigned j = 0; j < N; ++j)
	{
		unsigned m = 1;
		for(unsigned i = 0; i < j; ++i)
			m*=n_cells[i];
		idx += m*static_cast<unsigned>( floor( (x[j] - global_SW[j] )/deltas[j] ) );
	}

	return idx;
}

// per ora implemento con N=3
template<int M, int N>
std::unordered_set<unsigned> BoundingBoxes<M, N>::get_neighbouring_elements(const Element<M, N> & el) const
{
	
	// viene calcolato il punto medio del bounding box dell'elemento
	auto bounding_box = el.bounding_box();
	// vengono calcolate le coordinate dei due punti che definiscono il bounding box
	// PROBLEMA: QUALI TRA FIRST E SECOND SONO NE E SW??
	std::array<unsigned, N> coordinates_NE = compute_coordinates(bounding_box.second);
	std::array<unsigned, N> coordinates_SW = compute_coordinates(bounding_box.first);

	//
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
							if(boxes_intersection(boxes.find(id)->second, bounding_box))
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

#endif // __BOUNDING_BOXES_H__