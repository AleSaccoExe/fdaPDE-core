#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include <unordered_set>
#include <vector>
#include "element.h"

namespace fdapde{
namespace core{


template<typename Element_>
std::vector<std::pair<int, int>> project(const std::vector<Element_> & elems, 
 	DMatrix<double> & data, const std::unordered_set<unsigned> & ids)
{
	constexpr int N = Element_::embedding_dimension;
	std::vector<std::pair<int, int>> new_positions;
	for(unsigned id : ids) // loop sui dati
	{
		double opt_dist(std::numeric_limits<double>::max());
		bool inside_elem = false;
		unsigned opt_id; // id della posizione ottima del dato
		SVector<N> datum = data.row(id);
		SVector<N> opt_pos; // la posizione ottimale del dato
		for(const auto & elem : elems) // loop sugli elementi
		{
			SVector<N> p_datum = elem.hyperplane().project(datum); // viene proiettato il dato
			if(elem.contains(p_datum)) // se l'elemento contiene il dato proiettato sul piano del triangolo
			{
				double dist = (datum - p_datum).norm(); // distanza 
				if(dist < opt_dist)
				{
					opt_dist = dist;
					opt_pos = p_datum;
					inside_elem = true; // il dato cade in un elemento
					opt_id = elem.ID();
				}
			}
		} // loop sugli elementi
		if(inside_elem)
		{
			data.row(id) = opt_pos;
			new_positions.push_back({-2, opt_id}); // -2 indica che è nell'elemento opt_id
		}
		else // il dato non può essere proiettato dentro gli elementi
		{
			unsigned node_id1, node_id2;
			opt_dist = std::numeric_limits<double>::max();
			for(const auto & elem : elems) // loop sugli elementi
			{
				// vengono presi i nodi dell'elemento
				auto A = elem.coords()[0];
				auto B = elem.coords()[1];
				auto C = elem.coords()[2];
				// lato AB
				double t = -( datum - A ).dot( A - B )/( B - A ).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = A + t*(B - A);
					double dist = (p_datum - datum).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = elem.node_ids()[0];
						node_id2 = elem.node_ids()[1];
						opt_pos = p_datum;
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - A).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = -1; // vuol dire che la posizione ottima è sul vertice
						node_id2 = elem.node_ids()[0]; // id del nodo A
						opt_pos = A;
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - B).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = -1; // vuol dire che la posizione ottima è sul vertice
						node_id2 = elem.node_ids()[1]; // id del nodo B
						opt_pos = B;
					}
				}
				// lato BC
				t = -(datum - B).dot(B - C)/(C- B).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = B + t*(C - B);
					double dist = (p_datum - datum).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = elem.node_ids()[1];
						node_id2 = elem.node_ids()[2];
						opt_pos = p_datum;
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - B).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = -1; // vuol dire che la posizione ottima è sul vertice
						node_id2 = elem.node_ids()[1]; // id del nodo B
						opt_pos = B;
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - C).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = -1; // vuol dire che la posizione ottima è sul vertice
						node_id2 = elem.node_ids()[2]; // id del nodo C
						opt_pos = C;
					}
				}
				// lato CA
				t = -(datum - C).dot(C - A)/(A - C).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = C + t*(A - C);
					double dist = (p_datum - datum).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = elem.node_ids()[2];
						node_id2 = elem.node_ids()[0];
						opt_pos = p_datum;
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - C).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = -1; // vuol dire che la posizione ottima è sul vertice
						node_id2 = elem.node_ids()[0]; // id del nodo C
						opt_pos = C;
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - A).norm();
					if(opt_dist > dist)
					{
						opt_dist = dist;
						node_id1 = -1; // vuol dire che la posizione ottima è sul vertice
						node_id2 = elem.node_ids()[0]; // id del nodo A
						opt_pos = A;
					}
				}
				
				
			} // loop sugli elementi
			// aggiornata la posizione del dato
			data.row(id) = opt_pos;
			new_positions.push_back({node_id1, node_id2});
		} // else 

	} // loop sui dati
	return new_positions;
}


} // core
} // fdapde

#endif // __PROJECTION_H__