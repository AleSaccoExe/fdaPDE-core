#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include <unordered_set>
#include <vector>
#include "element.h"

namespace fdapde{
namespace core{

// proiezione dei dati sugli elementi forniti in input.
	/*
std::vector<std::pair<int, int>> project(const std::vector<Element<2, 3>> & elems, 
 	DMatrix<double> & data, const std::unordered_set<unsigned> & ids)
{
	constexpr int N = 3;
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
				inside_elem = true; // il dato cade in un elemento
				double dist = (datum - p_datum).norm(); // distanza
				else if(dist < opt_dist)
				{
					opt_dist = dist;
					opt_pos = p_datum;
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
*/

std::vector<std::unordered_set<unsigned>> project(const std::vector<Element<2, 3>> & elems, 
 	DMatrix<double> & data, const std::unordered_set<unsigned> & data_ids)
{
	std::vector<std::unordered_set<unsigned>> new_elem_to_data(elems.size());
	constexpr int N = 3;

	for(unsigned datum_id : data_ids) // loop sui dati
	{
		SVector<N> opt_pos; // posizione del dato dopo la proiezione
		double opt_dist(std::numeric_limits<double>::max()); // distanza della proiezione
		bool inside_elem = false;
		SVector<N> datum = data.row(datum_id); // dato da proiettare
		// il set degli elementi a cui il dato appartiene
		std::set<unsigned> elems_on_datum;
		for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
		{
			const auto & elem = elems[i];
			SVector<N> p_datum = elem.hyperplane().project(datum); // viene proiettato il dato
			if(elem.contains(p_datum)) // se l'elemento contiene il dato proiettato sul piano del triangolo
			{
				inside_elem = true;
				double dist = (datum - p_datum).norm(); // distanza
				// il dato potrebbe essere sul bordo:
				if( dist < DOUBLE_TOLERANCE )
				{
					opt_dist = 0.0;
					elems_on_datum.insert(i);
					opt_pos = datum; // il dato non è veramente proiettato
				}
				else if(dist < opt_dist)
				{
					opt_dist = dist;
					elems_on_datum.clear();
					elems_on_datum.insert(i);
					opt_pos = p_datum; 
				}
			}
		} // loop sugli elementi
		if(inside_elem) {
			for(unsigned elem_idx : elems_on_datum) {new_elem_to_data[elem_idx].insert(datum_id);}
		}
		else // il dato non può essere proiettato dentro gli elementi
		{
			assert(elems_on_datum.size() == 0);
			opt_dist = std::numeric_limits<double>::max();
			for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
			{
				const auto & elem = elems[i];
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
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = p_datum;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - A).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = A;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - B).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = B;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					}
				}
				// lato BC
				t = -(datum - B).dot(B - C)/(C- B).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = B + t*(C - B);
					double dist = (p_datum - datum).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = p_datum;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - B).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = B;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - C).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = C;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					}
				}
				// lato CA
				t = -(datum - C).dot(C - A)/(A - C).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = C + t*(A - C);
					double dist = (p_datum - datum).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = p_datum;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - C).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = C;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - A).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						opt_pos = A;
						elems_on_datum.clear();
						elems_on_datum.insert(i);
					}
				}
			} // loop sugli elementi
			// aggiornata la posizione del dato
			for(unsigned i : elems_on_datum) {new_elem_to_data[i].insert(datum_id);}
		} // else
		data.row(datum_id) = opt_pos; // cambiata la posizione del dato
	} // loop sui dati
	return new_elem_to_data;
}


// Calcola la distanza della proiezione sugli elementi in input dalla posizione iniziale del dato
SVector<3> project(const std::vector<Element<2, 3>> & elems, SVector<3> datum)
{
	double opt_dist(std::numeric_limits<double>::max());
	bool inside_elem = false;
	SVector<3> opt_pos; // la posizione ottimale del dato
	for(const auto & elem : elems) // loop sugli elementi
	{
		SVector<3> p_datum = elem.hyperplane().project(datum); // viene proiettato il dato
		if(elem.contains(p_datum)) // se l'elemento contiene il dato proiettato sul piano del triangolo
		{
			double dist = (datum - p_datum).norm(); // distanza 
			if(dist < opt_dist)
			{
				opt_dist = dist;
				opt_pos = p_datum;
				inside_elem = true; // il dato cade in un elemento
			}
		}
	} // loop sugli elementi
	if(inside_elem) {return opt_pos;}
	else // il dato non può essere proiettato dentro gli elementi
	{
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
				SVector<3> p_datum = A + t*(B - A);
				double dist = (p_datum - datum).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = p_datum;}
			}
			else if(t <=0)
			{
				double dist = (datum - A).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = A;}
			}
			else if(t >= 1)
			{
				double dist = (datum - B).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = B;}
			}
			// lato BC
			t = -(datum - B).dot(B - C)/(C- B).squaredNorm();
			if(t > 0 && t < 1) // la proiezione cade dentro al lato
			{
				SVector<3> p_datum = B + t*(C - B);
				double dist = (p_datum - datum).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = p_datum;}
			}
			else if(t <=0)
			{
				double dist = (datum - B).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = B;}
			}
			else if(t >= 1)
			{
				double dist = (datum - C).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = C;}
			}
			// lato CA
			t = -(datum - C).dot(C - A)/(A - C).squaredNorm();
			if(t > 0 && t < 1) // la proiezione cade dentro al lato
			{
				SVector<3> p_datum = C + t*(A - C);
				double dist = (p_datum - datum).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = p_datum;}
			}
			else if(t <=0)
			{
				double dist = (datum - C).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = C;}
			}
			else if(t >= 1)
			{
				double dist = (datum - A).norm();
				if(opt_dist > dist) {
					opt_dist = dist;
					opt_pos = C;}
			}
			
			
		} // loop sugli elementi
	} // else 
	return opt_pos;
} // dist_projection
/*
// questa funzione torna per ogni elemento un set degli id dei dati che ora appartengono all'elemento
std::vector<std::set<unsigned>> projection_info(const std::vector<Element<2, 3>> & elems, 
 	const DMatrix<double> & data, const std::unordered_set<unsigned> & data_ids)
{
	std::vector<std::set<unsigned>> new_elem_to_data(elems.size());
	constexpr int N = 3;
	for(unsigned datum_id : data_ids) // loop sui dati
	{
		double opt_dist(std::numeric_limits<double>::max());
		bool inside_elem = false;
		SVector<N> datum = data.row(datum_id);
		unsigned opt_pos; // la posizione ottimale del dato
		for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
		{
			const auto & elem = elems[i];
			SVector<N> p_datum = elem.hyperplane().project(datum); // viene proiettato il dato
			if(elem.contains(p_datum)) // se l'elemento contiene il dato proiettato sul piano del triangolo
			{
				double dist = (datum - p_datum).norm(); // distanza 
				if(dist < opt_dist)
				{
					opt_dist = dist;
					inside_elem = true; // il dato cade in un elemento
					opt_pos = i; 
				}
			}
		} // loop sugli elementi
		if(inside_elem) {new_elem_to_data[opt_pos].insert(datum_id);}
		else // il dato non può essere proiettato dentro gli elementi
		{
			// contiene gli id degli elementi a cui il dato appartiene dopo la proiezione
			std::set<unsigned> best_pos; 
			opt_dist = std::numeric_limits<double>::max();
			for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
			{
				const auto & elem = elems[i];
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
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - A).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - B).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				// lato BC
				t = -(datum - B).dot(B - C)/(C- B).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = B + t*(C - B);
					double dist = (p_datum - datum).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - B).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - C).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				// lato CA
				t = -(datum - C).dot(C - A)/(A - C).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = C + t*(A - C);
					double dist = (p_datum - datum).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - C).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - A).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
			} // loop sugli elementi
			// aggiornata la posizione del dato
			for(unsigned i : best_pos) {new_elem_to_data[i].insert(datum_id);}
		} // else 
	} // loop sui dati
	return new_elem_to_data;
}
*/

// questa funzione torna per ogni elemento un set degli id dei dati che ora appartengono all'elemento
std::vector<std::set<unsigned>> projection_info(const std::vector<Element<2, 3>> & elems, 
 	const DMatrix<double> & data, const std::unordered_set<unsigned> & data_ids)
{
	std::vector<std::set<unsigned>> new_elem_to_data(elems.size());
	constexpr int N = 3;
	for(unsigned datum_id : data_ids) // loop sui dati
	{
		double opt_dist(std::numeric_limits<double>::max());
		bool inside_elem = false;
		SVector<N> datum = data.row(datum_id);
		// il set degli elementi a cui il dato appartiene
		std::set<unsigned> elems_on_datum;
		unsigned opt_pos; // la posizione ottimale del dato (indice dell'elemento)
		for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
		{
			const auto & elem = elems[i];
			SVector<N> p_datum = elem.hyperplane().project(datum); // viene proiettato il dato
			if(elem.contains(p_datum)) // se l'elemento contiene il dato proiettato sul piano del triangolo
			{
				inside_elem = true;
				double dist = (datum - p_datum).norm(); // distanza
				// il dato potrebbe essere sul bordo:
				if( dist < DOUBLE_TOLERANCE )
				{
					opt_dist = 0.0;
					elems_on_datum.insert(i);
				}
				else if(dist < opt_dist)
				{
					opt_dist = dist;
					elems_on_datum.clear();
					elems_on_datum.insert(i); 
				}
			}
		} // loop sugli elementi
		if(inside_elem) {
			for(unsigned elem_idx : elems_on_datum) {new_elem_to_data[elem_idx].insert(datum_id);}
		}
		else // il dato non può essere proiettato dentro gli elementi
		{
			// contiene gli id degli elementi a cui il dato appartiene dopo la proiezione
			std::set<unsigned> best_pos; 
			opt_dist = std::numeric_limits<double>::max();
			for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
			{
				const auto & elem = elems[i];
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
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - A).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - B).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				// lato BC
				t = -(datum - B).dot(B - C)/(C- B).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = B + t*(C - B);
					double dist = (p_datum - datum).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - B).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - C).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				// lato CA
				t = -(datum - C).dot(C - A)/(A - C).squaredNorm();
				if(t > 0 && t < 1) // la proiezione cade dentro al lato
				{
					SVector<N> p_datum = C + t*(A - C);
					double dist = (p_datum - datum).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					} 
				}
				else if(t <=0)
				{
					double dist = (datum - C).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
				else if(t >= 1)
				{
					double dist = (datum - A).norm();
					if(std::abs(opt_dist-dist)<DOUBLE_TOLERANCE) {best_pos.insert(i);}
					else if(opt_dist > dist)
					{
						opt_dist = dist;
						best_pos.clear();
						best_pos.insert(i);
					}
				}
			} // loop sugli elementi
			// aggiornata la posizione del dato
			for(unsigned i : best_pos) {new_elem_to_data[i].insert(datum_id);}
		} // else 
	} // loop sui dati
	return new_elem_to_data;
}

} // core
} // fdapde

#endif // __PROJECTION_H__