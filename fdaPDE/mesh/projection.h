#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include <unordered_set>
#include <map>
#include <vector>
#include <chrono>
#include <iostream>

#include "element.h"

namespace fdapde{
namespace core{
// returns a map element -> data
// This implementation is valid for meshes that do not require projetion of data points (2D, 3D)
template<int M, int N>
std::vector<std::unordered_set<unsigned>> project(const std::vector<Element<M, N>>&  elems, DMatrix<double> & data, 
										const std::unordered_set<unsigned>& data_ids)
{
	// Initialize the connections:
	std::vector<std::unordered_set<unsigned>> new_elem_to_data(elems.size());
	for(unsigned datum_id : data_ids) // Loop over all the data points
	{
		bool ok = false;
		const SVector<N> datum = data.row(datum_id); // take the position of the current datum
		for(unsigned i = 0; i < elems.size(); ++i) // loop over all the input elements 
		{
			// if the datum belongs to the current element, update the connections:
			if(elems[i].contains(datum)) { new_elem_to_data[i].insert(datum_id); ok = true;} 
		}
		if(!ok)
		{
			std::cout<<"non va bene! Il dato non appartiene afìd alcun elemento\n";
			std::cout<<"elementi modificati:\n";
			for(const auto & elem : elems)
				std::cout<<elem.coords()[0]<<"\n"<<elem.coords()[1]<<"\n"<<elem.coords()[2]<<std::endl;
			std::cout<<"posizione del dato:\n";
			std::cout<<datum<<std::endl;
		}
	}
	return new_elem_to_data;
}


// returns a map element -> data
// valid implementation for surfaces (2.5D)
template<>
std::vector<std::unordered_set<unsigned>> project(const std::vector<Element<2, 3>>&  elems, DMatrix<double> & data, 
										const std::unordered_set<unsigned>& data_ids)
{ 
	std::vector<std::unordered_set<unsigned>> new_elem_to_data(elems.size());
	for(unsigned datum_id : data_ids) // loop over data points to project
	{
		unsigned opt_elem; // id of the element on which the datum will be projected
		SVector<3> datum = data.row(datum_id); // datum position
		SVector<3> opt_pos = datum; // optimal position of the projected point (minimizing the distance)
		double opt_dist = std::numeric_limits<double>::max(); // distance of projection
		bool done = false;
		bool inside = false;
		bool already_inside = false;
		for(unsigned i = 0; i < elems.size(); ++i) // loop over elements
		{
			const auto & elem = elems[i];
			// si vede se il dato coincide con un vertice
			/*for(SVector<3> vertex : elem.coords())
				if((vertex-datum).norm()<DOUBLE_TOLERANCE) // il dato coincide col vertice
				{
					opt_pos = datum;
					opt_dist = 0.0;
					new_elem_to_data[i].insert(datum_id);
					done = true;
					already_inside = true;
				}
			// si vede se il dato giace su uno spigolo del triangolo
			for(unsigned v1 = 0; v1 < 3; ++v1)
				for(unsigned v2 = v1+1; v2 < 3; ++v2)
				{
					SVector<3> A = elem.coords()[v1];
					SVector<3> B = elem.coords()[v2];
					if( std::abs( (A-datum).norm()+(B-datum).norm()-(A-B).norm() ) < DOUBLE_TOLERANCE ) // il punto è dentro il lato AB
					{
						opt_pos = datum;
						opt_dist = 0.0;
						new_elem_to_data[i].insert(datum_id);
						done = true;
						already_inside = true;
					}
				}*/
			if(elem.contains(datum)) // if the point is already inside an element ...
			{ 
				opt_pos = datum;
				opt_dist = 0.0;
				new_elem_to_data[i].insert(datum_id); // update the connections
				done = true; // will be useful later
				already_inside = true; // will be useful later
			}
			// if the datum is not already inside the current element, 
			// try to project the point onto it:
			if(!already_inside) 
			{
				SVector<3> p_datum = elem.hyperplane().project(datum); // project the datum on the element's plane
				double dist = (datum - p_datum).norm(); // distance of projection
				if(elem.contains(p_datum) && dist<opt_dist) // if the projection falls inside the element ...
				{
					opt_dist = dist;
					opt_pos = p_datum;
					opt_elem = i;
					done = true;
					inside = true; // will be useful later
				}
			}
		} // end loop on elements
		// if the point was projected inside an element, but it wasn't already inside of it,
		// update the connections:
		if(inside && !already_inside) {new_elem_to_data[opt_elem].insert(datum_id);}
		if(!done) // this is true if the datum couldn't be projected inside of an element
		{
			// Now the datum has to be projected onto the elements' edges.
			std::set<unsigned> elems_on_datum; // will contain the elements the datum will fall onto during the projection
			for(unsigned i = 0; i < elems.size(); ++i) // loop over the elements
			{
				const auto & elem = elems[i];
				SVector<3> A = elem.coords()[0];
				SVector<3> B = elem.coords()[1];
				SVector<3> C = elem.coords()[2];
				// projection onto edge AB
				double t = -(A-datum).dot(B-A)/(B-A).squaredNorm();
				SVector<3> p_datum = A + t*(B-A); // projections of the current datum
				if(0 <= t && t <= 1) // projection inside of the edge AB
				{
					double dist = (p_datum-datum).norm();
					if((p_datum-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // projection coincides with vertex A
				{
					double dist = (datum - A).norm();
					if((A-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = A;
						opt_dist = dist;
					}
				}
				else if(t > 1) // projection coincides with vertex B
				{
					double dist = (datum - B).norm();
					if((B-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = B;
						opt_dist = dist;
					}
				}
				// projection onto edge BC
				t = -(B-datum).dot(C-B)/(C-B).squaredNorm();
				p_datum = B + t*(C-B);
				if(0 <= t && t <= 1) // projection inside of the edge BC
				{
					double dist = (p_datum-datum).norm();
					if((p_datum-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // projection coincides with vertex B
				{
					double dist = (datum - B).norm();
					if((B-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = B;
						opt_dist = dist;
					}
				}
				else if(t > 1) // projection coincides with vertex C
				{
					double dist = (datum - C).norm();
					if((C-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = C;
						opt_dist = dist;
					}
				}
				// projection on edge CA
				t = -(C-datum).dot(A-C)/(A-C).squaredNorm();
				p_datum = C + t*(A-C);
				if(0 <= t && t <= 1) // projection inside of the edge CA
				{
					double dist = (p_datum-datum).norm();
					if((p_datum-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // projection coincides with vertex C
				{
					double dist = (datum - C).norm();
					if((C-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = C;
						opt_dist = dist;
					}
				}
				else if(t > 1) // projection coincides with vertex A
				{
					double dist = (datum - A).norm();
					if((A-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = A;
						opt_dist = dist;
					}
				}
 			} // end loop on elements
 			for(unsigned elem_idx : elems_on_datum) {new_elem_to_data[elem_idx].insert(datum_id);}
		} // if(!done)
		data.row(datum_id) = opt_pos; // update the position of the datum
	} // fine loop sui dati
	return new_elem_to_data;
} // project

// this function is used by collapsing cost DataDispCost
// This implementation is valid for meshes that do not require projetion of data points (2D, 3D)
template<int M, int N>
std::vector<std::unordered_set<unsigned>> projection_info(const std::vector<Element<M, N>>&  elems, const DMatrix<double> & data, 
										const std::unordered_set<unsigned>& data_ids)
{
	std::vector<std::unordered_set<unsigned>> new_elem_to_data(elems.size());
	for(unsigned datum_id : data_ids)
	{
		const SVector<N>& datum = data.row(datum_id);
		for(unsigned i = 0; i < elems.size(); ++i)
		{
			if(elems[i].contains(datum)) { new_elem_to_data[i].insert(datum_id); } 
		}
	}
	return new_elem_to_data;
}										

template<>
std::vector<std::unordered_set<unsigned>> projection_info(const std::vector<Element<2, 3>>&  elems, const DMatrix<double> & data, 
										const std::unordered_set<unsigned>& data_ids)
{
	std::vector<std::unordered_set<unsigned>> new_elem_to_data(elems.size());
	for(unsigned datum_id : data_ids) // loop sui dati da proiettare
	{
		unsigned opt_elem;
		SVector<3> datum = data.row(datum_id);
		SVector<3> opt_pos = datum; // posizione ottima del punto
		double opt_dist = std::numeric_limits<double>::max(); // distanza della proiezione
		bool done = false;
		bool already_inside = false;
		bool inside = false;
		for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
		{
			const auto & elem = elems[i];
			// si vede se il dato coincide con un vertice
			/*for(SVector<3> vertex : elem.coords())
				if((vertex-datum).norm()<DOUBLE_TOLERANCE) // il dato coincide col vertice
				{
					opt_pos = datum;
					opt_dist = 0.0;
					new_elem_to_data[i].insert(datum_id);
					done = true;
					already_inside = true;
				}
			// si vede se il dato giace su uno spigolo del triangolo
			for(unsigned v1 = 0; v1 < 3; ++v1)
				for(unsigned v2 = v1+1; v2 < 3; ++v2)
				{
					SVector<3> A = elem.coords()[v1];
					SVector<3> B = elem.coords()[v2];
					if( std::abs( (A-datum).norm()+(B-datum).norm()-(A-B).norm() ) < DOUBLE_TOLERANCE ) // il punto è dentro il lato AB
					{
						opt_pos = datum;
						opt_dist = 0.0;
						new_elem_to_data[i].insert(datum_id);
						done = true;
						already_inside = true;
					}
				}*/
			if(elem.contains(datum)) // si verifica se il punto è interno al triangolo
			{
				opt_pos = datum;
				opt_dist = 0.0;
				new_elem_to_data[i].insert(datum_id);
				done = true;
				already_inside = true;
			}
			if(!already_inside) // il dato deve essere proiettato
			{
				SVector<3> p_datum = elem.hyperplane().project(datum); // viene proiettato il dato
				double dist = (datum - p_datum).norm();
				if(elem.contains(p_datum) && dist<opt_dist)
				{
					opt_dist = dist;
					opt_pos = p_datum;
					done = true;
					opt_elem = i;
					inside = true;
				}
			}
		} // fine loop sugli elementi
		if(inside && !already_inside) {new_elem_to_data[opt_elem].insert(datum_id);}
		if(!done) // il dato non poteva essere proiettato all'interno di un elemento
		{
			std::unordered_set<unsigned> elems_on_datum; // contiene gli elementi a cui il dato appartiene
			for(unsigned i = 0; i < elems.size(); ++i)
			{
				const auto & elem = elems[i];
				const SVector<3>& A = elem.coords()[0];
				const SVector<3>& B = elem.coords()[1];
				const SVector<3>& C = elem.coords()[2];
				// proiezione sul lato AB
				double t = -(A-datum).dot(B-A)/(B-A).squaredNorm();
				SVector<3> p_datum = A + t*(B-A);
				if(0 <= t && t <= 1) // il punto cade nel lato AB
				{
					double dist = (p_datum-datum).norm();
					if((p_datum-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // la proiezione coincide con il punto A
				{
					double dist = (datum - A).norm();
					if((A-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = A;
						opt_dist = dist;
					}
				}
				else if(t > 1) // la proiezione coincide con il punto B
				{
					double dist = (datum - B).norm();
					if((B-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = B;
						opt_dist = dist;
					}
				}
				// proiezione sul lato BC
				t = -(B-datum).dot(C-B)/(C-B).squaredNorm();
				p_datum = B + t*(C-B);
				if(0 <= t && t <= 1) // il punto cade nel lato AB
				{
					double dist = (p_datum-datum).norm();
					if((p_datum-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // la proiezione coincide con il punto B
				{
					double dist = (datum - B).norm();
					if((B-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = B;
						opt_dist = dist;
					}
				}
				else if(t > 1) // la proiezione coincide con il punto C
				{
					double dist = (datum - C).norm();
					if((C-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = C;
						opt_dist = dist;
					}
				}
				// proiezione sul lato CA
				t = -(C-datum).dot(A-C)/(A-C).squaredNorm();
				p_datum = C + t*(A-C);
				if(0 <= t && t <= 1) // il punto cade nel lato AB
				{
					double dist = (p_datum-datum).norm();
					if((p_datum-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // la proiezione coincide con il punto C
				{
					double dist = (datum - C).norm();
					if((C-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = C;
						opt_dist = dist;
					}
				}
				else if(t > 1) // la proiezione coincide con il punto A
				{
					double dist = (datum - A).norm();
					if((A-opt_pos).norm()<DOUBLE_TOLERANCE) {elems_on_datum.insert(i);}
					else if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						elems_on_datum.clear();
						elems_on_datum.insert(i);
						opt_pos = A;
						opt_dist = dist;
					}
				}
 			} // fine loop su elementi
 			for(unsigned elem_idx : elems_on_datum) {new_elem_to_data[elem_idx].insert(datum_id);}
		}
	} // fine loop sui dati
	return new_elem_to_data;
} // project_info


// function used by collapsing cost DataDistCost
// returns the projection of the input datum onto the triangles contained in "elems"
SVector<3> project(const std::vector<Element<2, 3>>&  elems, SVector<3> datum)
{
		SVector<3> opt_pos = datum; // posizione ottima del punto
		double opt_dist = std::numeric_limits<double>::max(); // distanza della proiezione
		bool done = false;
		bool already_inside = false;
		for(unsigned i = 0; i < elems.size(); ++i) // loop sugli elementi
		{
			const auto & elem = elems[i];
			// si vede se il dato coincide con un vertice
			/*for(SVector<3> vertex : elem.coords())
				if((vertex-datum).norm()<DOUBLE_TOLERANCE) // il dato coincide col vertice
				{
					opt_pos = datum;
					opt_dist = 0.0;
					done = true;
					already_inside = true;
				}
			// si vede se il dato giace su uno spigolo del triangolo
			for(unsigned v1 = 0; v1 < 3; ++v1)
				for(unsigned v2 = v1+1; v2 < 3; ++v2)
				{
					SVector<3> A = elem.coords()[v1];
					SVector<3> B = elem.coords()[v2];
					if( std::abs( (A-datum).norm()+(B-datum).norm()-(A-B).norm() ) < DOUBLE_TOLERANCE ) // il punto è dentro il lato AB
					{
						opt_pos = datum;
						opt_dist = 0.0;
						done = true;
						already_inside = true;
					}
				}*/
			if(elem.contains(datum)) // si verifica se il punto è interno al triangolo
			{
				opt_pos = datum;
				opt_dist = 0.0;
				done = true;
				already_inside = true;
			}
			if(!already_inside) // il dato deve essere proiettato
			{
				unsigned opt_elem;
				SVector<3> p_datum = elem.hyperplane().project(datum); // viene proiettato il dato
				double dist = (datum - p_datum).norm();
				if(elem.contains(p_datum) && dist<opt_dist)
				{
					opt_dist = dist;
					opt_pos = p_datum;
					done = true;
				}
			}
		} // fine loop sugli elementi
		if(!done) // il dato non poteva essere proiettato all'interno di un elemento
		{
			for(unsigned i = 0; i < elems.size(); ++i)
			{
				const auto & elem = elems[i];
				SVector<3> A = elem.coords()[0];
				SVector<3> B = elem.coords()[1];
				SVector<3> C = elem.coords()[2];
				// proiezione sul lato AB
				double t = -(A-datum).dot(B-A)/(B-A).squaredNorm();
				SVector<3> p_datum = A + t*(B-A);
				if(0 <= t && t <= 1) // il punto cade nel lato AB
				{
					double dist = (p_datum-datum).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // la proiezione coincide con il punto A
				{
					double dist = (datum - A).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = A;
						opt_dist = dist;
					}
				}
				else if(t > 1) // la proiezione coincide con il punto B
				{
					double dist = (datum - B).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = B;
						opt_dist = dist;
					}
				}
				// proiezione sul lato BC
				t = -(B-datum).dot(C-B)/(C-B).squaredNorm();
				p_datum = B + t*(C-B);
				if(0 <= t && t <= 1) // il punto cade nel lato AB
				{
					double dist = (p_datum-datum).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				if(t < 0) // la proiezione coincide con il punto B
				{
					double dist = (datum - B).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{						opt_pos = B;
						opt_dist = dist;
					}
				}
				else if(t > 1) // la proiezione coincide con il punto C
				{
					double dist = (datum - C).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = C;
						opt_dist = dist;
					}
				}
				// proiezione sul lato CA
				t = -(C-datum).dot(A-C)/(A-C).squaredNorm();
				p_datum = C + t*(A-C);
				if(0 <= t && t <= 1) // il punto cade nel lato AB
				{
					double dist = (p_datum-datum).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = p_datum;
						opt_dist = dist;
					}
				}
				else if(t < 0) // la proiezione coincide con il punto C
				{
					double dist = (datum - C).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = C;
						opt_dist = dist;
					}
				}
				else if(t > 1) // la proiezione coincide con il punto A
				{
					double dist = (datum - A).norm();
					if(dist < opt_dist-DOUBLE_TOLERANCE)
					{
						opt_pos = A;
						opt_dist = dist;
					}
				}
 			} // fine loop su elementi
		}
	return opt_pos;
} // project


} // core
} // fdapde


#endif