#ifndef __SHARP_ELEMS_COST_H__
#define __SHARP_ELEMS_COST_H__

#include "../../mesh/element.h"
#include "../symbols.h"
#include "CostObjBase.h"
#include <cmath>
#include <iostream>

namespace fdapde{
namespace core{

template<int N, int K>
double circumradius(const std::array<SVector<N>, K> & points)
{
	return 0.0;
}
// specializzazione per triangoli
template<int N>
double circumradius(const std::array<SVector<N>, 3> & points)
{
	const SVector<N> A = points[0];
	const SVector<N> B = points[1];
	const SVector<N> C = points[2];

	SVector<N> a = A-C;
	SVector<N> b = B-C;
	return a.norm()*b.norm()*(a-b).norm()/(2.0*(a.cross(b)).norm());

}

// compute the max cos in a triangle identified by the coordinates of 
// its vertices
template<int N>
double get_max_cos(const SVector<N>& A, const SVector<N>& B, const SVector<N>& C) 
{
	double max = 0.0;
	double cos1 = std::abs((B-A).dot(C-A))/( (B-A).norm()*(C-A).norm() );
	if(max < cos1) {max = cos1;}
	double cos2 = std::abs((C-B).dot(A-B))/( (C-B).norm()*(A-B).norm() );
	if(max < cos2) {max = cos2;}
	double cos3 = std::abs((A-C).dot(B-C))/( (A-C).norm()*(B-C).norm() );
	if(max < cos3) {max = cos3;}
	return max;

}


/*template<int M, int N>
struct SharpElemsCost{
	double operator()(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete, 
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v, 
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		// return get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids)/max_;
		return get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids);
	}

	double get_cost(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete, 
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v, 
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		// all triangles in the element
		double max_cos1 = 0.0;
		// loop on elements to delete
		for(const auto & elem : elems_to_delete)
		{
			const auto& coords = elem.coords();
			for(unsigned i = 0; i < ct_binomial_coefficient(M+1, 3); ++i){
				double cos = get_max_cos(coords[triangles(i, 0)], coords[triangles(i, 1)], coords[triangles(i, 2)]);
			// double cos = get_max_cos(elem.coords());
				if(cos > max_cos1) {max_cos1 = cos;}
			}
		}
		// loop on elements to modify
		for(const auto & elem : elems_to_modify)
		{
			const auto& coords = elem.coords();
			for(unsigned i = 0; i < ct_binomial_coefficient(M+1, 3); ++i){
				double cos = get_max_cos(coords[triangles(i, 0)], coords[triangles(i, 1)], coords[triangles(i, 2)]);
			// double cos = get_max_cos(elem.coords());
				if(cos > max_cos1) {max_cos1 = cos;}
			}
		}
		double max_cos2 = 0.0;
		// loop on elements modified
		for(const auto & elem : elems_modified)
		{
			const auto& coords = elem.coords();
			for(unsigned i = 0; i < ct_binomial_coefficient(M+1, 3); ++i){
				double cos = get_max_cos(coords[triangles(i, 0)], coords[triangles(i, 1)], coords[triangles(i, 2)]);
			// double cos = get_max_cos(elem.coords());
				if(cos > max_cos2) {max_cos2 = cos;}
			}
		}
		double csi = (1.-max_cos1)*(1.-max_cos1)/( (1.-max_cos2)*(1.-max_cos2) );
		return 1./std::tanh(csi) - 1./csi;
	}
	
	void setup(Simplification<M, N> * p_simp){}

	void update_min(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete, 
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v, 
					  const std::unordered_set<unsigned> & data_ids )
	{
		double cost = get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids);
		if(cost < min_) {min_ = cost;}
	}

	void update_max()
	{
		if(max_ < min_) {max_ = min_;}
		min_  = std::numeric_limits<double>::max();
	}

	void update(const std::vector<Element<M, N>> & elems_to_delete,
				const std::vector<Element<M, N>> & elems_modified) {}
	void set_threshold(double new_threshold) {}
	bool check_update() {return false;}
	double max_ = 0.0;
	double min_ = std::numeric_limits<double>::max();

	SMatrix<ct_binomial_coefficient(M+1, 3), 3, int> triangles = combinations<3, M+1>();
};*/

template<int M, int N>
class SharpElemsCost: public CostObjBase<M, N, SharpElemsCost<M, N> >{
private:
	SMatrix<ct_binomial_coefficient(M+1, 3), 3, int> triangles = combinations<3, M+1>();
public:
	// some methods have to be empty since no normalization is needed
	void update_max() {}
	template<typename... UpdateArgs>
	void update_min(UpdateArgs&&... update_args) {}
	// rewrite the operetor() with no ormalization

	template<typename... CostArgs>
	double operator()(CostArgs&&... cost_args)
	{
		return get_cost(std::forward<CostArgs>(cost_args)...);
	}
	// since normalization is not needed, check_update always returns false
	bool check_update() {return false;}

	// method get_cost
	double get_cost(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete, 
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v, 
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		// all triangles in the element
		double max_cos1 = 0.0;
		// loop on elements to delete
		for(const auto & elem : elems_to_delete)
		{
			const auto& coords = elem.coords();
			for(unsigned i = 0; i < ct_binomial_coefficient(M+1, 3); ++i){
				double cos = get_max_cos(coords[triangles(i, 0)], coords[triangles(i, 1)], coords[triangles(i, 2)]);
			// double cos = get_max_cos(elem.coords());
				if(cos > max_cos1) {max_cos1 = cos;}
			}
		}
		// loop on elements to modify
		for(const auto & elem : elems_to_modify)
		{
			const auto& coords = elem.coords();
			for(unsigned i = 0; i < ct_binomial_coefficient(M+1, 3); ++i){
				double cos = get_max_cos(coords[triangles(i, 0)], coords[triangles(i, 1)], coords[triangles(i, 2)]);
			// double cos = get_max_cos(elem.coords());
				if(cos > max_cos1) {max_cos1 = cos;}
			}
		}
		double max_cos2 = 0.0;
		// loop on elements modified
		for(const auto & elem : elems_modified)
		{
			const auto& coords = elem.coords();
			for(unsigned i = 0; i < ct_binomial_coefficient(M+1, 3); ++i){
				double cos = get_max_cos(coords[triangles(i, 0)], coords[triangles(i, 1)], coords[triangles(i, 2)]);
			// double cos = get_max_cos(elem.coords());
				if(cos > max_cos2) {max_cos2 = cos;}
			}
		}
		double csi = (1.-max_cos1)*(1.-max_cos1)/( (1.-max_cos2)*(1.-max_cos2) );
		if(csi < DOUBLE_TOLERANCE) {return 0.;}
		return 1./std::tanh(csi) - 1./csi;
	}



};


} // core
} // fdapde

#endif // __SHARP_ELEMS_COST_H__