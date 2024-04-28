#ifndef __SHARP_ELEMS_COST_H__
#define __SHARP_ELEMS_COST_H__

#include "../../mesh/element.h"
#include "../symbols.h"
#include <cmath>

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

template<int N>
double get_max_cos(const std::array<SVector<N>, 3> & points)
{
	double max = 0.0;
	SVector<N> A = points[0];
	SVector<N> B = points[1];
	SVector<N> C = points[2];
	double cos1 = std::abs((B-A).dot(C-A))/( (B-A).norm()*(C-A).norm() );
	if(max < cos1) {max = cos1;}
	double cos2 = std::abs((C-B).dot(A-B))/( (C-B).norm()*(A-B).norm() );
	if(max < cos2) {max = cos2;}
	double cos3 = std::abs((A-C).dot(B-C))/( (A-C).norm()*(B-C).norm() );
	if(max < cos3) {max = cos3;}
	return max;

}


template<int M, int N>
struct SharpElemsCost{
	double operator()(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete, 
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v, 
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		return get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids)/max_;
	}

	double get_cost(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete, 
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v, 
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		/*
		double max_ratio = 0.0;
		for(const auto & elem : elems_modified)
		{
			double radius = circumradius(elem.coords());
			const auto & coords = elem.coords();
			double min_edge = std::numeric_limits<double>::max();
			for(unsigned i = 0; i < M+1; ++i)
				for(unsigned j = i+1; j < M+1; ++j)
				{
					double edge_length = (coords[i]-coords[j]).norm();
					if(edge_length < min_edge) {min_edge = edge_length;}
				}
			double ratio = radius/min_edge;
			if(ratio > max_ratio) {max_ratio = ratio;}
		}
		return max_ratio;
		*/
		double max_cos1 = 0.0;
		for(const auto & elem : elems_to_delete)
		{
			double cos = get_max_cos(elem.coords());
			if(cos > max_cos1) {max_cos1 = cos;}
		}
		for(const auto & elem : elems_to_modify)
		{
			double cos = get_max_cos(elem.coords());
			if(cos > max_cos1) {max_cos1 = cos;}
		}
		double max_cos2 = 0.0;
		for(const auto & elem : elems_modified)
		{
			double cos = get_max_cos(elem.coords());
			if(cos > max_cos2) {max_cos2 = cos;}
		}
		double csi = (1.-max_cos1)*(1.-max_cos1)/( (1.-max_cos2)*(1.-max_cos2) );
		return std::tanh(csi) - 1./csi;
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

};


} // core
} // fdapde

#endif // __SHARP_ELEMS_COST_H__