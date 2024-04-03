#ifndef __SHARP_ELEMS_COST_H__
#define __SHARP_ELEMS_COST_H__

#include "../../mesh/element.h"
#include "../symbols.h"

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
	}
	
	void setup(Simplification<M, N> * p_simp){}

	void update_max(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete, 
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v, 
					  const std::unordered_set<unsigned> & data_ids )
	{
		double cost = get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids);
		if(cost > max_) {max_ = cost;}
	}

	void update(const std::vector<Element<2, 3>> & elems_to_delete,
				const std::vector<Element<2, 3>> & elems_modified) {}

	double max_ = 0.0;
};


} // core
} // fdapde

#endif // __SHARP_ELEMS_COST_H__