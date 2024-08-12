#ifndef __GEOM_COST_H__
#define __GEOM_COST_H__

#include <vector>

#include "../../mesh/element.h"
#include "../symbols.h"

namespace fdapde{
namespace core{

// struct implementing the geometric cost for 2.5D meshes

struct GeomCost{
	// descrizione input:
	// - elems: elementi involved in collapse
	// - v: punto del collapse
	double operator()(const std::vector<Element<2, 3>> & elems_to_modify, 
					  const std::vector<Element<2, 3>> & elems_to_delete,
					  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v,
					  const std::unordered_set<unsigned> & data_ids)
	{  	
		double cost = get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids);
		// std::cout<<"cost disp not normalized: "<<cost<<"\n";
		if(cost < min_) {min_ = cost;}
		double ret = cost/max_;
		// std::cout<<"max disp: "<<max_<<"\n";
		return ret;
	}
	double get_cost(const std::vector<Element<2, 3>> & elems_to_modify, 
					  const std::vector<Element<2, 3>> & elems_to_delete,
					  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v,
					  const std::unordered_set<unsigned> & data_ids) const
	{
		SVector<10> Q;
		Q.setZero();
		for(const auto & elem : elems_to_modify)
		{
			SVector<3> n = elem.hyperplane().normal(); // normale all'elemento
			double d = -n.dot(elem.coords()[0]);
			// Construct matrix K
			SVector<10> K;
			K <<  n[0]*n[0], n[0]*n[1], n[0]*n[2], n[0]*d, n[1]*n[1], n[1]*n[2], n[1]*d, n[2]*n[2], n[2]*d, d*d;
			Q += K;
		}
		for(const auto & elem : elems_to_delete)
		{
			SVector<3> n = elem.hyperplane().normal(); // normale all'elemento
			double d = -n.dot(elem.coords()[0]);
			// Construct matrix K
			SVector<10> K;
			K <<  n[0]*n[0], n[0]*n[1], n[0]*n[2], n[0]*d, n[1]*n[1], n[1]*n[2], n[1]*d, n[2]*n[2], n[2]*d, d*d;
			Q += 2*K;
		}
		return Q[0]*v[0]*v[0] + Q[4]*v[1]*v[1] + Q[7]*v[2]*v[2]
			+ 2*Q[1]*v[0]*v[1] + 2*Q[2]*v[0]*v[2] + 2*Q[5]*v[1]*v[2]
			+ 2*Q[3]*v[0] + 2*Q[6]*v[1] + 2*Q[8]*v[2] + Q[9];
	}
	void setup(Simplification<2, 3> * p_simp){}

	void update_min(const std::vector<Element<2, 3>> & elems_to_modify, 
				  const std::vector<Element<2, 3>> & elems_to_delete,
				  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v,
				  const std::unordered_set<unsigned> & data_ids)
	{
		double cost = get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids);
		if(cost < min_) {min_ = cost;}
	}

	void update_max()
	{
		if(min_ > max_) {max_ = min_;}
		min_ = std::numeric_limits<double>::max();
	}

	void update(const std::vector<Element<2, 3>> & elems_to_delete,
				const std::vector<Element<2, 3>> & elems_modified) {}
	bool check_update()
	{
		// return false;
		double current_min = min_;
		min_ = std::numeric_limits<double>::max();
		if(current_min > threshold_*max_) { return true; }
		return (current_min > threshold_*max_);
	}
	void set_threshold(double new_threshold) {threshold_ = new_threshold;}

	double max_ = 0.0;
	double min_ = std::numeric_limits<double>::max();
	double threshold_ = 1.5;
};

}
}


#endif // __GEOM_COST_H__