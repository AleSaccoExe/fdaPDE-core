#ifndef __COST_OBJ_BASE_H__
#define __COST_OBJ_BASE_H__

#include <utility>

namespace fdapde{
namespace core{

template<int M, int N, typename ConcreteCost>
struct CostObjBase{
	// variable to be used for the max over the mesh
	double max_ = 0.0;
	double min_ = std::numeric_limits<double>::max();
	double threshold_ = 1.5;

	// methods

	void setup(Simplification<M, N> * p_simp){}

	template<typename... CostArgs>
	double operator()(CostArgs&&... cost_args)
	{  	
		double cost = static_cast<ConcreteCost&>(*this).get_cost(std::forward<CostArgs>(cost_args)...);
		if(cost < min_) {min_ = cost;}
		double ret = cost/max_;
		return ret;
	}

	template<typename... CostArgs>
	void update_min(CostArgs&&... cost_args)
	{
		double cost = static_cast<ConcreteCost&>(*this).get_cost(std::forward<CostArgs>(cost_args)...);
		if(cost < min_) {min_ = cost;}
	}

	void update_max()
	{
		if(min_ > max_) {max_ = min_;}
		min_ = std::numeric_limits<double>::max();
	}

	template<typename... UpdateArgs>
	void update(UpdateArgs&&...) {}

	bool check_update()
	{
		double current_min = min_;
		min_ = std::numeric_limits<double>::max();
		if(current_min > threshold_*max_) { return true; }
		return (current_min > threshold_*max_);
	}

	void set_threshold(double new_threshold) {threshold_ = new_threshold;}

};


} // namespace core
} // namespace fdapde

#endif