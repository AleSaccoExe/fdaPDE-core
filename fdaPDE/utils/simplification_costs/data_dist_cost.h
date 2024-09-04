#ifndef __DATA_DIST_COST_H__
#define __DATA_DIST_COST_H__

namespace fdapde{
namespace core{

#include "../../mesh/element.h"
#include "../../mesh/projection.h"
#include "../symbols.h"
#include "CostObjBase.h"


/*struct DataDistCost{

	double operator()(const std::vector<Element<2, 3>> & elems_to_modify, 
					  const std::vector<Element<2, 3>> & elems_to_delete, 
					  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v, 
					  const std::unordered_set<unsigned> & data_ids ) 
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
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		double max_dist = 0.0;
		for(unsigned datum_id : data_ids)
		{
			SVector<3> p_datum = project(elems_modified, p_simp_->get_data().row(datum_id));
			SVector<3> datum_original_position = data_.row(datum_id);
			double dist = (p_datum - datum_original_position).norm();
			if(dist > max_dist)
				max_dist = dist;
		}
		return max_dist;
	}
	void setup(Simplification<2, 3> * p_simp)
	{
		p_simp_ = p_simp;
		data_ = p_simp->get_data(); 
	}

	void update_min(const std::vector<Element<2, 3>> & elems_to_modify, 
					  const std::vector<Element<2, 3>> & elems_to_delete, 
					  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v, 
					  const std::unordered_set<unsigned> & data_ids )
	{
		double cost = get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids);
		if(cost < min_) {min_ = cost;}
	}

	void update_max()
	{
		if(min_>max_) {max_ = min_;}
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

	// Contiene le posizioni iniziali di tutti i dati
	DMatrix<double> data_;
	Simplification<2, 3>* p_simp_;
	double threshold_ = 1.5;
	double min_ = std::numeric_limits<double>::max();
	double max_ = 0.0;

};*/

class DataDistCost: public CostObjBase<2, 3, DataDistCost>{
private:
	DMatrix<double> data_;
	Simplification<2, 3>* p_simp_;
public:

	void setup(Simplification<2, 3> * p_simp)
	{
		p_simp_ = p_simp;
		data_ = p_simp->get_data(); 
	}

	double get_cost(const std::vector<Element<2, 3>> & elems_to_modify, 
					  const std::vector<Element<2, 3>> & elems_to_delete, 
					  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v, 
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		double max_dist = 0.0;
		for(unsigned datum_id : data_ids)
		{
			SVector<3> p_datum = project(elems_modified, p_simp_->get_data().row(datum_id));
			SVector<3> datum_original_position = data_.row(datum_id);
			double dist = (p_datum - datum_original_position).norm();
			if(dist > max_dist)
				max_dist = dist;
		}
		return max_dist;
	}

};


} // core
} // fdapde

#endif // __DATA_DIST_COST_H__