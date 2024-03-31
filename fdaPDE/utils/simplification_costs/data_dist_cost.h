#ifndef __DATA_DIST_COST_H__
#define __DATA_DIST_COST_H__

namespace fdapde{
namespace core{

#include "../../mesh/element.h"
#include "../../mesh/projection.h"
#include "../symbols.h"


struct DataDistCost{

	double operator()(const std::vector<Element<2, 3>> & elems_to_modify, 
					  const std::vector<Element<2, 3>> & elems_to_delete, 
					  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v, 
					  const std::unordered_set<unsigned> & data_ids ) const
	{
		return get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids)/max;
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

	void update_max(const std::vector<Element<2, 3>> & elems_to_modify, 
					  const std::vector<Element<2, 3>> & elems_to_delete, 
					  const std::vector<Element<2, 3>> & elems_modified, const SVector<3> & v, 
					  const std::unordered_set<unsigned> & data_ids )
	{
		double cost = get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids);
		if(cost > max) {max = cost;}
	}

	// Contiene le posizioni iniziali di tutti i dati
	DMatrix<double> data_;
	Simplification<2, 3>* p_simp_;
	double max = 0.0;

};


} // core
} // fdapde

#endif // __DATA_DIST_COST_H__