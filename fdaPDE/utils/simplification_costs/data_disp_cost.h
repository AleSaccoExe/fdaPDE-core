#ifndef __DATA_DISP_COST_H__
#define __DATA_DISP_COST_H__

#include "../../mesh/simplification.h"
#include "../symbols.h"
#include "CostObjBase.h"
#include <chrono>
#include <iostream>

namespace fdapde{
namespace core{

// struct implementing the cost on data distribution



template<int M, int N>
class DataEquiCost: public CostObjBase<M, N, DataEquiCost<M, N>>{
private:
	Simplification<M, N> * p_simp_;
	std::vector<double> qoi_;
	unsigned num_elems_;
	double mean_qoi_;

	double compute_qoi(unsigned elem_id) const
	{
		double Nt = 0.0;
		const auto & data_ids = p_simp_->elem_to_data(elem_id);
		for(unsigned datum_id : data_ids)
		{
			auto patch = p_simp_->data_to_elems(datum_id).size();
			
			if (patch == 1)
				Nt += 1.;
			else
				Nt += 1./patch;				
		}
		return Nt;
	}

public:

	void setup(Simplification<M, N> * p_simp)
	{
		p_simp_ = p_simp;
		num_elems_ = p_simp->active_elems().size();
		qoi_.reserve(num_elems_);
		double qoi_sum = 0.0;
		for(std::size_t elem_id = 0; elem_id < num_elems_; ++elem_id)
		{
			double qoi_elem = compute_qoi(elem_id);
			qoi_.push_back(qoi_elem);
			qoi_sum += qoi_elem;
		}
		mean_qoi_ = qoi_sum/num_elems_;

	} // setup

	double get_cost(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete,
					  const std::vector<Element<M, N>> & elems_modified, const SVector<N> & v,
					  const std::unordered_set<unsigned> & data_ids) const
	{
		double disp_cost = 0.0;
		// new element-data connections
		std::vector<std::unordered_set<unsigned>> new_elem_to_data = projection_info(elems_modified, p_simp_->get_data(),data_ids);
		std::unordered_map<unsigned, std::unordered_set<unsigned>> new_data_to_elem;
		for(unsigned i = 0; i < new_elem_to_data.size(); ++i)
		{
			auto& data_on_elem = new_elem_to_data[i];
			for(unsigned datum_id : data_on_elem) {new_data_to_elem[datum_id].insert(i);}
		}
		// for each element in elems_modified its qoi is computed
		for(unsigned i = 0; i < elems_modified.size(); ++i)
		{
			double Nt = 0.0;
			// take the data on the element
			const auto & data_on_elem = new_elem_to_data[i];
			for(unsigned datum_id : data_on_elem)
			{
				// the elements the data lied on before the collapse:
				auto patch = p_simp_->data_to_elems(datum_id);
				// erase from the patch the modified and the erased elements:
				for(const auto & elem : elems_to_modify) {patch.erase(elem.ID());}
				for(const auto & elem : elems_to_delete) {patch.erase(elem.ID());}
				// add the new connections:
				auto patch_size = patch.size() + new_data_to_elem.at(datum_id).size();
				if (patch_size == 1)
					Nt += 1.;
				else
					Nt += 1./patch_size;
				// assert(patch_size != 0);
			}
			disp_cost += (mean_qoi_ - Nt)*(mean_qoi_ - Nt);
		}
		return disp_cost/elems_modified.size();
	} // get_cost

	void update(const std::vector<Element<M, N>> & elems_to_delete,
				const std::vector<Element<M, N>> & elems_modified)
	{
		double sum_qoi = mean_qoi_*num_elems_;
		// erase from the total qoi the qoi of the deleted elements:
		for(const auto & elem : elems_to_delete) {sum_qoi = sum_qoi - qoi_[elem.ID()];}
		// moreover, erase the contribute of the modified elements ...
		for(const auto & elem : elems_modified) {
			sum_qoi = sum_qoi - qoi_[elem.ID()];
			// ... and now compute again their qoi:
			double new_qoi = compute_qoi(elem.ID());
			sum_qoi = sum_qoi + new_qoi;
			qoi_[elem.ID()] = new_qoi;
		}
		// update the number of elements:
		num_elems_ = p_simp_->active_elems().size();
		// update the mean value of the qoi:
		mean_qoi_ = sum_qoi/num_elems_;
	} // update

	double get_qoi(unsigned elem_id) const {return qoi_[elem_id];}

};

} // core
} // fdapde

#endif // __DATA_DISP_COST_H__