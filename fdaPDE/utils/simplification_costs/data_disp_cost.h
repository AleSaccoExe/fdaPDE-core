#ifndef __DATA_DISP_COST_H__
#define __DATA_DISP_COST_H__

#include "../../mesh/simplification.h"
#include "../symbols.h"

namespace fdapde{
namespace core{

template<int M, int N> 
struct DataDispCost{
	double operator()(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete,
					  const std::vector<Element<M, N>> & elems_modified, const SVector<3> & v,
					  const std::unordered_set<unsigned> & data_ids) const
	{ return get_cost(elems_to_modify, elems_to_delete,elems_modified, v, data_ids)/max_; }

	void update_min(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete,
					  const std::vector<Element<M, N>> & elems_modified, const SVector<3> & v,
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

	double get_cost(const std::vector<Element<M, N>> & elems_to_modify, 
					  const std::vector<Element<M, N>> & elems_to_delete,
					  const std::vector<Element<M, N>> & elems_modified, const SVector<3> & v,
					  const std::unordered_set<unsigned> & data_ids) const
	{
		double disp_cost = 0.0;
		// nuove connessioni elemento - dati
		std::vector<std::set<unsigned>> new_elem_to_data = projection_info(elems_modified, p_simp_->get_data(),data_ids);
		std::map<unsigned, std::set<unsigned>> new_data_to_elem;
		for(unsigned i = 0; i < new_elem_to_data.size(); ++i)
		{
			auto data_on_elem = new_elem_to_data[i];
			for(unsigned datum_id : data_on_elem) {new_data_to_elem[datum_id].insert(i);}
		}
		// adesso per ogni elemento in elems_modifed si calcola il qoi
		for(unsigned i = 0; i < elems_modified.size(); ++i)
		{
			double Nt = 0.0;
			// si prendono i dati sull'elemento
			const auto & data_on_elem = new_elem_to_data[i];
			for(unsigned datum_id : data_on_elem)
			{
				// gli elementi a cui il dato appartiene prima del collapse
				auto patch = p_simp_->data_to_elems(datum_id);
				// vengono eliminati dal patch gli elementi eliminati e modificati dal collapse
				for(const auto & elem : elems_to_modify) {patch.erase(elem.ID());}
				for(const auto & elem : elems_to_delete) {patch.erase(elem.ID());}
				// vengono aggiunte le nuove connessioni
				auto patch_size = patch.size() + new_data_to_elem.at(datum_id).size();
				if (patch_size == 1)
					Nt += 1.;
				else
					Nt += 1./patch_size;
				assert(patch_size != 0);
			}
			disp_cost += (mean_qoi_ - Nt)*(mean_qoi_ - Nt);
		}
		return disp_cost/elems_modified.size();
	} // get_cost

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

	void update(const std::vector<Element<M, N>> & elems_to_delete,
				const std::vector<Element<M, N>> & elems_modified)
	{
		double sum_qoi = mean_qoi_*num_elems_;
		// dalla somma dei qoi viene tolto il contributo degli elementi eliminati
		for(const auto & elem : elems_to_delete) {sum_qoi = sum_qoi - qoi_[elem.ID()];}
		// viene anche tolto il contributo degli elementi modificati ...
		for(const auto & elem : elems_modified) {
			sum_qoi = sum_qoi - qoi_[elem.ID()];
			// ... e viene quindi ricalcolato
			double new_qoi = compute_qoi(elem.ID());
			sum_qoi = sum_qoi + new_qoi;
			qoi_[elem.ID()] = new_qoi;
		}
		// viene aggiornato anche il numero degli elementi
		num_elems_ = p_simp_->active_elems().size();
		// infine è aggiornato il qoi medio
		mean_qoi_ = sum_qoi/num_elems_;
	} // update

	double compute_qoi(unsigned elem_id) const
	{
		double Nt = 0.0;
		const auto & data_ids = p_simp_->elem_to_data(elem_id);
		for(unsigned datum_id : data_ids)
		{
			auto patch = p_simp_->data_to_elems(datum_id).size();
			/*
			if (patch == 1)
				Nt += 1.;
			else
				Nt += 1./patch;
				*/
			Nt += 1./patch;
		}
		return Nt;
	}



	// membri
	double mean_qoi_;
	double max_ = 0.0;
	double min_ = std::numeric_limits<double>::max();
	Simplification<M, N> * p_simp_;
	std::vector<double> qoi_;
	unsigned num_elems_;
};

} // core
} // fdapde

#endif // __DATA_DISP_COST_H__