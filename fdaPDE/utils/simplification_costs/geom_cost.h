#ifndef __GEOM_COST_H__
#define __GEOM_COST_H__

#include <vector>

#include "CostObjBase.h"
#include "../../mesh/element.h"
#include "../symbols.h"

namespace fdapde{
namespace core{

// struct implementing the geometric cost for 2.5D meshes



struct GeomCost: public CostObjBase<2, 3, GeomCost>{
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
};

}
}


#endif // __GEOM_COST_H__