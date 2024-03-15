#ifndef __DATA_CONNECTIONS_H__
#define __DATA_CONNECTIONS_H__

#include <vector>
#include <unordered_set>
#include "../utils/symbols.h"

class DataConnections{
	using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
private:
	ConnectionsType data_to_elems_;
	ConnectionsType elem_to_data_;
	DMatrix<double> & data_;
	unsigned n_data_;
public:
	template<typename Mesh_>
	DataConnections(Mesh_&& mesh, DMatrix<double> & data);

};

// ===============
// IMPLEMENTAZIONE
// ===============

// COSTRUTTORE

template<typename Mesh_>
DataConnections::DataConnections(Mesh_&& mesh, DMatrix<double> & data):data_(data), n_data_(data.rows()), 
													  data_to_elems_(data.rows()), elem_to_data_(mesh.n_elements())
{
	unsigned n_elems = mesh.n_elements();
	auto locations = mesh.locate(data_);
	for(unsigned i = 0; i < n_data_; ++i)
	{
		data_to_elems_[i].insert(locations[i]);
		elem_to_data_[locations[i]].insert(i);
	}
}



#endif // __DATA_CONNECTIONS_H__