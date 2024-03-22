#ifndef __SIMPLIFICATION_H__
#define __SIMPLIFICATION_H__

#include "mesh.h"
#include "StructuredGridSearch.h"
#include "connections.h"

namespace fdapde{
namespace core{
template<int M, int N>
class Simplification{
	using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
private:
	StructuredGridSearch sgs_;
	Connections connections_;
	std::vector<Element<M, N>> elems_vec_;
	DMatrix<double> nodes_;
	DMatrix<int> elems_mat_;
	std::vector< std::array<int, Mesh<M, N>::n_vertices_per_facet> > facets;

	DMatrix<double> data_;
	ConnectionsType data_to_elems_;
    ConnectionsType elem_to_data_;

    std::multimap<double, std::pair<unsigned, SVector<N> >> costs_map;
    // mappa utile per l'update di collapse_costs dopo un collapse
    std::unordered_map<unsigned, double> facets_cost;

    unsigned n_nodes_; 
public:
	Simplification(const Mesh<M, N> & mesh);
};

// ===============
// implementazione
// ===============

// costruttore 

template<int M, int N>
Simplification<M, N>::Simplification<M, N>(const Mesh<M, N> & mesh):
	sgs_(mesh), connections_(mesh), elems_mat_(mesh.elements()), nodes_(mesh.nodes()),
	elem_to_data_(elems_mat_.rows()), n_nodes_(mesh.n_nodes())
{
	// vengono costruite le facce
	facets.reserve(mesh.n_facets());
    for(unsigned facet_id = 0; facet_id < mesh.n_facets(); ++facet_id)
        facets_.push_back(mesh.facet(facet_id).node_ids());
    // vengono costruite le connessioni elementi - dati
    data_ = mesh.nodes();
    data_to_elems_.resize(data_.rows());
    for(unsigned i = 0; i < data.rows(); ++i) // loop su tutti i dati
    {
        auto tmp_set = connections_.get_node_to_elems(i);
        data_to_elems_[i] = tmp_set;
        for(unsigned elem_id : tmp_set) // loop sugli elementi connessi al dato i
            elem_to_data_[elem_id].insert(i);
    }

}



} // fdapde
} // core

#endif // __SIMPLIFICATION_H__