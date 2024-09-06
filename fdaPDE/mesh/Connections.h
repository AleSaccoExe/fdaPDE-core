#ifndef __CONNECTIONS_H__
#define __CONNECTIONS_H__

#include "../utils/symbols.h"
#include <set>
#include <unordered_set>
#include <vector>
#include <algorithm>


// Class Connections

// This class is used by Simplification to modify and update
// the mesh' connectivity during the simplification process.

namespace fdapde{
namespace core{


class Connections{
	// The type used for all the connections:
	using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
	// A simple alias to reference an unordered set:
	using SetType = std::unordered_set<unsigned>;
private:
	// Private members:
	ConnectionsType node_to_nodes;
	ConnectionsType node_to_elems;
	// Member to identify a facets by the ids of its vertices:
	std::unordered_map<std::set<int>, unsigned, std_set_hash<int>> facets;
	std::set<unsigned> active_nodes;
	std::set<unsigned> active_elements;

	// Private methods:
	// Temporary method (not used by Simplification)
	std::unordered_set<unsigned> replace_node_in_node_to_nodes(unsigned , unsigned , const std::unordered_set<unsigned>& involved);
	// The nodes connected to old_id are taken from it and connected to new_id
	void replace_node_in_node_to_nodes(unsigned old_id, unsigned new_id)
	{
		active_nodes.erase(old_id);
		// Take the nodes connected to old_id:
		const auto & old_conn = node_to_nodes[old_id];
		// Now these nodes have to be connected to new_id:
		for(unsigned id_node : old_conn)
		{
			node_to_nodes[new_id].insert(id_node);
			node_to_nodes[id_node].erase(old_id);
			node_to_nodes[id_node].insert(new_id);
		}
		// erase the self connection new_id-new_id
		node_to_nodes[new_id].erase(new_id);
	} // replace_node_in_node_to_nodes

	// The elements connected to the node old_id are taken from it and connected to new_id
	// It is assumed that the elements to erase are already erased
	std::unordered_set<unsigned> replace_node_in_node_to_elems(unsigned old_id, unsigned new_id)
	{
		// Extract node-element connections for newId
		auto & new_id_conn = node_to_elems[new_id];
		auto & old_id_conn = node_to_elems[old_id];
		// Move the connected from oldId to newId
		for(unsigned id_elem : old_id_conn)
			node_to_elems[new_id].insert(id_elem);
		// Make oldId inactive
		active_nodes.erase(old_id);
		return new_id_conn;
	} // replace_node_in_node_to_elems

	// Erase from the connections the elements passed in input
	template<typename Element>
	void erase_elems_in_node_to_elems(const std::vector<Element> & to_remove)
	{
		// assert(to_remove.size() == 2);
		// Loop over all the elements in input:
		for (const Element & elem : to_remove)
		{
			// Take the ids of the current element's vertices 
			std::array<int, Element::n_vertices> node_ids = elem.node_ids();
			// Set element inactive in element-node connections
			active_elements.erase(elem.ID());
			// Erase element from node-element connections
			for (unsigned j = 0; j < Element::n_vertices; ++j)
				node_to_elems[node_ids[j]].erase(elem.ID());
		}
	} // erase_elems_in_node_to_elems

	// Update the facets assuming to collapse the input facet
	// In the std::pair returned as output:
	// first) ids of the facets erased by the collapse
	// second) ids of the facets modified by the collapse
	template<typename FacetType>
	std::pair<SetType, SetType> update_facets(const FacetType facet)
	{
		// If the facet in input is collaped, then some update of the neighbouring facets are needed:
		// teh facets connected only to the node facet[0] need no update, the facets connected to any other
		// vertex of the input facet have to be modified: one of their vertex has to be replaced by the node
		// with id facet[0].
		// The facets sharing more then one vertex with the input facet have to be deleted.

		SetType to_erase;
		SetType to_modify;
		unsigned collapsing_node = facet[0];
		// Case of tetrhaedra
		if constexpr(facet.size() == 3)
		{
			for(unsigned i = 0; i < 3; ++i)
			{
				const auto & conn_nodes_i = node_to_nodes[facet[i]]; // nodes connected to i
				for(unsigned j = i+1; j < 3; ++j)
				{
					std::set<int> tmp_facet = {static_cast<int>(facet[i]), static_cast<int>(facet[j])};
					std::unordered_set<unsigned> conn_nodes = node_to_nodes[facet[j]]; // nodes connected to j
					// Intersection between the tqo sets:
					for (auto it = conn_nodes.begin(); it != conn_nodes.end();) {
			    		if (!conn_nodes_i.count(*it)) { it = conn_nodes.erase(it); }
			    		else              { ++it; }
					}
					// Now the facets to erase are inserted: (i, j, ...)
					for(unsigned conn_node : conn_nodes)
					{
						tmp_facet.insert(conn_node);
						if(facets.find(tmp_facet)!=facets.end()){
							to_erase.insert(facets.at(tmp_facet));
							facets.erase(tmp_facet);	
						}
						tmp_facet.erase(conn_node);
					}

				}
			}
			for(unsigned i = 1; i < 3; ++i)
			{
				// questa è la faccia da controllare
				// i nodi connessi servono dentro un vettore, dato che ora devo controllare tutte le possibili
				// coppie di nodi
				SetType tmp_set(node_to_nodes[facet[i]].begin(), node_to_nodes[facet[i]].end());
				tmp_set.erase(collapsing_node);
				std::vector<unsigned> conn_nodes(tmp_set.begin(), tmp_set.end());
				for(unsigned j = 0; j < conn_nodes.size(); ++j)
				{
					for(unsigned k = j+1; k < conn_nodes.size(); ++k)
					{
						std::set<int> tmp_facet = {facet[i]};
						tmp_facet.insert(conn_nodes[j]);
						tmp_facet.insert(conn_nodes[k]);
						if(facets.find(tmp_facet)!=facets.end()) // If the three nodes in tmp_facet actually form an existing facet...
						{
							unsigned old_id = facets.at(tmp_facet);
							facets.erase(tmp_facet);
							tmp_facet.erase(facet[i]);
							tmp_facet.insert(collapsing_node);
							if(facets.find(tmp_facet)!=facets.end()) // The facet already exists, then old_id has to be deleted
								to_erase.insert(old_id);
							else{
								facets[tmp_facet] = old_id;
								to_modify.insert(old_id);
							}

						}
					}
				}
				
			}
		} // if constexpr(facet.size()==3)
		else // Case of triangles
		{
			std::set<unsigned> conn_nodes(node_to_nodes[facet[1]].begin(), node_to_nodes[facet[1]].end());
			assert(conn_nodes.find(collapsing_node)!=conn_nodes.end());
			conn_nodes.erase(collapsing_node);
			for(unsigned id_node : conn_nodes)
			{
				std::set<int> tmp_facet = {facet[1]};
				tmp_facet.insert(id_node);
				assert(facets.find(tmp_facet)!=facets.end());
				unsigned old_id = facets.at(tmp_facet);
				facets.erase(tmp_facet);
				tmp_facet.erase(facet[1]);
				tmp_facet.insert(collapsing_node);
				if(facets.find(tmp_facet)!=facets.end())
					to_erase.insert(old_id);
				else{
					facets[tmp_facet] = old_id;
					to_modify.insert(old_id);
				}
			}
			facets.erase({facet.begin(), facet.end()});
		}
		// unsigned facet_id = facets.at({facet.begin(), facet.end()});
		// to_erase.insert(facet_id);
		return {to_erase, to_modify};
	} // update_facets

public:
	// Public methods:
	Connections() = default; 
	template<typename Mesh_>
	Connections(Mesh_&& mesh)
	{
		using Mesh = std::decay_t<Mesh_>;
		constexpr unsigned M = Mesh::local_dimension;
		constexpr unsigned N = Mesh::embedding_dimension;
		constexpr unsigned n_vertices = Mesh::n_vertices;
		const auto & elements = mesh.elements();
		unsigned n_nodes = mesh.n_nodes();
		unsigned n_elements = mesh.n_elements();
		unsigned n_facets = mesh.n_facets();

		node_to_nodes.reserve(n_nodes);
		node_to_elems.reserve(n_elements);

		for(unsigned id_node = 0; id_node < n_nodes; ++id_node)
		{
			node_to_elems.push_back({});
			node_to_nodes.push_back({});
			active_nodes.insert(id_node);
		}			
		// The connections are built:
		// Loop over all the elements in the mesh:
		for (int id_elem = 0; id_elem < mesh.n_elements(); ++id_elem)
		{
			active_elements.insert(id_elem);
			// Extract element:
			const auto & elem = mesh.element(id_elem);
			// Extract the ids of the element's vertices:
			std::array<int, ct_nvertices(M)> node_ids = elem.node_ids();
			// Loop over the current element's vertices:
			for(unsigned i = 0; i < ct_nvertices(M); ++i)
			{
				// insert the connection node_ids[i] - id_elem in node_to_elems:
				node_to_elems[node_ids[i]].insert(id_elem);
				// Insert now the node-nodes connections:
				for(unsigned j = i+1; j < ct_nvertices(M); ++j)
				{
					// Insert the connections node_ids[i] - node_ids[j] in node_to_nodes:
					node_to_nodes[node_ids[i]].insert(node_ids[j]);
					node_to_nodes[node_ids[j]].insert(node_ids[i]);
				}
			}
		}

		// Add the facets with their id:
		for(unsigned id_facet = 0; id_facet < n_facets; ++id_facet){
			auto node_ids = mesh.facet(id_facet).node_ids();
			facets[{node_ids.begin(), node_ids.end()}] = id_facet;
		}
	} // constructor
	// Temporary method (not used by Simplification)
	void refresh();


	// Find the nodes that insist on the input facet
	template<typename FacetType>
	std::unordered_set<unsigned> nodes_on_facet(const FacetType & facet) const
	{
		// For each vertex of the facet, the ids of the nodes connected to 
		// it are taken, and then the set intersection is computed:
		std::unordered_set<unsigned> conn_nodes = node_to_nodes[facet[0]];
		for(unsigned i = 1; i < facet.size(); ++i){
			for (auto it = conn_nodes.begin(); it != conn_nodes.end();) {
			    if (!node_to_nodes[facet[i]].count(*it)) { it = conn_nodes.erase(it); }
			    else              { ++it; }
			}
		}
		return conn_nodes;
	}// nodes_on_facet

	// find all the nodes involved in the collapse of the input facet
	template<typename FacetType>
	std::unordered_set<unsigned> nodes_involved_in_collapse(const FacetType & facet) const
	{
		// The nodes involved in collapse are all the nodes connected to the input facet
		std::unordered_set<unsigned> conn_nodes;
		for(unsigned id_node : facet){
			const auto & tmp_conn = node_to_nodes[id_node];
			conn_nodes.insert(tmp_conn.begin(), tmp_conn.end());
		}
		for(unsigned id_node : conn_nodes)
			conn_nodes.erase(id_node);
		return conn_nodes;
	} // nodes_involved_in_collapse

	// Find all the elements that are erased by the collapse of the input facet
	template<typename FacetType>
	std::unordered_set<unsigned> elems_erased_in_collapse(const FacetType & facet) const
	{/*
		// intersezione degli elementi connessi ai nodi di facet
		std::unordered_set<unsigned> conn_elems = node_to_elems[facet[0]];
		for(unsigned i = 1; i < facet.size(); ++i){
			for (auto it = conn_elems.begin(); it != conn_elems.end();) {
			    if (!node_to_elems[facet[i]].count(*it)) { it = conn_elems.erase(it); }
			    else              { ++it; }
			}
		}*/

		// The elements erased by the collapse are all the elements that share at least two
		// vertices with the input facet
		std::unordered_set<unsigned> conn_elems; 
		for(int i = 0; i < facet.size(); ++i){
			for(int j = i+1; j < facet.size(); ++j){
				// intersezione tra gli elementi connessi a facet[i] e quelli connessi a facet[j]
				std::unordered_set<unsigned> tmp_set = node_to_elems[facet[i]];
				for (auto it = tmp_set.begin(); it != tmp_set.end();) {
		    		if (!node_to_elems[facet[j]].count(*it)) { it = tmp_set.erase(it); }
		    		else              { ++it; }
				}

				conn_elems.insert(tmp_set.begin(), tmp_set.end());

			}
		}
		return conn_elems;
	} // elems_erased_in_collapse

	// Find the elements that are modified by the collapse of the input facet
	template<typename FacetType>
	std::unordered_set<unsigned> elems_modified_in_collapse(const FacetType & facet) const
	{
		// The elements that are modified by the collapse are all the elements that share
		// only one veertex with th einput facet.
		// Take the elements erased:
		auto to_erase = elems_erased_in_collapse(facet);
		// Take all the elements connected to the input facet:
		std::unordered_set<unsigned> conn_elems;
		for(unsigned id_node : facet){
			const auto & tmp_conn = node_to_elems[id_node];
			conn_elems.insert(tmp_conn.begin(), tmp_conn.end());
		}
		// The elements modified are all the elements connected to the facet minus
		// the elements that are erased by the collapse:
		for(unsigned id_elem : to_erase)
			conn_elems.erase(id_elem);
		return conn_elems;
	}// elems_modified_in_collapse

	// Find the nodes connected on two depth levels to the input node
	std::unordered_set<unsigned> extended_node_patch(unsigned id_node) const;
	// Find all the facets whose cost has to be computed again
	std::set<unsigned> facets_to_update(unsigned node_id) const
	{
		std::set<unsigned> facet_ids;
		const auto & node_ids = node_to_nodes[node_id]; // noeds connected to node_id
		if(facets.begin()->first.size() == 2) // if the saved facets have two nodes
		{
			for(unsigned node_id2 : node_ids)
			{
				facet_ids.insert(facets.at({static_cast<int>(node_id), static_cast<int>(node_id2)}));
				const auto & node_ids2 = node_to_nodes[node_id2];
				for(unsigned node_id3 : node_ids2)
					facet_ids.insert(facets.at({static_cast<int>(node_id2), static_cast<int>(node_id3)}));
			}
			return facet_ids;
		}
		else // the saved facets have three nodes (3D mesh)
		{
			for(unsigned node_id2 : node_ids)
			{
				std::set<int> tmp_facet = {static_cast<int>(node_id), static_cast<int>(node_id2)};
				const auto & node_ids2 = node_to_nodes[node_id2];
				for(unsigned node_id3 : node_ids2) // for each node id3 connected to id2...
				{
					tmp_facet.insert(static_cast<unsigned>(node_id3));
					if(facets.find(tmp_facet)!=facets.end()) // ...check if facet (id, id2, id3) exists
						{ facet_ids.insert(facets.at(tmp_facet)); } // if it exists, it is added to the output set
					tmp_facet.erase(static_cast<int>(node_id3));
					const auto & node_ids3 = node_to_nodes[node_id3];
					for(unsigned node_id4 : node_ids3) // check the facets of the type: (id3, id4, ...)
					{
						std::set<int> tmp_facet2 = {static_cast<int>(node_id3), static_cast<int>(node_id4)};
						const auto & node_ids4 = node_to_nodes[node_id4];
						for(unsigned node_id5 : node_ids4)
						{
							tmp_facet2.insert(static_cast<int>(node_id5));
							if(facets.find(tmp_facet2)!=facets.end())
								{ facet_ids.insert(facets.at(tmp_facet2)); }
							tmp_facet2.erase(node_id5);
						}
					}
				}
			}
			return facet_ids;
		}
		return {}; 
	} // facets_to_update

	template<typename Element>
	std::unordered_set<unsigned> element_patch(const Element & elem) const
	{
		auto node_ids = elem.node_ids();
		std::unordered_set<unsigned> conn_elems;
		for(unsigned id_node : node_ids){
			const auto & tmp_set = node_to_elems[id_node];
			conn_elems.insert(tmp_set.begin(), tmp_set.end());
		}
		conn_elems.erase(elem.ID());
		return conn_elems;
	} // element_patch

	// getters
	std::unordered_set<unsigned> get_node_to_nodes(unsigned id_node) const{return node_to_nodes[id_node];}
	std::unordered_set<unsigned> get_node_to_elems(unsigned id_node) const{return node_to_elems[id_node];}
	std::set<unsigned> get_active_nodes() const{return active_nodes;}
	std::set<unsigned> get_active_elements() const{return active_elements;}


	// Modify all the connectivity needed after the collapse of the input facet
	template<typename Element, typename FacetType>
	std::pair<SetType, SetType> collapse_facet(const FacetType & facet, const std::vector<Element> & to_remove)
	{
		// The facet in input is collapsed into the node with id facet[0]
		// assert(facets.find({facet.begin(), facet.end()})!=facets.end());
		unsigned collapsing_node = facet[0];
		// update the facets:
		auto facets_info = update_facets(facet);
		// erase the elements:
		erase_elems_in_node_to_elems(to_remove);
		// replace the nodes facet[i] with the node facet[0]
		for(unsigned i = 1; i < facet.size(); ++i){
			replace_node_in_node_to_nodes(facet[i], collapsing_node);
			replace_node_in_node_to_elems(facet[i], collapsing_node);
		}
		for(unsigned i = 1; i < facet.size(); ++i)
			node_to_nodes[collapsing_node].erase(facet[i]);
		return facets_info;
	} // collapse_facet

};

// ======================
// IMPLEMENTATIVE DETAILS
// ======================



void Connections::refresh()
{
	// mappa dal vecchio id del nodo al nuovo id del nodo
	std::map<unsigned, unsigned> old_to_new_node_ids;
	// mappa dal vecchio id dell'elemento al nuovo id dell'elemento
	std::map<unsigned, unsigned> old_to_new_el_ids;
	unsigned new_id = 0;
	for(unsigned old_node_id : active_nodes){
		old_to_new_node_ids[old_node_id] = new_id;
		new_id++;
	}
	new_id = 0;
	for(unsigned old_el_id : active_elements){
		old_to_new_el_ids[old_el_id] = new_id;
		new_id++;
	}

	// qui salvo le nuove connessioni da sostituire poi nei membri della classe
	ConnectionsType new_node_to_nodes;
	ConnectionsType new_node_to_elems;

	new_id = 0; // il nuovo id del nodo comincia da 0
	for(unsigned i : active_nodes)
	{
		new_node_to_elems.push_back({});
		new_node_to_nodes.push_back({});
		// dentro a new_node_to_elems vengono sostituiti i vecchi id degli elementi
		std::unordered_set<unsigned> old_node_to_elems = node_to_elems[i];
		for(unsigned old_el_id : old_node_to_elems)
			new_node_to_elems[new_id].insert(old_to_new_el_ids[old_el_id]);
		
		// dentro a new_node_to_node vengono sostituiti i vecchi id dei nodi
		std::unordered_set<unsigned> old_node_to_nodes = node_to_nodes[i];
		for(unsigned old_node_id : old_node_to_nodes)
			new_node_to_nodes[new_id].insert(old_to_new_node_ids[old_node_id]);
		
		new_id++;
	}

	// ora manca solo da sostituire i nuovi vettori a quelli della classe
	node_to_nodes = new_node_to_nodes;
	node_to_elems = new_node_to_elems; 
} // refresh






} // namespace fdapde
} // namespace core
#endif // __CONNECTIONS_H__