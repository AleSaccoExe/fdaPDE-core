#ifndef __CONNECTIONS_H__
#define __CONNECTIONS_H__

#include "mesh.h"

namespace fdapde{
namespace core{


class Connections{
	using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
	using SetType = std::unordered_set<unsigned>;
private:
	// comincio l'implementazione utilizzando dei vettori. Un elemento inattivo semplicemente non
	// sarà contenuto nelle mappe
	ConnectionsType node_to_nodes;
	ConnectionsType node_to_elems;

	std::set<unsigned> active_nodes;
	std::set<unsigned> active_elements;

	// std::vector<unsigned> facets {};
	ConnectionsType facets {};
    std::unordered_map<unsigned, std::vector<int>> facet_to_elems {};   // map from facet id to elements insisting on it
    ConnectionsType elem_to_facets {}; // COME FARE PER CALCOLARE A COMPILE TIME IL NUMERO DI FACCE?
    int n_facets_ = 0;

public:
	// metodi privati
	// modifica le connessioni nodo-nodo in modo da rimpiazzare il nodo id_old con il nodo id_new:
	// i nodi connessi a id_old saranno ora connessi a id_new
	std::unordered_set<unsigned> replace_node_in_node_to_nodes(unsigned , unsigned , const std::unordered_set<unsigned>& involved);
	std::unordered_set<unsigned> replace_node_in_node_to_nodes(unsigned old_id, unsigned new_id);
	std::unordered_set<unsigned> replace_node_in_node_to_elems(unsigned old_id, unsigned new_id);
	void replace_node_in_facet(unsigned old_id, unsigned new_id, unsigned id_facet);
	template<typename Element>
	void erase_elems_in_node_to_elems(const std::vector<Element> & to_remove); // HA SENSO??

//public:
	Connections() = default;
	template<typename Mesh>
	Connections(Mesh&& mesh);
	// refresh rinumera tutti i dati contenuti nelle strutture node_to_nodes, node_to_elems
	void refresh();
	// non va bene: il numero di nodi che compongono una facet è noto a runtime, conviene utilizzare un array,
	// ma come posso capire la dimensione corretta?

	std::unordered_set<unsigned> nodes_involved_in_edge_collapse(const std::vector<unsigned> & facet) const;

	// getters
	std::unordered_set<unsigned> get_node_to_nodes(unsigned id_node) const{return node_to_nodes[id_node];}
	std::unordered_set<unsigned> get_node_to_elems(unsigned id_node) const{return node_to_elems[id_node];}
	std::set<unsigned> get_active_nodes() const{return active_nodes;}
	std::set<unsigned> get_active_elements() const{return active_elements;}


	// implementazione temporanea
	// template<typename Element>
	// void collapse_facet(const unsigned id_facet, const std::vector<Element> & to_remove);


};

// ===============
// IMPLEMENTAZIONE
// ===============

// Costruttore 
template<typename Mesh_>
Connections::Connections(Mesh_&& mesh) 
{
	using Mesh = std::decay_t<Mesh_>;
	constexpr unsigned M = Mesh::local_dimension;
	constexpr unsigned N = Mesh::embedding_dimension;
	constexpr unsigned n_vertices_per_facet = Mesh::n_vertices_per_facet;
	constexpr unsigned n_elements_per_facet = Mesh::n_elements_per_facet;
	constexpr unsigned n_vertices = Mesh::n_vertices;
	const auto & elements = mesh.elements();
	unsigned n_nodes = mesh.n_nodes();
	unsigned n_elements = mesh.n_elements();

	// da utilizzare per la costruzione delle connessioni facce-elemento
	auto facet_pattern = combinations<n_vertices_per_facet, n_vertices>();
    std::unordered_map<std::array<int, n_vertices_per_facet>, int, std_array_hash<int, n_vertices_per_facet>> visited;
    std::array<int, n_vertices_per_facet> facet;
    // vinene riservata memoria per i container
	node_to_nodes.reserve(n_nodes);
	node_to_elems.reserve(n_elements);
	elem_to_facets.reserve(n_elements);
	for(unsigned id_node = 0; id_node < n_nodes; ++id_node)
	{
		node_to_elems.push_back({});
		node_to_nodes.push_back({});
		active_nodes.insert(id_node);
	}			
	// Loop over all elements
	for (int id_elem = 0; id_elem < mesh.n_elements(); ++id_elem)
	{
		active_elements.insert(id_elem);
		// Extract element
		auto & elem = mesh.element(id_elem);
		std::array<int, ct_nvertices(M)> node_ids = elem.node_ids();
		// loop su tutte le coppie di vertici
		for(unsigned i = 0; i < ct_nvertices(M); ++i)
		{
			// inserisco la connessione node_ids[i] - id_elem in node_to_elems
			node_to_elems[node_ids[i]].insert(id_elem);
			for(unsigned j = i+1; j < ct_nvertices(M); ++j)
			{
				// inserisco la connessione node_ids[i] - node_ids[j] in node_to_nodes
				node_to_nodes[node_ids[i]].insert(node_ids[j]);
				node_to_nodes[node_ids[j]].insert(node_ids[i]);
			}
		}
		
	}
	// DOMANDA:
	// vanno inserite anche le facets o rimangono nella classe mesh? Se rimangono nella classe mesh
	// posso prendere in questa classe le informazioni necessarie quando serve? Oppure non sono informazioni importanti?
}


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

}

std::unordered_set<unsigned> Connections::replace_node_in_node_to_nodes(unsigned old_id, unsigned new_id, 
	const std::unordered_set<unsigned>& involved)
{
	std::unordered_set<unsigned> old_node_to_nodes = node_to_nodes[new_id];
	active_nodes.erase(old_id);
	node_to_nodes[new_id] = involved;

	// For all connected, replace oldId with newId
	for (unsigned id : involved){
		node_to_nodes[id].erase(old_id);
		node_to_nodes[id].insert(new_id);
	}
		
	return old_node_to_nodes;
}


std::unordered_set<unsigned> Connections::replace_node_in_node_to_elems(unsigned old_id, unsigned new_id)
{
	// Extract node-element connections for newId
	std::unordered_set<unsigned> old_node_to_elems = node_to_elems[new_id];
	
	// Move the connected from oldId to newId
	node_to_elems[new_id] = node_to_elems[old_id]; // ATTENZIONE: NELLA VECCHIA LIB QUI SI USAVA .insert(node2elem[oldId])
	
	// Make oldId inactive
	active_nodes.erase(old_id);
	
	return old_node_to_elems;
}

template<typename Element>
void Connections::erase_elems_in_node_to_elems(const std::vector<Element> & to_remove)
{
	for (const Element & elem : to_remove)
	{
		std::array<unsigned, Element::n_vertices> node_ids = elem.node_ids();
		// Set element inactive in element-node connections
		active_elements.erase(elem.ID());
		
		// Erase element from node-element connections
		for (unsigned j = 0; j < Element::n_vertices; ++j)
			node_to_elems[node_ids[j]].erase(elem.ID());
	}
}



std::unordered_set<unsigned> Connections::nodes_involved_in_edge_collapse(const std::vector<unsigned> & facet) const
{
	std::unordered_set<unsigned> conn_nodes = {};
	for(unsigned node_id : facet)
		std::set_union(conn_nodes.begin(), conn_nodes.end(), 
			node_to_nodes[node_id].begin(), node_to_nodes[node_id].end(),
			std::inserter(conn_nodes, conn_nodes.begin()));
	for(auto node_id : facet)
		conn_nodes.erase(node_id);
	return conn_nodes;
}

std::unordered_set<unsigned> Connections::replace_node_in_node_to_nodes(unsigned old_id, unsigned new_id)
{
	// old_id diventa inattivo
	active_nodes.erase(old_id);
	// vengono presi i nodi che sono connessi al nodo da sostituire
	auto old_connections = node_to_nodes[old_id];
	// un loop su tutte le connessioni
	for(unsigned id_node : old_connections)
	{
		// new_id deve essere connesso ai nodi che sono connessi al nodo da sostituire
		node_to_nodes[new_id].insert(id_node);
		// i nodi collegati al nodo da sostituire devono ora essere collegati al nodo new_id
		node_to_nodes[id_node].erase(old_id);
		node_to_nodes[id_node].insert(new_id);
	}
	// viene tolta l'autoconnessione new_id-new_id
	node_to_nodes[new_id].erase(new_id);
	// DEVO VEDERE COSA RITORNARE QUI
	return {};
}

/*
template<typename Element>
void Connections::collapse_facet(const unsigned id_facet, const std::vector<Element> & to_remove)
{
	SetType facet_nodes = facets[id_facet];
	// suppongo di far collassare tutta la faccia sul primo indice in facet_nodes
	unsigned id_main_node = *facet_nodes.begin()
	// gli elementi da eliminare devono essere solo 2
	assert(to_remove.size() == 2);
	// vengono eliminati gli elementi nelle connessioni
	erase_elems_in_node_to_elems(to_remove);
	auto it =facet_nodes.begin();
	++it;
	for(; it != facet_nodes.end(); ++it)
	{
		// *it è l'id da rimpiazzare con id_main_node
		replace_node_in_node_to_nodes(*it, id_main_node);
		replace_node_in_node_to_elems(*it, id_main_node);
	}

}
*/




}
}
#endif // __CONNECTIONS_H__