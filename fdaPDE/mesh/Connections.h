#ifndef __CONNECTIONS_H__
#define __CONNECTIONS_H__

#include "../utils/symbols.h"
#include <set>
#include <unordered_set>
#include <vector>
#include <algorithm>


/*
COMMENTI GENERALI:
Per ora una facet è un vwttore di nodi. Si può pensare di utilizzare un array, ma in questo caso 
devo prendere il numero di vertici per faccia dalla mesh.
Poi si può pensare di spostare il vettore di facce dalla classe mesh in questa classe. E allora tutti
i metodi che prendono in input una faccia rappresentata come vettore degli id dei suoi vertici possono essere
riscritti prendendo in input un solo unsigned (l'id della faccia).
*/

namespace fdapde{
namespace core{


class Connections{
	using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
	// using FacetType = std::vector<unsigned>;
	using SetType = std::unordered_set<unsigned>;
private:
	ConnectionsType node_to_nodes;
	ConnectionsType node_to_elems;
	// Da utilizzare per tenere conto delle facce che non possono essere più contratte
	// dopo aver utilizzato il metodo collapse_facet
	std::unordered_map<std::set<int>, unsigned, std_set_hash<int>> facets;

	// da eliminare e tenere dentro al metodo simplify?
	std::set<unsigned> active_nodes;
	std::set<unsigned> active_elements;

	// metodi privati
	// da rivedere
	std::unordered_set<unsigned> replace_node_in_node_to_nodes(unsigned , unsigned , const std::unordered_set<unsigned>& involved);
	// i nodi connessi a id_old vengono connessi invece a id_new
	void replace_node_in_node_to_nodes(unsigned old_id, unsigned new_id);
	// a new_id vengono aggiunti gli elementi che sono connessi a old_id
	// attenzione: si suppone che durante una contrazione, gli elementi da eliminare siano già stati eliminati
	std::unordered_set<unsigned> replace_node_in_node_to_elems(unsigned old_id, unsigned new_id);

	template<typename Element>
	void erase_elems_in_node_to_elems(const std::vector<Element> & to_remove); // HA SENSO??
	// il metodo aggiorna facets supponendo di contrarre l'input facet
	// nel std::pair di output vengono forniti: 
	// 1) gli id delle facce che vengono eliminate dalla contrazione
	// 2) gli id delle facce che vengono modificate dalla contrazione
	template<typename FacetType>
	std::pair<SetType, SetType> update_facets(const FacetType facet);

public:
	Connections() = default;
	template<typename Mesh>
	Connections(Mesh&& mesh);
	// refresh rinumera tutti i dati contenuti nelle strutture node_to_nodes, node_to_elems
	void refresh();

	// implementazione temporanea dei metodi che trovano connessioni di vario tipo. Questa implementazione 
	// suppone che le facce siano contenute nella classe mesh

	// intersezione dei nodi collegati ai nodi in facet. 
	template<typename FacetType>
	std::unordered_set<unsigned> nodes_on_facet(const FacetType & facet) const;
	// trova i nodi che sono legati ad almeno due vertici della faccia
	// template<typename FacetType> 
	// std::unordered_set<unsigned> nodes_attached_to_facet(const FacetType facet) const;
	// unione dei nodi collegati ai nodi in facet
	template<typename FacetType>
	std::unordered_set<unsigned> nodes_involved_in_collapse(const FacetType & facet) const; 
	template<typename FacetType>
	std::unordered_set<unsigned> elems_erased_in_collapse(const FacetType & facet) const;
	template<typename FacetType>
	std::unordered_set<unsigned> elems_modified_in_collapse(const FacetType & facet) const;
	// vengono trovati i due elementi che insistono sulla faccia
	// template<typename FacetType>
	// std::pair<unsigned, unsigned> elems_on_facet(const FacetType & facet) const;
	// vengono presi i nodi connessi a id_node su due livelli di profondità
	std::unordered_set<unsigned> extended_node_patch(unsigned id_node) const;
	// prende gli edge di cui calcolare di nuovo il costo
	// implementata solo la versione con M=2
	// NOTA: come implementare la versione M=3??
	std::set<unsigned> facets_to_update(unsigned node_id) const;

	template<typename Element>
	std::unordered_set<unsigned> element_patch(const Element & elem) const;

	// getters
	std::unordered_set<unsigned> get_node_to_nodes(unsigned id_node) const{return node_to_nodes[id_node];}
	std::unordered_set<unsigned> get_node_to_elems(unsigned id_node) const{return node_to_elems[id_node];}
	std::set<unsigned> get_active_nodes() const{return active_nodes;}
	std::set<unsigned> get_active_elements() const{return active_elements;}


	// implementazione temporanea
	template<typename Element, typename FacetType>
	std::pair<SetType, SetType> collapse_facet(const FacetType & facet, const std::vector<Element> & to_remove);

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
	constexpr unsigned n_vertices = Mesh::n_vertices;
	const auto & elements = mesh.elements();
	unsigned n_nodes = mesh.n_nodes();
	unsigned n_elements = mesh.n_elements();
	unsigned n_facets = mesh.n_facets();

    // vinene riservata memoria per i container
	node_to_nodes.reserve(n_nodes);
	node_to_elems.reserve(n_elements);
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
		const auto & elem = mesh.element(id_elem);
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

	// vengono aggiunte le facce con i loro id
	for(unsigned id_facet = 0; id_facet < n_facets; ++id_facet){
		// si può utilizzare l'iteratore definito per facets?
		auto node_ids = mesh.facet(id_facet).node_ids();
		facets[{node_ids.begin(), node_ids.end()}] = id_facet;
	}

	// DOMANDA:
	// vanno inserite anche le facets o rimangono nella classe mesh? Se rimangono nella classe mesh
	// posso prendere in questa classe le informazioni necessarie quando serve? Oppure non sono informazioni importanti?
}

template<typename FacetType>
std::unordered_set<unsigned> Connections::nodes_on_facet(const FacetType & facet) const
{
	// intersezione delle connessioni per ogni nodo in facet
	std::unordered_set<unsigned> conn_nodes = node_to_nodes[facet[0]];
	for(unsigned i = 1; i < facet.size(); ++i){
		for (auto it = conn_nodes.begin(); it != conn_nodes.end();) {
		    if (!node_to_nodes[facet[i]].count(*it)) { it = conn_nodes.erase(it); }
		    else              { ++it; }
		}
	}
	return conn_nodes;
}

template<typename FacetType>
std::unordered_set<unsigned> Connections::nodes_involved_in_collapse(const FacetType & facet) const
{
	std::unordered_set<unsigned> conn_nodes;
	for(unsigned id_node : facet){
		const auto & tmp_conn = node_to_nodes[id_node];
		conn_nodes.insert(tmp_conn.begin(), tmp_conn.end());
	}
	for(unsigned id_node : conn_nodes)
		conn_nodes.erase(id_node);
	return conn_nodes;
}

/*
template<typename FacetType>
std::pair<unsigned, unsigned> Connections::elems_on_facet(const FacetType & facet) const
{
	auto it_facet = facet.begin();
	auto conn_elems = node_to_elems[*it_facet];
	for(; it_facet!= facet.end(); ++it_facet)
		for (auto it = conn_elems.begin(); it != conn_elems.end();) {
		    if (!node_to_elems[*it_facet].count(*it)) { it = conn_elems.erase(it); }
		    else              { ++it; }
		}
	assert(conn_elems.size()==2);
}
*/

template<typename FacetType>
std::unordered_set<unsigned> Connections::elems_erased_in_collapse(const FacetType & facet) const
{/*
	// intersezione degli elementi connessi ai nodi di facet
	std::unordered_set<unsigned> conn_elems = node_to_elems[facet[0]];
	for(unsigned i = 1; i < facet.size(); ++i){
		for (auto it = conn_elems.begin(); it != conn_elems.end();) {
		    if (!node_to_elems[facet[i]].count(*it)) { it = conn_elems.erase(it); }
		    else              { ++it; }
		}
	}*/
	std::unordered_set<unsigned> conn_elems; // set di output
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
}

template<typename FacetType>
std::unordered_set<unsigned> Connections::elems_modified_in_collapse(const FacetType & facet) const
{
	auto to_erase = elems_erased_in_collapse(facet);
	// unione di tutti gli elementi collegati ai nodi di facet
	std::unordered_set<unsigned> conn_elems;
	for(unsigned id_node : facet){
		const auto & tmp_conn = node_to_elems[id_node];
		conn_elems.insert(tmp_conn.begin(), tmp_conn.end());
	}
	for(unsigned id_elem : to_erase)
		conn_elems.erase(id_elem);
	return conn_elems;
}


std::unordered_set<unsigned> Connections::extended_node_patch(unsigned id_node) const
{
	std::unordered_set<unsigned> conn_nodes = node_to_nodes[id_node];
	const auto & tmp_set1 = node_to_nodes[id_node];
	for(unsigned id : tmp_set1){
		const auto & tmp_set2 = node_to_nodes[id];
		conn_nodes.insert(tmp_set2.begin(), tmp_set2.end());
	}
	return conn_nodes;

}

template<typename Element>
std::unordered_set<unsigned> Connections::element_patch(const Element & elem) const
{
	auto node_ids = elem.node_ids();
	std::unordered_set<unsigned> conn_elems;
	for(unsigned id_node : node_ids){
		const auto & tmp_set = node_to_elems[id_node];
		conn_elems.insert(tmp_set.begin(), tmp_set.end());
	}
	conn_elems.erase(elem.ID());
	return conn_elems;
}

/*
template<typename FacetType> 
std::unordered_set<unsigned> Connections::nodes_attached_to_facet(const FacetType facet) const
{
	std::unordered_set<unsigned> conn_nodes;
	// sperando che il metodo size() di facet sia constexpr
	// auto facet_pattern = combinations<facet.size(), 2>();
	constexpr unsigned size = facet.size();
	auto facet_pattern = combinations<size, 2>();

	for(unsigned i = 0; i < facet_pattern.rows(); ++i)
	{
		std::unordered_set<unsigned> tmp_intersection;
		// vengono presi i set connessi a una coppia di nodi in facet
		const auto & node_set_1 = node_to_nodes[facet[facet_pattern(i, 0)]];
		const auto & node_set_2 = node_to_nodes[facet[facet_pattern(i, 1)]];
		// viene fatta l'intersezione dei due set
		for (unsigned id_node : node_set_1){
    		if (node_set_2.count(id_node)) { tmp_intersection.insert(id_node); }
		}
		// l'intersezione è posizionata in conn_nodes
		conn_nodes.insert(tmp_intersection.begin(), tmp_intersection.end());
	}
	
	return conn_nodes;
}
*/



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

void Connections::replace_node_in_node_to_nodes(unsigned old_id, unsigned new_id)
{
	active_nodes.erase(old_id);
	auto & old_conn = node_to_nodes[old_id];
	for(unsigned id_node : old_conn)
	{
		node_to_nodes[new_id].insert(id_node);
		node_to_nodes[id_node].erase(old_id);
		node_to_nodes[id_node].insert(new_id);
	}
	node_to_nodes[new_id].erase(new_id);
}


std::unordered_set<unsigned> Connections::replace_node_in_node_to_elems(unsigned old_id, unsigned new_id)
{
	// Extract node-element connections for newId
	auto & new_id_conn = node_to_elems[new_id];
	auto & old_id_conn = node_to_elems[old_id];
	// Move the connected from oldId to newId
	// node_to_elems[new_id] = node_to_elems[old_id]; // ATTENZIONE: NELLA VECCHIA LIB QUI SI USAVA .insert(node2elem[oldId])
	for(unsigned id_elem : old_id_conn)
		node_to_elems[new_id].insert(id_elem);
	// Make oldId inactive
	active_nodes.erase(old_id);
	return new_id_conn;
}

template<typename Element>
void Connections::erase_elems_in_node_to_elems(const std::vector<Element> & to_remove)
{
	assert(to_remove.size() == 2);
	for (const Element & elem : to_remove)
	{
		std::array<int, Element::n_vertices> node_ids = elem.node_ids();
		// Set element inactive in element-node connections
		active_elements.erase(elem.ID());
		
		// Erase element from node-element connections
		for (unsigned j = 0; j < Element::n_vertices; ++j)
			node_to_elems[node_ids[j]].erase(elem.ID());
	}
}


template<typename FacetType>
std::pair<std::unordered_set<unsigned>, std::unordered_set<unsigned>> Connections::update_facets(const FacetType facet)
{
	SetType to_erase;
	SetType to_modify;
	unsigned collapsing_node = facet[0];
	// caso di tetraedri
	if constexpr(facet.size() == 3)
	{
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
					tmp_facet.insert(conn_nodes[i]);
					if(facets.find(tmp_facet)!=facets.end()) // se i tre nodi in tmp_facet formano davvero una faccia
					{
						unsigned old_id = facets.at(tmp_facet);
						facets.erase(tmp_facet);
						tmp_facet.erase(facet[i]);
						tmp_facet.insert(collapsing_node);
						if(facets.find(tmp_facet)!=facets.end()) // la faccia già esiste, quindi old_id deve essere eliminato
							to_erase.insert(old_id);
						else{
							facets[tmp_facet] = old_id;
							to_modify.insert(old_id);
						}

					}
				}
			}

		}	
	}
	// caso di triangoli
	else
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
	}
	unsigned facet_id = facets.at({facet.begin(), facet.end()});
	// to_erase.insert(facet_id);
	facets.erase({facet.begin(), facet.end()});
	return {to_erase, to_modify};

} // update_facets

template<typename Element, typename FacetType>
std::pair<std::unordered_set<unsigned>, std::unordered_set<unsigned>> Connections::collapse_facet(const FacetType & facet, 
	const std::vector<Element> & to_remove)
{
	// la faccia viene contratta verso il primo nodo in facet
	assert(facets.find({facet.begin(), facet.end()})!=facets.end());
	unsigned collapsing_node = facet[0];
	auto facets_info = update_facets(facet);
	erase_elems_in_node_to_elems(to_remove);
	for(unsigned i = 1; i < facet.size(); ++i){
		replace_node_in_node_to_nodes(facet[i], collapsing_node);
		replace_node_in_node_to_elems(facet[i], collapsing_node);
	}
	for(unsigned i = 1; i < facet.size(); ++i)
		node_to_nodes[collapsing_node].erase(facet[i]);
	return facets_info;
}

std::set<unsigned> Connections::facets_to_update(unsigned node_id) const
{
	std::set<unsigned> facet_ids;
	const auto & node_ids = node_to_nodes[node_id];
	for(unsigned node_id2 : node_ids)
	{
		facet_ids.insert(facets.at({static_cast<int>(node_id), static_cast<int>(node_id2)}));
		const auto & node_ids2 = node_to_nodes[node_id2];
		for(unsigned node_id3 : node_ids2)
			facet_ids.insert(facets.at({static_cast<int>(node_id2), static_cast<int>(node_id3)}));
	}
	return facet_ids;
}

}
}
#endif // __CONNECTIONS_H__