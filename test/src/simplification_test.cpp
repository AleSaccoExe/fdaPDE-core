// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <fdaPDE/mesh.h>
#include <fdaPDE/utils.h>
#include <gtest/gtest.h>   // testing framework

#include <memory>
#include <random>
#include <set>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <string>
using fdapde::core::Element;
using fdapde::core::Mesh;

#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MESH_TYPE_LIST;
using fdapde::testing::MeshLoader;
using namespace fdapde::core;
using namespace std;


std::vector<int> generateUniqueRandomNumbers(int N, int M) {
    std::vector<int> numbers;
    // Riempire un vettore con numeri da 0 a M-1
    for (int i = 0; i < M; ++i) {
        numbers.push_back(i);
    }
    // Mescolare il vettore per avere ordine casuale
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(numbers.begin(), numbers.end(), g);
    // Ridimensionare il vettore per ottenere solo i primi N numeri
    numbers.resize(N);
    return numbers;
}


/*
TEST(simplification_test, test_with_intersection_check)
{
    using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
    constexpr unsigned N = 3;
    constexpr unsigned M = 2;
    // carica la mesh
    MeshLoader<Mesh<2, 3>> meshloader("surface");
    auto & mesh = meshloader.mesh;
    StructuredGridSearch sgs(mesh);
    // crea le connessioni
    Connections connections(mesh);
    // si sceglie l'insieme delle facce da contrarre
    auto vec = generateUniqueRandomNumbers(0.1*mesh.n_facets(), mesh.n_facets());
    std::set<unsigned> facets_to_collapse(vec.begin(), vec.end());
    // viene creato il vettore degli elementi della mesh
    std::vector<Element<2, 3>> elems_vec;
    for(unsigned id_elem = 0; id_elem < mesh.n_elements(); ++id_elem)
        elems_vec.push_back(mesh.element(id_elem));
    // copiata la matrice dei nodi
    auto nodes = mesh.nodes();
    // per ora suppongo di avere i dati sui vertici della mesh
    auto data = mesh.nodes();
    // copiata la matrice degli elementi
    auto elems_mat = mesh.elements();
    // copiato il vettore di facce
    std::vector< std::array<int, Mesh<M, N>::n_vertices_per_facet> > facets;
    facets.reserve(mesh.n_facets());
    for(unsigned facet_id = 0; facet_id < mesh.n_facets(); ++facet_id)
        facets.push_back(mesh.facet(facet_id).node_ids());

    // vengono inizializzate le connessioni dati-elementi
    ConnectionsType data_to_elems(data.rows());
    ConnectionsType elem_to_data(elems_vec.size());
    for(unsigned i = 0; i < data.rows(); ++i) // loop su tutti i dati
    {
        auto tmp_set = connections.get_node_to_elems(i);
        data_to_elems[i] = tmp_set;
        for(unsigned elem_id : tmp_set) // loop sugli elementi connessi al dato i
            elem_to_data[elem_id].insert(i);
    }
    // ===========================
    // comincia la semplificazione
    // ===========================

    std::vector<Element<M, N>> elems_tmp;

    for(unsigned facet_to_collapse : facets_to_collapse)
    {
        auto facet = facets[facet_to_collapse];
        if(connections.nodes_on_facet(facet).size()==2 && !mesh.is_on_boundary(facet[1]))
        {

        // std::cout<<"collapse di "<<facet_to_collapse<<"\n";
        elems_tmp.clear();
        // si controlla se la contrazione causa intersezioni
        // prima si vedono quali elementi vengono eliminati e quali vengono modificati
        auto elems_to_modify = connections.elems_modified_in_collapse(facet);
        auto elems_to_erase = connections.elems_erased_in_collapse(facet);
        // ora si modificano gli id dei vertici deli elementi da modificare
        for(unsigned elem_to_modify : elems_to_modify)
        {
            auto elem_tmp = elems_vec[elem_to_modify];
            std::array<int, ct_nvertices(M)> node_ids;
            // il vettore degli id viene riempito
            for(unsigned i = 0; i < ct_nvertices(M); ++i)
                node_ids[i] = elems_mat(elem_to_modify, i);
            for(int i = 0; i < Mesh<M, N>::n_vertices; ++i)
                for(int j = 1; j < Mesh<M, N>::local_dimension; ++j)
                    if(facet[j] == node_ids[i]){
                        node_ids[i] = facet[0];
            }
            std::array<SVector<N>, ct_nvertices(M)> coords; 
            for(unsigned j = 0; j < ct_nvertices(M); ++j)
                coords[j] = nodes.row(node_ids[j]);
            // si costruisce ora un nuovo elemento con l'id del nodo modificato
            elems_tmp.emplace_back(elem_tmp.ID(), node_ids, coords, elem_tmp.neighbors(), elem_tmp.is_on_boundary());
        }
        // std::cout<<"elementi modificati\n";
        // ora passo a sgs il vettore di elementi elems_tmp, che è formato di elementi da modificare
        // e il set elems_to_erase. La classe provvede quindi a capire le possibili intersezioni date dal collapse di facet
        // Una cosa è da vedere però: al metodo get_neighbouring_elements quali elementi vanno passati?

        sgs.update(elems_tmp, elems_to_erase);
        // std::cout<<"sgs updated\n";
        bool valid_collapse = true;
        
        for(auto & elem : elems_tmp){
            auto elems_to_check = sgs.get_neighbouring_elements(elem);
            for(unsigned elem_id : elems_to_check)
                if(elems_vec[elem_id].intersection(elem))
                {
                    // std::cout<<"trovata intersezione\n";
                    valid_collapse = false;
                    break;
                }  
            if(!valid_collapse)
                break;
        }
        

        // ora vengono modificate i dati in base al risultato dei test di intersezione
        if(valid_collapse) // non ci sono state intersezioni: collapse valido
        {
            // ========
            // collapse
            // ========

            // vengono presi i dati da proiettare
            std::unordered_set<unsigned> data_ids;
            for(unsigned elem_id : elems_to_erase){
                data_ids.insert(elem_to_data[elem_id].begin(), elem_to_data[elem_id].end());
            }
            for(unsigned elem_id : elems_to_modify)
                data_ids.insert(elem_to_data[elem_id].begin(), elem_to_data[elem_id].end());
            // vengono proiettati i dati sugli elementi già modificati
            auto new_data_positions = project(elems_tmp, data, data_ids);
            // ora si può eseguire il collapse e dopo si cambiano le informazioni di data_to_elems e elem_to_data

            // std::cout<<"collapse valido\n";
            for(auto & elem : elems_tmp){ // loop sugli elementi modificati dalla contrazione
                // viene sostituito nel vettore l'elemento modificato
                elems_vec[elem.ID()] = elem;
                // nella matrice degli elementi vengono modificati gli id dei nodi
                for(unsigned i = 0; i < ct_nvertices(M); ++i)
                    elems_mat(elem.ID(), i) = elem.node_ids()[i];
            }
            elems_tmp.clear(); // questo vettore viene ora riempito di elementi da eliminare
            for(unsigned elem_id : elems_to_erase) // loop sugli elementi da eliminare
                elems_tmp.push_back(elems_vec[elem_id]);
            // vengono modificate le connessioni e prese le informazioni necessarie a modificare
            // le facce
            auto facets_pair = connections.collapse_facet(facet, elems_tmp);
            // std::cout<<"facet collapse fatto\n";
            assert(facets_pair.first.size()==2);
            for(unsigned facet_id : facets_pair.first) // loop sulle facce da eliminare
                facets_to_collapse.erase(facet_id); // vengono eliminate le facce che non esistono più
            for(unsigned facet_id : facets_pair.second) // loop sulle facce da modificare
            {
                auto & facet_to_modify = facets[facet_id];
                for(int i = 0; i < Mesh<M, N>::n_vertices_per_facet; ++i)
                    for(int j = 1; j < Mesh<M, N>::n_vertices_per_facet; ++j)
                        if(facet_to_modify[i] == facet[j]){  facet_to_modify[i] = facet[0]; }
            }
            // std::cout<<"facce modificate\n";

            // vengono ora modificate le connessioni tra dati e elementi
            unsigned i = 0;
            for(unsigned data_id : data_ids)
            {
                // prendo i vecchi id degli elementi connessi al dato data_id
                auto old_data_to_elems = data_to_elems[data_id];
                // rimuovo le vecchie connessioni
                for(unsigned elem_id : old_data_to_elems) // loop sugli elementi
                    elem_to_data[elem_id].erase(data_id);

                auto data_info = new_data_positions[i];
                if(data_info.first == -2) // -2 vuol dire che il dato è sull'elemento
                {
                    data_to_elems[data_id] = {static_cast<unsigned>(data_info.second)};
                    elem_to_data[data_info.second].insert(data_id);
                }
                else if(data_info.first == -1) // -1 vuol dire che il dato è su un nodo
                {
                    auto conn_elems = connections.get_node_to_elems(data_info.second); // elementi connessi al dato
                    data_to_elems[data_id] = conn_elems;
                    for(unsigned elem_id : conn_elems)
                        elem_to_data[elem_id].insert(data_id);
                }
                else // ultimo caso: il dato è dentro ad un edge
                {
                    auto conn_elems1 = connections.get_node_to_elems(data_info.first);
                    auto conn_elems2 = connections.get_node_to_elems(data_info.second);
                    // intersezione tra i due set
                    for (auto it = conn_elems1.begin(); it != conn_elems1.end();) {
                        if (!conn_elems2.count(*it)) { it = conn_elems1.erase(it); }
                        else              { ++it; }
                    }
                    // vengono ora aggiornati le connessioni
                    // assert(conn_elems1.size()==2); devono essere per forza 2?
                    data_to_elems[data_id] = conn_elems1;
                    for(unsigned elem_id : conn_elems1) {elem_to_data[elem_id].insert(data_id);}
                }
                ++i;
            }

        } // if(valid_collapse)
        else // il collapse non è valido perchè causa intersezioni
        {
            // ===================
            // collapse non valido
            // ===================

            // std::cout<<"collapsse non valido\n";
            // viene ripristinata la classe sgs. Devono essere aggiunti di nuovo gli elementi che erano stati 
            // eliminati prima del controllo delle intersezioni, e devono essere riportati alla situazione originale 
            // gli elementi che erano stati modificati
            elems_tmp.clear(); // questo vettore viene ora riempito di elementi da aggiungere
            for(unsigned elem_id : elems_to_erase) // loop sugli elementi che erano stati eliminati
                elems_tmp.push_back(elems_vec[elem_id]);
            sgs.add_elements(elems_vec);
            elems_tmp.clear(); // ora viene riempito di elementi da modificare
            for(unsigned elem_id : elems_to_modify)
                elems_tmp.push_back(elems_vec[elem_id]);
            sgs.update_f(elems_tmp);
        } // else

        }
        else {std::cout<<"nodi sulla faccia: "<<connections.nodes_on_facet(facet).size()<<"\n";}
    }

    // =============================
    // viene ricostruita ora la mesh
    // =============================
    auto active_elems = connections.get_active_elements();
    auto active_nodes = connections.get_active_nodes();
    DMatrix<double> new_nodes(active_nodes.size(), N);
    DMatrix<int> new_elems(active_elems.size(), Mesh<M, N>::n_vertices);
    DMatrix<int> new_boundary(active_nodes.size(), 1);
    new_boundary.setZero();

    std::unordered_map<unsigned, unsigned> node_ids_map; // da utilizzare per mappare i vecchi id dei nodi ai nuovi id
    unsigned new_id = 0;
    for(unsigned old_id : active_nodes){
        for(unsigned j = 0; j < N; ++j)
            new_nodes(new_id, j) = nodes(old_id, j); 
        node_ids_map[old_id] = new_id;
        new_id++;
    }
    std::cout<<"matrice dei nodi creata\n";
    new_id = 0;
    for(unsigned old_id : active_elems){
        for(unsigned j = 0; j < Mesh<M, N>::n_vertices; ++j)
            new_elems(new_id, j) = node_ids_map.at(elems_mat(old_id, j));
        new_id++;
    }
    std::cout<<"matrice degli elementi creata\n";
    std::ofstream file_nodes("../../../meshes/nodes_simp.txt");
    std::ofstream file_elems("../../../meshes/elems_simp.txt");
    std::ofstream file_data("../../../meshes/data_simp.txt");
    file_nodes<<new_nodes;
    file_elems<<new_elems;
    file_data<<data;
    Mesh<M, N> mesh_simp(new_nodes, new_elems, new_boundary);

    // ==============================
    // test finali sulle intersezioni
    // ==============================
    std::cout<<"inizio test sulle intersezioni\n";
    bool found_intersection = false;
    for(unsigned i = 0; i < mesh_simp.n_elements(); ++i){
        for(unsigned j = 0; j < mesh_simp.n_elements(); ++j){
            if(i!=j)
                found_intersection = found_intersection || 
                                mesh_simp.element(i).intersection(mesh_simp.element(j));
            if(found_intersection)
                break;
        }
        if(found_intersection)
            break;
    }
    if(found_intersection)
        std::cout<<"intersezione trovata!!\n";
    else
        std::cout<<"elementi: "<<mesh_simp.n_elements()<<". Nessuna intersezione trovata\n";

}


TEST(simplification_test, only_geo)
{
    using ConnectionsType = std::vector<std::unordered_set<unsigned>>;
    constexpr unsigned N = 3;
    constexpr unsigned M = 2;
    // carica la mesh
    MeshLoader<Mesh<2, 3>> meshloader("surface");
    auto & mesh = meshloader.mesh;
    StructuredGridSearch sgs(mesh);
    // crea le connessioni
    Connections connections(mesh);
    // viene creato il vettore degli elementi della mesh
    std::vector<Element<2, 3>> elems_vec;
    for(unsigned id_elem = 0; id_elem < mesh.n_elements(); ++id_elem)
        elems_vec.push_back(mesh.element(id_elem));
    // copiata la matrice dei nodi
    auto nodes = mesh.nodes();
    // per ora suppongo di avere i dati sui vertici della mesh
    auto data = mesh.nodes();
    // copiata la matrice degli elementi
    auto elems_mat = mesh.elements();
    // copiato il vettore di facce
    std::vector< std::array<int, Mesh<M, N>::n_vertices_per_facet> > facets;
    facets.reserve(mesh.n_facets());
    for(unsigned facet_id = 0; facet_id < mesh.n_facets(); ++facet_id)
        facets.push_back(mesh.facet(facet_id).node_ids());

    // vengono inizializzate le connessioni dati-elementi
    ConnectionsType data_to_elems(data.rows());
    ConnectionsType elem_to_data(elems_vec.size());
    for(unsigned i = 0; i < data.rows(); ++i) // loop su tutti i dati
    {
        auto tmp_set = connections.get_node_to_elems(i);
        data_to_elems[i] = tmp_set;
        for(unsigned elem_id : tmp_set) // loop sugli elementi connessi al dato i
            elem_to_data[elem_id].insert(i);
    }

    // mappa per prendere la faccia da contrarre (la faccia con costo minore)
    // - key: costo della contrazione
    // - value: id della faccia da contrarre + punto di contrazione
    std::multimap<double, std::pair<unsigned, SVector<N> >> costs_map;
    // mappa utile per l'update di collapse_costs dopo un collapse
    std::unordered_map<unsigned, double> facets_cost;

    // ===============
    // INIZIO SIMPLIFY
    // ===============

    GeomCost geom_cost;

    // ===========================
    // calcolo costi per ogni lato
    // ===========================
    std::cout<<"inizializzazione: calcolo dei costi\n";
    unsigned facet_id = 0;
    for(const auto & facet : facets){
        if(connections.nodes_on_facet(facet).size()==2)
        {

        // SVector<N> collapse_points = get_collapse_points(facet, nodes);
        std::vector<SVector<N>> collapse_points = {nodes.row(facet[0])};
        // per ogni punto in collapse_points si calcola il costo della contrazione
        std::map<double, SVector<N>> tmp_costs_map;
        auto elems_to_modify = connections.elems_modified_in_collapse(facet);
        auto elems_to_erase = connections.elems_erased_in_collapse(facet);
        for(auto collapse_point : collapse_points)
        {
            // viene creato un vettore degli elementi involved nel collapse di facet
            std::vector<Element<M, N>> elems_tmp;
            for(unsigned elem_id : elems_to_modify)
                elems_tmp.push_back(elems_vec[elem_id]);
            for(unsigned elem_id : elems_to_erase)
                elems_tmp.push_back(elems_vec[elem_id]);
            // calcolato il costo del collapse
            tmp_costs_map[geom_cost(elems_tmp, collapse_point)] = collapse_point;
        }
        // ora si prende il costo minore e si controllano le intersezioni
        for(auto it = tmp_costs_map.begin(); it != tmp_costs_map.end(); ++it)
        {
            // viene creato un vettore degli elementi modificati nel collapse di facet
            std::vector<Element<M, N>> elems_tmp;
            for(unsigned elem_to_modify : elems_to_modify)
            {
                auto elem_tmp = elems_vec[elem_to_modify];
                std::array<int, ct_nvertices(M)> node_ids;
                std::array<SVector<N>, ct_nvertices(M)> coords; 
                // il vettore degli id viene riempito
                for(unsigned i = 0; i < ct_nvertices(M); ++i){
                    node_ids[i] = elems_mat(elem_to_modify, i);
                    coords[i] = nodes.row(node_ids[i]);
                }
                for(int i = 0; i < Mesh<M, N>::n_vertices; ++i)
                    for(int j = 1; j < Mesh<M, N>::local_dimension; ++j)
                        if(facet[j] == node_ids[i]){
                            node_ids[i] = facet[0];
                            coords[i] = it->second;
                }
                // si costruisce ora un nuovo elemento con l'id del nodo modificato
                elems_tmp.emplace_back(elem_tmp.ID(), node_ids, coords, elem_tmp.neighbors(), elem_tmp.is_on_boundary());
            }
            // ora passo a sgs il vettore di elementi elems_tmp, che è formato di elementi da modificare
            // e il set elems_to_erase. La classe provvede quindi a capire le possibili intersezioni date dal collapse di facet
            sgs.update(elems_tmp, elems_to_erase);
            bool valid_collapse = true;
            for(auto & elem : elems_tmp){
                auto elems_to_check = sgs.get_neighbouring_elements(elem);
                for(unsigned elem_id : elems_to_check)
                    if(elems_vec[elem_id].intersection(elem))
                    {
                        // std::cout<<"trovata intersezione\n";
                        valid_collapse = false;
                        break;
                    }  
                if(!valid_collapse)
                    break;
            }
            // se non ci sono intersezioni si procede ad aggiungere le informazioni nelle strutture adeguate
            if(valid_collapse)
            {
                costs_map.insert({it->first, {facet_id, it->second} });
                facets_cost[facet_id] = it->first;
                break;
            }

        }
        std::vector<Element<M, N>> elems_tmp; // questo vettore viene ora riempito di elementi da aggiungere
        // a questo punto si ripristina la classe sgs al suo stato iniziale
        for(unsigned elem_id : elems_to_erase) // loop sugli elementi che erano stati eliminati
            elems_tmp.push_back(elems_vec[elem_id]);
        sgs.add_elements(elems_vec);
        elems_tmp.clear(); // ora viene riempito di elementi da modificare
        for(unsigned elem_id : elems_to_modify)
            elems_tmp.push_back(elems_vec[elem_id]);
        sgs.update_f(elems_tmp);

        }
        facet_id++;
    }

    for(auto pair : costs_map)
        std::cout<<pair.first<<", "<<pair.second.first<<"\n";
    for(auto pair : facets_cost)
        std::cout<<pair.first<<", "<<pair.second<<"\n";
}
*/ 

TEST(simplification_test, only_geo)
{
    GeomCost geom_cost;
    DataDistCost data_dist_cost;
    SharpElemsCost<2, 3> sharp_elems_cost;

    // MeshLoader<Mesh<2, 3>> meshloader("surface");
    std::ifstream mesh_file("../../../meshes/sfera.inp");
    int n_nodes, n_elements;
    std::string line;
    getline(mesh_file, line);
    std::string s_num_elem, s_num_ver, num_data;
    std::istringstream ss(line);
    ss>>n_nodes;
    ss>>n_elements;
    DMatrix<double> nodes(n_nodes, 3);
    DMatrix<int> elements(n_elements, 3);
    for(unsigned i = 0; i<n_nodes; ++i)
    {
        getline(mesh_file, line);
        std::string x, y, z;
        std::istringstream ss(line);
        std::string useless;
        ss>>useless;
        ss>>x>>y>>z;
        nodes(i, 0) = std::stod(x);
        nodes(i, 1) = std::stod(y);
        nodes(i, 2) = std::stod(z);
    }
    for(unsigned i = 0; i< n_elements; ++i)
    {
        getline(mesh_file, line);
        std::string useless;
        std::istringstream ss(line);
        ss>>useless; ss>>useless; ss>>useless;
        std::string node_id1, node_id2, node_id3;
        ss>>node_id1>>node_id2>>node_id3;
        elements(i, 0) = std::stoi(node_id1)-1;
        elements(i, 1) = std::stoi(node_id2)-1;
        elements(i, 2) = std::stoi(node_id3)-1;
    }
    // Simplification simp(meshloader.mesh);
    MeshLoader<Mesh<2, 2>> meshloader("unit_square_128");
    DMatrix<int> boundary(n_nodes, 1);
    boundary.setZero();
    Mesh<2, 3> mesh(nodes, elements, boundary);
    std::cout<<"mesh creata\n";
    Simplification simp(mesh);
    std::cout<<"simp inizializzata\n";
    // std::cout<<"nodi mesh: "<<meshloader.mesh.n_nodes()<<"\nInserire il numero di nodi\n";
    std::cout<<"nodi mesh: "<<n_nodes<<"\nInserire il numero di nodi\n";
    unsigned target_nodes;
    std::cin>>target_nodes;
    std::array<double, 2> w = {0.5, 0.5};
    simp.simplify(target_nodes, w, geom_cost, data_dist_cost);
    std::cout<<"simplificazione finita\n";
    auto mesh_simp = simp.build_mesh();
    std::ofstream file_nodes("../../../meshes/nodes_simp.txt");
    std::ofstream file_elems("../../../meshes/elems_simp.txt");
    std::ofstream file_data("../../../meshes/data_simp.txt");
    file_nodes<<mesh_simp.nodes();
    file_elems<<mesh_simp.elements();
    file_data<<simp.get_data();
    file_nodes.close();
    file_data.close();
    file_elems.close();
    std::cout<<"check sull'intersezione tra elementi\n";
    bool do_intersect = false;
    for(unsigned i = 0; i < mesh.n_elements(); ++i){
        for(unsigned j = 0; j < mesh.n_elements(); ++j)
            if(i!=j)
            {
                do_intersect = do_intersect || mesh.element(i).intersection(mesh.element(j));
                if(do_intersect)
                {
                    std::cout<<"el: "<<i<<std::endl;
                    auto node_ids1 = mesh.element(i).node_ids();
                    std::cout<<node_ids1[0]<<" "<<node_ids1[1]<<" "<<node_ids1[2]<<"\n";
                    std::cout<<"el: "<<j<<std::endl;
                    auto node_ids2 = mesh.element(j).node_ids();
                    std::cout<<node_ids2[0]<<" "<<node_ids2[1]<<" "<<node_ids2[2]<<"\n";
                    break;
                }
            }
        if(do_intersect)
            break;
    }
    if(do_intersect)
        std::cout<<"trovata intersezione\n";
    
}