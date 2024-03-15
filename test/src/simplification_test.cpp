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
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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



TEST(simplification_test, test_with_intersection_check)
{
    constexpr unsigned N = 3;
    constexpr unsigned M = 2;
    // carica la mesh
    MeshLoader<Mesh<2, 3>> meshloader("surface");
    auto & mesh = meshloader.mesh;
    StructuredGridSearch sgs(mesh);
    // crea le connessioni
    Connections connections(mesh);
    // si sceglie l'insieme delle facce da contrarre
    auto vec = generateUniqueRandomNumbers(0.2*mesh.n_facets(), mesh.n_facets());
    std::set<unsigned> facets_to_collapse(vec.begin(), vec.end());
    // viene creato il vettore degli elementi da modificare
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

    // comincia la semplificazione

    std::vector<Element<M, N>> elems_tmp;

    for(unsigned facet_to_collapse : facets_to_collapse)
    {
        std::cout<<"collapse di "<<facet_to_collapse<<"\n";
        elems_tmp.clear();
        // si controlla se la contrazione causa intersezioni
        auto facet = facets[facet_to_collapse];
        // prima si vedono quali elementi vengono eliminati e quali vengono modificati
        auto elems_to_modify = connections.elems_modified_in_collapse(facet);
        auto elems_to_erase = connections.elems_on_facet(facet);
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
        std::cout<<"elementi modificati\n";
        // ora passo a sgs il vettore di elementi elems_tmp, che è formato di elementi da modificare
        // e il set elems_to_erase. La classe provvede quindi a capire le possibili intersezioni date dal collapse di facet
        // Una cosa è da vedere però: al metodo get_neighbouring_elements quali elementi vanno passati?

        sgs.update(elems_tmp, elems_to_erase);
        std::cout<<"sgs updated\n";
        bool valid_collapse = true;
        for(auto & elem : elems_tmp){
            auto elems_to_check = sgs.get_neighbouring_elements(elem);
            for(unsigned elem_id : elems_to_check)
                if(elems_vec[elem_id].intersection(elem))
                {
                    std::cout<<"trovata intersezione\n";
                    valid_collapse = false;
                    break;
                }  
            if(!valid_collapse)
                break;
        }

        // ora vengono modificate i dati in base al risultato dei test di intersezione
        if(valid_collapse) // non ci sono state intersezioni: collapse valido
        {
            std::cout<<"collapse valido\n";
            for(auto & elem : elems_tmp){ // loop sugli elementi modificati dalla contrazione
                // viene sostituito nel vettore l'elemento modificato
                elems_vec[elem.ID()] = elem;
                // nella matrice degli elementi vengono modificati gli id dei nodi
                for(unsigned i = 0; i < N; ++i)
                    elems_mat(elem.ID(), i) = elem.node_ids()[i];
            }
            elems_tmp.clear(); // questo vettore viene ora riempito di elementi da eliminare
            for(unsigned elem_id : elems_to_erase) // loop sugli elementi da eliminare
                elems_tmp.push_back(elems_vec[elem_id]);
            // vengono modificate le connessioni e prese le informazioni necessarie a modificare
            // le facce
            auto facets_pair = connections.collapse_facet(facet, elems_tmp);
            std::cout<<"facet collapse fatto\n";
            for(unsigned facet_id : facets_pair.first) // loop sulle facce da eliminare
                facets_to_collapse.erase(facet_id); // vengono eliminate le facce che non esistono più
            for(unsigned facet_id : facets_pair.second) // loop sulle facce da modificare
            {
                auto & facet_to_modify = facets[facet_id];
                for(int i = 0; i < Mesh<M, N>::n_vertices_per_facet; ++i)
                    for(int j = 1; j < Mesh<M, N>::n_vertices_per_facet; ++j)
                        if(facet_to_modify[i] == facet[j]){  facet_to_modify[i] = facet[0]; }
            }
            std::cout<<"facce modificate\n";
        } // if(valid_collapse)
        else // il collapse non è valido perchè causa intersezioni
        {
            std::cout<<"collapsse non valido\n";
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

    // viene ricostruita ora la mesh
    auto active_elems = connections.get_active_elements();
    auto active_nodes = connections.get_active_nodes();
    DMatrix<double> new_nodes(active_nodes.size(), N);
    DMatrix<int> new_elems(active_elems.size(), Mesh<M, N>::n_vertices);
    DMatrix<int> new_boundary(active_nodes.size(), 1);
    new_boundary.setZero();

    std::unordered_map<unsigned, unsigned> node_ids_map;
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
    file_nodes<<new_nodes;
    file_elems<<new_elems;
    Mesh<M, N> mesh_simp(new_nodes, new_elems, new_boundary);

}