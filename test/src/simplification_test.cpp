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

#include <set>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
using fdapde::core::Element;
using fdapde::core::Mesh;

#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MESH_TYPE_LIST;
using fdapde::testing::MeshLoader;
using namespace fdapde::core;
using namespace std;


/*
TEST(simplification_test, surface)
{
    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::duration_cast;

    GeomCost geom_cost;
    DataDistCost data_dist_cost;
    SharpElemsCost<2, 3> sharp_elems_cost;
    DataDispCost<2, 3> data_disp_cost;
    
    // MeshLoader<Mesh<2, 3>> meshloader("surface");
    std::ifstream mesh_file("../../../meshes/pawn.inp");
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
    }*/
    /*
    std::ifstream orig_elems_file("../../../meshes/simulation2_triangles.txt");
    std::ifstream orig_nodes_file("../../../meshes/simulation2_vertices.txt");
    std::ifstream orig_data_file("../../../meshes/simulation2_2500data.txt");
    unsigned n_elements = 0;
    unsigned n_nodes = 0;
    unsigned n_data = 0;
    std::string line;
    while(getline(orig_nodes_file, line)) {++n_nodes;}
    while(getline(orig_elems_file, line)) {++n_elements;}
    while(getline(orig_data_file, line)) {++n_data;}
    DMatrix<int> elements(n_elements, 3);
    DMatrix<double> nodes(n_nodes, 3);
    DMatrix<double> data(n_data, 3);
    orig_elems_file.clear();
    orig_nodes_file.clear();
    orig_data_file.clear();
    orig_nodes_file.seekg(0);
    orig_elems_file.seekg(0);
    orig_data_file.seekg(0);
    for (int i = 0; i < n_data; ++i) {
        for (int j = 0; j < 3; ++j) {
            orig_data_file >> data(i, j);
        }
    }
    for(unsigned i = 0; i<n_nodes; ++i)
    {
        getline(orig_nodes_file, line);
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
        getline(orig_elems_file, line);
        std::string useless;
        std::istringstream ss(line);
        ss>>useless;
        std::string node_id1, node_id2, node_id3;
        ss>>node_id1>>node_id2>>node_id3;
        elements(i, 0) = std::stoi(node_id1)-1;
        elements(i, 1) = std::stoi(node_id2)-1;
        elements(i, 2) = std::stoi(node_id3)-1;
    }*/
    /*
    // Simplification simp(meshloader.mesh);
    DMatrix<int> boundary(n_nodes, 1);
    boundary.setZero();
    Mesh<2, 3> mesh(nodes, elements, boundary);
    std::cout<<"mesh creata\n";

    Simplification simp(mesh);
    std::cout<<"simp inizializzata\n";
    // std::cout<<"nodi mesh: "<<meshloader.mesh.n_nodes()<<"\nInserire il numero di nodi\n";
    std::cout<<"nodi mesh: "<<n_nodes<<", numero di elementi: "<<n_elements<<"\nInserire il numero di nodi\n";
    unsigned target_nodes;
    std::cin>>target_nodes;
    std::array<double, 1> w = {0.5};
    // std::array<double, 3> w = {1./3., 1./3., 1./3.};
    auto start = Clock::now();
    // simp.simplify(target_nodes, w, geom_cost, data_disp_cost, data_dist_cost);
    simp.simplify(target_nodes, w, geom_cost);
    auto end = Clock::now();
    auto elapsed = duration_cast<duration<double>>(end - start);
    std::cout<<"simplificazione finita. Tempo impiegato: "<<elapsed.count()<<"\n";
    auto mesh_simp = simp.build_mesh();
    std::ofstream file_nodes("../../../meshes/nodes_simp.txt");
    std::ofstream file_elems("../../../meshes/elems_simp.txt");
    std::ofstream file_data("../../../meshes/data_simp.txt");
    file_nodes<<mesh_simp.nodes();
    file_elems<<mesh_simp.elements();
    file_data<<simp.get_data();
    file_nodes.close();
    file_data.close();
    file_elems.close();*/
    /*
    Connections conns(mesh);
    unsigned facet_id = 150;
    std::array<int, 2> facet = mesh.facet(facet_id).node_ids();
    auto elems_to_modify_ids = conns.elems_modified_in_collapse(facet);
    auto elems_modified = simp.modify_elements(elems_to_modify_ids, facet, 0.5*(mesh.node(facet[0]) + mesh.node(facet[1])) );
    std::vector<Element<2, 3>> elems_to_modify;
    for(unsigned elem_id : elems_to_modify_ids) {elems_to_modify.push_back(mesh.element(elem_id));}
    std::ofstream file_elems_before("../../../meshes/elems_before.txt");
    std::ofstream file_elems_after("../../../meshes/elems_after.txt");
    for(unsigned i = 0; i < elems_to_modify.size(); ++i)
    {
        auto coords_before = elems_to_modify[i].coords();
        auto coords_after = elems_modified[i].coords();
        for(unsigned j = 0; j < 3; ++j)
        {
            for(unsigned k = 0; k < 3; ++k){
                file_elems_after<<coords_after[j][k]<<" ";
                file_elems_before<<coords_before[j][k]<<" ";
            }
            file_elems_before<<std::endl;
            file_elems_after<<std::endl;
        }

    }
    file_elems_before.close();
    file_elems_after.close();
    */


    /*
    std::cout<<"scrivere il nome del file di qoi e dist\n";
    std::string nome_file;
    std::cin>>nome_file;

    std::ofstream file_dist("../../../meshes/dist_"+nome_file+".txt");
    // viene stampata la distanza dei dati
    for (int i = 0; i < mesh.n_nodes(); ++i) {
    // for (int i = 0; i < data.rows(); ++i) {
        SVector<3> new_pos = simp.get_data().row(i);
        SVector<3> orig_pos = mesh.node(i);
        // SVector<3> orig_pos = data.row(i);
        // Scriviamo il numero nel file di testo
        file_dist << (new_pos - orig_pos).norm() << std::endl;
    }
    std::cout<<"distanze scritte\n";

    // Chiudiamo il file dopo aver scritto tutti i numeri
    file_dist.close();
    
    std::ofstream file_qoi("../../../meshes/qoi_"+nome_file+".txt");
    auto active_elems = simp.active_elems();
    for(auto elem_id : active_elems)
        file_qoi<<data_disp_cost.qoi_[elem_id]<<std::endl;
    file_qoi.close();



    std::cout<<"check sull'intersezione tra elementi\n";
    unsigned intersezioni_trovate = 0;
    bool do_intersect = false;
    for(unsigned i = 0; i < mesh_simp.n_elements(); ++i){
        for(unsigned j = 0; j < mesh_simp.n_elements(); ++j)
            if(i!=j)
            {
                do_intersect = do_intersect || mesh_simp.element(i).intersection(mesh_simp.element(j));
                if(do_intersect)
                {
                    std::cout<<"el: "<<i<<std::endl;
                    auto node_ids1 = mesh_simp.element(i).node_ids();
                    std::cout<<node_ids1[0]<<" "<<node_ids1[1]<<" "<<node_ids1[2]<<"\n";
                    std::cout<<"el: "<<j<<std::endl;
                    auto node_ids2 = mesh_simp.element(j).node_ids();
                    std::cout<<node_ids2[0]<<" "<<node_ids2[1]<<" "<<node_ids2[2]<<"\n";
                    //break;
                    intersezioni_trovate++;
                }
            }
        // if(do_intersect)
            // break;
    }
    if(do_intersect)
        std::cout<<"trovata intersezione\n";
}*/

TEST(simplification_test, simplification_2D)
{
    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::duration_cast;

    DataDispCost<2, 2> data_disp_cost;
    MeshLoader<Mesh2D> CShaped("c_shaped");
    
    Simplification simp(CShaped.mesh);
    std::cout<<"simp inizializzata\n";
    unsigned n_nodes = CShaped.mesh.n_nodes();
    unsigned n_elements = CShaped.mesh.n_elements();
    std::cout<<"nodi mesh: "<<n_nodes<<", numero di elementi: "<<n_elements<<"\nInserire il numero di nodi\n";
    unsigned target_nodes;
    std::cin>>target_nodes;
    std::array<double, 1> w = {1.0};
    auto start = Clock::now();
    simp.simplify(target_nodes, w, data_disp_cost);
    auto end = Clock::now();
    auto elapsed = duration_cast<duration<double>>(end - start);
    std::cout<<"simplificazione finita. Tempo impiegato: "<<elapsed.count()<<"\n";
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

}