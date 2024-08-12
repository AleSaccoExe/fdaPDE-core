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
    std::ifstream mesh_file("../../../meshsimplification/mesh/pawn.inp");
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
    DMatrix<int> boundary(n_nodes, 1);
    boundary.setZero();
    Mesh<2, 3> mesh(nodes, elements, boundary);
    std::cout<<"mesh creata\n";
    StructuredGridSearch sgs(mesh);
    auto start_sgs = Clock::now();
    {
    auto ids = sgs.get_neighbouring_elements(mesh.element(10));
    for (auto id : ids)
            std::cout << id << " ";
    std::cout << std::endl;
    }
    
    /*{
    auto ids = sgs.get_neighbouring_elements(mesh.element(20));
    for (auto id : ids)
            std::cout << id << " ";
    std::cout << std::endl;
    }*/
    auto end_sgs = Clock::now();
    std::cout<<"tempo get_neighbouring_elements: "<<duration_cast<duration<double>>(end_sgs - start_sgs).count()<<"\n";
    Simplification simp(mesh);
    std::cout<<"simp inizializzata\n";
    // std::cout<<"nodi mesh: "<<meshloader.mesh.n_nodes()<<"\nInserire il numero di nodi\n";
    std::cout<<"nodi mesh: "<<n_nodes<<", numero di elementi: "<<n_elements<<"\nInserire il numero di nodi\n";
    unsigned target_nodes;
    std::cin>>target_nodes;
    // std::array<double, 3> w = {1./3., 1./3., 1./3.};
    // std::array<double, 4> w = {0.3, 0.3, 0.3, 0.1};
    // std::array<double, 3> w = {1./3., 0., 0.};
    auto start = Clock::now();
    // simp.simplify(target_nodes, w, geom_cost, data_disp_cost, data_dist_cost, sharp_elems_cost);
    simp.set_check_intersections(true);
    simp.simplify(target_nodes, geom_cost);
    auto end = Clock::now();
    auto elapsed = duration_cast<duration<double>>(end - start);
    std::cout<<"simplificazione finita. Tempo impiegato: "<<elapsed.count()<<"\n";
    auto mesh_simp = simp.build_mesh();
    StructuredGridSearch sgs_simp(mesh_simp);

}
