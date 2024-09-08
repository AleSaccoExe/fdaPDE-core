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
#include <string>
#include <chrono>
#include <fstream>

using fdapde::core::Element;
using fdapde::core::Mesh;

#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MESH_TYPE_LIST;
using fdapde::testing::MeshLoader;
using namespace fdapde::core;
using namespace std;


Mesh<3,3> read_3D(std::string file_elems, std::string file_nodes, std::string file_boundary)
{
    std::ifstream orig_elems_file(file_elems);
    std::ifstream orig_nodes_file(file_nodes);
    std::ifstream orig_boundary_file(file_boundary);
    unsigned n_elements = 0;
    unsigned n_nodes = 0;
    std::string line;
    while(getline(orig_nodes_file, line)) {++n_nodes;}
    while(getline(orig_elems_file, line)) {++n_elements;}
    DMatrix<int> elements(n_elements, 4);
    DMatrix<double> nodes(n_nodes, 3);
    DMatrix<int> boundary(n_nodes, 1);

    orig_elems_file.clear();
    orig_nodes_file.clear();
    orig_nodes_file.seekg(0);
    orig_elems_file.seekg(0);
    for(unsigned i = 0; i<n_nodes; ++i)
    {
        std::string boundary_line;
        getline(orig_nodes_file, line);
        getline(orig_boundary_file, boundary_line);
        std::string x, y, z;
        std::istringstream ss(line);
        ss>>x>>y>>z;
        nodes(i, 0) = std::stod(x);
        nodes(i, 1) = std::stod(y);
        nodes(i, 2) = std::stod(z);
        boundary(i, 0) = std::stoi(boundary_line);
    }
    for(unsigned i = 0; i< n_elements; ++i)
    {
        getline(orig_elems_file, line);
        std::istringstream ss(line);
        std::string node_id1, node_id2, node_id3, node_id4;
        ss>>node_id1>>node_id2>>node_id3>>node_id4;
        elements(i, 0) = std::stoi(node_id1)-1;
        elements(i, 1) = std::stoi(node_id2)-1;
        elements(i, 2) = std::stoi(node_id3)-1;
        elements(i, 3) = std::stoi(node_id4)-1;
    }
    Mesh<3, 3> mesh(nodes, elements, boundary);  
    return mesh;  
}


TEST(connections_test, connections_test_1)
{
    // build a 2D mesh
    MeshLoader<Mesh<2, 2>> meshloader("unit_square_16");
    const auto& mesh2D = meshloader.mesh;
    // add a third coordinate
    DMatrix<double> nodes(mesh2D.n_nodes(), 3);
    nodes.setZero();
    nodes.col(0) = meshloader.points_.col(0);
    nodes.col(1) = meshloader.points_.col(1);
    // build the test mesh
    Mesh<2, 3> mesh(nodes, meshloader.elements_, meshloader.boundary_);

    // TEST
    {
        Connections connections(mesh);
        const auto& node_to_nodes_0 = connections.get_node_to_nodes(18);
        std::cout<<"\nNodes connected to node 18:\n";
        for(auto node_id : node_to_nodes_0) {std::cout<<node_id<<"  ";}
        std::cout<<"\nNodes connected to node 19:\n";
        const auto& node_to_nodes_8 = connections.get_node_to_nodes(19);
        for(auto node_id : node_to_nodes_8) {std::cout<<node_id<<"  ";}
        const auto& node_to_elems_0 = connections.get_node_to_elems(18);
        std::cout<<"\nElements connected to node 18:\n";
        for(auto elem_id : node_to_elems_0) {std::cout<<elem_id<<"   ";}
        const auto& node_to_elems_8 = connections.get_node_to_elems(19);
        std::cout<<"\nElements connected to node 19:\n";
        for(auto elem_id : node_to_elems_8) {std::cout<<elem_id<<"   ";}

        // Define the edge vith vertices (18, 19)
        std::array<int, 2> facet({18, 19});
        const auto& elems_erased = connections.elems_erased_in_collapse(facet);
        const auto& elems_modified = connections.elems_modified_in_collapse(facet);
        std::cout<<"\nElems erased by the collapse of facet (18, 19):\n";
        for(auto elem_id : elems_erased) {std::cout<<elem_id<<"   ";}
        std::cout<<"\nElems modified by the collapse of facet (18, 19):\n";
        for(auto elem_id : elems_modified) {std::cout<<elem_id<<"   ";}   

        const auto& facets_to_update = connections.facets_to_update(38);
        std::cout<<"\nSuppoing that node 38 is the collapsing node, the facets to update are:\n";
        for(auto facet_id : facets_to_update) {std::cout<<facet_id<<"   ";}
    }            
}

TEST(connections_test, connectiont_test_3D)
{
    Mesh<3, 3> mesh_3D = read_3D("../data/mesh/meshes_3D/elems_cubo.txt", 
                                 "../data/mesh/meshes_3D/nodes_cubo.txt",
                                 "../data/mesh/meshes_3D/boundary_cubo.txt");
    Connections connections(mesh_3D);
    // TEST
    {
        unsigned id1=10, id2=11, id3=186;
        auto conn_nodes1 = connections.get_node_to_nodes(id1);
        auto conn_nodes2 = connections.get_node_to_nodes(id2);
        auto conn_nodes3 = connections.get_node_to_nodes(id3);
        cout<<"\nNodes connected to "<<id1<<":\n";
        for(unsigned id : conn_nodes1)
            cout<<id<<", ";
        cout<<"\nNodes connected to "<<id2<<":\n";
        for(unsigned id : conn_nodes2)
            cout<<id<<", ";
        cout<<"\nNodes connected to "<<id3<<":\n";
        for(unsigned id : conn_nodes3)
            cout<<id<<", ";
        cout<<endl;
        std::array<int, 3> facet({static_cast<int>(id1), static_cast<int>(id2), static_cast<int>(id3)});
        auto elems_to_erase = connections.elems_erased_in_collapse(facet);
        auto conn_elems1 = connections.get_node_to_elems(id1);
        auto conn_elems2 = connections.get_node_to_elems(id2);
        auto conn_elems3 = connections.get_node_to_elems(id3);
        cout<<"\nElements connected to "<<id1<<":\n";
        for(unsigned id : conn_elems1)
            cout<<id<<", ";
        cout<<"\nElements connected to "<<id2<<":\n";
        for(unsigned id : conn_elems2)
            cout<<id<<", ";
        cout<<"\nElements connected to "<<id3<<":\n";
        for(unsigned id : conn_elems3)
            cout<<id<<", ";
        cout<<"\nElements to erase when collapseing: "<<id1<<", "<<id2<<", "<<id3<<endl;
        for(unsigned id : elems_to_erase)
            cout<<id<<", ";
        cout<<endl;
    }
}  
