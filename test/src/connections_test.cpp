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
