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
OK:
- get_facet

FORSE OK:
- facets_connected_to_node
- get_elem_to_facet


DA DEBUGGARE:
- nodes_on_facet
- elems_on_facet. DA ELIMINARE??
*/

TEST(connections_test, test_1)
{
    unsigned id_facet_test = 25;
    unsigned id_node_test = 30;
    unsigned id_elem_test = 15;
    using namespace std::chrono;
    MeshLoader<Mesh<2, 3>> meshloader("surface");
    high_resolution_clock::time_point start = high_resolution_clock::now();
    Connections connections(meshloader.mesh);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end-start).count();
    cout << "tempo per le connessioni: " << duration << " ms" << endl;
/*
    bool ok = true;

    auto facets = connections.facets_connected_to_node(id_node_test);
    for(unsigned id_facet : facets)
    {
        auto node_ids = connections.get_facet(id_facet);
        ok = ok && (node_ids.find(id_node_test) != node_ids.end());
    }
    */
    
    auto facet = connections.get_facet(id_facet_test);
    auto it = facet.begin();
    unsigned old_id = *it;
    ++it;
    unsigned new_id = *it;

    connections.replace_node_in_node_to_elems(old_id, new_id); 

    EXPECT_TRUE(true);
}