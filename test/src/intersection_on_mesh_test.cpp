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
TEST(ProvaLoadingMesh, Prova_1)
{
    MeshLoader<Mesh<2, 3>> meshloader("surface");

    StructuredGridSearch bb(meshloader.mesh);
    /*
    auto global_bb = bb.get_global_bounding_box();
    std::cout<<global_bb.first<<std::endl;
    cout<<global_bb.second<<endl;
    */

    unsigned el_id_test = 15;

    auto neighbouring_elements = bb.get_neighbouring_elements(meshloader.mesh.element(el_id_test));
    const auto & el = meshloader.mesh.element(el_id_test);
    
    bool do_intersect = false;
    for(unsigned el_id : neighbouring_elements)
        do_intersect = do_intersect || el.intersection(meshloader.mesh.element(el_id));

    const auto & nodes = meshloader.mesh.nodes();
    const auto & elems = meshloader.mesh.elements();
    unsigned n_nodes = meshloader.mesh.n_nodes();
    unsigned n_elements = meshloader.mesh.n_elements();

    std::string nodes_file_name = "../../../meshes/nodes.txt";
    std::string elems_file_name = "../../../meshes/elements.txt";

    // Apri il file di output
    std::ofstream nodes_file(nodes_file_name);
    std::ofstream elems_file(elems_file_name);

    if (nodes_file.is_open()) {
        // Scrivi i dati della matrice nel file
        nodes_file << nodes << endl;

        // Chiudi il file
        nodes_file.close();
    }
    else{
        std::cout<<"file non aperto\n";
    }

    if (elems_file.is_open()) {
        // Scrivi i dati della matrice nel file
        elems_file << elems << endl;

        // Chiudi il file
        elems_file.close();
    }
    else{
        std::cout<<"file non aperto\n";
    }


    EXPECT_FALSE(do_intersect);
}