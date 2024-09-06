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


Mesh<2, 3> read_inp(std::string file)
{
    std::ifstream mesh_file(file);
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
    DMatrix<int> boundary(n_nodes, 1);
    boundary.setZero();
    return Mesh<2, 3>(nodes, elements, boundary);

}



TEST(StructuredGridSearch_test, StructuredGridSearch_test_1)
{
    MeshLoader<SurfaceMesh> mesh_loader("surface");
    Mesh<2, 3> pawn_mesh = read_inp("../data/mesh/pawn.inp");
    // StructureGridSearch instantiation
    StructuredGridSearch sgs(pawn_mesh);

    {
        auto cell_size = sgs.get_cell_size();
        auto n_cells = sgs.get_n_cells();
        cout<<"cells along x axis: "<<n_cells[0]<<", dimension: "<<cell_size[0]<<endl;
        cout<<"cells along y axis: "<<n_cells[1]<<", dimension: "<<cell_size[2]<<endl;
        cout<<"cells along z axis: "<<n_cells[1]<<", dimension: "<<cell_size[2]<<endl;
    }

    {
        auto element_test = mesh_loader.mesh.element(10);
        auto neighbouring_elems = sgs.get_neighbouring_elements(element_test);
        cout<<"neighbouring elements:"<<endl;
        for(auto id : neighbouring_elems)
            cout<<id<<",  ";
        cout<<endl;
    }
}
