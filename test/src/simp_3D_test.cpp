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


#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"




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


TEST(simplification_test, simplification_3D)
{
    Mesh<3, 3> mesh_3D = read_3D("../data/mesh/meshes_3D/elems_cubo.txt", 
                                 "../data/mesh/meshes_3D/nodes_cubo.txt",
                                 "../data/mesh/meshes_3D/boundary_cubo.txt");
    DataDispCost<3, 3> data_disp_cost;
    Simplification simp(mesh_3D);
    std::cout<<BLUE<<"\nStarting simplification of the 3D cube with data distribution costs\n"<<RESET;
    std::cout<<"Initial nodes: 2091, final nodes: 1500\n";
    simp.simplify(1500, data_disp_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}