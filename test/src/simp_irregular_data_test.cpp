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





Simplification<2, 3> read_surface_mesh(std::string nodes_file, std::string elems_file, std::string data_file)
{
    std::ifstream orig_elems_file(elems_file);
    std::ifstream orig_nodes_file(nodes_file);
    std::ifstream orig_data_file(data_file);
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
    }
    DMatrix<int> boundary(n_nodes, 1);
    boundary.setZero();
    return Simplification(Mesh<2, 3>(nodes, elements, boundary), data);
}


TEST(simplification_test, sphere_with_irregular_data)
{
    Simplification<2, 3> simp = read_surface_mesh("../data/mesh/simulation2_vertices.txt",
                                                  "../data/mesh/simulation2_triangles.txt",
                                                  "../data/mesh/simulation2_2500data.txt");
    GeomCost geom_cost;
    DataDispCost<2, 3> data_disp_cost;
    DataDistCost data_dist_cost;
    std::cout<<BLUE<<"\nStarting simplification of the sphere with data scattered irregularly with geometric, data distance and data distribution costs\n"<<RESET;
    std::cout<<"Initial nodes: 3097, final nodes: 2500\n";
    std::array<double, 3> w = {1./3., 1./3., 1./3.};
    simp.simplify(1500, w,geom_cost, data_dist_cost, data_disp_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}