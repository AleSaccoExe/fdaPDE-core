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






TEST(simplification_test, simplification_3D)
{
    Mesh<3, 3> mesh_3D = read_3D("../data/mesh/meshes_3D/elems_cubo.txt", 
                                 "../data/mesh/meshes_3D/nodes_cubo.txt",
                                 "../data/mesh/meshes_3D/boundary_cubo.txt");
    DataEquiCost<3, 3> data_equi_cost;
    Simplification simp(mesh_3D);
    std::cout<<BLUE<<"\nStarting simplification of the 3D cube with data distribution costs\n"<<RESET;
    std::cout<<"Initial nodes: 2091, final nodes: 1500\n";
    simp.simplify(1500, data_equi_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}