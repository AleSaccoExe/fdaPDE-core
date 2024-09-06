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



TEST(StructuredGridSearch_test, StructuredGridSearch_test_1)
{
    MeshLoader<SurfaceMesh> mesh_loader("surface");
    // StructureGridSearch instantiation
    StructuredGridSearch sgs(mesh_loader.mesh);

    {
        auto cell_size = sgs.get_cell_size();
        auto n_cells = sgs.get_n_cells();
        cout<<"cells along x axis: "<<n_cells[0]<<", dimension: "<<cell_size[0]<<endl;
        cout<<"cells along y axis: "<<n_cells[1]<<", dimension: "<<cell_size[1]<<endl;
        cout<<"cells along z axis: "<<n_cells[2]<<", dimension: "<<cell_size[2]<<endl;
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
