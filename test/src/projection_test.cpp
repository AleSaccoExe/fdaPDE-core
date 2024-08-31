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

TEST(projection_test, projection_test_1)
{
    // build a 2D square mesh: 
    MeshLoader<Mesh<2, 2>> meshloader("unit_square_16");
    const auto& mesh2D = meshloader.mesh;
    // add a third coordinate
    DMatrix<double> nodes(mesh2D.n_nodes(), 3);
    nodes.setZero();
    nodes.col(0) = meshloader.points_.col(0);
    nodes.col(1) = meshloader.points_.col(1);
    // build the coresponding 2.5D mesh
    Mesh<2, 3> mesh(nodes, meshloader.elements_, meshloader.boundary_);

    // TEST
    {
        std::vector<Element<2, 3>> elems;
        for(unsigned elem_id = 0; elem_id < mesh.n_elements(); ++elem_id) {elems.push_back(mesh.element(elem_id));}
        std::cout<<std::endl;
        // Point to project:
        SVector<3> A = {0.5, 0.5, 1.0};
        std::cout<<"\nOriginal datum position:\n";
        std::cout <<A;
        // build the vector of element of the mesh
        auto p_A = project(elems, A);
        std::cout<<"\nProjected datum:\n";
        std::cout<<p_A<<"\n";

        SVector<3> B = {1.5, 2.0, 3.0};
        std::cout<<"\nOriginal datum position:\n";
        std::cout <<B;
        // build the vector of element of the mesh
        auto p_B = project(elems, B);
        std::cout<<"\nProjected datum:\n";
        std::cout<<p_B<<"\n";
    }
}