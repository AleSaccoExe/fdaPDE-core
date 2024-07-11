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

TEST(test_vari, test_1)
{
    // mesh quadrata:
    /*
    nodes << 0, 0, 0, 
    		 1, 0, 0,
    		 1, 1, 0,
    		 0, 1, 0,
    		 0.5, 0, 0,
             0, 0.5, 0,
             0.5, 1, 0,
             0, 0.5, 0,
             0.5, 0.5, 0,
             0.25, 0, 0,
             0.5, 0.25, 0,
             0.25, 0.5, 0,
             0, 0.25, 0,
             0.25, 0.25, 0,
             0.75, 0, 0,
             1, 0.25, 0
             0.75, 0.5, 0,
             0.75, 0.25, 0,
             1, 0.75, 0,
             0.75, 1, 0,
             0.5, 0.75, 0,
             0.75, 0.75, 0,
             0.25, 1, 0,
             0, 0.75, 0,
             0.25, 0.75, 0;
    */
    // Costruzione di una mesh quadrata 2.5D

    // prima costruisco una mesh 2D quadrata: 
    MeshLoader<Mesh<2, 2>> meshloader("unit_square_16");
    const auto& mesh2D = meshloader.mesh;
    // ora i nodi diventano tridimensionali ponendo 0 la coordinata z
    DMatrix<double> nodes(mesh2D.n_nodes(), 3);
    nodes.setZero();
    nodes.col(0) = meshloader.points_.col(0);
    nodes.col(1) = meshloader.points_.col(1);

    Mesh<2, 3> mesh(nodes, meshloader.elements_, meshloader.boundary_);

    // test sulle connessioni
    {
        Connections connections(mesh);
        const auto& node_to_nodes_0 = connections.get_node_to_nodes(0);
        std::cout<<"\nNodes connected to node 0:\n";
        for(auto node_id : node_to_nodes_0) {std::cout<<node_id<<"  ";}
        std::cout<<"\nNodes connected to node 8:\n";
        const auto& node_to_nodes_8 = connections.get_node_to_nodes(8);
        for(auto node_id : node_to_nodes_8) {std::cout<<node_id<<"  ";}
        const auto& node_to_elems_0 = connections.get_node_to_elems(0);
        std::cout<<"\nElements connected to node 0:\n";
        for(auto elem_id : node_to_elems_0) {std::cout<<elem_id<<"   ";}
        const auto& node_to_elems_8 = connections.get_node_to_elems(8);
        std::cout<<"\nElements connected to node 8:\n";
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
    
    // test sulle proiezioni
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
    // test su simplification
    {
        
    }


    		 
}