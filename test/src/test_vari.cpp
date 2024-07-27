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
        Simplification<2, 3> simp(mesh);
           
    }
    		 
}


TEST(simplification_test, sphere_onlygeo)
{
    Mesh<2, 3> sphere_mesh = read_inp("../../../meshes/sfera.inp");
    GeomCost geom_cost;
    Simplification simp(sphere_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the sphere with only geometric cost\n"<<RESET;
    std::cout<<"Initial nodes: 5040, final nodes: 2500\n";
    simp.simplify(2500, geom_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

TEST(simplification_test, sphere_datageo)
{
    Mesh<2, 3> sphere_mesh = read_inp("../../../meshes/sfera.inp");
    GeomCost geom_cost;
    DataDispCost<2, 3> data_disp_cost;
    DataDistCost data_dist_cost;
    Simplification simp(sphere_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the sphere with geometric, data distance and data distribution costs\n"<<RESET;
    std::cout<<"Initial nodes: 5040, final nodes: 2500\n";
    std::array<double, 3> w = {1./3., 1./3., 1./3.};
    simp.simplify(2500, w,geom_cost, data_dist_cost, data_disp_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

TEST(simplification_test, pawn_onlygeo)
{
    Mesh<2, 3> pawn_mesh = read_inp("../../../meshes/pawn.inp");
    GeomCost geom_cost;
    Simplification simp(pawn_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the pawn with only geometric cost\n"<<RESET;
    std::cout<<"Initial nodes: 2522, final nodes: 1500\n";
    simp.simplify(1500, geom_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

TEST(simplification_test, pawn_datageo)
{
    Mesh<2, 3> pawn_mesh = read_inp("../../../meshes/pawn.inp");
    GeomCost geom_cost;
    DataDispCost<2, 3> data_disp_cost;
    DataDistCost data_dist_cost;
    Simplification simp(pawn_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the pawn with geometric, data distance and data distribution costs\n"<<RESET;
    std::cout<<"Initial nodes: 2522, final nodes: 1500\n";
    std::array<double, 3> w = {1./3., 1./3., 1./3.};
    simp.simplify(1500, w,geom_cost, data_dist_cost, data_disp_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

TEST(simplification_test, torus_onlygeo)
{
    Mesh<2, 3> pawn_mesh = read_inp("../../../meshes/toro.inp");
    GeomCost geom_cost;
    Simplification simp(pawn_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the pawn with only geometric cost\n"<<RESET;
    std::cout<<"Initial nodes: 7496, final nodes: 3500\n";
    simp.simplify(3500, geom_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

TEST(simplification_test, torus_sharp)
{
    Mesh<2, 3> pawn_mesh = read_inp("../../../meshes/toro.inp");
    GeomCost geom_cost;
    DataDispCost<2, 3> data_disp_cost;
    DataDistCost data_dist_cost;
    SharpElemsCost<2, 3> sharp_elems_cost;
    Simplification simp(pawn_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the torus with geometric, data distance, data distribution costs and cost on sharp elements\n"<<RESET;
    std::cout<<"Initial nodes: 7496, final nodes: 3500\n";
    std::array<double, 4> w = {0.3, 0.3, 0.3, 0.1};
    simp.simplify(1500, w,geom_cost, data_dist_cost, data_disp_cost, sharp_elems_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

TEST(simplification_test, sphere_with_irregular_data)
{
    Simplification<2, 3> simp = read_surface_mesh("../../../meshes/simulation2_vertices.txt",
                                                  "../../../meshes/simulation2_triangles.txt",
                                                  "../../../meshes/simulation2_2500data.txt");
    GeomCost geom_cost;
    DataDispCost<2, 3> data_disp_cost;
    DataDistCost data_dist_cost;
    std::cout<<BLUE<<"\nStarting simplification of the sphere with data scattered irregularly with geometric, data distance and data distribution costs\n"<<RESET;
    std::cout<<"Initial nodes: 2522, final nodes: 1500\n";
    std::array<double, 3> w = {1./3., 1./3., 1./3.};
    simp.simplify(1500, w,geom_cost, data_dist_cost, data_disp_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}