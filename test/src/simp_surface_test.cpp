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



TEST(simplification_test, sphere_onlygeo)
{
    Mesh<2, 3> sphere_mesh = read_inp("../data/mesh/sfera.inp");
    GeomCost geom_cost;
    Simplification simp(sphere_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the sphere with only geometric cost\n"<<RESET;
    std::cout<<"Initial nodes: 5040, final nodes: 2500\n";
    simp.simplify(2500, geom_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

/*TEST(simplification_test, sphere_datageo)
{
    Mesh<2, 3> sphere_mesh = read_inp("../data/mesh/sfera.inp");
    GeomCost geom_cost;
    DataDispCost<2, 3> data_disp_cost;
    DataDistCost data_dist_cost;
    Simplification simp(sphere_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the sphere with geometric, data distance and data distribution costs\n"<<RESET;
    std::cout<<"Initial nodes: 5040, final nodes: 2500\n";
    std::array<double, 3> w = {1./3., 1./3., 1./3.};
    simp.simplify(2500, w,geom_cost, data_dist_cost, data_disp_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}*/

TEST(simplification_test, pawn_onlygeo)
{
    Mesh<2, 3> pawn_mesh = read_inp("../data/mesh/pawn.inp");
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
    Mesh<2, 3> pawn_mesh = read_inp("../data/mesh/toro.inp");
    GeomCost geom_cost;
    Simplification simp(pawn_mesh);
    std::cout<<BLUE<<"\nStarting simplification of the torus with only geometric cost\n"<<RESET;
    std::cout<<"Initial nodes: 7496, final nodes: 3500\n";
    simp.simplify(3500, geom_cost);
    std::cout<<BLUE<<"Simplification completed\n"<<RESET;
}

/*TEST(simplification_test, torus_sharp)
{
    Mesh<2, 3> pawn_mesh = read_inp("../data/mesh/toro.inp");
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
}*/