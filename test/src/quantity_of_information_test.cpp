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


TEST(qoi_test, qoi_test_1)
{	

	MeshLoader<SurfaceMesh> mesh_loader("surface");
	DataEquiCost<2, 3> data_equi_cost;
	Simplification simp(mesh_loader.mesh);

	data_equi_cost.setup(&simp);

	unsigned elem_id = 30;
	double qoi = data_equi_cost.get_qoi(elem_id);
	cout<<"qoi of element "<< elem_id<<": "<<qoi<<endl;

	unsigned elem_id2 = 40;
	double qoi2 = data_equi_cost.get_qoi(elem_id2);
	cout<<"qoi of element "<< elem_id2<<": "<<qoi2<<endl; 
}