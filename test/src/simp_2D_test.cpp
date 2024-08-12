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


TEST(simplification_test, simplification_2D)
{
    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::duration_cast;

    DataDispCost<2, 2> data_disp_cost;
    SharpElemsCost<2, 2> sharp_elems_cost;
    MeshLoader<Mesh2D> CShaped("unit_square_64");
    
    Simplification simp(CShaped.mesh);
    std::cout<<"simp inizializzata\n";
    unsigned n_nodes = CShaped.mesh.n_nodes();
    unsigned n_elements = CShaped.mesh.n_elements();
    std::cout<<"nodi mesh: "<<n_nodes<<", numero di elementi: "<<n_elements<<"\nInserire il numero di nodi\n";
    unsigned target_nodes;
    std::cin>>target_nodes;
    std::array<double, 2> w = {0.9, 0.1};
    // std::array<double, 1> w = {0.9};
    auto start = Clock::now();
    simp.simplify(target_nodes, w, data_disp_cost, sharp_elems_cost);
    // simp.simplify(target_nodes, data_disp_cost);
    auto end = Clock::now();
    auto elapsed = duration_cast<duration<double>>(end - start);
    std::cout<<"simplificazione finita. Tempo impiegato: "<<elapsed.count()<<"\n";
    auto mesh_simp = simp.build_mesh();
    std::ofstream file_nodes("../../../meshes/nodes_simp.txt");
    std::ofstream file_elems("../../../meshes/elems_simp.txt");
    std::ofstream file_data("../../../meshes/data_simp.txt");
    file_nodes<<mesh_simp.nodes();
    file_elems<<mesh_simp.elements();
    file_data<<simp.get_data();
    file_nodes.close();
    file_data.close();
    file_elems.close();



    std::cout<<"scrivere il nome del file della qoi\n";
    std::string nome_file;
    std::cin>>nome_file;

    // ================
    // SCRITTURA ANGOLI
    // ================
    std::ofstream file_angles("../../../meshes/angles_"+nome_file+".txt");
    for(unsigned id_elem = 0; id_elem < mesh_simp.n_elements(); ++id_elem)
    {
        SVector<2> A = mesh_simp.element(id_elem).coords()[0];
        SVector<2> B = mesh_simp.element(id_elem).coords()[1];
        SVector<2> C = mesh_simp.element(id_elem).coords()[2];
        double ABC = std::acos(  (A-B).dot(C-B)/( (B-A).norm()*(B-C).norm() )  );
        double BCA = std::acos(  (B-C).dot(A-C)/( (B-C).norm()*(A-C).norm() )  );
        double CAB = std::acos(  (C-A).dot(B-A)/( (C-A).norm()*(B-A).norm() )  );
        file_angles<<ABC<<std::endl<<BCA<<std::endl<<CAB<<std::endl;
    }
    file_angles.close();

    std::ofstream file_qoi("../../../meshes/dati/2D/"+nome_file+".txt");
    auto active_elems = simp.active_elems();
    for(auto elem_id : active_elems)
        file_qoi<<data_disp_cost.qoi_[elem_id]<<std::endl;
    file_qoi.close();
}