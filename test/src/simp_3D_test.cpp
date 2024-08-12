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


TEST(simplification_test, simplification_3D)
{
    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::duration_cast;

    DataDispCost<3, 3> data_disp_cost;
    SharpElemsCost<3, 3> sharp_elems_cost;
    // estrazione della mesh
    std::string file_name = "cubo";
    // vengono estratti i file: nodes_file_name.txt, elems_file_name.txt, boundary_file_name.txt
    std::ifstream orig_elems_file("../../../meshes/meshes_3D/elems_"+file_name+".txt");
    std::ifstream orig_nodes_file("../../../meshes/meshes_3D/nodes_"+file_name+".txt");
    std::ifstream orig_boundary_file("../../../meshes/meshes_3D/boundary_"+file_name+".txt");
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
    Simplification simp(mesh);
    std::cout<<"simp inizializzata\n";
    std::cout<<"nodi mesh: "<<n_nodes<<", numero di elementi: "<<n_elements<<"\nInserire il numero di nodi\n";
    unsigned target_nodes;
    std::cin>>target_nodes;
    std::array<double, 2> w = {0.9, 0.1};
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

    std::cout<<"inserire il nome dei file\n";
    std::string nome_file;
    std::cin>>nome_file;

    std::ofstream file_qoi("../../../meshes/dati/3D/"+nome_file+".txt");
    auto active_elems = simp.active_elems();
    for(auto elem_id : active_elems)
        file_qoi<<data_disp_cost.qoi_[elem_id]<<std::endl;
    file_qoi.close();

    std::ofstream file_angles("../../../meshes/dati/3D/angles_"+nome_file+".txt");
    for(unsigned id_elem = 0; id_elem < mesh_simp.n_elements(); ++id_elem)
    {
        SVector<3> A = mesh_simp.element(id_elem).coords()[0];
        SVector<3> B = mesh_simp.element(id_elem).coords()[1];
        SVector<3> C = mesh_simp.element(id_elem).coords()[2];
        SVector<3> D = mesh_simp.element(id_elem).coords()[3];
        // triangolo ABC
        double ABC = std::acos(  (A-B).dot(C-B)/( (B-A).norm()*(B-C).norm() )  );
        double BCA = std::acos(  (B-C).dot(A-C)/( (B-C).norm()*(A-C).norm() )  );
        double CAB = std::acos(  (C-A).dot(B-A)/( (C-A).norm()*(B-A).norm() )  );
        file_angles<<ABC<<std::endl<<BCA<<std::endl<<CAB<<std::endl;
        // triangolo ABD
        double ABD = std::acos(  (A-B).dot(D-B)/( (A-B).norm()*(D-B).norm() )  );
        double BDA = std::acos(  (B-D).dot(A-D)/( (B-D).norm()*(A-D).norm() )  );
        double DAB = std::acos(  (D-A).dot(B-A)/( (D-A).norm()*(B-A).norm() )  );
        file_angles<<ABD<<std::endl<<BDA<<std::endl<<DAB<<std::endl;
        // triangolo BDC
        double BDC = std::acos(  (B-D).dot(C-D)/( (B-D).norm()*(C-D).norm() )  );
        double DCB = std::acos(  (D-C).dot(B-C)/( (D-C).norm()*(B-C).norm() )  );
        double CBD = std::acos(  (C-B).dot(D-B)/( (C-B).norm()*(D-B).norm() )  );
        file_angles<<ABC<<std::endl<<BCA<<std::endl<<CAB<<std::endl;
        // triangolo ACD
        double ACD = std::acos(  (A-C).dot(D-C)/( (A-C).norm()*(D-C).norm() )  );
        double CDA = std::acos(  (C-D).dot(A-D)/( (C-D).norm()*(A-D).norm() )  );
        double DAC = std::acos(  (D-A).dot(C-A)/( (D-A).norm()*(C-A).norm() )  );
        file_angles<<ABC<<std::endl<<BCA<<std::endl<<CAB<<std::endl;
    }
    file_angles.close();


    // la mesh viene scritta nel formato vtk
    std::cout<<"comincia la scrittura del file .vtk\n";
    std::string file_vtk_name = "cubo.vtu";
    n_elements = mesh_simp.n_elements();
    n_nodes = mesh_simp.n_nodes();
    constexpr int dim = 3;

    /*std::ofstream vtu;
    vtu.open("../../../meshes/"+file_vtk_name);

    vtu << "<?xml version=\"1.0\"?>\n";
    vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtu << "<UnstructuredGrid>\n";
    vtu << "<Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_elements << "\">\n";

    // Nodi della mesh (punti VTK)
    vtu << "<Points>\n";
    vtu << "<DataArray type=\"Float64\" Name=\"nodes\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">\n";
    for(unsigned i = 0; i < n_nodes; ++i){
        for(unsigned j = 0; j < 3; ++j){
            vtu << mesh_simp.nodes()(i, j) << " ";
        }
        vtu << "\n";
    }
    vtu << "</DataArray>\n";
    vtu << "</Points>\n";

    // Elementi della mesh (celle VTK)
    vtu << "<Cells>\n";
    // Specifica la connettività dei punti. Tutti gli elenchi dei punti delle celle sono concatenati insieme.
    vtu << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(unsigned i = 0; i < n_elements; ++i){
        for(unsigned j = 0; j < 4; ++j){
            vtu << mesh_simp.elements()(i, j) << " ";
        }
        vtu << "\n";
    }
    vtu << "</DataArray>\n";

    // Specifica l'offset nella matrice di connettività per la fine di ogni cella
    vtu << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset =  4;  // Numero di punti per elemento per mesh triangolari e tetraedriche
    for (int i = 1; i <= n_elements; ++i) {
        vtu << i * offset << " ";
    }
    vtu << "\n</DataArray>\n";

    // Specifica il tipo di ogni cella
    vtu << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    int vtk_cell_type =  10;  // Codice VTK per elementi triangolari e tetraedrici
    for (int i = 0; i < n_elements; ++i) {
        vtu << vtk_cell_type << " ";
    }
    vtu << "\n</DataArray>\n";
    vtu << "</Cells>\n";
    vtu << "</Piece>\n";
    vtu << "</UnstructuredGrid>\n";
    vtu << "</VTKFile>\n";

    vtu.close();

    std::cout<<"scrivere il nome del file della qoi\n";
    std::string nome_file;
    std::cin>>nome_file;
    std::ofstream file_qoi("../../../meshes/dati/3D/"+nome_file+".txt");
    auto active_elems = simp.active_elems();
    for(auto elem_id : active_elems)
        file_qoi<<data_disp_cost.qoi_[elem_id]<<std::endl;
    file_qoi.close();
    std::cout<<"numero totale di elementi: "<<active_elems.size()<<"\n";*/
}