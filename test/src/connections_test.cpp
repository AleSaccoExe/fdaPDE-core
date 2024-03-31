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
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <fdaPDE/mesh.h>
#include <fdaPDE/utils.h>
#include <gtest/gtest.h>   // testing framework

#include <memory>
#include <random>
#include <set>
#include <unordered_set>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>

using fdapde::core::Element;
using fdapde::core::Mesh;

#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MESH_TYPE_LIST;
using fdapde::testing::MeshLoader;
using namespace fdapde::core;
using namespace std;

/*
OK:
- get_facet

FORSE OK:
- facets_connected_to_node
- get_elem_to_facet


DA DEBUGGARE:
- nodes_on_facet
- elems_on_facet. DA ELIMINARE??
*/
std::vector<int> generateUniqueRandomNumbers(int N, int M) {
    std::vector<int> numbers;
    // Riempire un vettore con numeri da 0 a M-1
    for (int i = 0; i < M; ++i) {
        numbers.push_back(i);
    }
    // Mescolare il vettore per avere ordine casuale
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(numbers.begin(), numbers.end(), g);
    // Ridimensionare il vettore per ottenere solo i primi N numeri
    numbers.resize(N);
    return numbers;
}



/*
TEST(connections_test, test_1)
{
    unsigned id_facet_test = 25;
    unsigned id_node_test = 30;
    unsigned id_elem_test = 15;
    using namespace std::chrono;
    MeshLoader<Mesh<2, 3>> meshloader("surface");
    unsigned n_facets = meshloader.mesh.n_facets();
    high_resolution_clock::time_point start = high_resolution_clock::now();
    std::vector<int> facets_to_collapse = generateUniqueRandomNumbers(5, n_facets); 
    auto mesh_simp = meshloader.mesh.simplify({facets_to_collapse.begin(), facets_to_collapse.end()});
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end-start).count();
    cout << "tempo per simplify: " << duration << " ms" << endl;


    std::ofstream file_elems("../../../meshes/elems_simp.txt");
    std::ofstream file_nodes("../../../meshes/nodes_simp.txt");

    // Verifica che il file sia stato aperto correttamente
    if (file_elems.is_open()) {
        // Scrittura della matrice nel file
        file_elems << mesh_simp.elements() << std::endl;

        // Chiusura del file
        file_elems.close();

        std::cout << "Matrice degli elementi scritta con successo nel file." << std::endl;
    } else {
        std::cerr << "Impossibile aprire il file degli elementi per la scrittura." << std::endl;
    }

    if (file_nodes.is_open()) {
        // Scrittura della matrice nel file
        file_nodes << mesh_simp.nodes() << std::endl;

        // Chiusura del file
        file_nodes.close();

        std::cout << "Matrice dei nodi scritta con successo nel file." << std::endl;
    } else {
        std::cerr << "Impossibile aprire il file dei nodi per la scrittura." << std::endl;
    }
    EXPECT_TRUE(true);
}
*/

TEST(connections_test, test_2)
{
    using namespace std::chrono;
    MeshLoader<Mesh<3, 3>> meshloader("unit_sphere");
    unsigned n_facets = meshloader.mesh.n_facets();

    std::vector<int> facets_to_collapse = generateUniqueRandomNumbers(5, n_facets);

    auto mesh_simp = meshloader.mesh.simplify({facets_to_collapse.begin(), facets_to_collapse.end()});


    std::ofstream file_elems("../../../meshes/sphere_elems_simp.txt");
    std::ofstream file_nodes("../../../meshes/sphere_nodes_simp.txt");

    // Verifica che il file sia stato aperto correttamente
    if (file_elems.is_open()) {
        // Scrittura della matrice nel file
        file_elems << mesh_simp.elements() << std::endl;

        // Chiusura del file
        file_elems.close();

        std::cout << "Matrice degli elementi scritta con successo nel file." << std::endl;
    } else {
        std::cerr << "Impossibile aprire il file degli elementi per la scrittura." << std::endl;
    }

    if (file_nodes.is_open()) {
        // Scrittura della matrice nel file
        file_nodes << mesh_simp.nodes() << std::endl;

        // Chiusura del file
        file_nodes.close();

        std::cout << "Matrice dei nodi scritta con successo nel file." << std::endl;
    } else {
        std::cerr << "Impossibile aprire il file dei nodi per la scrittura." << std::endl;
    }

    EXPECT_TRUE(true);
}
