#ifndef __MESH_H__
#define __MESH_H__

#include <Eigen/Core>
#include <type_traits>
#include <vector>
#include <memory>
#include <array>

#include "Element.h"
#include "../utils/IO/CSVReader.h"
#include "../utils/Symbols.h"

namespace fdaPDE{
namespace core{
namespace MESH{

  // Mesh is the access point to mesh informations (triangulated domains). It offers an abstraction layer allowing to reason
  // on the mesh from a geometrical perspective, i.e. without considering its internal representation in memory. The class is able to transparently
  // handle manifold and non manifold meshes, exposing the same interface in any case. Two non-type template parameters are used:
  //     * M: local dimension of the mesh (dimension of the space to which a mesh element belongs)
  //     * N: embeddding dimension of the mesh (dimension of the space where the whole mesh lies)
  // if M != N the mesh is a manifold. Currently are implemented:
  //     * 1.5D meshes (linear newtorks,    M=1, N=2)
  //     * 2D meshes   (planar domains,     M=2, N=2)
  //     * 2.5D meshes (surfaces,           M=2, N=3)
  //     * 3D meshes   (volumetric domains, M=3, N=3)
  // A third non-type template parameter R is used to keep track of the order of the elements contained in the mesh object. Mesh are defaulted to
  // use finite linear elements (R = 1)
  
  // NB about internal implementaiton: special care is needed in the development of linear networks, since differently from any other case the number
  // of neighboing elements is not known at compile time. This implies the usage of specialized data structures wrt any other case

  // trait to detect if the mesh is a manifold
  template <unsigned int M, unsigned int N>
  struct is_manifold{
    static constexpr bool value = (M != N);
  };
  
  // trait to detect if the mesh is a linear network
  template <unsigned int M, unsigned int N>
  struct is_linear_network{
    static constexpr bool value = std::conditional<
      (M == 1 && N == 2), std::true_type, std::false_type
      >::type::value;
  };

  // trait to select a proper neighboring storage structure depending on the type of mesh. In case of linear networks this information is stored as
  // a sparse matrix where entry (i,j) is set to 1 if and only if elements i and j are neighbors
  template <unsigned int M, unsigned int N>
  struct neighboring_structure{
    using type = typename std::conditional<
      is_linear_network<M, N>::value, SpMatrix<int>, DMatrix<int>
      >::type;
  };
  
  template <unsigned int M, unsigned int N, unsigned int R = 1>
  class Mesh{
  private:
    // coordinates of points constituting the vertices of mesh elements
    DMatrix<double> points_{};
    unsigned int numNodes_ = 0;
    // matrix of edges. Each row of the matrix contains the row numbers in points_ matrix
    // of the points which form the edge
    DMatrix<int> edges_{};
    // matrix of triangles in the triangulation. Each row of the matrix contains the row
    // numbers in points_ matrix of the points which form the triangle
    DMatrix<int> elements_{};
    unsigned int numElements_ = 0;
    // in case of non linear-networks neighbors_ is a dense matrix where row i contains the indexes as row number in triangles_ matrix of the
    // neighboring triangles to triangle i (all triangles in the triangulation which share an edge with i). In case of linear-newtorks neighbors_
    // is a sparse matrix where entry (i,j) is set to 1 iff i and j are neighbors
    typename neighboring_structure<M, N>::type neighbors_{};
    // store boundary informations. This is a vector of binary coefficients such that, if element j is 1
    // then mesh node j is on boundary, otherwise 0
    DMatrix<int> boundary_{};
    
    // store min-max values for each dimension of the mesh
    std::array<std::pair<double, double>, N> range_{};
    // is often required to access just to the minimum value along each dimension and to the quantity
    // kk_[dim] = 1/(max[dim] - min[dim]) = 1/(meshRange[dim].second - meshRange[dim].first). Compute here once and cache results for efficiency
    std::array<double, N> minRange_{};
    std::array<double, N> kk_{}; // kk_[dim] = 1/(range[dim].second - range[dim].first)

    // elements informations are computed once and cached here for fast re-access
    std::vector<std::shared_ptr<Element<M,N,R>>> cache_{};
    void fill_cache();
  public:
    Mesh() = default;
    // construct from .csv files, strings are names of file where raw data is contained
    Mesh(const std::string& points,    const std::string& edges, const std::string& triangles,
	 const std::string& neighbors, const std::string& boundary);

    // construct directly from eigen matrices
    Mesh(const DMatrix<double>& points, const DMatrix<int>& edges, const DMatrix<int>& elements,
	 const typename neighboring_structure<M, N>::type& neighbors, const DMatrix<int>& boundary);
    
    // returns an element object given its ID (its row number in the triangles_ matrix) from raw (matrix-like) informations
    std::shared_ptr<Element<M,N,R>> element(unsigned int ID) const;
    // return the coordinate of a node given its ID (its row number in the points_matrix)
    SVector<N> node(unsigned int ID) const;

    // allow range-for loop over mesh elements
    struct iterator{
    private:
      friend Mesh;
      const Mesh* meshContainer; // pointer to mesh object
      int index;           // keep track of current iteration during for-loop
      // constructor
      iterator(const Mesh* container_, int index_) : meshContainer(container_), index(index_) {}; 
    public:
      // just increment the current iteration and return this iterator
      iterator& operator++() {
	++index;
	return *this;
      }
      // dereference the iterator means to create Element object at current index
      std::shared_ptr<Element<M,N,R>> operator*() {
	return meshContainer->element(index);
      }
      // two iterators are different when their indexes are different
      friend bool operator!=(const iterator& lhs, const iterator& rhs) {
	return lhs.index != rhs.index;
      }

      // const version to enable const auto& syntax
      std::shared_ptr<Element<M,N,R>> operator*() const { return meshContainer->element(index); }
    };
    // provide begin() and end() methods
    iterator begin() const { return iterator(this, 0); }
    iterator end()   const { return iterator(this, elements_.rows()); }

    // return true if the given node is on boundary, false otherwise
    bool isOnBoundary(size_t j) const { return boundary_(j) == 1; }
    
    // getters
    unsigned int elements() const { return numElements_; }
    unsigned int nodes() const { return numNodes_; }
    std::array<std::pair<double, double>, N> range() const { return range_; }
    std::array<double, N> lowerBound() const { return minRange_; }
    std::array<double, N> kk() const { return kk_; }
    
    // expose compile time informations to outside
    static constexpr bool manifold = is_manifold<M, N>::value;
    static constexpr unsigned int local_dimension = M;
    static constexpr unsigned int embedding_dimension = N;
    static constexpr unsigned int order = R;
  };

  // export some aliases
  template <unsigned int R=1> using Mesh2D = Mesh<2,2,R>;
  template <unsigned int R=1> using Mesh3D = Mesh<3,3,R>;
  // manifold cases
  template <unsigned int R=1> using SurfaceMesh = Mesh<2,3,R>;
  template <unsigned int R=1> using NetworkMesh = Mesh<1,2,R>;

#include "Mesh.tpp"
}}}
  
#endif // __MESH_H__
