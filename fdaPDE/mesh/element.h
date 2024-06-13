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

#ifndef __ELEMENT_H__
#define __ELEMENT_H__

#include <Eigen/LU>
#include <array>

#include "hyperplane.h"
#include "../linear_algebra/vector_space.h"
#include "../utils/assert.h"
#include "../utils/compile_time.h"
#include "../utils/symbols.h"
#include "../utils/intersection.h"
#include "mesh_utils.h"

namespace fdapde {
namespace core {

// M local dimension, N embedding dimension, R order of the element (defaulted to linear finite elements)
template <int M, int N> class Element;   // forward declaration

// A mesh element. M: local dimension, N: embedding dimension
template <int M, int N> class Element {
   private:
    typedef typename std::conditional<is_network<M, N>::value, std::vector<int>, std::array<int, M + 1>>::type
      NeighborsContainer;

    int ID_;                                              // ID of this element
    std::array<int, ct_nvertices(M)> node_ids_ {};        // ID of nodes composing the element
    std::array<SVector<N>, ct_nvertices(M)> coords_ {};   // nodes coordinates (1-1 mapped with node_ids_)
    NeighborsContainer neighbors_ {};                     // ID of neighboring elements
    bool boundary_;                                       // true if the element has at least one vertex on the boundary
    double measure_ = 0;                                  // measure of the element

    // affine transformations from cartesian to barycentric coordinates and viceversa
    SMatrix<N, M> J_;       // [J_]_ij = (coords_[j][i] - coords_[0][i])
    SMatrix<M, N> inv_J_;   // J^{-1} (Penrose pseudo-inverse for manifold elements)
   public:
    // constructor
    Element() = default;
    Element(int ID, const std::array<int, ct_nvertices(M)>& node_ids,
	    const std::array<SVector<N>, ct_nvertices(M)>& coords, const NeighborsContainer& neighbors, bool boundary);

    // getters
    const std::array<SVector<N>, ct_nvertices(M)>& coords() const { return coords_; }
    const NeighborsContainer& neighbors() const { return neighbors_; }
    const int ID() const { return ID_; }
    Eigen::Matrix<double, N, M> barycentric_matrix() const { return J_; }
    Eigen::Matrix<double, M, N> inv_barycentric_matrix() const { return inv_J_; }
    std::array<int, ct_nvertices(M)> node_ids() const { return node_ids_; }
    double measure() const { return measure_; }   // measure of the element

    // check if x is contained in the element
    template <bool is_manifold = (N != M)>
    typename std::enable_if<!is_manifold, bool>::type contains(const SVector<N>& x) const;
    template <bool is_manifold = (N != M)>   // specialization for manifold case
    typename std::enable_if<is_manifold, bool>::type contains(const SVector<N>& x) const;

    SVector<M + 1> to_barycentric_coords(const SVector<N>& x) const;   // move x into barycentric reference system
    SVector<N> mid_point() const;                                      // computes midpoint of the element
    std::pair<SVector<N>, SVector<N>> bounding_box() const;   // compute the smallest rectangle containing this element
    bool is_on_boundary() const { return boundary_; }   // true if the element has at least one node on the boundary
    VectorSpace<M, N> spanned_space() const;            // vector space passing throught this element

    // allow range for over element's coordinates
    typename std::array<SVector<N>, ct_nvertices(M)>::const_iterator begin() const { return coords_.cbegin(); }
    typename std::array<SVector<N>, ct_nvertices(M)>::const_iterator end() const { return coords_.cend(); }

    // compile time constants
    enum {
        n_vertices = ct_nvertices(M),
	local_dimension = M,
	embedding_dimension = N
    };

    // AGGIUNTE MIE 
    // ===========>
    bool intersection(const Element<M, N> & other_el) const;
    // il metodo trova il nodo dell'elemento in facet[1:end] e lo sostituisce con l'id facet[0]
    void find_and_change_node(const std::array<int, M> & facet);
    // hyperplane passing throught this element
    HyperPlane<M, N> hyperplane() const;
    // ================>
    // FINE AGGIUNTE MIE
};

// implementation details

template <int M, int N>
Element<M, N>::Element(
  int ID, const std::array<int, ct_nvertices(M)>& node_ids, const std::array<SVector<N>, ct_nvertices(M)>& coords,
  const NeighborsContainer& neighbors, bool boundary) :
    ID_(ID), node_ids_(node_ids), coords_(coords), neighbors_(neighbors), boundary_(boundary) {
    // precompute barycentric coordinate matrix for fast access
    // use first point as reference
    SVector<N> ref = coords_[0];
    for (std::size_t j = 0; j < M; ++j) { J_.col(j) = coords_[j + 1] - ref; }

    // precompute and cache inverse of barycentric matrix
    if constexpr (N == M) {
        inv_J_ = J_.inverse();
    } else {
        // for manifold domains compute the generalized Penrose inverse
        inv_J_ = (J_.transpose() * J_).inverse() * J_.transpose();
    }
    // precompute element measure
    if constexpr (M == N) {   // non-manifold case
        measure_ = std::abs(J_.determinant()) / (ct_factorial(M));
    } else {
        if constexpr (M == 2)   // surface element, compute area of 3D triangle
            measure_ = 0.5 * J_.col(0).cross(J_.col(1)).norm();
        if constexpr (M == 1)   // network element, compute length of 2D segment
            measure_ = J_.col(0).norm();
    }
}

template <int M, int N>
SVector<M + 1> Element<M, N>::to_barycentric_coords(const SVector<N>& x) const {
    // solve: barycenitrc_matrix_*z = (x-ref)
    SVector<M> z = inv_J_ * (x - coords_[0]);

    SVector<M + 1> result;
    result << SVector<1>(1 - z.sum()), z;
    return result;
}

template <int M, int N> SVector<N> Element<M, N>::mid_point() const {
    // The center of gravity of an element has all its barycentric coordinates equal to 1/(M+1).
    // The midpoint of an element can be computed by mapping it to a cartesian reference system
    SVector<M> barycentric_mid_point;
    barycentric_mid_point.fill(1.0 / (M + 1));   // avoid implicit conversion to int
    return J_ * barycentric_mid_point + coords_[0];
}

template <int M, int N>
std::pair<SVector<N>, SVector<N>> Element<M, N>::bounding_box() const {
    // define lower-left and upper-right corner of bounding box
    SVector<N> ll, ur;
    // project each vertex coordinate on the reference axis
    std::array<std::array<double, ct_nvertices(M)>, N> proj_coords;

    for (std::size_t j = 0; j < ct_nvertices(M); ++j) {
        for (std::size_t dim = 0; dim < N; ++dim) { proj_coords[dim][j] = coords_[j][dim]; }
    }
    // take minimum and maximum value along each dimension, i.e. the lower-left and upper-right corner of the box
    for (std::size_t dim = 0; dim < N; ++dim) {
        ll[dim] = *std::min_element(proj_coords[dim].begin(), proj_coords[dim].end());
        ur[dim] = *std::max_element(proj_coords[dim].begin(), proj_coords[dim].end());
    }
    return std::make_pair(ll, ur);
}

// check if a point is contained in the element
template <int M, int N>
template <bool is_manifold>
typename std::enable_if<!is_manifold, bool>::type Element<M, N>::contains(const SVector<N>& x) const {
    // A point is inside the element if all its barycentric coordinates are positive
    return (to_barycentric_coords(x).array() >= -100 * std::numeric_limits<double>::epsilon()).all();
}

template <int M, int N> VectorSpace<M, N> Element<M, N>::spanned_space() const {
    // build a basis for the space containing the element
    std::array<SVector<N>, M> basis;
    for (size_t i = 0; i < M; ++i) { basis[i] = (coords_[i + 1] - coords_[0]); }
    return VectorSpace<M, N>(basis, coords_[0]);
}

// specialization for manifold elements of contains() routine.
template <int M, int N>
template <bool is_manifold>
typename std::enable_if<is_manifold, bool>::type Element<M, N>::contains(const SVector<N>& x) const {
    // check if the point is contained in the affine space spanned by the element
    VectorSpace<M, N> vs = spanned_space();
    if (vs.distance(x) > 10 * std::numeric_limits<double>::epsilon()) {
        return false;   // distance too large, not contained
    }
    // check if its barycentric coordinates are all positive
    return (to_barycentric_coords(x).array() >= -10 * std::numeric_limits<double>::epsilon()).all();
}

// AGGIUNTE MIE
// ===========>
template<>
bool Element<2, 3>::intersection(const Element<2, 3>& other_el) const
{
    // Compute the normal to the triangle and the RHS 
    // of the equation of the plane the triangle lies in
    const SVector<3>& A = coords_[0];
    const SVector<3>& B = coords_[1];
    const SVector<3>& C = coords_[2];

    const SVector<3>& D = other_el.coords_[0];
    const SVector<3>& E = other_el.coords_[1];
    const SVector<3>& F = other_el.coords_[2];
    
    SVector<3> N = (B - A).cross(C - B);
    N.normalize();
    double RHS = N.dot(A);

    // Extract the maximum coordinate of N
    int z = getMaxCoord(N); // coordinata massima
    int x = (z+1) % 3;
    int y = (z+2) % 3;

    // Project the triangle ABC onto the "xy"-plane
    Eigen::Vector2d a(A[x],A[y]);
    Eigen::Vector2d b(B[x],B[y]);
    Eigen::Vector2d c(C[x],C[y]);
    
    auto de_abc = intSegTri(D,E,a,b,c,N,RHS,x,y);
    //utility::printIntersectionType(de_abc,"de and abc");
    auto ef_abc = intSegTri(E,F,a,b,c,N,RHS,x,y);
    //utility::printIntersectionType(ef_abc,"ef and abc");
    auto fd_abc = intSegTri(F,D,a,b,c,N,RHS,x,y);
    //utility::printIntersectionType(fd_abc,"fd and abc");
    
    // If at least one segment-triangle intersection
        // is not conformal, the triangle-triangle intersection is 
        // invalid too
        if ((de_abc == IntersectionType::INVALID) ||
            (ef_abc == IntersectionType::INVALID) ||
            (fd_abc == IntersectionType::INVALID))
            return true;
            
        // 
        // Triangle DEF as reference
        //
        
            // Compute the normal to the triangle and the RHS 
            // of the equation of the plane the triangle lies in
            N = (E - D).cross(F - E);
            N.normalize();
            RHS = N.dot(D);
        
            // Extract the maximum coordinate of N
            z = getMaxCoord(N);
            x = (z+1) % 3;
            y = (z+2) % 3;
        
            // Project the triangle ABC onto the "xy"-plane
            Eigen::Vector2d d(D[x],D[y]);
            Eigen::Vector2d e(E[x],E[y]);
            Eigen::Vector2d f(F[x],F[y]);
        
        // Test intersection of each edge of DEF with ABC
        auto ab_def = intSegTri(A,B,d,e,f,N,RHS,x,y);
        //utility::printIntersectionType(ab_def,"ab and def");
        auto bc_def = intSegTri(B,C,d,e,f,N,RHS,x,y);
        //utility::printIntersectionType(bc_def,"bc and def");
        auto ca_def = intSegTri(C,A,d,e,f,N,RHS,x,y);
        //utility::printIntersectionType(ca_def,"ca and def");
                
        // If at least one segment-triangle intersection
        // is not conformal, the triangle-triangle intersection is 
        // invalid too
        if ((ab_def == IntersectionType::INVALID) ||
            (bc_def == IntersectionType::INVALID) ||
            (ca_def == IntersectionType::INVALID))
            return true;
            
        // If here, the triangles do not intersection one each other
        // or they do in a conformal way
        return false;
    }
template<int M, int N>
void Element<M, N>::find_and_change_node(const std::array<int, M> & facet)
{
    for(int i = 0; i < n_vertices; ++i)
        for(int j = 1; j < M; ++j)
            if(facet[j] == node_ids_[i]){
                node_ids_[i] = facet[0];
                return;
            }   
}
template<>
HyperPlane<2, 3> Element<2, 3>::hyperplane() const
{
    return HyperPlane<2, 3>(coords_[0], coords_[1], coords_[2]);
}

template<>
HyperPlane<2, 2> Element<2, 2>::hyperplane() const
{
    return HyperPlane<2, 2>(coords_[0], coords_[1], coords_[2]);
}

template<>
HyperPlane<3, 3> Element<3, 3>::hyperplane() const
{
    SMatrix<3, 4> coords;
    for(int i = 0; i < 4; ++i)
        coords.col(i) = coords_[i];
    return HyperPlane<3, 3>(coords);
}

// ================>
// FINE AGGIUNTE MIE

}   // namespace core
}   // namespace fdapde

#endif   // __ELEMENT_H__
