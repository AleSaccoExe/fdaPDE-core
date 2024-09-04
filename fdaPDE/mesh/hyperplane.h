#ifndef __HYPERPLANE_H__
#define __HYPERPLANE_H__

#include <array>

#include "../utils/assert.h"
#include "../utils/symbols.h"

namespace fdapde {
namespace core {


// orthogonal projection of vector v over u
template <int N> constexpr SVector<N> orthogonal_project(const SVector<N>& v, const SVector<N>& u) {
    return (v.dot(u) / u.squaredNorm() * u.array()).matrix();
}

// a template class representing an M-dimensional plane embedded in an N-dimensional space
template <int M, int N> class HyperPlane {
    static_assert(M <= N);
   private:
    // let x_1, x_2, \ldots, x_{N+1} be a set of N+1 points through which the plane passes
    SMatrix<N, M> basis_;   // matrix [x2 - x1, x3 - x1, \ldots, x_{N+1} - x1] of vectors generating the plane
    double offset_;         // hyperplane's intercept (the d in the equation ax + by + cz + d = 0)
    SVector<N> normal_;     // normal vector to the hyperplane
    SVector<N> p_;          // point through which the plane is guaranteed to pass
   public:
    // constructor
    HyperPlane() = default;
    // constructs a hyperplane passing through 2 points, e.g. a line
    HyperPlane(const SVector<N>& x1, const SVector<N>& x2) : p_(x1) {
        fdapde_static_assert(M == 1, THIS_METHOD_IS_ONLY_FOR_LINES);
	SVector<N> tmp = x2 - x1;
        basis_.col(0) = tmp.normalized();
        normal_ << -tmp[0], tmp[1];
	normal_.normalized();
        offset_ = -x1.dot(normal_);
    }
    // constructs an hyperplane passing through 3 (non-collinear) points, e.g. a plane
    HyperPlane(const SVector<N>& x1, const SVector<N>& x2, const SVector<N> x3) : p_(x1) {
        fdapde_static_assert(M == 2, THIS_METHOD_IS_ONLY_FOR_PLANES);
        basis_.col(0) = (x2 - x1).normalized();
        basis_.col(1) = ((x3 - x1) - orthogonal_project(SVector<N>(x3 - x1), SVector<N>(basis_.col(0)))).normalized();
    if constexpr(N!=2)
	   normal_ = ((x2 - x1).cross(x3 - x1)).normalized();
    else{
        /*normal_ = basis_.fullPivLu().kernel();
	    offset_ = -x1.dot(normal_);*/
        static_assert(N==2);
        normal_[0] = 0.0;
        SVector<N> v = x2-x1;
        SVector<N> w = x3-x1;
        normal_[1] = v[0]*w[1]-v[1]*w[0];
        }
    }
    HyperPlane(const SMatrix<N, M + 1>& coords) : p_(coords.col(0)) {
        basis_ = coords.rightCols(M).colwise() - coords.col(0);
    // basis orthonormalization via modified Gram-Schmidt method
        basis_.col(0) /= basis_.col(0).norm();
        for (int i = 1; i < M; ++i) {
            for (int j = 0; j < i; ++j) {
                basis_.col(i) = basis_.col(i) - orthogonal_project<N>(basis_.col(i), basis_.col(j));
            }
            basis_.col(i) /= basis_.col(i).norm();
        }
        // normal_ = basis_.fullPivLu().kernel();   // normal to the hyperplane is any element in the null space of basis_
        SVector<N> u = coords.col(1) - coords.col(0);
        SVector<N> v = coords.col(2) - coords.col(0);
        SVector<N> w = coords.col(3) - coords.col(0);
        SVector<N> l1 = u.cross(v);
        SVector<N> l2 = v.cross(w);
        normal_ = l1.cross(l2);
        normal_/=normal_.norm();
        offset_ = -coords.col(0).dot(normal_);
    }
    // constructors from matrix coordinates
    /*HyperPlane(const SMatrix<N, 2>& coords) requires (M == 1) : HyperPlane(coords.col(0), coords.col(1)) { }
    HyperPlane(const SMatrix<N, 3>& coords) requires (M == 2) :
      HyperPlane(coords.col(0), coords.col(1), coords.col(2)) { }
    // general hyperplane constructor
    HyperPlane(const SMatrix<N, M + 1>& coords) requires(M > 2) : p_(coords.col(0)) {
        basis_ = coords.rightCols(M).colwise() - coords.col(0);
	// basis orthonormalization via modified Gram-Schmidt method
        basis_.col(0) /= basis_.col(0).norm();
        for (int i = 1; i < M; ++i) {
            for (int j = 0; j < i; ++j) {
                basis_.col(i) = basis_.col(i) - orthogonal_project(basis_.col(i), basis_.col(j));
            }
            basis_.col(i) /= basis_.col(i).norm();
        }
        normal_ = basis_.fullPivLu().kernel();   // normal to the hyperplane is any element in the null space of basis_
        offset_ = -coords.col(0).dot(normal_);
    }*/
    // projection
    SVector<M> project_onto(const SVector<N>& x) {
        if constexpr (M == N) {
            return x;
        } else {
            // build the projection onto the space spanned by basis_
            SVector<M> proj;
            for (int i = 0; i < M; ++i) { proj[i] = (x - p_).dot(basis_.col(i)); }
            return proj;
        }
    }
    SVector<N> project(const SVector<N>& x) const {
        if constexpr (M == N) {
            return x;
        } else {
            return basis_ * basis_.transpose() * (x - p_) + p_;
        }
    }
    // euclidean distance of a point from the space
    double distance(const SVector<N>& x) { return (x - project(x)).norm(); }
    SVector<N> operator()(const SVector<M>& coeffs) const {
        SVector<N> res = p_;
        for (int i = 0; i < M; ++i) res += coeffs[i] * basis_.col(i);
        return res;
    }
    // normal direction to the hyperplane
    const SVector<N>& normal() { return normal_; }  
};

}   // namespace core
}   // namespace fdapde

#endif   // __VECTOR_SPACE_H__