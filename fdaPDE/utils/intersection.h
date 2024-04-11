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

#ifndef __INTERSECTION_H__
#define __INTERSECTION_H__

#include "compile_time.h"
#include "symbols.h"

namespace fdapde {
namespace core {

/*! Relative position of a point with respect to a segment. */
enum class Point2Seg {EXTERN, INTERN, ONVERTEX};

/*! Relative position of a point with respect to a triangle. */
enum class Point2Tri {EXTERN, INTERN, ONEDGE, ONVERTEX};

/*! Relative position of a straight line with respect to a plane. */
enum class Line2Plane {PARALLEL, COMPLANAR, INCIDENT};

/*! Type of intersection:
        <ol>
        <li> NONE:      the elements do not intersect;
        <li> VALID:     the elements intersect in a conformal way,
                        i.e. they share a vertex or an entire edge;
        <li> UNVALID:   the elements intersect in a non-conformal way.
        <\ol> */
enum class IntersectionType {NONE, VALID, INVALID};

/*! Test if a segment intersect a plane.
    \param Q    querying end-point of the segment
    \param R    ray end-point of the segment 
    \param N    unit normal to the triangle 
    \param D    RHS term in the equation of the plane
                determined by the triangle
    \return     relative position line-plane
    \return     position of the intersection point \$ p \$
                (if any) in the segment
    \return     \$ t \$ s.t. \$ p = q + t \cdot (r - q) \$ */

template<int N>
unsigned getMaxCoord(const SVector<N> & v)
{
    unsigned max_coord = 0;
    for(int i = 1; i<N; ++i)
        if(abs(v[i]) > abs(v[max_coord]) )
            max_coord = i;
    return max_coord;

}

double getTriArea2d(const Eigen::Vector2d & A, 
        const Eigen::Vector2d & B, const Eigen::Vector2d & C)
    {
        Eigen::Vector2d l1 = B-A;
        Eigen::Vector2d l2 = C-B;
        return 0.5*( l1[0]*l2[1] - l2[0]*l1[1] );
       
       
    }

std::tuple<Line2Plane, Point2Seg, double> intSegPlane
        (const Eigen::Vector3d & Q, const Eigen::Vector3d & R, const Eigen::Vector3d & N, const double & D)
    {
        double TOLL = DOUBLE_TOLERANCE;
        auto l2p = Line2Plane::INCIDENT; 
        auto p2s = Point2Seg::ONVERTEX;
        double t;
        
        // Compute numerator and denumerator of Equation (7.7), p. 228
        double q_plane = D - Q.dot(N);
        double den = (R - Q).dot(N);
        
        //
        // The segment is parallel to the plane
        //
        // This happens iff the denumerator vanishes
        // If the numerator vanishes too: the segment belongs to the plane
        // If the numerator does not vanish: the segment does not belong to the plane
        if ((-TOLL <= den) && (den <= TOLL))
        {
            if ((-TOLL <= q_plane) && (q_plane <= TOLL))
                l2p = Line2Plane::COMPLANAR;
            else
                l2p = Line2Plane::PARALLEL;
                
            return std::make_tuple(l2p, p2s, t);
        }
        
        //
        // The segment is not parallel to the plane
        //
        
        // Check if the plane and the segment intersect
        // in the querying vertex of the segment
        if ((-TOLL <= q_plane) && (q_plane <= TOLL))
            return std::make_tuple(l2p, p2s, 0.);
            
        // Check if the plane and the segment intersect
        // in the ray vertex of the segment
        double r_plane = D - R.dot(N);
        if ((-TOLL <= r_plane) && (r_plane <= TOLL))
            return std::make_tuple(l2p, p2s, 1.);
        
        t = q_plane/den;
        
        // They intersect one each other iff \$ 0 \leq t \leq 1 \$
        if ((t < -TOLL) || (t > 1.+TOLL))
        {
            p2s = Point2Seg::EXTERN;
            return std::make_tuple(l2p, p2s, t);
        }
        
        // The intersection point is strictly inside the segment
        // iff \$ 0 < t < 1 \$
        if ((TOLL < t) && (t < 1.-TOLL))
        {
            p2s = Point2Seg::INTERN;
            return std::make_tuple(l2p, p2s, t);
        }
        
        // Due to floating point arithmetic, t may turn out to be 
        // close to 0 or 1, i.e. $-TOLL \leq t \leq TOLL$ or 
        // $1-TOLL \leq t \leq 1+TOLL$. To avoid any run time issue, 
        // we should handle this situation properly.
        if ((-TOLL <= t) && (t <= TOLL))
            return std::make_tuple(l2p, p2s, 0.);
        return std::make_tuple(l2p, p2s, 1.);
    }



IntersectionType intSegSeg2d(const Eigen::Vector2d & q1, 
        const Eigen::Vector2d & r1, const Eigen::Vector2d & q2, const Eigen::Vector2d & r2)
    {
        double TOLL = DOUBLE_TOLERANCE;
        // Compute signed area of the triangles q1r1q2, q1r1r2, q2r2q1, q2r2r1
        double q1r1q2 = getTriArea2d(q1,r1,q2);
        double q1r1r2 = getTriArea2d(q1,r1,r2);
        double q2r2q1 = getTriArea2d(q2,r2,q1);
        double q2r2r1 = getTriArea2d(q2,r2,r1);
        
        //
        // The segments are collinear
        //
        // We have to check if either q2 or r2 belongs to q1-r1
        // A point p belongs to a segment qr iff
        //  \$ p = q + t \cdot (r - q) \$
        // with \$ 0 \leq t \leq 1 \$
        //
        // To test if a point p is collinear with other two points q and r:
        // check if the (signed) area of the triangle qrp is zero
        if (((-TOLL <= q1r1q2) && (q1r1q2 <= TOLL))
            && ((-TOLL <= q1r1r2) && (q1r1r2 <= TOLL)))
        {
            double t_q2, t_r2;
            
            // Properly check for horizontal or vertical segments
            auto den = r1[0] - q1[0];
            if ((den < -TOLL) || (den > TOLL))
            {
                t_q2 = (q2[0] - q1[0])/den;
                t_r2 = (r2[0] - q1[0])/den;
            }
            else
            {
                den = r1[1] - q1[1];
                t_q2 = (q2[1] - q1[1])/den;
                t_r2 = (r2[1] - q1[1])/den;
            }
            
            // Check if q2 and/or r2 are internal to q1-r1
            if (((TOLL < t_q2) && (t_q2 < 1.-TOLL)) ||
                ((TOLL < t_r2) && (t_r2 < 1.-TOLL)))
                return IntersectionType::INVALID;
                
            // If two vertices coincide: we need to check if
            // either q1 or r1 are internal to q2-r2
            if (((-TOLL <= t_q2) && (t_q2 <= TOLL)) ||
                ((1.-TOLL <= t_q2) && (t_q2 <= 1.+TOLL)) ||
                ((-TOLL <= t_r2) && (t_r2 <= TOLL)) ||
                ((1.-TOLL <= t_r2) && (t_r2 <= 1.+TOLL)))
            {
                double t_q1, t_r1;
            
                // Properly check for horizontal or vertical segments
                auto den = r2[0] - q2[0];
                if ((den < -TOLL) || (den > TOLL))
                {
                    t_q1 = (q1[0] - q2[0])/den;
                    t_r1 = (r1[0] - q2[0])/den;
                }
                else
                {
                    den = r2[1] - q2[1];
                    t_q1 = (q1[1] - q2[1])/den;
                    t_r1 = (r1[1] - q2[1])/den;
                }
            
                // Check if either q1 or r1 are internal to q2-r2...
                if (((TOLL < t_q1) && (t_q1 < 1.-TOLL)) ||
                    ((TOLL < t_r1) && (t_r1 < 1.-TOLL)))
                    return IntersectionType::INVALID;
                
                // ...If not, the segments share a vertex
                return IntersectionType::VALID;
            }
            
            // Last remaining scenarios: the segments do not intersect
            return IntersectionType::NONE;
        }
        
        //
        // The segments do not intersect
        //
        // One of the following two statements must hold:
        //  - q2 and r2 are both on the left or on the right
        //    of the segment q1-r1 
        //  - q1 and r1 are both on the left or on the right
        //    of the segment q2-r2
        //
        // To test if a point p is on the left (right) of a line qr:
        // check if the signed area of the triangle qrp is positive (negative)
        if (((q1r1q2 > TOLL && q1r1r2 > TOLL) || (q1r1q2 < -TOLL && q1r1r2 < -TOLL))
            || ((q2r2q1 > TOLL && q2r2r1 > TOLL) || (q2r2q1 < -TOLL && q2r2r1 < -TOLL)))
            return IntersectionType::NONE; 
        
        //
        // The intersection point is internal to both segments
        //
        // Sufficient condition is that the following two statements
        // simultaneously hold:
        //  - q2 is on the left (right) of q1-r1 and r2 is on the right 
        //    (left) of q1-r1 
        //  - q1 is on the left (right) of q2-r2 and r1 is on the right 
        //    (left) of q2-r2
        if (((q1r1q2 > TOLL && q1r1r2 < -TOLL) || (q1r1q2 < -TOLL && q1r1r2 > TOLL))
            && ((q2r2q1 > TOLL && q2r2r1 < -TOLL) || (q2r2q1 < -TOLL && q2r2r1 > TOLL)))
            return IntersectionType::INVALID;
            
        //
        // The intersection point coincides with a vertex of a segment
        // but it is internal to the other one
        //
        // One of the following two statements must hold:
        //  - q2 is on the left (right) of q1-r1, r2 is on the right
        //    (left) of q1-r1 and either q1 or r1 is collinear with q2 and r2
        //  - q1 is on the left (right) of q2-r2, r1 is on the right
        //    (left) of q2-r2 and either q2 or r2 is collinear with q1 and r1
        if ((((q1r1q2 > TOLL && q1r1r2 < -TOLL) || (q1r1q2 < -TOLL && q1r1r2 > TOLL))
            && (((-TOLL <= q2r2q1) && (q2r2q1 <= TOLL)) || ((-TOLL <= q2r2r1) && (q2r2r1 <= TOLL)))) ||
            (((q2r2q1 > TOLL && q2r2r1 < -TOLL) || (q2r2q1 < -TOLL && q2r2r1 > TOLL))
            && (((-TOLL <= q1r1q2) && (q1r1q2 <= TOLL)) || ((-TOLL <= q1r1r2) && (q1r1r2 <= TOLL)))))
            return IntersectionType::INVALID;
            
        //
        // Only remaining scenarios: the segments share a vertex
        //
        return IntersectionType::VALID;
    }



Point2Tri inTri2d(const Eigen::Vector2d & p, 
        const Eigen::Vector2d & a, const Eigen::Vector2d & b, const Eigen::Vector2d & c)
    {
        double TOLL = DOUBLE_TOLERANCE;
        // Compute signed area of the triangle pab, pbc and pac
        double pab = getTriArea2d(p,a,b);
        double pbc = getTriArea2d(p,b,c);
        double pca = getTriArea2d(p,c,a);
                
        // If the areas are all positive or all negative:
        // the point is internal to the triangle
        if (((pab > TOLL) && (pbc > TOLL) && (pca > TOLL)) ||
            ((pab < -TOLL) && (pbc < -TOLL) && (pca < -TOLL)))
            return Point2Tri::INTERN;
            
        // If two areas are zero: the point coincides with
        // the vertex shared by the associated edges
        bool pab_iszero = (-TOLL <= pab) && (pab <= TOLL);
        bool pbc_iszero = (-TOLL <= pbc) && (pbc <= TOLL);
        bool pca_iszero = (-TOLL <= pca) && (pca <= TOLL);
        if ((pab_iszero && pbc_iszero) ||
            (pbc_iszero && pca_iszero) ||
            (pca_iszero && pab_iszero))
            return Point2Tri::ONVERTEX;
            
        // If one area is zero and the others are concorde:
        // the point lays on an edge
        if ((pab_iszero && (((pbc > 0) && (pca > 0)) || ((pbc < 0) && (pca < 0)))) ||
            (pbc_iszero && (((pab > 0) && (pca > 0)) || ((pab < 0) && (pca < 0)))) ||
            (pca_iszero && (((pab > 0) && (pbc > 0)) || ((pab < 0) && (pbc < 0)))))
            return Point2Tri::ONEDGE;
                        
        // Otherwise, the point does not belong to the triangle
        return Point2Tri::EXTERN;   
    }



IntersectionType intSegTri(const Eigen::Vector3d & Q, const Eigen::Vector3d & R,
        const Eigen::Vector2d & a, const Eigen::Vector2d & b, const Eigen::Vector2d & c, 
        const Eigen::Vector3d & N, const double & D, const unsigned & x, const unsigned & y)
    {
        //
        // Segment-plane intersection
        //
        
        Line2Plane l2p;
        Point2Seg p2s;
        double t;
        std::tie(l2p,p2s,t) = intSegPlane(Q,R,N,D);
        //utility::printLine2Plane(l2p);
        //utility::printPoint2Seg(p2s);
                
        // Necessary condition for the segment to intersect the triangle
        // is that the segment intersects the plane
        if ((l2p == Line2Plane::PARALLEL) || (p2s == Point2Seg::EXTERN))
            return IntersectionType::NONE;
            
        // Project Q and R onto the "xy"-plane
        Eigen::Vector2d q(Q[x],Q[y]);
        Eigen::Vector2d r(R[x],R[y]);
                
        //
        // The segment is parallel to the plane
        //
        
        if (l2p == Line2Plane::COMPLANAR)
        {
            // Test intersection of the segment with each edge
            auto qr_ab = intSegSeg2d(q,r,a,b);
            if (qr_ab == IntersectionType::INVALID)
                return qr_ab;
            
            auto qr_bc = intSegSeg2d(q,r,b,c);
            if (qr_bc == IntersectionType::INVALID)
                return qr_bc;
            
            auto qr_ca = intSegSeg2d(q,r,c,a);
            if (qr_ca == IntersectionType::INVALID)
                return qr_ca;
            
            // The segment may be completely within the triangle
            auto q_abc = inTri2d(q,a,b,c);
            auto r_abc = inTri2d(r,a,b,c);
            if ((q_abc == Point2Tri::INTERN) || (r_abc == Point2Tri::INTERN))
                return IntersectionType::INVALID;
                
            // If here, the triangles do not intersect or
            // they intersect in a conformal way
            // Since we do not really care about that,
            // just classify the intersection as VALID
            return IntersectionType::VALID;
        }
        
        //
        // The segment is not parallel to the plane
        //
        
        // Compute the intersection point
        Eigen::Vector2d p = q + t*(r - q);
        
        // Find the relative position between the intersection
        // point and the triangle
        auto p2t = inTri2d(p,a,b,c);
        //utility::printPoint2Tri(p2t);
        
        // The segment does not intersect the triangle
        if (p2t == Point2Tri::EXTERN)
            return IntersectionType::NONE;
            
        // The segment intersects the triangle in a conformal way
        if ((p2s == Point2Seg::ONVERTEX) && (p2t != Point2Tri::INTERN)) 
            return IntersectionType::VALID;
            
        // Only remaining scenario: the segment intersects the
        // triangle in a non-conformal way
        return IntersectionType::INVALID;
    }



Point2Tri inTri3d(const Eigen::Vector3d & P, 
        const Eigen::Vector3d & A, const Eigen::Vector3d & B, const Eigen::Vector3d & C)
    {
        double TOLL = DOUBLE_TOLERANCE;
        // Compute the normal to the triangle and the RHS 
        // of the equation of the plane the triangle lies in

        SVector<3> N = ((B - A).cross(C - B));
        N.normalize();
        double D = N.dot(A);
        
        //
        // Test if the point belongs to the plane of the triangle
        //
        
        double dist = abs(N.dot(P) - D);
        if (dist > TOLL)
            return Point2Tri::EXTERN;
            
        //
        // The point belongs to the plane of the triangle
        //
        
        // Extract the maximum coordinate of N
        unsigned z = getMaxCoord(N);
        unsigned x = (z+1) % 3;
        unsigned y = (z+2) % 3;
        
        // Project all points onto the "xy"-plane
        Eigen::Vector2d p(P[x],P[y]);
        Eigen::Vector2d a(A[x],A[y]);
        Eigen::Vector2d b(B[x],B[y]);
        Eigen::Vector2d c(C[x],C[y]);
        
        // Check belonging in 2D
        return inTri2d(p,a,b,c);
    }


    // check the intersection of two bounding boxex identified by their NE and SW points
    template<int N>
    bool boxes_intersection(const std::pair<SVector<N>, SVector<N>> & b1, const std::pair<SVector<N>, SVector<N>> & b2)
    {
        for (int i = 0; i < N; ++i) {
            if (b1.second[i] < b2.first[i] || 
                b2.second[i] < b1.first[i]) {
                return false;
            }
        }
        return true;   
    }


}
}

#endif // __INTERSECTION_H__
