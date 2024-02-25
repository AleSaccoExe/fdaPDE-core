// TEST MIO

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

#include <gtest/gtest.h>   // testing framework
#include <cstddef>

#include <fdaPDE/utils.h>
#include <fdaPDE/mesh.h>
#include <fdaPDE/linear_algebra.h>
using fdapde::core::Element;
using fdapde::core::VectorSpace;
//using fdapde::core::circumcenter;

#include "utils/mesh_loader.h"
#include "utils/constants.h"
#include "utils/utils.h"

using fdapde::testing::MESH_TYPE_LIST;
using fdapde::testing::MeshLoader;
using fdapde::testing::DOUBLE_TOLERANCE;
using fdapde::testing::almost_equal;



TEST(IntersectionTest, SegmentPlaneIntersection_1)
{
	// normale al triangolo
	SVector<3> N = {0.5, 0.5, 1.0/sqrt(2.0)};
	// parte destra dell'equazione del piano
	double d = 1.0;
	// estremi del segmento
	SVector<3> Q = {0.0, 0.0, 0.0};
	SVector<3> R = {2.0, 2.0, 2.0};
	fdapde::core::Line2Plane l2p;
	fdapde::core::Point2Seg p2s;
	double t;
	std::tie(l2p, p2s, t) = fdapde::core::intSegPlane(Q, R, N, d);

	//Element<2, 3> el1(0, {1, 2, 3}, {p1, p2, p3}, {}, false);
	//Element<2, 3> el2(0, {1, 2, 3}, {p4, p5, p6}, {}, false);
	EXPECT_TRUE(l2p == fdapde::core::Line2Plane::INCIDENT);
	EXPECT_TRUE(p2s == fdapde::core::Point2Seg::INTERN);

}

TEST(IntersectionTest, SegmentPlaneIntersection_2)
{
	// normale al triangolo
	SVector<3> N = {0.5, 0.5, 1.0/sqrt(2.0)};
	// parte destra dell'equazione del piano
	double d = 1.0;
	// estremi del segmento
	SVector<3> Q = {1.1, 1.1, 1.1};
	SVector<3> R = {2.0, 2.0, 2.0};
	fdapde::core::Line2Plane l2p;
	fdapde::core::Point2Seg p2s;
	double t;
	std::tie(l2p, p2s, t) = fdapde::core::intSegPlane(Q, R, N, d);

	EXPECT_TRUE(l2p == fdapde::core::Line2Plane::INCIDENT);
	EXPECT_TRUE(p2s == fdapde::core::Point2Seg::EXTERN);

}

TEST(IntersectionTest, SegmentPlaneIntersection_3)
{
	// normale al triangolo
	SVector<3> N = {0.5, 0.5, 1.0/sqrt(2.0)};
	// parte destra dell'equazione del piano
	double d = 1.0;
	// estremi del segmento
	SVector<3> Q = {1.0, 1.0, 1.0};
	SVector<3> R = {2.0, 0.0, 1.0};
	fdapde::core::Line2Plane l2p;
	fdapde::core::Point2Seg p2s;
	double t;
	std::tie(l2p, p2s, t) = fdapde::core::intSegPlane(Q, R, N, d);

	EXPECT_TRUE(l2p == fdapde::core::Line2Plane::PARALLEL);
	EXPECT_TRUE(p2s == fdapde::core::Point2Seg::ONVERTEX);
}

TEST(IntersectionTest, SegmentPlaneIntersection_4)
{
	// normale al triangolo
	SVector<3> N = {0.5, 0.5, 1.0/sqrt(2.0)};
	// parte destra dell'equazione del piano
	double d = 1.0;
	// estremi del segmento
	SVector<3> Q = {1.0, 1.0, 0.0};
	SVector<3> R = {2.0, 0.0, 0.0};
	fdapde::core::Line2Plane l2p;
	fdapde::core::Point2Seg p2s;
	double t;
	std::tie(l2p, p2s, t) = fdapde::core::intSegPlane(Q, R, N, d);

	EXPECT_TRUE(l2p == fdapde::core::Line2Plane::COMPLANAR);
	EXPECT_TRUE(p2s == fdapde::core::Point2Seg::ONVERTEX);
}


// TEST SU INTERSEZIONI DI SEGMENTI IN 2D
/*
// Segmenti che si intersecano danno IntersectionType == VALID
TEST(IntersectionTest, SegmentSegment2D_1)
{
	Eigen::Vector2d q1(0.0, 0.0);
	Eigen::Vector2d r1(1.0, 0.0);
	Eigen::Vector2d q2(1.0, 1.0);
	Eigen::Vector2d r2(0.0, 1.0);
	auto int_type = fdapde::core::intSegSeg2d(q1, r1, q2, r2);
	EXPECT_TRUE(int_type == fdapde::core::IntersectionType::VALID);	
}

// segmenti coincidenti danno IntersectionType == INVALID
TEST(IntersectionTest, SegmentSegment2D_2)
{
	Eigen::Vector2d q1(0.0, 0.0);
	Eigen::Vector2d r1(0.0, 0.0);
	Eigen::Vector2d q2(1.0, 1.0);
	Eigen::Vector2d r2(1.0, 1.0);
	auto int_type = fdapde::core::intSegSeg2d(q1, r1, q2, r2);
	EXPECT_TRUE(int_type == fdapde::core::IntersectionType::VALID);	
}
// Segmanti che condivisono un vertice danno IntersectionType == Valid
TEST(IntersectionTest, SegmentSegment2D_3)
{
	Eigen::Vector2d q1(0.0, 0.0);
	Eigen::Vector2d r1(0.0, 0.0);
	Eigen::Vector2d q2(1.0, 1.0);
	Eigen::Vector2d r2(0.0, 1.0);
	auto int_type = fdapde::core::intSegSeg2d(q1, r1, q2, r2);
	EXPECT_TRUE(int_type == fdapde::core::IntersectionType::INVALID);	
}
// Segmenti che non si intersecano danno IntersectionType == NONE
TEST(IntersectionTest, SegmentSegment2D_4)
{
	Eigen::Vector2d q1(20.0, 20.0);
	Eigen::Vector2d r1(1.0, 0.0);
	Eigen::Vector2d q2(21.0, 21.0);
	Eigen::Vector2d r2(0.0, 1.0);
	auto int_type = fdapde::core::intSegSeg2d(q1, r1, q2, r2);
	EXPECT_TRUE(int_type == fdapde::core::IntersectionType::NONE);	
}
*/

// INTERSEZIONI TRA UN PUNTO E UN TRIANGOLO IN 2D

// punto esterno al triangolo
TEST(IntersectionTest, PointTriangle2D_1)
{
	Eigen::Vector2d p(0.0, 0.0);
	Eigen::Vector2d a(1.0, 1.0); 
	Eigen::Vector2d b(1.0, 2.0);
	Eigen::Vector2d c(2.0, 1.0); 
	auto int_type = fdapde::core::inTri2d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::EXTERN);
}

// punto interno al triangolo
TEST(IntersectionTest, PointTriangle2D_2)
{
	Eigen::Vector2d p(1.1, 1.1);
	Eigen::Vector2d a(1.0, 1.0); 
	Eigen::Vector2d b(1.0, 2.0);
	Eigen::Vector2d c(2.0, 1.0); 
	auto int_type = fdapde::core::inTri2d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::INTERN);
}

// punto sul vertice di un triangolo
TEST(IntersectionTest, PointTriangle2D_3)
{
	Eigen::Vector2d p(1.0, 1.0);
	Eigen::Vector2d a(1.0, 1.0); 
	Eigen::Vector2d b(1.0, 2.0);
	Eigen::Vector2d c(2.0, 1.0); 
	auto int_type = fdapde::core::inTri2d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::ONVERTEX);
}

// punto sul lato del triangolo
TEST(IntersectionTest, PointTriangle2D_4)
{
	Eigen::Vector2d p(1.0, 1.5);
	Eigen::Vector2d a(1.0, 1.0); 
	Eigen::Vector2d b(1.0, 2.0);
	Eigen::Vector2d c(2.0, 1.0); 
	auto int_type = fdapde::core::inTri2d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::ONEDGE);
}


// INTERSEZIONI TRA PUNTO E TRIANGOLO IN 3D

// punto esterno al triangolo
TEST(IntersectionTest, PointTriangle3D_1)
{
	Eigen::Vector3d p(0.0, 0.0, 0.0);
	Eigen::Vector3d a(1.0, 1.0, 1.0); 
	Eigen::Vector3d b(1.0, 2.0, 1.0);
	Eigen::Vector3d c(2.0, 1.0, 1.0);
	auto int_type = fdapde::core::inTri3d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::EXTERN);
}

// punto interno al triangolo
TEST(IntersectionTest, PointTriangle3D_2)
{
	Eigen::Vector3d p(1.1, 1.1, 1.0);
	Eigen::Vector3d a(1.0, 1.0, 1.0); 
	Eigen::Vector3d b(1.0, 2.0, 1.0);
	Eigen::Vector3d c(2.0, 1.0, 1.0);
	auto int_type = fdapde::core::inTri3d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::INTERN);
}

// punto sul vertice
TEST(IntersectionTest, PointTriangle3D_3)
{
	Eigen::Vector3d p(1.0, 2.0, 1.0);
	Eigen::Vector3d a(1.0, 1.0, 1.0); 
	Eigen::Vector3d b(1.0, 2.0, 1.0);
	Eigen::Vector3d c(2.0, 1.0, 1.0);
	auto int_type = fdapde::core::inTri3d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::ONVERTEX);
}

// punto sul lato
TEST(IntersectionTest, PointTriangle3D_4)
{
	Eigen::Vector3d p(1.0, 1.5, 1.0);
	Eigen::Vector3d a(1.0, 1.0, 1.0); 
	Eigen::Vector3d b(1.0, 2.0, 1.0);
	Eigen::Vector3d c(2.0, 1.0, 1.0);
	auto int_type = fdapde::core::inTri3d(p, a, b, c);
	EXPECT_TRUE(int_type == fdapde::core::Point2Tri::ONEDGE);
}


// INTERSEZIONI TRA TRIANGOLI NELLO SPAZIO

// Triangoli che non si intersecano
TEST(IntersectionTest, Triangles_1)
{
	Eigen::Vector3d p1(0.0, 0.0, 0.0);
	Eigen::Vector3d p2(1.0, 1.0, 1.0);
	Eigen::Vector3d p3(2.0, 2.0, 2.0);
	Eigen::Vector3d p4(10.0, 10.0, 10.0);
	Eigen::Vector3d p5(11.0, 11.0, 11.0);
	Eigen::Vector3d p6(12.0, 12.0, 12.0);
	Element<2, 3> el1(0, {1, 2, 3}, {p1, p2, p3}, {}, false);
	Element<2, 3> el2(0, {1, 2, 3}, {p4, p5, p6}, {}, false);
	EXPECT_FALSE(el1.intersection(el2));
}

// Triangoli che si intersecano
TEST(IntersectionTest, Triangles_2)
{
	Eigen::Vector3d p1(0.0, 0.0, 0.0);
	Eigen::Vector3d p2(2.0, 2.0, 0.0);
	Eigen::Vector3d p3(2.0, 2.0, 4.0);
	Eigen::Vector3d p4(1.0, 0.0, 0.1);
	Eigen::Vector3d p5(0.0, 1.0, 0.1);
	Eigen::Vector3d p6(0.0, 0.0, 3.0);
	Element<2, 3> el1(0, {1, 2, 3}, {p1, p2, p3}, {}, false);
	Element<2, 3> el2(0, {1, 2, 3}, {p4, p5, p6}, {}, false);
	EXPECT_TRUE(el1.intersection(el2));
}
