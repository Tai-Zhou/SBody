/**
 * @file TestView.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-02-12
 *
 * @copyright Copyright (c) 2026
 *
 */

#define BOOST_TEST_MODULE TestView

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/test/included/unit_test.hpp>

#include "Utility.hpp"

using namespace std;
using namespace SBody;

#define PRECISION long double

BOOST_AUTO_TEST_CASE(TestPolySolveQuadratic) {
	namespace ublas = boost::numeric::ublas;
	ublas::bounded_matrix<double, 4, 4> gmunu;
	gmunu(0, 0) = -1.;
	gmunu(1, 1) = gmunu(2, 2) = gmunu(3, 3) = 1.0;
	ublas::bounded_matrix<double, 4, 4> coordinate;
	for (int i = 0; i < 4; ++i)
		coordinate(i, i) = 1.0;
	ublas::bounded_vector<double, 4> photon_transform;
	photon_transform(0) = 1.;
	photon_transform(1) = photon_transform(1) = photon_transform(1) = 0.5;
	ublas::bounded_vector<double, 4> photon_in_object_frame_cartesian, photon_in_object_frame_spherical;
	ublas::bounded_matrix<double, 4, 4> coordinate_gmunu = ublas::prod(gmunu, coordinate); // object local inertial coordinate frame
	photon_in_object_frame_cartesian = ublas::prod(coordinate_gmunu, photon_transform);
	photon_in_object_frame_cartesian /= photon_in_object_frame_cartesian(0);
	CartesianToSpherical(photon_in_object_frame_cartesian, photon_in_object_frame_spherical);
}
