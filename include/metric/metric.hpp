/**
 * @file metric.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <cerrno>
#include <cmath>
#include <functional>
#include <numeric>
#include <string>

// #include <boost/math/tools/quartic_roots.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "infra/infra_context.hpp"
#include "infra/status.hpp"
#include "math/scalar.hpp"

namespace sbody {
/**
 * @brief The time system used in the record of the object.
 *
 * - T: the time step of the integrator is the observer's time \f$t\f$, and the first element of the object vector is its
 * intrinsic time \f$\tau\f$.
 * - TAU: the time step of the integrator is the intrinsic time of the object, and the first element is the observer's time.
 *
 */
enum class TimeSystem { T, TAU, TIME_SYSTEM_NUM };

/**
 * @brief The dynamical system to describe the motion of the object.
 *
 * - LAGRANGIAN: the velocity part of the object vector is \f$\mathrm d{x^\mu}/\mathrm dt\f$ or \f$\mathrm d{x^\mu}/\mathrm
 * d\tau\f$.
 * - HAMILTONIAN: the velocity is recorded as \f$p_\mu\f$.
 *
 */
enum class DynamicalSystem { LAGRANGIAN, HAMILTONIAN, DYNAMICAL_SYSTEM_NUM };

/**
 * @brief The mode of the motion of the object.
 *
 * - GEODESIC: the motion is described by the geodesic equation.
 * - CIRCULAR: the orbit has fixed radius.
 * - CONICAL: the motion is described by the conical orbit equation.
 *
 */
enum class MotionMode {
    GEODESIC,
    CIRCULAR,
    CONICAL,
    CYLINDRICAL,
    CONICAL_FIXED_ENERGY,
    CYLINDRICAL_FIXED_ENERGY,
    MOTION_MODE_NUM
};

/**
 * @brief The base class of all kinds of metrics.
 *
 */
template <typename T>
class Metric {
  public:
    /**
     * @brief Return the name of the metric.
     *
     * @return name of the metric.
     */
    virtual std::string name() const = 0;

    /**
     * @brief Calculate the metric tensor at `position`, stored in `metric`.
     *
     * @param position 4 dimensional vector
     * @param metric matrix with size 4×4
     * @return status
     */
    virtual int metricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position,
                             boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const = 0;

    virtual T dotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2,
                           T y3) const = 0;
    virtual T dotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1,
                           T y2, T y3) const = 0;

    /**
     * @brief Dot product of vector `x` and `y` at `position`. \f$g_{\mu\nu}x^\mu y^\nu\f$
     *
     * @param position 4 dimensional vector, position to calcuate the dot product of `x` and `y`.
     * @param x 4 dimensional vector
     * @param y 4 dimensional vector
     * @param include_time include time related terms if true.
     * @return result
     */
    template <typename VecExp1, typename VecExp2>
    T dotProduct(const boost::numeric::ublas::bounded_vector<T, 8> &position, const VecExp1 &x, const VecExp2 &y,
                 bool include_time) const {
        if (include_time)
            return dotProduct4D(position, x(0), x(1), x(2), x(3), y(0), y(1), y(2), y(3));
        return dotProduct3D(position, x(1), x(2), x(3), y(1), y(2), y(3));
    }

    /**
     * @brief Calculate the square of the distance between `x` and `y` at `x`. \f$g_{\mu\nu}(x^\mu-y^\mu)(x^\nu-y^\nu)\f$
     *
     * @param x 8 dimensional vector
     * @param y 8 dimensional vector
     * @return result
     */
    virtual T distanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x,
                             const boost::numeric::ublas::bounded_vector<T, 8> &y, bool include_time) const = 0;

    /**
     * @brief Calculate the local inertial frame coordinate of the object at `position`, stored in `coordinate`.
     *
     * @param position 8 dimensional vector
     * @param time time systme of the object
     * @param coordinate matrix with size 4×4
     * @return status
     */
    int localInertialFrame(const boost::numeric::ublas::bounded_vector<T, 8> &position, TimeSystem time,
                           boost::numeric::ublas::bounded_matrix<T, 4, 4> &coordinate) const {
        namespace ublas = boost::numeric::ublas;
        coordinate = ublas::identity_matrix<T>(4);
        ublas::bounded_matrix<T, 4, 4> metric, product = ublas::identity_matrix<T>(4), product_LU;
        metricTensor(position, metric);
        ublas::permutation_matrix<std::size_t> permutation(4);
        for (int i = 0; i < 4; ++i) {
            auto coordinate_row = ublas::matrix_row<ublas::bounded_matrix<T, 4, 4>>(coordinate, i);
            if (i == 0) {
                if (time == TimeSystem::T)
                    std::copy(position.begin() + 5, position.end(), coordinate_row.begin() + 1);
                else
                    std::copy(position.begin() + 4, position.end(), coordinate_row.begin());
            } else {
                product_LU = product;
                std::iota(permutation.begin(), permutation.end(), 0);
                if (int status = ublas::lu_factorize(product_LU, permutation); status != StatusCode::SUCCESS)
                    return status;
                ublas::lu_substitute(product_LU, permutation, coordinate_row);
            }
            coordinate_row /= std::sqrt(std::abs(dotProduct(position, coordinate_row, coordinate_row, true)));
            ublas::matrix_row<ublas::bounded_matrix<T, 4, 4>>(product, i) = ublas::prod(metric, coordinate_row);
        }
        return std::isnan(coordinate(3, 0)) ? StatusCode::NUMERIC_ERROR : StatusCode::SUCCESS;
    }

    /**
     * @brief Convert the coordinate system from Lagrangian to Hamiltonian.
     *
     * @param y 8 dimensional vector
     * @return status
     */
    virtual int lagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const = 0;

    /**
     * @brief Convert the coordinate system from Hamiltonian to Lagrangian.
     *
     * @param y 8 dimensional vector
     * @return status
     */
    virtual int hamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const = 0;

    int initializePhoton(boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time, T alpha, T beta, T r, T r2,
                         T theta, T sin_theta) const {
        photon[0] = 0.;
        photon[1] = r;
        photon[4] = 1.;
        photon[5] = 1.;
        photon_time = r;
        if (sin_theta < boost::math::tools::root_epsilon<T>()) {
            const T k = std::hypot(alpha, beta);
            if (theta < boost::math::constants::half_pi<T>()) {
                photon[2] = 1e-15;
                photon[3] = atan2(alpha, -beta);
                photon[6] = -k / r2;
            } else {
                photon[2] = boost::math::constants::pi<T>() - 1e-15;
                photon[3] = atan2(alpha, beta);
                photon[6] = k / r2;
            }
            photon[7] = 0.;
        } else {
            photon[2] = theta;
            photon[3] = 0., photon[6] = beta / r2;
            photon[7] = -alpha / (r2 * sin_theta);
        }
        return normalizeNullGeodesic(photon, 1.);
    }

    int angularMomentumCarterConstantToAlphaBeta(T l, T q2, T cos_theta, T sin_theta, T &alpha, T &abs_beta) const {
        if (sin_theta == 0.)
            return StatusCode::FAILURE;
        alpha = -l / sin_theta;
        const T beta2 = q2 - power2(cos_theta * alpha);
        if (beta2 < 0.)
            return StatusCode::DOMAIN_ERROR;
        abs_beta = std::sqrt(beta2);
        return StatusCode::SUCCESS;
    }

    /**
     * @brief Trace the photon from the observer to the target, using elliptic integrals.
     *
     * @param r_observer radius of the observer, \f$r_\text{obs}\f$.
     * @param theta_observer theta of the observer, \f$\theta_\text{obs}\f$.
     * @param sin_theta_observer sin(theta_observer), \f$\sin\theta_\text{obs}\f$.
     * @param cos_theta_observer cos(theta_observer), \f$\cos\theta_\text{obs}\f$.
     * @param r_object radius of the target, \f$r_\text{tar}\f$.
     * @param theta_object theta of the target, \f$\theta_\text{tar}\f$.
     * @param phi_object phi of the target, \f$\phi_\text{tar}\f$.
     * @param alpha x position of the target in the observer's view.
     * @param beta y position of the target in the observer's view.
     * @param photon 8 dimensional vector, position and the velocity of the photon traced to the target.
     * @param photon_time look back time of the photon.
     * @param absolute_accuracy absolute accuracy for the trace.
     * @param relative_accuracy relative accuracy for the trace.
     * @param infra_context context for log and output.
     * @return status
     */
    virtual int fastTrace(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer,
                          const T r_object, const T theta_object, const T phi_object, T &alpha, T &beta,
                          boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time, T absolute_accuracy,
                          T relative_accuracy, const InfraContext &infra_context) const = 0;

    virtual int fastShadow(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer,
                           const T alpha, const T beta, const T r_min) const {
        return StatusCode::NOT_IMPLEMENTED;
    }

    /**
     * @brief Calculate the energy of the object.
     *
     * @param y 8 dimensional vector
     * @param time time system of the object
     * @param dynamics dynamical system of the object
     * @return result
     */
    virtual T energy(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const = 0;
    /**
     * @brief Calculate the angular momentum of the object.
     *
     * @param y 8 dimensional vector
     * @param time time system of the object
     * @param dynamics dynamical system of the object
     * @return result
     */
    virtual T angularMomentum(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time,
                              DynamicalSystem dynamics) const = 0;

    /**
     * @brief Calculate the Carter constant of the object.
     *
     * @param y 8 dimensional vector
     * @param time time system of the object
     * @param dynamics dynamical system of the object
     * @return result
     */
    virtual T carterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, const T mu2, TimeSystem time,
                             DynamicalSystem dynamics) const = 0;

    /**
     * @brief Calculate the redshift of the object, \f$1+z\f$
     *
     * @param y 8 dimensional vector of the object
     * @param photon photon traced to the object
     * @param object_time time system of the object
     * @param photon_time time system of the photon
     * @return result
     */
    virtual T redshift(const boost::numeric::ublas::bounded_vector<T, 8> &y,
                       const boost::numeric::ublas::bounded_vector<T, 8> &photon_position, TimeSystem object_time,
                       TimeSystem photon_time) const {
        namespace ublas = boost::numeric::ublas;
        ublas::bounded_vector<T, 4> u, v;
        u(0) = 1.;
        std::copy(y.begin() + 5, y.end(), u.begin() + 1);
        v(0) = 1.;
        std::copy(photon_position.begin() + 5, photon_position.end(), v.begin() + 1);
        if (object_time == TimeSystem::T) {
            if (photon_time == TimeSystem::T)
                return -dotProduct(y, u, v, true) / (y[4] * photon_position[4]);
            return -dotProduct(y, u,
                               ublas::vector_range<const ublas::bounded_vector<T, 8>>(photon_position, ublas::range(4, 8)),
                               true) /
                   y[4];
        }
        if (photon_time == TimeSystem::T)
            return -dotProduct(y, ublas::vector_range<const ublas::bounded_vector<T, 8>>(y, ublas::range(4, 8)), v, true) /
                   photon_position[4];
        return -dotProduct(y, ublas::vector_range<const ublas::bounded_vector<T, 8>>(y, ublas::range(4, 8)),
                           ublas::vector_range<const ublas::bounded_vector<T, 8>>(photon_position, ublas::range(4, 8)), true);
    }

    /**
     * @brief Normalize the timelike geodesic.
     *
     * @param y 8 dimensional vector
     * @return status
     */
    virtual int normalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const = 0;

    /**
     * @brief Normalize the null geodesic.
     *
     * @param y 8 dimensional vector
     * @return status
     */
    virtual int normalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const = 0;

    /**
     * @brief Get the integration system to calculate the motion of the object.
     *
     * @param time time system of the object
     * @param dynamics dynamical system of the object
     * @param motion motion mode of the object
     * @return integrator system
     */
    virtual int getIntegrationSystem(std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &,
                                                        boost::numeric::ublas::bounded_vector<T, 8> &, T)> &system,
                                     TimeSystem time, DynamicalSystem dynamics,
                                     MotionMode motion = MotionMode::GEODESIC) const = 0;

    virtual T innermostStableCircularOrbitRadius(bool prograde) const { return 6.; }
};
} // namespace sbody
