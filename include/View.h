/**
 * @file View.h
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_VIEW_H
#define SBODY_VIEW_H

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "IO.h"
#include "Metric.h"
#include "Object.h"

namespace SBody {
	class View {
	  protected:
		std::shared_ptr<Metric> metric_;
		/// Distance to the center black hole
		const double r_;
		/// Square of the distance to the center black hole
		const double r2_;
		/// Angle between the observer and the \f$z\f$ axis, \f$\theta\f$.
		const double theta_;
		/// \f$\sin\theta\f$
		const double sin_theta_;
		/// \f$\cos\theta\f$
		const double cos_theta_;
		/// Rotational angle of the coordiante of the view, \f$\iota\f$.
		const double iota_;
		/// \f$\sin\iota\f$
		const double sin_iota_;
		/// \f$\cos\iota\f$
		const double cos_iota_;
		/// Position and velocity of the observer.
		double position_[8];
		/// time limit for the integration of photons.
		const double t_final_;
		std::unique_ptr<indicators::BlockProgressBar> bar_;

	  public:
		/**
		 * @brief Constructor
		 *
		 * @param r Distance to the black hole
		 * @param theta Angle between the observer and the \f$z\f$ axis, \f$\theta\f$.
		 * @param iota Rotational angle of the coordiante of the view, \f$\iota\f$.
		 */
		View(std::shared_ptr<Metric> metric, double r, double theta, double iota, double v_alpha = 0.0, double v_beta = 0.0);

		/**
		 * @brief Initialize the position and velocity of a trace back photon.
		 *
		 * @param photon 9 dimensional vector, photon[8] is used to store the look back time.
		 * @param alpha x position of the target in the observer's view.
		 * @param beta y position of the target in the observer's view.
		 * @return status
		 */
		int
		InitializePhoton(double photon[], double alpha, double beta);

		/**
		 * @brief
		 *
		 * @param position 8 dimensional vector
		 * @param object_time time system of the object
		 * @param record position to save the observational values
		 * @param calculate_magnification switch of the magnification calculation
		 * @param fast_trace switch of the fast trace
		 * @return status
		 */
		int Trace(const double position[], TimeSystem object_time, double record[], bool calculate_magnification, bool fast_trace = true);

		static int TraceToPlane(const gsl_vector *alpha_beta, void *params, gsl_vector *delta_apparent_alpha_beta);

		int PhotonInformation(const double position[], TimeSystem object_time, double record[], const double photon[], double alpha, double beta);

		/**
		 * @brief Calculate the magnification of the object.
		 *
		 * @param position 8 dimensional vector
		 * @param object_time time system of the object
		 * @param magnification position to save the magnification
		 * @param photon 8 dimensional vector of the photon traced to the object
		 * @param redshift redshift of the photon
		 * @return int
		 */
		int Magnification(const double position[], TimeSystem object_time, double &magnification, const double photon[], double redshift);

		/**
		 * @brief Calculate the radius of the black hole shadow at different direction, saved in the file.
		 *
		 * @param file_name file name.
		 * @return status
		 */
		int Shadow(std::string file_name);

		int OmegaTest();
	};
	struct TraceParameters {
		std::shared_ptr<Metric> metric;
		std::shared_ptr<Integrator> integrator;
		double *photon;
		const double r, r2;
		const double theta_obs, sin_theta_obs, cos_theta_obs;
		const double r_obj, x_obj, y_obj, z_obj;
		const double t_final;
		TraceParameters(std::shared_ptr<Metric> metric, std::shared_ptr<Integrator> integrator, double photon[], double r, double r2, double theta_obs, double sin_theta_obs, double cos_theta_obs, double r_obj, double sin_theta_obj, double cos_theta_obj, double sin_phi_obj, double cos_phi_obj, double t_final);
	};
	class Camera : public View {
	  protected:
		const size_t pixel_;
		const double half_angle_;
		std::vector<std::array<double, 10>> initials_;
		std::vector<std::vector<double>> screen_;

	  public:
		/**
		 * @brief Construct a new camera object
		 *
		 * @param pixel
		 * @param half_angle
		 * @param r
		 * @param theta
		 * @param file_name
		 */
		Camera(std::shared_ptr<Metric> metric, size_t pixel, double half_angle, double r, double theta, double iota);

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Trace();

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Lens();

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Save(std::string file_name);
	};
} // namespace SBody

#endif
