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
		const double r_;
		const double theta_;
		const double sin_theta_observer_;
		const double cos_theta_observer_;
		const double t_final_;
		std::unique_ptr<indicators::BlockProgressBar> bar_;

	  public:
		/**
		 * @brief Construct a new view object
		 *
		 * @param r
		 * @param theta
		 * @param file_name
		 */
		View(std::shared_ptr<Metric> metric, double r, double theta);

		/**
		 * @brief
		 *
		 * @param photon
		 * @param alpha
		 * @param beta
		 * @return int
		 */
		int InitializePhoton(double *photon, double alpha, double beta);

		/**
		 * @brief
		 *
		 * @param star
		 * @param ray_number
		 * @return int
		 */

		int TraceStar(Star &star, double position[], int ray_number, double record[]);

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Shadow(std::string file_name);

		int OmegaTest();
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
		Camera(std::shared_ptr<Metric> metric, size_t pixel, double half_angle, double r, double theta);

		/**
		 * @brief
		 *
		 * @return int
		 */
		int TraceStar();

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
