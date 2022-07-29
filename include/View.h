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
#include "Object.h"

namespace SBody {
	class View {
	  protected:
		const double r_;
		const double theta_;
		const double sin_theta_observer_;
		const double cos_theta_observer_;
		const double t_final_;
		std::unique_ptr<File> output_;

	  public:
		/**
		 * @brief Construct a new view object
		 *
		 * @param r
		 * @param theta
		 * @param file_name
		 */
		View(double r, double theta, std::string file_name);
		/**
		 * @brief
		 *
		 * @return int
		 */
		int RayInit();
		/**
		 * @brief
		 *
		 * @return int
		 */
		int RayTrace();
		/**
		 * @brief
		 *
		 * @param star
		 * @param ray_number
		 * @return int
		 */
		int TraceStar(Star &star, int ray_number);
		/**
		 * @brief
		 *
		 * @param n
		 * @return int
		 */
		int Shadow(int n);
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
		Camera(size_t pixel, double half_angle, double r, double theta, std::string file_name);
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
		int Save();
	};
} // namespace SBody

#endif
