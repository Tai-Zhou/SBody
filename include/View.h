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
	class view {
	  protected:
		const double r;
		const double theta;
		const double sinto;
		const double costo;
		const double tFinal;
		std::unique_ptr<file> output;

	  public:
		/**
		 * @brief Construct a new view object
		 *
		 * @param r
		 * @param theta
		 * @param fileName
		 */
		view(double r, double theta, std::string fileName);
		/**
		 * @brief
		 *
		 * @return int
		 */
		int rayInit();
		/**
		 * @brief
		 *
		 * @return int
		 */
		int rayTrace();
		/**
		 * @brief
		 *
		 * @param s
		 * @param rayNO
		 * @return int
		 */
		int traceStar(Object::star &s, int rayNO);
		/**
		 * @brief
		 *
		 * @param n
		 * @return int
		 */
		int shadow(int n);
	};
	class camera : public view {
	  protected:
		const size_t pixel;
		const double halfAngle;
		std::vector<std::array<double, 10>> initials;
		std::vector<std::vector<double>> screen;

	  public:
		/**
		 * @brief Construct a new camera object
		 *
		 * @param pixel
		 * @param halfAngle
		 * @param r
		 * @param theta
		 * @param fileName
		 */
		camera(size_t pixel, double halfAngle, double r, double theta, std::string fileName);
		/**
		 * @brief
		 *
		 * @return int
		 */
		int traceStar();
		/**
		 * @brief
		 *
		 * @return int
		 */
		int lens();
		/**
		 * @brief
		 *
		 * @return int
		 */
		int save();
	};
} // namespace SBody

#endif
