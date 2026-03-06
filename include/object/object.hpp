/**
 * @file object.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_OBJECT_HPP
#define SBODY_OBJECT_HPP

#include <array>
#include <functional>
#include <memory>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "metric/metric.hpp"

template <typename T, std::size_t N>
using bounded_vector = boost::numeric::ublas::bounded_vector<T, N>;
template <typename T, std::size_t M, std::size_t N>
using bounded_matrix = boost::numeric::ublas::bounded_matrix<T, M, N>;

namespace sbody {
	/**
	 * @brief
	 *
	 */
	template <typename T>
	class Object {
	  protected:
		const Metric<T> &metric_;

	  public:
		/**
		 * @brief Construct a new Object object
		 *
		 * @param metric
		 */
		Object(const Metric<T> &metric) : metric_(metric) {}
		/**
		 * @brief
		 *
		 * @param current
		 * @param last
		 * @return bool
		 */
		virtual bool Hit(const bounded_vector<T, 8> &current, const bounded_vector<T, 8> &last) const = 0;

		/**
		 * @brief
		 * \f[z=\frac{E_\mathrm{obj}}{E_\mathrm{obs}}-1\f]
		 * where
		 * \f[E_\mathrm{obj}=-g_{\mu\nu}u_\mathrm{ph}^\mu u_\mathrm{obj}^\nu\f]
		 * For the distant observer, the equation is the same. With the tetrad velocity \f$u_\mathrm{obs}=(1,0,0,0)\f$ and metric
		 * \f[
		 * 	g_\mathrm{obs}=
		 * 	\begin{pmatrix}
		 * 		-1 & 0 & 0 & 0 \\
		 * 		0  & 1 & 0 & 0 \\
		 * 		0  & 0 & 1 & 0 \\
		 * 		0  & 0 & 0 & 1
		 * 	\end{pmatrix},
		 * \f]
		 * the energy of the photon can be simplified to
		 * \f[
		 * 	E_\mathrm{obs}=u_\mathrm{ph}^0.
		 * \f]
		 *
		 * @param photon 8 dimensional information of photon
		 * @param photon_time
		 * @return T
		 */
		virtual T Redshift(const boost::numeric::ublas::bounded_vector<T, 8> &photon, TimeSystem photon_time) const = 0;
	};
} // namespace sbody

#endif
