/**
 * @file flux.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_FLUX_HPP
#define SBODY_FLUX_HPP

namespace sbody {
	/**
	 * @brief Calculate the bolometric flux. The result is \f[F_\text{obs}=\frac{L_\text{src}}{(1+z)^2}\frac{\Omega_\text{src}}{\Omega_\text{obs}}.\f] The additional \f$(1+z)^{-1}\f$ comes from the expansion of the frequency that \f$(1+z)^{-1}=\frac{\nu_\text{obs}}{\nu_\text{src}}\f$.
	 *
	 * @param luminosity intrinsic luminosity of the source, \f$L_\text{src}\f$
	 * @param magnification magnification of the source, including the relativistic redshift and beaming, and gravitational lensing. \f$\frac{\Omega_\text{src}}{(1+z)\Omega_\text{obs}}\f$
	 * @param redshift relativistic redshift of the source, \f$1+z\f$
	 * @return result
	 */
	template <typename T>
	T Flux(T luminosity, T magnification, T redshift) {
		return luminosity * magnification / redshift;
	}

	/**
	 * @brief Calculate the flux density. The result is \f[F_{\nu,\text{obs}}=\frac{L_{\nu,\text{src}}}{1+z}\frac{\Omega_\text{src}}{\Omega_\text{obs}}.\f]
	 *
	 * @param spectral_density intrinsic spectral density of the source, \f$L_{\nu,\text{src}}\f$
	 * @param magnification magnification of the source, including the relativistic redshift and beaming, and gravitational lensing. \f$\frac{\Omega_\text{src}}{(1+z)\Omega_\text{obs}}\f$
	 * @return result
	 */
	template <typename T>
	T FluxDensity(T spectral_density, T magnification) {
		return spectral_density * magnification;
	}
} // namespace sbody

#endif
