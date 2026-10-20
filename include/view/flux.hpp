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

#pragma once

namespace sbody {
/**
 * @brief Calculate the bolometric flux. The result is
 * \f[F_\text{obs}=\frac{L_\text{src}}{(1+z)^2}\frac{\Omega_\text{src}}{\Omega_\text{obs}}.\f] The additional \f$(1+z)^{-1}\f$
 * comes from the expansion of the frequency that \f$(1+z)^{-1}=\frac{\nu_\text{obs}}{\nu_\text{src}}\f$.
 *
 * @param luminosity intrinsic luminosity of the source, \f$L_\text{src}\f$
 * @param magnification magnification of the source, including the relativistic redshift and beaming, and gravitational lensing.
 * \f$\frac{\Omega_\text{src}}{(1+z)\Omega_\text{obs}}\f$
 * @param doppler_factor Doppler factor of the source, \f$g=1/(1+z)\f$
 * @return result
 */
template <typename T>
T flux(T luminosity, T magnification, T doppler_factor) {
    // [OBS-03] formulas/obs/obs.py
    return luminosity * magnification * doppler_factor;
}

/**
 * @brief Calculate the flux density. The result is
 * \f[F_{\nu,\text{obs}}=\frac{L_{\nu,\text{src}}}{1+z}\frac{\Omega_\text{src}}{\Omega_\text{obs}}.\f]
 *
 * @param spectral_density frequency-shifted spectral density
 * \f$L_\nu(\nu_\text{obs}/g)=L_\nu g^\alpha\f$ for \f$L_\nu\propto\nu^{-\alpha}\f$
 * @param magnification magnification of the source, including the relativistic redshift and beaming, and gravitational lensing.
 * \f$\frac{\Omega_\text{src}\,g}{\Omega_\text{obs}}\f$; the solid-angle ratio already contains
 * the beaming factor \f$g^2\f$, so this is \f$g^3\f$ times the geometric lensing ratio.
 * @return result
 */
template <typename T>
T fluxDensity(T spectral_density, T magnification) {
    // [OBS-03] formulas/obs/obs.py
    return spectral_density * magnification;
}
} // namespace sbody
