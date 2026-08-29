#pragma once
#include "es/types.hpp"

namespace es::constants {

// CODATA 2018 / SI-2019 exact values where applicable.
inline constexpr Real e         = 1.602176634e-19;   ///< elementary charge [C]
inline constexpr Real kB        = 1.380649e-23;      ///< Boltzmann constant [J/K]
inline constexpr Real h_planck  = 6.62607015e-34;    ///< Planck constant [J s]
inline constexpr Real N_A       = 6.02214076e23;     ///< Avogadro constant [1/mol]
inline constexpr Real eps0      = 8.8541878128e-12;  ///< vacuum permittivity [F/m]
inline constexpr Real amu       = 1.66053906660e-27; ///< atomic mass unit [kg]
inline constexpr Real g0        = 9.80665;           ///< standard gravity [m/s^2]
inline constexpr Real pi        = 3.14159265358979323846;

/// Taylor's equilibrium cone half-angle [rad] (49.2915...deg).
/// Root of P_{1/2}(cos theta) = 0, the only angle for which the electrostatic
/// pull ~ r^{-1/2} balances the capillary pressure ~ 1/r on a conical surface.
inline constexpr Real taylor_angle = 0.8603358;

/// eV -> J
inline constexpr Real eV = e;

}  // namespace es::constants
