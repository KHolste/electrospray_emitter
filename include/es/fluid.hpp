#pragma once
#include <string>
#include <vector>

#include "es/types.hpp"

namespace es {

// ---------------------------------------------------------------------------
// Working fluid
// ---------------------------------------------------------------------------
//
// IMPORTANT, please read before trusting a number out of this table.
//
// Ionic-liquid transport properties scatter substantially between sources --
// conductivity in particular is very sensitive to water content and to the
// measurement method, and published values for the same liquid can differ by
// tens of percent.  The entries below are literature-typical room-temperature
// values meant as sane starting points for design sweeps, NOT as a substitute
// for characterising the batch you actually spray.  Every field is public and
// every one can be overridden from the config file.
//
// The ion-evaporation activation energy dG_solvation is the least certain
// quantity here by far: reported values for EMI-BF4 span roughly 1.0-1.4 eV
// depending on how they were extracted, and the emitted current depends on it
// exponentially.  Treat it as the primary fitting parameter when matching
// measured I-V curves, not as a known constant.

struct Fluid {
  std::string name{"custom"};

  // --- reference state ------------------------------------------------------
  Real T_ref{298.15};    ///< temperature at which the values below apply [K]

  // --- bulk properties at T_ref --------------------------------------------
  Real rho{1279.0};      ///< density [kg/m^3]
  Real gamma{0.0452};    ///< surface tension [N/m]
  Real K{1.36};          ///< electrical conductivity [S/m]
  Real mu{0.0371};       ///< dynamic viscosity [Pa s]
  Real eps_r{12.8};      ///< relative permittivity [-]

  // --- temperature coefficients --------------------------------------------
  /// Vogel-Fulcher-Tammann form for conductivity and viscosity:
  ///   K(T)  = K(T_ref)  * exp(-B_K /(T-T0)) / exp(-B_K /(T_ref-T0))
  ///   mu(T) = mu(T_ref) * exp( B_mu/(T-T0)) / exp( B_mu/(T_ref-T0))
  /// The defaults are generic values for imidazolium ionic liquids; they
  /// reproduce the right order of magnitude for the temperature sensitivity but
  /// are not fluid-specific.  Override if you have your own fit.
  Real vft_B_K{700.0};   ///< [K]
  Real vft_B_mu{700.0};  ///< [K]
  Real vft_T0{165.0};    ///< [K]
  Real dgamma_dT{-5.0e-5};  ///< [N/m/K]
  Real drho_dT{-0.6};       ///< [kg/m^3/K]

  // --- species --------------------------------------------------------------
  Real M_cation{111.17e-3 / 6.02214076e23};   ///< mass of the bare cation [kg]
  Real M_anion{86.81e-3 / 6.02214076e23};     ///< mass of the bare anion [kg]
  /// Mean solvation number of the emitted ion cluster.  Field emission from
  /// ionic liquids produces monomers, dimers (n=1) and higher clusters; the
  /// mixture sets the mean q/m and therefore Isp.  Fitted from time-of-flight
  /// data in practice.
  Real mean_solvation_n{0.4};

  // --- ion evaporation ------------------------------------------------------
  Real dG_solvation{1.09 * 1.602176634e-19};  ///< activation energy [J]
  /// Multiplier on the Iribarne-Thomson prefactor, absorbing the unknown
  /// attempt frequency / surface-site density.  1.0 = textbook form.
  Real evap_prefactor{1.0};

  // --- derived --------------------------------------------------------------
  Fluid at_temperature(Real T) const;

  /// Charge relaxation time eps0 eps_r / K [s].  Compare against the
  /// hydrodynamic time scale: if it is much shorter (the usual case for ionic
  /// liquids, ~1e-10 s), the free surface may be treated as an equipotential.
  Real charge_relaxation_time() const;

  /// Electrohydrodynamic length scale r* = (gamma eps0^2 eps_r^2 / (rho K^2))^(1/3)
  /// -- the cone-jet radius at the minimum stable flow rate.  A few nanometres
  /// for ionic liquids, which is why their jets break into such small droplets.
  Real ehd_length() const;

  /// Minimum flow rate for a stable cone-jet, Q_min ~ gamma eps0 eps_r/(rho K).
  /// An order-of-magnitude criterion, not a sharp threshold.
  Real q_min() const;

  /// q/m of a bare cation and of the mean emitted cluster [C/kg].
  Real qm_bare_cation() const;
  Real qm_cluster() const;

  void print(std::FILE* f) const;
};

/// Built-in fluids.  Names are case-insensitive; underscores and hyphens are
/// interchangeable.  Throws std::runtime_error if unknown.
Fluid fluid_by_name(const std::string& name);
std::vector<std::string> fluid_names();

}  // namespace es
