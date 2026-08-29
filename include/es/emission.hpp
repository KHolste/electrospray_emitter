#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/fluid.hpp"
#include "es/status.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// 1.  Closed-form reference values from the literature
// ===========================================================================
//
// These are NOT results of this model.  They are published closed forms, kept
// so that a computed static fold can be put next to a number of the right order
// of magnitude.  Printing them next to computed values without saying so was
// one of the defects of the prototype, so every caller must label them.

/// Field at which a hemispherical meniscus of radius r ceases to be balanced by
/// capillary pressure: eps0 E^2 / 2 = 2 gamma / r  =>  E = 2 sqrt(gamma/(eps0 r)).
/// A static balance estimate, not a stability criterion.
Real hemisphere_balance_field(Real r, Real gamma);

/// Taylor/Smith closed form for the emission onset voltage of a capillary of
/// bore radius r_c a distance d from a counter-electrode plane:
///
///     V = sqrt( 2 gamma r_c cos(49.3 deg) / eps0 ) * ln(4 d / r_c)
///
/// Source: Smith (1986), IEEE Trans. Ind. Appl. IA-22(3), 527-535.  Derived for
/// a perfectly conducting liquid on an infinitely thin tube; ignores wall
/// thickness, extractor aperture and flow.  It is an EMISSION ONSET in the
/// experimental sense, which is a different quantity from the static fold this
/// code computes -- the two must not be presented as the same thing.
Real literature_onset_voltage_smith(Real r_c, Real d, Real gamma);

/// Rayleigh charge limit of a droplet of diameter d [C].
Real rayleigh_charge(Real d, Real gamma);

// ===========================================================================
// 2.  Ion evaporation (Iribarne-Thomson)
// ===========================================================================

/// Schottky-type barrier lowering by the local field [J]:
///     G(E) = sqrt(e^3 E / (4 pi eps0))
Real schottky_lowering(Real E);

/// Emitted ion current density [A/m^2] from a surface at normal field E:
///     j = (kT/h) * sigma_s * exp[ -(dG - G(E)) / (kT) ],   sigma_s = eps0 E
/// Source: Iribarne & Thomson (1976), J. Chem. Phys. 64(6), 2287-2294.
Real ion_current_density(Real E, const Fluid& f, Real T);

/// Inverse: the field required to reach a given current density [V/m].
Real field_for_current_density(Real j_target, const Fluid& f, Real T);

/// Field at which ion evaporation becomes energetically free (G(E) = dG).
Real characteristic_evaporation_field(const Fluid& f);

/// Throws NotImplementedInThisPhase when the applied polarity would emit
/// anions.  Only the cation series is modelled: masses, solvation energies and
/// the species mix all differ between polarities, and the prototype produced
/// identical currents for both signs because it used |E_n| and cation masses
/// throughout.  Rather than return a number that is wrong for one polarity, the
/// negative case fails closed.
void require_modelled_polarity(const BemSolver& bem);
/// Same check on a raw emitter-to-extractor voltage, so that an application
/// can refuse before it spends minutes on a meniscus continuation -- and so
/// that the reason it reports is the polarity, not some later symptom.
void require_modelled_polarity(Real U_emitter_minus_extractor);

struct IonEmission {
  Real current{0};   ///< integrated emitted ion current [A]
  /// Effective emitting area  A_eff = (int j dA)^2 / int j^2 dA.
  ///
  /// This replaces the "smallest area carrying 99% of the current" of the
  /// prototype, which was a sum over whole elements and therefore quantised:
  /// it did not converge smoothly under mesh refinement (measured 3.198 /
  /// 3.150 / 3.103 / 3.150 e-10 m^2 at 61/81/121/161 nodes).  A_eff is a ratio
  /// of two integrals of j, converges like any other integrated quantity, and
  /// equals the true area for a uniform current density.
  Real effective_area{0};
  Real peak_j{0};    ///< peak current density [A/m^2]
  Real peak_E{0};    ///< field at the peak [V/m]
  Real mdot{0};      ///< ion mass flow implied by the mean cluster q/m [kg/s]
};

/// Integrate the Iribarne-Thomson rate over the solved surface.
///
/// IMPORTANT -- this is a DIAGNOSTIC ESTIMATE, not a predicted current.  It
/// applies the evaporation rate to the field of a static, perfectly conducting,
/// non-emitting meniscus.  In the pure ionic regime the flow is
/// viscosity-dominated and the current is controlled by the finite conductivity
/// (Higuera 2008, Phys. Rev. E 77, 026308), so the perfect-conductor field is
/// not the field that would actually be present.  The self-consistent model is
/// due in phase P5.
IonEmission integrate_ion_emission(const BemSolver& bem, const Fluid& f, Real T,
                                   bool include_wetted_metal = false);

// ===========================================================================
// 3.  Cone-jet (droplet) mode -- empirical, not coupled to anything
// ===========================================================================

/// CALIBRATION WARNING.  Both prefactors are experimental fits.  The current
/// scaling I = f(eps_r) sqrt(gamma K Q/eps_r) was established for
/// high-permittivity liquids where f saturates near 18; ionic liquids have
/// eps_r ~ 12, well outside that range.  Source: Fernandez de la Mora &
/// Loscertales (1994), J. Fluid Mech. 260, 155-184.
struct ConeJetModel {
  Real f_current{18.0};
  Real jet_to_drop{1.89};
  Real d_prefactor{1.0};
};

struct ConeJetState {
  Real Q{0};
  Real current{0};
  Real d_jet{0};
  Real d_droplet{0};
  Real qm{0};
  Real mdot{0};
  Real q_droplet{0};
  Real fissility{0};
  Real q_over_qmin{0};
  bool extrapolated{false};
  /// Always true: this is a correlation, not a computed operating point.  It is
  /// carried into every output so that it cannot be mistaken for one.
  static constexpr bool empirical = true;
};

ConeJetState cone_jet(const Fluid& f, Real Q, const ConeJetModel& model = {});

// ===========================================================================
// 4.  Beam figures of merit
// ===========================================================================

struct Species {
  std::string name{};
  Real mdot{0};  ///< [kg/s]
  Real qm{0};    ///< [C/kg]
};

struct BeamFigures {
  Real current{0};
  Real mdot{0};
  Real thrust{0};
  Real Isp{0};
  Real beam_power{0};
  Real mean_qm{0};
  Real eta_polydispersity{1.0};
};

BeamFigures beam_figures(const std::vector<Species>& species, Real V_accel);

/// Print the non-coupled diagnostic estimate.  The wording is deliberate: this
/// is not an operating point and must not be presented as one.
void print_diagnostic_estimate(std::FILE* out, const Fluid& f, const IonEmission* ion,
                               const BeamFigures* fig);

/// Print the cone-jet correlation in its own block, marked empirical.
void print_cone_jet_correlation(std::FILE* out, const Fluid& f, const ConeJetState& cj);

}  // namespace es
