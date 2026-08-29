#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/fluid.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// 1.  Onset of emission
// ===========================================================================

/// Field at which a hemispherical meniscus of radius r ceases to be balanced by
/// capillary pressure: eps0 E^2 / 2 = 2 gamma / r  =>  E = 2 sqrt(gamma/(eps0 r)).
/// A crude but useful upper-bound estimate; the meniscus solver gives the real
/// turning point, which is lower because the meniscus elongates first.
Real onset_field_hemisphere(Real r, Real gamma);

/// Classical Taylor/Smith onset voltage for a capillary of bore radius r_c held
/// a distance d from a counter-electrode plane:
///
///     V_on = sqrt( 2 gamma r_c cos(49.3 deg) / eps0 ) * ln(4 d / r_c)
///
/// Derived for a perfectly conducting liquid forming an equilibrium Taylor cone
/// on an infinitely thin tube.  It ignores tube wall thickness, extractor
/// aperture geometry and flow, so expect tens of percent of disagreement with
/// a real emitter -- use it as a sanity check on the meniscus continuation,
/// not as a prediction.
Real onset_voltage_taylor(Real r_c, Real d, Real gamma);

/// Rayleigh charge limit of a droplet of diameter d [C].
Real rayleigh_charge(Real d, Real gamma);

// ===========================================================================
// 2.  Ion evaporation (Iribarne-Thomson)
// ===========================================================================

/// Schottky-type barrier lowering by the local field [J]:
///     G(E) = sqrt(e^3 E / (4 pi eps0))
/// At E = 1 V/nm this is ~1.2 eV, which is why field emission from ionic
/// liquids switches on so abruptly around that field.
Real schottky_lowering(Real E);

/// Emitted ion current density [A/m^2] from a surface at normal field E:
///     j = (kT/h) * sigma_s * exp[ -(dG - G(E)) / (kT) ],   sigma_s = eps0 E
/// The exponential makes this extraordinarily sensitive to dG and to E; see the
/// warning in fluid.hpp.
Real ion_current_density(Real E, const Fluid& f, Real T);

/// Inverse: the field required to reach a given current density [V/m].
/// Solved by bisection; returns 0 if unreachable below 1e11 V/m.
Real field_for_current_density(Real j_target, const Fluid& f, Real T);

/// Field at which ion evaporation becomes energetically free (G(E) = dG).
Real characteristic_evaporation_field(const Fluid& f);

struct IonEmission {
  Real current{0};        ///< total emitted ion current [A]
  Real emitting_area{0};  ///< area contributing 99% of the current [m^2]
  Real peak_j{0};         ///< peak current density [A/m^2]
  Real peak_E{0};         ///< field at the peak [V/m]
  Real mdot{0};           ///< ion mass flow [kg/s]
};

/// Integrate the Iribarne-Thomson rate over the solved surface.  Only elements
/// tagged FreeSurface emit unless `include_wetted_metal` is set, which is the
/// right choice for externally wetted and porous emitters where the whole tip
/// carries a liquid film.
IonEmission integrate_ion_emission(const BemSolver& bem, const Fluid& f, Real T,
                                   bool include_wetted_metal = false);

// ===========================================================================
// 3.  Cone-jet (droplet) mode
// ===========================================================================

/// Empirical correlations for the cone-jet regime.
///
/// CALIBRATION WARNING.  Both prefactors below are experimental fits, not
/// derived constants.  The current scaling I = f(eps_r) sqrt(gamma K Q/eps_r)
/// was established for high-permittivity liquids, where f saturates near 18.
/// Ionic liquids have eps_r ~ 12, well outside that range, and f is known to be
/// smaller there -- but there is no consensus closed form.  Rather than invent
/// one, the code keeps f as an explicit input, defaults it to the
/// high-permittivity value, and warns when it is being extrapolated.  Fit it to
/// your own I(Q) data before trusting absolute currents.
struct ConeJetModel {
  Real f_current{18.0};    ///< prefactor in I = f sqrt(gamma K Q / eps_r)
  Real jet_to_drop{1.89};  ///< droplet diameter / jet diameter (Rayleigh breakup)
  Real d_prefactor{1.0};   ///< prefactor in d_jet = C (Q eps0 eps_r / K)^(1/3)
};

struct ConeJetState {
  Real Q{0};            ///< volumetric flow [m^3/s]
  Real current{0};      ///< emitted current [A]
  Real d_jet{0};        ///< jet diameter [m]
  Real d_droplet{0};    ///< primary droplet diameter [m]
  Real qm{0};           ///< beam charge-to-mass ratio I/(rho Q) [C/kg]
  Real mdot{0};         ///< mass flow [kg/s]
  Real q_droplet{0};    ///< charge per primary droplet [C]
  Real fissility{0};    ///< q_droplet / Rayleigh limit; > ~0.8 means it fissions
  Real q_over_qmin{0};  ///< Q / Q_min -- below ~1 the cone-jet is not stable
  bool extrapolated{false};  ///< true when eps_r < 40, i.e. f_current is a guess
};

ConeJetState cone_jet(const Fluid& f, Real Q, const ConeJetModel& model = {});

// ===========================================================================
// 4.  Beam figures of merit
// ===========================================================================

/// One population in the beam: a mass flow at a given charge-to-mass ratio.
/// Droplets and each ion cluster species are separate entries.
struct Species {
  std::string name{};
  Real mdot{0};  ///< [kg/s]
  Real qm{0};    ///< [C/kg]
};

struct BeamFigures {
  Real current{0};    ///< [A]
  Real mdot{0};       ///< [kg/s]
  Real thrust{0};     ///< [N]
  Real Isp{0};        ///< [s]
  Real beam_power{0}; ///< I * V_accel [W]
  Real mean_qm{0};    ///< [C/kg]
  /// Polydispersive efficiency F^2 / (2 mdot P).  Unity for a monoenergetic,
  /// single-q/m beam; the dominant loss in mixed droplet/ion operation.
  Real eta_polydispersity{1.0};
};

/// Figures of merit for a beam accelerated through V_accel [V].
BeamFigures beam_figures(const std::vector<Species>& species, Real V_accel);

void print_operating_point(std::FILE* out, const Fluid& f, const ConeJetState* cj,
                           const IonEmission* ion, const BeamFigures* fig);

}  // namespace es
