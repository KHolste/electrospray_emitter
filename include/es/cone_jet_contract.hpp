#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/material_data.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P8 -- the cone-jet / droplet mode: its own contract, its own status
// ===========================================================================
//
// THIS MODE IS STRICTLY SEPARATE from pure-ion emission (P5) and from the
// static meniscus (P3a/P3b).  Mixing them was one of the prototype's defects:
// the cone-jet correlation produced a droplet current that was printed next to
// an ion current as if the two were alternatives of one model.  They are not.
// A cone-jet has a JET -- a slender charged filament that breaks into drops --
// and a static pinned meniscus has none; the two are different states of the
// same liquid and cannot be computed from one another.
//
// STATUS: blocked.  A cone-jet model needs a two-phase free-surface flow with
// finite conductivity, and P3/P4 established that none of that exists here.
// What this file provides instead is the CONTRACT: the equations that would
// have to be solved, the inputs they would need, the dimensionless numbers that
// ARE unambiguously defined, and a closed failure while the sub-models are
// missing.
//
// WHY THE EMPIRICAL SCALING IS NOT ADOPTED EITHER
//
// The standard scaling of Ganan-Calvo (Phys. Rev. Lett. 79, 217 (1997)), with
// its erratum (Phys. Rev. Lett. 85, 4193 (2000)), gives the emitted current and
// the jet radius as powers of the flow rate.  The task allows adopting such a
// correlation only after the ORIGINAL equation, the ERRATUM, the PREFACTOR and
// the VALIDITY RANGE have each been checked at the source.
//
// None of the four was checked in this run: neither the 1997 paper nor the 2000
// erratum was reachable in full text, and the numbers that a search returns are
// snippet values, which this project does not accept as sources.  The
// correlation is therefore NOT implemented, and the old ConeJetModel in
// emission.hpp stays as it is -- unused by anything here and never a prediction.
//
// WHAT IS PROVIDED: dimensionless numbers that are DEFINITIONS, so that they
// need no literature source.  Each one is derived in its own comment from the
// balance it expresses.  They diagnose which physics would dominate; they do
// NOT predict a regime, and the figure that shows them says so.

// ---------------------------------------------------------------------------

enum class ConeJetStatus {
  /// The default and the shipped state.
  Blocked = 0,
  /// A material property the diagnosis needs is MissingMaterialData.
  MissingMaterialData,
  /// A flow rate or a voltage was not given.
  MissingOperatingPoint,
  /// The dimensionless diagnosis could be evaluated.  It is a DIAGNOSIS.
  DiagnosisAvailable,
};
const char* to_string(ConeJetStatus s);
const char* explain(ConeJetStatus s);

/// Which sub-model a cone-jet computation would need, and whether this project
/// has it.  The list is the contract; the flags are measured against the other
/// phases rather than asserted.
struct ConeJetRequirement {
  const char* name{""};
  const char* why{""};
  bool available{false};
  const char* provided_by{""};
};

/// The full list, in the order a solver would need them.
std::vector<ConeJetRequirement> cone_jet_requirements();

// ---------------------------------------------------------------------------
// Dimensionless numbers that are DEFINITIONS
// ---------------------------------------------------------------------------
//
// Every one below is derived from a balance of two terms and needs no source.
// The derivations are in the comments; the tests check the algebra and the
// dimensions, not a literature value.

/// Charge relaxation time  tau_e = eps0 eps_r / K  [s].
/// From  d rho/dt + (K/eps) rho = 0  (see P3, docs/14).
Real charge_relaxation_time_cj(Real eps_r, Real K);

/// Capillary-inertial time of a structure of size a:  t_c = sqrt(rho a^3 / gamma) [s].
/// From balancing inertia rho a/t^2 against the capillary pressure gamma/a^2
/// over the length a.
Real capillary_inertial_time(Real rho, Real gamma, Real a);

/// Visco-capillary time  t_v = mu a / gamma  [s].
/// From balancing the viscous stress mu/t against gamma/a.
Real visco_capillary_time(Real mu, Real gamma, Real a);

/// Electrohydrodynamic length, DERIVED HERE and not quoted:
/// set the charge relaxation time equal to the capillary-inertial time of a
/// structure of size r,
///     eps0 eps_r / K = sqrt(rho r^3 / gamma)   =>
///     r* = ( gamma eps0^2 eps_r^2 / (rho K^2) )^(1/3)  [m].
/// It is the size at which a structure can no longer outrun charge relaxation.
/// It is NOT "the jet radius" and it is NOT the Q0 of any published scaling --
/// those carry prefactors this project has not checked.
Real electrohydrodynamic_length(Real gamma, Real eps_r, Real rho, Real K);

/// Ohnesorge number  Oh = mu / sqrt(rho gamma a)  [-].
/// The ratio of the visco-capillary to the capillary-inertial time.
Real ohnesorge(Real mu, Real rho, Real gamma, Real a);

/// Electric Bond number  Bo_E = eps0 E^2 a / (2 gamma)  [-].
/// The ratio of the Maxwell pressure eps0 E^2 / 2 to the capillary pressure
/// gamma / a -- the SAME two terms P3b balances, so this number is the one
/// place where the cone-jet diagnosis and the static equilibrium meet.
Real electric_bond(Real E, Real a, Real gamma);

/// Reynolds number of the feed flow  Re = 2 rho Q / (pi R mu)  [-].
Real feed_reynolds(Real rho, Real Q, Real R, Real mu);

/// Capillary number of the feed flow  Ca = mu u / gamma with u = Q/(pi R^2).
Real feed_capillary(Real mu, Real Q, Real R, Real gamma);

// ---------------------------------------------------------------------------

struct ConeJetDiagnosis {
  ConeJetStatus status{ConeJetStatus::Blocked};
  std::string message;

  // inputs, echoed with their provenance status
  Real gamma{0}, rho{0}, mu{0}, K{0}, eps_r{0};
  MaterialDataStatus gamma_status{MaterialDataStatus::MissingMaterialData};
  MaterialDataStatus rho_status{MaterialDataStatus::MissingMaterialData};
  MaterialDataStatus mu_status{MaterialDataStatus::MissingMaterialData};
  MaterialDataStatus K_status{MaterialDataStatus::MissingMaterialData};
  MaterialDataStatus eps_r_status{MaterialDataStatus::MissingMaterialData};

  Real a{0}, Q{0}, R_channel{0}, E_surface{0};

  // the derived numbers; NaN where an input is missing
  Real tau_e{0}, t_capillary{0}, t_viscous{0}, t_residence{0};
  Real r_star{0}, Oh{0}, Bo_E{0}, Re{0}, Ca{0};

  void print(std::FILE* out) const;
};

/// Evaluate the diagnosis from a MATERIAL DATASET, so that a missing property
/// fails closed instead of being replaced.  `E_surface` is the surface field a
/// caller measured elsewhere (P3b); it enters only the electric Bond number.
ConeJetDiagnosis diagnose_cone_jet(const MaterialDataset& d, Real T, Real a, Real Q,
                                   Real R_channel, Real E_surface);

/// The correlation itself.  ALWAYS fails: see the header.  It exists so that
/// the refusal has a place and a reason, and so that no caller can drift into
/// believing a droplet current came out of this project.
[[noreturn]] void cone_jet_current();

}  // namespace es
