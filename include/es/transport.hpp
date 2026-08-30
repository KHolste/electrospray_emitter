#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/axisym_fem.hpp"
#include "es/material_data.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P3 -- finite conductivity and liquid flow, in the smallest honest step
// ===========================================================================
//
// P3b treats the liquid as an IDEAL CONDUCTOR: the free surface is an
// equipotential, all charge sits on it, and no current flows anywhere.  That is
// a limit, and like every limit it is worth what its small parameter is worth.
// This file makes the small parameter computable and states the contract of the
// finite-conductivity problem it is a limit of.  It does NOT solve that problem
// on a moving meniscus.
//
// TWO SEPARATE THINGS LIVE HERE, and they are separate on purpose:
//
//   1. FULLY DEVELOPED PIPE FLOW.  The steady axisymmetric Stokes problem in a
//      straight circular pipe, solved on the existing axisymmetric FEM and
//      validated against Hagen-Poiseuille.  It is the flow counterpart of the
//      hydraulic resistance of P1: P1 asserts the closed form, this solves the
//      field and shows the closed form is what comes out.
//
//   2. CHARGE TRANSPORT.  The contract for conduction in the liquid --
//      j = sigma E, charge conservation, the relaxation time tau = eps/sigma --
//      with the steady conduction problem solved on the same FEM and validated
//      against the closed form, and with the perfect-conductor limit stated as
//      a measurable criterion rather than an assumption.
//
// WHAT IS NOT HERE, and must not be read into any number:
//
//   * no general Stokes solver.  No pressure-velocity coupling, no entrance
//     flow, no curved geometry, no free surface.  The reduction below is exact
//     ONLY for a straight pipe with a fully developed profile.
//   * no coupled finite-conductivity meniscus.  The boundary conditions on a
//     deforming free surface with surface charge convection are written down in
//     docs/14 and are NOT implemented.  Nothing here returns a meniscus shape.
//   * no emission.  Consequently NO steady normal current may cross the free
//     surface: with no emission there is nowhere for the charge to go, and a
//     solver that produced one would be producing charge from nothing.  That is
//     a testable statement and it is tested.
//   * no time stepping of the free surface.  tau is computed and compared to
//     other time scales; nothing is integrated in time except the closed-form
//     relaxation of a uniform charge density, which is an ODE with a solution.

// ---------------------------------------------------------------------------
// 1.  Fully developed axisymmetric pipe flow
// ---------------------------------------------------------------------------
//
// WHY THE REDUCTION IS EXACT, AND WHERE IT STOPS
//
// For steady incompressible flow in a straight circular pipe, take
// u = u_z(r, z) e_z.  Continuity div u = 0 gives du_z/dz = 0 immediately, so
// u_z = u_z(r) and the convective term (u . grad) u vanishes IDENTICALLY --
// not by a small-Reynolds argument, but because u is orthogonal to its own
// gradient.  The z momentum equation is then
//
//     0 = -dp/dz + mu * (1/r) d/dr ( r du_z/dr ) ,
//
// and the r momentum equation reduces to dp/dr = 0.  The right-hand side is
// exactly the axisymmetric Laplacian of u_z, so
//
//     -div( mu grad u_z ) = -dp/dz
//
// is a scalar problem the existing axisymmetric FEM solves without any change
// of formulation: AxisymProblem with coefficient_scale = 1, the coefficient
// carrying mu, and a constant source -dp/dz.
//
// THIS IS NOT A STOKES SOLVER.  It cannot do an entrance region, a contraction,
// a bend or a free surface, because in all of those du_z/dz is not zero and the
// pressure is not a function of z alone.  What it IS is an independent check of
// the P1 hydraulic resistance: P1 writes down 8 mu L Q / (pi R^4), and this
// solves the field and integrates the profile to see whether that is what comes
// out.  The two share no code.

struct PipeFlowSolution {
  Index n_nodes{0}, n_cells{0};
  Real radius{0}, length{0}, mu{0}, dpdz{0};
  std::vector<Real> u;         ///< nodal axial velocity [m/s]
  Real centreline_velocity{0}; ///< [m/s]
  Real flow_rate{0};           ///< integral of u over the cross section [m^3/s]
  Real mean_velocity{0};       ///< [m/s]
  Real wall_shear_stress{0};   ///< from the recovered gradient at the wall [Pa]
  /// The closed forms the solution is judged against.
  Real flow_rate_closed_form{0};
  Real centreline_closed_form{0};
  /// Largest relative deviation of the nodal profile from
  /// u(r) = -(dp/dz)(R^2 - r^2)/(4 mu), over the nodes with u > 0.
  Real max_profile_error{0};
  Real flow_rate_error{0};
  Real fem_residual{0};
};

/// Solve the reduced problem on a straight cylinder of radius `R` and length
/// `L`, discretised `nr` x `nz`, with no slip on the wall and a prescribed
/// pressure gradient.  Throws on an invalid request.
PipeFlowSolution solve_pipe_flow(Real R, Real L, Real mu, Real dpdz, Index nr, Index nz);

// ---------------------------------------------------------------------------
// 2.  Charge transport: the contract
// ---------------------------------------------------------------------------
//
// THE EQUATIONS, stated so that what is implemented and what is not can be told
// apart.
//
//   conduction current   j = sigma E = -sigma grad phi        [A/m^2]
//   charge conservation  d rho / dt + div j = 0
//   Gauss                div( eps grad phi ) = -rho
//
// Combining the three for a homogeneous liquid (sigma, eps constant) gives
//
//     d rho / dt + (sigma / eps) rho = 0   =>   rho(t) = rho(0) exp(-t / tau),
//     tau = eps / sigma = eps0 eps_r / sigma .
//
// So free charge cannot live in the BULK of a homogeneous conductor: it decays
// with the relaxation time and ends up on the boundary.  That is the whole
// content of the perfect-conductor limit, and it is a statement about a RATIO
// of times, not about sigma being large.
//
//   STEADY STATE.  With d rho/dt = 0 and sigma constant, charge conservation
//   reduces to div( sigma grad phi ) = 0 -- the SAME operator as electrostatics
//   with eps replaced by sigma.  The nodal reaction of that operator is then a
//   CURRENT in ampere instead of a charge in coulomb.  That is why the existing
//   FEM solves it unchanged.
//
//   THE PERFECT-CONDUCTOR LIMIT is tau << t_process.  Which t_process is the
//   right one is a modelling decision, not a fact: for a static equilibrium it
//   is the time the shape takes to settle; for emission it is the transit time
//   of the charge across the emitting region.  This file therefore takes
//   t_process as an EXPLICIT INPUT and reports the ratio; it does not pick one.
//
//   WITHOUT EMISSION THERE IS NO STEADY NORMAL CURRENT through the free
//   surface.  Charge that crossed it would have nowhere to go, so conservation
//   forbids it: the correct boundary condition on a non-emitting free surface
//   is j . n = 0, i.e. the natural (zero-flux) condition of the same operator.
//   A solver that imposed a potential there instead would be driving a current
//   out of the liquid and calling it physics.  This is checked, not assumed.

/// Charge relaxation time tau = eps0 eps_r / sigma [s].
/// Returns NaN if either input is non-positive -- there is no default.
Real charge_relaxation_time(Real eps_r, Real sigma);

/// Closed-form decay of a uniform free charge density [C/m^3].
Real relaxed_charge_density(Real rho0, Real t, Real tau);

enum class ConductorLimit {
  NotAttempted = 0,
  /// tau is more than kPerfectConductorMargin times shorter than the process
  /// time: the equipotential treatment is justified for THAT process time.
  PerfectConductorJustified,
  /// tau is comparable to or longer than the process time.  The equipotential
  /// treatment is then an assumption, not a limit.
  FiniteConductivityRequired,
  /// eps_r or sigma is missing.  tau is not computable and NOTHING is claimed.
  MissingMaterialData,
};
const char* to_string(ConductorLimit v);
const char* explain(ConductorLimit v);

namespace transport {
/// How much shorter tau must be than the process time before the
/// perfect-conductor treatment is called justified.  Fixed before the
/// measurement; it is a declared convention, not a derived number.
inline constexpr Real kPerfectConductorMargin = 100.0;
}  // namespace transport

struct RelaxationVerdict {
  ConductorLimit limit{ConductorLimit::NotAttempted};
  Real eps_r{0}, sigma{0};
  Real tau{0};             ///< [s], NaN when not computable
  Real process_time{0};    ///< [s], the input it is compared against
  Real ratio{0};           ///< process_time / tau
  MaterialDataStatus eps_status{MaterialDataStatus::MissingMaterialData};
  MaterialDataStatus sigma_status{MaterialDataStatus::MissingMaterialData};
  std::string message;
  void print(std::FILE* out) const;
};

/// Judge the perfect-conductor limit from a MATERIAL DATASET, so that a missing
/// permittivity fails closed instead of silently using a plausible number.
///
/// This asks for ONE selected permittivity and therefore still reports
/// MissingMaterialData for EMI-BF4: no permittivity source states purity and
/// water content, so the selection rule of P2 selects none.  That is correct
/// and stays.  It is NOT the end of the question -- see the band and the
/// self-consistent relaxation below, which answer it without inventing a
/// single value.
RelaxationVerdict judge_conductor_limit(const MaterialDataset& d, Real T, Real process_time);

/// The same from explicit numbers, for tests and for a deliberate what-if.
RelaxationVerdict judge_conductor_limit_explicit(Real eps_r, Real sigma, Real process_time);

// ---------------------------------------------------------------------------
// 2b.  WHICH permittivity belongs in tau_q -- and why "not DC" is not a reason
//      to throw a measurement away
// ---------------------------------------------------------------------------
//
// An earlier version of this file demanded a single "DC permittivity" and
// rejected the 1-18 GHz measurement in the data set on the ground that it was
// "not DC".  Both halves of that were wrong, and in opposite directions.
//
// FIVE DIFFERENT THINGS GET CALLED "THE PERMITTIVITY" OF A CONDUCTIVE IONIC
// LIQUID.  They are not variants of one number:
//
//   (1) STATIC / LOW-FREQUENCY APPARENT PERMITTIVITY.  What a capacitance
//       bridge reports at kHz.  For a liquid with K ~ 1 S/m this is dominated
//       by the DC conduction term and by electrode polarisation, and reaches
//       values of 10^4 and beyond.  It is a property of the CELL, not of the
//       liquid, and no model may use it.
//
//   (2) INTRINSIC STATIC (DIELECTRIC) PERMITTIVITY eps_s.  The zero-frequency
//       limit of the DIELECTRIC dispersion after the conduction contribution
//       has been separated out.  For ionic liquids it is not measured directly
//       at zero frequency at all -- it is obtained by fitting a relaxation
//       model to microwave spectra and extrapolating.  The three "static"
//       values in this data set (12.8, 12.9, 13.6) are of this kind.
//
//   (3) FREQUENCY-RESOLVED COMPLEX PERMITTIVITY eps*(f) = eps'(f) - i eps''(f).
//       The Bennett et al. (2019) set in this data set is eps'(f) at 18
//       frequencies from 1 to 18 GHz, at 298.15 K, with stated uncertainties.
//
//   (4) ELECTRODE POLARISATION.  A measurement artefact of the cell: ions pile
//       up at the electrodes and screen the applied field, which inflates the
//       apparent capacitance.  It is what makes (1) useless and is the reason
//       microwave spectroscopy is used for these liquids in the first place.
//
//   (5) DC CONDUCTIVITY K.  A separate property with its own selected source.
//       Writing it as an imaginary permittivity K/(eps0 omega) is a bookkeeping
//       choice, not a second dielectric constant.
//
// WHICH ONE ENTERS tau_q.  This is not a matter of taste; it follows from the
// derivation.  Charge conservation plus Gauss give
//
//     d rho_f / dt + (K / (eps0 eps_r)) rho_f = 0 ,
//
// where eps_r is the BOUND-charge response -- the polarisation that follows the
// field while the free charge is decaying.  The free charge decays on the time
// scale tau itself, so its spectral content sits at angular frequency
// omega ~ 1/tau, i.e. at f* = 1/(2 pi tau).  The permittivity that belongs in
// the formula is therefore eps_r evaluated AT THAT FREQUENCY -- and since eps_r
// is dispersive, the equation for tau is IMPLICIT:
//
//     tau = eps0 * eps_r( f* ) / K ,      f* = 1 / (2 pi tau) .
//
// For this liquid K ~ 1.6 S/m and eps_r ~ 10 put f* in the low GHz range --
// exactly where Bennett et al. measured.  The 1-18 GHz data are therefore not
// "the wrong frequency"; they are the RIGHT frequency, and the reason the
// earlier version dismissed them was a misreading of the formula.  Equally,
// none of them may be quoted as a DC value: the loop below evaluates the
// measured curve at the frequency the physics selects, and reports whether
// that frequency lies inside the measured range or had to be extrapolated to.
//
// WHAT IS STILL MISSING, and stays missing.  Not one permittivity source states
// purity and water content, so the P2 selection rule selects none and
// material_value() keeps reporting MissingMaterialData.  No single released
// number rests on an unsourced eps_r.  What CAN be established is a justified
// BAND together with a sensitivity computation over it -- which is the second
// of the three outcomes P3 is allowed to reach.

enum class PermittivityConcept {
  StaticApparentLowFrequency = 0,
  IntrinsicStatic,
  FrequencyResolved,
  ElectrodePolarisation,
  DcConductivity,
};
const char* to_string(PermittivityConcept c);
const char* explain(PermittivityConcept c);

namespace transport {
/// Below this frequency a permittivity of a liquid with K of order 1 S/m is
/// dominated by electrode polarisation and DC conduction.  A point below it is
/// concept (1) and is refused as a dielectric datum.  Declared before looking
/// at the data; the data set happens to contain no such point at all, which is
/// checked rather than assumed.
inline constexpr Real kElectrodePolarisationFloor = 1.0e6;   ///< [Hz]
}  // namespace transport

/// The justified range of the dielectric permittivity at T, over every source
/// that is admissible as a DIELECTRIC datum.  No selection, no average, no
/// single value -- a band, and what it is made of.
struct PermittivityBand {
  bool ok{false};
  Real T{0};
  Real lo{0}, hi{0};                  ///< [-]
  std::size_t n_sources{0};           ///< admissible sources
  std::size_t n_static_points{0};     ///< concept (2), reported as static
  std::size_t n_frequency_points{0};  ///< concept (3)
  Real f_min{0}, f_max{0};            ///< [Hz] over the frequency-resolved points
  /// True if ANY point sits below kElectrodePolarisationFloor.  Such a point
  /// would be concept (1) and is excluded from the band; the flag exists so
  /// that the exclusion is reported instead of being invisible.
  bool any_below_polarisation_floor{false};
  std::string basis;
  std::string blocker;
  void print(std::FILE* out) const;
};
PermittivityBand permittivity_band(const MaterialDataset& d, Real T);

/// eps_r'(f) at temperature T from the MEASURED dispersion points, log-linear
/// between neighbouring measured frequencies.  Outside the measured range it
/// returns the nearest measured value and sets `*extrapolated` -- it never
/// fits a dispersion model, because inventing one would be inventing data.
/// Returns NaN if the data set has no frequency-resolved point at T.
Real dielectric_permittivity_at(const MaterialDataset& d, Real T, Real frequency_Hz,
                                bool* extrapolated);

/// tau from the IMPLICIT equation tau = eps0 eps_r(1/(2 pi tau)) / K, solved on
/// the measured dispersion curve by fixed-point iteration.  This is the value
/// the physics asks for; it is reported together with the frequency it was
/// taken at and with whether that frequency is inside the measured range.
struct SelfConsistentRelaxation {
  bool ok{false};
  Real T{0};
  Real sigma{0};                  ///< [S/m], the SELECTED conductivity
  Real eps_r{0};                  ///< [-] at f_star
  Real tau{0};                    ///< [s]
  Real f_star{0};                 ///< [Hz] = 1/(2 pi tau)
  int iterations{0};
  /// |eps_r(f*) - eps_r| / eps_r at the returned solution.  tau = eps0 eps_r/K
  /// and f* = 1/(2 pi tau) hold EXACTLY by construction; this is the one
  /// identity that a fixed point can only satisfy to within its convergence,
  /// and it is reported rather than being asserted away.
  Real residual{0};
  bool f_star_inside_measured{false};
  Real f_measured_min{0}, f_measured_max{0};
  /// tau computed with the intrinsic STATIC permittivity instead, for
  /// comparison.  The two differ by the dispersion between 0 and f_star and
  /// the difference is reported rather than hidden.
  Real tau_static{0}, eps_static{0};
  std::string blocker;
  void print(std::FILE* out) const;
};
SelfConsistentRelaxation self_consistent_relaxation(const MaterialDataset& d, Real T);

/// The perfect-conductor verdict evaluated over the WHOLE justified band of
/// eps_r AND the whole documented band of K -- i.e. at the corner that is worst
/// for the approximation.  `limit` is PerfectConductorJustified only if the
/// worst corner passes, so the verdict never rests on a convenient value.
struct BandedRelaxationVerdict {
  ConductorLimit limit{ConductorLimit::NotAttempted};
  Real T{0};
  Real process_time{0};
  Real eps_lo{0}, eps_hi{0};
  Real sigma_lo{0}, sigma_hi{0};
  Real tau_min{0}, tau_max{0};      ///< [s] over the four corners
  Real ratio_min{0}, ratio_max{0};  ///< process_time / tau, worst and best
  /// The self-consistent value, which lies inside the band and is the one the
  /// physics actually selects.
  Real tau_self_consistent{0};
  std::string message;
  void print(std::FILE* out) const;
};
BandedRelaxationVerdict judge_conductor_limit_over_band(const MaterialDataset& d, Real T,
                                                        Real process_time);

// ---------------------------------------------------------------------------
// 3.  Steady conduction in a cylinder
// ---------------------------------------------------------------------------

struct ConductionSolution {
  Index n_nodes{0};
  Real radius{0}, length{0}, sigma{0}, voltage{0};
  std::vector<Real> phi;       ///< nodal potential [V]
  Real current{0};             ///< total current through the driven end [A]
  Real current_closed_form{0}; ///< sigma * pi R^2 * V / L
  Real current_error{0};
  Real resistance{0};          ///< V / I [Ohm]
  Real resistance_closed_form{0};
  /// Largest |j . n| on the LATERAL surface, divided by the mean axial current
  /// density.  The lateral surface carries the natural zero-flux condition --
  /// this is the "no current leaves where it may not" check, and it is what a
  /// non-emitting free surface must also satisfy.
  Real lateral_leakage{0};
  Real max_potential_error{0};  ///< against the linear closed form
  Real fem_residual{0};
};

/// Steady conduction in a cylinder of radius R and length L with the two flat
/// ends held at 0 and V and the lateral surface insulating.  The exact solution
/// is phi = V z / L, which makes every quantity above closed form.
ConductionSolution solve_cylinder_conduction(Real R, Real L, Real sigma, Real V, Index nr,
                                             Index nz);

}  // namespace es
