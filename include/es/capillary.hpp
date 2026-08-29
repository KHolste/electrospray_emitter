#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/device_geometry.hpp"
#include "es/liquid.hpp"
#include "es/status.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P3a -- static capillary meniscus, NO electric field
// ===========================================================================
//
// WHAT IS SOLVED
//
// The shape of the free liquid surface that closes the bore, held by surface
// tension alone against a prescribed liquid pressure.  Axisymmetric
// Young-Laplace:
//
//     gamma * kappa = delta_p_exit,     kappa = dpsi/ds + sin(psi)/r
//
// with kappa the SUM of the two principal curvatures (twice the mean
// curvature), s the arclength measured from the apex on the axis, and psi the
// angle of the tangent.  Nothing else is on the right-hand side: no Maxwell
// stress, no hydrostatic term, no viscous term.  See "WHAT IS NOT IN HERE".
//
// SIGN CONVENTION -- fixed here, tested, and printed on every figure
//
//   delta_p_exit = p_liquid - p_vacuum  at the exit plane z = 0  [Pa]
//
//   delta_p_exit > 0  ->  the surface bulges OUT of the bore, towards the
//                         extractor, i.e. towards +z.  Apex height h > 0.
//                         The centre of curvature lies inside the liquid.
//   delta_p_exit = 0  ->  the surface is EXACTLY flat, z == 0 everywhere.
//   delta_p_exit < 0  ->  the surface is drawn INTO the bore, towards -z.
//                         Apex height h < 0.
//
// The outward normal points from the liquid into the vacuum, n = (sin psi,
// cos psi) with the parametrisation below, so kappa > 0 means "convex seen from
// the vacuum" and gamma*kappa is the pressure jump the surface can carry.
//
// PARAMETRISATION -- and why not z(r)
//
// The state is (r(s), z(s), psi(s)) with
//
//     dr/ds = cos psi,   dz/ds = -sin psi,   dpsi/ds = kappa - sin(psi)/r
//
// integrated from the apex outwards.  A height function z(r) has an infinite
// derivative as soon as the surface turns vertical, which happens already at
// the hemispherical limit of this very problem and always in a Taylor cone.
// Building on z(r) now would force a rewrite later, so it is not built.
//
// THE AXIS IS TREATED ANALYTICALLY, NOT GUESSED
//
// At r = 0 the term sin(psi)/r is 0/0.  Regularity of a smooth surface makes
// the two principal curvatures equal there, so
//
//     lim_{r->0} sin(psi)/r = dpsi/ds   and therefore   dpsi/ds|_apex = kappa/2.
//
// That limit is what the integrator uses at r = 0; there is no division by zero
// and no small-r fudge factor anywhere.
//
// BOUNDARY CONDITIONS
//
//   * axis r = 0: regularity, psi(0) = 0 (horizontal tangent at the apex).
//   * contact line: PINNED at the sharp exit edge, r = phi_2/2, z = 0, taken
//     directly from DeviceParameters.  No radius is hard-coded here.
//   * NO contact angle.  A pinned contact line and a prescribed Young angle are
//     two mutually exclusive descriptions of the same edge; imposing both
//     over-determines the problem.  Asking for both is refused, it is not
//     silently resolved in favour of one of them.
//
// WHAT IS NOT IN HERE -- and must not be read into any result
//
//   * no electric field, no Maxwell pressure, no coupling to the electrostatic
//     solver, and therefore no operating voltage;
//   * no emission, no space charge, no beam;
//   * no flow: no viscous pressure drop, no feed impedance, no Poiseuille term;
//   * no time dependence and no stability analysis.  A solution existing here
//     says nothing about whether it is stable;
//   * no Taylor cone and no cone-jet.  The shapes here are spherical caps, and
//     the range in which they exist at all is bounded (see kMaxAbsPi);
//   * no gravity.  Left out on purpose, and the Bond number that justifies it
//     is computed and reported rather than assumed -- see
//     LiquidProperties::bond_number();
//   * the bore is ASSUMED to be full up to the exit edge.  P3a does not compute
//     the capillary rise from the reservoir; the upstream state enters solely
//     through delta_p_exit.
//   * the P2c reservoir geometry is untouched and is not solved volumetrically.

// ---------------------------------------------------------------------------
// Status contract.  Deliberately not a bool: "did not converge", "the pressure
// admits no smooth shape at all" and "the material data are unusable" are three
// different answers and a caller must be able to tell them apart.
// ---------------------------------------------------------------------------

enum class CapillaryStatus {
  Solved = 0,
  /// No smooth surface pinned at this radius exists for this pressure: the
  /// meridian turns vertical before it reaches the contact radius.  |Pi| > 2.
  /// NOTHING is returned in this case -- no last iterate, no clipped shape.
  PressureOutsideCapillaryRange,
  /// |Pi| is at the hemispherical limit within round-off.  The hemisphere is a
  /// legitimate solution and its closed form is known, but the shooting
  /// condition r(L) = a becomes tangential there (a double root), so the
  /// numerical shape is NOT reported as if it had been solved.
  HemisphericalLimit,
  /// The pinning radius from the device geometry is not usable.
  InvalidGeometry,
  /// The material data set is unusable; LiquidProperties::why_unusable() says why.
  InvalidLiquid,
  /// A contact angle was prescribed in addition to the pinned contact line.
  ContactAngleAndPinningBothPrescribed,
  /// The refinement ran into its cap without reaching the requested accuracy.
  /// The shape is present but must not be used at the requested accuracy.
  AccuracyNotReached,
  /// The arclength search found no bracket within its bounds -- a numerical
  /// failure, distinct from the physical one above.
  ArclengthNotBracketed,
  NotAttempted,
};
const char* to_string(CapillaryStatus s);
const char* explain(CapillaryStatus s);
inline bool is_usable(CapillaryStatus s) { return s == CapillaryStatus::Solved; }

// ---------------------------------------------------------------------------

namespace capillary {

/// Dimensionless pressure  Pi = delta_p * a / gamma,  a = phi_2/2.
///
/// A smooth cap pinned at a exists exactly for |Pi| <= 2: the sphere radius is
/// R = 2*gamma/delta_p, and it must not be smaller than the pinning radius.
/// Pi = +-2 is the hemisphere, where the surface meets the edge vertically.
inline constexpr Real kMaxAbsPi = 2.0;

/// Distance from kMaxAbsPi within which the shooting condition is tangential.
/// Chosen from the conditioning of the root, not from a result: the residual
/// derivative is cos(psi_contact) = sqrt(1 - (Pi/2)^2), so at |Pi| = 2 - 1e-9
/// the arclength is still recoverable to about 1e-13 of a bore radius.
inline constexpr Real kHemisphereBand = 1.0e-9;

Real pi_from_pressure(Real delta_p, Real a, Real gamma);
Real pressure_from_pi(Real Pi, Real a, Real gamma);

}  // namespace capillary

// ---------------------------------------------------------------------------
// Closed-form reference: the spherical cap.  Independent of the solver -- it
// shares no code with it, on purpose, so that agreement means something.
// ---------------------------------------------------------------------------

struct SphericalCap {
  Real a{0};             ///< pinning radius [m]
  Real gamma{0};         ///< surface tension [N/m]
  Real delta_p{0};       ///< pressure difference [Pa]
  Real Pi{0};            ///< dimensionless pressure
  Real sphere_radius{0}; ///< signed, 2*gamma/delta_p; infinite for delta_p = 0
  Real curvature{0};     ///< kappa = delta_p/gamma, constant over the whole cap
  Real apex_height{0};   ///< signed, apex above the contact plane [m]
  Real arclength{0};     ///< meridian length apex -> contact line [m]
  Real revolved_area{0}; ///< 2*pi*int r ds [m^2]
  Real revolved_volume{0};        ///< signed volume between cap and z = 0 [m^3]
  Real contact_tangent_angle{0};  ///< psi at the contact line [rad], signed

  /// z of the cap at radius r in [0, a], with z(a) = 0.
  Real z_at_radius(Real r) const;
};

/// Throws std::runtime_error if |Pi| > 2 -- the closed form does not exist
/// there either, and must not be evaluated as if it did.
SphericalCap spherical_cap(Real a, Real delta_p, Real gamma);
bool spherical_cap_exists(Real a, Real delta_p, Real gamma);

// ---------------------------------------------------------------------------

struct CapillaryRequest {
  /// p_liquid - p_vacuum at the exit plane [Pa].
  Real delta_p_exit{0.0};

  /// Requested relative accuracy of the profile, in units of the bore radius.
  /// This is the ONLY discretisation input.  The number of arclength intervals
  /// is chosen automatically from it by refinement until the change between two
  /// successive resolutions falls below it; there is no mesh size, no h_tip and
  /// no h_far anywhere in this module.
  Real target_relative_accuracy{1.0e-10};

  /// Upper bound of the automatic refinement.  Reaching it without meeting the
  /// accuracy is reported as AccuracyNotReached, never accepted silently.
  int max_intervals{131072};

  /// STUDY CONTROL, NOT A USER INPUT.  0 = automatic refinement as described
  /// above.  A positive value fixes the number of arclength intervals and is
  /// used by the mesh-convergence study, which has to see a chosen resolution
  /// by construction.  Applications must not expose this through a config key.
  int forced_intervals{0};

  /// A Young contact angle, if a caller insists on prescribing one.  Any value
  /// is REFUSED together with the pinned contact line
  /// (ContactAngleAndPinningBothPrescribed); it exists so that the refusal is
  /// explicit rather than the option being absent and silently ignored.
  Real prescribed_contact_angle_deg{0.0};
  bool contact_angle_prescribed{false};
};

/// Point-wise Young-Laplace residual, evaluated from the NODE COORDINATES only.
/// It shares nothing with the integrator: the tangent angles come from the
/// chords and the meridional curvature from their turning, so a solver that
/// merely satisfies its own discrete equations cannot pass it.
struct ResidualProfile {
  std::vector<Real> s;          ///< arclength of the evaluated node [m]
  std::vector<Real> residual;   ///< (gamma*kappa - delta_p) * a / gamma  [-]
  Real max_abs{0};
  Real rms{0};
};

struct CapillaryMeniscus {
  CapillaryStatus status{CapillaryStatus::NotAttempted};
  std::string message;   ///< non-empty whenever status != Solved

  // --- the problem that was posed ------------------------------------------
  Real contact_radius{0};   ///< a, from DeviceParameters -- never hard-coded
  Real contact_z{0};        ///< z of the pinned edge, from the device geometry
  Real delta_p_exit{0};
  Real gamma{0};
  Real Pi{0};

  // --- the discretisation that was chosen ----------------------------------
  int n_intervals{0};
  Real estimated_relative_error{0};  ///< difference to the half-resolution run
  bool discretisation_was_forced{false};

  // --- the solution ---------------------------------------------------------
  std::vector<Vec2> nodes;  ///< apex -> contact line, uniform in arclength
  std::vector<Real> psi;    ///< tangent angle at each node [rad]
  Real apex_height{0};      ///< signed [m]
  Real arclength{0};        ///< [m]
  Real revolved_area{0};    ///< [m^2], integrated with the ODE
  Real revolved_volume{0};  ///< signed [m^3], integrated with the ODE
  Real contact_tangent_angle{0};  ///< psi at the contact line [rad]

  /// The same two measures taken from the polyline instead of from the ODE,
  /// using the tested helpers in device_geometry.hpp.  Second-order accurate,
  /// so they differ from the ODE values by the discretisation error -- which is
  /// the point: two independent evaluations of the same quantity.
  Real polyline_area{0};
  Real polyline_volume{0};

  const Vec2& apex() const { return nodes.front(); }
  const Vec2& contact() const { return nodes.back(); }
};

/// Solve the static meniscus for the exit edge of `geometry`.
///
/// The pinning radius and the axial position of the contact line are taken from
/// the device geometry's PinnedContactEdge feature.  No geometry is described a
/// second time in this module.
CapillaryMeniscus solve_capillary_meniscus(const DeviceGeometry& geometry,
                                           const LiquidProperties& liquid,
                                           const CapillaryRequest& request);

/// Same solver on a bare pinning radius, for the dimensionless tests and for
/// the scaling study.  The device path above is the one applications use.
CapillaryMeniscus solve_capillary_meniscus(Real contact_radius, Real contact_z,
                                           const LiquidProperties& liquid,
                                           const CapillaryRequest& request);

ResidualProfile young_laplace_residual(const CapillaryMeniscus& m);

/// Largest distance of a node from the analytic cap SURFACE, measured normal to
/// it and divided by the pinning radius.  The normal distance is used instead
/// of a vertical difference because dz/dr diverges as the hemispherical limit
/// is approached, where a vertical difference would report the conditioning of
/// the comparison rather than the error of the solution.  Throws if the cap
/// does not exist.
Real profile_error_against_cap(const CapillaryMeniscus& m);

/// Largest vertical difference |z_numeric - z_cap(r)| over the nodes, divided
/// by the pinning radius.  Reported alongside the normal distance because it is
/// the difference a reader sees in a profile plot.  Throws if the cap does not
/// exist.
Real profile_z_error_against_cap(const CapillaryMeniscus& m);

void print(const CapillaryMeniscus& m, std::FILE* out);

}  // namespace es
