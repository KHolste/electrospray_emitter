#pragma once
#include <cstdio>
#include <string>

#include "es/liquid.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P1 -- the pressure budget at the exit plane
// ===========================================================================
//
// WHAT THIS REPLACES, AND WHAT IT DOES NOT
//
// P3a and P3b take delta_p_exit as a free input.  That is honest but empty: the
// number that decides the meniscus shape has no relation to anything one could
// set on a bench.  This file makes it a BUDGET,
//
//     delta_p_exit = p_reservoir - p_vacuum
//                    - delta_p_hydrostatic
//                    - delta_p_viscous ,
//
// in which every term is either an explicit input or a closed-form expression
// of explicit inputs.  The direct input stays available as a controlled mode --
// PressureMode::Direct -- because every existing result was computed with it
// and must remain reproducible.
//
// WHAT THIS IS NOT.  It is NOT a flow solver.  It is a one-dimensional pressure
// bookkeeping with ONE closed-form hydraulic resistance, and it is worth
// exactly what its assumptions are worth:
//
//   * The feed channel is a straight circular tube, COMPLETELY FILLED with the
//     liquid, of constant radius R and length L.
//   * The flow in it is steady, laminar, fully developed and incompressible.
//   * The reservoir is a boundary condition, not a volume: its geometry does
//     not enter, and no reservoir emptying, no free surface inside it and no
//     gas cushion is modelled.
//   * The volumetric flow rate Q is an INPUT.  Nothing here computes it -- that
//     would require the emission that this phase does not have.
//   * The reservoir pressure is an INPUT.  Nothing here computes it either.
//
// NO CAPILLARY RISE IS CLAIMED.  A meniscus with a moving contact line inside
// the feed channel is not modelled and must not be read into any number here.
// In a completely filled channel there is no second free surface at all, so
// there is no Young angle to insert: the ONLY free surface in this model is the
// pinned meniscus at the exit, and its capillary pressure is what the P3a/P3b
// equilibrium computes -- it is not added here a second time.  Setting a
// contact angle is therefore refused rather than ignored.

// ---------------------------------------------------------------------------
// Sign convention, reference height and flow direction -- fixed here, once
// ---------------------------------------------------------------------------
//
// GEOMETRY.  z runs along the emitter axis and increases towards the extractor.
// The exit plane -- the plane the contact line is pinned in -- is the reference
// height z_exit.  The reservoir reference level z_reservoir is the height at
// which p_reservoir is stated; for the P1 device it lies BELOW the exit plane,
// so z_exit - z_reservoir > 0 is the length of liquid column that has to be
// lifted.
//
// GRAVITY.  gravity_axial is the component of the gravitational acceleration
// along +z, in m/s^2.  It is ZERO by default, because the device this code
// models is a spacecraft thruster and there is no gravity to account for; on a
// bench with the emitter pointing up it is -9.80665.  The hydrostatic term is
//
//     delta_p_hydrostatic = -rho * gravity_axial * (z_exit - z_reservoir) ,
//
// which for an upward-pointing emitter on Earth is +rho g H, a LOSS, as it must
// be: pressure at the top of a raised column is lower than at its foot.
//
// FLOW DIRECTION.  Q > 0 means liquid flows in +z, i.e. from the reservoir
// towards the exit -- the direction it flows when the emitter is fed.  The
// viscous term is then a LOSS and enters with a minus sign.  Q < 0 (suck-back)
// is allowed and simply reverses that sign; it is not refused, because nothing
// in the Poiseuille relation cares, but it is reported.
//
// EVERY TERM IS A SUBTRACTION FROM THE DRIVING PRESSURE.  A budget that added
// one of them would be a sign error and the tests pin each one separately.

enum class PressureMode {
  /// delta_p_exit is given directly.  The mode every earlier result used; kept
  /// so that they stay reproducible.
  Direct = 0,
  /// delta_p_exit is computed from the budget below.
  Budget,
};
const char* to_string(PressureMode m);

enum class FeedStatus {
  NotAttempted = 0,
  Ok,
  /// A required input of the chosen mode is missing or non-physical.
  MissingFeedInput,
  /// The channel radius or length is not positive.
  ChannelGeometryInvalid,
  /// The liquid data set does not carry a viscosity or a density.
  MissingLiquidProperty,
  /// The Reynolds number exceeds the laminar bound, so the Poiseuille relation
  /// does not apply.  Reported as a status, never silently used anyway.
  NotLaminar,
  /// The hydrodynamic entrance length is not short against the channel, so the
  /// flow is not fully developed and the closed form is not the whole loss.
  EntranceLengthNotShort,
};
const char* to_string(FeedStatus s);
const char* explain(FeedStatus s);
inline bool is_usable(FeedStatus s) { return s == FeedStatus::Ok; }

// ---------------------------------------------------------------------------

namespace feed {
/// Above this Reynolds number a straight circular pipe is no longer reliably
/// laminar.  The classical value; it is a bound on the validity of the closed
/// form, not a model parameter.
inline constexpr Real kReynoldsLaminar = 2300.0;
/// The hydrodynamic entrance length of laminar pipe flow is about
/// 0.06 Re d.  The fully developed profile is assumed over the WHOLE channel,
/// so that length must be small against L; this is the fraction allowed.
inline constexpr Real kEntranceFraction = 0.05;
}  // namespace feed

struct FeedChannel {
  Real radius{0.0};   ///< R [m]
  Real length{0.0};   ///< L [m]

  /// Hagen-Poiseuille hydraulic resistance of a straight circular tube,
  ///
  ///     R_h = 8 mu L / (pi R^4)      [Pa s / m^3] ,
  ///
  /// so that delta_p = R_h Q.  Derived from the parabolic profile
  /// u(r) = (delta_p / (4 mu L)) (R^2 - r^2), whose integral over the cross
  /// section is Q = pi R^4 delta_p / (8 mu L).  tests/test_feed.cpp checks the
  /// resistance against that integral computed independently.
  Real hydraulic_resistance(Real mu) const;

  /// Mean velocity Q / (pi R^2) [m/s].
  Real mean_velocity(Real Q) const;
  /// Centreline velocity of the developed profile, 2 * mean [m/s].
  Real centreline_velocity(Real Q) const;
  /// The developed axial velocity at radius r [m/s].  Only a diagnostic and a
  /// test reference -- the budget uses the resistance, not the profile.
  Real velocity_at(Real r, Real Q) const;
  /// Wall shear stress, tau_w = 4 mu u_mean / R = delta_p R / (2 L) [Pa].
  Real wall_shear_stress(Real mu, Real Q) const;
  /// Reynolds number 2 rho Q / (pi R mu) = rho u_mean d / mu [-].
  Real reynolds(Real rho, Real mu, Real Q) const;
  /// Hydrodynamic entrance length, 0.06 Re d [m].
  Real entrance_length(Real rho, Real mu, Real Q) const;
};

struct FeedRequest {
  PressureMode mode{PressureMode::Direct};

  /// Used in PressureMode::Direct, and reported in both.
  Real delta_p_exit_direct{0.0};

  // --- inputs of the budget -------------------------------------------------
  Real p_reservoir{0.0};   ///< absolute pressure of the liquid at z_reservoir [Pa]
  Real p_vacuum{0.0};      ///< ambient pressure on the vacuum side [Pa]
  Real z_exit{0.0};        ///< height of the exit plane [m]
  Real z_reservoir{0.0};   ///< height at which p_reservoir is stated [m]
  Real gravity_axial{0.0}; ///< component of g along +z [m/s^2]; 0 in orbit
  Real Q{0.0};             ///< volumetric flow rate, positive towards the exit [m^3/s]
  FeedChannel channel;

  /// A contact angle has no place here; see the header.  Kept as a field only
  /// so that a caller that sets it gets a clear refusal instead of silence.
  bool contact_angle_requested{false};
};

struct PressureBudget {
  FeedStatus status{FeedStatus::NotAttempted};
  std::string message;
  PressureMode mode{PressureMode::Direct};

  Real delta_p_exit{0.0};        ///< the result [Pa]

  // --- the terms, each with its own sign already applied --------------------
  Real driving{0.0};             ///< p_reservoir - p_vacuum [Pa]
  Real hydrostatic{0.0};         ///< delta_p_hydrostatic [Pa], subtracted
  Real viscous{0.0};             ///< delta_p_viscous [Pa], subtracted

  // --- what the channel is doing --------------------------------------------
  Real hydraulic_resistance{0.0};  ///< [Pa s / m^3]
  Real mean_velocity{0.0};         ///< [m/s]
  Real wall_shear_stress{0.0};     ///< [Pa]
  Real reynolds{0.0};
  Real entrance_length{0.0};       ///< [m]
  Real entrance_fraction{0.0};     ///< entrance_length / L

  // --- where the result sits on the capillary scale -------------------------
  Real gamma_over_a{0.0};          ///< [Pa]
  Real Pi{0.0};                    ///< delta_p_exit / (gamma/a), the P3a number
  /// |delta_p_exit| <= 2 gamma / a is the range in which a pinned static shape
  /// exists at all (P3a).  Reported, so that a budget outside it is visible
  /// before the meniscus solver is asked.
  bool within_capillary_range{false};

  void print(std::FILE* out) const;
};

/// Evaluate the budget.  Fails closed: a missing input gives a status, never a
/// silent default.
PressureBudget solve_pressure_budget(const FeedRequest& q, const LiquidProperties& liquid,
                                     Real contact_radius);

}  // namespace es
