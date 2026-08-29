#pragma once
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "es/status.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P4 -- the time-dependent free surface: what is implemented and what is not
// ===========================================================================
//
// THE PRECONDITION, CHECKED FIRST
//
// A dynamic free surface may be implemented only on a closed flow and
// charge-transport foundation.  P3 does NOT provide one:
//
//   * its flow model is fully developed pipe flow in a STRAIGHT tube.  It has
//     no free surface at all, and its exactness rests on du_z/dz = 0, which a
//     deforming surface destroys;
//   * its charge model is STEADY conduction with no surface charge, no surface
//     charge convection and -- the structural gap -- no tangential traction
//     q_s E_t.  That traction is exactly what drives the liquid in a finitely
//     conducting meniscus and it does not exist in the perfect-conductor limit;
//   * eps_r is MissingMaterialData, so even the relaxation time is not
//     computable from this project's own sourced data.
//
// So a dynamic meniscus solver is NOT implemented.  solve_dynamic_meniscus()
// below throws NotImplementedInThisPhase and names what is missing.  The full
// contract -- state variables, boundary conditions, time-step rule, required
// discretisation, stability and energy checks -- is written down in
// docs/15_free_surface_dynamics.md so that the gap is specified rather than
// merely admitted.
//
// WHAT IS IMPLEMENTED AND VALIDATED HERE
//
// Exactly one piece, and it is the one that can be validated against an exact
// solution: the KINEMATIC BOUNDARY CONDITION, integrated in time on a
// PRESCRIBED velocity field.  That is kinematics, not dynamics: the velocity
// comes from the caller, nothing is solved for it, and no force balance is
// evaluated anywhere in this file.
//
// THREE PROHIBITIONS, enforced by there being no code for them:
//
//   * NO MOBILITY.  There is no coefficient anywhere that relaxes a shape
//     towards equilibrium at a chosen rate.  Such a coefficient would be a
//     free parameter dressed as physics.
//   * NO ARTIFICIAL DAMPING.  Nothing is smoothed, filtered or relaxed.  The
//     tangential node redistribution below moves nodes ALONG the surface and
//     therefore cannot change the surface; that it does not is measured.
//   * NO STABILITY CLAIM.  The end of a static branch is not a dynamic
//     instability and nothing here turns it into one.

// ---------------------------------------------------------------------------
// The state
// ---------------------------------------------------------------------------

/// The free surface as a polyline in the meridian half-plane, apex first.
/// This is the STATE VARIABLE of the kinematic problem and nothing else is.
struct SurfacePolyline {
  std::vector<Vec2> nodes;
  Real time{0};

  /// Outward normal at node k, from the neighbouring segments.  At the apex
  /// (k = 0, on the axis) axisymmetry forces the normal to be axial.
  Vec2 normal_at(std::size_t k) const;
  Vec2 tangent_at(std::size_t k) const;
  Real arclength() const;
  /// Volume enclosed between the surface and the plane z = z_base [m^3],
  /// by the exact revolved integral of the polyline.
  Real revolved_volume(Real z_base) const;
  Real revolved_area() const;
  /// Largest deviation of the segment lengths from their mean, relative.  A
  /// measure of how badly the nodes have bunched, and the reason the tangential
  /// redistribution exists.
  Real node_spacing_nonuniformity() const;
};

/// A prescribed velocity field in the meridian half-plane [m/s].
using VelocityField = std::function<Vec2(Vec2, Real)>;

// ---------------------------------------------------------------------------
// The kinematic boundary condition
// ---------------------------------------------------------------------------
//
//     dx/dt . n = u . n            on the free surface.
//
// ONLY THE NORMAL COMPONENT IS PHYSICS.  The tangential motion of a node is a
// property of the parametrisation, not of the surface: a node that slides along
// the surface leaves the surface unchanged.  Two modes are therefore offered
// and they are NOT interchangeable:
//
//   Lagrangian     dx/dt = u.  Nodes are material points.  This is the mode in
//                  which an exact solution can be written down, because the
//                  node trajectories are the characteristics of the field.
//   NormalOnly     dx/dt = (u . n) n, plus an explicit tangential
//                  redistribution that keeps the nodes equally spaced in
//                  arclength.  The surface is the same; the parametrisation is
//                  not.  The redistribution is MESH MOTION and is labelled so.
//
// The difference between the two after the same time is a measure of the
// redistribution error and is reported, not hidden.
enum class KinematicMode {
  Lagrangian = 0,
  NormalOnly,
};
const char* to_string(KinematicMode m);

/// Contact-line treatment.  Pinned is the P3a/P3b condition; Free is used by
/// the validation cases, whose exact solutions move the whole surface.
enum class ContactLine { Pinned = 0, Free };
const char* to_string(ContactLine c);

// ---------------------------------------------------------------------------
// The time-step contract
// ---------------------------------------------------------------------------
//
// A step is admissible only if the surface it produces is still a surface.  The
// three conditions below are checked AFTER every step and a failing step is
// rejected, never accepted with a warning:
//
//   1. no node crosses the axis (r < 0);
//   2. no segment reverses or collapses (the polyline stays simple in the
//      sense that every segment keeps positive length);
//   3. the CFL-like bound: no node moves further than kMaxNodeMotion times the
//      shortest segment length in one step.
//
// The last is a discretisation bound, not a physical one, and it is stated as
// such.  Nothing here integrates a force, so no acoustic or capillary time-step
// limit applies -- those belong to the dynamic solver that does not exist.
namespace kinematics {
inline constexpr Real kMaxNodeMotion = 0.25;
inline constexpr Real kMinSegmentFraction = 1.0e-6;
}  // namespace kinematics

enum class StepStatus {
  Ok = 0,
  NodeCrossedAxis,
  SegmentCollapsed,
  StepTooLarge,
  NotAttempted,
};
const char* to_string(StepStatus s);
const char* explain(StepStatus s);

struct AdvectionResult {
  StepStatus status{StepStatus::NotAttempted};
  SurfacePolyline surface;
  int steps{0};
  Real dt{0};
  Real max_node_motion_fraction{0};  ///< worst node motion / shortest segment
  /// Volume at the start and at the end, and their relative change.  For a
  /// DIVERGENCE-FREE field the exact change is zero, so this is a pure measure
  /// of the time integration.
  Real volume_initial{0}, volume_final{0}, volume_change{0};
  Real spacing_nonuniformity{0};
  std::string message;
  void print(std::FILE* out) const;
};

/// Integrate the kinematic condition over `n_steps` steps of `dt` with
/// classical RK4 on the node trajectories.  `u` is PRESCRIBED; nothing is
/// solved for it.
AdvectionResult advect_surface(const SurfacePolyline& initial, const VelocityField& u, Real dt,
                               int n_steps, KinematicMode mode, ContactLine contact,
                               Real z_base);

// ---------------------------------------------------------------------------
// The analytically controllable validation cases
// ---------------------------------------------------------------------------
//
// Both are prescribed fields whose Lagrangian solution is a closed-form map, so
// the exact node positions at any time are known -- not merely the volume.
//
//   Dilation      u = alpha * x.  A sphere of radius R stays a sphere of radius
//                 R e^{alpha t}; the volume grows as e^{3 alpha t}.  The field
//                 has div u = 3 alpha, so this case tests that the volume
//                 follows a KNOWN NON-ZERO change rather than being conserved
//                 by accident.
//   Squeeze       u = (-alpha r / 2, alpha z).  Its divergence is
//                 (1/r) d(r u_r)/dr + du_z/dz = -alpha + alpha = 0 EXACTLY, so
//                 the enclosed volume is exactly conserved while the shape
//                 changes a lot.  A sphere becomes a spheroid with semi-axes
//                 R e^{-alpha t/2} and R e^{alpha t}.
//
// Between them they separate "the volume is right" from "the shape is right",
// which a single case cannot.

VelocityField dilation_field(Real alpha);
VelocityField squeeze_field(Real alpha);

/// Exact Lagrangian map of the two fields above.
Vec2 dilation_exact(Vec2 x0, Real alpha, Real t);
Vec2 squeeze_exact(Vec2 x0, Real alpha, Real t);

/// A half sphere of radius R centred at the origin, apex on the axis at +R,
/// discretised with `n` nodes uniform in polar angle.  Its revolved volume
/// above z = 0 is exactly (2/3) pi R^3.
SurfacePolyline hemisphere(Real R, std::size_t n);

// ---------------------------------------------------------------------------
// What is NOT implemented
// ---------------------------------------------------------------------------

/// Fails closed, always.  See the header comment and
/// docs/15_free_surface_dynamics.md.  It exists so that the gap has a name and
/// a place, and so that no caller can drift into believing it works.
[[noreturn]] void solve_dynamic_meniscus();

}  // namespace es
