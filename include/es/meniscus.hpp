#pragma once
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/fluid.hpp"
#include "es/geometry.hpp"
#include "es/status.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// Static electrified meniscus
// ===========================================================================
//
// The free surface of a perfectly conducting liquid pinned at a contact line of
// radius r_c satisfies the Young-Laplace equation with a Maxwell traction:
//
//     gamma ( dphi/ds + sin(phi)/r ) = dp - rho g z + (eps0/2) E_n^2
//
// with s the arclength from the apex, phi the angle of the surface tangent
// below the horizontal, and the two terms on the left the meridional and
// azimuthal curvatures.  The electric term always pulls outward, so raising the
// voltage inflates the meniscus.
//
// TERMINOLOGY -- READ THIS BEFORE USING ANY RESULT
//
// What this solver can locate is the TURNING POINT of a branch of static
// solutions: the maximum of U(h) at fixed feed pressure.  It is called a
// "static fold" throughout, never an "onset".  Three different things were
// previously conflated under that word:
//
//   (a) loss of static stability -- a bifurcation.  The fold is a CANDIDATE for
//       it, and only under assumptions about which parameter is held fixed and
//       which perturbation mode is considered.  No stability analysis is
//       performed here, so even (a) is not established.
//   (b) onset of emission in the sense of a measurable current -- an instrument
//       threshold.  In the pure ionic regime the evaporation rate is a smooth
//       exponential of the field with no threshold at all.
//   (c) transition to the cone-jet regime -- involves the flow and cannot be
//       determined from a static model.
//
// An emission onset may only be reported once a physical stability criterion
// has been implemented and checked; see docs/02_model_specification.md 2.4.
//
// WHAT IS AN INDEPENDENT VARIABLE, AND WHY IT MATTERS NUMERICALLY
//
// dp (the feed pressure) and U (the applied voltage) are the physical inputs;
// the apex height h comes out.  But h(U) has a vertical tangent at the onset of
// emission -- that IS the onset: a saddle-node bifurcation beyond which no
// static meniscus exists and the liquid erupts into a jet.  Marching in U can
// therefore never reach it; the iteration simply stops converging somewhere
// short of the fold and the result looks like a numerical failure rather than
// the physics it is.
//
// So the solver inverts the roles: h is prescribed and U is solved for.  That
// parameterisation is regular through the fold, and the whole branch including
// its turning point is traced.  The maximum of U(h) is reported as the static
// fold voltage -- not as an onset.
//
// ASSUMPTIONS
//  * static: no flow.  Valid near onset and at the low flow rates of the pure
//    ionic regime.  At cone-jet flow rates the viscous pressure drop along the
//    cone is not negligible and this model overestimates the apex height.
//  * perfectly conducting liquid.  Justified by the charge relaxation time,
//    ~1e-10 s for ionic liquids -- see Fluid::charge_relaxation_time().
//  * pinned contact line at r_c.  Correct for a capillary with a sharp edge.
//    For externally wetted or porous emitters the contact line is not pinned;
//    treating the tip as pinned at its radius is a first-order stand-in only.
//  * the free surface is monotone in r.  Overhanging menisci -- a pendant drop
//    bulged past a hemisphere -- are legitimate Young-Laplace solutions but are
//    not on the path to a Taylor cone, and the solver reports them as
//    unreachable rather than tracking them.  In practice this bounds the
//    reachable apex height at roughly r_c in the field-free limit; with a field
//    the cone-like branch goes far higher without ever turning vertical.

struct MeniscusParams {
  Real r_contact{1.0e-5};  ///< pinning radius [m]
  Real z_contact{0.0};     ///< axial position of the contact line [m]
  Real delta_p{0.0};       ///< feed pressure at the apex, above ambient [Pa]
  Real gamma{0.0452};      ///< surface tension [N/m]
  Real rho{0.0};           ///< density, for the hydrostatic term; 0 disables it

  int n_nodes{81};             ///< nodes on the free surface
  Real apex_clustering{1.8};   ///< >1 clusters nodes toward the apex
  Real h_far{0};               ///< element size on the electrodes (0 = builder default)

  int max_outer{40};       ///< shape <-> field iterations
  Real relax{0.6};         ///< under-relaxation of the shape update
  Real tol{2.0e-4};        ///< convergence: max node motion / r_contact
  /// Relative tolerance on the voltage for solve_at_voltage().  A result whose
  /// voltage misses the request by more than this is reported as
  /// SolveStatus::VoltageMismatch, never as converged.
  Real voltage_tol{1.0e-3};
  bool verbose{false};
};

struct MeniscusShape {
  std::vector<Vec2> nodes;  ///< contact line -> apex (BEM traversal order)
  Real height{0};           ///< apex above the contact line [m]
  Real arclength{0};        ///< meridian length of the free surface [m]
  Real apex_radius{0};      ///< radius of curvature at the apex [m]
  Real half_angle{0};       ///< local cone half-angle at mid-arc [rad]

  /// Free-surface mesh in BEM orientation (outward normal into the vacuum).
  Mesh to_mesh(Real potential) const;
};

/// Initial guess: an ellipse arc of height h pinned at r_c.  Reduces to a
/// spherical cap when h == r_c.
MeniscusShape initial_shape(Real r_c, Real z_c, Real h, int n_nodes, Real clustering);

struct MeniscusSolution {
  MeniscusShape shape;
  Real voltage{0};        ///< emitter-to-extractor voltage sustaining this shape [V]
  /// Voltage that was requested by solve_at_voltage(); 0 if the call was
  /// solve_at_height(), where no voltage target exists.
  Real target_voltage{0};
  Real apex_field{0};     ///< |E_n| at the apex [V/m]
  Real peak_field{0};     ///< max |E_n| on the free surface [V/m]
  Real delta_p{0};
  SolveStatus status{SolveStatus::NotAttempted};
  int iterations{0};
  Real residual{0};       ///< final max node motion / r_contact

  /// The only admissible success test.  Anything else must be treated as a
  /// failure, however plausible the numbers look.
  bool ok() const { return is_usable(status); }
};

/// Couples the Young-Laplace shape solver to the BEM field solver.
class MeniscusSolver {
 public:
  /// `electrodes` must contain everything except the free surface: the emitter
  /// body (see make_capillary_open) and any extractor.  The emitter is held at
  /// U and the extractor at 0, so U is the emitter-to-extractor voltage.
  MeniscusSolver(Mesh electrodes, MeniscusParams params);

  /// Solve for the voltage that sustains an apex height h.  Warm-starts from
  /// the shape left by the previous call, which is what makes continuation
  /// cheap; pass a shape explicitly to override.
  MeniscusSolution solve_at_height(Real h, const MeniscusShape* start = nullptr);

  /// Solve for the meniscus on the rising side of the branch that a given
  /// voltage sustains.  Runs a coarse continuation to locate the static fold,
  /// then bisects on apex height below it.
  ///
  /// Returns SolveStatus::Converged ONLY if the voltage of the returned shape
  /// agrees with U to within params().voltage_tol (relative).  If U lies above
  /// the fold voltage or below the lowest voltage on the traced branch, the
  /// status is VoltageNotBracketed and the shape must not be used.
  MeniscusSolution solve_at_voltage(Real U, Real h_max, int scout_steps = 14);

  /// Put the solver's BEM into the state described by `sol`, so that any
  /// surface or mesh dump provably belongs to that state rather than to
  /// whatever the last internal iteration happened to leave behind.
  void realize(const MeniscusSolution& sol);

  /// Trace the branch h = h_min ... h_max.  The onset voltage is the maximum of
  /// `voltage` over the returned branch; everything beyond it is unstable.
  std::vector<MeniscusSolution> continuation(Real h_min, Real h_max, int n_steps);

  /// Turning point of the traced branch: the maximum of U(h).
  ///
  /// This is a STATIC FOLD, not an emission onset -- see the terminology note
  /// at the top of this header.  It is only reported when it is a genuine
  /// INTERIOR maximum: at least three converged points, strictly rising before
  /// and strictly falling after.  A branch with a single point, or one that is
  /// monotone throughout, has no turning point and says so.
  struct StaticFold {
    FoldStatus status{FoldStatus::NotAttempted};
    Real voltage{0};      ///< fold voltage [V], parabolically refined
    Real height{0};
    Real apex_field{0};
    Real apex_radius{0};
    Real half_angle{0};
    std::size_t index{0}; ///< index of the sampled maximum in the branch
    bool found() const { return status == FoldStatus::Found; }
  };
  static StaticFold find_static_fold(const std::vector<MeniscusSolution>& branch);

  /// The solver keeps the last assembled system, so callers can post-process
  /// (ion emission, beam launch) without rebuilding it.
  const BemSolver& bem() const { return bem_; }
  BemSolver& bem() { return bem_; }

  const MeniscusParams& params() const { return params_; }

  static void write_branch_csv(const std::string& path,
                               const std::vector<MeniscusSolution>& branch,
                               const std::string& header = {});

 private:
  Mesh electrodes_;
  MeniscusParams params_;
  BemSolver bem_;
  MeniscusShape last_;
  bool have_last_{false};
};

}  // namespace es
