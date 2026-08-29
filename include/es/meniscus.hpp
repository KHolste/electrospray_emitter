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

/// Which intersection of the branch with the target voltage is meant.
///
/// The names refer to the APEX HEIGHT and to nothing else.  In particular they
/// are not "stable" and "unstable": no stability analysis is implemented, so
/// no such claim may be encoded in an identifier.  Which of the two, if either,
/// is dynamically stable is open until phase P3.
enum class BranchSide {
  /// Refuse to choose.  A target voltage with more than one solution returns
  /// SolveStatus::AmbiguousBranch.
  Unspecified = 0,
  LowerHeight,   ///< the intersection at the smaller apex height
  UpperHeight,   ///< the intersection at the larger apex height
};

const char* to_string(BranchSide s);

/// Why the traced branch ended where it did.  Needed to tell "the branch really
/// stops here" from "we stopped looking here".
enum class BranchTermination {
  NotTraced = 0,
  /// The continuation covered the whole requested apex-height range.
  ReachedRequestedRange,
  /// The continuation aborted early because the shape solver stopped
  /// converging.  Nothing is known about the branch past that point.
  SolverStopped,
};

const char* to_string(BranchTermination s);

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
  /// How many times the branch crosses the requested voltage WITHIN THE RANGE
  /// THAT WAS ACTUALLY TRACED.  This says nothing about the rest of the branch;
  /// read coverage_complete before drawing any conclusion about uniqueness.
  int crossings_in_range{0};
  /// True only when the traced range is sufficient to decide the question for
  /// this target: the branch has passed a turning point and has moved away
  /// below the target at the far end, so no further crossing follows within the
  /// monotone segment that was examined.  False means "not investigated far
  /// enough", not "no further solution".
  bool coverage_complete{false};
  /// Complement of coverage_complete, stated positively for callers that want
  /// to warn.  A further crossing MAY exist outside the traced range; whether
  /// it does is unknown.
  bool additional_crossing_possible{true};
  /// Why the traced branch ended.
  BranchTermination termination_reason{BranchTermination::NotTraced};
  /// Which side was delivered.  When branch_crossings == 1 the choice was
  /// vacuous -- there was only one solution, and it is returned whichever
  /// side was asked for.  Read branch_crossings to tell the two cases apart.
  BranchSide side{BranchSide::Unspecified};
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
  /// U(h) rises to the turning point and falls again, so a target voltage below
  /// the fold generally has TWO solutions of very different shape -- measured on
  /// the reference geometry at 1154 V: apex heights 0.30 and 0.63 r_c, apex
  /// radii differing by a factor 3.3 and apex fields by a factor 1.8.  Picking
  /// one of them silently is not admissible.
  ///
  /// Crossings are counted ONLY within the apex-height range that was actually
  /// traced.  A single crossing there does not prove global uniqueness: the
  /// branch may cross again beyond h_max, or the continuation may have stopped
  /// before reaching a turning point.  `coverage_complete` records whether the
  /// traced range was sufficient to decide.
  ///
  /// Contract:
  ///
  ///  * `Unspecified`
  ///      - two or more crossings in range           -> AmbiguousBranch
  ///      - exactly one, coverage incomplete         -> BranchCoverageIncomplete
  ///      - exactly one, coverage complete           -> Converged
  ///      - none, coverage complete                  -> VoltageNotBracketed
  ///      - none, coverage incomplete                -> BranchCoverageIncomplete
  ///    Converged therefore only ever follows from demonstrated uniqueness.
  ///
  ///  * `LowerHeight`
  ///      - at least one crossing in range           -> that solution, Converged
  ///        (coverage_complete / additional_crossing_possible still reported, so
  ///         the caller can see whether more solutions may exist)
  ///      - none, coverage complete                  -> VoltageNotBracketed
  ///      - none, coverage incomplete                -> BranchCoverageIncomplete
  ///
  ///  * `UpperHeight`
  ///      - two or more crossings in range           -> the highest, Converged
  ///      - exactly one and coverage complete        -> that unique solution
  ///      - otherwise                                -> BranchCoverageIncomplete.
  ///        Never the lower crossing as a substitute.
  ///
  ///  * `VoltageNotBracketed` is used only when a fully investigated branch
  ///    provably contains no matching solution.
  ///
  ///  * `Converged` additionally requires the delivered voltage to agree with U
  ///    to within params().voltage_tol (relative).
  ///
  /// `scout_steps` sets the SAMPLING RESOLUTION of the continuation and nothing
  /// else.  Raising it can reveal a pair of crossings closer together than the
  /// previous spacing; it never extends the traced range.  Use `h_max` for that.
  MeniscusSolution solve_at_voltage(Real U, Real h_max,
                                    BranchSide side = BranchSide::Unspecified,
                                    int scout_steps = 14);

  /// Put the solver's BEM into the state described by `sol`, so that any
  /// surface or mesh dump provably belongs to that state rather than to
  /// whatever the last internal iteration happened to leave behind.
  void realize(const MeniscusSolution& sol);

  /// Trace the branch h = h_min ... h_max.  The onset voltage is the maximum of
  /// `voltage` over the returned branch; everything beyond it is unstable.
  std::vector<MeniscusSolution> continuation(Real h_min, Real h_max, int n_steps);

  /// CANDIDATE turning point of the traced branch: the maximum of U(h) as seen
  /// on a finite set of sampled branch points.
  ///
  /// What it is: a discrete maximum, refined parabolically through its two
  /// neighbours, reported only when it is a genuine INTERIOR maximum -- at
  /// least three converged points, strictly rising before and strictly falling
  /// after.  A single point, or a monotone branch, has no turning point.
  ///
  /// What it is NOT, and what may NOT be derived from it:
  ///  * dynamic stability, or the loss of it.  No eigenvalue analysis is
  ///    performed.  Whether the fold coincides with the stability limit depends
  ///    on the control parameter held fixed and on the perturbation mode
  ///    considered; neither is examined here.
  ///  * an emission onset.  In the pure ionic regime the evaporation rate is a
  ///    smooth exponential of the field with no threshold; what experiments
  ///    call onset is a detection limit.
  ///  * the transition to the cone-jet regime, which involves the flow and is
  ///    not determined by a static model.
  ///
  /// Being a discrete maximum, its position and value depend on the sampling of
  /// the continuation.  Refine the continuation before quoting the number.
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
