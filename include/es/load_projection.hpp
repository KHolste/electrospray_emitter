#pragma once
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "es/electrocapillary.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P0 -- the load projection, made visible and testable
// ===========================================================================
//
// WHY THIS FILE EXISTS
//
// P3b hands the capillary solver a surface load.  Three different objects were
// all called "the load" in that phase, and only one of them is what the ODE
// actually integrates:
//
//   1. the NODAL load        p_M = eps0 E_n^2 / 2 at the mesh nodes.  It is the
//      raw field quantity.  At the contact line it does not converge, and it is
//      used for nothing.
//   2. the SEGMENT load      one constant per surface segment, equal to the
//      integrated normal force of that segment divided by its revolved area.
//      This is the conservative projection and the force bookkeeping.  It is a
//      STAIRCASE and must never be an ODE right-hand side.
//   3. the HANDED load       the continuous reconstruction p(tau) = G'/A' that
//      solve_coupled() actually passes to solve_capillary_meniscus().
//
// The three are now separate, named types, all of them measurable from outside,
// so the claim "the load handed over is continuous and carries the integrated
// Maxwell force" can be checked instead of read.
//
// WHAT IS NOT CLAIMED HERE.  Nothing in this file computes a field, a shape or
// a physical result.  It audits a projection, and the manufactured loads below
// are prescribed functions -- they are not physics and must never be reported
// as one.

// ---------------------------------------------------------------------------
// The load actually handed to the capillary solver
// ---------------------------------------------------------------------------
//
// This was an anonymous-namespace type inside src/electrocapillary.cpp, which
// is why it could not be shown to be continuous.  It is unchanged in substance:
// the bins, the Fritsch-Carlson slopes and at() are the same arithmetic.
//
// THE CONSTRUCTION.  The segment data are binned onto a fixed uniform grid in
// the normalised arclength tau = s/L, giving a bin force and a bin area.  Their
// running sums are the CUMULATIVE force G and area A, which are monotone by
// construction.  Both are interpolated by a monotone cubic and the load is
//
//     p(tau) = G'(tau) / A'(tau) .
//
// Two properties follow, and both are measured in tests/test_load_projection:
//
//   * CONTINUITY.  A cubic Hermite has the prescribed slope at both ends of
//     every interval, and neighbouring intervals share it, so G' and A' are
//     continuous across the bin boundaries and so is their quotient wherever
//     A' > 0.  This is the whole reason the reconstruction exists: the
//     staircase made the adaptive integrator of the capillary solver run to its
//     refinement cap at every voltage.
//   * CONSERVATION.  Over one bin,  int p A' dtau = int G' dtau = G(b+1)-G(b),
//     which is exactly the force that bin carried.  The integrated Maxwell
//     force therefore survives the reconstruction bin by bin -- with respect to
//     the reconstructed area measure A'.  Against the TRUE surface element
//     2 pi r ds the same integral agrees only to the interpolation order, and
//     the audit below reports BOTH numbers rather than the flattering one.
struct ProjectedLoad {
  static constexpr int kBins = 128;

  std::vector<Real> p;     ///< bin pressure [Pa] -- the conservative bin mean
  std::vector<Real> area;  ///< bin revolved area [m^2]

  ProjectedLoad();

  /// Bin the segment data of `L`.  No clipping, no exclusion zone, no maximum:
  /// every segment contributes its full force to the bins it overlaps, split by
  /// the fraction of its tau span.
  static ProjectedLoad from(const MaxwellLoad& L);

  /// THE HANDED LOAD.  p(tau) = G'/A' from the monotone cubic interpolants.
  /// Continuous in tau; this is the function solve_coupled() passes on.
  Real at(Real tau) const;

  /// The staircase the bins are, evaluated at the same tau.  It is NOT what is
  /// handed over; it exists so that a test can show the two differ and by how
  /// much, and so that a figure can draw both.
  Real bin_pressure_at(Real tau) const;

  /// sum p_b * area_b [N] -- the force content of the bins.
  Real integrated_force() const;
  /// sum area_b [m^2].
  Real total_area() const;
  /// Number of bins that received no segment area at all.  Must be zero for a
  /// surface whose segments cover [0,1]; a non-zero count is a coverage gap and
  /// is reported, never smoothed away.
  Index empty_bins() const;

  static Real difference(const ProjectedLoad& a, const ProjectedLoad& b);
  static ProjectedLoad blend(const ProjectedLoad& old_load, const ProjectedLoad& fresh, Real w);

  /// True when the table carries no load at all.  Then the problem IS the P3a
  /// problem and is handed to the capillary solver as such -- same right-hand
  /// side, same requested accuracy, same answer bit for bit.
  bool is_zero() const;

  /// The cumulative force and area the reconstruction is built on, exposed so
  /// that a figure can show them next to the segment means instead of a
  /// summary of them.  Size kBins + 1, starting at zero.
  const std::vector<Real>& cumulative_force() const;
  const std::vector<Real>& cumulative_area() const;

 private:
  mutable std::vector<Real> cum_force_, cum_area_, slope_force_, slope_area_;
  void build_cumulative() const;
  Real slope_of(const std::vector<Real>& y, const std::vector<Real>& d, int b, Real u) const;
};

// ---------------------------------------------------------------------------
// Manufactured loads: a prescribed pressure on a prescribed surface
// ---------------------------------------------------------------------------
//
// The projection is a piece of quadrature and bookkeeping.  To test it, the
// input must be one whose integral is known in closed form -- so the load is
// PRESCRIBED and no field is solved.  The segment quadrature used is literally
// the one maxwell_load() uses, because the same helper assembles it; otherwise
// the audit would test a copy of the code rather than the code.
//
// THE EDGE NODE OF A SINGULAR LOAD.  For p = C d^beta with -1 < beta < 0 the
// pointwise value at the contact line is infinite, and no discretisation can
// carry it.  The manufactured input therefore gives that ONE node the local
// mean over the last half segment,
//
//     p_edge = (2/h) int_0^{h/2} C d^beta dd = C (h/2)^beta / (1 + beta) ,
//
// which is finite, well defined and stated.  It mirrors what the recovered FEM
// field does there -- a cell average, not a point value -- and it is a declared
// property of the MANUFACTURED INPUT, not a regularisation of the projection.
// Nothing else is regularised anywhere.

/// Prescribed pressure as a function of the distance d from the contact line,
/// measured along the surface [m] -> [Pa].
using PrescribedPressure = std::function<Real(Real)>;

/// Build a MaxwellLoad on `fs` whose nodes sit at the given radii (ascending,
/// starting at 0 and ending at the contact radius) and whose nodal pressure is
/// `p_of_d`, with the edge node treated as documented above.
///
/// `node_En` is filled with the E_n that would produce this pressure with the
/// positive sign, so that the record is self-consistent; it is not used by the
/// projection.
MaxwellLoad manufactured_load(const FreeSurface& fs, const std::vector<Real>& node_r,
                              const PrescribedPressure& p_of_d, Real gamma_over_a);

/// Node radii of a uniform-in-radius surface discretisation with `n_segments`
/// segments, which is the layout the P2c mesh row has inside the bore.
std::vector<Real> uniform_radius_nodes(Real contact_radius, Index n_segments);

/// Closed-form integrated normal force of p = C d^beta on a FLAT disc of radius
/// a, with d = a - r:
///
///     F = int_0^a C (a-r)^beta 2 pi r dr = 2 pi C a^(2+beta)
///                                          / ((1+beta)(2+beta)) .
///
/// Valid for beta > -1; the integral does not exist otherwise.
Real flat_disc_power_law_force(Real C, Real beta, Real a);

/// Closed-form integrated normal force of p = p0 (1 + (r/a)^2) on a flat disc:
/// F = (3/2) pi p0 a^2.  A smooth load with no edge structure at all.
Real flat_disc_smooth_force(Real p0, Real a);

// ---------------------------------------------------------------------------
// The audit
// ---------------------------------------------------------------------------

struct LoadProjectionAudit {
  std::string tag;
  Index n_nodes{0}, n_segments{0};

  // --- force, at every stage of the chain ----------------------------------
  /// Known in closed form, or NaN when there is none (a solved field).
  Real analytic_force{0};
  /// Sum of the segment forces: the conservative projection of the nodal load.
  Real segment_force{0};
  /// Force content of the bins: sum p_b area_b.
  Real bin_force{0};
  /// Integral of the HANDED load against the reconstructed area measure A'.
  /// This is the quantity the reconstruction preserves exactly.
  Real handed_force_reconstructed{0};
  /// Integral of the HANDED load against the TRUE surface element 2 pi r ds.
  /// This is what the capillary solver really integrates, and it agrees only to
  /// the interpolation order.  Both are reported; neither is hidden.
  Real handed_force_true{0};

  Real error_segment_vs_analytic{0};   ///< relative, NaN without an analytic value
  Real error_bin_vs_segment{0};
  Real error_handed_reconstructed{0};  ///< must be at machine level
  Real error_handed_true{0};

  // --- continuity -----------------------------------------------------------
  //
  // HOW CONTINUITY IS MEASURED, AND WHY NOT BY A THRESHOLD.  The load is
  // sampled a probe offset delta on either side of every bin boundary.  A
  // continuous function with a kink gives a difference proportional to delta;
  // a staircase gives one that does not depend on delta at all.  So the jump is
  // measured at delta and again at delta/10, and the DECAY of the two is the
  // statement: about 0.1 for a continuous load, about 1 for a staircase.  No
  // tolerance has to be invented for it.
  static constexpr Real kProbeOffset = 1.0e-6;  ///< in units of the bin width
  /// Largest jump of the HANDED load across a bin boundary, at the probe offset
  /// and at a tenth of it [Pa].
  Real handed_max_jump{0}, handed_max_jump_tenth{0};
  /// jump(delta/10) / jump(delta).  About 0.1 for a continuous load.
  Real handed_jump_decay{0};
  /// Largest jump of the STAIRCASE across the same boundaries [Pa], and its
  /// decay.  It is the control: without it a small number above would mean
  /// nothing.
  Real staircase_max_jump{0}, staircase_max_jump_tenth{0};
  Real staircase_jump_decay{0};
  /// Spread of the load over the surface [Pa], the scale both jumps are judged
  /// against.
  Real load_span{0};
  Real handed_jump_ratio{0};      ///< handed_max_jump / load_span
  Real staircase_jump_ratio{0};   ///< staircase_max_jump / load_span

  // --- edge and coverage ----------------------------------------------------
  Real tau_first{0}, tau_last{0};   ///< tau range the segments actually cover
  Real segment_area{0}, surface_area{0}, bin_area{0};
  Real area_gap{0};                 ///< |bin_area - segment_area| / segment_area
  Index empty_bins{0};
  /// Force carried by the last bin, i.e. by the surface nearest the contact
  /// line, as a fraction of the total.  A hidden exclusion zone would make it
  /// zero; it is measured so that it cannot be.
  Real last_bin_force_fraction{0};
  /// Largest handed pressure and largest bin pressure [Pa].  A cap anywhere in
  /// the chain would show as the first being pinned below the second.
  Real max_handed_pressure{0}, max_bin_pressure{0};
  Real max_node_pressure{0};        ///< the raw nodal peak, for the contrast

  void print(std::FILE* out) const;
};

/// Audit a load and its projection.  `analytic_force` may be NaN.
LoadProjectionAudit audit_projection(const MaxwellLoad& raw, const FreeSurface& fs,
                                     const std::string& tag, Real analytic_force);

// ---------------------------------------------------------------------------
// Richardson extrapolation and the discretisation verdict
// ---------------------------------------------------------------------------
//
// The mesh levels of this code scale the whole size field by 2^(-level/2), so
// consecutive levels have a refinement ratio of exactly sqrt(2) in h.  With
// three values f1, f2, f3 on such a sequence,
//
//     p     = ln((f1-f2)/(f2-f3)) / ln(ratio)          observed order
//     f_inf = f3 + (f3-f2)/(ratio^p - 1)               extrapolated limit
//
// and the estimated relative discretisation error of the finest value is
// |f_inf - f3| / |f_inf|.
//
// WHEN THIS IS NOT USABLE.  If the three differences do not have the same sign,
// or one of them vanishes, the sequence is not in its asymptotic range and no
// order can be observed.  The estimate is then NOT produced; `usable` is false
// and the honest fallback is the relative change between the two finest levels,
// which is a change and not an error bound.  That distinction is kept in the
// field names.
struct RichardsonEstimate {
  int n_levels{0};
  Real ratio{0};                   ///< h_coarse / h_fine between consecutive levels
  Real observed_order{0};          ///< NaN when not observable
  Real extrapolated{0};            ///< NaN when not observable
  Real relative_error_finest{0};   ///< NaN when not observable
  Real last_relative_change{0};    ///< always available; a change, not an error
  bool monotone{false};
  bool usable{false};
  std::string note;
};

/// `values` in order of increasing refinement (coarse first).  At least three.
RichardsonEstimate richardson(const std::vector<Real>& values, Real ratio);

/// The refinement ratio in h between two consecutive mesh levels of this code.
inline constexpr Real kMeshLevelRatio = 1.4142135623730951;  // sqrt(2)

/// The target fixed BEFORE the measurement for P0: the estimated
/// discretisation error of the finest level must be below one per cent.
inline constexpr Real kDiscretizationTarget = 1.0e-2;

enum class DiscretizationVerdict {
  NotAttempted = 0,
  Converged,                 ///< estimated error below kDiscretizationTarget
  DiscretizationNotConverged,///< an estimate exists and misses the target
  NotInAsymptoticRange,      ///< no order observable; only a change is known
  InsufficientLevels,
};
const char* to_string(DiscretizationVerdict v);
const char* explain(DiscretizationVerdict v);

/// Verdict from an estimate, against kDiscretizationTarget.
DiscretizationVerdict verdict_of(const RichardsonEstimate& e);

// ---------------------------------------------------------------------------
// Why the exclusion-halving criterion of P3b cannot be met -- in closed form
// ---------------------------------------------------------------------------
//
// P3b declared edge_gate::kTolExclusion: halving the exclusion distance must
// change the integrated force by less than five per cent.  Concretely it forms
//
//     |F(d0/2) - F(d0)| / |F(d0/2)| ,   d0 = edge_gate::kExclusionMid = 0.05 a ,
//
// on the finest mesh.  It measured 11.3 % for the flat shape, reported the
// miss, and let a different criterion decide.  P0 has to settle that instead of
// leaving two criteria standing.
//
// THE MATHEMATICS.  Near the contact line the load behaves as p = C d^beta with
// beta > -1, and the revolved radius tends to a, so the force inside a distance
// d0 of the edge is
//
//     F_total - F(d0) = 2 pi a C d0^(1+beta) / (1+beta) + O(d0^(2+beta)) .
//
// Halving d0 therefore changes the integrated force by
//
//     [F(d0/2) - F(d0)] / F(d0/2) = K d0^(1+beta) (1 - 2^-(1+beta)) / F(d0/2) ,
//
// which is a property of the SURFACE AND THE EXPONENT and does not contain the
// mesh size at all.  It has a finite NON-ZERO limit under refinement.  A
// criterion that demands it become small is therefore not a convergence test;
// it asks a geometric quantity to vanish, and for beta near -0.44 with
// d0 = 0.05 a its limit is of the order of ten per cent -- which is the number
// P3b reported and could not have avoided.
//
// The function below is that closed-form limit.  tests/test_load_projection.cpp
// pins it against the measured value of a manufactured load with a KNOWN beta,
// and pins that the measured value stops changing under refinement while the
// discretisation error of the same load keeps falling.  That is what makes the
// replacement of the criterion a derivation rather than a swap.
//
// The bound itself is not touched.  It is still declared, still measured and
// still reported as not met.
struct ExclusionHalvingLimit {
  Real beta{0};
  Real d0_over_a{0};
  /// The closed-form limit of [F(d0/2) - F(d0)] / F_total for p = C d^beta on
  /// a flat disc of radius a.  Independent of C.
  Real limit_change{0};
  /// The same quantity measured on a load, for comparison.
  Real measured_change{0};
  Real agreement{0};   ///< |measured - limit| / limit
};

/// Closed-form limit for the flat disc.  With the exact force
/// F(d0) = 2 pi C [a (a^(1+b) - d0^(1+b))/(1+b) - (a^(2+b) - d0^(2+b))/(2+b)]
/// this is evaluated exactly rather than through the leading term alone.
Real exclusion_halving_limit(Real beta, Real d0_over_a);

/// Force of p = C d^beta on a flat disc of radius a beyond a distance d0 of
/// the edge, in closed form [N].
Real flat_disc_power_law_force_beyond(Real C, Real beta, Real a, Real d0);

}  // namespace es
