#pragma once
#include <array>
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/boundary_mesh.hpp"
#include "es/device_geometry.hpp"
#include "es/geometry.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P2a -- adapter from the P1 boundary mesh to the existing vacuum BEM
// ===========================================================================
//
// WHAT THE EXISTING BEM CORE ACTUALLY IS  (read off src/bem.cpp, not assumed)
//
//   * Indirect single-layer formulation.  The unknown is a surface charge
//     density sigma on each element; the potential anywhere is
//         V(x) = sum_j sigma_j \int_{elem j} G(x; x') ds'
//     with G the azimuthally integrated free-space Green's function.
//     Collocation is at element midpoints, so the linear system is
//         A sigma = V_applied,   A_ij = \int_{elem j} G(x_i; x') ds'.
//
//   * BOUNDARY CONDITION AT INFINITY: V -> 0.  It is built into the kernel --
//     G is the free-space Green's function, which decays like 1/|x|.  There is
//     no truncation surface, no far-field box and nothing to impose there.
//     Consequence: applied potentials are absolute against infinity, and the
//     total charge of the system is generally NOT zero.
//
//   * SURFACES: the single-layer representation is well posed on any surface,
//     closed or open; on an open sheet sigma is the SUM of the densities on the
//     two faces, so En() = sigma/eps0 is then not a one-sided field.  Closed
//     conductors are therefore the case in which En() is exact.  This adapter
//     documents, and the P2a tests verify, which of the two applies here.
//
//   * SEVERAL CONDUCTORS: elements are grouped into es::Electrode by their
//     es::Tag (bem.cpp, electrode_of).  solve_basis() computes one unit
//     potential basis vector per electrode, and every voltage combination is a
//     superposition.  Only three electrodes exist, and -- a trap --
//     Tag::Other maps to Electrode::Collector, so an untagged element would
//     silently become a third, grounded conductor.  This adapter never emits
//     Tag::Other or Tag::Collector.
//
// WHAT MUST NOT ENTER THE VACUUM BEM
//
//   The open rectangular computational domain is not a conductor.  Feeding its
//   edges to a Dirichlet solver would impose a grounded box that does not
//   exist.  The symmetry axis is not a surface either (it carries no area),
//   and the bore wall lies inside the conductor union.  The selection below is
//   made from region and boundary tags only -- never from coordinates and
//   never from element order.
//
//   The numerical back closure is the one surface that IS a conductor without
//   being a device part.  It enters the solve, because a conductor with a hole
//   in it is a different mathematical object; it does not enter any result,
//   because it is not a piece of hardware.  The two roles are kept apart by
//   VacuumPanel::numerical and VacuumPanel::evaluable, not by a comment.

/// Which conductor an accepted element belongs to, and why it was accepted.
struct VacuumPanel {
  int mesh_element{-1};   ///< index into BoundaryMesh::elements()
  int bem_element{-1};    ///< index into the produced es::Mesh
  BoundaryId boundary{BoundaryId::OpenBoundary};
  Region material{Region::Vacuum};  ///< the non-vacuum side
  Tag tag{Tag::Other};
  /// Part of the NUMERICAL rearward closure rather than of the device.  It
  /// carries charge and must be solved for -- leaving it out would reopen the
  /// conductor -- but no field value from it is a result.
  bool numerical{false};
  /// May sigma/eps0 on this panel be read as a physical one-sided vacuum field?
  /// False on the numerical closure, on the arbitrary-length shank behind the
  /// taper foot, and on panels inside a marked edge zone.
  bool evaluable{false};
};

/// Full account of the selection: what went in, what stayed out and why.
struct VacuumSelectionReport {
  struct CurveDecision {
    std::string curve;      ///< unique BoundaryCurve name
    BoundaryId boundary{BoundaryId::OpenBoundary};
    Region side_a{Region::Vacuum}, side_b{Region::Outside};
    int n_elements{0};
    bool accepted{false};
    Tag tag{Tag::Other};
    std::string reason;
  };
  std::vector<CurveDecision> curves;
  std::vector<VacuumPanel> panels;
  int n_mesh_elements{0};
  int n_selected{0};
  int n_emitter{0}, n_free_surface{0}, n_extractor{0};
  int n_numerical_closure{0};   ///< of the emitter panels, how many are the cap

  /// Meridian-length-weighted totals, for cross-checking against the mesher.
  Real revolved_area_emitter{0.0};
  Real revolved_area_free_surface{0.0};
  Real revolved_area_extractor{0.0};
  Real revolved_area_numerical_closure{0.0};

  void print(std::FILE* out) const;
  void write_csv(const std::string& path) const;
};

/// Build the vacuum-BEM mesh from the P1 boundary mesh.
///
/// Accepts exactly those elements that separate vacuum from a conductor:
/// emitter solid, extractor solid, or the liquid column (a perfect conductor in
/// P2a).  Everything else is rejected with a recorded reason.  Throws if a
/// vacuum-facing interface carries a boundary identifier that this phase does
/// not know how to charge -- silence would be worse than a stop.
///
/// THE CONDUCTOR MUST BE CLOSED.  Every accepted arc has to close, either onto
/// itself or through the symmetry axis, which for the emitter conductor means
/// that DeviceParameters::emitter_back_length must be set.  An arc
/// that stops in mid-vacuum is an open sheet: the single-layer density there is
/// the SUM over both faces, sigma/eps0 is then not a one-sided vacuum field
/// anywhere on that conductor, and the free edge carries a 1/sqrt(d) density
/// singularity that dominates every peak-field number on it.  Earlier phases
/// left the emitter shank cut open at the domain floor and relied on the cavity
/// behind it being screened; that argument bounds the error but does not remove
/// it, and it produced exactly the artificial field maximum this closure
/// exists to eliminate.  vacuum_bem_mesh() now REFUSES an open conductor
/// instead of documenting one.
Mesh vacuum_bem_mesh(const BoundaryMesh& bm, const DeviceGeometry& g,
                     VacuumSelectionReport* report = nullptr);

/// Endpoints of the conductor contour that exactly one panel touches and that
/// do not lie on r = 0.  An endpoint on the axis is where the surface of
/// revolution closes on itself and is not an edge at all.  Empty means the
/// conductor is closed.  Found from panel connectivity, never from coordinates,
/// so it cannot be fooled by a geometry whose numbers happen to look right.
std::vector<Vec2> open_arc_ends(const Mesh& bem_mesh);

// ---------------------------------------------------------------------------
// Capacitance -- named, not just "C"
// ---------------------------------------------------------------------------
//
// With V = 0 at infinity a two-electrode system has a 2x2 Maxwell capacitance
// matrix, not one capacitance:
//
//     Q_E = c_EE V_E + c_EX V_X
//     Q_X = c_XE V_E + c_XX V_X
//
// c_EE and c_XX are positive coefficients of capacitance, c_EX = c_XE < 0 are
// coefficients of induction.  Derived quantities that are often all called "C":
//
//     mutual capacitance             C_m       = -c_EX
//     emitter self-capacitance       C_E_inf   = c_EE + c_EX
//     two-terminal ratio at a point  Q_E / (V_E - V_X)
//
// The last one equals C_m only when the system is charge neutral.  Every
// number this module reports is labelled with which of these it is.
//
// WHAT DEPENDS ON WHICH DIMENSION
//
//   c_EE is the charge on the emitter conductor at unit potential.  A conductor
//   that is made longer holds more charge, full stop.  It does NOT converge as
//   the emitter is extended rearwards, and it must not be expected to: quoting
//   c_EE without stating emitter_back_length would be inventing a dimension.
//
//   c_EX is the charge induced on the EXTRACTOR by the emitter.  The intention
//   was that this one converges, because the far end of the shank points away
//   from the electrode.  It does not: it grows by about a third per doubling of
//   emitter_back_length over the whole range that can be meshed here.  With no
//   return electrode anywhere in the model, the emitter's charge grows nearly
//   linearly with its length and the electrode intercepts a roughly constant
//   share of it.
//
//   The local field quantities -- E_z at the axial reference point, the potential
//   at fixed points between the electrodes -- move more slowly (a 1/L tail) but
//   likewise do not reach the tolerance fixed in advance.
//
//   So: every one of these numbers is a function of the stated geometry, and the
//   run reports it as such.  See truncation::kTol* below for the tolerances, and
//   results/.../truncation.csv for the measurement they failed.
struct CapacitanceMatrix {
  Real c_EE{0.0}, c_EX{0.0}, c_XE{0.0}, c_XX{0.0};  ///< [F]

  Real mutual() const { return -0.5 * (c_EX + c_XE); }
  Real emitter_to_infinity() const { return c_EE + c_EX; }
  Real extractor_to_infinity() const { return c_XX + c_XE; }
  /// |c_EX - c_XE| / max(|c_EX|, |c_XE|); reciprocity demands this be ~0.
  Real reciprocity_error() const;
};

/// Maxwell coefficients from the unit-potential basis.  Calls solve_basis() if
/// needed; leaves the solver's active solution untouched apart from that.
CapacitanceMatrix maxwell_capacitance(BemSolver& bem);

// ---------------------------------------------------------------------------
// Truncation tolerances -- fixed BEFORE the measurement, not after it
// ---------------------------------------------------------------------------
//
// A quantity counts as independent of emitter_back_length when doubling that
// length changes it by less than the bound below.  The bounds were chosen from
// what the numbers are for and what else limits them, before the study was run:
//
//   * the mesh discretisation error of the reference level is a few times 1e-5,
//     so a bound at 1e-3 keeps the two error sources an order apart instead of
//     chasing round-off;
//   * the axial reference field is wanted to about a per cent, because that is
//     the accuracy at which a field-dependent emission law is worth evaluating,
//     so 1e-3 is an order of magnitude below what the answer is used for;
//   * the same bound is applied to c_EX, and to the potential at the fixed probe
//     points relative to the applied span |V_E - V_X|.
//
// The study reports pass or fail against these numbers.  It fails.  That is a
// result about the model, not a licence to move the bound.
namespace truncation {
inline constexpr Real kTolEzRef = 1.0e-3;
inline constexpr Real kTolCEX = 1.0e-3;
inline constexpr Real kTolVProbe = 1.0e-3;
}  // namespace truncation

/// Charge on one electrode for a given applied voltage pair [C].
Real electrode_charge(const BemSolver& bem, Electrode e);

/// Name of an electrode in P2a output.  es::tag_name() calls Tag::FreeSurface a
/// "meniscus"; in this phase the plane at z = 0 is precisely not one, so P2a
/// output never uses that name.
const char* electrode_label(Tag t);

// ---------------------------------------------------------------------------
// Residuals
// ---------------------------------------------------------------------------

/// How well the Dirichlet condition is met AWAY from the collocation points.
/// Sampling at the midpoints would return round-off by construction, so the
/// samples sit at t = 1/4 and 3/4 along every element.
///
/// Reported three times, because the three answer different questions:
///
///   * over the whole conductor -- dominated by the integrable density
///     singularities at the unrounded edges, where a piecewise-constant density
///     cannot represent sigma at all;
///   * outside the marked edge zones;
///   * over the EVALUABLE panels only -- the device surfaces in front of the
///     taper foot.  This is the error of the solution where the answer is
///     actually read off, and it is the only one that has to fall under
///     refinement for the reported numbers to mean anything.  The numerical
///     closure is excluded from it: its own discretisation error is a property
///     of a surface that does not exist on the device.
struct PotentialResidual {
  Real max_emitter{0.0};      ///< [V], everywhere
  Real max_extractor{0.0};
  Real max_emitter_clear{0.0};    ///< [V], outside every marked edge zone
  Real max_extractor_clear{0.0};
  Real max_physical{0.0};     ///< [V], over the evaluable panels only
  Real rms_physical{0.0};     ///< [V], likewise -- the one that must fall with h
  Real rms_emitter{0.0};
  Real rms_extractor{0.0};
  Vec2 worst_position{};      ///< where max over both conductors sits
  Real reference_span{0.0};   ///< |V_E - V_X| used for normalising [V]

  Real relative() const {
    const Real m = std::max(max_emitter_clear, max_extractor_clear);
    return reference_span > 0.0 ? m / reference_span : 0.0;
  }
  Real relative_including_edges() const {
    const Real m = std::max(max_emitter, max_extractor);
    return reference_span > 0.0 ? m / reference_span : 0.0;
  }
  Real relative_physical() const {
    return reference_span > 0.0 ? max_physical / reference_span : 0.0;
  }
  /// The maximum sits on whichever element happens to lie just outside an edge
  /// zone, so it is a max over a set that changes with the mesh and need not
  /// fall monotonically.  The root mean square over the evaluable surfaces does,
  /// and it is the quantity to judge refinement by.
  Real relative_rms_physical() const {
    return reference_span > 0.0 ? rms_physical / reference_span : 0.0;
  }
};

/// `evaluable` is the mask from mark_evaluable_panels(); pass an empty vector to
/// leave max_physical at zero.
PotentialResidual potential_residual(const BemSolver& bem,
                                     const std::vector<struct EdgeZone>& zones,
                                     const std::vector<char>& evaluable = {});

// ---------------------------------------------------------------------------
// Sharp edges
// ---------------------------------------------------------------------------
//
// The exit edge, the outer land edge and both aperture edges are unrounded in
// P1.  The field of a re-entrant wedge diverges there, so |E| at such an edge
// is set by the local element size and converges to nothing.  It is marked, not
// reported as a peak field.  The marking radius is purely geometric: a quarter
// of the local feature size the mesher already computed for that edge.

///
/// The rim of the numerical back closure is marked the same way, for a
/// different reason: it is a perfectly ordinary convex conducting edge, but it
/// belongs to a surface that was invented to close the conductor, so its field
/// is not a device field at any mesh density.
///
/// TruncationEnd is the third kind and no longer occurs in a solved problem.
/// It flags the free edge of an OPEN conductor arc, where the single-layer
/// density diverges like 1/sqrt(distance).  It is integrable, so charges and
/// capacitances still converge, but |E_n| beside it does not, and it used to be
/// the largest value anywhere on the emitter.  edge_zones() still detects it,
/// and vacuum_bem_mesh() refuses any mesh in which one appears -- detecting the
/// defect is what proves the closure did its job.
enum class EdgeKind {
  SharpFeature = 0,     ///< unrounded geometric edge of the device
  TruncationEnd,        ///< open end of a conductor arc -- rejected, never solved
  NumericalClosure,     ///< rim of the numerical rearward closure
};
const char* to_string(EdgeKind k);

struct EdgeZone {
  EdgeKind kind{EdgeKind::SharpFeature};
  std::string name;
  Vec2 position;
  Real radius{0.0};              ///< marking radius [m]
  Real local_feature_size{0.0};
};

/// One zone per named feature, plus one per open end of the conductor arcs.
/// `bem_mesh` is the output of vacuum_bem_mesh for the same boundary mesh; the
/// open ends are found from panel connectivity, not from coordinates.  Ends
/// that lie on the symmetry axis are closures of the surface of revolution, not
/// edges, and are not marked.
std::vector<EdgeZone> edge_zones(const DeviceGeometry& g, const BoundaryMesh& bm,
                                 const Mesh& bem_mesh);

/// Is x inside any marked edge zone?
bool in_edge_zone(const std::vector<EdgeZone>& zones, Vec2 x);

/// Largest |E_n| over elements of the given electrode whose midpoint lies
/// OUTSIDE every edge zone.  Returns the element index through `which`.
Real peak_field_outside_edges(const BemSolver& bem, Electrode e,
                              const std::vector<EdgeZone>& zones, Index* which = nullptr);

// ---------------------------------------------------------------------------
// Where a surface field may be read off at all
// ---------------------------------------------------------------------------
//
// After the closure, sigma/eps0 is the one-sided vacuum normal field exactly
// where the panel bounds the closed perfect-conductor region against vacuum.
// That is every device panel -- and it is NOT the closing disc, whose other
// side is the conductor's own hollow interior, nor the shank behind the taper
// foot, whose length was chosen numerically and whose surface field therefore
// answers a question about the closure rather than about the emitter, nor any
// panel inside a marked edge zone.
//
// Fills `report.panels[k].evaluable` and returns the same information indexed by
// BEM element, so a caller can mask a field sweep without re-deriving the rule.

/// Mark which panels may yield a physical field value.  Modifies `report`.
std::vector<char> mark_evaluable_panels(const Mesh& bem_mesh, VacuumSelectionReport& report,
                                        const DeviceGeometry& g,
                                        const std::vector<EdgeZone>& zones);

/// Largest |E_n| over the evaluable panels of one electrode.  `ok` is the mask
/// returned by mark_evaluable_panels().
Real peak_field_evaluable(const BemSolver& bem, Electrode e, const std::vector<char>& ok,
                          Index* which = nullptr);

}  // namespace es
