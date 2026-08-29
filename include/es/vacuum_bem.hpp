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
//   and the bore wall and the liquid inlet lie inside the conductor union.
//   The selection below is made from region and boundary tags only -- never
//   from coordinates and never from element order.

/// Which conductor an accepted element belongs to, and why it was accepted.
struct VacuumPanel {
  int mesh_element{-1};   ///< index into BoundaryMesh::elements()
  int bem_element{-1};    ///< index into the produced es::Mesh
  BoundaryId boundary{BoundaryId::OpenBoundary};
  Region material{Region::Vacuum};  ///< the non-vacuum side
  Tag tag{Tag::Other};
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

  /// Meridian-length-weighted totals, for cross-checking against the mesher.
  Real revolved_area_emitter{0.0};
  Real revolved_area_free_surface{0.0};
  Real revolved_area_extractor{0.0};

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
/// The emitter panels form an arc that is closed through the symmetry axis at
/// one end and OPEN at the other, where the domain floor cuts the shank.  That
/// open end is a truncation of a feed-through, not a surface; it is left open
/// rather than capped with the domain edge.  The cavity behind it is a tube of
/// radius phi_3/2 and length emitter_height, so the interior is screened many
/// e-foldings over; test_vacuum_bem checks that numerically, which is what
/// makes En() = sigma/eps0 usable on the tip.
Mesh vacuum_bem_mesh(const BoundaryMesh& bm, const DeviceGeometry& g,
                     VacuumSelectionReport* report = nullptr);

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
/// Reported twice: over the whole conductor, and over the part outside the
/// marked edge zones.  The first is dominated by the integrable density
/// singularities at the unrounded edges and at the open truncation rim, where a
/// piecewise-constant density cannot represent sigma at all; the second is the
/// error of the solution in the region where the answer is meant to be used.
struct PotentialResidual {
  Real max_emitter{0.0};      ///< [V], everywhere
  Real max_extractor{0.0};
  Real max_emitter_clear{0.0};    ///< [V], outside every marked edge zone
  Real max_extractor_clear{0.0};
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
};
PotentialResidual potential_residual(const BemSolver& bem, const std::vector<struct EdgeZone>& zones);

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
/// A second, different kind of place has to be marked as well.  The emitter
/// panels form an arc that is OPEN where the domain floor cuts the shank.  In a
/// single-layer formulation the density at the free edge of an open sheet
/// diverges like 1/sqrt(distance) -- an artefact of truncating a feed-through,
/// not a property of the device.  It is integrable, so charges and capacitances
/// converge, but |E_n| next to it does not, and it is by far the largest value
/// on the emitter.  It is marked exactly like a sharp edge, and for the same
/// reason: nothing there may be reported as a field.
enum class EdgeKind {
  SharpFeature = 0,  ///< unrounded geometric edge of the device
  TruncationEnd,     ///< open end of a modelled conductor arc
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

}  // namespace es
