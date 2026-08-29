#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/device_geometry.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P1 boundary mesher -- automatic, axisymmetric, meridian half-plane (r, z)
// ===========================================================================
//
// Turns a DeviceGeometry into a connected discretisation of its boundary.  It
// is a STANDALONE component: it depends on device_geometry.hpp and nothing
// else, it knows no solver, and it deliberately offers no conversion into the
// P0 BEM `es::Mesh`.  The P0 BEM keeps its own meshes; switching it over is a
// separate, explicit step.
//
// NOT a volume mesh.  Only the boundary is discretised, which is what a BEM
// needs.  The volume-mesher decision is documented in docs/07_mesher_decision.md
// and is due before P3, because the coupled liquid model needs a volume
// discretisation.
//
// AXISYMMETRIC, NOT PLANAR.  Every element is a segment in the meridian
// half-plane; the object it represents is the surface swept by rotating that
// segment about r = 0.  The two measures never share a name:
//
//   meridian_length  [m]    length of the segment in the (r,z) plane
//   revolved_area    [m^2]  2*pi * int r ds  -- the actual surface area
//
// Elements lying on r = 0 are SYMMETRY elements, not ring elements: they sweep
// no area at all.  They are tagged ElementKind::AxisSymmetry and their
// revolved_area is exactly zero.  Treating one as a ring, or dividing by its
// radius, is the mistake this tagging exists to prevent.

// ---------------------------------------------------------------------------

enum class ElementKind {
  Ring = 0,      ///< sweeps a surface of revolution of positive area
  AxisSymmetry,  ///< both endpoints exactly on r = 0; sweeps nothing
};
const char* to_string(ElementKind k);

/// Long, unambiguous label for a boundary identifier.  Used in reports and
/// figures.  In particular the plane at z = 0 is named for what it is -- the
/// initial flat liquid surface, not a computed meniscus.
const char* boundary_long_name(BoundaryId b);

// ---------------------------------------------------------------------------

struct MeshNode {
  Vec2 p;
  bool on_axis{false};    ///< exactly r = 0
  bool is_corner{false};  ///< coincides with a vertex of a boundary polyline
  int feature{-1};        ///< index into DeviceGeometry::features(), -1 if none
};

/// One boundary element.  Carries everything the axisymmetric contract needs.
struct BoundaryElement {
  // --- identity -------------------------------------------------------------
  BoundaryId id{BoundaryId::OpenBoundary};
  int curve{-1};    ///< index into DeviceGeometry::boundaries(); names are unique
  int segment{-1};  ///< index of the straight geometric segment inside that curve

  // --- endpoints ------------------------------------------------------------
  int node_a{-1}, node_b{-1};
  Vec2 a{}, b{};

  // --- orientation ----------------------------------------------------------
  /// Unit tangent a -> b.
  Vec2 tangent{};
  /// Unit normal, equal to perp(tangent) = (t_z, -t_r).  For an interface it
  /// points AWAY from side_a and INTO side_b, which is verified geometrically
  /// rather than inherited from the order of the input polyline.  For an
  /// AxisSymmetry element it is (-1, 0), points out of the meridian half-plane
  /// and carries no material meaning.
  Vec2 normal{};

  // --- adjacency ------------------------------------------------------------
  Region side_a{Region::Vacuum};
  Region side_b{Region::Outside};

  // --- measures -------------------------------------------------------------
  Real meridian_length{0.0};
  Real revolved_area{0.0};  ///< pi*(r_a + r_b)*L, exact for a straight segment

  ElementKind kind{ElementKind::Ring};

  Vec2 midpoint() const { return 0.5 * (a + b); }
  Real mid_radius() const { return 0.5 * (a.r + b.r); }
  bool is_axis() const { return kind == ElementKind::AxisSymmetry; }
};

// ---------------------------------------------------------------------------
// Size function
// ---------------------------------------------------------------------------
//
// There are NO user-facing mesh parameters -- no h_tip, no h_far.  The element
// size is a pure function of the geometry:
//
//     h(x) = min( min_s [ h_s + G * |x - x_s| ],  h_max )
//
// a minimum over cones seated on a finite set of refinement sources.  A minimum
// of G-Lipschitz functions is G-Lipschitz, so the size field itself bounds the
// growth rate; no post-smoothing is needed and adjacent element sizes cannot
// differ by more than (1 + G/2)/(1 - G/2).
//
// Sources, both seated on geometric points, both scaled by the LOCAL FEATURE
// SIZE there -- the smaller of (a) the shortest boundary segment meeting at
// that point and (b) the distance to the nearest boundary segment not meeting
// there:
//
//   * every named zero-dimensional feature (pinned exit edge, outer land edge,
//     the two extractor aperture edges):  h = lfs / kFeatureDivisions
//   * every remaining polyline vertex:    h = lfs / kCornerDivisions
//
// and one global ceiling h_max = domain diagonal / kDomainDivisions, which is
// what makes the far field coarse.
//
// A sharp edge MAY be refined, and is.  That does not make an electric field
// computed at such an edge a mesh-convergent quantity: the corner field of an
// unrounded edge diverges, so its value tracks the local element size and must
// never later be reported as converged.  Rounding radii are a P3 parameter.

namespace mesher {
/// Element size at a named feature = local feature size / this.
inline constexpr Real kFeatureDivisions = 32.0;
/// Element size at an ordinary polyline vertex = local feature size / this.
inline constexpr Real kCornerDivisions = 8.0;
/// Coarsest element = domain diagonal / this.
inline constexpr Real kDomainDivisions = 40.0;
/// Lipschitz constant of the size field: growth of h per unit distance.
inline constexpr Real kGradation = 0.25;
/// Every straight geometric segment gets at least this many elements, so a
/// short segment can never collapse into a single element.
inline constexpr int kMinElementsPerSegment = 4;
/// Two mesh points closer than this fraction of the domain diagonal are the
/// same node.  Orders of magnitude below the smallest element this geometry
/// produces, so it can only ever merge points that are meant to coincide.
inline constexpr Real kNodeSnapRelative = 1e-12;
/// Offset, in units of the element length, used to probe which region lies on
/// either side of an element.
inline constexpr Real kSideProbeRelative = 1e-3;
/// Bound on the size ratio of two elements sharing a node.  The size field
/// gives (1+G/2)/(1-G/2) = 1.286 for G = 0.25; the check allows a margin for
/// the integer rounding of the element count per segment.
inline constexpr Real kMaxNeighbourRatio = 1.5;
}  // namespace mesher

class SizeField {
 public:
  struct Source {
    Vec2 x;
    Real h{0.0};
    Real local_feature_size{0.0};
    bool is_named_feature{false};
    std::string origin;  ///< human-readable provenance
  };

  /// `size_scale` multiplies every target size uniformly.  It is a
  /// CONVERGENCE-STUDY knob and nothing else: 1.0 is the automatic mesh, and
  /// production code must not pass anything else.  It exists because a mesh
  /// with no user parameters cannot otherwise be refined to demonstrate
  /// discretisation convergence.  It changes no shape, no source and no
  /// gradation -- only the overall element size.
  static SizeField from_geometry(const DeviceGeometry& g, Real size_scale = 1.0);

  /// Target element size at x [m].  Strictly positive.
  Real operator()(Vec2 x) const;

  Real h_max() const { return h_max_; }
  Real h_min() const { return h_min_; }  ///< smallest source size
  Real gradation() const { return mesher::kGradation; }
  const std::vector<Source>& sources() const { return sources_; }
  void print(std::FILE* out) const;

 private:
  std::vector<Source> sources_;
  Real h_max_{0.0};
  Real h_min_{0.0};
};

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

struct MeshCheck {
  std::string name;
  bool passed{false};
  std::string detail;
};

struct MeshReport {
  std::vector<MeshCheck> checks;
  bool all_passed() const;
  int failures() const;
  void print(std::FILE* out) const;
};

/// Element-length statistics for one group of elements.
struct LengthStats {
  int count{0};
  Real min{0.0}, median{0.0}, max{0.0};
  Real total_meridian_length{0.0};
  Real total_revolved_area{0.0};
};

// ---------------------------------------------------------------------------

class BoundaryMesh {
 public:
  /// Deterministic: identical DeviceParameters and size_scale give a bitwise
  /// identical mesh.  Throws std::runtime_error on a degenerate element (zero
  /// length, negative radius) or an ambiguous material assignment.
  /// See SizeField::from_geometry for what `size_scale` is and is not for.
  static BoundaryMesh generate(const DeviceGeometry& g, Real size_scale = 1.0);

  /// The uniform size multiplier this mesh was generated with; 1.0 = automatic.
  Real size_scale() const { return size_scale_; }

  const std::vector<MeshNode>& nodes() const { return nodes_; }
  const std::vector<BoundaryElement>& elements() const { return elements_; }
  const SizeField& size_field() const { return size_; }

  std::vector<const BoundaryElement*> elements_with(BoundaryId id) const;
  std::vector<const BoundaryElement*> elements_of_curve(int curve) const;

  LengthStats stats_of(BoundaryId id) const;
  LengthStats stats_of_curve(int curve) const;
  LengthStats stats_total() const;

  /// Largest size ratio between two elements sharing a node.
  Real max_neighbour_ratio() const;

  /// Runs every topological, geometric and reproducibility check.  `g` must be
  /// the geometry this mesh was generated from.
  MeshReport validate(const DeviceGeometry& g) const;

  void print(std::FILE* out, const DeviceGeometry& g) const;
  /// Writes mesh_nodes.csv, mesh_elements.csv, mesh_boundaries.csv and
  /// mesh_size_field.csv into `dir`.
  void write_csv(const std::string& dir, const DeviceGeometry& g) const;

 private:
  std::vector<MeshNode> nodes_;
  std::vector<BoundaryElement> elements_;
  SizeField size_;
  Real size_scale_{1.0};
};

}  // namespace es
