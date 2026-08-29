#pragma once
#include <string>
#include <vector>

#include "es/types.hpp"

namespace es {

/// What a boundary element belongs to.  Drives boundary conditions, the
/// meniscus solver (which only moves FreeSurface elements) and the beam code
/// (which reports interception per tag).
enum class Tag : int {
  Emitter = 0,     ///< metal of the emitter / capillary, held at V_emitter
  FreeSurface = 1, ///< liquid meniscus, equipotential with the emitter
  Extractor = 2,   ///< extractor / accelerator electrode
  Collector = 3,   ///< downstream target or beam dump
  Other = 4,
};

const char* tag_name(Tag t);

/// Straight constant-density boundary element in the meridian half-plane.
/// The generating curve is traversed counter-clockwise in (r,z), so that the
/// outward normal (into the vacuum) is perp(tangent) = (t_z, -t_r).
struct Element {
  Vec2 a{}, b{};      ///< endpoints
  Vec2 mid{};         ///< collocation point (midpoint)
  Vec2 tangent{};     ///< unit tangent a->b
  Vec2 normal{};      ///< unit outward normal
  Real len{0};        ///< meridian length [m]
  Real area{0};       ///< 2*pi*r_mid*len, the actual surface area [m^2]
  Tag tag{Tag::Other};
  Real potential{0};  ///< Dirichlet value [V]
  int body{0};        ///< index of the closed contour this element came from
};

class Mesh {
 public:
  std::vector<Element> elems;

  Index size() const { return static_cast<Index>(elems.size()); }
  void clear() { elems.clear(); nodes_pending_.clear(); }

  /// --- contour construction -------------------------------------------------
  /// Start a new closed contour.  Nodes are accumulated and turned into
  /// elements by end_body().
  void begin_body(Tag tag, Real potential);
  void add_node(Vec2 p);
  /// Graded straight run from the current node to `to`, with element size
  /// varying linearly along the run from h_start to h_end.
  void line_to(Vec2 to, Real h_start, Real h_end);
  void line_to(Vec2 to, Real h) { line_to(to, h, h); }
  /// Graded circular arc from the current node, about `center`, sweeping to
  /// angle `a_end` (radians, measured from +r axis in the (r,z) plane).
  void arc_to(Vec2 center, Real radius, Real a_start, Real a_end, Real h);
  /// Close the contour (optionally back to the first node) and emit elements.
  void end_body(bool close_loop);

  /// Recompute midpoints, normals, lengths.  Called by end_body(); call again
  /// after moving nodes by hand (meniscus iteration).
  void finalize();

  /// Diagnostics
  Real total_area() const;
  /// `header` is written as leading '#'-prefixed comment lines, so that a
  /// stray file still says which run and which state it belongs to.
  void write_csv(const std::string& path, const std::string& header = {}) const;

 private:
  std::vector<Vec2> nodes_pending_;
  Tag tag_pending_{Tag::Other};
  Real pot_pending_{0};
  int body_counter_{0};
};

/// Number of nodes for a linearly varying target size h0 -> h1 over length L,
/// and the corresponding node positions in [0,1].  Exposed for testing.
std::vector<Real> graded_parameters(Real length, Real h0, Real h1);

// ---------------------------------------------------------------------------
// Ready-made geometries
// ---------------------------------------------------------------------------

/// Isolated conducting sphere -- BEM verification case (C = 4 pi eps0 R).
Mesh make_sphere(Real R, Real V, int n_elem);

/// Prolate spheroid, semi-axis `a` along z, `b` transverse, centred at origin.
/// Verification case with an analytic tip field -- the closest analytic
/// stand-in for a sharpened needle.
Mesh make_prolate_spheroid(Real a, Real b, Real V, int n_elem);

/// Analytic reference values for the prolate spheroid at potential V.
Real spheroid_capacitance(Real a, Real b);
Real spheroid_tip_field(Real a, Real b, Real V);

struct CapillaryParams {
  Real r_inner{1.0e-5};   ///< bore radius [m]
  Real r_outer{2.0e-5};   ///< outer radius [m]
  Real rim_radius{0};     ///< rounding of the rim; 0 -> (r_outer-r_inner)/2
  Real shank_length{5.0e-4};
  Real z_tip{0.0};        ///< axial position of the rim
  Real h_tip{0};          ///< element size at the rim; 0 -> r_inner/12
  Real h_far{0};          ///< element size far from the tip; 0 -> shank/25
  Real potential{0.0};
};
/// Hollow capillary (annular cross-section) with a rounded rim, axis along z,
/// opening toward +z.  The bore is left open -- add the meniscus separately.
Mesh make_capillary(const CapillaryParams& p);

struct NeedleParams {
  Real tip_radius{1.0e-6};       ///< radius of curvature of the apex [m]
  Real half_angle{10.0 * 3.14159265358979323846 / 180.0};  ///< cone half-angle [rad]
  Real shank_radius{2.5e-4};     ///< radius where the cone meets the shank
  Real length{2.0e-3};           ///< total axial extent below the apex
  Real z_tip{0.0};               ///< apex position
  Real h_tip{0};                 ///< element size at the apex; 0 -> tip_radius/8
  Real h_far{0};
  Real potential{0.0};
};
/// Solid externally-wetted needle: spherical cap + cone + cylindrical shank.
Mesh make_needle(const NeedleParams& p);

struct OpenCapillaryParams {
  Real r_bore{1.0e-5};     ///< bore radius = contact-line pinning radius [m]
  Real r_outer{2.0e-5};    ///< outer radius [m]
  Real shank_length{5.0e-4};
  Real z_rim{0.0};         ///< axial position of the rim face
  Real h_rim{0};           ///< element size near the rim; 0 -> r_bore/12
  Real h_far{0};
  Real potential{0.0};
};
/// Capillary with a SHARP inner edge, meshed as an OPEN chain that runs
///     (0, z_bot) -> (r_outer, z_bot) -> (r_outer, z_rim) -> (r_bore, z_rim)
/// i.e. fictitious bottom cap, outer wall, annular rim face -- stopping exactly
/// at the contact line.  The meniscus solver appends the free surface, which
/// runs on to the apex on the axis and closes the body.  The bore is never
/// meshed: it is full of liquid, and liquid plus metal form one conductor.
///
/// A sharp edge is deliberate here: it makes the pinning radius unambiguous,
/// which is what the Young-Laplace boundary condition needs.  Use
/// make_capillary() instead when you want a rounded rim and no free surface.
Mesh make_capillary_open(const OpenCapillaryParams& p);

struct ExtractorParams {
  Real aperture_radius{2.0e-4};  ///< hole radius [m]
  Real outer_radius{3.0e-3};     ///< outer radius of the modelled plate [m]
  Real thickness{1.0e-4};
  Real z_plate{5.0e-4};          ///< axial position of the downstream face
  Real edge_radius{0};           ///< rounding of the aperture lip; 0 -> thickness/2
  Real h_edge{0};
  Real h_far{0};
  Real potential{-1000.0};
};
/// Annular extractor electrode with a rounded aperture lip.
Mesh make_extractor(const ExtractorParams& p);

/// Concatenate meshes (bodies keep distinct `body` indices).
Mesh merge(const std::vector<Mesh>& parts);

}  // namespace es
