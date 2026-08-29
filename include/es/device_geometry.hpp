#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/status.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P1 device geometry -- parametric, axisymmetric, meridian half-plane (r, z)
// ===========================================================================
//
// THIS IS NOT A PLANAR 2D GEOMETRY.  Everything lives in the meridian
// half-plane and every area or volume is what results from revolving it about
// r = 0, so every measure carries the Jacobian 2*pi*r.  The vocabulary keeps
// the two apart on purpose:
//
//   meridian_length   [m]    a length in the (r,z) plane
//   meridian_area     [m^2]  an area in the (r,z) plane -- NOT a surface
//   revolved_area     [m^2]  the surface of revolution, 2*pi * int r ds
//   revolved_volume   [m^3]  the volume of revolution, int 2*pi*r dA
//
// Confusing meridian_area with revolved_area is a factor 2*pi*r error, which is
// why they never share a name.  The analytic tests in tests/test_device_geometry
// pin both against a cylinder and a truncated cone.
//
// SI units throughout.  Diameters are inputs because that is how the drawing is
// dimensioned; radii are used internally.
//
// Coordinate convention: the emitter face lies at z = 0 and z increases towards
// the extractor.  The axis of symmetry is r = 0.

// ---------------------------------------------------------------------------

enum class Region {
  Vacuum = 0,
  Liquid,          ///< ionic liquid filling the bore
  EmitterSolid,    ///< the massive capillary body
  ExtractorSolid,  ///< the annular extraction electrode
  Outside,         ///< beyond the computational domain; not a material
};
const char* to_string(Region r);

/// Labelled pieces of the boundary.  No physical boundary condition is attached
/// to any of them in this phase -- these are identifiers, nothing more.
enum class BoundaryId {
  SymmetryAxis = 0,       ///< r = 0; a symmetry plane, not an interface
  EmitterOuterSurface,    ///< cone frustum flank and cylindrical shank
  EmitterTipLand,         ///< annular face at z = 0 between the bore and the flank
  BoreWall,               ///< cylindrical wall of the bore, liquid against solid
  FreeSurfaceReference,   ///< flat plane closing the liquid column at z = 0
  LiquidInlet,            ///< upstream cut through the liquid column
  ExtractorSurface,       ///< aperture wall and the two faces of the electrode
  OpenBoundary,           ///< outer edges of the computational domain
};
const char* to_string(BoundaryId b);

/// Zero-dimensional features of the meridian plane -- circles in three
/// dimensions.  They are the places a mesher has to resolve.
enum class FeatureId {
  PinnedContactEdge = 0,      ///< sharp exit edge; the contact line is pinned here
  EmitterOuterEdge,           ///< outer edge of the tip land
  ExtractorApertureEdgeFront, ///< aperture edge on the emitter-facing side
  ExtractorApertureEdgeBack,
};
const char* to_string(FeatureId f);

// ---------------------------------------------------------------------------

struct DeviceParameters {
  // --- mandatory ------------------------------------------------------------
  Real phi_3{4.0e-5};                       ///< outer diameter at the emitter foot [m]
  Real phi_1{2.0e-5};                       ///< outer diameter of the tip land [m]
  Real phi_2{1.0e-5};                       ///< diameter of the exit bore [m]
  Real emitter_height{6.0e-5};              ///< axial length of the outer taper [m]
  Real extraction_distance{5.0e-4};         ///< tip plane to extractor face [m]
  Real extractor_aperture_diameter{4.0e-4}; ///< [m]
  Real extractor_thickness{1.0e-4};         ///< [m]

  // --- open computational domain -------------------------------------------
  Real domain_radius{3.0e-3};
  Real domain_z_min{-4.0e-4};
  Real domain_z_max{1.5e-3};

  /// Radial extent of the electrode; 0 means "out to the domain boundary",
  /// which is how the drawing shows it.
  Real extractor_outer_radius{0.0};

  // --- reserved for later phases -------------------------------------------
  //
  // Present so that the parameter set does not have to be reshaped later, and
  // rejected if used, so that nothing is pretended.  build() throws
  // NotImplementedInThisPhase when any of these leaves its default.
  struct Reserved {
    Real edge_radius_inner{0.0};      ///< rounding of the exit edge (P3)
    Real edge_radius_outer{0.0};      ///< rounding of the outer tip edge (P3)
    Real contact_angle_deg{0.0};      ///< Young condition instead of pinning (P3)
    Real bore_diameter_at_inlet{0.0}; ///< tapered bore (later)
    bool porous_emitter{false};       ///< Darcy feed instead of Poiseuille (later)
    bool collector_enabled{false};    ///< downstream collector electrode (later)
  } reserved;
};

// ---------------------------------------------------------------------------

/// One labelled piece of boundary, as a meridian polyline.
struct BoundaryCurve {
  BoundaryId id{BoundaryId::OpenBoundary};
  std::string name;                  ///< unique, e.g. "open_boundary.z_max"
  std::vector<Vec2> points;
  Region side_a{Region::Vacuum};     ///< region the normal points away from
  Region side_b{Region::Outside};    ///< region on the other side

  Real meridian_length() const;
  Real revolved_area() const;
};

/// One material region, as a closed meridian contour with optional holes.
struct RegionBody {
  Region region{Region::Vacuum};
  std::vector<Vec2> outer_loop;                 ///< closed, counter-clockwise
  std::vector<std::vector<Vec2>> holes;

  Real meridian_area() const;   ///< area in the (r,z) plane -- NOT a surface
  Real revolved_volume() const; ///< volume of revolution
};

struct NamedFeature {
  FeatureId id{FeatureId::PinnedContactEdge};
  Vec2 position;
};

// ---------------------------------------------------------------------------

class DeviceGeometry {
 public:
  /// Build and validate.  Throws std::runtime_error on an invalid parameter set
  /// and NotImplementedInThisPhase if a reserved parameter was used.
  static DeviceGeometry build(const DeviceParameters& p);

  const DeviceParameters& parameters() const { return p_; }
  const std::vector<BoundaryCurve>& boundaries() const { return boundaries_; }
  const std::vector<RegionBody>& regions() const { return regions_; }
  const std::vector<NamedFeature>& features() const { return features_; }

  const RegionBody& region(Region r) const;
  Vec2 feature(FeatureId f) const;
  /// All boundary pieces carrying the given identifier.
  std::vector<const BoundaryCurve*> boundaries_with(BoundaryId id) const;

  // --- derived quantities ---------------------------------------------------
  Real contact_radius() const { return 0.5 * p_.phi_2; }
  Real land_width() const { return 0.5 * (p_.phi_1 - p_.phi_2); }
  Real cone_half_angle() const;      ///< of the outer taper [rad]
  Real extractor_outer_radius() const;
  Real domain_revolved_volume() const;

  void print(std::FILE* out) const;
  /// Write the geometry as CSV into `dir` for plotting.  Nothing is read back.
  void write_csv(const std::string& dir) const;

 private:
  DeviceParameters p_;
  std::vector<BoundaryCurve> boundaries_;
  std::vector<RegionBody> regions_;
  std::vector<NamedFeature> features_;
};

// ---------------------------------------------------------------------------
// Measures of revolution -- exposed because they are what the analytic tests
// pin down.
// ---------------------------------------------------------------------------

/// Surface of revolution of an open meridian polyline: 2*pi * int r ds.
/// Exact for straight segments.
Real revolved_area(const std::vector<Vec2>& polyline);

/// Volume of revolution enclosed by a closed meridian loop, from
/// int_A 2*pi*r dA = oint pi r^2 dz (Green).  Exact for polygonal loops.
/// The sign follows the orientation; the magnitude is what callers want.
Real revolved_volume(const std::vector<Vec2>& closed_loop);

/// Signed area of a closed meridian polygon; positive for counter-clockwise.
Real meridian_signed_area(const std::vector<Vec2>& closed_loop);

/// Total meridian length of a polyline.
Real meridian_length(const std::vector<Vec2>& polyline);

}  // namespace es
