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
  /// The dielectric body of the liquid reservoir behind the base plate.  A
  /// SEPARATE region from EmitterSolid on purpose: it is a different part, made
  /// of a different polymer (PEEK in the built device), and keeping it apart is
  /// what lets the boundary audit name it and a later phase give it its own
  /// permittivity.  It is a DIELECTRIC.  It is never an electrode.
  ReservoirSolid,
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
  /// NOT a part of the device.  The disc that closes the emitter conductor at
  /// the rear end of its NUMERICAL continuation, so that the boundary-integral
  /// formulation sees a closed perfect conductor instead of an open sheet.
  /// It is never a computational-domain edge and never a real component.
  NumericalEmitterBackClosure,
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
  /// Rim where the numerical rearward continuation meets its closing disc.
  /// An artefact of the closure, not a device edge: it exists only so that the
  /// conductor is closed, it lies far behind the physically evaluated region,
  /// and no field value from it is ever reported.
  NumericalBackClosureEdge,
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

  /// Radial extent of the modelled electrode [m].  MANDATORY.
  ///
  /// It used to be optional, with 0 meaning "out to the domain boundary".  That
  /// silently identified a conductor with the open edge of the computational
  /// box, which is not a conductor.  A boundary-integral vacuum solve has no
  /// truncation boundary at all (the free-space kernel decays on its own), so
  /// the electrode needs a real outer radius or it has none.  Zero is now
  /// rejected, and domain_radius must be STRICTLY larger.
  ///
  /// The value below is an EXAMPLE VALUE, not a measured dimension.
  Real extractor_outer_radius{2.0e-3};

  /// Axial length of the modelled emitter conductor behind the tip plane z = 0
  /// [m].  A DIMENSION OF THE DEVICE, and mandatory for every P2a result.
  ///
  ///   0  -- off.  The emitter solid and the liquid column then run down to
  ///         domain_z_min and are cut there, exactly as in P1: the conductor
  ///         ends as an OPEN sheet.  That is a legitimate geometric sketch, but
  ///         it is not a boundary-integral conductor, and vacuum_bem_mesh()
  ///         refuses it.
  ///   > 0 -- the emitter solid and the liquid column end at z = -this value and
  ///         are closed there by a full conducting disc from r = 0 out to the
  ///         foot radius, tagged NumericalEmitterBackClosure.
  ///
  /// WHY IT IS A GEOMETRY PARAMETER AND NOT A CONVERGENCE KNOB.  The intention
  /// was the opposite: extend the conductor backwards far enough that the local
  /// tip field stops moving, quote the local field, and admit only the total
  /// capacitance as length dependent.  The measurement says that is not
  /// available in this model.  Doubling the length from 200 um to 3.2 mm moves
  /// E_z at the axial reference point by 5.6, 4.0, 2.5 and 1.5 per cent per
  /// doubling -- a 1/L tail that would need centimetres of emitter to reach one
  /// part in a thousand -- and moves c_EX by about a third per doubling with no
  /// sign of settling at all.
  ///
  /// The reason is in the boundary condition, not in the discretisation.  With
  /// V -> 0 at infinity and only two conductors, the system carries net charge:
  /// there is no return electrode and no enclosure anywhere in the model, so a
  /// longer emitter simply holds more charge, and that charge is felt at the tip
  /// as well as at the extractor.  Truncating the conductor is therefore not a
  /// numerical detail that can be made to disappear; it fixes a dimension of the
  /// device.  Whether the local field becomes insensitive once a grounded
  /// enclosure is present is a question for the phase that adds one.
  ///
  /// Consequence for the results: NOTHING in P2a is truncation-converged.
  /// E_z(ref), c_EX, C_m and c_EE are all reported as functions of this
  /// parameter, and the value used is an EXAMPLE VALUE, not a measured
  /// dimension -- exactly like extractor_outer_radius above.
  ///
  /// The closing disc itself remains numerical and is tagged as such: a flat cap
  /// is a stand-in for whatever really terminates the emitter, its rim is marked
  /// as a non-evaluable edge, and no field is ever read off it.
  Real emitter_back_length{0.0};

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

  /// Is the emitter conductor closed by the numerical rearward continuation?
  bool has_back_closure() const { return p_.emitter_back_length > 0.0; }
  /// z of the closing disc [m].  Throws if there is no closure.
  Real back_closure_z() const;
  /// Rear limit of the region in which a field value is a physical result.
  /// The taper foot: everything behind it is a shank of arbitrary, numerically
  /// chosen length and carries no device meaning.
  Real evaluation_z_min() const { return -p_.emitter_height; }
  /// Axial distance between the closing disc and that region [m].
  Real back_closure_clearance() const;

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
