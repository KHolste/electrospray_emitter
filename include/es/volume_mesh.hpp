#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/axisym_fem.hpp"
#include "es/device_geometry.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P2b -- automatic axisymmetric VOLUME mesh of the capillary emitter
// ===========================================================================
//
// docs/07_mesher_decision.md, section 7.3, is binding here and says two things:
// no general unstructured volume mesher in-house, and no heavyweight external
// dependency introduced casually.  What is built here is neither.
//
// It is a BLOCK-STRUCTURED, LOGICALLY RECTANGULAR mesh that follows from the
// device parameters in closed form.  The meridian half-plane of this device is
// a stack of z-strips, each of which is a stack of radial blocks, and the only
// non-axis-aligned boundary in the whole geometry is the straight taper flank.
// A tensor-product grid in (r_ref, z), radially WARPED so that one grid line
// lies exactly on the flank, therefore represents every material interface
// exactly, with no hanging nodes, no staircasing and no Delaunay machinery:
//
//     r(i, j) = W( r_ref[i], z[j] )
//
// where W is a continuous, piecewise linear, strictly increasing map of the
// radius that pins r = r_bore and r = r_anchor and sends the reference radius
// r_land to the local outer radius of the emitter, rho(z).  Because rho(z) is
// piecewise linear in z and the grid rows are the breakpoints, the warped grid
// lines are straight inside every cell -- so the mesh reproduces every region
// volume EXACTLY, which is check 5 of validate().
//
// WHAT THIS DOES NOT DO.  It does not adapt, it does not refine a posteriori,
// and it does not handle a moving free surface.  P3 needs all three; when that
// comes, the decision in 7.3 has to be made properly.  This mesher is the
// static-geometry answer, it is a pure function of the parameters, and it is
// cheap to throw away.
//
// WHAT DOES NOT COME FROM THE USER.  No element count, no element size, no
// refinement zone.  The single control is an integer `mesh_level`, which scales
// the whole size field by 2^(-level/2) and exists so that convergence can be
// measured.  It is not a quality knob.

/// Which replacement the model uses for everything behind the emitter base.
///
/// THIS ENUM IS THE POINT OF P2c.  Until now a single number, `liquid_feed_z`,
/// moved three different things at once -- the length of the CONDUCTING liquid
/// column, the rearward extent of the DIELECTRIC body, and with it the whole
/// rear device geometry.  A study that varied it was therefore not moving a
/// boundary; it was building a different device each time.  The two entries
/// below separate the two questions that were tangled up in it.
enum class ReservoirModel {
  /// SUPERSEDED, kept as a DIAGNOSIS.  The dielectric body and the liquid
  /// column both stop at the rear face of the base plate and the liquid is cut
  /// there; the emitter potential is imposed on the liquid cross section only.
  /// Growing `base_plate_thickness` in this mode lengthens the conductor and
  /// the dielectric together -- which is exactly the coupling that made the old
  /// "convergence against the position of the feed boundary" misleading.  No
  /// result of this mode is a statement about a reservoir.
  TruncatedColumn = 0,
  /// The replacement model.  The device geometry in front of the base plate is
  /// FIXED; the bore continues through a fixed feed channel into a coaxial
  /// liquid plenum inside a dielectric (PEEK) body.  All connected liquid --
  /// bore, channel and plenum -- is one equipotential at V_emitter.  There is
  /// no metal holder, no rear disc and no conducting back plate anywhere.
  AxisymmetricPlenum,
};
const char* to_string(ReservoirModel m);

/// Parameters of the dielectric device model.
///
/// EVERY LENGTH BELOW IS ONE THING.  The grouping is the deliverable: the front
/// emitter, the base body, the feed channel and the reservoir are four separate
/// parameter groups, and changing one of them may not move any of the others.
struct DielectricDeviceParameters {
  /// P1 geometry of the FRONT of the device -- taper, bore, extractor, domain.
  /// emitter_back_length MUST be zero here: the conducting rear closure is a
  /// P2a construction for a metallic emitter and has no place in a dielectric
  /// model.
  DeviceParameters device;

  // --- the dielectric base body --------------------------------------------

  /// Axial thickness of the dielectric base body behind the taper foot [m], so
  /// that the printed emitter ends at z = -(emitter_height + this).  Its radius
  /// is the foot radius phi_3/2 -- the cylindrical continuation of the printed
  /// structure, which is a dimension the drawing gives.
  ///
  /// PROVISIONAL EXAMPLE VALUE, not a measured dimension.  The Kunze
  /// dissertation dimensions the tapered structure, not what carries it, and
  /// docs/04_geometry_model.md read a wide "base plate at emitter potential"
  /// out of a hatched bar in a sketch -- a reading that is neither dimensioned
  /// nor, as it turns out, electrically right (see docs/08, 8.9).  What is
  /// modelled here is a stated thickness of DIELECTRIC, flagged as provisional
  /// wherever it is reported.
  ///
  /// It is FIXED across every reservoir comparison.  In AxisymmetricPlenum it
  /// is a dimension of the fixed front geometry and nothing else.
  Real base_plate_thickness{1.4e-4};

  // --- what stands behind it ------------------------------------------------

  ReservoirModel reservoir{ReservoirModel::TruncatedColumn};

  /// Radius of the liquid feed channel behind the base plate [m].  0 means "the
  /// bore radius", which is the reference case.  It may not EXCEED the bore
  /// radius: r > r_bore is inside the radially warped zone of the mesher and
  /// would no longer fall on a grid line exactly.
  Real feed_channel_radius{0.0};

  /// Axial length of the fixed feed channel [m] -- equivalently, the thickness
  /// of the reservoir body's top wall, which the channel passes through.  It is
  /// the stand-off between the emitter base and the plenum roof.  FIXED across
  /// every reservoir comparison.  PROVISIONAL EXAMPLE VALUE.
  Real feed_channel_length{3.0e-4};

  /// Radius of the coaxial liquid plenum [m].
  ///
  /// It must lie radially OUTSIDE the modelled extractor
  /// (> device.extractor_outer_radius).  That is not a convenience: the radial
  /// node list is a tensor factor, so a level inside the near field would move
  /// near-field nodes when the plenum is resized and the comparison would be
  /// measuring the mesh as well.  Outside the extractor, only nodes beyond the
  /// extractor rim move, and every node of the front device stays bitwise
  /// identical -- which tests/test_reservoir.cpp checks rather than assumes.
  ///
  /// ERSATZGEOMETRIE.  Kunze's Abb. A.5 shows an improved reservoir that is NOT
  /// a solid of revolution, roughly 25 mm high with about 10 mm outer and 8 mm
  /// inner width.  Those figures are used HERE ONLY to pick a sensible order of
  /// magnitude -- millimetres, not micrometres.  The axisymmetric plenum is a
  /// replacement geometry and is never claimed to be a reconstruction of it.
  Real plenum_radius{2.5e-3};
  /// Axial depth of the liquid plenum cavity [m].
  Real plenum_depth{1.0e-3};
  /// Thickness of the dielectric wall and floor around the cavity [m].
  Real plenum_wall_thickness{5.0e-4};
  /// Filled fraction of the cavity, measured DOWN from the roof, in (0, 1].
  /// 1 means completely filled, which is the reference case.  Gravity, wetting
  /// and the shape of a free surface inside the reservoir are NOT modelled; the
  /// fill level is a flat plane and a parameter of the replacement model.
  Real plenum_fill_fraction{1.0};

  /// Size-field scale exponent.  h -> h * 2^(-level/2).  Level 0 is the base
  /// mesh; each level is a factor sqrt(2) in element size in both directions.
  int mesh_level{2};

  /// Rear face of the dielectric base body [m].  In TruncatedColumn this is
  /// also where the liquid column is cut.
  Real base_z() const { return -(device.emitter_height + base_plate_thickness); }
};

/// Fixed, documented constants of the size field.  Same shape as the boundary
/// mesher in docs/07_mesher_decision.md, section 7.2, and for the same reason:
/// every number that controls the mesh is here, named, and none of them is a
/// user input.
namespace volume_mesher {
inline constexpr Real kGradation = 0.25;        ///< Lipschitz bound on h(x)
inline constexpr int kFeatureDivisions = 16;    ///< h at a named feature = lfs/16
inline constexpr int kKeyDivisions = 8;         ///< h at any other key level = lfs/8
inline constexpr int kDomainDivisions = 40;     ///< h_max = domain diagonal / 40
inline constexpr Index kMinCellsPerInterval = 3;
inline constexpr Real kMaxNeighbourRatio = 1.6; ///< checked, not enforced by smoothing
}  // namespace volume_mesher

// ---------------------------------------------------------------------------

struct VolumeMeshCheck {
  std::string name;
  bool passed{false};
  Real measured{0.0}, bound{0.0};
  std::string note;
};

struct VolumeMeshReport {
  std::vector<VolumeMeshCheck> checks;
  bool all_passed() const;
  void print(std::FILE* out) const;
  void write_csv(const std::string& path) const;
};

// ---------------------------------------------------------------------------

struct DeviceVolumeMesh {
  DielectricDeviceParameters p;
  QuadMesh grid;
  std::vector<Region> cell_region;  ///< one per cell, size grid.n_cells()

  /// Unwarped radial node coordinates and the axial ones.  The key levels are
  /// exactly members of these lists; that is checked, not assumed.
  std::vector<Real> r_ref, z_ref;

  // Landmark indices.  Every boundary condition is expressed through these and
  // never through a coordinate comparison, so no surface can be misidentified
  // by a rounding accident.
  Index i_axis{0};        ///< r = 0
  Index i_channel{0};     ///< r = r_channel (feed-channel wall; = i_bore by default)
  Index i_bore{0};        ///< r = r_bore   (bore wall / liquid edge)
  Index i_land{0};        ///< reference radius of the tip land; warps to rho(z)
  Index i_foot{0};        ///< r = r_foot   (reference position only)
  Index i_aperture{0};    ///< r = aperture radius
  Index i_ext_outer{0};   ///< r = extractor outer radius
  Index i_plenum{-1};     ///< r = plenum_radius        (-1 without a plenum)
  Index i_plenum_outer{-1};  ///< r = plenum outer radius (-1 without a plenum)
  Index i_far{0};         ///< r = domain radius
  Index j_min{0};         ///< z = domain z_min
  Index j_block_bottom{-1};  ///< z = underside of the reservoir body (-1 without)
  Index j_cav_bottom{-1};    ///< z = floor of the plenum cavity      (-1 without)
  Index j_fill{-1};          ///< z = liquid level in the plenum      (-1 without)
  Index j_roof{-1};          ///< z = roof of the plenum cavity       (-1 without)
  Index j_base{0};        ///< z = rear face of the dielectric base body
  Index j_foot{0};        ///< z = -emitter_height
  Index j_tip{0};         ///< z = 0, the flat free-surface reference plane
  Index j_ex_front{0};    ///< z = extraction_distance
  Index j_ex_back{0};     ///< z = extraction_distance + thickness
  Index j_max{0};         ///< z = domain z_max

  Real r_bore{0.0}, r_land{0.0}, r_foot{0.0}, r_aperture{0.0}, r_ext_outer{0.0};
  Real r_channel{0.0};    ///< feed-channel radius; <= r_bore
  Real r_plenum{0.0}, r_plenum_outer{0.0};
  Real r_anchor{0.0};     ///< outer fixed point of the radial warp
  /// Axial landmarks, in decreasing z: base body rear face, plenum roof, liquid
  /// level, cavity floor, underside of the reservoir body.  The last four are
  /// only meaningful with ReservoirModel::AxisymmetricPlenum.
  Real z_base{0.0}, z_roof{0.0}, z_fill{0.0}, z_cav_bottom{0.0}, z_block_bottom{0.0};

  bool has_plenum() const { return p.reservoir == ReservoirModel::AxisymmetricPlenum; }
  /// Rear end of the modelled device: the underside of the reservoir body with
  /// a plenum, the rear face of the base body without one.
  Real device_z_min() const { return has_plenum() ? z_block_bottom : z_base; }

  /// Outer radius of the emitter body at height z [m]; r_foot below the taper,
  /// r_land at the tip plane, and 0 where there is no body.
  Real emitter_outer_radius_at(Real z) const;
  /// The radial warp: reference radius -> physical radius at height z.
  Real warp(Real r_reference, Real z) const;

  /// Which material occupies the point, from the closed-form geometry alone --
  /// no indices, no mesh.  validate() uses it to check the index-based region
  /// assignment, so the two can never quietly agree by construction.
  Region region_at(Vec2 x) const;

  Index n_cells_of(Region r) const;
  Real revolved_volume_of(Region r) const;
  /// Closed-form volume of a region, independent of the mesh.
  Real analytic_volume_of(Region r) const;

  VolumeMeshReport validate() const;
  void print(std::FILE* out) const;
  /// nodes.csv, cells.csv, mesh_checks.csv into `dir`.
  void write_csv(const std::string& dir) const;
};

/// Build the volume mesh.  Pure function of the parameters: same input, bitwise
/// same output.  Throws std::runtime_error on a parameter set this mesher
/// cannot represent, and NotImplementedInThisPhase if a P2a-only parameter
/// (the conducting rear closure) is still set.
DeviceVolumeMesh build_volume_mesh(const DielectricDeviceParameters& p);

/// Size-field scale for a level: 2^(-level/2).
Real mesh_level_scale(int level);

}  // namespace es
