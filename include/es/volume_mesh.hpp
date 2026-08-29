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

/// Parameters of the P2b dielectric device model.
struct DielectricDeviceParameters {
  /// P1 geometry.  emitter_back_length MUST be zero here: the conducting rear
  /// closure is a P2a construction for a metallic emitter and has no place in a
  /// dielectric model.  The liquid column is terminated by the feed boundary
  /// below instead.
  DeviceParameters device;

  /// z of the liquid feed boundary [m].  The ionic liquid is cut here and the
  /// EMITTER POTENTIAL IS IMPOSED ON THE LIQUID CROSS SECTION ONLY, r <= r_bore.
  /// The rest of that plane is the polymer's rear face and is NOT an electrode.
  ///
  /// The reservoir is not meshed.  This boundary stands in for it, and its
  /// position is a modelling choice whose influence is measured, not assumed:
  /// moving it further back must stop changing the field at the meniscus.
  ///
  /// It is also, by construction, the place where the HYDRAULIC inlet condition
  /// will go in P3 -- reservoir pressure, volume flow, or a hydraulic impedance.
  /// Nothing hydraulic is implemented here; the boundary is named and located so
  /// that adding one is a new condition on an existing entity, not a new entity.
  Real liquid_feed_z{-2.0e-4};

  /// Size-field scale exponent.  h -> h * 2^(-level/2).  Level 0 is the base
  /// mesh; each level is a factor sqrt(2) in element size in both directions.
  int mesh_level{2};
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
  Index i_bore{0};        ///< r = r_bore   (bore wall / liquid edge)
  Index i_land{0};        ///< reference radius of the tip land; warps to rho(z)
  Index i_foot{0};        ///< r = r_foot   (reference position only)
  Index i_aperture{0};    ///< r = aperture radius
  Index i_ext_outer{0};   ///< r = extractor outer radius
  Index i_far{0};         ///< r = domain radius
  Index j_min{0};         ///< z = domain z_min
  Index j_feed{0};        ///< z = liquid_feed_z
  Index j_foot{0};        ///< z = -emitter_height
  Index j_tip{0};         ///< z = 0, the flat free-surface reference plane
  Index j_ex_front{0};    ///< z = extraction_distance
  Index j_ex_back{0};     ///< z = extraction_distance + thickness
  Index j_max{0};         ///< z = domain z_max

  Real r_bore{0.0}, r_land{0.0}, r_foot{0.0}, r_aperture{0.0}, r_ext_outer{0.0};
  Real r_anchor{0.0};     ///< outer fixed point of the radial warp

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
