#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/axisym_fem.hpp"
#include "es/materials.hpp"
#include "es/types.hpp"
#include "es/volume_mesh.hpp"

namespace es {

// ===========================================================================
// P2b -- the dielectric electrostatic model of the capillary Kunze emitter
// ===========================================================================
//
// THE CORRECTED PHYSICAL CONTRACT.  What is on potential and what is not:
//
//   ionic liquid           IDEAL CONDUCTOR at V_emitter.  It fills the bore and
//                          is cut at the feed boundary; the whole liquid volume
//                          is one equipotential, so its entire surface -- bore
//                          wall, flat reference plane at z = 0, and the feed
//                          cross section -- carries the Dirichlet condition.
//   emitter body (SU-8)    DIELECTRIC.  eps_r from materials.hpp.  It is NOT an
//                          electrode: it carries no free charge, it only
//                          polarises, and no surface of it is fixed to any
//                          potential.  In P2a it was metal.  That was wrong.
//   extractor carrier      DIELECTRIC, same treatment.
//   metallised face        IDEAL CONDUCTOR at V_extractor, ZERO THICKNESS.
//                          Which faces are coated is a parameter, because it is
//                          a manufacturing fact and not a law.
//   vacuum                 eps_r = 1.
//
//   liquid feed boundary   phi = V_emitter on the LIQUID CROSS SECTION ONLY,
//                          r <= r_bore at z = liquid_feed_z.  The remainder of
//                          that plane is the rear face of the polymer and is an
//                          ordinary dielectric/vacuum interface.  Treating the
//                          whole cut plane as an electrode would put a grounded
//                          -- or worse, energised -- disc behind the emitter
//                          that does not exist, and would dominate the field it
//                          is supposed to leave alone.
//
// Solved:  div(eps grad phi) = 0, axisymmetric, static, no free charge.
// Not solved, and not pretended: space charge, emission, meniscus motion,
// flow, finite liquid conductivity, time dependence.

// ---------------------------------------------------------------------------
// Tolerances for the feed-boundary study -- fixed BEFORE the measurement
// ---------------------------------------------------------------------------
//
// The requirement on the feed boundary is that pushing it further back stops
// changing the field at the meniscus.  A quantity counts as independent of
// liquid_feed_z when DOUBLING the modelled column length changes it by less
// than the bound below.  The bounds are the same ones P2a fixed for its own
// truncation study, and for the same reasons: the mesh discretisation error at
// the reference level is a few times 1e-4, so a bound at 1e-3 keeps the two
// error sources apart; and a field-dependent emission law is worth evaluating
// at about a per cent, so 1e-3 is an order below what the answer is used for.
//
// THE STUDY FAILS THESE BOUNDS, AND THAT IS A RESULT ABOUT THE MODEL.
//
// Measured on the reference geometry, doubling the column length from 400 um to
// 800 um moves E_z on the axis above the reference plane by about 3.6 per cent
// and the potential at mid gap by about 1.5 per cent of the applied span; the
// emitter charge nearly doubles.  The changes do shrink -- roughly like
// 1/ln(L/a) -- but nowhere near fast enough, and centimetres of column would be
// needed to reach one part in a thousand.
//
// THE REASON, and why it is not a numerical defect.  The liquid column is an
// energised conductor.  A thin cylinder of radius a and length L has a
// self-capacitance of about 2 pi eps0 L / (ln(2L/a) - 1); making it longer makes
// it hold proportionally more charge, and that charge is felt at the tip.  The
// measured Q_emitter follows that formula to within the contribution of the tip
// and the electrode.  Nothing about the far field is responsible: repeating the
// study inside a GROUNDED enclosure at 25 mm changes every number by less than
// 0.1 per cent, so it is not an artefact of the open boundary either.
//
// WHAT WOULD FIX IT, and why P2b does not do it.  The real emitter does not end
// in mid vacuum: docs/04_geometry_model.md, section 4.1, already provides for a
// base plate at emitter potential from which the tapered structure protrudes.
// A real conductor of stated dimensions terminating the column would make the
// feed position a numerical parameter again.  Adding one is a DEVICE GEOMETRY
// decision, and this phase was explicitly instructed not to treat the rear cut
// plane as an electrode -- which is right, because a disc invented to make a
// convergence study succeed is not a device.
//
// CONSEQUENCE FOR THE RESULTS.  Every P2b number is reported as a function of
// liquid_feed_z, and the value used is an EXAMPLE VALUE, not a measured
// dimension -- exactly like extractor_outer_radius.
namespace feed_truncation {
inline constexpr Real kTolPhiOverSpan = 1.0e-3;
inline constexpr Real kTolFieldRelative = 1.0e-3;
}  // namespace feed_truncation

/// Which bodies are conductors.
enum class ConductorModel {
  /// The P2b physics above.
  Dielectric = 0,
  /// SUPERSEDED P2a ARRANGEMENT, kept for ONE purpose: it is the configuration
  /// the independent boundary-element solver can also solve, so it is how the
  /// finite element machinery is cross-checked against the BEM at eps_r = 1.
  /// Emitter body and liquid are one perfect conductor, the whole extractor is
  /// another, all dielectrics are forced to eps_r = 1.  It is a solver test.
  /// It is NOT a model of this device and no physical claim may rest on it.
  MetallicReference,
};
const char* to_string(ConductorModel c);

/// Which extractor surfaces carry the metal film.
enum class Metallisation {
  FrontOnly = 0,     ///< emitter-facing face only
  FrontAndAperture,  ///< front face and the aperture wall -- the P2b reference
  AllSurfaces,       ///< fully coated electrode; also what MetallicReference uses
};
const char* to_string(Metallisation m);

/// What a node is, electrically.  Everything downstream -- the Dirichlet set,
/// the charge sums, the audit -- is expressed through this and never through a
/// coordinate test.
enum class NodeRole {
  Free = 0,
  LiquidConductor,        ///< liquid surface or interior; phi = V_emitter
  LiquidFeedBoundary,     ///< the liquid cross section at z = liquid_feed_z
  ExtractorMetallisation, ///< the metal film; phi = V_extractor
  FarFieldDirichlet,      ///< only with FarField::Grounded
  EmitterMetalReference,  ///< only in ConductorModel::MetallicReference
};
const char* to_string(NodeRole r);

// ---------------------------------------------------------------------------

struct DielectricSetup {
  DielectricDeviceParameters geometry;
  DielectricMaterials materials;
  ConductorModel conductor_model{ConductorModel::Dielectric};
  Metallisation metallisation{Metallisation::FrontAndAperture};
  FarField far_field{FarField::Asymptotic};
  Real V_emitter{1500.0};
  Real V_extractor{0.0};

  Real applied_span() const { return V_emitter - V_extractor; }
};

// ---------------------------------------------------------------------------
// The audit that no polymer surface became a conductor
// ---------------------------------------------------------------------------
//
// This is a required check, not a nicety: the single defect P2b exists to
// repair is a dielectric that was treated as an electrode.  The audit walks the
// NAMED polymer surfaces of the device -- outer flank, tip land, rear face,
// body interior, and the uncoated extractor faces -- and reports every node on
// them that carries a Dirichlet condition.  The count must be zero.
struct BoundaryAudit {
  Index n_nodes{0};
  Index n_dirichlet{0};
  Index n_liquid{0}, n_feed{0}, n_metal{0}, n_far{0}, n_emitter_metal_reference{0};
  /// Nodes on a polymer surface that are nevertheless fixed.  MUST be zero in
  /// ConductorModel::Dielectric.
  Index n_polymer_dirichlet{0};
  /// Nodes at z = liquid_feed_z that are fixed but lie OUTSIDE the liquid cross
  /// section.  MUST be zero: the cut plane is not an electrode.
  Index n_feed_plane_outside_liquid{0};
  /// Nodes where a coated face meets an uncoated one.  They belong to the metal
  /// film's rim, are legitimately fixed, and are counted rather than silently
  /// skipped so that the exclusion is visible.
  Index n_film_rim_shared{0};
  std::vector<std::string> violations;

  bool ok() const { return violations.empty(); }
  void print(std::FILE* out) const;
  void write_csv(const std::string& path) const;
};

// ---------------------------------------------------------------------------

struct Probe {
  std::string name;
  Vec2 x;
  Real clearance{0.0};  ///< distance to the nearest unrounded edge [m]
  std::string note;
};

/// Fixed evaluation points, placed from the geometry so that they are the same
/// physical points at every mesh level, every feed position and every eps_r.
/// All of them lie in vacuum and away from every unrounded edge; `clearance`
/// records how far away, because "edge-far" is a measurement, not an adjective.
std::vector<Probe> dielectric_probes(const DeviceVolumeMesh& m);

// ---------------------------------------------------------------------------

struct DielectricSolution {
  DeviceVolumeMesh mesh;
  std::vector<NodeRole> role;
  std::vector<char> emitter_mask, extractor_mask;
  AxisymSolution fem;
  BoundaryAudit audit;
  /// Per-cell permittivity and activity, kept so that a caller can recover
  /// fields (field_recovered_at) without rebuilding the problem.
  std::vector<Real> cell_eps_r;
  std::vector<char> cell_active;

  Real Q_emitter{0.0}, Q_extractor{0.0}, Q_net{0.0};   ///< [C]

  std::vector<Probe> probes;
  /// All three from the RECOVERED field, so that a mesh study measures the
  /// discretisation error and not the cell-to-cell jitter of a raw Q1 gradient.
  std::vector<Real> phi_probe;    ///< [V]
  std::vector<Real> Ez_probe;     ///< [V/m]
  std::vector<Real> Emag_probe;   ///< [V/m]

  /// One-sided normal field just above the flat liquid reference plane, per
  /// radial cell, r < r_bore.  It is the quantity a meniscus model will need,
  /// and it is edge-limited: the last cell before the pinned edge is inside the
  /// singularity and is flagged.
  std::vector<Real> surface_r, surface_Ez;
  Index surface_edge_cells{0};    ///< cells excluded as edge-affected

  /// Normal electric flux density either side of the polymer/vacuum interface
  /// on the taper flank; the two must agree.  Sampled at mid height.
  Real Dn_polymer_side{0.0}, Dn_vacuum_side{0.0};
  Real phi_interface_jump{0.0};   ///< identically zero by construction; measured anyway

  Real relative_interface_error() const;
  void print(std::FILE* out) const;
  void write_csv(const std::string& dir) const;
};

/// Assemble and solve the device problem.  Throws on an inconsistent setup --
/// most importantly, MetallicReference with anything other than eps_r = 1.
DielectricSolution solve_dielectric(const DielectricSetup& s);

/// Node roles alone, without solving.  Used by the audit test.
std::vector<NodeRole> node_roles(const DeviceVolumeMesh& m, const DielectricSetup& s);
BoundaryAudit audit_boundaries(const DeviceVolumeMesh& m, const std::vector<NodeRole>& role,
                               const DielectricSetup& s);

}  // namespace es
