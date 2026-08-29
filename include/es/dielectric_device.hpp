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
//   ionic liquid           IDEAL CONDUCTOR at V_emitter.  ALL connected liquid
//                          is ONE equipotential -- bore, feed channel and, with
//                          a plenum, the reservoir -- so its entire surface
//                          carries the Dirichlet condition and nothing else
//                          does.
//   emitter body (SU-8)    DIELECTRIC.  eps_r from materials.hpp.  It is NOT an
//                          electrode: it carries no free charge, it only
//                          polarises, and no surface of it is fixed to any
//                          potential.  In P2a it was metal.  That was wrong.
//   reservoir body (PEEK)  DIELECTRIC, same treatment.  There is NO conducting
//                          holder, NO rear metal disc and NO base plate on
//                          emitter potential anywhere in this model.
//   extractor carrier      DIELECTRIC, same treatment.
//   metallised face        IDEAL CONDUCTOR at V_extractor, ZERO THICKNESS.
//                          Which faces are coated is a parameter, because it is
//                          a manufacturing fact and not a law.
//   vacuum                 eps_r = 1.
//
//   rear cut plane         ONLY in ReservoirModel::TruncatedColumn, which is
//                          kept as a diagnosis of the superseded arrangement:
//                          phi = V_emitter on the LIQUID CROSS SECTION ONLY.
//                          The remainder of that plane is the rear face of the
//                          polymer and is an ordinary dielectric interface.
//                          Treating the whole cut plane as an electrode would
//                          put a disc behind the emitter that does not exist
//                          and would dominate the field it is supposed to leave
//                          alone.
//
// Solved:  div(eps grad phi) = 0, axisymmetric, static, no free charge.
// Not solved, and not pretended: space charge, emission, meniscus motion,
// flow, finite liquid conductivity, time dependence.

// ---------------------------------------------------------------------------
// Tolerances for the reservoir study -- fixed BEFORE the measurement
// ---------------------------------------------------------------------------
//
// The requirement is that making the modelled liquid reservoir LARGER stops
// changing the field at the meniscus.  A quantity counts as independent of the
// reservoir when growing it by one step changes the quantity by less than the
// bound below.  The bounds are the ones P2a and P2b fixed for their own
// truncation studies, and for the same reasons: the mesh discretisation error
// at the reference level is a few times 1e-4, so a bound at 1e-3 keeps the two
// error sources apart; and a field-dependent emission law is worth evaluating
// at about a per cent, so 1e-3 is an order below what the answer is used for.
//
// THEY ARE NOT LOOSENED AFTERWARDS.  Whatever the measurement says, it is
// reported against these numbers.
//
// WHAT THE OLD "FEED BOUNDARY" STUDY ACTUALLY MEASURED -- the correction this
// phase exists for.  Until now one parameter, liquid_feed_z, set the rear end
// of the CONDUCTING liquid column, the rear end of the DIELECTRIC emitter body,
// and with them the whole rear geometry of the device, all at once.  Varying it
// was therefore never "moving a boundary condition": it built a longer
// energised conductor inside a longer dielectric each time, and the strong
// field change that followed is what a longer conductor does.  Reporting that
// as a failed convergence "against the position of the feed boundary" was
// misleading, and section 8.9 of docs/08_dielectric_model.md says so.
//
// The position of a fully immersed electrical contact, by contrast, is
// irrelevant in a model whose liquid is one ideal conductor, and it is not
// represented geometrically at all.
//
// LOCAL AND GLOBAL QUANTITIES ARE JUDGED SEPARATELY.  The total charge on a
// conductor of finite size in an open domain is a function of that size --
// there is no return electrode, so a larger vessel simply holds more charge,
// and no tolerance is placed on it.  What has to converge is the LOCAL
// extraction field: phi at edge-far probe points and E just above the flat
// liquid surface.
namespace reservoir_convergence {
inline constexpr Real kTolPhiOverSpan = 1.0e-3;
inline constexpr Real kTolFieldRelative = 1.0e-3;
}  // namespace reservoir_convergence

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
  /// ONLY in ReservoirModel::TruncatedColumn: the liquid cross section at the
  /// rear cut, where the column is chopped off.  With a plenum the liquid is
  /// closed and there is no cut, so this role does not occur.
  LiquidFeedBoundary,
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
// repair is a dielectric that was treated as an electrode, and the reservoir
// added afterwards is exactly the place where a conducting holder would be
// easiest to slip in by accident.
//
// The audit therefore runs TWO checks, not one.
//
//   STRUCTURAL, and the one that decides.  A node may carry the emitter
//   potential only if it touches a LIQUID cell.  Every other fixed node -- on
//   the emitter body, on the reservoir body, on an uncoated extractor face, or
//   anywhere else -- is a violation.  This does not depend on anybody having
//   named the surface, so a new part cannot be added with a conducting face
//   that nobody thought to check.
//
//   NAMED, and the one that explains.  The same test on an explicit list of
//   surfaces -- taper flank, tip land, rear face, emitter interior, reservoir
//   top face, outer rim, underside, reservoir interior, uncoated extractor
//   faces -- so that a violation can be reported by the surface it sits on.
//
// Both counts must be zero.
struct BoundaryAudit {
  Index n_nodes{0};
  Index n_dirichlet{0};
  Index n_liquid{0}, n_feed{0}, n_metal{0}, n_far{0}, n_emitter_metal_reference{0};
  /// Fixed nodes that touch NO liquid cell and are neither the extractor film
  /// nor the far-field box.  This is the STRUCTURAL form of the check: it does
  /// not depend on a list of surface names, so a surface nobody thought to name
  /// cannot slip through.  MUST be zero in ConductorModel::Dielectric.
  Index n_polymer_dirichlet{0};
  /// Same count restricted to the NAMED polymer surfaces below, kept so that a
  /// violation can be reported by the surface it sits on.
  Index n_named_surface_dirichlet{0};
  /// Nodes on the rear cut plane of a truncated column that are fixed but lie
  /// OUTSIDE the liquid cross section.  MUST be zero: the cut plane is not an
  /// electrode.  Always zero with a plenum, which has no cut plane.
  Index n_feed_plane_outside_liquid{0};
  /// Fixed nodes on the surface of the reservoir body that are liquid surface
  /// (channel wall, plenum roof, plenum wall, plenum floor).  Reported so that
  /// the liquid really is one connected equipotential and not silently split.
  Index n_reservoir_liquid_surface{0};
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

/// Which post-processing a solve performs.
enum class DielectricDiagnostics {
  /// Everything: the probe points, the one-sided normal field on the FLAT
  /// liquid reference plane and the interface flux check on the taper flank.
  /// This is what P2b and P2c report.
  Full = 0,
  /// Field, charges and audit only.  For a mesh whose free surface has been
  /// DEFORMED (P3b): the three diagnostics above assume the P2c mesh -- level
  /// rows, so that point location is exact, and a flat liquid surface at z = 0
  /// -- and on a deformed mesh they would either locate a point in the wrong
  /// cell or describe a plane that is no longer there.  They are therefore not
  /// computed rather than computed wrongly; P3b brings its own.
  FieldOnly,
};

/// Assemble and solve the device problem.  Throws on an inconsistent setup --
/// most importantly, MetallicReference with anything other than eps_r = 1.
DielectricSolution solve_dielectric(const DielectricSetup& s);

/// The same solve on a mesh the caller has already built.  `solve_dielectric`
/// is exactly this with `build_volume_mesh(s.geometry)`; it exists so that P3b
/// can hand in the same mesh with its free surface deformed onto a prescribed
/// meniscus, WITHOUT this file learning anything about menisci and without the
/// boundary conditions, the audit or the assembly being written a second time.
DielectricSolution solve_dielectric_on(DeviceVolumeMesh mesh, const DielectricSetup& s,
                                       DielectricDiagnostics diag = DielectricDiagnostics::Full);

/// Node roles alone, without solving.  Used by the audit test.
std::vector<NodeRole> node_roles(const DeviceVolumeMesh& m, const DielectricSetup& s);
BoundaryAudit audit_boundaries(const DeviceVolumeMesh& m, const std::vector<NodeRole>& role,
                               const DielectricSetup& s);

}  // namespace es
