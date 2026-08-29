#pragma once
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "es/capillary.hpp"
#include "es/dielectric_device.hpp"
#include "es/liquid.hpp"
#include "es/types.hpp"
#include "es/volume_mesh.hpp"

namespace es {

// ===========================================================================
// P3b -- self-consistent static electro-capillary equilibrium
// ===========================================================================
//
// WHAT IS SOLVED
//
//     gamma * kappa(s) = delta_p_exit + p_M(s),     p_M = eps0 E_n^2 / 2,
//
// with kappa the sum of the principal curvatures as in P3a, and E_n the
// VACUUM-SIDE normal field on the free surface of the liquid, which is treated
// as an ideal conductor at V_emitter.  The electrostatics is the P2c dielectric
// problem, unchanged: div(eps grad phi) = 0, polymer bodies polarise and are
// never electrodes, the metal film on the extractor carrier is at V_extractor.
//
// WHERE THE SIGN COMES FROM, AND WHY p_M IS ALWAYS OUTWARD
//
// On the surface of a conductor the field is normal, E = E_n n with n the
// outward (vacuum-pointing) normal, and the surface charge is sigma = eps0 E_n.
// The traction on the surface is sigma * <E>, with <E> the mean of the field on
// the two sides -- zero inside, E_n n outside -- hence
//
//     t = sigma * (E_n/2) n = (eps0 E_n^2 / 2) n = p_M n,   p_M >= 0.
//
// It pulls the surface OUTWARD whatever the sign of E_n, so it enters the
// curvature equation exactly like an increase of the internal pressure, and it
// cannot change sign with the polarity.  Two consequences are tested rather
// than asserted: at delta_p_exit = 0 and V != 0 the meniscus must bulge OUT
// (h > 0), and reversing every applied potential must leave the shape identical
// while reversing E_n.
//
// WHAT IS NOT IN HERE -- and must not be read into any number
//
//   * no emission and no emitted current;
//   * no finite liquid conductivity: the liquid is an ideal equipotential
//     conductor, the P2b/P2c contract, and no charge relaxation is modelled;
//   * no flow, no viscosity, no feed impedance;
//   * no space charge -- the vacuum carries no free charge;
//   * no time dependence and NO DYNAMIC STABILITY.  A branch that ends is a
//     branch this solver could not continue; it is not a stability limit;
//   * no Taylor-cone onset and no cone-jet.  The word "onset" does not appear
//     in any status, figure or report of this phase;
//   * gravity stays out, on the Bond number computed in P3a;
//   * delta_p_exit remains an INPUT.  P3b determines neither the reservoir
//     pressure nor the feed pressure.
//
// THE EDGE.  The contact line sits on an unrounded conductor/dielectric/vacuum
// edge, where the field is singular and no pointwise value converges.  Whether
// the load can be used at all is therefore decided BEFORE any coupled solve, by
// the gate below, on measurements: the exponent is fitted, the integrated force
// is followed under refinement and under a shrinking exclusion zone, and the
// weak (segment-averaged) projection is checked for convergence.  Nothing is
// clipped, no maximum is hard-coded, no exclusion zone is chosen freely and no
// edge radius is invented.

// ---------------------------------------------------------------------------
// The prescribed free surface
// ---------------------------------------------------------------------------

/// A meridian curve from the apex on the axis to the pinned contact line, in
/// the P3a convention: nodes uniform in arclength, tangent angle psi per node.
///
/// It is kept in ARCLENGTH form, never as a height function z(r).  The mesher
/// below needs z at a given radius and gets it by inverting r(s) with the cubic
/// Hermite that the tangent angles make exact to fourth order -- which is a
/// local evaluation of an arclength curve, not a change of parametrisation.
struct FreeSurface {
  std::vector<Vec2> nodes;   ///< apex -> contact line
  std::vector<Real> psi;     ///< tangent angle [rad]
  Real contact_radius{0}, contact_z{0};
  Real apex_height{0}, arclength{0};
  bool flat{true};

  static FreeSurface flat_surface(Real a, Real z_contact);
  static FreeSurface from(const CapillaryMeniscus& m);

  /// z of the surface at radius r in [0, contact_radius].
  Real z_at_radius(Real r) const;
  /// Tangent angle at radius r; the outward normal is (sin psi, cos psi).
  Real psi_at_radius(Real r) const;
  /// Arclength from the apex to radius r.
  Real s_at_radius(Real r) const;
  /// Revolved volume between the surface and the plane z = contact_z [m^3].
  Real revolved_volume() const;
  Real revolved_area() const;
};

// ---------------------------------------------------------------------------
// The moving mesh
// ---------------------------------------------------------------------------
//
// HOW THE FREE SURFACE IS PUT INTO THE P2c MESH, AND WHY IT STAYS VALID
//
// The P2c mesh is a tensor grid in (r_reference, z) with a purely RADIAL warp
// that is the identity inside the bore, so the columns i <= i_bore are vertical
// lines and the row j_tip is the flat liquid surface.  P3b adds an AXIAL warp
// on those columns only:
//
//   * row j_tip is carried from z = 0 to z = z_meniscus(r_i);
//   * two rows of the existing grid, the first below -kBandFactor*a and the
//     first above +kBandFactor*a, are held FIXED and act as the band edges;
//   * between them the map is piecewise linear in z and continuous.
//
// Outside the bore nothing moves at all, and at r = r_bore the meniscus height
// is zero by the pinning condition, so the two halves of the mesh agree there
// exactly.  Cells cannot invert: r depends on i alone in the deformed columns,
// so the Jacobian determinant is (dr/dxi)(dz/deta) with both factors positive
// by construction.  It is measured all the same.
//
// The CELL REGIONS are unchanged -- P2c assigns them from indices alone, so the
// cells below row j_tip inside the bore are liquid wherever that row is put.
// The boundary conditions, the audit and the assembly are therefore the P2c
// ones, untouched.
//
// The band edges are EXISTING grid rows, not new ones.  That keeps the flat
// case bitwise identical to the P2c mesh -- which is what makes the zero-field
// cross-check against P2c a check and not a comparison of two different meshes.

namespace meniscus_mesh {
/// The axial warp band reaches at least this multiple of the pinning radius
/// above and below the exit plane.  The largest apex height a pinned cap can
/// have is exactly a, so 1.5 leaves half a bore radius of margin on each side;
/// nothing is clipped and no shape is refused for lack of band.
inline constexpr Real kBandFactor = 1.5;
}  // namespace meniscus_mesh

struct MeniscusMeshQuality {
  Real min_jacobian{0};          ///< over all cells and Gauss points; must be > 0
  Real max_cell_aspect{0};       ///< worst edge-length ratio in the deformed band
  Real max_shear{0};             ///< worst |dz across a cell| / cell height
  Index inverted_cells{0};       ///< must be 0
  Index deformed_cells{0};
  Real band_z_lo{0}, band_z_hi{0};
  Real contact_radius_error{0};  ///< |r(node at the contact) - a| / a
  Real contact_z_error{0};       ///< |z(node at the contact) - z_contact| / a
  Real apex_error{0};            ///< |z(apex node) - h| / a
  Real surface_error{0};         ///< max |z(node) - z_surface(r)| / a over the row
  /// Liquid volume of the mesh against the closed form (bore column plus the
  /// revolved meniscus volume).  The mesh volume is a sum of cell volumes; the
  /// reference comes from the surface, so the two are independent.
  Real liquid_volume_mesh{0}, liquid_volume_reference{0}, liquid_volume_error{0};
  bool ok() const;
  void print(std::FILE* out) const;
};

struct MeniscusMesh {
  DeviceVolumeMesh device;   ///< the P2c mesh with the free surface moved
  FreeSurface surface;
  MeniscusMeshQuality quality;
  Index j_surface{0};        ///< the row that carries the free surface
  Index i_contact{0};        ///< the column of the pinned contact line
};

/// Build the P2c mesh and move its free surface onto `surface`.
/// Throws std::runtime_error if the surface does not fit the band or if a cell
/// would invert -- it does not return a mesh it knows to be broken.
MeniscusMesh build_meniscus_mesh(const DielectricDeviceParameters& p, const FreeSurface& surface);

// --- evaluation on the deformed mesh ---------------------------------------
//
// The generic locate() in axisym_fem.hpp assumes level rows, which the deformed
// mesh does not have inside the bore.  These three do the same job for it and
// are exact: outside the bore the rows ARE level, and inside it the columns are
// vertical, so both halves reduce to a bisection plus a linear map.

bool locate_meniscus(const MeniscusMesh& m, Vec2 x, Index* i, Index* j, Real* xi, Real* eta);
Real potential_at_meniscus(const MeniscusMesh& m, const std::vector<Real>& phi, Vec2 x);
Vec2 field_recovered_at_meniscus(const MeniscusMesh& m, const std::vector<Real>& phi,
                                 const std::vector<Real>& eps_r, const std::vector<char>& active,
                                 Vec2 x);

// ---------------------------------------------------------------------------
// The Maxwell load on the free surface
// ---------------------------------------------------------------------------
//
// Nodal values are the one-sided vacuum limit, recovered from the vacuum cells
// only.  They are reported and plotted, and near the edge they are NOT a
// converged quantity -- that is the point of the gate.
//
// The quantity the coupling uses is the SEGMENT one: for every surface segment
// the electric normal force is integrated with the one-sided field,
//
//     dF_k = int_segment (eps0 E_n^2 / 2) 2 pi r ds ,
//
// and the projected pressure of that segment is dF_k divided by its revolved
// area.  Summing the segments returns the total force exactly, so the
// projection is conservative by construction; near the edge the segment mean of
// an integrable singularity stays finite while the pointwise value does not.

struct MaxwellLoad {
  // --- nodal, one-sided from the vacuum ------------------------------------
  std::vector<Real> node_r, node_z, node_s, node_tau;
  std::vector<Real> node_En;        ///< [V/m], signed
  std::vector<Real> node_pM;        ///< [Pa], NOT converged in the last cells
  std::vector<Real> node_d_edge;    ///< distance to the contact line along the surface [m]
  /// |E_t| / |E| at the node.  The surface is an equipotential, so this is zero
  /// in the exact solution; it is a discretisation error measure, not physics.
  std::vector<Real> node_tangential_fraction;

  // --- per segment, conservative -------------------------------------------
  std::vector<Real> seg_tau0, seg_tau1, seg_tau_mid, seg_d_mid;
  std::vector<Real> seg_area;       ///< revolved area of the segment [m^2]
  std::vector<Real> seg_force;      ///< integrated normal force on it [N]
  std::vector<Real> seg_pressure;   ///< seg_force / seg_area [Pa]

  Real total_force{0};              ///< sum of seg_force [N]
  Real axial_force{0};              ///< integral of p_M n_z dA [N]
  Real gamma_over_a{0};             ///< the pressure scale the load is judged against

  /// Integrated force from the segments whose midpoint lies farther than
  /// `d_exclude` from the contact line [N].
  Real force_beyond(Real d_exclude) const;
  /// The conservative piecewise-constant load at normalised arclength tau [Pa].
  Real pressure_at_tau(Real tau) const;
  /// Largest segment pressure whose midpoint is farther than d from the edge.
  Real pressure_beyond(Real d_exclude) const;
  /// Segment pressure at a distance d from the contact line, linearly
  /// interpolated between segment midpoints.  It is how two different meshes
  /// are compared at the SAME physical place, since their segments do not line
  /// up.  Returns 0 outside the sampled range.
  Real pressure_at_distance(Real d) const;
  bool empty() const { return seg_pressure.empty(); }
};

/// Evaluate the load on the free surface of `m` from the solved field.
MaxwellLoad maxwell_load(const MeniscusMesh& m, const DielectricSolution& sol,
                         Real gamma_over_a);

/// Assemble the conservative segment quantities of `out` from its already
/// filled nodal arrays (node_r, node_z, node_tau, node_d_edge, node_pM) and the
/// surface geometry.  maxwell_load() calls it; so do the manufactured loads of
/// load_projection.hpp, so that an audit exercises the production quadrature
/// rather than a second copy of it.  It knows nothing of where the nodal values
/// came from.
void assemble_load_segments(MaxwellLoad& out, const FreeSurface& fs);

// ---------------------------------------------------------------------------
// The coupling gate
// ---------------------------------------------------------------------------
//
// Tolerances fixed BEFORE the measurement, and not loosened afterwards.
namespace edge_gate {
/// Exclusion distances, in units of the pinning radius.  Three of them, halving,
/// so that the contribution of the edge zone can be seen to converge instead of
/// being asserted to.
inline constexpr Real kExclusionCoarse = 0.10;
inline constexpr Real kExclusionMid = 0.05;
inline constexpr Real kExclusionFine = 0.025;
/// The segment-averaged load away from the edge must be mesh converged: largest
/// relative change over the segments with d >= kExclusionCoarse * a, between the
/// two finest levels.
inline constexpr Real kTolEdgeFarLoad = 2.0e-2;
/// The total integrated normal force must be mesh converged to this.
inline constexpr Real kTolTotalForce = 5.0e-2;
/// DECLARED BEFORE THE MEASUREMENT, MEASURED, AND FOUND TO BE THE WRONG TEST.
///
/// It was meant to show that the singular contribution is summable: halve the
/// exclusion distance on the finest mesh and see the integrated force settle.
/// What the number actually compares is the force in two DIFFERENT regions of
/// the surface, and the force contained in the annulus between two exclusion
/// radii is a physical quantity that does not become small under refinement and
/// must not.  With the measured exponent the missing force behaves as
/// d^(1+beta) = d^0.56, so the change per halving cannot fall below about a
/// third of the previous one no matter how fine the mesh is.
///
/// The number is still measured and still reported -- it is
/// EdgeGateResult::measured_exclusion_change, and it does not meet this bound.
/// The decision is carried by kTolLimitAgreement below, which tests summability
/// the way it can be tested.  Nothing was loosened: this bound stands, it is
/// reported as not met, and the defect is written down in
/// docs/10_electrocapillary_model.md rather than repaired in silence.
inline constexpr Real kTolExclusion = 5.0e-2;

/// The summability test that replaces it, at the SAME tolerance as the total
/// force, so that no new number was invented after the fact.
///
/// Two extrapolations of the same limit, from data that share nothing:
///   * over the MESH -- Aitken on the total force of the three finest levels;
///   * over the EXCLUSION DISTANCE -- a least-squares fit of
///     F(d0) = F_inf - C d0^(1+beta) with the independently fitted beta.
/// If the edge contribution were not summable the second would not exist and
/// the two could not agree.  They must agree to this.
inline constexpr Real kTolLimitAgreement = 5.0e-2;
/// p_M ~ d^beta must be integrable against the revolved area element, i.e.
/// beta > -1.  Measured by a fit, never assumed.
inline constexpr Real kMinExponent = -1.0;
}  // namespace edge_gate

struct EdgeStudyPoint {
  int mesh_level{0};
  Index n_nodes{0};
  Index n_surface_segments{0};
  Real smallest_d{0};          ///< distance of the innermost segment midpoint [m]
  Real total_force{0};
  Real force_coarse{0}, force_mid{0}, force_fine{0};  ///< with the three exclusions
  Real peak_node_pM{0};        ///< the pointwise edge value: reported, never used
  Real max_tangential_fraction{0};
  Real fit_exponent{0};        ///< beta in p_M ~ d^beta
  Real fit_r2{0};
  Real fit_d_lo{0}, fit_d_hi{0};
  Index fit_points{0};
};

enum class GateVerdict {
  NotAttempted = 0,
  Passed,                 ///< the weak projection may be used
  FailedTotalForce,       ///< the integrated force does not converge with the mesh
  FailedEdgeFarLoad,      ///< the load away from the edge does not converge
  FailedExclusion,        ///< the two extrapolations of the limit force disagree
  FailedNotIntegrable,    ///< the fitted exponent is <= -1
  FailedNoData,
};
const char* to_string(GateVerdict v);
const char* explain(GateVerdict v);

struct EdgeGateResult {
  std::string shape_tag;
  Real Pi{0};
  std::vector<EdgeStudyPoint> levels;
  /// The full load of every level, so that the figures can show E_n(d) and
  /// p_M(d) themselves instead of a summary of them.
  std::vector<MaxwellLoad> loads;
  GateVerdict verdict{GateVerdict::NotAttempted};
  Real measured_total_force_change{0};
  Real measured_edge_far_change{0};
  /// Reported, not decisive -- see edge_gate::kTolExclusion.
  Real measured_exclusion_change{0};
  /// The two independent extrapolations of the limit force [N] and their
  /// relative disagreement.
  Real limit_force_mesh{0};
  Real limit_force_exclusion{0};
  Real measured_limit_agreement{0};
  Real fitted_exponent{0};
  /// The conducting-wedge value for the measured contact tangent angle, quoted
  /// as a reference only: it ignores the dielectric on the other side of the
  /// edge, so it is not a prediction of the number that was fitted.
  Real wedge_reference_exponent{0};
  std::string note;
};

/// Run the gate for ONE prescribed shape over a list of mesh levels.
///
/// Nothing is coupled here and no shape is solved for: the surface is
/// prescribed, the field is solved on it, and the question asked is only
/// whether the resulting surface load is a usable quantity at all.
EdgeGateResult run_edge_gate(const DielectricDeviceParameters& base,
                             const DielectricMaterials& materials, Real V_emitter,
                             Real V_extractor, Metallisation metallisation, FarField far_field,
                             const FreeSurface& surface, const std::string& tag, Real Pi,
                             const std::vector<int>& levels, Real gamma_over_a,
                             std::size_t memory_cap_bytes = (2ull << 30));

// ---------------------------------------------------------------------------
// The coupled solve
// ---------------------------------------------------------------------------

enum class CouplingStatus {
  NotAttempted = 0,
  Converged,
  ElectrostaticFailure,            ///< the field solve threw or did not assemble
  MeshInvalid,                     ///< the free surface could not be meshed
  MechanicalResidualNotConverged,  ///< the fixed point did not settle
  ContinuationStepTooSmall,        ///< the branch could not be continued further
  CapillaryRangeExceeded,          ///< no shape pinned at a exists for this load
  EdgeLoadNotWellPosed,            ///< the gate did not pass; nothing was coupled
  MultipleSolutionsDetected,       ///< several shapes satisfy the same data
};
const char* to_string(CouplingStatus s);
const char* explain(CouplingStatus s);
inline bool is_usable(CouplingStatus s) { return s == CouplingStatus::Converged; }

/// Convergence bounds of the fixed point.  Fixed before the measurement.
namespace coupling {
inline constexpr Real kTolShape = 1.0e-6;      ///< max node motion / a
inline constexpr Real kTolLoad = 1.0e-4;       ///< max |dp_M| / (gamma/a)
inline constexpr Real kTolContact = 1.0e-12;   ///< |r_contact - a| / a
/// Young-Laplace residual * a/gamma, judged AWAY FROM THE EDGE.
///
/// DEFECT FOUND AND CORRECTED.  It was first applied to the maximum over the
/// whole surface.  That contradicts what the gate had already established: the
/// pointwise load at the contact line does not converge, so a pointwise
/// residual there cannot either, and the continuation stopped at 733 V for a
/// reason that was a property of the load representation rather than of the
/// solution.  The bound is unchanged; what changed is that it is applied to the
/// quantity it was meant for.  The whole-surface value is still measured and
/// still reported.
inline constexpr Real kTolMechanical = 1.0e-3;
inline constexpr int kMaxIterations = 60;
inline constexpr Real kRelaxation = 0.5;       ///< under-relaxation of the load
}  // namespace coupling

struct CoupledPoint {
  CouplingStatus status{CouplingStatus::NotAttempted};
  std::string message;

  Real V_emitter{0}, V_extractor{0}, delta_p_exit{0};
  int mesh_level{0};

  CapillaryMeniscus shape;
  MaxwellLoad load;
  Real apex_height{0};
  Real max_curvature{0};          ///< max |kappa| along the surface [1/m]
  Real E_edge_far{0};             ///< |E| at the surface point d = 0.25 a from the edge
  Real E_apex{0};                 ///< normal field on the axis [V/m]
  Real total_force{0};
  Real liquid_volume{0}, surface_area{0};

  // --- what the iteration did -----------------------------------------------
  int iterations{0};
  Real final_shape_change{0};
  Real final_load_change{0};
  /// Young-Laplace residual * a/gamma, from the node coordinates alone.
  ///
  /// TWO NUMBERS, because they answer different questions.  `mechanical_residual`
  /// is the maximum over the WHOLE surface, including the last nodes before the
  /// pinned edge; there the applied load is the projection of a field the gate
  /// declared non-convergent pointwise, and the geometric curvature estimate
  /// straddles the steepest part of it, so that number measures the load
  /// representation and not the solution.  `mechanical_residual_edge_far` is the
  /// maximum over the nodes farther than edge_gate::kExclusionCoarse * a from
  /// the contact line, and it is the one the convergence bound is applied to.
  /// Both are always reported.
  Real mechanical_residual{0};
  Real mechanical_residual_edge_far{0};
  Real contact_error{0};
  Real fem_residual{0};           ///< max |K phi - f| over the free equations [C]
  Real fem_residual_relative{0};  ///< the same against |Q_emitter|
  Real min_jacobian{0};
  int crossings{1};
};

struct CoupledRequest {
  DielectricDeviceParameters geometry;
  DielectricMaterials materials;
  LiquidProperties liquid;
  Real V_emitter{0}, V_extractor{0};
  Real delta_p_exit{0};
  Metallisation metallisation{Metallisation::FrontAndAperture};
  FarField far_field{FarField::Asymptotic};
  Real relaxation{coupling::kRelaxation};
  /// Handed to the field solve; see DielectricSetup::memory_cap_bytes.
  std::size_t memory_cap_bytes{2ull << 30};
  /// Shape to start from.  Empty = the P3a solution for delta_p_exit, which is
  /// the field-free member of this branch.
  const CapillaryMeniscus* initial_shape{nullptr};
  bool verbose{false};
};

/// One operating point, solved to the fixed point of shape and load.
CoupledPoint solve_coupled(const CoupledRequest& q);

// ---------------------------------------------------------------------------
// Continuation over the voltage
// ---------------------------------------------------------------------------
//
// The branch followed is the one that contains the field-free P3a solution.  It
// is followed by raising |V| in steps, each started from the previous converged
// shape; a step that fails is halved, and the continuation stops when the step
// would fall below the floor below.  A point that did not converge is NEVER
// stored as a solution.
//
// The end of the branch is where THIS SOLVER stopped.  It is not an emission
// onset, not a Taylor-cone onset and not a stability limit; none of those is
// computed anywhere in this phase.
namespace continuation {
inline constexpr Real kFirstStep = 250.0;     ///< [V]
inline constexpr Real kMinStep = 15.625;      ///< [V]; kFirstStep / 2^4
inline constexpr Real kGrowth = 1.5;          ///< step growth after an easy point
}  // namespace continuation

struct ContinuationResult {
  std::vector<CoupledPoint> points;   ///< converged points only, in order of |V|
  CouplingStatus end_status{CouplingStatus::NotAttempted};
  std::string end_message;
  Real last_converged_voltage{0};
  Real first_failed_voltage{0};
  int steps_attempted{0}, steps_rejected{0};
};

ContinuationResult continue_over_voltage(CoupledRequest q, Real V_max);

}  // namespace es
