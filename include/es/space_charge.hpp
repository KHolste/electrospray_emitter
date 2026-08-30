#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/axisym_fem.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P6 -- Poisson with free space charge, and the foundation for a PIC loop
// ===========================================================================
//
// WHAT THIS REPLACES.  The prototype's ring macroparticle model is DISABLED
// because its self-field diverges: the analytic potential of a charged ring is
// logarithmically singular on the ring, so a particle sitting on its own ring
// sees an infinite field.  Refining the mesh made it worse, not better.
//
// WHAT IS DONE INSTEAD, and why it is not the same mistake with a new name.
// The charge is not a singularity in the field solution at all.  It is
// DEPOSITED onto the mesh nodes with the element shape functions and enters the
// FINITE ELEMENT load vector,
//
//     f_a = sum_p q_p N_a(x_p) ,
//
// which is exactly the weak form of div(eps grad phi) = -rho for
// rho = sum_p q_p delta(x - x_p).  The solution of that discrete system is a
// piecewise-bilinear function: it is FINITE everywhere, including at the
// particle, and its gradient is finite too.  There is no singularity to
// regularise because the discretisation never creates one.
//
// THE PRICE, stated exactly, because an earlier version of this comment stated
// it wrongly.  The finite element solution is bounded ON A FIXED MESH.  It is
// NOT regularised: as h -> 0 the self-potential of a macroparticle grows
// without bound, and how fast depends on where the particle sits.
//
//   * OFF THE AXIS a macroparticle is a charged RING, whose exact potential is
//     logarithmically singular on the ring itself.  Smearing it over a cell of
//     size h cuts the logarithm off at h, so the self-potential grows like
//     ln(1/h) -- slowly, and without limit.
//   * ON THE AXIS the ring degenerates to a point charge, whose exact potential
//     is 1/d singular.  Cutting that off at h gives a self-potential that grows
//     like 1/h.
//
// Both are measured here (see self_potential_scaling), and an earlier comment
// in this file claimed the growth "must be like 1/h in the potential"
// everywhere.  That is only the on-axis law; off the axis the measurement shows
// the logarithm, which is the ring behaving like a ring.
//
// So "no divergence" is a statement about a FIXED MESH and nothing else.  Every
// output of this module that says a quantity is bounded now says at which mesh,
// and every statement about h -> 0 says whether the quantity converges.
//
// WHAT IS DONE ABOUT IT.  The self-force a macroparticle exerts on itself is
// exactly zero in the physics and is a pure artefact of the deposition in the
// discretisation.  Because the discrete problem is LINEAR, the artefact can be
// removed exactly rather than damped: see exclude_self_field() below.  Its cost
// is one extra solve per particle and that cost is reported, so nothing here
// pretends the mechanism scales to a production PIC loop.
//
// WHAT IS NOT DONE, and why: see pic_options() -- a typed record of the three
// candidate treatments with a verdict and the measurement behind each.
//
// WHAT IS NOT HERE:
//   * no self-consistent emission-PIC loop.  It is BLOCKED, and by two
//     independent reasons at once -- see pic_loop_status().  P5 is blocked, so
//     there is no physical source of particles; and the only self-field
//     treatment that survives scrutiny costs a solve per particle.  Any charge
//     distribution used here is PRESCRIBED and is labelled a test source;
//   * no particle motion -- that is P7;
//   * no magnetic field, no collisions, no space-charge-limited emission.

// ---------------------------------------------------------------------------

/// A macroparticle: a position in the meridian half-plane and the charge it
/// carries.  In an axisymmetric problem a particle off the axis represents a
/// RING of charge; `charge` is the total charge of that ring in coulomb, which
/// is what the weak form needs and what is conserved.
struct Macroparticle {
  Vec2 x;
  Real charge{0};   ///< [C], signed
};

struct DepositionResult {
  std::vector<Real> node_charge;   ///< [C] per node
  Real total_deposited{0};         ///< [C]
  Real total_particles{0};         ///< [C], the sum over the particles
  Real conservation_error{0};      ///< |deposited - particles| / |particles|
  Index n_outside{0};              ///< particles that fell outside the mesh
  Index n_on_axis{0};
  /// Largest |sum_a N_a - 1| over the deposited particles.  The partition of
  /// unity is what makes the deposition conservative; it is measured.
  Real partition_of_unity_error{0};
};

/// Deposit macroparticles onto the nodes of `m` with the bilinear shape
/// functions.  A particle outside the mesh is COUNTED and not deposited --
/// silently dropping charge would break the conservation statement.
DepositionResult deposit(const QuadMesh& m, const std::vector<Macroparticle>& p);

// ---------------------------------------------------------------------------

struct SpaceChargeSolution {
  Index n_nodes{0};
  std::vector<Real> phi;          ///< [V]
  std::vector<Real> phi_no_charge;///< [V], the same problem with rho = 0
  Real fem_residual{0};
  Real deposited_charge{0};
  Real conservation_error{0};
  /// max |phi - phi_no_charge| and the same for the recovered field.
  Real max_potential_shift{0};
  Real max_field_shift{0};
};

/// Solve the electrostatic problem on `m` with the given Dirichlet data and a
/// PRESCRIBED charge distribution, and solve it again with rho = 0 on the same
/// mesh and the same boundary data.  Returning both is deliberate: the
/// difference is the space-charge effect and nothing else, and it is what the
/// figures show.
///
/// `node_charge` may be empty (then only the Laplace problem is solved) and
/// `node_density` may be empty (then there is no volumetric source).  Both may
/// be given; they add.
SpaceChargeSolution solve_with_space_charge(const QuadMesh& m, const std::vector<Real>& eps_r,
                                           const std::vector<char>& active,
                                           const std::vector<char>& fixed,
                                           const std::vector<Real>& fixed_value,
                                           const std::vector<Real>& node_charge,
                                           const std::vector<Real>& node_density);

// ---------------------------------------------------------------------------
// The manufactured solution
// ---------------------------------------------------------------------------
//
// On the cylinder r in [0, R], z in [-L, L]:
//
//     phi(r,z) = phi0 (R^2 - r^2)(L^2 - z^2) / (R^2 L^2)
//
// vanishes on r = R and on z = +-L, is regular on the axis, and has
//
//     lap phi = -(2 phi0 / (R^2 L^2)) [ 2(L^2 - z^2) + (R^2 - r^2) ] ,
//     rho = -eps0 lap phi
//         =  (2 eps0 phi0 / (R^2 L^2)) [ 2(L^2 - z^2) + (R^2 - r^2) ] .
//
// Both are written out below so that the test compares the solver against an
// independently evaluated closed form rather than against itself.

Real manufactured_potential(Vec2 x, Real R, Real L, Real phi0);
Real manufactured_charge_density(Vec2 x, Real R, Real L, Real phi0);

/// A uniform cylinder mesh with `nr` x `nz` nodes over r in [0,R], z in [-L,L].
QuadMesh cylinder_mesh_symmetric(Real R, Real L, Index nr, Index nz);

// ---------------------------------------------------------------------------
// Field interpolation for particles
// ---------------------------------------------------------------------------
//
// The field a particle feels is the RECOVERED field of axisym_fem, evaluated by
// the same bilinear interpolation the deposition uses.  Using the same basis
// for deposition and interpolation is not a convenience: it is what keeps the
// self-force from acquiring a spurious direction on a uniform mesh.

/// E = -grad(phi), recovered at the nodes and interpolated bilinearly.
/// Returns {0,0} outside the mesh; the caller must check containment itself if
/// that matters.
Vec2 interpolated_field(const QuadMesh& m, const std::vector<Real>& phi,
                        const std::vector<Real>& eps_r, const std::vector<char>& active,
                        Vec2 x);

// ===========================================================================
// The self-field: measured honestly, and removed exactly
// ===========================================================================
//
// THE FIVE THINGS THAT MUST BE LOOKED AT SEPARATELY, because they behave
// differently and an aggregate statement about "the self-field" hides that:
//
//   1. CHARGE CONSERVATION of the deposition.  Exact to rounding, at every h.
//      It is a property of the partition of unity and has nothing to do with
//      the singularity.  -> DepositionResult.
//   2. THE FOREIGN FIELD of a particle, i.e. what a DIFFERENT particle at a
//      fixed distance d feels.  This converges under refinement as long as
//      d stays fixed, because away from the charge the discrete solution is
//      just a finite element approximation of a smooth function.
//      -> foreign_field_convergence().
//   3. THE SELF-POTENTIAL AND SELF-FIELD at the particle itself.  This does NOT
//      converge; see the header comment.  -> self_potential_scaling().
//   4. THE MESH DEPENDENCE OF THE SHAPE WIDTH.  The deposited cloud is one cell
//      wide, so its width IS h.  There is no physical length in it at all,
//      which is the reason 2. and 3. differ.  -> deposition_width().
//   5. THE SCALING WITH PARTICLES PER CELL.  The self-potential is proportional
//      to the charge of ONE macroparticle, so at fixed total charge it falls
//      like 1/N_p while the collective field stays put.
//      -> self_to_total_ratio().

/// The width of the deposited charge cloud of a single particle, measured as
/// the charge-weighted RMS distance of the receiving nodes from the particle.
/// It has no physical content: it is proportional to h and to nothing else,
/// which is exactly what makes a "smoothing width" a free parameter rather
/// than a model.
struct DepositionWidth {
  Real h{0};             ///< [m] radial cell size -- the reporting scale
  Real h_z{0};           ///< [m] axial cell size; these meshes are NOT square
  Real diagonal{0};      ///< [m] sqrt(h^2 + h_z^2), the largest distance inside a cell
  Real rms{0};           ///< [m], charge-weighted RMS node distance
  Real rms_over_h{0};    ///< the width in units of the ONLY length there is
  Real max_distance{0};  ///< [m], farthest node that received charge
  Real max_over_diagonal{0};  ///< never exceeds 1: the support is one cell
  Index n_nodes_receiving{0};
  /// True when the particle sits on a node to within rounding.  Then all the
  /// charge lands on that one node and the width is exactly zero -- a real and
  /// degenerate sub-case, reported rather than averaged away.
  bool on_node{false};
};
DepositionWidth deposition_width(const QuadMesh& m, const Macroparticle& p);

/// How the self-potential of ONE macroparticle behaves under refinement.
/// Both candidate laws are fitted and both residuals are reported, so that the
/// answer is read off the data instead of being assumed.
struct SelfPotentialScaling {
  std::vector<Real> h;              ///< cell size, one per level
  std::vector<Real> phi_self;       ///< [V] at the particle itself
  std::vector<Real> field_self;     ///< [V/m] magnitude at the particle
  /// phi ~ C h^-p : the fitted p and the residual of the fit in log-log.
  Real power_exponent{0}, power_residual{0};
  /// phi ~ a ln(1/h) + b : the fitted a and the residual.
  Real log_slope{0}, log_residual{0};
  bool prefers_logarithmic{false};  ///< which of the two fits the data better
  bool grows_under_refinement{false};
  Real growth_factor{0};            ///< last / first
  std::string verdict;
};
SelfPotentialScaling self_potential_scaling(Real R, Real L, Vec2 x_particle, Real charge,
                                            const std::vector<Index>& levels);

/// What a DIFFERENT particle at a fixed distance feels from a deposited one.
/// Unlike the self-field this converges, and the measured order is reported.
struct ForeignFieldConvergence {
  std::vector<Real> h;
  std::vector<Real> phi;       ///< [V] at the probe point
  std::vector<Real> field;     ///< [V/m] at the probe point
  Real distance{0};            ///< [m] probe distance, held fixed
  Real order_phi{0};           ///< observed order from the last three levels
  Real order_field{0};
  Real relative_change_last{0};///< |phi_n - phi_{n-1}| / |phi_n|
  bool converges{false};
  std::string verdict;
};
ForeignFieldConvergence foreign_field_convergence(Real R, Real L, Vec2 x_particle, Real charge,
                                                  Real distance,
                                                  const std::vector<Index>& levels);

/// The field at one particle with the contribution of its OWN deposited charge
/// removed.  The discrete problem is linear, so this is not a filter and not a
/// smoothing: it is a subtraction that is exact up to the linear solve, and the
/// superposition error is measured and returned rather than assumed.
struct SelfFieldExclusion {
  bool ok{false};
  Vec2 field_total{};      ///< what a naive PIC loop would hand the particle
  Vec2 field_self{};       ///< the particle's own deposited charge, alone
  Vec2 field_external{};   ///< total - self: the only one that may be used
  Real phi_total{0}, phi_self{0}, phi_external{0};
  /// |(boundary-only + sum of single-particle solutions) - total| / |total|,
  /// over the nodes.  Linearity is CHECKED, not trusted.
  Real superposition_error{0};
  Index solves{0};         ///< the cost of this call, stated
  std::string message;
};
SelfFieldExclusion exclude_self_field(const QuadMesh& m, const std::vector<Real>& eps_r,
                                      const std::vector<char>& active,
                                      const std::vector<char>& fixed,
                                      const std::vector<Real>& fixed_value,
                                      const std::vector<Macroparticle>& parts,
                                      std::size_t index);

/// The self-potential of one macroparticle divided by the potential the whole
/// prescribed charge produces, as a function of the number of macroparticles
/// the same total charge is split into.  This is the quantity a PIC convergence
/// argument would have to drive to zero.
struct SelfToTotalRatio {
  std::vector<Index> n_particles;
  std::vector<Real> phi_self;
  std::vector<Real> phi_total;
  std::vector<Real> ratio;
  /// ratio ~ N^-q over ALL points.  At small N the ratio is still near one --
  /// a single macroparticle IS its own total field -- so this underestimates
  /// the asymptotic law and is reported alongside it rather than instead.
  Real fitted_exponent{0};
  /// the same fit over the last three points, where the law has set in.
  Real fitted_exponent_asymptotic{0};
};
SelfToTotalRatio self_to_total_ratio(Real R, Real L, Real total_charge, Index level,
                                     const std::vector<Index>& counts);

// ---------------------------------------------------------------------------
// The three candidate treatments, with a verdict each
// ---------------------------------------------------------------------------

enum class PicOption {
  /// Subtract the particle's own deposited field from the field it feels.
  SelfFieldExclusion = 0,
  /// Give the macroparticle a physical width independent of h.
  FixedPhysicalShapeWidth,
  /// Let the macroparticle number grow with refinement so that the charge per
  /// macroparticle, and with it the self-field, goes to zero.
  ScaledMacroparticleNumber,
};
const char* to_string(PicOption o);

enum class PicOptionVerdict {
  /// Justified and implemented here, with tests that pin the justification.
  Implemented = 0,
  /// Rejected: it would introduce a free parameter with no physical origin.
  RejectedFreeParameter,
  /// Sound, and measured here, but not implemented -- it is a property of a
  /// PIC loop, and there is no PIC loop to run.
  MeasuredNotImplemented,
};
const char* to_string(PicOptionVerdict v);

struct PicOptionAssessment {
  PicOption option;
  PicOptionVerdict verdict;
  const char* what_it_does;
  const char* why;
  const char* evidence;   ///< which measurement in this file backs the verdict
};
const PicOptionAssessment* pic_options(std::size_t& n);

/// Why the self-consistent PIC loop is blocked.  Two independent reasons, so
/// that removing one does not silently unblock it.
struct PicLoopStatus {
  bool blocked{true};
  const char* reason_source{""};
  const char* reason_cost{""};
};
PicLoopStatus pic_loop_status();

}  // namespace es
