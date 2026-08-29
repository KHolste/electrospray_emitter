#pragma once
#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "es/types.hpp"

namespace es {

// ===========================================================================
// Axisymmetric finite elements for  div( eps(x) grad phi ) = 0
// ===========================================================================
//
// WHY A VOLUME METHOD AT ALL, GIVEN THAT THERE IS A WORKING BEM
//
// The boundary-integral solver in src/bem.cpp is an indirect single-layer
// formulation on conductors in vacuum.  It cannot represent a piecewise
// permittivity without becoming a different method (a multi-domain
// direct/indirect BEM with interface unknowns), and the P2b device is exactly
// that: a polymer body with eps_r != 1 between the liquid conductor and the
// electrode.  The BEM is therefore kept, untouched, as the INDEPENDENT
// reference for the one case both methods can solve -- eps_r = 1 everywhere,
// rho_f = 0 -- and the dielectric problem gets a volume method.
//
// WHAT THIS FILE IS AND IS NOT
//
// It is a generic axisymmetric Q1 finite element solver on a STRUCTURED
// quadrilateral grid of the meridian half-plane.  It knows nothing about
// emitters, extractors or materials; it takes a grid, a permittivity per cell,
// a Dirichlet mask per node and a far-field rule, and returns nodal potentials
// plus nodal reactions.  The device-specific part -- which cell is SU-8, which
// node is the metallisation -- lives in dielectric_device.hpp.
//
// It is NOT a general unstructured mesher or solver, and does not pretend to
// be.  The grid comes from volume_mesh.hpp, which builds it deterministically
// from the device parameters.
//
// THE AXISYMMETRIC WEIGHT
//
// Everything is an integral over the solid of revolution, so every element
// integral carries 2*pi*r:
//
//     a(phi, v) = \int_Omega eps grad(phi) . grad(v)  2 pi r  dr dz
//
// with grad = (d/dr, d/dz).  There is no 1/r term to add by hand and no special
// treatment at r = 0: the weight vanishes there, which IS the natural
// (symmetry) condition dphi/dr = 0 on the axis.  Getting the 2*pi*r wrong is
// the classic axisymmetric mistake, so tests/test_axisym_fem.cpp pins it
// against a coaxial capacitor, whose capacitance is a closed-form function of
// the radii and is wrong by a factor of order r if the weight is missing.

// ---------------------------------------------------------------------------
// Structured quadrilateral mesh of the meridian half-plane
// ---------------------------------------------------------------------------
//
// Logically a (nr x nz) grid: node (i, j) has index j*nr + i.  Cells are
// general quadrilaterals -- the r coordinate of a node may depend on BOTH
// indices, which is how a slanted material interface is made to lie exactly on
// a grid line.  The z coordinate depends on j ONLY ("level rows"); that is what
// makes point location exact and cheap, and validate_level_rows() checks it
// rather than assuming it.
struct QuadMesh {
  Index nr{0}, nz{0};
  std::vector<Vec2> nodes;  ///< size nr*nz

  Index node(Index i, Index j) const { return j * nr + i; }
  Index n_nodes() const { return nr * nz; }
  Index n_cells() const { return (nr - 1) * (nz - 1); }
  Index cell(Index i, Index j) const { return j * (nr - 1) + i; }
  Vec2 at(Index i, Index j) const { return nodes[static_cast<std::size_t>(node(i, j))]; }

  /// Corner nodes of cell (i, j), counter-clockwise in the (r, z) plane:
  /// (i,j), (i+1,j), (i+1,j+1), (i,j+1).
  std::array<Index, 4> cell_nodes(Index i, Index j) const;

  Real z_of_row(Index j) const { return nodes[static_cast<std::size_t>(node(0, j))].z; }

  Real cell_meridian_area(Index i, Index j) const;
  Real cell_revolved_volume(Index i, Index j) const;
  Real total_revolved_volume() const;

  /// Largest |ratio - 1| by which the z coordinate varies within a row.
  /// Must be zero for point location to be exact.
  bool validate_level_rows(Real tol_rel = 1e-12) const;
  /// Smallest cell Jacobian determinant over all Gauss points; must be > 0.
  Real min_jacobian() const;
};

// ---------------------------------------------------------------------------
// Far-field treatment of the open problem
// ---------------------------------------------------------------------------
//
// The device is not enclosed.  The physical condition is phi -> 0 at infinity
// with a system that carries NET charge (the emitter is at V, the extractor at
// ground, and nothing returns the emitter's charge).  A finite box therefore
// needs a condition that is consistent with the leading far-field term.
//
//   Asymptotic  the monopole condition.  For any bounded charge distribution,
//               phi = Q/(4 pi eps0 R) + O(R^-2) with R = |x - x0| measured from
//               a point inside the device, hence
//                     dphi/dn = -phi * (n . (x - x0)) / R^2 .
//               Imposed as a Robin term; the residual is dipole order, i.e. it
//               falls as R^-1 relative to the leading term rather than being
//               O(1) as a truncated Dirichlet condition is.  This is the
//               reference treatment.
//   Grounded    phi = 0 on the box.  A real, different device: one inside a
//               grounded enclosure at that distance.  Kept because it is the
//               honest way to show what the truncation is worth -- the two
//               bracket the answer, and the box-size study reports both.
//   Insulated   natural condition dphi/dn = 0.  A mirror plane, not an open
//               boundary.  Available for tests whose exact solution has it.
enum class FarField { Asymptotic = 0, Grounded, Insulated };
const char* to_string(FarField f);

// ---------------------------------------------------------------------------

struct AxisymProblem {
  const QuadMesh* mesh{nullptr};

  /// Relative permittivity per cell.  Cells with active = 0 are not part of the
  /// field domain at all (the interior of an ideal conductor).
  std::vector<Real> eps_r;
  std::vector<char> active;

  std::vector<char> fixed;        ///< per node: Dirichlet?
  std::vector<Real> fixed_value;  ///< per node: the value [V]

  FarField far_field{FarField::Asymptotic};
  /// Origin from which the monopole decay is measured.  Any point inside the
  /// device works to leading order; the box-size study reports the sensitivity.
  Vec2 far_field_origin{0.0, 0.0};
  /// Boundary segments carrying the far-field condition, as node pairs.
  /// Empty with FarField::Insulated.
  std::vector<std::array<Index, 2>> far_edges;

  /// Whether check() insists on level rows.
  ///
  /// The ASSEMBLY never needs them: every element integral is isoparametric and
  /// a general quadrilateral is fine.  What needs them is locate(), and through
  /// it every evaluation helper in this file.  A caller that has deformed the
  /// rows -- P3b moves the free surface off the plane z = 0 -- therefore turns
  /// the requirement off HERE and must then locate points itself; nothing else
  /// changes, and the default keeps every existing caller exactly as it was.
  bool require_level_rows{true};

  void check() const;  ///< throws on an inconsistent problem
};

// ---------------------------------------------------------------------------

enum class LinearSolver {
  Band = 0,  ///< exact symmetric band Cholesky; memory ~ N * nr * 8 byte
  Dense,     ///< dense LU via es::Matrix; only for tiny test systems
};

struct AxisymSolution {
  std::vector<Real> phi;       ///< nodal potential [V]
  /// Nodal reaction  R_i = sum_b K_ib phi_b  over the UNREDUCED system, in
  /// coulomb.  Zero at a free node (that is its equation).  At a Dirichlet node
  /// it is the charge that the constraint holds, so the total charge of a
  /// conductor is the sum of R over its nodes -- exact, and far more accurate
  /// than differentiating phi and integrating the result.
  std::vector<Real> reaction;
  Index n_nodes{0}, n_free{0};
  Index half_bandwidth{0};
  std::size_t factor_bytes{0};
  Real residual_inf{0.0};  ///< max |K phi - f| over the free equations [C]
};

/// Assemble and solve.  Throws std::runtime_error if the factorisation would
/// need more than `memory_cap_bytes` (default 2 GiB) -- a clear stop beats
/// swapping for ten minutes and then failing.
AxisymSolution solve_axisym(const AxisymProblem& p, LinearSolver which = LinearSolver::Band,
                            std::size_t memory_cap_bytes = 2ull << 30);

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

/// Cell containing `x`, plus its local coordinates xi, eta in [0,1]^2.
/// Returns false if `x` lies outside the mesh.  Exact for level-row meshes.
bool locate(const QuadMesh& m, Vec2 x, Index* i, Index* j, Real* xi, Real* eta);

/// Potential at a point (bilinear interpolation).  Throws if outside.
Real potential_at(const QuadMesh& m, const std::vector<Real>& phi, Vec2 x);

/// E = -grad(phi) at a point, evaluated inside the cell that contains it.
/// Q1 gradients jump across cell faces, so at a material interface the value
/// depends on which side the point falls -- which is the physically correct
/// behaviour and is what the interface tests exploit.
Vec2 field_at(const QuadMesh& m, const std::vector<Real>& phi, Vec2 x);

/// E inside a named cell at named local coordinates.  This is the one-sided
/// evaluator: it is how a field is read off ON a boundary from a chosen side.
Vec2 field_in_cell(const QuadMesh& m, const std::vector<Real>& phi, Index i, Index j, Real xi,
                   Real eta);

/// Potential inside a named cell at named local coordinates.  Evaluating the
/// same interface point from both sides must give the same number -- which it
/// does by construction, because the interface carries ONE nodal unknown.  The
/// measurement is kept because "by construction" is a claim about the code.
Real potential_in_cell(const QuadMesh& m, const std::vector<Real>& phi, Index i, Index j,
                       Real xi, Real eta);

// ---------------------------------------------------------------------------
// Recovered field
// ---------------------------------------------------------------------------
//
// A Q1 gradient taken at an arbitrary point inside a cell is only first-order
// accurate and jumps from cell to cell, so a mesh study built on it measures
// the jitter rather than the error -- on this device it moves the axial field
// by ten per cent between neighbouring levels while the potential is already
// converged to three digits.  The gradient at the CELL CENTRE is
// superconvergent for Q1; averaging those centre values into the nodes, with
// the cell volume as weight, gives a second-order, continuous field.
//
// Cells are averaged ONLY with cells of the same permittivity, and inactive
// cells are skipped.  Otherwise the recovery would smear the deliberate jump of
// E across a dielectric interface, which is precisely the quantity the
// interface test is about.  For a one-sided value AT an interface, use
// field_in_cell() -- these two answer different questions and are not
// interchangeable.

/// Recovered E at node (i, j), averaging only cells whose relative permittivity
/// equals `eps_select` (compare exactly; the values come from the same table).
Vec2 field_recovered_at_node(const QuadMesh& m, const std::vector<Real>& phi,
                             const std::vector<Real>& eps_r, const std::vector<char>& active,
                             Index i, Index j, Real eps_select);

/// Recovered E at a point: nodal recovery on the four corners of the containing
/// cell, using that cell's permittivity as the selector, then bilinear
/// interpolation.  Throws if the point is outside or in an inactive cell.
Vec2 field_recovered_at(const QuadMesh& m, const std::vector<Real>& phi,
                        const std::vector<Real>& eps_r, const std::vector<char>& active, Vec2 x);

/// Total charge held by the Dirichlet nodes flagged in `mask` [C].
Real charge_of(const AxisymSolution& s, const std::vector<char>& mask);

}  // namespace es
