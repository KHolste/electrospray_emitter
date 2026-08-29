#pragma once
#include <array>
#include <vector>

#include "es/geometry.hpp"
#include "es/linalg.hpp"
#include "es/types.hpp"

namespace es {

// ---------------------------------------------------------------------------
// Axisymmetric Laplace kernel
// ---------------------------------------------------------------------------
// A surface charge density sigma [C/m^2] on a body of revolution produces
//
//     V(r,z) = \int sigma(s') G(r,z; r',z') ds'
//     G      = r' K(m) / (pi eps0 S),
//     S      = sqrt((r+r')^2 + (z-z')^2),   m = 4 r r' / S^2
//
// where K is the complete elliptic integral of the first kind (parameter
// convention).  G is the free-space Green's function integrated once around the
// azimuth, so the far field decays on its own -- no truncation boundary is
// needed.  Close to the ring, G -> ln(8 r / d) / (2 pi eps0): an integrable
// logarithmic singularity, handled by geometric panel refinement.

/// Potential kernel G at field point x due to a unit ring density at xp.
Real kernel_G(Vec2 x, Vec2 xp);

/// Gradient of G with respect to the FIELD point x (needed for E = -grad V).
Vec2 kernel_gradG(Vec2 x, Vec2 xp);

// ---------------------------------------------------------------------------
// Solver
// ---------------------------------------------------------------------------

/// Electrodes are groups of elements that share one applied potential.  The
/// meniscus is electrically part of the emitter: for ionic liquids the charge
/// relaxation time eps eps_r / K is ~1e-10 s, orders of magnitude below every
/// hydrodynamic time scale, so the free surface is an equipotential.
enum class Electrode : int { Emitter = 0, Extractor = 1, Collector = 2, Count = 3 };

Electrode electrode_of(Tag t);

class BemSolver {
 public:
  BemSolver() = default;
  explicit BemSolver(Mesh mesh) : mesh_(std::move(mesh)) {}

  /// Replacing the mesh invalidates EVERYTHING derived from the old one.
  /// Failing to clear the basis here meant that a later solve() silently
  /// superposed the previous mesh's basis vectors onto the new geometry --
  /// wrong by a factor of 30 when the two shapes differ substantially, and
  /// quietly wrong by a small amount when they nearly agree.
  void set_mesh(Mesh mesh) {
    mesh_ = std::move(mesh);
    factored_ = false;
    basis_.clear();
    sigma_.clear();
    applied_ = {{0.0, 0.0, 0.0}};
  }
  Mesh& mesh() { return mesh_; }
  const Mesh& mesh() const { return mesh_; }

  /// Build and LU-factor the influence matrix.  O(N^2) kernel work, O(N^3)
  /// factorisation; both are cheap at the N ~ 10^3 an axisymmetric mesh needs.
  void assemble();

  /// Solve for the unit-potential basis: sigma_k is the density produced by
  /// holding electrode k at 1 V and all others at 0 V.  Because the problem is
  /// linear, every subsequent voltage combination is then a vector sum -- a
  /// voltage sweep costs nothing beyond the first assembly.
  void solve_basis();

  /// Superpose the basis for a given set of electrode potentials [V].
  std::vector<Real> sigma_for(const std::array<Real, 3>& V) const;

  /// One-shot convenience: assemble (if needed), solve, store the result as the
  /// active density, and write the applied potentials back into the mesh so
  /// every consumer sees one consistent state.
  void solve(const std::array<Real, 3>& V);

  /// Potentials passed to the last solve() [V], indexed by Electrode.
  const std::array<Real, 3>& applied() const { return applied_; }

  /// An externally imposed potential sampled at the collocation points (beam
  /// space charge).  Subtracted from the right-hand side.  Invalidates the
  /// stored basis but not the factorisation.
  void set_external_potential(std::vector<Real> phi_ext);
  void clear_external_potential();

  // --- results --------------------------------------------------------------
  const std::vector<Real>& sigma() const { return sigma_; }

  /// Outward normal field on element i.  For a closed conductor the interior
  /// field vanishes, so E_n = sigma / eps0 exactly -- far more accurate than
  /// differentiating the potential numerically at the surface.
  Real En(Index i) const;

  /// Total charge on all elements carrying a given tag [C].
  Real charge_on(Tag t) const;
  Real total_charge() const;

  /// Self-capacitance of electrode `e` from the current solution: Q_e / V_e.
  /// Only meaningful for a single-conductor problem.
  Real capacitance(Electrode e = Electrode::Emitter) const;

  /// Field-point evaluation (valid anywhere in the vacuum, including close to
  /// the surface: the quadrature refines automatically).
  Real potential_at(Vec2 x) const;
  Vec2 field_at(Vec2 x) const;

  /// Peak |E_n| over elements with the given tag, and the element index.
  Real peak_field(Tag t, Index* which = nullptr) const;
  /// Peak |E_n| over the emitter-side electrode (metal + meniscus).
  Real peak_emitter_field(Index* which = nullptr) const;

  Index size() const { return mesh_.size(); }
  bool is_factored() const { return factored_; }

  /// Dump element-wise sigma, E_n and Maxwell stress.
  void write_surface_csv(const std::string& path, const std::string& header = {}) const;

 private:
  Mesh mesh_;
  Matrix lu_;
  std::vector<Index> piv_;
  bool factored_{false};
  std::vector<std::vector<Real>> basis_;  // per electrode
  std::vector<Real> sigma_;
  std::vector<Real> phi_ext_;
  std::array<Real, 3> applied_{{0.0, 0.0, 0.0}};
};

// ---------------------------------------------------------------------------
// Quadrature over one element (exposed for testing)
// ---------------------------------------------------------------------------

/// \int_{element} G(x; x') ds', with panel refinement toward the point on the
/// element closest to x.  `self` forces the singular treatment.
Real integrate_G(const Element& el, Vec2 x, bool self);

/// \int_{element} grad_x G(x; x') ds'.
Vec2 integrate_gradG(const Element& el, Vec2 x);

}  // namespace es
