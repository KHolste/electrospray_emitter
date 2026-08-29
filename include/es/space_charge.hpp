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
// The price is stated rather than hidden: the self-field of a macroparticle
// depends on the mesh -- it is the field of a charge smeared over the
// neighbouring cells, not of a ring.  The peak grows as the mesh is refined,
// because a finer mesh smears it less.  That is measured (it must grow like
// 1/h in the potential, not like log(1/d) in the distance), and it is the
// reason a PIC loop must never let a particle feel its own deposited charge
// without a correction.  This file does not implement such a correction; it
// makes the effect measurable and says so.
//
// WHAT IS NOT HERE:
//   * no self-consistent emission-PIC loop.  P5 is blocked, so there is no
//     physical source of particles.  Any charge distribution used here is
//     PRESCRIBED and is labelled a test source;
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

}  // namespace es
