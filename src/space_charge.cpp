#include "es/space_charge.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::eps0;

// ---------------------------------------------------------------------------

QuadMesh cylinder_mesh_symmetric(Real R, Real L, Index nr, Index nz) {
  if (!(R > 0.0) || !(L > 0.0) || nr < 3 || nz < 3)
    throw std::runtime_error("cylinder_mesh_symmetric: ungueltige Vorgaben");
  QuadMesh m;
  m.nr = nr;
  m.nz = nz;
  m.nodes.resize(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nz));
  for (Index j = 0; j < nz; ++j)
    for (Index i = 0; i < nr; ++i)
      m.nodes[static_cast<std::size_t>(j * nr + i)] =
          Vec2{R * static_cast<Real>(i) / static_cast<Real>(nr - 1),
               -L + 2.0 * L * static_cast<Real>(j) / static_cast<Real>(nz - 1)};
  return m;
}

Real manufactured_potential(Vec2 x, Real R, Real L, Real phi0) {
  return phi0 * (R * R - x.r * x.r) * (L * L - x.z * x.z) / (R * R * L * L);
}

Real manufactured_charge_density(Vec2 x, Real R, Real L, Real phi0) {
  // rho = -eps0 lap phi, written from the closed form, not from a difference.
  return (2.0 * eps0 * phi0 / (R * R * L * L)) *
         (2.0 * (L * L - x.z * x.z) + (R * R - x.r * x.r));
}

// ---------------------------------------------------------------------------

DepositionResult deposit(const QuadMesh& m, const std::vector<Macroparticle>& parts) {
  DepositionResult d;
  d.node_charge.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  for (const Macroparticle& p : parts) {
    d.total_particles += p.charge;
    Index i, j;
    Real xi, eta;
    if (!locate(m, p.x, &i, &j, &xi, &eta)) {
      ++d.n_outside;
      continue;   // NOT deposited, and counted -- dropping it silently would
                  // break the conservation statement below
    }
    if (std::abs(p.x.r) < 1e-15) ++d.n_on_axis;
    const Real N[4] = {(1 - xi) * (1 - eta), xi * (1 - eta), xi * eta, (1 - xi) * eta};
    d.partition_of_unity_error =
        std::max(d.partition_of_unity_error, std::abs(N[0] + N[1] + N[2] + N[3] - 1.0));
    const std::array<Index, 4> cn = m.cell_nodes(i, j);
    for (int a = 0; a < 4; ++a)
      d.node_charge[static_cast<std::size_t>(cn[a])] += p.charge * N[a];
  }
  for (Real q : d.node_charge) d.total_deposited += q;
  d.conservation_error =
      (std::abs(d.total_particles) > 0.0)
          ? std::abs(d.total_deposited - d.total_particles) / std::abs(d.total_particles)
          : std::abs(d.total_deposited);
  return d;
}

// ---------------------------------------------------------------------------

namespace {

AxisymSolution solve_one(const QuadMesh& m, const std::vector<Real>& eps_r,
                         const std::vector<char>& active, const std::vector<char>& fixed,
                         const std::vector<Real>& fixed_value,
                         const std::vector<Real>& node_charge,
                         const std::vector<Real>& node_density) {
  AxisymProblem p;
  p.mesh = &m;
  p.eps_r = eps_r;
  p.active = active;
  p.fixed = fixed;
  p.fixed_value = fixed_value;
  p.node_charge = node_charge;
  p.node_source_density = node_density;
  p.far_field = FarField::Insulated;
  p.check();
  return solve_axisym(p);
}

}  // namespace

SpaceChargeSolution solve_with_space_charge(const QuadMesh& m, const std::vector<Real>& eps_r,
                                           const std::vector<char>& active,
                                           const std::vector<char>& fixed,
                                           const std::vector<Real>& fixed_value,
                                           const std::vector<Real>& node_charge,
                                           const std::vector<Real>& node_density) {
  SpaceChargeSolution out;
  out.n_nodes = m.n_nodes();

  const AxisymSolution with =
      solve_one(m, eps_r, active, fixed, fixed_value, node_charge, node_density);
  const AxisymSolution without =
      solve_one(m, eps_r, active, fixed, fixed_value, {}, {});

  out.phi = with.phi;
  out.phi_no_charge = without.phi;
  out.fem_residual = with.residual_inf;
  for (Real q : node_charge) out.deposited_charge += q;

  for (std::size_t k = 0; k < out.phi.size(); ++k)
    out.max_potential_shift =
        std::max(out.max_potential_shift, std::abs(out.phi[k] - out.phi_no_charge[k]));

  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i) {
      const Vec2 a = field_recovered_at_node(m, out.phi, eps_r, active, i, j, 1.0);
      const Vec2 b = field_recovered_at_node(m, out.phi_no_charge, eps_r, active, i, j, 1.0);
      out.max_field_shift = std::max(out.max_field_shift, norm(a - b));
    }
  return out;
}

// ---------------------------------------------------------------------------

Vec2 interpolated_field(const QuadMesh& m, const std::vector<Real>& phi,
                        const std::vector<Real>& eps_r, const std::vector<char>& active,
                        Vec2 x) {
  Index i, j;
  Real xi, eta;
  if (!locate(m, x, &i, &j, &xi, &eta)) return {0.0, 0.0};
  const Vec2 e00 = field_recovered_at_node(m, phi, eps_r, active, i, j, 1.0);
  const Vec2 e10 = field_recovered_at_node(m, phi, eps_r, active, i + 1, j, 1.0);
  const Vec2 e11 = field_recovered_at_node(m, phi, eps_r, active, i + 1, j + 1, 1.0);
  const Vec2 e01 = field_recovered_at_node(m, phi, eps_r, active, i, j + 1, 1.0);
  const Real n0 = (1 - xi) * (1 - eta), n1 = xi * (1 - eta), n2 = xi * eta,
             n3 = (1 - xi) * eta;
  return n0 * e00 + n1 * e10 + n2 * e11 + n3 * e01;
}

}  // namespace es
