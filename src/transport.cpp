#include "es/transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::eps0;
using constants::pi;

namespace {

constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();

/// A uniform cylinder mesh: nr columns in r, nz rows in z, level rows.
QuadMesh cylinder_mesh(Real R, Real L, Index nr, Index nz) {
  QuadMesh m;
  m.nr = nr;
  m.nz = nz;
  m.nodes.resize(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nz));
  for (Index j = 0; j < nz; ++j)
    for (Index i = 0; i < nr; ++i)
      m.nodes[static_cast<std::size_t>(j * nr + i)] =
          Vec2{R * static_cast<Real>(i) / static_cast<Real>(nr - 1),
               L * static_cast<Real>(j) / static_cast<Real>(nz - 1)};
  return m;
}

}  // namespace

// ===========================================================================
// 1.  Fully developed pipe flow
// ===========================================================================

PipeFlowSolution solve_pipe_flow(Real R, Real L, Real mu, Real dpdz, Index nr, Index nz) {
  if (!(R > 0.0) || !(L > 0.0) || !(mu > 0.0))
    throw std::runtime_error("solve_pipe_flow: R, L und mu muessen positiv sein");
  if (nr < 3 || nz < 3) throw std::runtime_error("solve_pipe_flow: Netz zu grob");

  PipeFlowSolution out;
  out.radius = R;
  out.length = L;
  out.mu = mu;
  out.dpdz = dpdz;

  const QuadMesh m = cylinder_mesh(R, L, nr, nz);
  out.n_nodes = m.n_nodes();
  out.n_cells = m.n_cells();

  AxisymProblem p;
  p.mesh = &m;
  // THE COEFFICIENT IS THE VISCOSITY and the scale is one: see the header of
  // AxisymProblem.  Nothing about the assembly changes; only what the numbers
  // mean does.
  p.coefficient_scale = 1.0;
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), mu);
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  // -div(mu grad u) = -dp/dz
  p.cell_source.assign(static_cast<std::size_t>(m.n_cells()), -dpdz);
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  // No slip on the wall.  The two ends and the axis carry the NATURAL condition
  // du/dn = 0, which is exactly right here: the profile is z-independent, so
  // du/dz = 0 on the ends, and regularity gives du/dr = 0 on the axis.  Nothing
  // is imposed that the exact solution does not satisfy.
  for (Index j = 0; j < m.nz; ++j) {
    const Index n = m.node(m.nr - 1, j);
    p.fixed[static_cast<std::size_t>(n)] = 1;
    p.fixed_value[static_cast<std::size_t>(n)] = 0.0;
  }
  p.far_field = FarField::Insulated;
  p.check();

  const AxisymSolution sol = solve_axisym(p);
  out.u = sol.phi;
  out.fem_residual = sol.residual_inf;

  // --- the closed forms this is judged against ------------------------------
  const Real u0 = -dpdz * R * R / (4.0 * mu);
  out.centreline_closed_form = u0;
  out.flow_rate_closed_form = -dpdz * pi * std::pow(R, 4.0) / (8.0 * mu);

  out.centreline_velocity = out.u[static_cast<std::size_t>(m.node(0, m.nz / 2))];

  // Flow rate by integrating the nodal profile over one cross section with the
  // trapezoidal rule against 2 pi r dr.  It is a different operation from the
  // closed form and from the FEM assembly.
  {
    const Index j = m.nz / 2;
    Real q = 0.0;
    for (Index i = 0; i + 1 < m.nr; ++i) {
      const Real r0 = m.at(i, j).r, r1 = m.at(i + 1, j).r;
      const Real u0n = out.u[static_cast<std::size_t>(m.node(i, j))];
      const Real u1n = out.u[static_cast<std::size_t>(m.node(i + 1, j))];
      // int_{r0}^{r1} u(r) 2 pi r dr with u linear in r.
      q += 2.0 * pi * (r1 - r0) * ((2.0 * r0 * u0n + r0 * u1n + r1 * u0n + 2.0 * r1 * u1n) / 6.0);
    }
    out.flow_rate = q;
  }
  out.mean_velocity = out.flow_rate / (pi * R * R);

  // Wall shear stress from the last cell's gradient: tau_w = mu |du/dr| at r=R.
  {
    const Index j = m.nz / 2;
    const Real r0 = m.at(m.nr - 2, j).r, r1 = m.at(m.nr - 1, j).r;
    const Real ua = out.u[static_cast<std::size_t>(m.node(m.nr - 2, j))];
    const Real ub = out.u[static_cast<std::size_t>(m.node(m.nr - 1, j))];
    out.wall_shear_stress = mu * std::abs((ub - ua) / (r1 - r0));
  }

  // Profile error against the parabola, on the mid plane.
  {
    const Index j = m.nz / 2;
    Real worst = 0.0;
    for (Index i = 0; i < m.nr; ++i) {
      const Real r = m.at(i, j).r;
      const Real want = u0 * (1.0 - (r / R) * (r / R));
      const Real got = out.u[static_cast<std::size_t>(m.node(i, j))];
      if (std::abs(want) > 1e-12 * std::abs(u0))
        worst = std::max(worst, std::abs(got - want) / std::abs(u0));
    }
    out.max_profile_error = worst;
  }
  out.flow_rate_error = (std::abs(out.flow_rate_closed_form) > 0.0)
                            ? std::abs(out.flow_rate - out.flow_rate_closed_form) /
                                  std::abs(out.flow_rate_closed_form)
                            : kNaN;
  return out;
}

// ===========================================================================
// 2.  Charge transport
// ===========================================================================

Real charge_relaxation_time(Real eps_r, Real sigma) {
  if (!(eps_r > 0.0) || !(sigma > 0.0)) return kNaN;
  return eps0 * eps_r / sigma;
}

Real relaxed_charge_density(Real rho0, Real t, Real tau) {
  if (!(tau > 0.0)) return kNaN;
  return rho0 * std::exp(-t / tau);
}

const char* to_string(ConductorLimit v) {
  switch (v) {
    case ConductorLimit::NotAttempted: return "NotAttempted";
    case ConductorLimit::PerfectConductorJustified: return "PerfectConductorJustified";
    case ConductorLimit::FiniteConductivityRequired: return "FiniteConductivityRequired";
    case ConductorLimit::MissingMaterialData: return "MissingMaterialData";
  }
  return "?";
}

const char* explain(ConductorLimit v) {
  switch (v) {
    case ConductorLimit::NotAttempted:
      return "Der Grenzfall wurde nicht beurteilt.";
    case ConductorLimit::PerfectConductorJustified:
      return "Die Ladungsrelaxationszeit ist um mehr als den festgelegten Faktor kuerzer als "
             "die angegebene Prozesszeit.  Die Behandlung der Fluessigkeit als "
             "Aequipotentialflaeche ist FUER DIESE Prozesszeit gerechtfertigt -- und fuer "
             "keine andere.";
    case ConductorLimit::FiniteConductivityRequired:
      return "Die Ladungsrelaxationszeit ist vergleichbar mit der Prozesszeit oder laenger.  "
             "Die Behandlung als idealer Leiter ist dann eine Annahme und kein Grenzfall; "
             "das Problem braucht endliche Leitfaehigkeit.";
    case ConductorLimit::MissingMaterialData:
      return "Die Permittivitaet oder die Leitfaehigkeit fehlt.  tau ist nicht berechenbar, "
             "und es wird KEIN Ersatzwert gesetzt: ohne beide Zahlen ist der Grenzfall "
             "weder gerechtfertigt noch widerlegt, sondern unbekannt.";
  }
  return "?";
}

void RelaxationVerdict::print(std::FILE* out) const {
  std::fprintf(out, "  Ladungsrelaxation: %s\n", to_string(limit));
  std::fprintf(out, "    eps_r = %.6g [%s], sigma = %.6g S/m [%s]\n", eps_r,
               to_string(eps_status), sigma, to_string(sigma_status));
  if (std::isfinite(tau)) {
    std::fprintf(out, "    tau = eps0 eps_r / sigma = %.6e s\n", tau);
    std::fprintf(out, "    Prozesszeit %.6e s, Verhaeltnis t/tau = %.4g (verlangt > %.0f)\n",
                 process_time, ratio, transport::kPerfectConductorMargin);
  } else {
    std::fprintf(out, "    tau ist nicht berechenbar.\n");
  }
  std::fprintf(out, "    %s\n", message.c_str());
}

RelaxationVerdict judge_conductor_limit_explicit(Real eps_r, Real sigma, Real process_time) {
  RelaxationVerdict v;
  v.eps_r = eps_r;
  v.sigma = sigma;
  v.process_time = process_time;
  v.eps_status = (eps_r > 0.0) ? MaterialDataStatus::Literature
                               : MaterialDataStatus::MissingMaterialData;
  v.sigma_status = (sigma > 0.0) ? MaterialDataStatus::Literature
                                 : MaterialDataStatus::MissingMaterialData;
  v.tau = charge_relaxation_time(eps_r, sigma);
  if (!std::isfinite(v.tau)) {
    v.limit = ConductorLimit::MissingMaterialData;
    v.message = explain(v.limit);
    v.ratio = kNaN;
    return v;
  }
  v.ratio = process_time / v.tau;
  v.limit = (v.ratio > transport::kPerfectConductorMargin)
                ? ConductorLimit::PerfectConductorJustified
                : ConductorLimit::FiniteConductivityRequired;
  v.message = explain(v.limit);
  return v;
}

RelaxationVerdict judge_conductor_limit(const MaterialDataset& d, Real T, Real process_time) {
  const MaterialValue e = material_value(d, PropertyKind::RelativePermittivity, T);
  const MaterialValue s = material_value(d, PropertyKind::ElectricalConductivity, T);
  RelaxationVerdict v;
  v.process_time = process_time;
  v.eps_status = e.status;
  v.sigma_status = s.status;
  v.eps_r = e.usable() ? e.value : 0.0;
  v.sigma = s.usable() ? s.value : 0.0;
  if (!e.usable() || !s.usable()) {
    v.limit = ConductorLimit::MissingMaterialData;
    v.tau = kNaN;
    v.ratio = kNaN;
    v.message = std::string(explain(v.limit)) + "  ";
    if (!e.usable()) v.message += "eps_r: " + e.message + "  ";
    if (!s.usable()) v.message += "sigma: " + s.message;
    return v;
  }
  return judge_conductor_limit_explicit(v.eps_r, v.sigma, process_time);
}

// ===========================================================================
// 3.  Steady conduction in a cylinder
// ===========================================================================

ConductionSolution solve_cylinder_conduction(Real R, Real L, Real sigma, Real V, Index nr,
                                             Index nz) {
  if (!(R > 0.0) || !(L > 0.0) || !(sigma > 0.0))
    throw std::runtime_error("solve_cylinder_conduction: R, L und sigma muessen positiv sein");
  if (nr < 3 || nz < 3) throw std::runtime_error("solve_cylinder_conduction: Netz zu grob");

  ConductionSolution out;
  out.radius = R;
  out.length = L;
  out.sigma = sigma;
  out.voltage = V;

  const QuadMesh m = cylinder_mesh(R, L, nr, nz);
  out.n_nodes = m.n_nodes();

  AxisymProblem p;
  p.mesh = &m;
  // THE COEFFICIENT IS THE CONDUCTIVITY.  div(sigma grad phi) = 0 is the same
  // equation as div(eps grad phi) = 0; the nodal reaction is then a current in
  // ampere rather than a charge in coulomb.
  p.coefficient_scale = 1.0;
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), sigma);
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  for (Index i = 0; i < m.nr; ++i) {
    const Index a = m.node(i, 0), b = m.node(i, m.nz - 1);
    p.fixed[static_cast<std::size_t>(a)] = 1;
    p.fixed_value[static_cast<std::size_t>(a)] = 0.0;
    p.fixed[static_cast<std::size_t>(b)] = 1;
    p.fixed_value[static_cast<std::size_t>(b)] = V;
  }
  // The LATERAL surface carries the natural condition j . n = 0.  That is the
  // condition a non-emitting free surface must also carry: charge that crossed
  // it would have nowhere to go.
  p.far_field = FarField::Insulated;
  p.check();

  const AxisymSolution sol = solve_axisym(p);
  out.phi = sol.phi;
  out.fem_residual = sol.residual_inf;

  // The current through the driven end, from the OPERATOR reaction -- the same
  // route by which the electrostatic problem gets a conductor charge, and far
  // more accurate than differentiating phi and integrating the result.
  {
    std::vector<char> mask(static_cast<std::size_t>(m.n_nodes()), 0);
    for (Index i = 0; i < m.nr; ++i) mask[static_cast<std::size_t>(m.node(i, m.nz - 1))] = 1;
    out.current = std::abs(charge_of(sol, mask));
  }
  out.current_closed_form = sigma * pi * R * R * std::abs(V) / L;
  out.current_error = (out.current_closed_form > 0.0)
                          ? std::abs(out.current - out.current_closed_form) /
                                out.current_closed_form
                          : kNaN;
  out.resistance = (out.current > 0.0) ? std::abs(V) / out.current : kNaN;
  out.resistance_closed_form = L / (sigma * pi * R * R);

  // The exact solution is linear in z; check it node by node.
  {
    Real worst = 0.0;
    for (Index j = 0; j < m.nz; ++j)
      for (Index i = 0; i < m.nr; ++i) {
        const Real want = V * m.at(i, j).z / L;
        const Real got = out.phi[static_cast<std::size_t>(m.node(i, j))];
        worst = std::max(worst, std::abs(got - want));
      }
    out.max_potential_error = (std::abs(V) > 0.0) ? worst / std::abs(V) : worst;
  }

  // Lateral leakage: the radial current density on the outer column, against
  // the mean axial one.  In the exact solution it is identically zero.
  {
    const Real j_axial = sigma * std::abs(V) / L;
    Real worst = 0.0;
    for (Index j = 1; j + 1 < m.nz; ++j) {
      const Real r0 = m.at(m.nr - 2, j).r, r1 = m.at(m.nr - 1, j).r;
      const Real p0 = out.phi[static_cast<std::size_t>(m.node(m.nr - 2, j))];
      const Real p1 = out.phi[static_cast<std::size_t>(m.node(m.nr - 1, j))];
      worst = std::max(worst, sigma * std::abs((p1 - p0) / (r1 - r0)));
    }
    out.lateral_leakage = (j_axial > 0.0) ? worst / j_axial : kNaN;
  }
  return out;
}

}  // namespace es
