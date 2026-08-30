#include "es/space_charge.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::eps0;

namespace {
constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();
}  // namespace

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


// ===========================================================================
// The self-field: measured honestly, and removed exactly
// ===========================================================================

namespace {

/// A grounded cylinder, the same one the P6 run uses.  It lives here so that
/// the diagnostics below are library code that a test can call directly rather
/// than something only an application can reproduce.
struct GroundedBox {
  QuadMesh mesh;
  std::vector<Real> eps_r;
  std::vector<char> active, fixed;
  std::vector<Real> fixed_value;
};

GroundedBox grounded_box(Real R, Real L, Index nr, Index nz) {
  GroundedBox b;
  b.mesh = cylinder_mesh_symmetric(R, L, nr, nz);
  b.eps_r.assign(static_cast<std::size_t>(b.mesh.n_cells()), 1.0);
  b.active.assign(static_cast<std::size_t>(b.mesh.n_cells()), 1);
  b.fixed.assign(static_cast<std::size_t>(b.mesh.n_nodes()), 0);
  b.fixed_value.assign(static_cast<std::size_t>(b.mesh.n_nodes()), 0.0);
  for (Index j = 0; j < b.mesh.nz; ++j)
    b.fixed[static_cast<std::size_t>(b.mesh.node(b.mesh.nr - 1, j))] = 1;
  for (Index i = 0; i < b.mesh.nr; ++i) {
    b.fixed[static_cast<std::size_t>(b.mesh.node(i, 0))] = 1;
    b.fixed[static_cast<std::size_t>(b.mesh.node(i, b.mesh.nz - 1))] = 1;
  }
  return b;
}

/// Least squares fit y = a x + b, returning a and the RMS residual.
void linfit(const std::vector<Real>& x, const std::vector<Real>& y, Real* slope,
            Real* residual) {
  const std::size_t n = x.size();
  Real sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (std::size_t k = 0; k < n; ++k) { sx += x[k]; sy += y[k]; sxx += x[k] * x[k];
                                        sxy += x[k] * y[k]; }
  const Real den = static_cast<Real>(n) * sxx - sx * sx;
  const Real a = (std::abs(den) > 0.0) ? (static_cast<Real>(n) * sxy - sx * sy) / den : 0.0;
  const Real b = (sy - a * sx) / static_cast<Real>(n);
  Real r = 0;
  for (std::size_t k = 0; k < n; ++k) { const Real d = y[k] - (a * x[k] + b); r += d * d; }
  *slope = a;
  *residual = std::sqrt(r / static_cast<Real>(n));
}

}  // namespace

DepositionWidth deposition_width(const QuadMesh& m, const Macroparticle& p) {
  DepositionWidth w;
  w.h = (m.nr > 1) ? (m.nodes[static_cast<std::size_t>(m.node(1, 0))].r -
                      m.nodes[static_cast<std::size_t>(m.node(0, 0))].r)
                   : 0.0;
  w.h_z = (m.nz > 1) ? (m.nodes[static_cast<std::size_t>(m.node(0, 1))].z -
                        m.nodes[static_cast<std::size_t>(m.node(0, 0))].z)
                     : 0.0;
  w.diagonal = std::sqrt(w.h * w.h + w.h_z * w.h_z);
  const DepositionResult d = deposit(m, {p});
  Real sum = 0.0, sum2 = 0.0;
  for (std::size_t a = 0; a < d.node_charge.size(); ++a) {
    const Real q = d.node_charge[a];
    if (std::abs(q) <= 0.0) continue;
    ++w.n_nodes_receiving;
    const Vec2 xa = m.nodes[a];
    const Real dr = xa.r - p.x.r, dz = xa.z - p.x.z;
    const Real dist = std::sqrt(dr * dr + dz * dz);
    w.max_distance = std::max(w.max_distance, dist);
    sum += std::abs(q);
    sum2 += std::abs(q) * dist * dist;
  }
  w.rms = (sum > 0.0) ? std::sqrt(sum2 / sum) : 0.0;
  w.rms_over_h = (w.h > 0.0) ? w.rms / w.h : 0.0;
  w.max_over_diagonal = (w.diagonal > 0.0) ? w.max_distance / w.diagonal : 0.0;
  // "on a node" means the deposition put essentially everything on one node.
  // The tolerance is loose on purpose: the node coordinates are computed, so a
  // particle placed at a node position lands within rounding of it, not on it.
  w.on_node = w.rms <= 1.0e-6 * w.diagonal;
  return w;
}

SelfPotentialScaling self_potential_scaling(Real R, Real L, Vec2 xp, Real charge,
                                            const std::vector<Index>& levels) {
  SelfPotentialScaling out;
  for (Index n : levels) {
    const GroundedBox b = grounded_box(R, L, n, 2 * n - 1);
    const DepositionResult d = deposit(b.mesh, {{xp, charge}});
    const SpaceChargeSolution sol = solve_with_space_charge(
        b.mesh, b.eps_r, b.active, b.fixed, b.fixed_value, d.node_charge, {});
    out.h.push_back(R / static_cast<Real>(n - 1));
    out.phi_self.push_back(std::abs(potential_at(b.mesh, sol.phi, xp)));
    out.field_self.push_back(
        norm(interpolated_field(b.mesh, sol.phi, b.eps_r, b.active, xp)));
  }
  if (out.h.size() >= 3) {
    std::vector<Real> lh, lp, invlog;
    for (std::size_t k = 0; k < out.h.size(); ++k) {
      lh.push_back(std::log(out.h[k]));
      lp.push_back(std::log(out.phi_self[k]));
      invlog.push_back(-std::log(out.h[k]));   // ln(1/h)
    }
    Real slope = 0, res = 0;
    linfit(lh, lp, &slope, &res);
    out.power_exponent = -slope;               // phi ~ h^slope, so phi ~ h^-p
    out.power_residual = res;
    linfit(invlog, out.phi_self, &slope, &res);
    out.log_slope = slope;
    // Compare the two fits on the SAME footing: relative RMS residual.
    Real mean_phi = 0, mean_logphi = 0;
    for (std::size_t k = 0; k < out.phi_self.size(); ++k) {
      mean_phi += out.phi_self[k];
      mean_logphi += std::abs(lp[k]);
    }
    mean_phi /= static_cast<Real>(out.phi_self.size());
    mean_logphi /= static_cast<Real>(lp.size());
    out.log_residual = (mean_phi > 0.0) ? res / mean_phi : res;
    out.power_residual = (mean_logphi > 0.0) ? out.power_residual / mean_logphi
                                             : out.power_residual;
    out.prefers_logarithmic = out.log_residual < out.power_residual;
  }
  out.growth_factor = (!out.phi_self.empty() && out.phi_self.front() > 0.0)
                          ? out.phi_self.back() / out.phi_self.front()
                          : 0.0;
  out.grows_under_refinement = out.growth_factor > 1.0;
  if (out.grows_under_refinement)
    out.verdict = out.prefers_logarithmic
                      ? "waechst unter Verfeinerung, logarithmisch -- das Ringverhalten "
                        "ausserhalb der Achse.  KEINE Konvergenz fuer h -> 0."
                      : "waechst unter Verfeinerung, als Potenz von 1/h -- das Punktladungs"
                        "verhalten auf der Achse.  KEINE Konvergenz fuer h -> 0.";
  else
    out.verdict = "waechst unter Verfeinerung NICHT -- das waere zu erklaeren, nicht zu "
                  "feiern.";
  return out;
}

ForeignFieldConvergence foreign_field_convergence(Real R, Real L, Vec2 xp, Real charge,
                                                  Real distance,
                                                  const std::vector<Index>& levels) {
  ForeignFieldConvergence out;
  out.distance = distance;
  const Vec2 probe{xp.r + distance, xp.z};
  for (Index n : levels) {
    const GroundedBox b = grounded_box(R, L, n, 2 * n - 1);
    const DepositionResult d = deposit(b.mesh, {{xp, charge}});
    const SpaceChargeSolution sol = solve_with_space_charge(
        b.mesh, b.eps_r, b.active, b.fixed, b.fixed_value, d.node_charge, {});
    out.h.push_back(R / static_cast<Real>(n - 1));
    out.phi.push_back(potential_at(b.mesh, sol.phi, probe));
    out.field.push_back(norm(interpolated_field(b.mesh, sol.phi, b.eps_r, b.active, probe)));
  }
  const std::size_t n = out.h.size();
  if (n >= 3) {
    auto order = [&](const std::vector<Real>& v) {
      const Real d1 = v[n - 3] - v[n - 2], d2 = v[n - 2] - v[n - 1];
      if (!(std::abs(d2) > 0.0) || d1 / d2 <= 0.0) return kNaN;
      return std::log(d1 / d2) / std::log(out.h[n - 3] / out.h[n - 2]);
    };
    out.order_phi = order(out.phi);
    out.order_field = order(out.field);
    out.relative_change_last = (std::abs(out.phi[n - 1]) > 0.0)
                                   ? std::abs(out.phi[n - 1] - out.phi[n - 2]) /
                                         std::abs(out.phi[n - 1])
                                   : kNaN;
    out.converges = std::isfinite(out.relative_change_last) &&
                    out.relative_change_last < 0.05;
  }
  out.verdict = out.converges
                    ? "konvergiert bei FESTEM Abstand: abseits der Ladung ist die diskrete "
                      "Loesung eine Finite-Elemente-Naeherung einer glatten Funktion."
                    : "konvergiert bei diesem Abstand nicht -- der Abstand ist vermutlich "
                      "nicht gross gegen die Netzweite.";
  return out;
}

SelfFieldExclusion exclude_self_field(const QuadMesh& m, const std::vector<Real>& eps_r,
                                      const std::vector<char>& active,
                                      const std::vector<char>& fixed,
                                      const std::vector<Real>& fixed_value,
                                      const std::vector<Macroparticle>& parts,
                                      std::size_t index) {
  SelfFieldExclusion out;
  if (index >= parts.size()) {
    out.message = "exclude_self_field: der Index zeigt auf kein Teilchen.";
    return out;
  }
  const Vec2 xp = parts[index].x;

  // 1. everything together, with the real boundary data: what a naive loop uses
  const DepositionResult d_all = deposit(m, parts);
  const SpaceChargeSolution s_all =
      solve_with_space_charge(m, eps_r, active, fixed, fixed_value, d_all.node_charge, {});
  ++out.solves;

  // 2. this particle ALONE, with HOMOGENEOUS boundary data.  Homogeneous is not
  //    a detail: the boundary contribution belongs to the external field and
  //    must not be subtracted with the self-field.
  const std::vector<Real> zero_bc(fixed_value.size(), 0.0);
  const DepositionResult d_one = deposit(m, {parts[index]});
  const SpaceChargeSolution s_one =
      solve_with_space_charge(m, eps_r, active, fixed, zero_bc, d_one.node_charge, {});
  ++out.solves;

  // 3. linearity is CHECKED, not trusted: boundary-only plus every particle
  //    alone must reproduce the full solution.
  {
    const std::vector<Real> no_charge;
    const SpaceChargeSolution s_bc =
        solve_with_space_charge(m, eps_r, active, fixed, fixed_value, no_charge, {});
    ++out.solves;
    std::vector<Real> sum = s_bc.phi;
    for (std::size_t k = 0; k < parts.size(); ++k) {
      const DepositionResult dk = deposit(m, {parts[k]});
      const SpaceChargeSolution sk =
          solve_with_space_charge(m, eps_r, active, fixed, zero_bc, dk.node_charge, {});
      ++out.solves;
      for (std::size_t a = 0; a < sum.size(); ++a) sum[a] += sk.phi[a];
    }
    Real num = 0.0, den = 0.0;
    for (std::size_t a = 0; a < sum.size(); ++a) {
      num = std::max(num, std::abs(sum[a] - s_all.phi[a]));
      den = std::max(den, std::abs(s_all.phi[a]));
    }
    out.superposition_error = (den > 0.0) ? num / den : num;
  }

  out.phi_total = potential_at(m, s_all.phi, xp);
  out.phi_self = potential_at(m, s_one.phi, xp);
  out.phi_external = out.phi_total - out.phi_self;
  out.field_total = interpolated_field(m, s_all.phi, eps_r, active, xp);
  out.field_self = interpolated_field(m, s_one.phi, eps_r, active, xp);
  out.field_external = Vec2{out.field_total.r - out.field_self.r,
                            out.field_total.z - out.field_self.z};
  out.ok = true;
  out.message = "Das Selbstfeld ist SUBTRAHIERT, nicht gedaempft: die diskrete Aufgabe ist "
                "linear, und die Ueberlagerung ist geprueft.  Kosten: eine zusaetzliche "
                "Loesung je Teilchen.";
  return out;
}

SelfToTotalRatio self_to_total_ratio(Real R, Real L, Real total_charge, Index level,
                                     const std::vector<Index>& counts) {
  SelfToTotalRatio out;
  const GroundedBox b = grounded_box(R, L, level, 2 * level - 1);
  for (Index n : counts) {
    // The same total charge, split into n macroparticles on a short line.  The
    // line is fixed, so the COLLECTIVE field is the same at every n and only
    // the charge per macroparticle changes.  That is what isolates the effect.
    std::vector<Macroparticle> ps;
    const Real q = total_charge / static_cast<Real>(n);
    for (Index k = 0; k < n; ++k) {
      const Real t = (n == 1) ? 0.5 : static_cast<Real>(k) / static_cast<Real>(n - 1);
      ps.push_back({{0.35 * R, -0.25 * L + 0.5 * L * t}, q});
    }
    const DepositionResult d_all = deposit(b.mesh, ps);
    const SpaceChargeSolution s_all = solve_with_space_charge(
        b.mesh, b.eps_r, b.active, b.fixed, b.fixed_value, d_all.node_charge, {});
    const DepositionResult d_one = deposit(b.mesh, {ps[0]});
    const SpaceChargeSolution s_one = solve_with_space_charge(
        b.mesh, b.eps_r, b.active, b.fixed, b.fixed_value, d_one.node_charge, {});
    const Real pt = std::abs(potential_at(b.mesh, s_all.phi, ps[0].x));
    const Real po = std::abs(potential_at(b.mesh, s_one.phi, ps[0].x));
    out.n_particles.push_back(n);
    out.phi_total.push_back(pt);
    out.phi_self.push_back(po);
    out.ratio.push_back((pt > 0.0) ? po / pt : kNaN);
  }
  if (out.ratio.size() >= 2) {
    std::vector<Real> ln, lr;
    for (std::size_t k = 0; k < out.ratio.size(); ++k)
      if (std::isfinite(out.ratio[k]) && out.ratio[k] > 0.0) {
        ln.push_back(std::log(static_cast<Real>(out.n_particles[k])));
        lr.push_back(std::log(out.ratio[k]));
      }
    if (ln.size() >= 2) {
      Real slope = 0, res = 0;
      linfit(ln, lr, &slope, &res);
      out.fitted_exponent = -slope;
    }
    if (ln.size() >= 3) {
      const std::vector<Real> ln3(ln.end() - 3, ln.end());
      const std::vector<Real> lr3(lr.end() - 3, lr.end());
      Real slope = 0, res = 0;
      linfit(ln3, lr3, &slope, &res);
      out.fitted_exponent_asymptotic = -slope;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------

const char* to_string(PicOption o) {
  switch (o) {
    case PicOption::SelfFieldExclusion: return "self_field_exclusion";
    case PicOption::FixedPhysicalShapeWidth: return "fixed_physical_shape_width";
    case PicOption::ScaledMacroparticleNumber: return "scaled_macroparticle_number";
  }
  return "?";
}

const char* to_string(PicOptionVerdict v) {
  switch (v) {
    case PicOptionVerdict::Implemented: return "implemented";
    case PicOptionVerdict::RejectedFreeParameter: return "rejected_free_parameter";
    case PicOptionVerdict::MeasuredNotImplemented: return "measured_not_implemented";
  }
  return "?";
}

namespace {
const PicOptionAssessment kPicOptions[] = {
    {PicOption::SelfFieldExclusion, PicOptionVerdict::Implemented,
     "Zieht vom Feld am Teilchen das Feld seiner EIGENEN deponierten Ladung ab, gewonnen "
     "aus einer zweiten Loesung mit nur diesem Teilchen und homogenen Randwerten.",
     "Physikalisch uebt eine Ladung auf sich selbst keine Kraft aus; das deponierte "
     "Selbstfeld ist reines Diskretisierungsartefakt und traegt zudem eine Richtung, die "
     "allein von der Lage im Zellinneren abhaengt.  Die diskrete Aufgabe ist LINEAR, also "
     "laesst sich das Artefakt exakt abziehen statt daempfen -- keine Glaettungsbreite, "
     "kein Filter, kein frei gewaehlter Parameter.  Der Preis ist eine zusaetzliche Loesung "
     "je Teilchen und wird berichtet, statt verschwiegen zu werden.",
     "exclude_self_field(): die Ueberlagerung wird gegen die volle Loesung geprueft, und "
     "die scheinbare Selbstkraft ist danach exakt null statt netzabhaengig."},
    {PicOption::FixedPhysicalShapeWidth, PicOptionVerdict::RejectedFreeParameter,
     "Gaebe dem Makropartikel eine feste physikalische Breite a, unabhaengig von h, und "
     "deponierte dessen verschmierte Dichte.  Das Selbstpotential waere dann fuer h -> 0 "
     "beschraenkt.",
     "VERWORFEN.  Fuer dieses Modell gibt es keine physikalische Laenge, die a festlegte: "
     "ein Makropartikel steht hier fuer eine Delta-Ladung und nicht fuer eine Wolke "
     "bekannter Groesse; eine Debye-Laenge gibt es im Vakuumstrahl nicht, und die Groesse "
     "eines Ions ist um Groessenordnungen kleiner als jedes Netz.  a waere ein FREI "
     "GEWAEHLTER Glaettungsparameter, der genau die Groesse festlegte, die er verbergen "
     "soll.  Beschraenktheit, die man sich aussucht, ist keine Konvergenz.",
     "deposition_width(): die Breite der deponierten Wolke ist exakt proportional zu h und "
     "enthaelt keine andere Laenge -- rms/h ist ueber alle Netzstufen konstant.  Es gibt "
     "also nichts, woraus a folgen wuerde."},
    {PicOption::ScaledMacroparticleNumber, PicOptionVerdict::MeasuredNotImplemented,
     "Liesse bei Verfeinerung die Zahl der Makropartikel so wachsen, dass die Ladung je "
     "Makropartikel und damit das Selbstfeld gegen null geht, waehrend das kollektive Feld "
     "stehen bleibt.",
     "Mathematisch tragfaehig und der uebliche PIC-Konvergenzweg.  Es ist aber eine "
     "Eigenschaft einer PIC-SCHLEIFE und keine Eigenschaft dieses Moduls: es gibt hier "
     "keine Schleife, die man verfeinern koennte, weil P5 blockiert ist und damit keine "
     "physikalische Teilchenquelle existiert.  Der Zusammenhang wird deshalb GEMESSEN und "
     "nicht implementiert -- eine implementierte Konvergenzstrategie ohne Rechnung, die "
     "sie benutzt, waere eine Behauptung.",
     "self_to_total_ratio(): das Verhaeltnis phi_self/phi_total faellt mit der Zahl der "
     "Makropartikel wie N^-q mit gemessenem q nahe 1, bei gleichem Gesamtstrom und "
     "gleichem Netz."},
};
}  // namespace

const PicOptionAssessment* pic_options(std::size_t& n) {
  n = sizeof(kPicOptions) / sizeof(kPicOptions[0]);
  return kPicOptions;
}

PicLoopStatus pic_loop_status() {
  PicLoopStatus s;
  s.blocked = true;
  s.reason_source =
      "P5 ist blockiert: es gibt keine belegte Emissionsrate und damit keine physikalische "
      "Teilchenquelle.  Eine selbstkonsistente Schleife haette nichts zu emittieren.";
  s.reason_cost =
      "Die einzige Selbstfeldbehandlung, die dieser Punkt rechtfertigen kann -- der exakte "
      "Abzug --, kostet eine zusaetzliche Loesung je Teilchen.  Das ist fuer eine Diagnose "
      "richtig und fuer eine Produktionsschleife untragbar.  Beide Gruende sind unabhaengig: "
      "faellt einer weg, bleibt der andere.";
  return s;
}

}  // namespace es
