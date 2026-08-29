// Axisymmetric finite elements for div(eps grad phi) = 0.
//
// Everything here is checked against a CLOSED-FORM solution or against a second,
// independent computation.  Two properties get the most attention, because they
// are the two that a plausible-looking wrong implementation gets wrong:
//
//   * the 2*pi*r weight.  Dropping it leaves a solver that still converges, still
//     looks smooth and is wrong by a factor of order r.  The coaxial capacitor
//     pins it: its capacitance is a closed-form function of the radii, and no
//     planar formulation reproduces it.
//   * the jump condition at a material interface.  phi is continuous because the
//     interface carries one nodal unknown; eps*dphi/dn is continuous only if the
//     assembly is right, so that is the check with content.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/axisym_fem.hpp"
#include "es/constants.hpp"

using namespace es;
using constants::eps0;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-56s got=%-14.8g want=%-14.8g rel=%-10.3g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

/// For quantities whose exact value is zero, a relative bound is meaningless.
static void check_small(const char* what, Real got, Real bound) {
  const bool ok = std::abs(got) <= bound;
  if (!ok) ++failures;
  std::printf("  %-56s |%.3e| <= %.3e %s\n", what, got, bound, ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-56s %s\n", what, ok ? "OK" : "FAIL");
}

// ---------------------------------------------------------------------------

static QuadMesh rect_mesh(const std::vector<Real>& r, const std::vector<Real>& z) {
  QuadMesh m;
  m.nr = static_cast<Index>(r.size());
  m.nz = static_cast<Index>(z.size());
  m.nodes.resize(static_cast<std::size_t>(m.n_nodes()));
  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i)
      m.nodes[static_cast<std::size_t>(m.node(i, j))] = {r[static_cast<std::size_t>(i)],
                                                         z[static_cast<std::size_t>(j)]};
  return m;
}

static std::vector<Real> linspace(Real a, Real b, Index n) {
  std::vector<Real> v(static_cast<std::size_t>(n));
  for (Index k = 0; k < n; ++k)
    v[static_cast<std::size_t>(k)] = a + (b - a) * static_cast<Real>(k) / (n - 1);
  v.front() = a;
  v.back() = b;
  return v;
}

/// Logarithmically spaced radii from a to b with `c` hit exactly.
static std::vector<Real> logspace_through(Real a, Real c, Real b, Index n1, Index n2) {
  std::vector<Real> v;
  for (Index k = 0; k <= n1; ++k)
    v.push_back(a * std::pow(c / a, static_cast<Real>(k) / n1));
  v.back() = c;
  for (Index k = 1; k <= n2; ++k)
    v.push_back(c * std::pow(b / c, static_cast<Real>(k) / n2));
  v.back() = b;
  return v;
}

// ==========================================================================
// 1.  Coaxial capacitor with two radial dielectric layers
// ==========================================================================
//
//   phi(a) = V, phi(b) = 0, dphi/dz = 0 on both end faces (the natural
//   condition), eps_r = e1 on a<r<c and e2 on c<r<b.
//
//   D_r r = K = const,  E_r = K/(eps r),  so
//     V = K [ ln(c/a)/(eps0 e1) + ln(b/c)/(eps0 e2) ]
//     Q = 2 pi height K,   C = Q/V.
//
// This one test covers three of the required checks at once: the layered
// dielectric, the 2*pi*r weighting against a known axisymmetric geometry, and
// the continuity of phi and of the normal flux density at the interface.
struct Coax {
  Real a{1.0e-3}, c{3.0e-3}, b{1.0e-2}, height{2.0e-3};
  Real e1{2.5}, e2{6.0}, V{100.0};

  Real K() const {
    return V / (std::log(c / a) / (eps0 * e1) + std::log(b / c) / (eps0 * e2));
  }
  Real phi(Real r) const {
    if (r <= c) return V - K() * std::log(r / a) / (eps0 * e1);
    return V - K() * (std::log(c / a) / (eps0 * e1) + std::log(r / c) / (eps0 * e2));
  }
  Real capacitance() const { return 2.0 * pi * height * K() / V; }
  Real Er(Real r, Real e) const { return K() / (eps0 * e * r); }
};

static void solve_coax(const Coax& cx, Index n1, Index n2, AxisymProblem& p, QuadMesh& m,
                       AxisymSolution& s, std::vector<char>& inner_mask) {
  const std::vector<Real> r = logspace_through(cx.a, cx.c, cx.b, n1, n2);
  const std::vector<Real> z = linspace(0.0, cx.height, 5);
  m = rect_mesh(r, z);
  p = AxisymProblem{};
  p.mesh = &m;
  p.far_field = FarField::Insulated;
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), 1.0);
  for (Index j = 0; j + 1 < m.nz; ++j)
    for (Index i = 0; i + 1 < m.nr; ++i)
      p.eps_r[static_cast<std::size_t>(m.cell(i, j))] =
          (0.5 * (r[static_cast<std::size_t>(i)] + r[static_cast<std::size_t>(i + 1)]) < cx.c)
              ? cx.e1
              : cx.e2;
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  inner_mask.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  for (Index j = 0; j < m.nz; ++j) {
    p.fixed[static_cast<std::size_t>(m.node(0, j))] = 1;
    p.fixed_value[static_cast<std::size_t>(m.node(0, j))] = cx.V;
    inner_mask[static_cast<std::size_t>(m.node(0, j))] = 1;
    p.fixed[static_cast<std::size_t>(m.node(m.nr - 1, j))] = 1;
  }
  s = solve_axisym(p);
}

static void test_coaxial_layered() {
  std::printf("\n=== 1. Koaxialkondensator mit zwei Dielektrikumsschichten ===\n");
  const Coax cx;
  QuadMesh m;
  AxisymProblem p;
  AxisymSolution s;
  std::vector<char> inner;
  solve_coax(cx, 120, 160, p, m, s, inner);

  // Potential at three interior radii.
  Real worst = 0.0;
  for (Index i = 1; i + 1 < m.nr; ++i) {
    const Real r = m.at(i, 2).r;
    const Real got = s.phi[static_cast<std::size_t>(m.node(i, 2))];
    worst = std::max(worst, std::abs(got - cx.phi(r)) / cx.V);
  }
  check_small("groesster Potentialfehler, relativ zu V", worst, 2.0e-5);

  // Capacitance from the nodal reactions -- this is the 2*pi*r check.
  const Real Q = charge_of(s, inner);
  check("Kapazitaet aus den Knotenreaktionen", Q / cx.V, cx.capacitance(), 5.0e-5);

  // A planar (Cartesian) assembly would give this instead; the test is only
  // meaningful if the two differ by a lot, so assert that too.
  const Real planar = eps0 * cx.height / ((cx.c - cx.a) / cx.e1 + (cx.b - cx.c) / cx.e2);
  expect("die ebene Formel unterscheidet sich um mehr als Faktor 2",
         std::abs(planar - cx.capacitance()) / cx.capacitance() > 1.0);

  // Continuity across the material interface, evaluated one cell either side.
  const Index ic = [&] {
    for (Index i = 0; i < m.nr; ++i)
      if (m.at(i, 0).r == cx.c) return i;
    return Index{-1};
  }();
  expect("die Materialgrenze ist exakt ein Netzknoten", ic > 0);
  // At the CELL CENTRE the Q1 gradient is superconvergent, so it may be pinned
  // against the closed form tightly.
  const Real rc_in = 0.5 * (m.at(ic - 1, 0).r + m.at(ic, 0).r);
  const Real rc_out = 0.5 * (m.at(ic, 0).r + m.at(ic + 1, 0).r);
  check("E_r im Zellmittelpunkt innen", field_in_cell(m, s.phi, ic - 1, 2, 0.5, 0.5).r,
        cx.Er(rc_in, cx.e1), 5.0e-5);
  check("E_r im Zellmittelpunkt aussen", field_in_cell(m, s.phi, ic, 2, 0.5, 0.5).r,
        cx.Er(rc_out, cx.e2), 5.0e-5);

  // ON the interface the value is a one-sided edge extrapolation of a Q1
  // gradient and therefore only first order.  The physical statement is that the
  // two sides agree; the numerical statement is that the residual mismatch is an
  // O(h) discretisation error, so it is bounded AND shown to fall.
  const Vec2 E_in = field_in_cell(m, s.phi, ic - 1, 2, 1.0, 0.5);
  const Vec2 E_out = field_in_cell(m, s.phi, ic, 2, 0.0, 0.5);
  const Real Dn_in = eps0 * cx.e1 * E_in.r, Dn_out = eps0 * cx.e2 * E_out.r;
  const Real dn_err = std::abs(Dn_in - Dn_out) / std::abs(Dn_in);
  check_small("D_n-Sprung an der Materialgrenze (einseitig, O(h))", dn_err, 1.5e-2);
  check("E_r springt um das Permittivitaetsverhaeltnis", E_in.r / E_out.r, cx.e2 / cx.e1,
        1.5e-2);
  check_small("phi stetig ueber die Materialgrenze (Konstruktion)",
              (potential_in_cell(m, s.phi, ic - 1, 2, 1.0, 0.5) -
               potential_in_cell(m, s.phi, ic, 2, 0.0, 0.5)) / cx.V,
              1.0e-15);

  // Refinement halves nothing if the answer is already exact, so check the
  // error falls instead of pinning a single number.
  QuadMesh m2;
  AxisymProblem p2;
  AxisymSolution s2;
  std::vector<char> inner2;
  solve_coax(cx, 240, 320, p2, m2, s2, inner2);
  const Real e1 = std::abs(charge_of(s, inner) / cx.V - cx.capacitance());
  const Real e2 = std::abs(charge_of(s2, inner2) / cx.V - cx.capacitance());
  std::printf("  Kapazitaetsfehler grob %.3e, fein %.3e, Verhaeltnis %.2f\n", e1, e2,
              e2 > 0 ? e1 / e2 : 0.0);
  expect("der Kapazitaetsfehler faellt bei Verfeinerung", e2 < 0.6 * e1);

  const Index ic2 = [&] {
    for (Index i = 0; i < m2.nr; ++i)
      if (m2.at(i, 0).r == cx.c) return i;
    return Index{-1};
  }();
  const Real din2 = eps0 * cx.e1 * field_in_cell(m2, s2.phi, ic2 - 1, 2, 1.0, 0.5).r;
  const Real dout2 = eps0 * cx.e2 * field_in_cell(m2, s2.phi, ic2, 2, 0.0, 0.5).r;
  const Real dn2 = std::abs(din2 - dout2) / std::abs(din2);
  std::printf("  D_n-Sprung grob %.3e, fein %.3e, Verhaeltnis %.2f\n", dn_err, dn2,
              dn2 > 0 ? dn_err / dn2 : 0.0);
  expect("der D_n-Sprung faellt bei Verfeinerung", dn2 < 0.7 * dn_err);
}

// ==========================================================================
// 2.  Axially layered dielectric -- Q1 is EXACT here
// ==========================================================================
//
// A cylinder with phi fixed on both end faces and a material interface at a
// mesh level.  The exact solution is piecewise linear in z, which Q1 represents
// exactly, so the tolerance is round-off and not a discretisation bound.
static void test_axial_layers() {
  std::printf("\n=== 2. Axial geschichtetes Dielektrikum (Q1 exakt) ===\n");
  const Real R = 2.0e-3, h = 1.0e-3, c = 4.0e-4, e1 = 3.3, e2 = 1.0, V = 250.0;
  std::vector<Real> z;
  for (Index k = 0; k <= 8; ++k) z.push_back(c * static_cast<Real>(k) / 8);
  z.back() = c;
  for (Index k = 1; k <= 12; ++k) z.push_back(c + (h - c) * static_cast<Real>(k) / 12);
  z.back() = h;
  const std::vector<Real> r = linspace(0.0, R, 7);
  QuadMesh m = rect_mesh(r, z);

  AxisymProblem p;
  p.mesh = &m;
  p.far_field = FarField::Insulated;
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), 1.0);
  for (Index j = 0; j + 1 < m.nz; ++j)
    for (Index i = 0; i + 1 < m.nr; ++i)
      p.eps_r[static_cast<std::size_t>(m.cell(i, j))] =
          (0.5 * (z[static_cast<std::size_t>(j)] + z[static_cast<std::size_t>(j + 1)]) < c) ? e1
                                                                                           : e2;
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  std::vector<char> bottom(static_cast<std::size_t>(m.n_nodes()), 0);
  for (Index i = 0; i < m.nr; ++i) {
    p.fixed[static_cast<std::size_t>(m.node(i, 0))] = 1;
    p.fixed_value[static_cast<std::size_t>(m.node(i, 0))] = V;
    bottom[static_cast<std::size_t>(m.node(i, 0))] = 1;
    p.fixed[static_cast<std::size_t>(m.node(i, m.nz - 1))] = 1;
  }
  const AxisymSolution s = solve_axisym(p);

  // D_z is constant; phi is piecewise linear with slopes in the ratio e2:e1.
  const Real D = V / (c / (eps0 * e1) + (h - c) / (eps0 * e2));
  auto exact = [&](Real zz) {
    return (zz <= c) ? V - D * zz / (eps0 * e1)
                     : V - D * (c / (eps0 * e1) + (zz - c) / (eps0 * e2));
  };
  Real worst = 0.0;
  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i)
      worst = std::max(worst,
                       std::abs(s.phi[static_cast<std::size_t>(m.node(i, j))] -
                                exact(m.at(i, j).z)) / V);
  check_small("Potential exakt reproduziert (stueckweise linear)", worst, 1.0e-12);

  const Real C = eps0 * pi * R * R / (c / e1 + (h - c) / e2);
  check("Kapazitaet der geschichteten Scheibe", charge_of(s, bottom) / V, C, 1.0e-12);

  // The interface flux jump: E_z jumps by e1/e2, D_z does not.
  const Index jc = [&] {
    for (Index j = 0; j < m.nz; ++j)
      if (m.at(0, j).z == c) return j;
    return Index{-1};
  }();
  expect("die Schichtgrenze ist exakt eine Netzzeile", jc > 0);
  const Vec2 Ea = field_in_cell(m, s.phi, 3, jc - 1, 0.5, 1.0);
  const Vec2 Eb = field_in_cell(m, s.phi, 3, jc, 0.5, 0.0);
  check("D_z stetig ueber die Schichtgrenze", eps0 * e1 * Ea.z, eps0 * e2 * Eb.z, 1.0e-12);
  check("E_z springt um e1/e2", Eb.z / Ea.z, e1 / e2, 1.0e-12);
}

// ==========================================================================
// 3.  Harmonic patch tests
// ==========================================================================
//
// phi = A + B z is harmonic and lies in the Q1 space: it must be reproduced to
// round-off, which tests the assembly and the gradients but NOT the 2*pi*r
// weight (a wrong weight leaves it in the kernel too).
//
// phi = z^2 - r^2/2 IS axisymmetric-harmonic and is not in the Q1 space, so it
// tests the weight: with the 2*pi*r dropped it is no longer a solution at all
// and the error stays O(1) instead of falling like h^2.
static void run_patch(const std::vector<Real>& r, const std::vector<Real>& z,
                      Real (*exact)(Vec2), Real* worst_out) {
  QuadMesh m = rect_mesh(r, z);
  AxisymProblem p;
  p.mesh = &m;
  p.far_field = FarField::Insulated;
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), 2.7);
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i) {
      if (i != 0 && i != m.nr - 1 && j != 0 && j != m.nz - 1) continue;
      p.fixed[static_cast<std::size_t>(m.node(i, j))] = 1;
      p.fixed_value[static_cast<std::size_t>(m.node(i, j))] = exact(m.at(i, j));
    }
  const AxisymSolution s = solve_axisym(p);
  Real worst = 0.0, scale = 0.0;
  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i) {
      const Real e = exact(m.at(i, j));
      worst = std::max(worst, std::abs(s.phi[static_cast<std::size_t>(m.node(i, j))] - e));
      scale = std::max(scale, std::abs(e));
    }
  *worst_out = worst / std::max(scale, 1e-300);
}

static Real linear_exact(Vec2 x) { return 12.0 + 3400.0 * x.z; }
static Real harmonic_exact(Vec2 x) { return x.z * x.z - 0.5 * x.r * x.r; }

static void test_patch() {
  std::printf("\n=== 3. Patch-Tests auf einem ungleichmaessigen Gitter ===\n");
  // Deliberately non-uniform, so that a mesh-dependent assembly error shows.
  std::vector<Real> r, z;
  for (Index k = 0; k <= 14; ++k) r.push_back(1.0e-3 * std::pow(static_cast<Real>(k) / 14, 1.7));
  for (Index k = 0; k <= 11; ++k)
    z.push_back(-5.0e-4 + 1.5e-3 * std::pow(static_cast<Real>(k) / 11, 0.8));
  Real w = 0.0;
  run_patch(r, z, linear_exact, &w);
  check_small("phi = A + B z exakt reproduziert", w, 1.0e-12);

  // The axisymmetric-harmonic quadratic, on two levels.
  std::vector<Real> r2, z2;
  for (Index k = 0; k <= 28; ++k) r2.push_back(1.0e-3 * std::pow(static_cast<Real>(k) / 28, 1.7));
  for (Index k = 0; k <= 22; ++k)
    z2.push_back(-5.0e-4 + 1.5e-3 * std::pow(static_cast<Real>(k) / 22, 0.8));
  Real wa = 0.0, wb = 0.0;
  run_patch(r, z, harmonic_exact, &wa);
  run_patch(r2, z2, harmonic_exact, &wb);
  std::printf("  phi = z^2 - r^2/2 : grob %.3e, fein %.3e, Verhaeltnis %.2f\n", wa, wb,
              wb > 0 ? wa / wb : 0.0);
  expect("achsensymmetrisch harmonisch: Fehler klein", wa < 5.0e-3);
  expect("achsensymmetrisch harmonisch: Fehler faellt etwa wie h^2", wa / wb > 2.5);
}

// ==========================================================================
// 4.  Band solver against a dense LU of the same system
// ==========================================================================
static void test_band_vs_dense() {
  std::printf("\n=== 4. Bandloeser gegen dichte LU-Zerlegung ===\n");
  const std::vector<Real> r = linspace(0.0, 1.0e-3, 11);
  const std::vector<Real> z = linspace(0.0, 2.0e-3, 9);
  QuadMesh m = rect_mesh(r, z);
  AxisymProblem p;
  p.mesh = &m;
  p.far_field = FarField::Asymptotic;
  p.far_field_origin = {0.0, 1.0e-3};
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), 1.0);
  for (Index j = 0; j + 1 < m.nz; ++j)
    for (Index i = 0; i + 1 < m.nr; ++i)
      if (i < 3 && j > 2 && j < 6) p.eps_r[static_cast<std::size_t>(m.cell(i, j))] = 4.0;
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  for (Index i = 0; i < 4; ++i) {
    p.fixed[static_cast<std::size_t>(m.node(i, 4))] = 1;
    p.fixed_value[static_cast<std::size_t>(m.node(i, 4))] = 500.0;
  }
  for (Index i = 0; i + 1 < m.nr; ++i) {
    p.far_edges.push_back({m.node(i, 0), m.node(i + 1, 0)});
    p.far_edges.push_back({m.node(i, m.nz - 1), m.node(i + 1, m.nz - 1)});
  }
  for (Index j = 0; j + 1 < m.nz; ++j)
    p.far_edges.push_back({m.node(m.nr - 1, j), m.node(m.nr - 1, j + 1)});

  const AxisymSolution a = solve_axisym(p, LinearSolver::Band);
  const AxisymSolution b = solve_axisym(p, LinearSolver::Dense);
  Real worst = 0.0, scale = 0.0;
  for (std::size_t k = 0; k < a.phi.size(); ++k) {
    worst = std::max(worst, std::abs(a.phi[k] - b.phi[k]));
    scale = std::max(scale, std::abs(b.phi[k]));
  }
  check_small("Bandloeser gegen dichte LU", worst / scale, 1.0e-10);
}

// ==========================================================================
// 5.  Linearity and polarity reversal
// ==========================================================================
static void test_linearity() {
  std::printf("\n=== 5. Linearitaet und Polaritaetsumkehr ===\n");
  const Coax cx;
  QuadMesh m;
  AxisymProblem p;
  AxisymSolution s;
  std::vector<char> inner;
  solve_coax(cx, 40, 50, p, m, s, inner);

  // Same problem with the applied potential scaled and reversed.
  auto resolve = [&](Real factor) {
    AxisymProblem q = p;
    q.mesh = &m;
    for (std::size_t k = 0; k < q.fixed_value.size(); ++k) q.fixed_value[k] *= factor;
    return solve_axisym(q);
  };
  const AxisymSolution s3 = resolve(3.0);
  const AxisymSolution sm = resolve(-1.0);
  Real w3 = 0.0, wm = 0.0, scale = 0.0;
  for (std::size_t k = 0; k < s.phi.size(); ++k) {
    w3 = std::max(w3, std::abs(s3.phi[k] - 3.0 * s.phi[k]));
    wm = std::max(wm, std::abs(sm.phi[k] + s.phi[k]));
    scale = std::max(scale, std::abs(s.phi[k]));
  }
  check_small("phi(3V) = 3 phi(V)", w3 / scale, 1.0e-12);
  check_small("phi(-V) = -phi(V)", wm / scale, 1.0e-12);
  check("Q(-V) = -Q(V)", charge_of(sm, inner), -charge_of(s, inner), 1.0e-12);
}

// ==========================================================================
// 6.  Failure modes that must not pass silently
// ==========================================================================
static void test_refusals() {
  std::printf("\n=== 6. Was abgelehnt werden muss ===\n");
  const std::vector<Real> r = linspace(0.0, 1.0e-3, 6);
  const std::vector<Real> z = linspace(0.0, 1.0e-3, 6);
  QuadMesh m = rect_mesh(r, z);
  AxisymProblem p;
  p.mesh = &m;
  p.far_field = FarField::Insulated;
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), 1.0);
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);

  bool threw = false;
  try {
    (void)solve_axisym(p);
  } catch (const std::exception&) {
    threw = true;
  }
  expect("ohne Dirichlet-Bedingung wird abgelehnt", threw);

  // An inactive region whose nodes are not fixed leaves empty equations.
  p.fixed[static_cast<std::size_t>(m.node(0, 0))] = 1;
  for (Index j = 1; j + 2 < m.nz; ++j)
    for (Index i = 1; i + 2 < m.nr; ++i) p.active[static_cast<std::size_t>(m.cell(i, j))] = 0;
  threw = false;
  try {
    (void)solve_axisym(p);
  } catch (const std::exception&) {
    threw = true;
  }
  expect("freie Knoten ohne aktive Nachbarzelle werden abgelehnt", threw);

  // A mesh whose rows are not level breaks exact point location.
  QuadMesh bad = rect_mesh(r, z);
  bad.nodes[static_cast<std::size_t>(bad.node(2, 2))].z += 1.0e-6;
  expect("nicht-ebene Netzzeilen werden erkannt", !bad.validate_level_rows());
}

// ==========================================================================
int main() {
  std::printf("Achsensymmetrische FEM fuer div(eps grad phi) = 0\n");
  test_coaxial_layered();
  test_axial_layers();
  test_patch();
  test_band_vs_dense();
  test_linearity();
  test_refusals();
  std::printf("\n%s (%d Fehler)\n", failures == 0 ? "ALLE TESTS BESTANDEN" : "TESTS FEHLGESCHLAGEN",
              failures);
  return failures == 0 ? 0 : 1;
}
