#include <cmath>
#include <cstdio>

#include "es/bem.hpp"
#include "es/constants.hpp"
#include "es/geometry.hpp"

using namespace es;
using constants::eps0;
using constants::pi;

static int failures = 0;

static bool check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-46s got=%-16.9g want=%-16.9g relerr=%-9.2g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
  return ok;
}

static void check_abs(const char* what, Real got, Real want, Real atol) {
  const Real err = std::abs(got - want);
  const bool ok = err <= atol;
  if (!ok) ++failures;
  std::printf("  %-46s got=%-16.9g want=%-16.9g abserr=%-9.2g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

// ---------------------------------------------------------------------------
static void test_kernel_gradient() {
  std::printf("\n=== kernel: grad G vs finite differences ===\n");
  const Vec2 src{1.3e-3, -0.4e-3};
  const Vec2 pts[] = {{2.0e-3, 1.0e-3}, {0.2e-3, 0.1e-3}, {1.35e-3, -0.39e-3}, {1e-9, 2e-3}};
  for (const Vec2& x : pts) {
    const Vec2 g = kernel_gradG(x, src);
    const Real hr = 1e-7 * std::max(x.r, 1e-6);
    const Real hz = 1e-7 * std::max(std::abs(x.z), 1e-6);
    const Real gr = (kernel_G({x.r + hr, x.z}, src) - kernel_G({x.r - hr, x.z}, src)) / (2 * hr);
    const Real gz = (kernel_G({x.r, x.z + hz}, src) - kernel_G({x.r, x.z - hz}, src)) / (2 * hz);
    // On the axis dG/dr vanishes identically, so a purely relative tolerance is
    // meaningless there; compare against the magnitude of the full gradient.
    const Real scale = std::max(std::abs(gr), std::abs(gz));
    char buf[96];
    std::snprintf(buf, sizeof buf, "dG/dr at (%.2e,%.2e)", x.r, x.z);
    check_abs(buf, g.r, gr, 2e-5 * scale);
    std::snprintf(buf, sizeof buf, "dG/dz at (%.2e,%.2e)", x.r, x.z);
    check_abs(buf, g.z, gz, 2e-5 * scale);
  }
}

// ---------------------------------------------------------------------------
static void test_sphere() {
  std::printf("\n=== isolated conducting sphere (R = 1 mm, V = 1000 V) ===\n");
  const Real R = 1e-3, V = 1000.0;

  for (int n : {40, 160, 640}) {
    BemSolver s(make_sphere(R, V, n));
    s.solve({V, 0.0, 0.0});

    const Real C = s.capacitance();
    const Real Cex = 4.0 * pi * eps0 * R;

    // sigma must be uniform: E_n = V/R everywhere
    Real emin = 1e300, emax = -1e300;
    for (Index i = 0; i < s.size(); ++i) {
      const Real e = s.En(i);
      emin = std::min(emin, e);
      emax = std::max(emax, e);
    }

    std::printf(" N = %d\n", static_cast<int>(s.size()));
    char buf[96];
    std::snprintf(buf, sizeof buf, "capacitance (4 pi eps0 R)");
    check(buf, C, Cex, n >= 160 ? 2e-4 : 3e-3);
    check("E_n min", emin, V / R, n >= 160 ? 2e-3 : 2e-2);
    check("E_n max", emax, V / R, n >= 160 ? 2e-3 : 2e-2);

    if (n == 640) {
      // Exterior potential V R / d, and the field V R / d^2.
      for (Real d : {1.5e-3, 3e-3, 1e-2}) {
        const Vec2 x{0.0, d};
        std::snprintf(buf, sizeof buf, "V at z = %.1f mm", d * 1e3);
        check(buf, s.potential_at(x), V * R / d, 1e-4);
        const Vec2 E = s.field_at(x);
        std::snprintf(buf, sizeof buf, "E_z at z = %.1f mm", d * 1e3);
        check(buf, E.z, V * R / (d * d), 1e-4);
        std::snprintf(buf, sizeof buf, "E_r at z = %.1f mm (must vanish)", d * 1e3);
        check(buf, E.r + 1.0, 1.0, 1e-6);
      }
      // Off-axis point, 45 degrees
      const Real d = 4e-3;
      const Vec2 x{d / std::sqrt(2.0), d / std::sqrt(2.0)};
      check("V off-axis", s.potential_at(x), V * R / d, 1e-4);
      const Vec2 E = s.field_at(x);
      check("|E| off-axis", std::sqrt(E.r * E.r + E.z * E.z), V * R / (d * d), 1e-4);
    }
  }
}

// ---------------------------------------------------------------------------
static void test_spheroid() {
  // The closest analytic stand-in for a sharpened needle: high aspect ratio,
  // strong tip-field enhancement, and a known surface-field distribution.
  std::printf("\n=== prolate spheroid a = 200 um, b = 10 um (aspect 20), V = 1000 V ===\n");
  const Real a = 200e-6, b = 10e-6, V = 1000.0;
  const Real c = std::sqrt(a * a - b * b);
  const Real Q0 = 0.5 * std::log((a + c) / (a - c));

  BemSolver s(make_prolate_spheroid(a, b, V, 600));
  s.solve({V, 0.0, 0.0});
  std::printf(" N = %d,  field enhancement E_tip/(V/a) = %.2f\n", static_cast<int>(s.size()),
              spheroid_tip_field(a, b, V) / (V / a));

  check("capacitance", s.capacitance(), spheroid_capacitance(a, b), 3e-3);

  // Element-wise surface field:  E_n = V / ( Q0(xi0) c sqrt((xi0^2-1)(xi0^2-eta^2)) )
  const Real xi0 = a / c;
  Real worst = 0.0;
  Real tip_num = 0.0, tip_ana = 0.0;
  for (Index i = 0; i < s.size(); ++i) {
    const Element& e = s.mesh().elems[static_cast<std::size_t>(i)];
    const Real eta = e.mid.z / a;
    const Real ex = V / (Q0 * c * std::sqrt((xi0 * xi0 - 1.0) * (xi0 * xi0 - eta * eta)));
    const Real rel = std::abs(s.En(i) - ex) / ex;
    worst = std::max(worst, rel);
    if (i == s.size() - 1) { tip_num = s.En(i); tip_ana = ex; }
  }
  std::printf("  worst element-wise E_n error over the whole surface: %.3g\n", worst);
  if (worst > 5e-3) ++failures;
  std::printf("  %s\n", worst <= 5e-3 ? "OK" : "FAIL");
  check("E_n on the last (apex) element", tip_num, tip_ana, 5e-3);
}

// ---------------------------------------------------------------------------
static void test_bc_residual() {
  // The strongest self-consistency check: re-evaluate the potential integral at
  // every collocation point and confirm it reproduces the Dirichlet data.
  std::printf("\n=== boundary-condition residual, needle + extractor ===\n");
  NeedleParams np;
  np.tip_radius = 5e-6;
  np.half_angle = 15.0 * pi / 180.0;
  np.shank_radius = 1.5e-4;
  np.length = 1.0e-3;
  np.z_tip = 0.0;
  ExtractorParams ep;
  ep.aperture_radius = 2.5e-4;
  ep.outer_radius = 2.0e-3;
  ep.thickness = 1.0e-4;
  ep.z_plate = 3.0e-4;
  ep.potential = -1500.0;

  Mesh m = merge({make_needle(np), make_extractor(ep)});
  BemSolver s(m);
  const Real Ve = 0.0, Vx = -1500.0;
  s.solve({Ve, Vx, 0.0});
  std::printf(" N = %d elements\n", static_cast<int>(s.size()));

  Real worst = 0.0;
  const Real scale = std::max(std::abs(Ve), std::abs(Vx));
  for (Index i = 0; i < s.size(); ++i) {
    const Element& e = s.mesh().elems[static_cast<std::size_t>(i)];
    const Real v = s.potential_at(e.mid);
    worst = std::max(worst, std::abs(v - e.potential) / scale);
  }
  std::printf("  worst |V_computed - V_bc| / V_scale = %.3g\n", worst);
  if (worst > 1e-6) ++failures;
  std::printf("  %s\n", worst <= 1e-6 ? "OK" : "FAIL");

  // Charge neutrality: an isolated two-electrode system in free space carries
  // net charge, but E far away must fall off as a monopole with exactly that
  // charge.  Check Gauss' law on a large sphere.
  s.solve({Ve, Vx, 0.0});
  const Real Q = s.total_charge();
  const Real Rg = 0.5;  // 0.5 m -- far outside the 2 mm structure
  const Real Vfar = s.potential_at({0.0, Rg});
  check("far-field potential = Q/(4 pi eps0 R)", Vfar, Q / (4 * pi * eps0 * Rg), 5e-3);

  std::printf("  peak |E_n| on the needle: %.4g V/m at V_ext = %.0f V\n",
              s.peak_emitter_field(), Vx);
}

// ---------------------------------------------------------------------------
static void test_linearity() {
  std::printf("\n=== superposition (basis solve) ===\n");
  BemSolver s(make_sphere(1e-3, 0.0, 120));
  s.solve_basis();
  const std::vector<Real> s1 = s.sigma_for({1000.0, 0.0, 0.0});
  const std::vector<Real> s2 = s.sigma_for({2500.0, 0.0, 0.0});
  Real worst = 0.0;
  for (std::size_t i = 0; i < s1.size(); ++i)
    worst = std::max(worst, std::abs(s2[i] - 2.5 * s1[i]) / std::abs(2.5 * s1[i]));
  std::printf("  worst deviation from exact linearity: %.3g\n", worst);
  if (worst > 1e-12) ++failures;
  std::printf("  %s\n", worst <= 1e-12 ? "OK" : "FAIL");
}

int main() {
  test_kernel_gradient();
  test_sphere();
  test_spheroid();
  test_bc_residual();
  test_linearity();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
