#include <cmath>
#include <cstdio>

#include "es/constants.hpp"
#include "es/emission.hpp"
#include "es/fluid.hpp"
#include "es/meniscus.hpp"

using namespace es;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-46s got=%-15.7g want=%-15.7g relerr=%-9.2g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-46s %s\n", what, ok ? "OK" : "FAIL");
}

// ---------------------------------------------------------------------------
// 1.  Field-free limit.  With no voltage, the Young-Laplace equation reduces to
//     a surface of constant mean curvature -- a spherical cap of radius
//     R = 2 gamma / dp.  This exercises the ODE, the apex series start and the
//     resampling completely independently of the electrostatics.
// ---------------------------------------------------------------------------
static void test_zero_field_spherical_cap() {
  std::printf("\n=== field-free limit: spherical cap of radius 2 gamma / dp ===\n");
  const Real gamma = 0.0452;
  const Real r_c = 1.0e-5;

  // Kept below a hemisphere: past h = r_c the cap overhangs, which the
  // r-monotone shape solver deliberately does not track (see meniscus.hpp).
  for (Real h_over_rc : {0.3, 0.6, 0.85}) {
    const Real h = h_over_rc * r_c;
    const Real R = (h * h + r_c * r_c) / (2.0 * h);  // sphere through apex and rim
    const Real dp = 2.0 * gamma / R;

    OpenCapillaryParams cp;
    cp.r_bore = r_c;
    cp.r_outer = 2.0 * r_c;
    cp.shank_length = 40.0 * r_c;
    cp.z_rim = 0.0;

    MeniscusParams mp;
    mp.r_contact = r_c;
    mp.z_contact = 0.0;
    mp.gamma = gamma;
    mp.delta_p = dp;
    mp.n_nodes = 121;
    mp.tol = 1e-5;

    MeniscusSolver solver(make_capillary_open(cp), mp);
    MeniscusSolution s = solver.solve_at_height(h);

    char buf[96];
    std::printf(" h/r_c = %.2f  (R_sphere = %.4g m)\n", h_over_rc, R);
    expect("converged", s.converged);
    std::snprintf(buf, sizeof buf, "voltage required (must be 0)");
    check(buf, s.voltage + 1.0, 1.0, 1e-9);
    std::snprintf(buf, sizeof buf, "apex radius of curvature");
    check(buf, s.shape.apex_radius, R, 1e-6);
    // Every node must lie on the sphere centred at (0, h - R).
    Real worst = 0.0;
    for (const Vec2& p : s.shape.nodes) {
      const Real d = std::sqrt(p.r * p.r + (p.z - (h - R)) * (p.z - (h - R)));
      worst = std::max(worst, std::abs(d - R) / R);
    }
    std::printf("  %-46s %.3g %s\n", "worst node deviation from the sphere", worst,
                worst < 2e-5 ? "OK" : "FAIL");
    if (worst >= 2e-5) ++failures;
  }
}

// ---------------------------------------------------------------------------
// 2.  Onset of emission: continuation in apex height, onset = max of U(h).
// ---------------------------------------------------------------------------
static void test_onset() {
  std::printf("\n=== onset of emission, capillary opposite an extractor ===\n");
  const Fluid f = fluid_by_name("EMI-BF4");
  const Real r_c = 1.0e-5;   // 10 um bore
  const Real gap = 5.0e-4;   // 500 um to the extractor

  OpenCapillaryParams cp;
  cp.r_bore = r_c;
  cp.r_outer = 2.0e-5;
  cp.shank_length = 1.0e-3;
  cp.z_rim = 0.0;
  cp.h_rim = r_c / 14.0;

  ExtractorParams ep;
  ep.aperture_radius = 2.0e-4;
  ep.outer_radius = 3.0e-3;
  ep.thickness = 1.0e-4;
  ep.z_plate = gap;
  ep.h_edge = 1.0e-5;

  MeniscusParams mp;
  mp.r_contact = r_c;
  mp.z_contact = 0.0;
  mp.gamma = f.gamma;
  mp.delta_p = 0.0;
  mp.n_nodes = 81;
  mp.max_outer = 60;
  mp.relax = 0.5;
  mp.tol = 3e-4;

  MeniscusSolver solver(merge({make_capillary_open(cp), make_extractor(ep)}), mp);
  const std::vector<MeniscusSolution> branch = solver.continuation(0.15 * r_c, 2.2 * r_c, 22);
  MeniscusSolver::write_branch_csv("branch_test.csv", branch);

  std::printf("   h/r_c      U [V]     E_apex [V/m]   R_apex/r_c   half-angle   conv\n");
  for (const MeniscusSolution& s : branch)
    std::printf("  %6.3f  %10.1f  %13.4g  %10.4f  %8.2f deg  %s (%d it)\n",
                s.shape.height / r_c, s.voltage, s.apex_field, s.shape.apex_radius / r_c,
                s.shape.half_angle * 180.0 / pi, s.converged ? "yes" : "NO ", s.iterations);

  const MeniscusSolver::Onset on = MeniscusSolver::find_onset(branch);
  expect("a turning point was found", on.found);
  if (!on.found) return;

  const Real V_taylor = onset_voltage_taylor(r_c, gap, f.gamma);
  std::printf("\n  onset voltage (BEM + Young-Laplace) : %10.1f V\n", on.voltage);
  std::printf("  Taylor/Smith closed form            : %10.1f V\n", V_taylor);
  std::printf("  ratio                               : %10.3f\n", on.voltage / V_taylor);
  std::printf("  apex height at onset                : %10.4f r_c\n", on.height / r_c);
  std::printf("  apex field at onset                 : %10.4g V/m\n", on.apex_field);
  std::printf("  cone half-angle at mid-arc          : %10.2f deg  (Taylor: 49.29)\n",
              on.half_angle * 180.0 / pi);

  // The closed form ignores wall thickness, aperture geometry and the finite
  // extractor, so agreement is expected only to within a factor of order one.
  expect("onset within a factor 2 of the closed form",
         on.voltage > 0.5 * V_taylor && on.voltage < 2.0 * V_taylor);
  // The turning point must be an interior maximum, not the end of the sweep.
  expect("turning point is interior to the branch",
         on.height > 1.05 * branch.front().shape.height &&
             on.height < 0.95 * branch.back().shape.height);
  // Physical signature of the branch: as the meniscus elongates it must sharpen
  // monotonically, and the far end of the branch must approach Taylor's
  // equilibrium cone angle.  Nothing in the solver knows about 49.3 deg -- it
  // has to come out of the coupled Young-Laplace / BEM problem.
  {
    bool monotone = true;
    Real last_angle = 1e9;
    Real min_angle = 1e9;
    for (const MeniscusSolution& s : branch) {
      if (!s.converged) continue;
      if (s.shape.half_angle > last_angle + 1e-9) monotone = false;
      last_angle = s.shape.half_angle;
      min_angle = std::min(min_angle, s.shape.half_angle);
    }
    expect("cone half-angle decreases monotonically along the branch", monotone);
    std::printf("  smallest half-angle on the branch   : %10.2f deg\n",
                min_angle * 180.0 / pi);
    expect("branch approaches Taylor's 49.3 deg within 3 deg",
           std::abs(min_angle - constants::taylor_angle) < 3.0 * pi / 180.0);
  }

  // Mesh convergence of the onset voltage.
  mp.n_nodes = 141;
  MeniscusSolver s2(merge({make_capillary_open(cp), make_extractor(ep)}), mp);
  const MeniscusSolver::Onset on2 =
      MeniscusSolver::find_onset(s2.continuation(0.15 * r_c, 2.2 * r_c, 22));
  std::printf("\n  onset with 141 free-surface nodes   : %10.1f V\n", on2.voltage);
  check("onset voltage, 81 vs 141 nodes", on2.voltage, on.voltage, 0.02);
}

int main() {
  test_zero_field_spherical_cap();
  test_onset();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
