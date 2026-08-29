#include <cmath>
#include <cstdio>

#include "es/constants.hpp"
#include "es/emission.hpp"
#include "es/fluid.hpp"
#include <vector>
#include <algorithm>
#include "es/meniscus.hpp"
#include "es/status.hpp"

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
    expect("converged", s.ok());
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
// 2.  Static fold: continuation in apex height, turning point of U(h).
//     This is NOT an emission onset -- see docs/02_model_specification.md 2.4.
// ---------------------------------------------------------------------------
static Mesh fold_electrodes() {
  OpenCapillaryParams cp;
  cp.r_bore = 1.0e-5;
  cp.r_outer = 2.0e-5;
  cp.shank_length = 1.0e-3;
  cp.z_rim = 0.0;
  cp.h_rim = 1.0e-5 / 14.0;
  ExtractorParams ep;
  ep.aperture_radius = 2.0e-4;
  ep.outer_radius = 3.0e-3;
  ep.thickness = 1.0e-4;
  ep.z_plate = 5.0e-4;
  ep.h_edge = 1.0e-5;
  return merge({make_capillary_open(cp), make_extractor(ep)});
}

static void test_static_fold() {
  std::printf("\n=== static fold of the branch, capillary opposite an extractor ===\n");
  const Fluid f = fluid_by_name("EMI-BF4");
  const Real r_c = 1.0e-5;
  const Real gap = 5.0e-4;

  MeniscusParams mp;
  mp.r_contact = r_c;
  mp.z_contact = 0.0;
  mp.gamma = f.gamma;
  mp.delta_p = 0.0;
  mp.n_nodes = 81;
  mp.max_outer = 60;
  mp.relax = 0.5;
  mp.tol = 3e-4;

  MeniscusSolver solver(fold_electrodes(), mp);
  const std::vector<MeniscusSolution> branch = solver.continuation(0.15 * r_c, 2.2 * r_c, 22);

  std::printf("   h/r_c      U [V]     E_apex [V/m]   R_apex/r_c   half-angle   status\n");
  for (const MeniscusSolution& s : branch)
    std::printf("  %6.3f  %10.1f  %13.4g  %10.4f  %8.2f deg  %s (%d it)\n",
                s.shape.height / r_c, s.voltage, s.apex_field, s.shape.apex_radius / r_c,
                s.shape.half_angle * 180.0 / pi, to_string(s.status), s.iterations);

  const MeniscusSolver::StaticFold fold = MeniscusSolver::find_static_fold(branch);
  expect("an interior turning point was demonstrated", fold.found());
  if (!fold.found()) {
    std::printf("  %s\n", explain(fold.status));
    return;
  }

  const Real V_lit = literature_onset_voltage_smith(r_c, gap, f.gamma);
  std::printf("\n  static fold voltage                 : %10.1f V\n", fold.voltage);
  std::printf("  Smith (1986) onset closed form      : %10.1f V   (ratio %.3f)\n", V_lit,
              fold.voltage / V_lit);
  std::printf("  apex height at the fold             : %10.4f r_c\n", fold.height / r_c);
  std::printf("  apex field at the fold              : %10.4g V/m\n", fold.apex_field);
  std::printf("  cone half-angle at mid-arc          : %10.2f deg  (Taylor: 49.29)\n",
              fold.half_angle * 180.0 / pi);

  expect("fold within a factor 2 of the closed form",
         fold.voltage > 0.5 * V_lit && fold.voltage < 2.0 * V_lit);
  expect("turning point is interior to the branch",
         fold.height > 1.05 * branch.front().shape.height &&
             fold.height < 0.95 * branch.back().shape.height);
  {
    bool monotone = true;
    Real last_angle = 1e9, min_angle = 1e9;
    for (const MeniscusSolution& s : branch) {
      if (!s.ok()) continue;
      if (s.shape.half_angle > last_angle + 1e-9) monotone = false;
      last_angle = s.shape.half_angle;
      min_angle = std::min(min_angle, s.shape.half_angle);
    }
    expect("cone half-angle decreases monotonically along the branch", monotone);
    std::printf("  smallest half-angle on the branch   : %10.2f deg\n", min_angle * 180.0 / pi);
    expect("branch approaches Taylor's 49.3 deg within 3 deg",
           std::abs(min_angle - constants::taylor_angle) < 3.0 * pi / 180.0);
  }
}

// ---------------------------------------------------------------------------
// 3.  Mesh convergence.  Finding 9: the prototype checked only the fold
//     voltage.  Every quantity that is still reported must be checked.
// ---------------------------------------------------------------------------
static void test_mesh_convergence() {
  std::printf("\n=== mesh convergence of every reported quantity ===\n");
  const Fluid f = fluid_by_name("EMI-BF4");
  const Real r_c = 1.0e-5;

  struct Row { int n; Real V, E, R, I, A; };
  std::vector<Row> rows;

  for (int n : {61, 81, 121, 161}) {
    MeniscusParams mp;
    mp.r_contact = r_c;
    mp.gamma = f.gamma;
    mp.n_nodes = n;
    mp.max_outer = 60;
    mp.relax = 0.5;
    mp.tol = 3e-4;
    MeniscusSolver s(fold_electrodes(), mp);
    const std::vector<MeniscusSolution> br = s.continuation(0.15 * r_c, 2.2 * r_c, 22);
    const auto fold = MeniscusSolver::find_static_fold(br);
    if (!fold.found()) { ++failures; std::printf("  n = %d: kein Umkehrpunkt\n", n); continue; }
    const MeniscusSolution& at = br[fold.index];
    s.realize(at);
    const IonEmission ie = integrate_ion_emission(s.bem(), f, 298.15, false);
    rows.push_back({n, fold.voltage, at.apex_field, at.shape.apex_radius, ie.current,
                    ie.effective_area});
  }

  std::printf("  %6s %12s %14s %13s %13s %13s\n", "nodes", "U_fold [V]", "E_apex [V/m]",
              "R_apex [m]", "I_ion [A]", "A_eff [m^2]");
  for (const Row& r : rows)
    std::printf("  %6d %12.3f %14.6g %13.6g %13.6g %13.6g\n", r.n, r.V, r.E, r.R, r.I, r.A);
  if (rows.size() < 2) { expect("at least two refinement levels", false); return; }

  const Row& a = rows.front();
  const Row& b = rows.back();
  auto rel = [](Real x, Real y) { return std::abs(x - y) / std::max(std::abs(y), 1e-300); };
  check("fold voltage, coarsest vs finest", a.V, b.V, 2e-3);
  check("apex field, coarsest vs finest", a.E, b.E, 5e-3);
  check("apex radius, coarsest vs finest", a.R, b.R, 1e-2);
  check("ion current estimate, coarsest vs finest", a.I, b.I, 2e-2);
  check("effective emission area, coarsest vs finest", a.A, b.A, 5e-2);

  // The effective area must converge, not wander: successive differences shrink.
  if (rows.size() >= 3) {
    bool shrinking = true;
    for (std::size_t i = 2; i + 1 < rows.size(); ++i) {
      const Real d0 = rel(rows[i - 1].A, rows[i - 2].A);
      const Real d1 = rel(rows[i].A, rows[i - 1].A);
      if (d1 > d0 * 1.5 + 1e-4) shrinking = false;
    }
    expect("A_eff differences shrink under refinement (no quantisation)", shrinking);
  }
}

int main() {
  test_zero_field_spherical_cap();
  test_static_fold();
  test_mesh_convergence();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
