// P3a -- static capillary meniscus without an electric field.
//
// Every tolerance in this file is declared at the top, BEFORE any number is
// looked at, and each one carries the reason for its size.  None of them was
// adjusted after a result was seen.
//
// The reference is the closed-form spherical cap, which is what
// gamma*(dpsi/ds + sin psi / r) = delta_p must produce when nothing else acts
// on the surface.  The solver never evaluates that closed form: it integrates
// the Young-Laplace ODE and finds the arclength at which the meridian reaches
// the pinning radius, so agreement is a statement about the solver and not a
// tautology.

#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "es/capillary.hpp"
#include "es/constants.hpp"
#include "es/device_geometry.hpp"
#include "es/liquid.hpp"

using namespace es;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-56s got=%-14.7g want=%-14.7g rel=%-9.2e %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

static void check_abs(const char* what, Real got, Real want, Real atol) {
  const Real err = std::abs(got - want);
  const bool ok = err <= atol;
  if (!ok) ++failures;
  std::printf("  %-56s got=%-14.7g want=%-14.7g abs=%-9.2e %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-56s %s\n", what, ok ? "OK" : "FAIL");
}

// ===========================================================================
// Tolerances -- fixed before the first run
// ===========================================================================
namespace tol {

/// The flat surface is not an approximation: with delta_p = 0 the right-hand
/// side and psi are identically zero in floating point, so every z must be
/// exactly the contact plane.  Zero, not "small".
constexpr Real kFlatExact = 0.0;

/// The arclength search terminates at |r_end - a| <= 8 eps a, so the pinned
/// edge must be hit to a few tens of eps.  1e-14 is about 45 eps.
constexpr Real kPinning = 1.0e-14;

/// Requested accuracy handed to the solver in the analytic comparisons.
constexpr Real kRequestedAccuracy = 1.0e-12;

/// Agreement with the closed form.  The refinement is asked for 1e-12 and RK4
/// is fourth order, so 1e-10 leaves two decades of head-room and is still far
/// below anything a wrong equation could survive.
constexpr Real kAnalytic = 1.0e-10;

/// Young-Laplace residual from the node coordinates, normalised by gamma/a.
/// The geometric curvature estimate is second order in the node spacing, so its
/// floor is set by the automatically chosen resolution (>= 128 intervals ->
/// ~1e-5); 1e-4 is a deliberately loose bound on an independent check whose
/// convergence order is verified separately below.
constexpr Real kResidual = 1.0e-4;

/// Observed convergence orders.  RK4 gives four, the geometric residual two;
/// the bounds are set well below the theoretical values so that a
/// slowly-converging scheme fails but round-off noise does not.
constexpr Real kMinOrderProfile = 3.5;
constexpr Real kMinOrderResidual = 1.5;

/// Geometric similarity: dimensionless quantities of a scaled problem must
/// agree to round-off accumulated over the integration, not to a physical
/// tolerance.
constexpr Real kScaling = 1.0e-11;

}  // namespace tol

// ---------------------------------------------------------------------------

static LiquidProperties test_liquid() { return unit_liquid(); }

static CapillaryMeniscus solve_pi(Real a, Real Pi, const LiquidProperties& L,
                                  int forced = 0) {
  CapillaryRequest q;
  q.delta_p_exit = capillary::pressure_from_pi(Pi, a, L.gamma);
  q.target_relative_accuracy = tol::kRequestedAccuracy;
  q.forced_intervals = forced;
  return solve_capillary_meniscus(a, 0.0, L, q);
}

// ===========================================================================
// 1.  delta_p_exit = 0 gives an exactly flat surface
// ===========================================================================
static void test_flat() {
  std::printf("\n=== 1. delta_p_exit = 0 -> exakt ebene Oberflaeche ===\n");
  const LiquidProperties L = test_liquid();
  const Real a = 5.0e-6;
  const CapillaryMeniscus m = solve_pi(a, 0.0, L);
  expect("Status Solved", is_usable(m.status));

  Real worst_z = 0.0, worst_psi = 0.0;
  for (std::size_t i = 0; i < m.nodes.size(); ++i) {
    worst_z = std::max(worst_z, std::abs(m.nodes[i].z - m.contact_z));
    worst_psi = std::max(worst_psi, std::abs(m.psi[i]));
  }
  check_abs("groesste Abweichung von der Ebene z = 0", worst_z, 0.0, tol::kFlatExact);
  check_abs("groesster Tangentenwinkel psi", worst_psi, 0.0, tol::kFlatExact);
  check_abs("Apexhoehe", m.apex_height, 0.0, tol::kFlatExact);
  check("Bogenlaenge = a", m.arclength, a, tol::kAnalytic);
  check("Rotationsflaeche = pi a^2", m.revolved_area, pi * a * a, tol::kAnalytic);
  check_abs("Rotationsvolumen = 0", m.revolved_volume, 0.0, 1e-30);

  // The same statement against the closed form, which is the flat disc.
  const SphericalCap c = spherical_cap(a, 0.0, L.gamma);
  check_abs("geschlossene Form: Apexhoehe = 0", c.apex_height, 0.0, tol::kFlatExact);
  check("geschlossene Form: Flaeche = pi a^2", c.revolved_area, pi * a * a, 1e-15);
}

// ===========================================================================
// 2.  The pinned exit edge is hit exactly
// ===========================================================================
static void test_pinning() {
  std::printf("\n=== 2. Gepinnte Austrittskante wird exakt getroffen ===\n");
  const LiquidProperties L = test_liquid();
  const Real a = 5.0e-6;
  for (Real Pi : {-1.7, -0.9, -0.2, 0.0, 0.2, 0.9, 1.7}) {
    const CapillaryMeniscus m = solve_pi(a, Pi, L);
    char what[96];
    std::snprintf(what, sizeof what, "Pi = %+5.2f: |r_Kontakt - a|/a", Pi);
    check_abs(what, std::abs(m.contact().r - a) / a, 0.0, tol::kPinning);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: z_Kontakt - z_Kante", Pi);
    check_abs(what, m.contact().z - m.contact_z, 0.0, 0.0);
  }
}

// ===========================================================================
// 3.  Symmetry and the regular curvature limit on the axis
// ===========================================================================
static void test_axis_regularity() {
  std::printf("\n=== 3. Symmetrie und regulaerer Kruemmungsgrenzwert bei r = 0 ===\n");
  const LiquidProperties L = test_liquid();
  const Real a = 5.0e-6;
  for (Real Pi : {-1.5, -0.5, 0.5, 1.5}) {
    const CapillaryMeniscus m = solve_pi(a, Pi, L);
    const Real dp = m.delta_p_exit;
    char what[96];

    std::snprintf(what, sizeof what, "Pi = %+5.2f: r(Apex) = 0 exakt", Pi);
    check_abs(what, m.nodes.front().r, 0.0, 0.0);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: psi(Apex) = 0 exakt", Pi);
    check_abs(what, m.psi.front(), 0.0, 0.0);

    // dpsi/ds at the apex must be kappa/2 -- the analytic limit of sin(psi)/r.
    // Read off the numerical solution, not from the formula the solver uses.
    const Real ds = norm(m.nodes[1] - m.nodes[0]);
    const Real dpsi_ds = m.psi[1] / ds;
    std::snprintf(what, sizeof what, "Pi = %+5.2f: dpsi/ds(Apex) = kappa/2", Pi);
    check(what, dpsi_ds, 0.5 * dp / m.gamma, 1e-6);

    // Both principal curvatures are equal there, so the total is kappa.
    std::snprintf(what, sizeof what, "Pi = %+5.2f: 2*dpsi/ds(Apex) = dp/gamma", Pi);
    check(what, 2.0 * dpsi_ds, dp / m.gamma, 1e-6);

    // Nothing on the axis is NaN or infinite -- the 0/0 term is handled.
    bool finite = true;
    for (const Vec2& p : m.nodes) finite = finite && std::isfinite(p.r) && std::isfinite(p.z);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: alle Knoten endlich", Pi);
    expect(what, finite);
  }
}

// ===========================================================================
// 4./5.  Profile, apex height, area and volume against the closed form
// ===========================================================================
static void test_against_spherical_cap() {
  std::printf("\n=== 4./5. Profil, Apexhoehe, Flaeche, Volumen gegen die Kugelkappe ===\n");
  const LiquidProperties L = test_liquid();
  const Real a = 5.0e-6;
  for (Real Pi : {-1.98, -1.5, -1.0, -0.5, -0.05, 0.05, 0.5, 1.0, 1.5, 1.98}) {
    const CapillaryMeniscus m = solve_pi(a, Pi, L);
    const SphericalCap c = spherical_cap(a, m.delta_p_exit, L.gamma);
    char what[96];

    std::snprintf(what, sizeof what, "Pi = %+5.2f: Status Solved", Pi);
    expect(what, is_usable(m.status));

    std::snprintf(what, sizeof what, "Pi = %+5.2f: Profil (Normalabstand)/a", Pi);
    check_abs(what, profile_error_against_cap(m), 0.0, tol::kAnalytic);

    std::snprintf(what, sizeof what, "Pi = %+5.2f: Profil (|dz|)/a", Pi);
    check_abs(what, profile_z_error_against_cap(m), 0.0, tol::kAnalytic);

    std::snprintf(what, sizeof what, "Pi = %+5.2f: Apexhoehe", Pi);
    check(what, m.apex_height, c.apex_height, tol::kAnalytic);

    std::snprintf(what, sizeof what, "Pi = %+5.2f: Bogenlaenge", Pi);
    check(what, m.arclength, c.arclength, tol::kAnalytic);

    std::snprintf(what, sizeof what, "Pi = %+5.2f: Rotationsflaeche", Pi);
    check(what, m.revolved_area, c.revolved_area, tol::kAnalytic);

    std::snprintf(what, sizeof what, "Pi = %+5.2f: Rotationsvolumen", Pi);
    check(what, m.revolved_volume, c.revolved_volume, tol::kAnalytic);

    std::snprintf(what, sizeof what, "Pi = %+5.2f: psi(Kontaktlinie)", Pi);
    check(what, m.contact_tangent_angle, c.contact_tangent_angle, tol::kAnalytic);

    // The sign convention, stated as a test and not only in a comment.
    std::snprintf(what, sizeof what, "Pi = %+5.2f: Vorzeichen der Woelbung", Pi);
    expect(what, (Pi > 0.0 && m.apex_height > 0.0 && m.revolved_volume > 0.0) ||
                     (Pi < 0.0 && m.apex_height < 0.0 && m.revolved_volume < 0.0));
  }

  // The hemisphere, where the closed form is h = a and A = 2 pi a^2.
  const SphericalCap hemi = spherical_cap(a, capillary::pressure_from_pi(2.0, a, L.gamma),
                                          L.gamma);
  check("Halbkugel: h = a", hemi.apex_height, a, 1e-15);
  check("Halbkugel: Flaeche = 2 pi a^2", hemi.revolved_area, 2.0 * pi * a * a, 1e-15);
  check("Halbkugel: Volumen = 2/3 pi a^3", hemi.revolved_volume, 2.0 / 3.0 * pi * a * a * a,
        1e-15);
  check("Halbkugel: psi(Kontakt) = 90 deg", hemi.contact_tangent_angle, 0.5 * pi, 1e-15);
}

// ===========================================================================
// 6.  Young-Laplace residual along the whole surface
// ===========================================================================
static void test_residual() {
  std::printf("\n=== 6. Young-Laplace-Residuum entlang der gesamten Oberflaeche ===\n");
  const LiquidProperties L = test_liquid();
  const Real a = 5.0e-6;
  for (Real Pi : {-1.9, -1.0, -0.3, 0.3, 1.0, 1.9}) {
    const CapillaryMeniscus m = solve_pi(a, Pi, L);
    const ResidualProfile R = young_laplace_residual(m);
    char what[96];
    std::snprintf(what, sizeof what, "Pi = %+5.2f: max|Residuum| (%zu Knoten)", Pi,
                  R.residual.size());
    check_abs(what, R.max_abs, 0.0, tol::kResidual);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: Residuum an allen Knoten ausgewertet", Pi);
    expect(what, R.residual.size() == m.nodes.size());
  }
}

// ===========================================================================
// 7.  Mesh convergence over four resolutions
// ===========================================================================
static Real order_of(Real e_coarse, Real e_fine) {
  if (!(e_coarse > 0.0) || !(e_fine > 0.0)) return 0.0;
  return std::log(e_coarse / e_fine) / std::log(2.0);
}

static void test_mesh_convergence() {
  std::printf("\n=== 7. Netzkonvergenz ueber vier Aufloesungen ===\n");
  const LiquidProperties L = test_liquid();
  const Real a = 5.0e-6;
  const int levels[5] = {16, 32, 64, 128, 256};

  for (Real Pi : {-1.2, 0.8}) {
    std::printf("  Pi = %+5.2f\n", Pi);
    std::vector<Real> e_profile, e_res, e_h, e_A, e_V;
    for (int n : levels) {
      const CapillaryMeniscus m = solve_pi(a, Pi, L, n);
      const SphericalCap c = spherical_cap(a, m.delta_p_exit, L.gamma);
      const ResidualProfile R = young_laplace_residual(m);
      e_profile.push_back(profile_error_against_cap(m));
      e_res.push_back(R.max_abs);
      e_h.push_back(std::abs(m.apex_height - c.apex_height) / a);
      e_A.push_back(std::abs(m.revolved_area - c.revolved_area) / (a * a));
      e_V.push_back(std::abs(m.revolved_volume - c.revolved_volume) / (a * a * a));
      std::printf("    n=%4d  Profil %.3e  h %.3e  A %.3e  V %.3e  Residuum %.3e\n", n,
                  e_profile.back(), e_h.back(), e_A.back(), e_V.back(), e_res.back());
    }
    for (std::size_t k = 0; k + 1 < e_profile.size(); ++k) {
      char what[96];
      std::snprintf(what, sizeof what, "    Ordnung Profil %d->%d", levels[k], levels[k + 1]);
      const Real p = order_of(e_profile[k], e_profile[k + 1]);
      if (p < tol::kMinOrderProfile) ++failures;
      std::printf("  %-56s %.2f %s\n", what, p, p >= tol::kMinOrderProfile ? "OK" : "FAIL");
    }
    for (std::size_t k = 0; k + 1 < e_res.size(); ++k) {
      char what[96];
      std::snprintf(what, sizeof what, "    Ordnung Residuum %d->%d", levels[k], levels[k + 1]);
      const Real p = order_of(e_res[k], e_res[k + 1]);
      if (p < tol::kMinOrderResidual) ++failures;
      std::printf("  %-56s %.2f %s\n", what, p, p >= tol::kMinOrderResidual ? "OK" : "FAIL");
    }
    expect("    Apexhoehenfehler faellt monoton", e_h[0] > e_h[1] && e_h[1] > e_h[2] &&
                                                      e_h[2] > e_h[3] && e_h[3] > e_h[4]);
    expect("    Flaechenfehler faellt monoton", e_A[0] > e_A[1] && e_A[1] > e_A[2] &&
                                                    e_A[2] > e_A[3] && e_A[3] > e_A[4]);
    expect("    Volumenfehler faellt monoton", e_V[0] > e_V[1] && e_V[1] > e_V[2] &&
                                                   e_V[2] > e_V[3] && e_V[3] > e_V[4]);
  }

  // The automatic discretisation must actually meet what it was asked for.
  CapillaryRequest q;
  q.delta_p_exit = capillary::pressure_from_pi(1.3, a, L.gamma);
  q.target_relative_accuracy = 1.0e-9;
  const CapillaryMeniscus m = solve_capillary_meniscus(a, 0.0, L, q);
  expect("automatische Aufloesung: Status Solved", is_usable(m.status));
  expect("automatische Aufloesung: Schaetzfehler <= Vorgabe",
         m.estimated_relative_error <= 1.0e-9);
  expect("automatische Aufloesung: Intervalle selbst gewaehlt",
         m.n_intervals >= 128 && !m.discretisation_was_forced);
  std::printf("    gewaehlt: %d Intervalle, Schaetzfehler %.2e\n", m.n_intervals,
              m.estimated_relative_error);
}

// ===========================================================================
// 8.  Scaling: geometrically similar problems give similar shapes
// ===========================================================================
static void test_scaling() {
  std::printf("\n=== 8. Skalierungspruefung bei proportional veraenderter Geometrie ===\n");
  LiquidProperties L1 = test_liquid();
  LiquidProperties L2 = test_liquid();
  L2.gamma = 0.037;                       // a different surface tension as well
  const Real a1 = 5.0e-6, a2 = 3.7e-4;    // 74x larger bore

  for (Real Pi : {-1.6, -0.4, 0.4, 1.6}) {
    const CapillaryMeniscus m1 = solve_pi(a1, Pi, L1);
    const CapillaryMeniscus m2 = solve_pi(a2, Pi, L2);
    char what[96];
    std::snprintf(what, sizeof what, "Pi = %+5.2f: h/a gleich", Pi);
    check(what, m2.apex_height / a2, m1.apex_height / a1, tol::kScaling);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: s/a gleich", Pi);
    check(what, m2.arclength / a2, m1.arclength / a1, tol::kScaling);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: A/a^2 gleich", Pi);
    check(what, m2.revolved_area / (a2 * a2), m1.revolved_area / (a1 * a1), tol::kScaling);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: V/a^3 gleich", Pi);
    check(what, m2.revolved_volume / (a2 * a2 * a2), m1.revolved_volume / (a1 * a1 * a1),
          tol::kScaling);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: psi(Kontakt) gleich", Pi);
    check(what, m2.contact_tangent_angle, m1.contact_tangent_angle, tol::kScaling);
  }
}

// ===========================================================================
// 9.  Invalid material data and an unrepresentable pressure fail closed
// ===========================================================================
static void test_fail_closed() {
  std::printf("\n=== 9. Ungueltige Stoffdaten und nicht darstellbarer Druck ===\n");
  const Real a = 5.0e-6;
  const LiquidProperties L = test_liquid();

  // --- pressure beyond the hemispherical limit ------------------------------
  for (Real Pi : {-4.0, -2.05, 2.05, 4.0}) {
    const CapillaryMeniscus m = solve_pi(a, Pi, L);
    char what[96];
    std::snprintf(what, sizeof what, "Pi = %+5.2f: PressureOutsideCapillaryRange", Pi);
    expect(what, m.status == CapillaryStatus::PressureOutsideCapillaryRange);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: keine Form zurueckgegeben", Pi);
    expect(what, m.nodes.empty() && m.apex_height == 0.0 && m.arclength == 0.0);
    std::snprintf(what, sizeof what, "Pi = %+5.2f: Begruendung vorhanden", Pi);
    expect(what, !m.message.empty());
  }

  // The physical criterion and the geometric detection must agree: the solver
  // decides from the integrated shape, the closed form from |Pi| <= 2.
  for (Real Pi : {-2.5, -1.99, 1.99, 2.5}) {
    const Real dp = capillary::pressure_from_pi(Pi, a, L.gamma);
    const CapillaryMeniscus m = solve_pi(a, Pi, L);
    char what[96];
    std::snprintf(what, sizeof what, "Pi = %+5.2f: Loeser und Existenzkriterium einig", Pi);
    expect(what, is_usable(m.status) == spherical_cap_exists(a, dp, L.gamma));
  }

  // The hemispherical limit itself is its own answer, not a silent shape.
  {
    const CapillaryMeniscus m = solve_pi(a, 2.0, L);
    expect("Pi = 2 exakt: HemisphericalLimit", m.status == CapillaryStatus::HemisphericalLimit);
    expect("Pi = 2 exakt: keine Form zurueckgegeben", m.nodes.empty());
  }

  // --- unusable material data ----------------------------------------------
  {
    LiquidProperties bad = test_liquid();
    bad.gamma = 0.0;
    const CapillaryMeniscus m = solve_pi(a, 1.0, bad);
    expect("gamma = 0: InvalidLiquid", m.status == CapillaryStatus::InvalidLiquid);
    expect("gamma = 0: keine Form zurueckgegeben", m.nodes.empty());
  }
  {
    LiquidProperties bad = test_liquid();
    bad.rho = -1.0;
    CapillaryRequest q;
    q.delta_p_exit = 1.0;
    const CapillaryMeniscus m = solve_capillary_meniscus(a, 0.0, bad, q);
    expect("rho < 0: InvalidLiquid", m.status == CapillaryStatus::InvalidLiquid);
  }
  {
    LiquidProperties bad;             // default: status Unknown, no numbers
    CapillaryRequest q;
    q.delta_p_exit = 1.0;
    const CapillaryMeniscus m = solve_capillary_meniscus(a, 0.0, bad, q);
    expect("Status unknown: InvalidLiquid", m.status == CapillaryStatus::InvalidLiquid);
    expect("Status unknown: why_unusable() benennt den Grund", !bad.usable() &&
                                                                  !m.message.empty());
  }
  {
    LiquidProperties bad = test_liquid();
    bad.status = LiquidDataStatus::Provisional;
    bad.source.clear();
    expect("provisional ohne Fundstelle ist unbrauchbar", !bad.usable());
  }

  // --- invalid geometry ------------------------------------------------------
  {
    CapillaryRequest q;
    q.delta_p_exit = 1.0;
    const CapillaryMeniscus m = solve_capillary_meniscus(0.0, 0.0, test_liquid(), q);
    expect("a = 0: InvalidGeometry", m.status == CapillaryStatus::InvalidGeometry);
  }

  // --- contact angle and pinning together ------------------------------------
  {
    CapillaryRequest q;
    q.delta_p_exit = capillary::pressure_from_pi(1.0, a, L.gamma);
    q.contact_angle_prescribed = true;
    q.prescribed_contact_angle_deg = 40.0;
    const CapillaryMeniscus m = solve_capillary_meniscus(a, 0.0, L, q);
    expect("Kontaktwinkel + Pinning: abgelehnt",
           m.status == CapillaryStatus::ContactAngleAndPinningBothPrescribed);
    expect("Kontaktwinkel + Pinning: keine Form zurueckgegeben", m.nodes.empty());
  }

  // --- the closed form refuses the same pressures ----------------------------
  {
    bool threw = false;
    try {
      spherical_cap(a, capillary::pressure_from_pi(2.5, a, L.gamma), L.gamma);
    } catch (const std::exception&) {
      threw = true;
    }
    expect("geschlossene Form wirft bei |Pi| > 2", threw);
  }
}

// ===========================================================================
// 10.  The exit edge comes from the device geometry, and the liquid record
// ===========================================================================
static void test_device_coupling_and_liquid() {
  std::printf("\n=== 10. Austrittskante aus DeviceParameters, Stoffdatensatz ===\n");
  DeviceParameters p;                       // phi_2 = 10 um by default
  p.extractor_outer_radius = 2.0e-3;
  const DeviceGeometry g = DeviceGeometry::build(p);
  const LiquidProperties L = emibf4_illustrative();

  CapillaryRequest q;
  q.delta_p_exit = capillary::pressure_from_pi(1.0, g.contact_radius(), L.gamma);
  q.target_relative_accuracy = tol::kRequestedAccuracy;
  const CapillaryMeniscus m = solve_capillary_meniscus(g, L, q);

  expect("Status Solved", is_usable(m.status));
  check("Pinningradius = phi_2/2", m.contact_radius, 0.5 * p.phi_2, 0.0);
  check_abs("Kontaktlinie bei z = 0", m.contact_z, 0.0, 0.0);
  check("Pi = 1 wie gefordert", m.Pi, 1.0, 1e-14);

  const SphericalCap c = spherical_cap(m.contact_radius, m.delta_p_exit, L.gamma);
  check("Apexhoehe gegen Kugelkappe", m.apex_height, c.apex_height, tol::kAnalytic);

  // A scaled bore must move the answer, and proportionally.
  DeviceParameters p2 = p;
  p2.phi_2 = 2.0 * p.phi_2;
  p2.phi_1 = 2.0 * p.phi_1;
  const DeviceGeometry g2 = DeviceGeometry::build(p2);
  CapillaryRequest q2 = q;
  q2.delta_p_exit = capillary::pressure_from_pi(1.0, g2.contact_radius(), L.gamma);
  const CapillaryMeniscus m2 = solve_capillary_meniscus(g2, L, q2);
  check("doppelte Bohrung, gleiches Pi: h verdoppelt sich", m2.apex_height,
        2.0 * m.apex_height, tol::kScaling);

  // --- the material record --------------------------------------------------
  expect("EMI-BF4 ist ausdruecklich illustrative",
         L.status == LiquidDataStatus::Illustrative);
  expect("EMI-BF4 nennt die Fundstelle der Stoffidentitaet",
         L.source.find("2.3.2") != std::string::npos);
  expect("EMI-BF4 traegt keine quantitative Aussage",
         !carries_quantitative_claim(L.status));
  expect("mu, K, eps_r sind nur vorgemerkt",
         L.documented_only.mu > 0.0 && L.documented_only.K > 0.0 &&
             L.documented_only.eps_r > 0.0);

  // Bond number of the reference geometry -- the justification for leaving
  // gravity out, computed rather than assumed.
  const Real Bo = L.bond_number(g.contact_radius());
  std::printf("    Bond-Zahl der Referenzgeometrie: %.4e (a = %.3g m)\n", Bo,
              g.contact_radius());
  expect("Bond-Zahl << 1", Bo < 1.0e-4);
  check("Bond-Zahl = rho g a^2 / gamma",
        Bo, L.rho * constants::g0 * g.contact_radius() * g.contact_radius() / L.gamma, 1e-15);
}

// ===========================================================================
int main() {
  std::printf("P3a -- statischer Kapillarmeniskus ohne elektrisches Feld\n");
  std::printf("Alle Toleranzen sind vor der Auswertung festgelegt (siehe namespace tol).\n");
  test_flat();
  test_pinning();
  test_axis_regularity();
  test_against_spherical_cap();
  test_residual();
  test_mesh_convergence();
  test_scaling();
  test_fail_closed();
  test_device_coupling_and_liquid();
  std::printf("\n%s: %d Fehler\n", failures ? "FEHLGESCHLAGEN" : "BESTANDEN", failures);
  return failures ? 1 : 0;
}
