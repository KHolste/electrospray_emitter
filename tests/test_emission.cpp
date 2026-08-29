#include <cmath>
#include <cstdio>

#include "es/constants.hpp"
#include "es/emission.hpp"
#include "es/fluid.hpp"

using namespace es;
using constants::eV;

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

static void test_fluids() {
  std::printf("\n=== fluid database ===\n");
  for (const std::string& n : fluid_names()) {
    const Fluid f = fluid_by_name(n);
    std::printf("  %-16s rho=%7.1f gamma=%.4f K=%6.3f eps_r=%6.1f  tau_e=%.2e s  r*=%.2e m\n",
                f.name.c_str(), f.rho, f.gamma, f.K, f.eps_r, f.charge_relaxation_time(),
                f.ehd_length());
    expect("charge relaxation is fast enough for the equipotential assumption",
           f.charge_relaxation_time() < 1e-6);
    // The cone-jet radius scale must come out nanometric, not micrometric --
    // a good tripwire for a dropped factor of eps0.
    expect("EHD length is in the nm-to-100nm range",
           f.ehd_length() > 1e-10 && f.ehd_length() < 1e-6);
  }
  // Name canonicalisation
  expect("EMI_BF4 / emibf4 / EMI-BF4 all resolve",
         fluid_by_name("EMI_BF4").name == fluid_by_name("emibf4").name &&
             fluid_by_name("emibf4").name == fluid_by_name("EMI-BF4").name);

  // Temperature dependence: conductivity must rise and viscosity fall with T.
  const Fluid a = fluid_by_name("EMI-BF4");
  const Fluid b = a.at_temperature(333.15);
  std::printf("  EMI-BF4 at 298 K: K = %.3f S/m, mu = %.4f Pa s\n", a.K, a.mu);
  std::printf("  EMI-BF4 at 333 K: K = %.3f S/m, mu = %.4f Pa s\n", b.K, b.mu);
  expect("conductivity increases with temperature", b.K > a.K);
  expect("viscosity decreases with temperature", b.mu < a.mu);
  expect("surface tension decreases with temperature", b.gamma < a.gamma);
}

static void test_schottky() {
  std::printf("\n=== ion evaporation ===\n");
  // The barrier lowering at 1 V/nm is the number that sets the ~1 V/nm rule of
  // thumb for the pure ionic regime.
  const Real G = schottky_lowering(1.0e9);
  std::printf("  Schottky lowering at 1 V/nm: %.4f eV\n", G / eV);
  check("G(1 V/nm)", G / eV, 1.1998, 1e-3);
  expect("G grows as sqrt(E)",
         std::abs(schottky_lowering(4.0e9) / schottky_lowering(1.0e9) - 2.0) < 1e-12);

  const Fluid f = fluid_by_name("EMI-BF4");
  const Real Echar = characteristic_evaporation_field(f);
  std::printf("  field at which G(E) = dG (%.2f eV): %.4g V/m = %.3f V/nm\n",
              f.dG_solvation / eV, Echar, Echar * 1e-9);
  expect("characteristic field is of order 1 V/nm", Echar > 5e8 && Echar < 2e9);

  // Round-trip the inverse.
  for (Real E : {6e8, 8e8, 1.0e9, 1.2e9}) {
    const Real j = ion_current_density(E, f, 300.0);
    const Real Eback = field_for_current_density(j, f, 300.0);
    char buf[96];
    std::snprintf(buf, sizeof buf, "field <-> current density round trip at %.1e V/m", E);
    check(buf, Eback, E, 1e-6);
  }

  // Steepness: the current density must climb many decades over a modest field
  // range.  This is the property that localises emission to the apex.
  const Real j1 = ion_current_density(0.9e9, f, 300.0);
  const Real j2 = ion_current_density(1.1e9, f, 300.0);
  std::printf("  j(0.9 V/nm) = %.3e A/m^2, j(1.1 V/nm) = %.3e A/m^2  (x%.1f)\n", j1, j2, j2 / j1);
  expect("j rises by more than an order of magnitude from 0.9 to 1.1 V/nm", j2 / j1 > 10.0);
}

static void test_cone_jet() {
  std::printf("\n=== cone-jet scaling ===\n");
  const Fluid f = fluid_by_name("formamide+NaI");  // high eps_r: inside the fit range
  const ConeJetState a = cone_jet(f, 1.0e-12);
  const ConeJetState b = cone_jet(f, 4.0e-12);
  std::printf("  Q = 1 nL/s: I = %.4g A, d = %.4g nm, q/m = %.4g C/kg\n", a.current,
              a.d_droplet * 1e9, a.qm);
  std::printf("  Q = 4 nL/s: I = %.4g A, d = %.4g nm, q/m = %.4g C/kg\n", b.current,
              b.d_droplet * 1e9, b.qm);
  check("I scales as sqrt(Q)", b.current / a.current, 2.0, 1e-12);
  check("d scales as Q^(1/3)", b.d_droplet / a.d_droplet, std::cbrt(4.0), 1e-12);
  check("q/m scales as Q^(-1/2)", b.qm / a.qm, 0.5, 1e-12);
  expect("high-permittivity liquid is not flagged as extrapolated", !a.extrapolated);
  expect("ionic liquid IS flagged as extrapolated",
         cone_jet(fluid_by_name("EMI-BF4"), 1e-13).extrapolated);

  // Rayleigh limit consistency: a droplet at its Rayleigh charge has fissility 1.
  const Real d = 100e-9;
  const Real qR = rayleigh_charge(d, f.gamma);
  const Real R = 0.5 * d;
  check("Rayleigh charge = 8 pi sqrt(eps0 gamma R^3)", qR,
        8.0 * constants::pi * std::sqrt(constants::eps0 * f.gamma * R * R * R), 1e-14);
}

static void test_beam_figures() {
  std::printf("\n=== beam figures of merit ===\n");
  const Real V = 2000.0;
  // Single species: polydispersive efficiency must be exactly 1.
  {
    const BeamFigures b = beam_figures({{"ions", 1e-12, 5.0e5}}, V);
    const Real v = std::sqrt(2.0 * 5.0e5 * V);
    check("thrust = mdot * v", b.thrust, 1e-12 * v, 1e-14);
    check("Isp = v / g0", b.Isp, v / constants::g0, 1e-14);
    check("eta_pol = 1 for one species", b.eta_polydispersity, 1.0, 1e-14);
  }
  // Mixed droplet/ion beam: efficiency must drop below 1.
  {
    const BeamFigures b =
        beam_figures({{"droplets", 5e-12, 1.0e4}, {"ions", 1e-12, 5.0e5}}, V);
    std::printf("  mixed beam: I = %.4g A, F = %.4g uN, Isp = %.1f s, eta_pol = %.4f\n",
                b.current, b.thrust * 1e6, b.Isp, b.eta_polydispersity);
    expect("mixed beam loses polydispersive efficiency",
           b.eta_polydispersity > 0.0 && b.eta_polydispersity < 1.0);
    check("current = sum mdot * q/m", b.current, 5e-12 * 1e4 + 1e-12 * 5e5, 1e-14);
  }

  // A pure-ionic EMI-BF4 beam should land in the few-thousand-second Isp range
  // that colloid thrusters are known for.
  const Fluid f = fluid_by_name("EMI-BF4");
  const BeamFigures pir = beam_figures({{"cluster", 1e-12, f.qm_cluster()}}, 1500.0);
  std::printf("  EMI-BF4 PIR at 1500 V: q/m = %.4g C/kg, Isp = %.0f s\n", f.qm_cluster(), pir.Isp);
  expect("PIR Isp is in the 1000-6000 s range", pir.Isp > 1000.0 && pir.Isp < 6000.0);
}

static void test_onset_formulas() {
  std::printf("\n=== onset closed forms ===\n");
  const Fluid f = fluid_by_name("EMI-BF4");
  for (Real rc : {5e-6, 1e-5, 5e-5}) {
    std::printf("  r_c = %5.1f um: E_hemisphere = %.3e V/m, V_Taylor(d=0.5mm) = %7.1f V\n",
                rc * 1e6, onset_field_hemisphere(rc, f.gamma),
                onset_voltage_taylor(rc, 5e-4, f.gamma));
  }
  // Both must scale as sqrt(r_c) in the appropriate sense.
  check("E_hemisphere scales as 1/sqrt(r)",
        onset_field_hemisphere(4e-6, f.gamma) / onset_field_hemisphere(1e-6, f.gamma), 0.5, 1e-12);
}

int main() {
  test_fluids();
  test_schottky();
  test_cone_jet();
  test_beam_figures();
  test_onset_formulas();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
