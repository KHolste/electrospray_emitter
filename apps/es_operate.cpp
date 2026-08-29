// es_operate -- complete operating point of a wetted emitter at a given voltage.
//
// Solves the meniscus that the applied voltage sustains, integrates the
// Iribarne-Thomson rate over that surface, compares against the cone-jet
// correlation at the requested flow rate, and reports thrust / Isp / efficiency.

#include <cstdio>
#include <stdexcept>

#include "es/constants.hpp"
#include "es/emission.hpp"
#include "es/io.hpp"
#include "es/meniscus.hpp"

using namespace es;
using constants::pi;

int main(int argc, char** argv) try {
  Config cfg;
  // Files first, command line second: an override on the command line must win.
  const std::vector<std::string> rest = Config::positional_args(argc, argv);
  for (const std::string& a : rest) {
    if (a == "--help" || a == "-h") {
      std::printf("usage: es_operate [file.cfg] [key=value ...]\n\n"
                  "Reports the full operating point at solve.voltage.\n"
                  "Extra keys: operate.flow_rate (e.g. 0.5nL/s), operate.sweep_max,\n"
                  "            operate.sweep_steps\n\n");
      print_key_reference(stdout);
      return 0;
    }
    cfg.load(a);
  }
  cfg.apply_cli(argc, argv);

  Setup s = build_setup(cfg);
  MeniscusParams mp = meniscus_params_from(cfg, s);
  s.print(stdout);

  const Real h_max = cfg.num("meniscus.h_max", 2.5) * mp.r_contact;
  MeniscusSolver solver(s.electrodes, mp);

  // ---- locate the onset first: it bounds everything below ------------------
  std::printf("\nlocating the onset ...\n");
  const std::vector<MeniscusSolution> branch =
      solver.continuation(0.15 * mp.r_contact, h_max, cfg.integer("meniscus.steps", 20));
  const MeniscusSolver::Onset on = MeniscusSolver::find_onset(branch);
  if (!on.found) throw std::runtime_error("no onset found -- widen meniscus.h_max");
  std::printf("  onset voltage         : %10.1f V\n", on.voltage);

  auto report_at = [&](Real U) {
    std::printf("\n=====================================================================\n");
    std::printf("operating point at U = %.1f V   (%.1f %% of onset)\n", U, 100.0 * U / on.voltage);
    std::printf("=====================================================================\n");
    if (U > on.voltage) {
      std::printf("  U exceeds the onset voltage: no static meniscus exists here.  The\n"
                  "  liquid is emitting; use the cone-jet block below for the droplet\n"
                  "  mode, and note that the static model cannot supply the apex field.\n");
      return;
    }
    MeniscusSolution m = solver.solve_at_voltage(U, h_max);
    if (!m.converged) {
      std::printf("  meniscus did not converge at this voltage.\n");
      return;
    }
    std::printf("\nmeniscus\n");
    std::printf("  apex height           : %10.4g m  (= %.3f r_c)\n", m.shape.height,
                m.shape.height / mp.r_contact);
    std::printf("  apex radius           : %10.4g m  (= %.3f r_c)\n", m.shape.apex_radius,
                m.shape.apex_radius / mp.r_contact);
    std::printf("  apex field            : %10.4g V/m (= %.4f V/nm)\n", m.apex_field,
                m.apex_field * 1e-9);
    std::printf("  cone half-angle       : %10.2f deg\n", m.shape.half_angle * 180.0 / pi);

    const IonEmission ion = integrate_ion_emission(solver.bem(), s.fluid, s.temperature,
                                                  cfg.flag("operate.wetted_metal", false));

    // Cone-jet comparison at the requested flow rate.
    const Real Q = cfg.num("operate.flow_rate", 0.0);
    ConeJetState cj;
    const bool have_flow = (Q > 0.0);
    if (have_flow) cj = cone_jet(s.fluid, Q);

    std::vector<Species> mix;
    if (ion.mdot > 0.0) mix.push_back({"ion clusters", ion.mdot, s.fluid.qm_cluster()});
    if (have_flow && cj.qm > 0.0) mix.push_back({"droplets", cj.mdot, cj.qm});
    BeamFigures fig;
    if (!mix.empty()) fig = beam_figures(mix, U);

    print_operating_point(stdout, s.fluid, have_flow ? &cj : nullptr, &ion,
                          mix.empty() ? nullptr : &fig);

    if (!mix.empty()) {
      std::printf("\n  beam composition\n");
      for (const Species& sp : mix)
        std::printf("    %-14s mdot = %10.4g kg/s, q/m = %10.4g C/kg, I = %10.4g A\n",
                    sp.name.c_str(), sp.mdot, sp.qm, sp.mdot * sp.qm);
    }
  };

  report_at(s.voltage);

  // ---- optional voltage sweep ---------------------------------------------
  if (cfg.has("operate.sweep_max")) {
    const Real v1 = cfg.num("operate.sweep_max", on.voltage);
    const Real v0 = cfg.num("operate.sweep_min", 0.6 * on.voltage);
    const int n = cfg.integer("operate.sweep_steps", 10);
    const std::string path = cfg.str("output.prefix", "operate") + "_iv.csv";
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot open " + path);
    std::fprintf(f, "voltage,apex_height,apex_radius,apex_field,ion_current,emitting_area,"
                    "thrust,Isp,converged\n");
    std::printf("\nI-V sweep\n");
    std::printf("  %10s %14s %14s %12s %10s\n", "U [V]", "E_apex [V/m]", "I_ion [A]", "F [uN]",
                "Isp [s]");
    for (int i = 0; i < n; ++i) {
      const Real U = v0 + (v1 - v0) * i / std::max(1, n - 1);
      MeniscusSolution m = solver.solve_at_voltage(U, h_max);
      IonEmission ion{};
      BeamFigures fig{};
      if (m.converged) {
        ion = integrate_ion_emission(solver.bem(), s.fluid, s.temperature,
                                     cfg.flag("operate.wetted_metal", false));
        if (ion.mdot > 0.0) fig = beam_figures({{"ion", ion.mdot, s.fluid.qm_cluster()}}, U);
      }
      std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%d\n", U, m.shape.height,
                   m.shape.apex_radius, m.apex_field, ion.current, ion.emitting_area, fig.thrust,
                   fig.Isp, m.converged ? 1 : 0);
      std::printf("  %10.1f %14.4g %14.4g %12.4g %10.1f%s\n", U, m.apex_field, ion.current,
                  fig.thrust * 1e6, fig.Isp, m.converged ? "" : "   (not converged)");
    }
    std::fclose(f);
    std::printf("wrote %s\n", path.c_str());
  }

  cfg.warn_about_unused(stdout, {"beam."});
  return 0;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_operate: %s\n", e.what());
  return 1;
}
