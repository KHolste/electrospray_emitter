// es_meniscus -- electrified meniscus shape and the onset of emission.
//
// Traces the equilibrium branch by prescribing the apex height and solving for
// the voltage that sustains it.  The maximum of U(h) is the onset voltage: the
// saddle-node beyond which no static meniscus exists.

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
      std::printf("usage: es_meniscus [file.cfg] [key=value ...]\n\n"
                  "Traces the meniscus equilibrium branch and reports the onset voltage.\n"
                  "Extra keys: meniscus.h_min, meniscus.h_max (in units of r_contact),\n"
                  "            meniscus.steps\n\n");
      print_key_reference(stdout);
      return 0;
    }
    cfg.load(a);
  }
  cfg.apply_cli(argc, argv);

  Setup s = build_setup(cfg);
  MeniscusParams mp = meniscus_params_from(cfg, s);
  s.print(stdout);

  const Real h_min = cfg.num("meniscus.h_min", 0.15) * mp.r_contact;
  const Real h_max = cfg.num("meniscus.h_max", 2.5) * mp.r_contact;
  const int steps = cfg.integer("meniscus.steps", 24);

  std::printf("\ntracing the equilibrium branch, h = %.3g .. %.3g m (%d steps)\n", h_min, h_max,
              steps);
  MeniscusSolver solver(s.electrodes, mp);
  const std::vector<MeniscusSolution> branch = solver.continuation(h_min, h_max, steps);

  std::printf("\n   h/r_c        U [V]    E_apex [V/m]   R_apex/r_c  half-angle   I_ion [A]  conv\n");
  for (const MeniscusSolution& m : branch) {
    // Ion current the apex field would drive, if this were a wetted emitter.
    Real Iion = 0.0;
    if (m.converged) {
      // Rough: apply the Iribarne-Thomson rate to the peak field over the
      // apex cap.  The full integral is done by es_operate.
      Iion = ion_current_density(m.apex_field, s.fluid, s.temperature) *
             (2.0 * pi * m.shape.apex_radius * m.shape.apex_radius);
    }
    std::printf("  %6.3f  %11.1f  %14.4g  %11.4f  %8.2f  %11.3g  %s\n", m.shape.height / mp.r_contact,
                m.voltage, m.apex_field, m.shape.apex_radius / mp.r_contact,
                m.shape.half_angle * 180.0 / pi, Iion, m.converged ? "yes" : "NO");
  }

  const MeniscusSolver::Onset on = MeniscusSolver::find_onset(branch);
  if (!on.found) {
    std::fprintf(stderr, "\nno turning point found -- widen meniscus.h_max or loosen "
                         "meniscus.tol\n");
    return 2;
  }

  std::printf("\nonset of emission\n");
  std::printf("  onset voltage         : %10.1f V\n", on.voltage);
  std::printf("  apex height at onset  : %10.4g m  (= %.3f r_c)\n", on.height,
              on.height / mp.r_contact);
  std::printf("  apex radius at onset  : %10.4g m  (= %.3f r_c)\n", on.apex_radius,
              on.apex_radius / mp.r_contact);
  std::printf("  apex field at onset   : %10.4g V/m (= %.4f V/nm)\n", on.apex_field,
              on.apex_field * 1e-9);
  std::printf("  cone half-angle       : %10.2f deg   (Taylor equilibrium: 49.29)\n",
              on.half_angle * 180.0 / pi);
  if (s.gap > 0.0) {
    const Real vt = onset_voltage_taylor(mp.r_contact, s.gap, s.fluid.gamma);
    std::printf("  Taylor/Smith estimate : %10.1f V   (ratio %.3f)\n", vt, on.voltage / vt);
  }
  std::printf("\n  Everything past the maximum of U(h) is the UNSTABLE branch: real\n"
              "  solutions of the equations, but not reachable at fixed voltage.\n");

  const std::string prefix = cfg.str("output.prefix", "meniscus");
  MeniscusSolver::write_branch_csv(prefix + "_branch.csv", branch);
  // Dump the shape at the turning point and at the last converged point.
  for (const MeniscusSolution& m : branch)
    if (m.converged && std::abs(m.shape.height - on.height) < 1e-12 + 0.5 * (h_max - h_min) / steps) {
      write_shape_csv(m.shape, prefix + "_shape_onset.csv");
      break;
    }
  for (auto it = branch.rbegin(); it != branch.rend(); ++it)
    if (it->converged) { write_shape_csv(it->shape, prefix + "_shape_last.csv"); break; }
  solver.bem().write_surface_csv(prefix + "_surface.csv");
  solver.bem().mesh().write_csv(prefix + "_mesh.csv");
  std::printf("\nwrote %s_branch.csv, %s_shape_onset.csv, %s_shape_last.csv,\n"
              "      %s_surface.csv, %s_mesh.csv\n",
              prefix.c_str(), prefix.c_str(), prefix.c_str(), prefix.c_str(), prefix.c_str());

  cfg.warn_about_unused(stdout, {"beam.", "operate."});
  return 0;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_meniscus: %s\n", e.what());
  return 1;
}
