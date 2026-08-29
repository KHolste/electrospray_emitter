// es_meniscus -- static meniscus branch and its turning point.
//
// Traces the equilibrium branch by prescribing the apex height and solving for
// the voltage that sustains it.  The maximum of U(h) is reported as the STATIC
// FOLD.  It is not an emission onset: see docs/02_model_specification.md 2.4.

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
  const std::vector<std::string> rest = Config::positional_args(argc, argv);
  for (const std::string& a : rest) {
    if (a == "--help" || a == "-h") {
      std::printf("usage: es_meniscus [file.cfg] [key=value ...]\n\n"
                  "Traces the static meniscus branch and reports its turning point\n"
                  "(static fold).  This is NOT an emission onset.\n"
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

  std::printf("\ntracing the static branch, h = %.3g .. %.3g m (%d steps)\n", h_min, h_max, steps);
  MeniscusSolver solver(s.electrodes, mp);
  const std::vector<MeniscusSolution> branch = solver.continuation(h_min, h_max, steps);

  std::printf("\n   h/r_c        U [V]    E_apex [V/m]   R_apex/r_c  half-angle  status\n");
  for (const MeniscusSolution& m : branch)
    std::printf("  %6.3f  %11.1f  %14.4g  %11.4f  %8.2f  %s (%d it)\n",
                m.shape.height / mp.r_contact, m.voltage, m.apex_field,
                m.shape.apex_radius / mp.r_contact, m.shape.half_angle * 180.0 / pi,
                to_string(m.status), m.iterations);

  const MeniscusSolver::StaticFold fold = MeniscusSolver::find_static_fold(branch);
  const std::string prefix = cfg.str("output.prefix", "out");

  if (!fold.found()) {
    std::fprintf(stderr, "\nKEIN Umkehrpunkt nachgewiesen: %s\n  %s\n", to_string(fold.status),
                 explain(fold.status));
    std::fprintf(stderr, "  Es wird keine Faltenspannung ausgegeben. Bereich ueber "
                         "meniscus.h_max erweitern oder meniscus.tol lockern.\n");
    MeniscusSolver::write_branch_csv(
        output_path(prefix, "meniscus", "branch-nofold", s.voltage, "branch"), branch,
        meta_header("es_meniscus", "branch without demonstrated turning point", s.voltage,
                    to_string(fold.status)));
    return 2;
  }

  std::printf("\nstatischer Umkehrpunkt (static fold)\n");
  std::printf("  Faltenspannung        : %10.1f V\n", fold.voltage);
  std::printf("  Apexhoehe             : %10.4g m  (= %.3f r_c)\n", fold.height,
              fold.height / mp.r_contact);
  std::printf("  Apexkruemmungsradius  : %10.4g m  (= %.3f r_c)\n", fold.apex_radius,
              fold.apex_radius / mp.r_contact);
  std::printf("  Apexfeld              : %10.4g V/m (= %.4f V/nm)\n", fold.apex_field,
              fold.apex_field * 1e-9);
  std::printf("  Kegelhalbwinkel       : %10.2f deg\n", fold.half_angle * 180.0 / pi);
  if (s.gap > 0.0) {
    const Real vt = literature_onset_voltage_smith(mp.r_contact, s.gap, s.fluid.gamma);
    std::printf("\n  Literaturformel Smith (1986), NICHT aus diesem Modell:\n");
    std::printf("    Emissions-Onset     : %10.1f V   (Verhaeltnis %.3f)\n", vt,
                fold.voltage / vt);
    std::printf("    Das ist eine andere Groesse als die Faltenspannung. Der Vergleich\n"
                "    zeigt nur die Groessenordnung, er ist kein Nachweis.\n");
  }
  std::printf("\n  WAS DIESE ZAHL NICHT IST: kein Emissions-Onset und kein nachgewiesener\n"
              "  Stabilitaetsverlust. Eine Stabilitaetsanalyse findet nicht statt; sie ist\n"
              "  fuer Phase P3 vorgesehen. Alles jenseits des Maximums ist der zweite Ast.\n");

  // --- output, bound to states that were actually realised -----------------
  const MeniscusSolution& at_fold = branch[fold.index];
  MeniscusSolver::write_branch_csv(
      output_path(prefix, "meniscus", "branch", fold.voltage, "branch"), branch,
      meta_header("es_meniscus", "full traced branch", fold.voltage,
                  "maximum over converged rows is the static fold"));

  solver.realize(at_fold);
  const std::string hdr = meta_header("es_meniscus", "static_fold", at_fold.voltage,
                                      "sampled branch point at the turning point");
  write_shape_csv(at_fold.shape,
                  output_path(prefix, "meniscus", "fold", at_fold.voltage, "shape"), hdr);
  solver.bem().write_surface_csv(
      output_path(prefix, "meniscus", "fold", at_fold.voltage, "surface"), hdr);
  solver.bem().mesh().write_csv(
      output_path(prefix, "meniscus", "fold", at_fold.voltage, "mesh"), hdr);

  // Last converged point of the branch, as a separate, separately named state.
  for (auto it = branch.rbegin(); it != branch.rend(); ++it) {
    if (!it->ok()) continue;
    solver.realize(*it);
    const std::string h2 = meta_header("es_meniscus", "last_converged", it->voltage,
                                       "beyond the fold: second branch, not realisable "
                                       "at fixed voltage");
    write_shape_csv(it->shape,
                    output_path(prefix, "meniscus", "last", it->voltage, "shape"), h2);
    solver.bem().write_surface_csv(
        output_path(prefix, "meniscus", "last", it->voltage, "surface"), h2);
    break;
  }

  std::printf("\ngeschrieben mit Anwendung, Zustand und Spannung im Dateinamen, Praefix '%s'\n",
              prefix.c_str());

  cfg.warn_about_unused(stdout, {"beam.", "operate."});
  return 0;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\nes_meniscus: %s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_meniscus: %s\n", e.what());
  return 1;
}
