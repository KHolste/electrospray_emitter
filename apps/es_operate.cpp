// es_operate -- non-coupled diagnostic estimate for a wetted emitter.
//
// This program does NOT compute an operating point.  It solves a static,
// perfectly conducting, non-emitting meniscus at the requested voltage and then
// applies the Iribarne-Thomson rate to that field.  In the pure ionic regime
// the flow is viscosity-dominated and the current is controlled by the finite
// conductivity (Higuera 2008), so the field used here is not the field that
// would actually be present.  The self-consistent model is due in phase P5.

#include <cstdio>
#include <stdexcept>

#include "es/constants.hpp"
#include "es/emission.hpp"
#include "es/io.hpp"
#include "es/meniscus.hpp"

using namespace es;
using constants::pi;

/// Read the branch selection from the configuration.  The application must make
/// this choice explicitly and say so in its output: the solver refuses to pick
/// when a voltage admits more than one meniscus.
static BranchSide branch_side_from(const Config& cfg) {
  const std::string v = cfg.str("meniscus.branch", "lower");
  if (v == "lower") return BranchSide::LowerHeight;
  if (v == "upper") return BranchSide::UpperHeight;
  if (v == "none") return BranchSide::Unspecified;
  throw std::runtime_error("meniscus.branch must be lower, upper or none");
}

/// Report the outcome of a voltage solve, including how many solutions the
/// traced branch offered.
static void report_branch_choice(const MeniscusSolution& m, BranchSide side) {
  std::printf("\n  Astabdeckung: %d Schnittpunkt(e) im verfolgten Bereich, "
              "Abdeckung %s,\n  Ende des Bereichs: %s. Gewaehlt: '%s' (meniscus.branch).\n",
              m.crossings_in_range, m.coverage_complete ? "vollstaendig" : "UNVOLLSTAENDIG",
              to_string(m.termination_reason), to_string(side));
  if (m.additional_crossing_possible)
    std::printf("  WARNUNG: der Ast wurde nicht weit genug verfolgt, um weitere Loesungen\n"
                "  auszuschliessen. meniscus.h_max vergroessern; scout_steps aendert nur\n"
                "  die Abtastung, nicht die Abdeckung.\n");
  if (m.crossings_in_range > 1)
    std::printf("  Die Bezeichnungen betreffen nur die Apexhoehe, nicht Stabilitaet --\n"
                "  eine Stabilitaetsanalyse gibt es nicht.\n");
}

int main(int argc, char** argv) try {
  Config cfg;
  const std::vector<std::string> rest = Config::positional_args(argc, argv);
  for (const std::string& a : rest) {
    if (a == "--help" || a == "-h") {
      std::printf("usage: es_operate [file.cfg] [key=value ...]\n\n"
                  "Reports a NON-COUPLED DIAGNOSTIC ESTIMATE at solve.voltage.\n"
                  "Not an operating point, not a current prediction.\n"
                  "Extra keys: operate.flow_rate (cone-jet correlation, printed\n"
                  "            separately and marked empirical)\n\n");
      print_key_reference(stdout);
      return 0;
    }
    cfg.load(a);
  }
  cfg.apply_cli(argc, argv);

  Setup s = build_setup(cfg);
  MeniscusParams mp = meniscus_params_from(cfg, s);
  s.print(stdout);

  // Refuse the unmodelled polarity here, with its own reason, rather than
  // letting it surface later as an unrelated symptom.
  require_modelled_polarity(s.voltage);

  const Real h_max = cfg.num("meniscus.h_max", 2.5) * mp.r_contact;
  const std::string prefix = cfg.str("output.prefix", "out");
  MeniscusSolver solver(s.electrodes, mp);

  std::printf("\nsuche den statischen Umkehrpunkt ...\n");
  const std::vector<MeniscusSolution> branch =
      solver.continuation(0.15 * mp.r_contact, h_max, cfg.integer("meniscus.steps", 20));
  const MeniscusSolver::StaticFold fold = MeniscusSolver::find_static_fold(branch);
  if (!fold.found()) {
    std::fprintf(stderr, "\nKEIN Umkehrpunkt nachgewiesen: %s\n  %s\n", to_string(fold.status),
                 explain(fold.status));
    return 2;
  }
  std::printf("  Faltenspannung        : %10.1f V  (kein Emissions-Onset)\n", fold.voltage);

  // --- meniscus at the requested voltage -----------------------------------
  const Real U = s.voltage;
  std::printf("\nloese den Meniskus bei U = %.1f V (%.1f %% der Faltenspannung) ...\n", U,
              100.0 * U / fold.voltage);
  const BranchSide side = branch_side_from(cfg);
  MeniscusSolution m = solver.solve_at_voltage(U, h_max, side);
  if (!m.ok()) {
    std::fprintf(stderr, "\nKEINE verwertbare Loesung: %s\n  %s\n", to_string(m.status),
                 explain(m.status));
    if (m.status == SolveStatus::VoltageMismatch)
      std::fprintf(stderr, "  angefordert %.1f V, erreicht %.1f V\n", m.target_voltage, m.voltage);
    std::fprintf(stderr, "  Es werden keine Zahlen ausgegeben.\n");
    return 2;
  }
  solver.realize(m);
  report_branch_choice(m, side);

  std::printf("\nMeniskus (statisch, perfekt leitend, nicht emittierend)\n");
  std::printf("  Apexhoehe             : %10.4g m  (= %.3f r_c)\n", m.shape.height,
              m.shape.height / mp.r_contact);
  std::printf("  Apexkruemmungsradius  : %10.4g m  (= %.3f r_c)\n", m.shape.apex_radius,
              m.shape.apex_radius / mp.r_contact);
  std::printf("  Apexfeld              : %10.4g V/m (= %.4f V/nm)\n", m.apex_field,
              m.apex_field * 1e-9);
  std::printf("  Kegelhalbwinkel       : %10.2f deg\n", m.shape.half_angle * 180.0 / pi);

  // --- diagnostic estimate --------------------------------------------------
  const IonEmission ion = integrate_ion_emission(solver.bem(), s.fluid, s.temperature,
                                                 cfg.flag("operate.wetted_metal", false));
  BeamFigures fig;
  std::vector<Species> mix;
  if (ion.mdot > 0.0) {
    mix.push_back({"ion clusters", ion.mdot, s.fluid.qm_cluster()});
    fig = beam_figures(mix, U);
  }
  print_diagnostic_estimate(stdout, s.fluid, &ion, mix.empty() ? nullptr : &fig);

  const Real E_evap = characteristic_evaporation_field(s.fluid);
  std::printf("\n  Einordnung: das Apexfeld betraegt %.3f V/nm, das Feld fuer G(E) = dG\n"
              "  liegt bei %.3f V/nm (Verhaeltnis %.3f). Die statische Rechnung loest die\n"
              "  nanometrische Emissionsstruktur nicht auf, aus der im PIR tatsaechlich\n"
              "  emittiert wird.\n",
              m.apex_field * 1e-9, E_evap * 1e-9, m.apex_field / E_evap);

  // --- cone-jet correlation, entirely separate ------------------------------
  const Real Q = cfg.num("operate.flow_rate", 0.0);
  if (Q > 0.0) print_cone_jet_correlation(stdout, s.fluid, cone_jet(s.fluid, Q));

  const std::string hdr =
      meta_header("es_operate", "static_meniscus_at_requested_voltage", m.voltage,
                  "diagnostic estimate only -- not an operating point");
  write_shape_csv(m.shape, output_path(prefix, "operate", "meniscus", m.voltage, "shape"), hdr);
  solver.bem().write_surface_csv(
      output_path(prefix, "operate", "meniscus", m.voltage, "surface"), hdr);
  std::printf("\ngeschrieben mit Praefix '%s'\n", prefix.c_str());

  cfg.warn_about_unused(stdout, {"beam."});
  return 0;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\nes_operate: %s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_operate: %s\n", e.what());
  return 1;
}
