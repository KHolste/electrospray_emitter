// es_beam -- beam transport from the emitting surface through the extractor.
//
// Laplace field only.  Space charge is disabled (phase P4), droplets are locked
// until the cone-jet coupling (phase P6), and only the cation polarity is
// modelled -- all three fail closed rather than producing numbers.

#include <cstdio>
#include <stdexcept>

#include "es/beam.hpp"
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
      std::printf("usage: es_beam [file.cfg] [key=value ...]\n\n"
                  "Traces the beam at solve.voltage in the Laplace field.\n"
                  "beam.space_charge_iters > 0, beam.species = droplet|both and a\n"
                  "negative polarity are rejected with an explanatory message.\n\n");
      print_key_reference(stdout);
      return 0;
    }
    cfg.load(a);
  }
  cfg.apply_cli(argc, argv);

  Setup s = build_setup(cfg);
  s.print(stdout);

  // Refuse the unmodelled polarity here, with its own reason, rather than
  // letting it surface later as an unrelated symptom.
  require_modelled_polarity(s.voltage);
  const std::string prefix = cfg.str("output.prefix", "out");

  BemSolver dry;
  MeniscusSolver mstore(s.electrodes, MeniscusParams{});
  BemSolver* bem = nullptr;
  std::string state;
  Real state_voltage = s.voltage;

  const bool use_meniscus = s.wetted && cfg.flag("beam.use_meniscus", true);
  if (use_meniscus) {
    MeniscusParams mp = meniscus_params_from(cfg, s);
    mstore = MeniscusSolver(s.electrodes, mp);
    const Real h_max = cfg.num("meniscus.h_max", 2.5) * mp.r_contact;
    std::printf("\nloese den Meniskus bei U = %.1f V ...\n", s.voltage);
    MeniscusSolution m = mstore.solve_at_voltage(s.voltage, h_max);
    if (!m.ok()) {
      std::fprintf(stderr, "\nKEINE verwertbare Meniskusloesung: %s\n  %s\n",
                   to_string(m.status), explain(m.status));
      if (m.status == SolveStatus::VoltageMismatch)
        std::fprintf(stderr, "  angefordert %.1f V, erreicht %.1f V\n", m.target_voltage,
                     m.voltage);
      std::fprintf(stderr, "  Es wird kein Strahl gerechnet. Zuerst es_meniscus ausfuehren.\n");
      return 2;
    }
    mstore.realize(m);
    bem = &mstore.bem();
    state = "meniscus";
    state_voltage = m.voltage;
    std::printf("  Apexhoehe %.4g m, Apexfeld %.4g V/m (%.4f V/nm)\n", m.shape.height,
                m.apex_field, m.apex_field * 1e-9);
    write_shape_csv(m.shape, output_path(prefix, "beam", state, state_voltage, "shape"),
                    meta_header("es_beam", "meniscus the beam was launched from", state_voltage));
  } else {
    dry = BemSolver(s.electrodes);
    dry.solve({s.voltage, 0.0, 0.0});
    bem = &dry;
    state = "drysurface";
    std::printf("\ntrockene Geometrie: Start von der Metalloberflaeche, "
                "Spitzenfeld %.4g V/m\n", bem->peak_emitter_field());
  }

  // --- species -------------------------------------------------------------
  std::vector<BeamSpecies> species;
  const std::string which = cfg.str("beam.species", "ion");
  if (which == "ion") {
    species.push_back({"ion", s.fluid.qm_cluster(), 1.0, SpeciesKind::IonEvaporated});
  } else if (which == "droplet" || which == "both") {
    // Deliberately constructed so that trace_beam refuses; the message names
    // the reason and the phase rather than silently dropping the request.
    if (which == "both")
      species.push_back({"ion", s.fluid.qm_cluster(), 0.5, SpeciesKind::IonEvaporated});
    species.push_back({"droplet", cfg.num("beam.droplet_qm", 1.0e4), 1.0, SpeciesKind::Droplet});
  } else {
    throw std::runtime_error("beam.species must be ion, droplet or both");
  }

  BeamParams bp;
  bp.z_end = cfg.num("beam.z_end", s.gap > 0 ? 4.0 * s.gap : 2.0e-3);
  bp.r_max = cfg.num("beam.r_max", 1.0e-2);
  bp.cfl = cfg.num("beam.cfl", 0.04);
  bp.max_steps = cfg.integer("beam.max_steps", 40000);
  bp.path_samples = cfg.integer("beam.path_samples", 80);
  bp.space_charge_iters = cfg.integer("beam.space_charge_iters", 0);
  bp.include_wetted_metal = cfg.flag("beam.wetted_metal", !use_meniscus);
  bp.verbose = cfg.flag("beam.verbose", true);

  std::printf("\nverfolge bis z = %.4g m (Laplace-Feld, keine Raumladung) ...\n", bp.z_end);
  BeamResult res = trace_beam(*bem, s.fluid, s.temperature, species, bp);
  res.print(stdout);

  Real on_extractor = 0.0, on_emitter = 0.0;
  for (const Ray& r : res.rays) {
    if (r.status != RayStatus::Intercepted) continue;
    if (r.hit_tag == Tag::Extractor) on_extractor += r.current; else on_emitter += r.current;
  }
  std::printf("  davon auf den Extraktor: %10.4g A\n", on_extractor);
  std::printf("  davon zurueck auf den Emitter: %10.4g A\n", on_emitter);

  std::printf("\n  Die Startgewichte stammen aus der Iribarne-Thomson-Rate auf dem Feld\n"
              "  eines nicht emittierenden Meniskus. Der absolute Strom ist deshalb eine\n"
              "  Abschaetzung, keine Vorhersage; die Strahlgeometrie ist davon unberuehrt.\n");

  const std::string hdr = meta_header("es_beam", state + " / Laplace, no space charge",
                                      state_voltage, "ion species only");
  res.write_rays_csv(output_path(prefix, "beam", state, state_voltage, "rays"));
  res.write_paths_csv(output_path(prefix, "beam", state, state_voltage, "paths"));
  bem->mesh().write_csv(output_path(prefix, "beam", state, state_voltage, "mesh"), hdr);
  bem->write_surface_csv(output_path(prefix, "beam", state, state_voltage, "surface"), hdr);
  std::printf("\ngeschrieben mit Praefix '%s'\n", prefix.c_str());

  cfg.warn_about_unused(stdout, {"operate.", "sweep."});
  return 0;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\nes_beam: %s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_beam: %s\n", e.what());
  return 1;
}
