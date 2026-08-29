// es_beam -- beam transport from the emitting surface through the extractor.
//
// Launches ring macroparticles weighted by the local Iribarne-Thomson emission
// rate, tracks them through the solved field (optionally with self-consistent
// space charge) and reports divergence, interception and beam energy.

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
  // Files first, command line second: an override on the command line must win.
  const std::vector<std::string> rest = Config::positional_args(argc, argv);
  for (const std::string& a : rest) {
    if (a == "--help" || a == "-h") {
      std::printf("usage: es_beam [file.cfg] [key=value ...]\n\n"
                  "Traces the beam at solve.voltage.  With emitter.type = capillary the\n"
                  "meniscus is solved first and the beam is launched from it; with a dry\n"
                  "geometry the beam is launched from the metal surface.\n\n");
      print_key_reference(stdout);
      return 0;
    }
    cfg.load(a);
  }
  cfg.apply_cli(argc, argv);

  Setup s = build_setup(cfg);
  s.print(stdout);

  BemSolver* bem = nullptr;
  BemSolver dry;
  MeniscusSolver* msolver = nullptr;
  MeniscusSolver mstore(s.electrodes, MeniscusParams{});

  const bool use_meniscus = s.wetted && cfg.flag("beam.use_meniscus", true);
  if (use_meniscus) {
    MeniscusParams mp = meniscus_params_from(cfg, s);
    mstore = MeniscusSolver(s.electrodes, mp);
    msolver = &mstore;
    const Real h_max = cfg.num("meniscus.h_max", 2.5) * mp.r_contact;
    std::printf("\nsolving the meniscus at U = %.1f V ...\n", s.voltage);
    MeniscusSolution m = msolver->solve_at_voltage(s.voltage, h_max);
    if (!m.converged)
      throw std::runtime_error("no static meniscus at this voltage -- it is above onset, or the "
                               "continuation did not converge.  Run es_meniscus first.");
    std::printf("  apex height %.4g m, apex field %.4g V/m (%.4f V/nm)\n", m.shape.height,
                m.apex_field, m.apex_field * 1e-9);
    bem = &msolver->bem();
    write_shape_csv(m.shape, cfg.str("output.prefix", "beam") + "_shape.csv");
  } else {
    dry = BemSolver(s.electrodes);
    dry.solve({s.voltage, 0.0, 0.0});
    bem = &dry;
    std::printf("\ndry geometry: launching from the metal surface, peak |E_n| = %.4g V/m\n",
                bem->peak_emitter_field());
  }

  // --- species mix ---------------------------------------------------------
  std::vector<BeamSpecies> species;
  const std::string which = cfg.str("beam.species", "ion");
  const Real drop_frac = cfg.num("beam.droplet_fraction", 0.0);
  if (which == "ion" || which == "both")
    species.push_back({"ion", s.fluid.qm_cluster(), which == "both" ? 1.0 - drop_frac : 1.0});
  if (which == "droplet" || which == "both")
    species.push_back({"droplet", cfg.num("beam.droplet_qm", 1.0e4),
                       which == "droplet" ? 1.0 : drop_frac});
  if (species.empty()) throw std::runtime_error("beam.species must be ion, droplet or both");

  BeamParams bp;
  bp.z_end = cfg.num("beam.z_end", s.gap > 0 ? 4.0 * s.gap : 2.0e-3);
  bp.r_max = cfg.num("beam.r_max", 1.0e-2);
  bp.cfl = cfg.num("beam.cfl", 0.04);
  bp.max_steps = cfg.integer("beam.max_steps", 40000);
  bp.path_samples = cfg.integer("beam.path_samples", 80);
  bp.space_charge_iters = cfg.integer("beam.space_charge_iters", 0);
  bp.space_charge_relax = cfg.num("beam.space_charge_relax", 0.5);
  bp.include_wetted_metal = cfg.flag("beam.wetted_metal", !use_meniscus);
  bp.verbose = cfg.flag("beam.verbose", true);

  std::printf("\ntracing to z = %.4g m ...\n", bp.z_end);
  BeamResult res = trace_beam(*bem, s.fluid, s.temperature, species, bp);
  res.print(stdout);

  // Per-surface interception breakdown.
  Real on_extractor = 0.0, on_emitter = 0.0;
  for (const Ray& r : res.rays) {
    if (r.status != RayStatus::Intercepted) continue;
    if (r.hit_tag == Tag::Extractor) on_extractor += r.current;
    else on_emitter += r.current;
  }
  std::printf("  hit the extractor   : %10.4g A\n", on_extractor);
  std::printf("  fell back on emitter: %10.4g A\n", on_emitter);

  if (res.current_launched > 0.0) {
    std::vector<Species> mix;
    for (const BeamSpecies& sp : species) {
      Real I = 0.0;
      for (const Ray& r : res.rays)
        if (r.species == sp.name && r.status == RayStatus::Transmitted) I += r.current;
      if (I > 0.0 && sp.qm > 0.0) mix.push_back({sp.name, I / sp.qm, sp.qm});
    }
    if (!mix.empty()) {
      const BeamFigures fig = beam_figures(mix, s.voltage);
      std::printf("\ntransmitted beam\n");
      std::printf("  thrust              : %10.4g N  (= %.4g uN)\n", fig.thrust,
                  fig.thrust * 1e6);
      std::printf("  specific impulse    : %10.1f s\n", fig.Isp);
      std::printf("  eta_polydispersity  : %10.4f\n", fig.eta_polydispersity);
      std::printf("  eta_transmission    : %10.4f\n",
                  res.current_transmitted / res.current_launched);
    }
  }

  const std::string prefix = cfg.str("output.prefix", "beam");
  res.write_rays_csv(prefix + "_rays.csv");
  res.write_paths_csv(prefix + "_paths.csv");
  bem->mesh().write_csv(prefix + "_mesh.csv");
  bem->write_surface_csv(prefix + "_surface.csv");
  std::printf("\nwrote %s_rays.csv, %s_paths.csv, %s_mesh.csv, %s_surface.csv\n", prefix.c_str(),
              prefix.c_str(), prefix.c_str(), prefix.c_str());

  cfg.warn_about_unused(stdout, {"operate.", "sweep."});
  return 0;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_beam: %s\n", e.what());
  return 1;
}
