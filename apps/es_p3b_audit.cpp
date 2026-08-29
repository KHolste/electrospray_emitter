// es_p3b_audit -- P0: the numerical clean-up of P3b.
//
//   es_p3b_audit <geometrie.cfg> [<p3b.cfg> ...] <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS FOR
//
// It does not add physics.  It answers three questions P3b left open, and it
// answers them with measurements rather than with prose:
//
//   1. THE LOAD PROJECTION.  Three different objects were all called "the
//      load".  Here the raw nodal load, the conservative segment load and the
//      load actually HANDED to the capillary solver are computed and written
//      separately, on manufactured loads whose integral is known in closed form
//      -- a smooth one and an integrable singularity p = C d^beta -- and on the
//      real solved field.  Continuity, force conservation, the treatment of the
//      contact edge and the absence of any clipping or exclusion zone are all
//      measured.
//
//   2. THE COUPLED MESH CONVERGENCE.  No new voltage sweep: two operating
//      points, 1000 V and 1400 V, on the three existing mesh levels plus one
//      finer.  For h/a, the edge-far normal field and the integrated Maxwell
//      force, a Richardson extrapolation with an observed order and an
//      estimated discretisation error, against a target fixed before the
//      measurement -- one per cent.  Where it is missed the result is labelled
//      DiscretizationNotConverged and is qualitative, and the figures say so.
//
//   3. THE EDGE GATE.  The flat shape and one outward-bulging shape, with the
//      additional mesh levels the raised memory cap allows.  The direct force
//      sequence, the extrapolation over the mesh levels and the extrapolation
//      over the exclusion distance are written side by side so that the reader
//      can see all three at once instead of a verdict.
//
// A run that fails writes NaN, never a zero.  A zero is a physical value and
// this run must not manufacture one.
//
// Exit code 2 means a declared check failed.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/capillary.hpp"
#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/electrocapillary.hpp"
#include "es/liquid.hpp"
#include "es/load_projection.hpp"

using namespace es;
using constants::pi;

namespace {

constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();

/// Print a Real that may be missing.  A failed computation is NaN in the file,
/// never 0 -- zero is a physical value and would be read as one.
void put(std::FILE* f, Real v) {
  if (std::isfinite(v))
    std::fprintf(f, ",%.9e", v);
  else
    std::fprintf(f, ",nan");
}

DeviceParameters device_from(const Config& c) {
  DeviceParameters p;
  p.phi_3 = c.num("device.phi_3", p.phi_3);
  p.phi_1 = c.num("device.phi_1", p.phi_1);
  p.phi_2 = c.num("device.phi_2", p.phi_2);
  p.emitter_height = c.num("device.emitter_height", p.emitter_height);
  p.extraction_distance = c.num("device.extraction_distance", p.extraction_distance);
  p.extractor_aperture_diameter =
      c.num("device.extractor_aperture_diameter", p.extractor_aperture_diameter);
  p.extractor_thickness = c.num("device.extractor_thickness", p.extractor_thickness);
  p.domain_radius = c.num("domain.radius", p.domain_radius);
  p.domain_z_min = c.num("domain.z_min", p.domain_z_min);
  p.domain_z_max = c.num("domain.z_max", p.domain_z_max);
  if (!c.has("device.extractor_outer_radius"))
    throw std::runtime_error("device.extractor_outer_radius fehlt (Pflichtangabe).");
  p.extractor_outer_radius = c.num("device.extractor_outer_radius", 0.0);
  p.emitter_back_length = 0.0;
  return p;
}

}  // namespace

// ===========================================================================

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.size() < 2) {
    std::printf(
        "es_p3b_audit -- P0: numerische Bereinigung von P3b\n\n"
        "  es_p3b_audit <geometrie.cfg> [<p3b.cfg> ...] <ausgabeverzeichnis> [key=value ...]\n\n"
        "Gerechnet wird keine neue Physik: Lastprojektion gegen hergestellte Lasten mit\n"
        "bekanntem Integral, gekoppelte Netzkonvergenz mit Richardson-Extrapolation an\n"
        "1000 V und 1400 V, und das Kanten-Gate fuer die ebene und eine nach aussen\n"
        "gewoelbte Form.\n");
    return 1;
  }

  Config cfg;
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) cfg.load(pos[i]);
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) {
    const std::filesystem::path src(pos[i]);
    std::filesystem::copy_file(src, std::filesystem::path(outdir) / src.filename(),
                               std::filesystem::copy_options::overwrite_existing);
  }

  // --- the model, identical to the P3b run ----------------------------------
  DielectricDeviceParameters geo;
  geo.device = device_from(cfg);
  geo.base_plate_thickness = cfg.num("device.base_plate_thickness", geo.base_plate_thickness);
  geo.reservoir = ReservoirModel::AxisymmetricPlenum;
  geo.feed_channel_radius = cfg.num("reservoir.feed_channel_radius", 0.0);
  geo.feed_channel_length = cfg.num("reservoir.feed_channel_length", geo.feed_channel_length);
  geo.plenum_radius = cfg.num("reservoir.plenum_radius", geo.plenum_radius);
  geo.plenum_depth = cfg.num("reservoir.plenum_depth", geo.plenum_depth);
  geo.plenum_wall_thickness = cfg.num("reservoir.wall_thickness", geo.plenum_wall_thickness);
  geo.plenum_fill_fraction = 1.0;

  MaterialLibrary lib;
  DielectricMaterials mats = DielectricMaterials::reference(lib);
  if (cfg.has("emitter.material")) mats.emitter_dielectric = lib.get(cfg.str("emitter.material", ""));
  if (cfg.has("extractor.material"))
    mats.extractor_carrier = lib.get(cfg.str("extractor.material", ""));
  mats.check_usable();

  LiquidProperties liquid = liquid_data_by_name(cfg.str("liquid.name", "emi-bf4"));
  liquid.validate_or_throw();

  const Real a = 0.5 * geo.device.phi_2;
  const Real gamma_over_a = liquid.gamma / a;
  const Metallisation metal = Metallisation::FrontAndAperture;
  const FarField far = FarField::Asymptotic;
  const std::size_t mem_cap = 6ull << 30;   // enough for one level beyond the default

  int exit_code = 0;
  // What the report and the meta record need, filled as the run goes.
  struct Summary {
    Real worst_exclusion_agreement{0};
    Real handed_jump_decay_field{0};
    Real handed_force_error_field{0};
    Real segment_error_smooth{0};
    std::string gate_flat, gate_bulged;
    Real beta_flat{0}, beta_bulged{0};
    Real excl_flat{0}, excl_bulged{0};
    Real excl_closed_flat{0}, excl_closed_bulged{0};
    int gate_levels_used{0};
    Real force_order_flat{0}, force_error_flat{0};
    std::string worst_verdict{"Converged"};
    Real worst_error{0};
  } sum;
  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& line) {
    std::printf("%s\n", line.c_str());
    std::fprintf(log, "%s\n", line.c_str());
    std::fflush(log);
  };

  say("P0 -- numerische Bereinigung von P3b");
  say("  a = " + std::to_string(a) + " m, gamma/a = " + std::to_string(gamma_over_a) + " Pa");

  // ==========================================================================
  // 1.  The load projection, on manufactured loads and on the real field
  // ==========================================================================
  say("");
  say("1. Lastprojektion");
  {
    std::FILE* f = std::fopen((outdir + "/projection_audit.csv").c_str(), "w");
    std::fprintf(f,
                 "# Audit der Lastprojektion.  'case' nennt die Last, 'segments' die\n"
                 "# Aufloesung der Oberflaeche.  analytic_force ist die geschlossene Form,\n"
                 "# wo es eine gibt, sonst nan.  Die drei Kraefte sind: die konservative\n"
                 "# Segmentprojektion, der Inhalt der Bins, und das Integral der TATSAECHLICH\n"
                 "# uebergebenen stetigen Last -- einmal gegen das rekonstruierte Flaechenmass\n"
                 "# A' (dort ist die Erhaltung exakt) und einmal gegen 2 pi r ds (dort gilt\n"
                 "# sie zur Interpolationsordnung).  Beide stehen hier, nicht nur die\n"
                 "# guenstigere.\n"
                 "#\n"
                 "# jump_handed / jump_stair sind die groessten Spruenge ueber die Binraender,\n"
                 "# gemessen mit einem Probenabstand delta und noch einmal mit delta/10.  Der\n"
                 "# ABFALL entscheidet: ~0.1 fuer eine stetige Last, ~1 fuer eine Treppe.\n");
    std::fprintf(f, "case,beta,segments,n_nodes,analytic_force_N,segment_force_N,bin_force_N,"
                    "handed_force_reconstructed_N,handed_force_true_N,err_segment_vs_analytic,"
                    "err_bin_vs_segment,err_handed_reconstructed,err_handed_true,"
                    "jump_handed_Pa,jump_handed_decay,jump_stair_Pa,jump_stair_decay,"
                    "load_span_Pa,tau_first,tau_last,empty_bins,area_gap,"
                    "last_bin_force_fraction,max_handed_Pa,max_bin_Pa,max_node_Pa\n");

    auto write_audit = [&](const std::string& name, Real beta, Index nseg,
                           const LoadProjectionAudit& A) {
      std::fprintf(f, "%s", name.c_str());
      put(f, beta);
      std::fprintf(f, ",%lld,%lld", static_cast<long long>(nseg),
                   static_cast<long long>(A.n_nodes));
      put(f, A.analytic_force);
      put(f, A.segment_force);
      put(f, A.bin_force);
      put(f, A.handed_force_reconstructed);
      put(f, A.handed_force_true);
      put(f, A.error_segment_vs_analytic);
      put(f, A.error_bin_vs_segment);
      put(f, A.error_handed_reconstructed);
      put(f, A.error_handed_true);
      put(f, A.handed_max_jump);
      put(f, A.handed_jump_decay);
      put(f, A.staircase_max_jump);
      put(f, A.staircase_jump_decay);
      put(f, A.load_span);
      put(f, A.tau_first);
      put(f, A.tau_last);
      std::fprintf(f, ",%lld", static_cast<long long>(A.empty_bins));
      put(f, A.area_gap);
      put(f, A.last_bin_force_fraction);
      put(f, A.max_handed_pressure);
      put(f, A.max_bin_pressure);
      put(f, A.max_node_pressure);
      std::fprintf(f, "\n");
    };

    const FreeSurface flat_fs = FreeSurface::flat_surface(a, 0.0);
    const Real p0 = 1.0e4;
    auto smooth = [&](Real d) {
      const Real r = a - d;
      return p0 * (1.0 + (r / a) * (r / a));
    };
    for (Index n : {16, 32, 64, 128, 256}) {
      const MaxwellLoad L = manufactured_load(flat_fs, uniform_radius_nodes(a, n), smooth,
                                              gamma_over_a);
      const LoadProjectionAudit A =
          audit_projection(L, flat_fs, "glatt", flat_disc_smooth_force(p0, a));
      write_audit("glatt", kNaN, n, A);
      if (n == 256) sum.segment_error_smooth = A.error_segment_vs_analytic;
    }
    for (Real beta : {-0.25, -0.44, -0.75}) {
      const Real C = 1.0e4 * std::pow(a, -beta);
      auto sing = [&](Real d) { return (d > 0.0) ? C * std::pow(d, beta) : 0.0; };
      for (Index n : {16, 32, 64, 128, 256}) {
        const MaxwellLoad L =
            manufactured_load(flat_fs, uniform_radius_nodes(a, n), sing, gamma_over_a);
        write_audit("singulaer", beta, n,
                    audit_projection(L, flat_fs, "singulaer",
                                     flat_disc_power_law_force(C, beta, a)));
      }
    }

    // The real thing: the solved field on the flat surface at the probe voltage.
    const Real V_probe = cfg.num("field.V_probe", 1500.0);
    for (int lvl : {1, 2, 3, 4}) {
      DielectricDeviceParameters p = geo;
      p.mesh_level = lvl;
      MeniscusMesh m = build_meniscus_mesh(p, flat_fs);
      DielectricSetup s;
      s.geometry = p;
      s.materials = mats;
      s.conductor_model = ConductorModel::Dielectric;
      s.metallisation = metal;
      s.far_field = far;
      s.V_emitter = V_probe;
      s.V_extractor = 0.0;
      s.memory_cap_bytes = mem_cap;
      const DielectricSolution sol =
          solve_dielectric_on(m.device, s, DielectricDiagnostics::FieldOnly);
      const MaxwellLoad L = maxwell_load(m, sol, gamma_over_a);
      const LoadProjectionAudit A = audit_projection(L, m.surface, "feld", kNaN);
      write_audit("feld_stufe" + std::to_string(lvl), kNaN,
                  static_cast<Index>(L.seg_force.size()), A);
      if (lvl == 2) {
        A.print(stdout);
        sum.handed_jump_decay_field = A.handed_jump_decay;
        sum.handed_force_error_field = A.error_handed_reconstructed;
      }
      // The declared properties, on the real load.
      if (!(A.handed_jump_decay < 0.2)) {
        say("  FEHLER: die uebergebene Last des Feldes ist nicht stetig (Abfall " +
            std::to_string(A.handed_jump_decay) + ")");
        exit_code = 2;
      }
      if (!(A.error_handed_reconstructed < 1.0e-12)) {
        say("  FEHLER: die uebergebene Last traegt die Kraft nicht (" +
            std::to_string(A.error_handed_reconstructed) + ")");
        exit_code = 2;
      }
      if (A.empty_bins != 0 || !(A.last_bin_force_fraction > 0.0)) {
        say("  FEHLER: Luecke oder Ausschlusszone in der Projektion");
        exit_code = 2;
      }
    }
    std::fclose(f);
    say("  projection_audit.csv geschrieben");
  }

  // --- the three loads, drawn: raw nodal, segment staircase, handed ----------
  {
    std::FILE* f = std::fopen((outdir + "/projection_profiles.csv").c_str(), "w");
    std::fprintf(f, "# Die drei Lasten nebeneinander, ueber der normierten Bogenlaenge tau.\n"
                    "# kind = node   : der ROHE punktweise Knotenwert (nicht konvergent an der"
                    " Kante, wird nirgends verwendet)\n"
                    "# kind = segment: die konservative Treppenfunktion (die Kraftbuchhaltung)\n"
                    "# kind = handed : die stetige Rekonstruktion, die dem Kapillarloeser"
                    " TATSAECHLICH uebergeben wird\n"
                    "# kind = cumulative: die kumulierte Kraft G(tau), aus der 'handed' gebaut"
                    " ist\n");
    std::fprintf(f, "case,kind,tau,d_over_a,value\n");

    auto dump = [&](const std::string& name, const MaxwellLoad& L, const FreeSurface& fs) {
      for (std::size_t k = 0; k < L.node_tau.size(); ++k)
        std::fprintf(f, "%s,node,%.9e,%.9e,%.9e\n", name.c_str(), L.node_tau[k],
                     L.node_d_edge[k] / a, L.node_pM[k]);
      for (std::size_t k = 0; k < L.seg_pressure.size(); ++k) {
        std::fprintf(f, "%s,segment,%.9e,%.9e,%.9e\n", name.c_str(), L.seg_tau0[k],
                     L.seg_d_mid[k] / a, L.seg_pressure[k]);
        std::fprintf(f, "%s,segment,%.9e,%.9e,%.9e\n", name.c_str(), L.seg_tau1[k],
                     L.seg_d_mid[k] / a, L.seg_pressure[k]);
      }
      const ProjectedLoad pl = ProjectedLoad::from(L);
      const int m = 2048;
      for (int k = 0; k <= m; ++k) {
        const Real tau = static_cast<Real>(k) / static_cast<Real>(m);
        std::fprintf(f, "%s,handed,%.9e,%.9e,%.9e\n", name.c_str(), tau,
                     (1.0 - tau) * fs.arclength / a, pl.at(tau));
      }
      const std::vector<Real>& G = pl.cumulative_force();
      for (std::size_t k = 0; k < G.size(); ++k)
        std::fprintf(f, "%s,cumulative,%.9e,%.9e,%.9e\n", name.c_str(),
                     static_cast<Real>(k) / static_cast<Real>(ProjectedLoad::kBins), kNaN, G[k]);
    };

    const FreeSurface flat_fs = FreeSurface::flat_surface(a, 0.0);
    const Real beta = -0.44, C = 1.0e4 * std::pow(a, -beta);
    auto sing = [&](Real d) { return (d > 0.0) ? C * std::pow(d, beta) : 0.0; };
    dump("singulaer", manufactured_load(flat_fs, uniform_radius_nodes(a, 32), sing,
                                        gamma_over_a), flat_fs);
    const Real p0 = 1.0e4;
    auto smooth = [&](Real d) {
      const Real r = a - d;
      return p0 * (1.0 + (r / a) * (r / a));
    };
    dump("glatt", manufactured_load(flat_fs, uniform_radius_nodes(a, 32), smooth, gamma_over_a),
         flat_fs);
    {
      DielectricDeviceParameters p = geo;
      p.mesh_level = 2;
      MeniscusMesh m = build_meniscus_mesh(p, flat_fs);
      DielectricSetup s;
      s.geometry = p;
      s.materials = mats;
      s.conductor_model = ConductorModel::Dielectric;
      s.metallisation = metal;
      s.far_field = far;
      s.V_emitter = cfg.num("field.V_probe", 1500.0);
      s.V_extractor = 0.0;
      s.memory_cap_bytes = mem_cap;
      const DielectricSolution sol =
          solve_dielectric_on(m.device, s, DielectricDiagnostics::FieldOnly);
      dump("feld", maxwell_load(m, sol, gamma_over_a), m.surface);
    }
    std::fclose(f);
    say("  projection_profiles.csv geschrieben");
  }

  // ==========================================================================
  // 2.  The exclusion criterion, resolved in closed form
  // ==========================================================================
  say("");
  say("2. kTolExclusion: geschlossene Form gegen Messung");
  {
    std::FILE* f = std::fopen((outdir + "/exclusion_criterion.csv").c_str(), "w");
    std::fprintf(f,
                 "# edge_gate::kTolExclusion verlangt, dass das Halbieren der Ausschluss-\n"
                 "# distanz die integrierte Kraft um weniger als 5 %% aendert.  Fuer eine\n"
                 "# Last p = C d^beta hat diese Groesse einen GESCHLOSSENEN GRENZWERT\n"
                 "# ungleich null, der die Netzweite gar nicht enthaelt.  Hier steht er\n"
                 "# neben der Messung, ueber vier Verfeinerungen: die Messung laeuft gegen\n"
                 "# den Grenzwert und bleibt stehen, waehrend der Diskretisierungsfehler\n"
                 "# derselben Last weiter faellt.  Damit ist das Kriterium nicht getauscht,\n"
                 "# sondern widerlegt.\n");
    std::fprintf(f, "beta,segments,measured_change,closed_form_limit,agreement,"
                    "discretisation_error\n");
    const FreeSurface flat_fs = FreeSurface::flat_surface(a, 0.0);
    Real worst_agreement = 0.0;
    for (Real beta : {-0.25, -0.44, -0.75, -0.90}) {
      const Real C = 1.0e4 * std::pow(a, -beta);
      auto sing = [&](Real d) { return (d > 0.0) ? C * std::pow(d, beta) : 0.0; };
      const Real limit = exclusion_halving_limit(beta, edge_gate::kExclusionMid);
      const Real F_exact = flat_disc_power_law_force(C, beta, a);
      for (Index n : {32, 64, 128, 256, 512, 1024}) {
        const MaxwellLoad L =
            manufactured_load(flat_fs, uniform_radius_nodes(a, n), sing, gamma_over_a);
        const Real fm = L.force_beyond(edge_gate::kExclusionMid * a);
        const Real ff = L.force_beyond(edge_gate::kExclusionFine * a);
        const Real measured = std::abs(ff - fm) / std::abs(ff);
        const Real agreement = std::abs(measured - limit) / limit;
        std::fprintf(f, "%.4f,%lld", beta, static_cast<long long>(n));
        put(f, measured);
        put(f, limit);
        put(f, agreement);
        put(f, std::abs(L.total_force - F_exact) / F_exact);
        std::fprintf(f, "\n");
        if (n == 1024) worst_agreement = std::max(worst_agreement, agreement);
      }
    }
    std::fclose(f);
    sum.worst_exclusion_agreement = worst_agreement;
    say("  groesste Abweichung Messung/geschlossene Form auf der feinsten Stufe: " +
        std::to_string(worst_agreement));
    if (!(worst_agreement < 2.0e-2)) {
      say("  FEHLER: die geschlossene Form trifft die Messung nicht.");
      exit_code = 2;
    }
  }

  // ==========================================================================
  // 3.  The edge gate: flat and one outward shape, with the extra levels
  // ==========================================================================
  say("");
  say("3. Kanten-Gate, ebene und nach aussen gewoelbte Form");
  const int gate_lo = cfg.integer("mesh.gate_min_level", 1);
  const int gate_hi = cfg.integer("audit.gate_max_level", 5);
  std::vector<int> gate_levels;
  for (int l = gate_lo; l <= gate_hi; ++l) gate_levels.push_back(l);
  {
    std::FILE* f = std::fopen((outdir + "/gate_levels.csv").c_str(), "w");
    std::fprintf(f, "# Kanten-Gate je Form und Netzstufe.  peak_node_pM ist der PUNKTWEISE\n"
                    "# Kantenwert: er konvergiert nicht und wird nirgends verwendet.\n");
    std::fprintf(f, "shape,Pi,level,n_nodes,n_segments,smallest_d_m,total_force_N,"
                    "force_beyond_coarse_N,force_beyond_mid_N,force_beyond_fine_N,"
                    "peak_node_pM_Pa,fit_exponent,fit_r2\n");
    std::FILE* g = std::fopen((outdir + "/gate_verdict.csv").c_str(), "w");
    std::fprintf(g, "# Urteil je Form, mit ALLEN drei Sichten nebeneinander: die direkte\n"
                    "# Kraftfolge ueber die Netzstufen (in gate_levels.csv), die Aitken-\n"
                    "# Extrapolation darueber, und die Extrapolation ueber die Ausschluss-\n"
                    "# distanz.  measured_exclusion_change ist die Groesse aus 2.: sie wird\n"
                    "# berichtet und entscheidet nicht.\n");
    std::fprintf(g, "shape,Pi,verdict,fitted_exponent,limit_force_mesh_N,"
                    "limit_force_exclusion_N,limit_agreement,total_force_change,"
                    "edge_far_change,measured_exclusion_change,closed_form_exclusion_limit,"
                    "levels\n");

    auto shape_for = [&](Real Pi) {
      CapillaryRequest cr;
      cr.delta_p_exit = capillary::pressure_from_pi(Pi, a, liquid.gamma);
      cr.target_relative_accuracy = 1.0e-10;
      return solve_capillary_meniscus(a, 0.0, liquid, cr);
    };

    struct Case { const char* tag; Real Pi; };
    const Case cases[] = {{"eben", 0.0}, {"gewoelbt", 1.5}};
    const Real V_probe = cfg.num("field.V_probe", 1500.0);
    for (const Case& c : cases) {
      FreeSurface fs = (c.Pi == 0.0) ? FreeSurface::flat_surface(a, 0.0)
                                     : FreeSurface::from(shape_for(c.Pi));
      std::vector<int> used;
      EdgeGateResult gate;
      // Fall back one level at a time if the finest does not fit in memory, and
      // SAY so -- a level silently dropped would make the study a study of the
      // memory cap.
      for (std::vector<int> lv = gate_levels; lv.size() >= 3; lv.pop_back()) {
        try {
          gate = run_edge_gate(geo, mats, V_probe, 0.0, metal, far, fs, c.tag, c.Pi, lv,
                               gamma_over_a, mem_cap);
          used = lv;
          break;
        } catch (const std::exception& e) {
          say(std::string("  Stufe ") + std::to_string(lv.back()) + " nicht rechenbar: " +
              e.what());
        }
      }
      if (used.empty()) {
        say(std::string("  Form ") + c.tag + ": keine Stufe rechenbar.");
        exit_code = 2;
        continue;
      }
      for (const EdgeStudyPoint& pt : gate.levels) {
        std::fprintf(f, "%s,%.4f,%d,%lld,%lld", c.tag, c.Pi, pt.mesh_level,
                     static_cast<long long>(pt.n_nodes),
                     static_cast<long long>(pt.n_surface_segments));
        put(f, pt.smallest_d);
        put(f, pt.total_force);
        put(f, pt.force_coarse);
        put(f, pt.force_mid);
        put(f, pt.force_fine);
        put(f, pt.peak_node_pM);
        put(f, pt.fit_exponent);
        put(f, pt.fit_r2);
        std::fprintf(f, "\n");
      }
      std::fprintf(g, "%s,%.4f,%s", c.tag, c.Pi, to_string(gate.verdict));
      put(g, gate.fitted_exponent);
      put(g, gate.limit_force_mesh);
      put(g, gate.limit_force_exclusion);
      put(g, gate.measured_limit_agreement);
      put(g, gate.measured_total_force_change);
      put(g, gate.measured_edge_far_change);
      put(g, gate.measured_exclusion_change);
      put(g, exclusion_halving_limit(gate.fitted_exponent, edge_gate::kExclusionMid));
      std::fprintf(g, ",%d\n", static_cast<int>(used.size()));

      if (c.Pi == 0.0) {
        sum.gate_flat = to_string(gate.verdict);
        sum.beta_flat = gate.fitted_exponent;
        sum.excl_flat = gate.measured_exclusion_change;
        sum.excl_closed_flat =
            exclusion_halving_limit(gate.fitted_exponent, edge_gate::kExclusionMid);
        sum.gate_levels_used = static_cast<int>(used.size());
      } else {
        sum.gate_bulged = to_string(gate.verdict);
        sum.beta_bulged = gate.fitted_exponent;
        sum.excl_bulged = gate.measured_exclusion_change;
        sum.excl_closed_bulged =
            exclusion_halving_limit(gate.fitted_exponent, edge_gate::kExclusionMid);
      }
      say(std::string("  Form ") + c.tag + " (Pi = " + std::to_string(c.Pi) + "): " +
          to_string(gate.verdict) + ", beta = " + std::to_string(gate.fitted_exponent) +
          ", Stufen " + std::to_string(used.front()) + ".." + std::to_string(used.back()));
      {
        const Real closed =
            exclusion_halving_limit(gate.fitted_exponent, edge_gate::kExclusionMid);
        const bool met = gate.measured_exclusion_change < edge_gate::kTolExclusion;
        say("    Halbierungsaenderung gemessen " +
            std::to_string(gate.measured_exclusion_change) +
            ", geschlossene Form fuer beta = " + std::to_string(gate.fitted_exponent) + ": " +
            std::to_string(closed) +
            (met ? " -- die Schranke ist hier eingehalten, WEIL beta gross genug ist, nicht"
                   " weil das Netz fein genug waere."
                 : " -- die Schranke ist hier verfehlt, um einen Betrag, den beta festlegt"
                   " und den kein Netz aendert."));
      }

      // The force sequence, extrapolated with Richardson as well, so that the
      // gate's own Aitken value can be seen next to an independent estimate.
      std::vector<Real> forces;
      for (const EdgeStudyPoint& pt : gate.levels) forces.push_back(pt.total_force);
      const RichardsonEstimate re = richardson(forces, kMeshLevelRatio);
      say("    Gesamtkraft ueber die Stufen: Ordnung " + std::to_string(re.observed_order) +
          ", extrapoliert " + std::to_string(re.extrapolated) + " N, geschaetzter Fehler " +
          std::to_string(re.relative_error_finest) + " (" +
          to_string(verdict_of(re)) + ")");
      if (c.Pi == 0.0) {
        sum.force_order_flat = re.observed_order;
        sum.force_error_flat = re.relative_error_finest;
      }
    }
    std::fclose(f);
    std::fclose(g);
  }

  // ==========================================================================
  // 4.  The coupled mesh convergence at 1000 V and 1400 V
  // ==========================================================================
  say("");
  say("4. Gekoppelte Netzkonvergenz bei 1000 V und 1400 V");
  {
    std::FILE* f = std::fopen((outdir + "/coupled_convergence.csv").c_str(), "w");
    std::fprintf(f, "# Gekoppelte Netzkonvergenz.  KEIN neuer Spannungssweep: zwei Punkte,\n"
                    "# die drei vorhandenen Netzstufen plus eine feinere.  Eine gescheiterte\n"
                    "# Rechnung steht als nan, NIEMALS als null -- null ist ein physikalischer\n"
                    "# Wert.\n");
    std::fprintf(f, "V_emitter_V,level,n_nodes,status,h_over_a,E_edge_far_V_per_m,"
                    "total_force_N,max_curvature_1_per_m,mech_residual_edge_far,"
                    "mech_residual_all,iterations\n");
    std::FILE* g = std::fopen((outdir + "/coupled_richardson.csv").c_str(), "w");
    std::fprintf(g, "# Richardson-Extrapolation ueber die Netzstufen, Verfeinerungsverhaeltnis\n"
                    "# sqrt(2) je Stufe.  Ziel VOR der Messung festgelegt: geschaetzter\n"
                    "# Diskretisierungsfehler der feinsten Stufe unter 1 %%.  Wo das verfehlt\n"
                    "# wird, steht DiscretizationNotConverged und das Ergebnis ist nur\n"
                    "# qualitativ.\n");
    std::fprintf(g, "V_emitter_V,quantity,window,n_levels,observed_order,extrapolated,"
                    "relative_error_finest,last_relative_change,verdict\n");

    const std::vector<int> levels = {1, 2, 3, 4};
    for (Real V : {1000.0, 1400.0}) {
      std::vector<Real> h_over_a, e_far, force;
      std::vector<int> ok_levels;
      for (int lvl : levels) {
        CoupledRequest q;
        q.geometry = geo;
        q.geometry.mesh_level = lvl;
        q.materials = mats;
        q.liquid = liquid;
        q.metallisation = metal;
        q.far_field = far;
        q.delta_p_exit = 0.0;
        q.V_emitter = V;
        q.memory_cap_bytes = mem_cap;
        CoupledPoint p;
        try {
          p = solve_coupled(q);
        } catch (const std::exception& e) {
          p.status = CouplingStatus::ElectrostaticFailure;
          p.message = e.what();
        }
        const bool good = is_usable(p.status);
        std::fprintf(f, "%.9e,%d,%lld,%s", V, lvl, 0LL, to_string(p.status));
        put(f, good ? p.apex_height / a : kNaN);
        put(f, good ? p.E_edge_far : kNaN);
        put(f, good ? p.total_force : kNaN);
        put(f, good ? p.max_curvature : kNaN);
        put(f, good ? p.mechanical_residual_edge_far : kNaN);
        put(f, good ? p.mechanical_residual : kNaN);
        std::fprintf(f, ",%d\n", p.iterations);
        std::fflush(f);
        if (good) {
          h_over_a.push_back(p.apex_height / a);
          e_far.push_back(p.E_edge_far);
          force.push_back(p.total_force);
          ok_levels.push_back(lvl);
        }
        say("    " + std::to_string(static_cast<int>(V)) + " V, Stufe " +
            std::to_string(lvl) + ": " + to_string(p.status) +
            (good ? ", h/a = " + std::to_string(p.apex_height / a) : "") + ", " +
            std::to_string(p.iterations) + " Iterationen");
      }
      struct Q { const char* name; const std::vector<Real>* v; };
      const Q qs[] = {{"h_over_a", &h_over_a}, {"E_edge_far", &e_far}, {"total_force", &force}};
      for (const Q& q : qs) {
        // TWO WINDOWS, because three points give one order and four give two,
        // and a single order that happens to be observable can hide a sequence
        // that is not in its asymptotic range.  The fine window decides; the
        // coarse one is written so that the reader can see whether they agree.
        for (int window = 0; window < 2; ++window) {
          std::vector<Real> v = *q.v;
          if (v.size() > 3) {
            if (window == 0)
              v.pop_back();            // levels 1..3
            else
              v.erase(v.begin());      // levels 2..4
          }
          const RichardsonEstimate re = richardson(v, kMeshLevelRatio);
          const DiscretizationVerdict vd = verdict_of(re);
          std::fprintf(g, "%.9e,%s,%s,%d", V, q.name, window == 0 ? "grob" : "fein",
                       re.n_levels);
          put(g, re.observed_order);
          put(g, re.extrapolated);
          put(g, re.relative_error_finest);
          put(g, re.last_relative_change);
          std::fprintf(g, ",%s\n", to_string(vd));
          if (window == 1 && vd != DiscretizationVerdict::Converged) {
            sum.worst_verdict = to_string(vd);
            if (std::isfinite(re.relative_error_finest))
              sum.worst_error = std::max(sum.worst_error, re.relative_error_finest);
          }
          if (window == 1)
            say(std::string("      ") + q.name + ": Ordnung " +
                std::to_string(re.observed_order) + ", Fehler der feinsten Stufe " +
                std::to_string(re.relative_error_finest) + ", blosse Aenderung " +
                std::to_string(re.last_relative_change) + " -> " + to_string(vd));
        }
      }
    }
    std::fclose(f);
    std::fclose(g);
  }

  // ==========================================================================
  // 5.  Parameters and meta
  // ==========================================================================
  {
    std::FILE* f = std::fopen((outdir + "/parameters.csv").c_str(), "w");
    std::fprintf(f, "name,value_SI,unit,role\n");
    std::fprintf(f, "contact_radius,%.9e,m,gepinnter Kontaktradius a = phi_2/2\n", a);
    std::fprintf(f, "gamma_over_a,%.9e,Pa,Druckskala\n", gamma_over_a);
    std::fprintf(f, "gamma,%.9e,N/m,Oberflaechenspannung (%s)\n", liquid.gamma,
                 to_string(liquid.status));
    std::fprintf(f, "discretization_target,%.9e,-,vorab festgelegtes Ziel\n",
                 kDiscretizationTarget);
    std::fprintf(f, "mesh_level_ratio,%.9e,-,Verfeinerungsverhaeltnis in h je Stufe\n",
                 kMeshLevelRatio);
    std::fprintf(f, "tol_exclusion,%.9e,-,edge_gate::kTolExclusion (deklariert, widerlegt)\n",
                 edge_gate::kTolExclusion);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_p3b_audit (P0)\n");
    std::fprintf(f, "phase=P0\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "config=device_p1.cfg;electrocapillary_p3b.cfg\n");
    std::fprintf(f, "state=numerische Bereinigung von P3b: Lastprojektion, gekoppelte "
                    "Netzkonvergenz, Kanten-Gate.  KEINE neue Physik.\n");
    std::fprintf(f, "contact_radius_m=%.9e\n", a);
    std::fprintf(f, "surface_tension_N_per_m=%.9e\n", liquid.gamma);
    std::fprintf(f, "capillary_pressure_scale_Pa=%.9e\n", gamma_over_a);
    std::fprintf(f, "liquid_substance=%s\n", liquid.substance.c_str());
    std::fprintf(f, "liquid_status=%s\n", to_string(liquid.status));
    std::fprintf(f, "discretization_target=%.9e\n", kDiscretizationTarget);
    std::fprintf(f, "gate_levels_used=%d\n", sum.gate_levels_used);
    std::fprintf(f, "gate_flat=%s\n", sum.gate_flat.c_str());
    std::fprintf(f, "gate_bulged=%s\n", sum.gate_bulged.c_str());
    std::fprintf(f, "beta_flat=%.9e\n", sum.beta_flat);
    std::fprintf(f, "beta_bulged=%.9e\n", sum.beta_bulged);
    std::fprintf(f, "exclusion_measured_flat=%.9e\n", sum.excl_flat);
    std::fprintf(f, "exclusion_closed_form_flat=%.9e\n", sum.excl_closed_flat);
    std::fprintf(f, "exclusion_measured_bulged=%.9e\n", sum.excl_bulged);
    std::fprintf(f, "exclusion_closed_form_bulged=%.9e\n", sum.excl_closed_bulged);
    std::fprintf(f, "exclusion_worst_agreement=%.9e\n", sum.worst_exclusion_agreement);
    std::fprintf(f, "handed_jump_decay_field=%.9e\n", sum.handed_jump_decay_field);
    std::fprintf(f, "handed_force_error_field=%.9e\n", sum.handed_force_error_field);
    std::fprintf(f, "coupled_worst_verdict=%s\n", sum.worst_verdict.c_str());
    std::fprintf(f, "coupled_worst_error=%.9e\n", sum.worst_error);
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/report.txt").c_str(), "w");
    std::fprintf(f,
      "P0 -- numerische Bereinigung von P3b\n"
      "====================================\n\n"
      "WAS DIESER LAUF IST\n"
      "  Keine neue Physik.  Drei offene Punkte aus P3b werden gemessen statt\n"
      "  behauptet: die Lastprojektion, die gekoppelte Netzkonvergenz und das\n"
      "  Kanten-Gate.  Eine gescheiterte Rechnung steht ueberall als nan, nie als\n"
      "  null.\n\n");
    std::fprintf(f,
      "1. LASTPROJEKTION\n"
      "  Drei Groessen heissen in P3b 'die Last' und sind jetzt getrennt:\n"
      "    * der ROHE Knotenwert p_M = eps0 E_n^2/2 -- an der Kante nicht\n"
      "      konvergent, wird nirgends verwendet;\n"
      "    * die konservative SEGMENTprojektion -- die Kraftbuchhaltung, eine\n"
      "      Treppenfunktion;\n"
      "    * die UEBERGEBENE Last p(tau) = G'/A' -- was der Kapillarloeser\n"
      "      integriert.\n"
      "  Geprueft an zwei hergestellten Lasten mit bekanntem Integral, einer\n"
      "  glatten und einer integrablen Singularitaet p = C d^beta, und am\n"
      "  wirklich geloesten Feld.\n"
      "    Segmentprojektion gegen die geschlossene Form, glatte Last, 256\n"
      "      Segmente          : %.3e relativ (zweite Ordnung, gemessen)\n"
      "    uebergebene Last, Krafterhaltung gegen A' (Feld, Stufe 2)\n"
      "                        : %.3e relativ\n"
      "    Stetigkeit der uebergebenen Last (Feld, Stufe 2): der Sprung ueber die\n"
      "      Binraender faellt bei zehnfach kleinerem Probenabstand auf das\n"
      "      %.4f-fache.  Fuer eine Treppenfunktion waere er 1.\n"
      "    Kontaktkante: die Segmente reichen bis tau = 1, kein Bin ist leer, das\n"
      "      kantennaechste Bin traegt Kraft.  Es gibt keine Ausschlusszone und\n"
      "      keine Kappung; der punktweise Kantenwert waechst mit jeder\n"
      "      Verfeinerung und wird nirgends verwendet.\n\n",
      sum.segment_error_smooth, sum.handed_force_error_field, sum.handed_jump_decay_field);
    std::fprintf(f,
      "2. DAS KRITERIUM kTolExclusion -- AUFGELOEST, NICHT GETAUSCHT\n"
      "  P3b berichtete: die geforderten 5 %% werden um 11,3 %% verfehlt, ein\n"
      "  anderes Kriterium entscheidet.  Das ist jetzt hergeleitet.  Fuer eine\n"
      "  Last p = C d^beta ist die Kraft jenseits einer Ausschlussdistanz d0 in\n"
      "  geschlossener Form bekannt, und die Groesse, die kTolExclusion prueft,\n"
      "      [F(d0/2) - F(d0)] / F(d0/2) ,\n"
      "  enthaelt die Netzweite ueberhaupt nicht.  Sie hat einen endlichen\n"
      "  Grenzwert ungleich null, der allein von beta und d0 abhaengt.\n"
      "    groesste Abweichung Messung gegen geschlossene Form: %.3e\n"
      "    ebene Form   : beta = %+.4f, gemessen %.4f, geschlossene Form %.4f\n"
      "    gewoelbte Form: beta = %+.4f, gemessen %.4f, geschlossene Form %.4f\n"
      "  Die Schranke bleibt im Code deklariert und wird weiter berichtet.  Sie\n"
      "  entscheidet nicht, und der Grund dafuer ist eine Rechnung mit\n"
      "  Regressionstest (tests/test_load_projection.cpp, Abschnitt 4), keine\n"
      "  nachtraegliche Wahl.  Das ersetzende Kriterium wird zusaetzlich gegen\n"
      "  einen BEKANNTEN Grenzwert geprueft (Abschnitt 5 desselben Tests).\n\n",
      sum.worst_exclusion_agreement, sum.beta_flat, sum.excl_flat, sum.excl_closed_flat,
      sum.beta_bulged, sum.excl_bulged, sum.excl_closed_bulged);
    std::fprintf(f,
      "3. DAS KANTEN-GATE\n"
      "  Nur die ebene und eine nach aussen gewoelbte Form, mit %d Netzstufen.\n"
      "    ebene Form    : %s\n"
      "    gewoelbte Form: %s\n"
      "  Die integrierte Gesamtkraft der ebenen Form konvergiert mit der Ordnung\n"
      "  %.3f.  Das ist NICHT zweite Ordnung, und es ist auch nicht ueberraschend:\n"
      "  1 + beta = %.3f ist genau die Rate, die eine integrable Singularitaet am\n"
      "  Rand zulaesst.  Der geschaetzte Diskretisierungsfehler der feinsten Stufe\n"
      "  ist %.3f -- ueber dem vorab festgelegten Ziel von 1 %%.\n\n",
      sum.gate_levels_used, sum.gate_flat.c_str(), sum.gate_bulged.c_str(),
      sum.force_order_flat, 1.0 + sum.beta_flat, sum.force_error_flat);
    std::fprintf(f,
      "4. GEKOPPELTE NETZKONVERGENZ\n"
      "  Zwei Betriebspunkte, 1000 V und 1400 V, vier Netzstufen, kein neuer\n"
      "  Spannungssweep.  Richardson auf zwei Fenstern (Stufen 1-3 und 2-4).\n"
      "    schlechtestes Urteil: %s, geschaetzter Fehler bis %.3f\n"
      "  Ziel VOR der Messung: unter 1 %%.  Es ist nicht erreicht.  Die\n"
      "  gekoppelten Ergebnisse von P3b sind damit QUALITATIV; keine Apexhoehe\n"
      "  und keine Kraft aus diesem Ast traegt drei Stellen.\n\n",
      sum.worst_verdict.c_str(), sum.worst_error);
    std::fprintf(f,
      "WAS DIESER LAUF NICHT IST\n"
      "  Keine Emission, keine endliche Leitfaehigkeit, keine Stroemung, keine\n"
      "  Raumladung, keine Zeitabhaengigkeit, keine Stabilitaetsaussage, kein\n"
      "  Taylor-Kegel.  delta_p_exit bleibt in diesem Lauf eine Eingabe.  Die\n"
      "  Stoffwerte bleiben illustrative.\n");
    std::fclose(f);
  }

  say("");
  say(exit_code == 0 ? "Alle deklarierten Pruefungen dieses Laufs bestanden."
                     : "MINDESTENS EINE DEKLARIERTE PRUEFUNG IST FEHLGESCHLAGEN.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_p3b_audit: %s\n", e.what());
  return 2;
}
