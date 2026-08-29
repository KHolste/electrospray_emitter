// es_dielectric -- P2b: dielectric axisymmetric electrostatics of the capillary
// Kunze emitter.
//
//   es_dielectric <geometry.cfg> [<p2b.cfg> ...] <output-directory> [key=value ...]
//
// Solves  div(eps(x) grad phi) = 0  on the axisymmetric device, with
//
//   * the ionic liquid as an ideal equipotential conductor at V_emitter,
//     terminated at a named feed boundary whose position is a parameter;
//   * the 3D-printed emitter body as a DIELECTRIC (SU-8 by default) -- not an
//     electrode, which is the correction this phase exists for;
//   * the extractor as a polymer carrier with a metallised face at V_extractor;
//   * an asymptotic (monopole) condition on the open far boundary.
//
// NOT in this phase: space charge, emission, meniscus motion, flow, finite
// liquid conductivity, time dependence.  The plane at z = 0 is the initial flat
// liquid surface, not a computed meniscus.
//
// Exit code 2 means a check failed; 3 means an unimplemented option was asked
// for.  Anything the run cannot back up, it refuses to print.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/boundary_mesh.hpp"
#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/dielectric_device.hpp"
#include "es/io.hpp"
#include "es/vacuum_bem.hpp"

using namespace es;
using constants::eps0;

namespace {

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
  // P2a's conducting rear closure has no place in a dielectric model; refusing
  // it here is what stops the superseded arrangement from creeping back.
  if (c.has("device.emitter_back_length") && c.num("device.emitter_back_length", 0.0) != 0.0)
    throw std::runtime_error(
        "device.emitter_back_length ist gesetzt.  Das ist der metallische P2a-Emitter mit "
        "leitender Abschlussscheibe und fuer den dielektrischen Emitter physikalisch falsch. "
        "In P2b wird die Fluessigkeitssaeule an feed.liquid_feed_z abgeschnitten.");
  p.emitter_back_length = 0.0;
  p.reserved.edge_radius_inner = c.num("reserved.edge_radius_inner", 0.0);
  p.reserved.edge_radius_outer = c.num("reserved.edge_radius_outer", 0.0);
  p.reserved.contact_angle_deg = c.num("reserved.contact_angle_deg", 0.0);
  p.reserved.bore_diameter_at_inlet = c.num("reserved.bore_diameter_at_inlet", 0.0);
  p.reserved.porous_emitter = c.flag("reserved.porous_emitter", false);
  p.reserved.collector_enabled = c.flag("reserved.collector_enabled", false);
  return p;
}

MaterialStatus status_from(const std::string& s) {
  if (s == "measured") return MaterialStatus::Measured;
  if (s == "manufacturer_spec") return MaterialStatus::ManufacturerSpec;
  if (s == "literature") return MaterialStatus::Literature;
  if (s == "provisional") return MaterialStatus::Provisional;
  throw std::runtime_error(
      "material.<name>.status muss measured, manufacturer_spec, literature oder provisional "
      "sein.  Ein Wert ohne Herkunftsangabe wird nicht angenommen; angegeben war: '" + s + "'");
}

/// Apply every material.<name>.* override.  This is the whole mechanism by
/// which IP-Q, IPx-Q or a better SU-8 measurement enters -- no code change.
void apply_material_overrides(const Config& c, MaterialLibrary& lib) {
  for (const char* raw : {"su8", "ip-q", "ipx-q"}) {
    const std::string name(raw);
    const std::string key = "material." + name + ".relative_permittivity";
    if (!c.has(key)) continue;
    const std::string src = c.str("material." + name + ".source", "");
    if (src.empty())
      throw std::runtime_error(
          "material." + name +
          ".relative_permittivity wurde gesetzt, aber material." + name +
          ".source fehlt.  Ein Materialwert ohne Herkunft ist in diesem Projekt keine "
          "Eingabe, sondern eine Behauptung.");
    lib.override_permittivity(name, c.num(key, 0.0),
                              status_from(c.str("material." + name + ".status", "provisional")),
                              src);
  }
}

Metallisation metallisation_from(const std::string& s) {
  if (s == "front_only") return Metallisation::FrontOnly;
  if (s == "front_and_aperture") return Metallisation::FrontAndAperture;
  if (s == "all_surfaces") return Metallisation::AllSurfaces;
  throw std::runtime_error("extractor.metallisation muss front_only, front_and_aperture oder "
                           "all_surfaces sein, war: '" + s + "'");
}

// ---------------------------------------------------------------------------

struct Window {
  std::string name;
  Real r0, r1, z0, z1;
  Index nr, nz;
};

/// Sample phi, |E| and the material region on a uniform window.  Written for
/// plotting only; nothing is read back into the physics.
void write_grid_csv(const std::string& path, const DielectricSolution& s, const Window& w) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# Abtastung fuer die Abbildungen: %s\n", w.name.c_str());
  std::fprintf(f, "# |E| ist im Inneren der Fluessigkeit nicht definiert (idealer Leiter) "
                  "und dort 0\n");
  std::fprintf(f, "r_m,z_m,phi_V,E_V_per_m,region\n");
  for (Index j = 0; j < w.nz; ++j) {
    const Real z = w.z0 + (w.z1 - w.z0) * static_cast<Real>(j) / (w.nz - 1);
    for (Index i = 0; i < w.nr; ++i) {
      const Real r = w.r0 + (w.r1 - w.r0) * static_cast<Real>(i) / (w.nr - 1);
      const Vec2 x{r, z};
      Real phi = 0.0, E = 0.0;
      try {
        phi = potential_at(s.mesh.grid, s.fem.phi, x);
      } catch (const std::exception&) {
        continue;
      }
      try {
        E = norm(field_recovered_at(s.mesh.grid, s.fem.phi, s.cell_eps_r, s.cell_active, x));
      } catch (const std::exception&) {
        E = 0.0;
      }
      std::fprintf(f, "%.7e,%.7e,%.7e,%.7e,%s\n", r, z, phi, E,
                   to_string(s.mesh.region_at(x)));
    }
  }
  std::fclose(f);
}

Real rel_change(Real a, Real b) {
  const Real d = std::max(std::abs(a), std::abs(b));
  return d > 0.0 ? std::abs(a - b) / d : 0.0;
}

Index probe_index(const std::vector<Probe>& p, const char* name) {
  for (std::size_t k = 0; k < p.size(); ++k)
    if (p[k].name == name) return static_cast<Index>(k);
  return -1;
}

}  // namespace

// ---------------------------------------------------------------------------

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.size() < 2) {
    std::fprintf(stderr,
                 "es_dielectric <geometrie.cfg> [<p2b.cfg> ...] <ausgabeverzeichnis> "
                 "[key=value ...]\n");
    return 1;
  }
  Config cfg;
  for (std::size_t k = 0; k + 1 < pos.size(); ++k) cfg.load(pos[k]);
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  // --- setup ---------------------------------------------------------------
  MaterialLibrary lib;
  apply_material_overrides(cfg, lib);

  DielectricSetup S;
  S.geometry.device = device_from(cfg);
  S.geometry.liquid_feed_z = cfg.num("feed.liquid_feed_z", -2.0e-4);
  const int ref_level = cfg.integer("mesh.reference_level", 3);
  const int max_level = cfg.integer("mesh.max_level", 4);
  S.geometry.mesh_level = ref_level;
  S.materials = DielectricMaterials::reference(lib);
  S.materials.emitter_dielectric = lib.get(cfg.str("emitter.material", "su8"));
  S.materials.extractor_carrier = lib.get(cfg.str("extractor.material", "su8"));
  S.metallisation = metallisation_from(cfg.str("extractor.metallisation", "front_and_aperture"));
  S.far_field = FarField::Asymptotic;
  if (!cfg.has("field.V_emitter") || !cfg.has("field.V_extractor"))
    throw std::runtime_error("field.V_emitter und field.V_extractor sind Pflichtangaben.");
  S.V_emitter = cfg.num("field.V_emitter", 0.0);
  S.V_extractor = cfg.num("field.V_extractor", 0.0);

  std::printf("P2b -- dielektrische achsensymmetrische Elektrostatik\n\n");
  lib.print(stdout);
  std::printf("\n");
  S.materials.print(stdout);
  std::printf("\n");

  int exit_code = 0;
  auto fail = [&exit_code](const char* what) {
    std::fprintf(stderr, "PRUEFUNG FEHLGESCHLAGEN: %s\n", what);
    exit_code = 2;
  };

  // --- reference solution ---------------------------------------------------
  const DielectricSolution R = solve_dielectric(S);
  R.print(stdout);
  if (!R.mesh.validate().all_passed()) fail("Volumennetz");
  if (!R.audit.ok()) fail("Randbedingungs-Audit");

  R.mesh.write_csv(outdir);
  R.write_csv(outdir);
  lib.write_csv(outdir + "/materials_library.csv");
  S.materials.write_csv(outdir + "/materials.csv");
  // The outlines the figures are drawn from come from the VOLUME mesh
  // (device_outline.csv, written above).  DeviceGeometry is deliberately not
  // used for that: it is the P1 description, in which the emitter runs down to
  // the domain floor, and drawing it over a P2b field would show a body that
  // was never solved.

  // --- 1. mesh convergence --------------------------------------------------
  std::printf("\n=== Netzkonvergenz ===\n");
  std::vector<DielectricSolution> levels;
  for (int L = 0; L <= max_level; ++L) {
    DielectricSetup s = S;
    s.geometry.mesh_level = L;
    levels.push_back(solve_dielectric(s));
    std::printf("  Stufe %d: %lld Knoten, Q_E = %.8e C\n", L,
                static_cast<long long>(levels.back().fem.n_nodes), levels.back().Q_emitter);
  }
  {
    std::FILE* f = std::fopen((outdir + "/convergence_mesh.csv").c_str(), "w");
    std::fprintf(f, "# Netzkonvergenz: das Groessenfeld wird als Ganzes mit 2^(-Stufe/2) "
                    "skaliert\n");
    std::fprintf(f, "level,size_scale,nodes,nr,nz,half_bandwidth,Q_emitter_C,Q_extractor_C,"
                    "residual_C,interface_Dn_rel");
    for (const Probe& p : levels[0].probes)
      std::fprintf(f, ",phi_%s_V,E_%s_V_per_m", p.name.c_str(), p.name.c_str());
    std::fprintf(f, "\n");
    for (std::size_t k = 0; k < levels.size(); ++k) {
      const DielectricSolution& s = levels[k];
      std::fprintf(f, "%zu,%.9g,%lld,%lld,%lld,%lld,%.9e,%.9e,%.9e,%.9e", k,
                   mesh_level_scale(static_cast<int>(k)),
                   static_cast<long long>(s.fem.n_nodes), static_cast<long long>(s.mesh.grid.nr),
                   static_cast<long long>(s.mesh.grid.nz),
                   static_cast<long long>(s.fem.half_bandwidth), s.Q_emitter, s.Q_extractor,
                   s.fem.residual_inf, s.relative_interface_error());
      for (std::size_t p = 0; p < s.probes.size(); ++p)
        std::fprintf(f, ",%.9e,%.9e", s.phi_probe[p], s.Emag_probe[p]);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }
  Real worst_mesh_phi = 0.0, worst_mesh_E = 0.0;
  for (std::size_t p = 0; p < levels[0].probes.size(); ++p) {
    const std::size_t a = levels.size() - 2, b = levels.size() - 1;
    worst_mesh_phi =
        std::max(worst_mesh_phi,
                 std::abs(levels[b].phi_probe[p] - levels[a].phi_probe[p]) /
                     std::abs(S.applied_span()));
    worst_mesh_E = std::max(worst_mesh_E, rel_change(levels[a].Emag_probe[p],
                                                     levels[b].Emag_probe[p]));
  }
  const Real mesh_dQ = rel_change(levels[levels.size() - 2].Q_emitter, levels.back().Q_emitter);

  // --- 2. feed boundary position -------------------------------------------
  std::printf("\n=== Lage der Zulaufgrenze ===\n");
  const Real zf0 = S.geometry.liquid_feed_z;
  const std::vector<Real> feed_z{0.5 * zf0, zf0, 2.0 * zf0, 4.0 * zf0};
  std::vector<DielectricSolution> feeds;
  for (Real z : feed_z) {
    DielectricSetup s = S;
    s.geometry.mesh_level = std::min(ref_level, 2);
    s.geometry.liquid_feed_z = z;
    feeds.push_back(solve_dielectric(s));
    std::printf("  z_feed = %10.4g m: Q_E = %.8e C\n", z, feeds.back().Q_emitter);
  }
  {
    std::FILE* f = std::fopen((outdir + "/convergence_feed.csv").c_str(), "w");
    std::fprintf(f, "# Konvergenz gegen die Lage der Zulaufgrenze; das Netz ab dem Kegelfuss "
                    "ist dabei unveraendert\n");
    std::fprintf(f, "liquid_feed_z_m,nodes,Q_emitter_C");
    for (const Probe& p : feeds[0].probes)
      std::fprintf(f, ",phi_%s_V,E_%s_V_per_m", p.name.c_str(), p.name.c_str());
    std::fprintf(f, "\n");
    for (std::size_t k = 0; k < feeds.size(); ++k) {
      std::fprintf(f, "%.9e,%lld,%.9e", feed_z[k],
                   static_cast<long long>(feeds[k].fem.n_nodes), feeds[k].Q_emitter);
      for (std::size_t p = 0; p < feeds[k].probes.size(); ++p)
        std::fprintf(f, ",%.9e,%.9e", feeds[k].phi_probe[p], feeds[k].Emag_probe[p]);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }
  Real worst_feed_phi = 0.0, worst_feed_E = 0.0;
  for (std::size_t p = 0; p < feeds[0].probes.size(); ++p) {
    worst_feed_phi = std::max(worst_feed_phi,
                              std::abs(feeds[3].phi_probe[p] - feeds[2].phi_probe[p]) /
                                  std::abs(S.applied_span()));
    worst_feed_E = std::max(worst_feed_E, rel_change(feeds[2].Emag_probe[p],
                                                     feeds[3].Emag_probe[p]));
  }
  const bool feed_converged = worst_feed_phi < feed_truncation::kTolPhiOverSpan &&
                              worst_feed_E < feed_truncation::kTolFieldRelative;
  std::printf("  Trunkierungsgrenzen (vorab festgelegt): phi %.1e, |E| %.1e der Spannweite\n",
              feed_truncation::kTolPhiOverSpan, feed_truncation::kTolFieldRelative);
  std::printf("  gemessen bei der letzten Verdopplung   : phi %.3e, |E| %.3e  -> %s\n",
              worst_feed_phi, worst_feed_E, feed_converged ? "konvergiert" : "NICHT KONVERGIERT");
  if (!feed_converged)
    std::printf("  BEFUND: die Lage der Zulaufgrenze ist eine Abmessung des Modells, kein\n"
                "  Konvergenzparameter.  Die Fluessigkeitssaeule ist ein Leiter auf "
                "V_emitter;\n  eine laengere Saeule traegt mehr Ladung, und die wird an der "
                "Spitze gespuert.\n  Jede berichtete Zahl gilt fuer liquid_feed_z = %.6g m.\n",
                S.geometry.liquid_feed_z);

  // --- 3. permittivity sensitivity -----------------------------------------
  std::printf("\n=== Empfindlichkeit gegenueber eps_r ===\n");
  const Material su8 = S.materials.emitter_dielectric;
  std::vector<Real> eps_values{1.0};
  if (su8.has_range()) {
    for (int k = 0; k <= 6; ++k)
      eps_values.push_back(su8.eps_r_low +
                           (su8.eps_r_high - su8.eps_r_low) * static_cast<Real>(k) / 6);
  }
  eps_values.push_back(su8.relative_permittivity);
  std::sort(eps_values.begin(), eps_values.end());
  eps_values.erase(std::unique(eps_values.begin(), eps_values.end()), eps_values.end());
  std::vector<DielectricSolution> eps_runs;
  for (Real e : eps_values) {
    DielectricSetup s = S;
    s.geometry.mesh_level = std::min(ref_level, 2);
    s.materials.emitter_dielectric.relative_permittivity = e;
    s.materials.emitter_dielectric.status = MaterialStatus::Provisional;
    s.materials.extractor_carrier.relative_permittivity = e;
    s.materials.extractor_carrier.status = MaterialStatus::Provisional;
    eps_runs.push_back(solve_dielectric(s));
    std::printf("  eps_r = %5.3f: Q_E = %.8e C\n", e, eps_runs.back().Q_emitter);
  }
  {
    std::FILE* f = std::fopen((outdir + "/sensitivity_permittivity.csv").c_str(), "w");
    std::fprintf(f, "# Sensitivitaet gegenueber der relativen Permittivitaet von SU-8.\n");
    std::fprintf(f, "# Der Nominalwert ist VORLAEUFIG; aus dieser Studie darf keine "
                    "Validierungsaussage folgen.\n");
    std::fprintf(f, "eps_r,nodes,Q_emitter_C,Q_extractor_C");
    for (const Probe& p : eps_runs[0].probes)
      std::fprintf(f, ",phi_%s_V,E_%s_V_per_m", p.name.c_str(), p.name.c_str());
    std::fprintf(f, "\n");
    for (std::size_t k = 0; k < eps_runs.size(); ++k) {
      std::fprintf(f, "%.9g,%lld,%.9e,%.9e", eps_values[k],
                   static_cast<long long>(eps_runs[k].fem.n_nodes), eps_runs[k].Q_emitter,
                   eps_runs[k].Q_extractor);
      for (std::size_t p = 0; p < eps_runs[k].probes.size(); ++p)
        std::fprintf(f, ",%.9e,%.9e", eps_runs[k].phi_probe[p], eps_runs[k].Emag_probe[p]);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }
  const Index k_tip = probe_index(eps_runs[0].probes, "axis_2_bore_radii");
  Real eps_lo_E = 0.0, eps_hi_E = 0.0, eps_nom_E = 0.0;
  for (std::size_t k = 0; k < eps_values.size(); ++k) {
    const Real E = eps_runs[k].Emag_probe[static_cast<std::size_t>(k_tip)];
    if (eps_values[k] == su8.eps_r_low) eps_lo_E = E;
    if (eps_values[k] == su8.eps_r_high) eps_hi_E = E;
    if (eps_values[k] == su8.relative_permittivity) eps_nom_E = E;
  }

  // --- 4. far-field treatment ----------------------------------------------
  std::printf("\n=== Fernrandbehandlung ===\n");
  const Real box_scale[4] = {0.25, 1.0, 2.0, 4.0};
  std::vector<std::array<Real, 5>> box_rows;  // R, phi_asym, phi_gnd, Q_asym, Q_gnd
  for (Real bs : box_scale) {
    std::array<Real, 5> row{};
    row[0] = S.geometry.device.domain_radius * bs;
    for (int mode = 0; mode < 2; ++mode) {
      DielectricSetup s = S;
      s.geometry.mesh_level = std::min(ref_level, 2);
      s.geometry.device.domain_radius = S.geometry.device.domain_radius * bs;
      s.geometry.device.domain_z_min = S.geometry.device.domain_z_min * bs;
      s.geometry.device.domain_z_max = S.geometry.device.domain_z_max * bs;
      s.far_field = mode ? FarField::Grounded : FarField::Asymptotic;
      const DielectricSolution sol = solve_dielectric(s);
      const Index i = probe_index(sol.probes, "axis_gap_mid");
      row[1 + mode] = sol.phi_probe[static_cast<std::size_t>(i)];
      row[3 + mode] = sol.Q_emitter;
    }
    box_rows.push_back(row);
    std::printf("  R = %8.3f mm: asymptotisch %11.5f V, geerdet %11.5f V, Differenz %.3e V\n",
                row[0] * 1e3, row[1], row[2], std::abs(row[1] - row[2]));
  }
  {
    std::FILE* f = std::fopen((outdir + "/farfield_study.csv").c_str(), "w");
    std::fprintf(f, "# Der offene Rand: asymptotische Monopolbedingung gegen geerdete Huelle.\n");
    std::fprintf(f, "# Die Differenz der beiden ist ein direktes Mass des "
                    "Trunkierungsfehlers.\n");
    std::fprintf(f, "domain_radius_m,phi_mid_asymptotic_V,phi_mid_grounded_V,"
                    "Q_emitter_asymptotic_C,Q_emitter_grounded_C\n");
    for (const auto& r : box_rows)
      std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,%.9e\n", r[0], r[1], r[2], r[3], r[4]);
    std::fclose(f);
  }
  const Real farfield_gap = std::abs(box_rows[1][1] - box_rows[1][2]) / std::abs(S.applied_span());

  // --- 5. FEM against the independent BEM at eps_r = 1 ----------------------
  std::printf("\n=== FEM gegen BEM (eps_r = 1, metallische Referenzanordnung) ===\n");
  Real bem_worst_phi = 0.0, bem_worst_E = 0.0;
  {
    DielectricSetup s = S;
    s.geometry.mesh_level = std::min(ref_level, 3);
    s.conductor_model = ConductorModel::MetallicReference;
    s.metallisation = Metallisation::AllSurfaces;
    s.materials = DielectricMaterials::all_vacuum(lib);
    const DielectricSolution fem = solve_dielectric(s);

    DeviceParameters dp = S.geometry.device;
    dp.emitter_back_length = -S.geometry.liquid_feed_z;
    const DeviceGeometry g = DeviceGeometry::build(dp);
    const BoundaryMesh bm = BoundaryMesh::generate(g, 0.5);
    BemSolver bem(vacuum_bem_mesh(bm, g));
    bem.solve({{S.V_emitter, S.V_extractor, 0.0}});

    std::FILE* f = std::fopen((outdir + "/fem_vs_bem.csv").c_str(), "w");
    std::fprintf(f, "# Querpruefung der FEM gegen die unabhaengige BEM fuer eps_r = 1.\n");
    std::fprintf(f, "# Beide loesen die SUPERSEDED metallische Anordnung -- das ist die "
                    "einzige, die eine\n# Einfachschicht-BEM darstellen kann.  Das prueft den "
                    "Loeser, NICHT die Physik von P2b.\n");
    std::fprintf(f, "probe,r_m,z_m,phi_fem_V,phi_bem_V,dphi_over_span,E_fem_V_per_m,"
                    "E_bem_V_per_m,E_rel\n");
    for (std::size_t k = 0; k < fem.probes.size(); ++k) {
      const Vec2 x = fem.probes[k].x;
      const Real pb = bem.potential_at(x);
      const Real eb = norm(bem.field_at(x));
      const Real dphi = std::abs(fem.phi_probe[k] - pb) / std::abs(S.applied_span());
      const Real de = rel_change(fem.Emag_probe[k], eb);
      bem_worst_phi = std::max(bem_worst_phi, dphi);
      bem_worst_E = std::max(bem_worst_E, de);
      std::fprintf(f, "%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", fem.probes[k].name.c_str(),
                   x.r, x.z, fem.phi_probe[k], pb, dphi, fem.Emag_probe[k], eb, de);
      std::printf("  %-26s phi FEM %11.5g BEM %11.5g (%.2e)  |E| FEM %10.4g BEM %10.4g "
                  "(%.2e)\n",
                  fem.probes[k].name.c_str(), fem.phi_probe[k], pb, dphi, fem.Emag_probe[k], eb,
                  de);
    }
    std::fclose(f);
  }

  // --- 6. linearity and polarity -------------------------------------------
  Real lin_err = 0.0, pol_err = 0.0;
  {
    DielectricSetup a = S;
    a.geometry.mesh_level = std::min(ref_level, 2);
    const DielectricSolution s1 = solve_dielectric(a);
    DielectricSetup b = a;
    b.V_emitter *= 2.5;
    b.V_extractor *= 2.5;
    const DielectricSolution s2 = solve_dielectric(b);
    DielectricSetup c = a;
    c.V_emitter = -a.V_emitter;
    c.V_extractor = -a.V_extractor;
    const DielectricSolution s3 = solve_dielectric(c);
    Real scale = 0.0;
    for (std::size_t k = 0; k < s1.fem.phi.size(); ++k) {
      lin_err = std::max(lin_err, std::abs(s2.fem.phi[k] - 2.5 * s1.fem.phi[k]));
      pol_err = std::max(pol_err, std::abs(s3.fem.phi[k] + s1.fem.phi[k]));
      scale = std::max(scale, std::abs(s1.fem.phi[k]));
    }
    lin_err /= scale;
    pol_err /= scale;
    std::FILE* f = std::fopen((outdir + "/linearity_polarity.csv").c_str(), "w");
    std::fprintf(f, "# Linearitaet und Polaritaetsumkehr des linearen Randwertproblems\n");
    std::fprintf(f, "quantity,value\n");
    std::fprintf(f, "max_rel_error_phi_scaled_2.5,%.9e\n", lin_err);
    std::fprintf(f, "max_rel_error_phi_reversed,%.9e\n", pol_err);
    std::fprintf(f, "Q_emitter_V,%.9e\n", s1.Q_emitter);
    std::fprintf(f, "Q_emitter_2.5V,%.9e\n", s2.Q_emitter);
    std::fprintf(f, "Q_emitter_minus_V,%.9e\n", s3.Q_emitter);
    std::fclose(f);
    std::printf("\nLinearitaet %.2e, Polaritaetsumkehr %.2e (relativ)\n", lin_err, pol_err);
  }

  // --- 7. field maps for the figures ---------------------------------------
  {
    const DeviceVolumeMesh& m = R.mesh;
    const Real L = S.geometry.device.extraction_distance;
    // Two windows chosen so that both come out close to square at equal aspect
    // -- the meridian half-plane is drawn to scale, and a window with a ten to
    // one aspect ratio would be honest but unreadable.
    const std::vector<Window> ws{
        {"uebersicht", 0.0, 1.6 * L, 1.15 * S.geometry.liquid_feed_z,
         1.25 * (L + S.geometry.device.extractor_thickness), 360, 360},
        {"spitze", 0.0, 4.0 * m.r_foot, -4.0 * m.r_foot, 4.0 * m.r_foot, 360, 360}};
    for (const Window& w : ws)
      write_grid_csv(outdir + "/field_" + w.name + ".csv", R, w);
  }

  // --- report ---------------------------------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/report.txt").c_str(), "w");
    std::fprintf(f, "P2b -- dielektrische achsensymmetrische Elektrostatik des kapillaren "
                    "Kunze-Emitters\n");
    std::fprintf(f, "============================================================================"
                    "====\n\n");
    std::fprintf(f, "KORRIGIERTER MODELLVERTRAG\n");
    std::fprintf(f, "  Auf Hochspannung liegt die ionische Fluessigkeit, nicht der Emitter.\n");
    std::fprintf(f, "  Der Emitterkoerper ist ein Dielektrikum (%s, eps_r = %.4g, Status %s).\n",
                 S.materials.emitter_dielectric.name.c_str(),
                 S.materials.emitter_dielectric.relative_permittivity,
                 to_string(S.materials.emitter_dielectric.status));
    std::fprintf(f, "  Der Extraktor ist ein Polymertraeger mit metallisierter Flaeche "
                    "(%s) auf V_extractor.\n", to_string(S.metallisation));
    std::fprintf(f, "  Das Reservoir ist nicht vernetzt: die Fluessigkeit endet bei "
                    "z = %.6g m,\n  und NUR ihr Querschnitt dort traegt V_emitter.  Die "
                    "uebrige Schnittebene ist\n  die Rueckflaeche des Polymers und keine "
                    "Elektrode.\n", S.geometry.liquid_feed_z);
    std::fprintf(f, "  P2a behandelte den Emitterkoerper als Metall; alle davon abhaengigen "
                    "P2a-Zahlen\n  sind damit ueberholt, nicht nur ungenauer.\n\n");

    std::fprintf(f, "VERFAHREN\n");
    std::fprintf(f, "  Achsensymmetrische Q1-Finite-Elemente auf einem blockstrukturierten, "
                    "radial\n  gewarpten Volumennetz, das aus den Geraeteparametern folgt. "
                    "2*pi*r-Gewichtung,\n  stueckweise Permittivitaet, symmetrischer "
                    "Banddirektloeser.\n");
    std::fprintf(f, "  Fernrand: asymptotische Monopolbedingung als Robin-Term.\n");
    std::fprintf(f, "  Referenzstufe %d: %lld Knoten, Halbbandbreite %lld, Residuum %.2e C.\n\n",
                 ref_level, static_cast<long long>(R.fem.n_nodes),
                 static_cast<long long>(R.fem.half_bandwidth), R.fem.residual_inf);

    std::fprintf(f, "MATERIALWERT SU-8\n");
    std::fprintf(f, "  eps_r = %.4g  [%s]\n", S.materials.emitter_dielectric.relative_permittivity,
                 to_string(S.materials.emitter_dielectric.status));
    std::fprintf(f, "  Quelle    : %s\n", S.materials.emitter_dielectric.source.c_str());
    std::fprintf(f, "  Bedingung : %s\n", S.materials.emitter_dielectric.conditions.c_str());
    std::fprintf(f, "  Vorbehalt : %s\n\n", S.materials.emitter_dielectric.caveat.c_str());

    std::fprintf(f, "ERGEBNISSE AUF DER REFERENZSTUFE\n");
    std::fprintf(f, "  V_emitter = %.6g V, V_extractor = %.6g V\n", S.V_emitter, S.V_extractor);
    std::fprintf(f, "  Q_emitter = %.9e C, Q_extractor = %.9e C, Summe = %.9e C\n", R.Q_emitter,
                 R.Q_extractor, R.Q_net);
    std::fprintf(f, "  Die Summe ist nicht null: bei phi -> 0 im Unendlichen traegt das "
                    "System Nettoladung.\n");
    std::fprintf(f, "  %-28s %11s %11s %11s %13s %13s\n", "Punkt", "r [m]", "z [m]",
                 "Kantenabstand", "phi [V]", "|E| [V/m]");
    for (std::size_t k = 0; k < R.probes.size(); ++k)
      std::fprintf(f, "  %-28s %11.4g %11.4g %11.4g %13.6g %13.6g\n", R.probes[k].name.c_str(),
                   R.probes[k].x.r, R.probes[k].x.z, R.probes[k].clearance, R.phi_probe[k],
                   R.Emag_probe[k]);
    std::fprintf(f, "\n  Auf der ebenen Fluessigkeitsreferenz bei z = 0+ liegt E_z bei r = 0 "
                    "bei %.6g V/m;\n  %lld kantennahe Zellen sind ausgeschlossen.  Die "
                    "Austrittskante ist unverrundet:\n  ihr Feld divergiert und folgt der "
                    "Elementgroesse.  Ein \"konvergiertes Spitzenfeld\"\n  gibt es dort nicht "
                    "und wird nicht berichtet.\n",
                 R.surface_Ez.empty() ? 0.0 : R.surface_Ez.front(),
                 static_cast<long long>(R.surface_edge_cells));

    std::fprintf(f, "\nKONVERGENZ UND EMPFINDLICHKEIT\n");
    std::fprintf(f, "  Netz, letzte Verfeinerung : phi %.2e der Spannweite, |E| %.2e relativ, "
                    "Q_E %.2e relativ\n", worst_mesh_phi, worst_mesh_E, mesh_dQ);
    std::fprintf(f, "  Zulaufgrenze, letzte Verdopplung: phi %.2e der Spannweite, |E| %.2e "
                    "relativ\n", worst_feed_phi, worst_feed_E);
    std::fprintf(f, "    vorab festgelegte Grenzen: phi %.1e, |E| %.1e  -> %s\n",
                 feed_truncation::kTolPhiOverSpan, feed_truncation::kTolFieldRelative,
                 feed_converged ? "eingehalten" : "NICHT EINGEHALTEN");
    std::fprintf(f, "  Fernrand bei R = %.4g m: asymptotisch gegen geerdet %.2e der "
                    "Spannweite\n", S.geometry.device.domain_radius, farfield_gap);
    std::fprintf(f, "  D_n-Sprung an der Dielektrikumsgrenze: %.2e relativ (faellt mit der "
                    "Verfeinerung)\n", R.relative_interface_error());
    std::fprintf(f, "  FEM gegen BEM bei eps_r = 1: phi %.2e der Spannweite, |E| %.2e "
                    "relativ\n", bem_worst_phi, bem_worst_E);
    std::fprintf(f, "  Linearitaet %.2e, Polaritaetsumkehr %.2e\n", lin_err, pol_err);
    if (su8.has_range())
      std::fprintf(f, "  eps_r von %.3g bis %.3g aendert |E| bei zwei Bohrungsradien um "
                      "%.2f %% des Nominalwerts\n    (%.6g bis %.6g V/m, nominal %.6g V/m). "
                      "KEINE Validierungsaussage.\n",
                   su8.eps_r_low, su8.eps_r_high,
                   eps_nom_E > 0.0 ? 100.0 * (eps_hi_E - eps_lo_E) / eps_nom_E : 0.0, eps_lo_E,
                   eps_hi_E, eps_nom_E);

    if (!feed_converged) {
      std::fprintf(f, "\nBEFUND: DIE LAGE DER ZULAUFGRENZE IST NICHT AUSKONVERGIERT\n");
      std::fprintf(f, "  Verlangt war der Nachweis, dass eine weitere Rueckverlagerung der\n"
                      "  Zulaufgrenze Potential und Feld am Meniskus nicht mehr aendert.  Sie\n"
                      "  aendert sie.  Gemessen (convergence_feed.csv): eine Verdopplung der\n"
                      "  modellierten Saeulenlaenge von %.4g m auf %.4g m verschiebt das\n"
                      "  Potential um %.2e der Spannweite und das Feld um %.2e relativ.  Die\n"
                      "  vorab festgelegten Grenzen von %.1e werden um mehr als eine\n"
                      "  Groessenordnung verfehlt.  Die Grenze wird nicht verschoben.\n",
                   -feed_z[2], -feed_z[3], worst_feed_phi, worst_feed_E,
                   feed_truncation::kTolPhiOverSpan);
      std::fprintf(f, "  URSACHE, und warum das kein numerischer Mangel ist: die\n"
                      "  Fluessigkeitssaeule ist ein Leiter auf V_emitter.  Eine duenne Saeule\n"
                      "  der Laenge L und des Radius a hat eine Selbstkapazitaet von etwa\n"
                      "  2 pi eps0 L / (ln(2L/a) - 1); eine laengere Saeule traegt "
                      "proportional\n"
                      "  mehr Ladung, und die wird an der Spitze gespuert.  Die gemessene\n"
                      "  Emitterladung folgt dieser Formel auf ein festes Verhaeltnis genau.\n"
                      "  Am offenen Rand liegt es nicht: dieselbe Studie in einer GEERDETEN\n"
                      "  Huelle bei 25 mm liefert dieselben Zahlen auf besser als 0.1 Prozent.\n"
                      "  Es ist derselbe Mechanismus, den P2a fuer emitter_back_length fand.\n");
      std::fprintf(f, "  WAS ES BEHEBEN WUERDE: der Emitterhalter.  04_geometry_model.md,\n"
                      "  Abschnitt 4.1, sieht eine Basisplatte auf Emitterpotential vor, aus\n"
                      "  der die verjuengte Struktur herausragt.  Ein echter Leiter mit\n"
                      "  angegebenen Abmessungen wuerde die Zulaufposition wieder zu einem\n"
                      "  numerischen Parameter machen.  Das ist eine GEOMETRIEENTSCHEIDUNG "
                      "fuer\n"
                      "  eine spaetere Phase.  Die hintere Schnittebene als solche zur\n"
                      "  Elektrode zu erklaeren waere keine -- und ist hier ausgeschlossen.\n");
      std::fprintf(f, "  KONSEQUENZ: liquid_feed_z = %.6g m ist ein BEISPIELWERT, keine\n"
                      "  gemessene Abmessung.  Jede Zahl in diesem Bericht gilt fuer diesen\n"
                      "  Wert.  Nichts in P2b ist bezueglich der Zulauftrunkierung "
                      "konvergiert.\n",
                   S.geometry.liquid_feed_z);
    }

    std::fprintf(f, "\nWAS HIER NICHT MODELLIERT IST\n");
    std::fprintf(f, "  %s\n", liquid_model::why_ideal_conductor_is_admissible_in_p2b());
    std::fprintf(f, "  Die Zulaufgrenze ist bereits als HYDRAULISCHER Zulauf vorgesehen: dort "
                    "wird in P3\n  Reservoirdruck, Volumenstrom oder eine hydraulische Impedanz "
                    "vorgegeben.  In\n  diesem Lauf ist nichts davon implementiert.\n");
    std::fprintf(f, "  Keine Raumladung, keine Emission, keine Meniskusbewegung, keine "
                    "Stroemung,\n  keine Zeitabhaengigkeit, keine 3D-Rechnung.\n");
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_dielectric (P2b)\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "V_emitter_V=%.9e\n", S.V_emitter);
    std::fprintf(f, "V_extractor_V=%.9e\n", S.V_extractor);
    std::fprintf(f, "reference_level=%d\n", ref_level);
    std::fprintf(f, "max_level=%d\n", max_level);
    std::fprintf(f, "reference_size_scale=%.9g\n", mesh_level_scale(ref_level));
    std::fprintf(f, "nodes=%lld\n", static_cast<long long>(R.fem.n_nodes));
    std::fprintf(f, "nr=%lld\n", static_cast<long long>(R.mesh.grid.nr));
    std::fprintf(f, "nz=%lld\n", static_cast<long long>(R.mesh.grid.nz));
    std::fprintf(f, "liquid_feed_z_m=%.9e\n", S.geometry.liquid_feed_z);
    std::fprintf(f, "conductor_model=%s\n", to_string(S.conductor_model));
    std::fprintf(f, "metallisation=%s\n", to_string(S.metallisation));
    std::fprintf(f, "far_field=%s\n", to_string(S.far_field));
    std::fprintf(f, "emitter_material=%s\n", S.materials.emitter_dielectric.name.c_str());
    std::fprintf(f, "emitter_eps_r=%.9g\n", S.materials.emitter_dielectric.relative_permittivity);
    std::fprintf(f, "emitter_eps_r_status=%s\n",
                 to_string(S.materials.emitter_dielectric.status));
    std::fprintf(f, "emitter_eps_r_low=%.9g\n", S.materials.emitter_dielectric.eps_r_low);
    std::fprintf(f, "emitter_eps_r_high=%.9g\n", S.materials.emitter_dielectric.eps_r_high);
    std::fprintf(f, "Q_emitter_C=%.9e\n", R.Q_emitter);
    std::fprintf(f, "Q_extractor_C=%.9e\n", R.Q_extractor);
    std::fprintf(f, "Q_net_C=%.9e\n", R.Q_net);
    std::fprintf(f, "mesh_change_phi_over_span=%.9e\n", worst_mesh_phi);
    std::fprintf(f, "mesh_change_E_rel=%.9e\n", worst_mesh_E);
    std::fprintf(f, "mesh_change_Q_rel=%.9e\n", mesh_dQ);
    std::fprintf(f, "feed_change_phi_over_span=%.9e\n", worst_feed_phi);
    std::fprintf(f, "feed_change_E_rel=%.9e\n", worst_feed_E);
    std::fprintf(f, "feed_tol_phi_over_span=%.9e\n", feed_truncation::kTolPhiOverSpan);
    std::fprintf(f, "feed_tol_E_rel=%.9e\n", feed_truncation::kTolFieldRelative);
    std::fprintf(f, "feed_converged=%s\n", feed_converged ? "yes" : "no");
    std::fprintf(f, "farfield_gap_over_span=%.9e\n", farfield_gap);
    std::fprintf(f, "interface_Dn_rel=%.9e\n", R.relative_interface_error());
    std::fprintf(f, "fem_vs_bem_phi_over_span=%.9e\n", bem_worst_phi);
    std::fprintf(f, "fem_vs_bem_E_rel=%.9e\n", bem_worst_E);
    std::fprintf(f, "linearity_rel=%.9e\n", lin_err);
    std::fprintf(f, "polarity_rel=%.9e\n", pol_err);
    std::fprintf(f, "polymer_dirichlet_nodes=%lld\n",
                 static_cast<long long>(R.audit.n_polymer_dirichlet));
    std::fprintf(f, "feed_plane_outside_liquid_nodes=%lld\n",
                 static_cast<long long>(R.audit.n_feed_plane_outside_liquid));
    std::fprintf(f, "state=dielectric electrostatics, rho_f=0, liquid an ideal conductor, "
                    "no emission, no flow\n");
    std::fclose(f);
  }

  cfg.warn_about_unused(stdout, {"meta.", "fluid.", "beam.", "output.", "bem.", "wetting.",
                                 "species.", "feed.mode", "feed.delta_p", "feed.impedance",
                                 "feed.Q", "feed.model"});
  std::printf("\ngeschrieben nach %s\n", outdir.c_str());
  (void)eps0;
  return exit_code;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\n%s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "\nFehler: %s\n", e.what());
  return 2;
}
