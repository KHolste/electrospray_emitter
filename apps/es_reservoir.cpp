// es_reservoir -- P2c: the liquid reservoir as a modelled body instead of a
// cut through the liquid column.
//
//   es_reservoir <geometrie.cfg> [<p2c.cfg> ...] <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS FOR
//
// Until now one number, the position of the "liquid feed boundary", moved three
// different things at the same time: the length of the CONDUCTING liquid
// column, the rearward extent of the DIELECTRIC body, and with them the whole
// rear geometry of the device.  A study that varied it was therefore never
// moving a boundary condition -- it was building a different high-voltage
// electrode each time -- and reporting the resulting field change as "does not
// converge against the position of the feed boundary" was misleading.
//
// This run separates the two questions:
//
//   * the FRONT of the device -- taper, bore, base body, extractor and the
//     near-field mesh -- is held fixed, and that is proved bitwise, not
//     asserted;
//   * only the modelled liquid reservoir behind it is varied, as an
//     axisymmetric plenum inside a DIELECTRIC body.
//
// No conducting holder, no rearward metal disc and no base plate on emitter
// potential is introduced.  The whole connected liquid -- bore, feed channel,
// plenum -- is one equipotential at V_emitter; the surrounding polymer stays a
// dielectric, and the boundary audit is what proves it.
//
// The old truncated-column arrangement is computed as well, and labelled as the
// DIAGNOSIS it is.
//
// NOT in this phase: surface tension, meniscus computation, flow, emission,
// space charge, 3D.  The plane at z = 0 is the initial flat liquid surface.
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

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/dielectric_device.hpp"
#include "es/io.hpp"

using namespace es;

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
  if (c.has("device.emitter_back_length") && c.num("device.emitter_back_length", 0.0) != 0.0)
    throw std::runtime_error(
        "device.emitter_back_length ist gesetzt.  Das ist die leitende P2a-Abschlussscheibe. "
        "Sie bleibt abgelehnt: ein Fluessigkeitsvorrat wird hier als dielektrisch "
        "umschlossener Fluessigkeitsraum modelliert, nicht als rueckwaertige Metallscheibe.");
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
      "sein; angegeben war: '" + s + "'");
}

void apply_material_overrides(const Config& c, MaterialLibrary& lib) {
  for (const char* raw : {"su8", "ip-q", "ipx-q", "peek"}) {
    const std::string name(raw);
    const std::string key = "material." + name + ".relative_permittivity";
    if (!c.has(key)) continue;
    const std::string src = c.str("material." + name + ".source", "");
    if (src.empty())
      throw std::runtime_error(
          "material." + name + ".relative_permittivity wurde gesetzt, aber material." + name +
          ".source fehlt.  Ein Materialwert ohne Herkunft ist keine Eingabe, sondern eine "
          "Behauptung.");
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

Real rel_change(Real a, Real b) {
  const Real d = std::max(std::abs(a), std::abs(b));
  return d > 0.0 ? std::abs(a - b) / d : 0.0;
}

Index probe_index(const std::vector<Probe>& p, const char* name) {
  for (std::size_t k = 0; k < p.size(); ++k)
    if (p[k].name == name) return static_cast<Index>(k);
  return -1;
}

// ---------------------------------------------------------------------------

/// One entry of the comparison.  `label` is what the figures show.
struct Variant {
  std::string tag, label;
  ReservoirModel model{ReservoirModel::AxisymmetricPlenum};
  Real plenum_radius{0.0}, plenum_depth{0.0}, fill{1.0};
  bool is_diagnosis{false};
};

struct Window {
  Real r0, r1, z0, z1;
  Index nr, nz;
};

/// Sample phi, |E| and the material region on a uniform window.  The window is
/// the SAME for every variant, which is what makes the figures comparable.
void write_grid_csv(const std::string& path, const DielectricSolution& s, const Window& w) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# Nahfeldausschnitt, fuer alle Vorratsvarianten identisch\n");
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
      std::fprintf(f, "%.7e,%.7e,%.7e,%.7e,%s\n", r, z, phi, E, to_string(s.mesh.region_at(x)));
    }
  }
  std::fclose(f);
}

/// Meridian outlines of one variant, appended with a variant column, so that
/// the geometry figure draws every variant to the same scale from one file.
void append_outline(std::FILE* f, const std::string& tag, const DeviceVolumeMesh& m) {
  auto poly = [&](const char* name, const std::vector<Vec2>& pts) {
    for (std::size_t k = 0; k < pts.size(); ++k)
      std::fprintf(f, "%s,%s,%zu,%.9e,%.9e\n", tag.c_str(), name, k, pts[k].r, pts[k].z);
  };
  const Real H = m.p.device.emitter_height;
  const Real zb = m.z_base;
  const Real ze = m.p.device.extraction_distance;
  const Real zt = ze + m.p.device.extractor_thickness;
  poly("emitter_dielectric",
       {{m.r_bore, zb}, {m.r_foot, zb}, {m.r_foot, -H}, {m.r_land, 0.0}, {m.r_bore, 0.0}});
  poly("extractor_carrier",
       {{m.r_aperture, ze}, {m.r_ext_outer, ze}, {m.r_ext_outer, zt}, {m.r_aperture, zt}});
  if (m.has_plenum()) {
    poly("liquid", {{0.0, m.z_fill},
                    {m.r_plenum, m.z_fill},
                    {m.r_plenum, m.z_roof},
                    {m.r_channel, m.z_roof},
                    {m.r_channel, zb},
                    {m.r_bore, zb},
                    {m.r_bore, 0.0},
                    {0.0, 0.0}});
    poly("reservoir_dielectric", {{m.r_channel, m.z_block_bottom},
                                  {m.r_plenum_outer, m.z_block_bottom},
                                  {m.r_plenum_outer, zb},
                                  {m.r_channel, zb},
                                  {m.r_channel, m.z_roof},
                                  {m.r_plenum, m.z_roof},
                                  {m.r_plenum, m.z_cav_bottom},
                                  {m.r_channel, m.z_cav_bottom}});
  } else {
    poly("liquid", {{0.0, zb}, {m.r_bore, zb}, {m.r_bore, 0.0}, {0.0, 0.0}});
    poly("liquid_column_cut", {{0.0, zb}, {m.r_bore, zb}});
  }
}

}  // namespace

// ---------------------------------------------------------------------------

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.size() < 2) {
    std::fprintf(stderr, "es_reservoir <geometrie.cfg> [<p2c.cfg> ...] <ausgabeverzeichnis> "
                         "[key=value ...]\n");
    return 1;
  }
  Config cfg;
  for (std::size_t k = 0; k + 1 < pos.size(); ++k) cfg.load(pos[k]);
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  MaterialLibrary lib;
  apply_material_overrides(cfg, lib);

  DielectricSetup S;
  S.geometry.device = device_from(cfg);
  S.geometry.base_plate_thickness =
      cfg.num("device.base_plate_thickness", S.geometry.base_plate_thickness);
  S.geometry.feed_channel_length =
      cfg.num("reservoir.feed_channel_length", S.geometry.feed_channel_length);
  S.geometry.feed_channel_radius = cfg.num("reservoir.feed_channel_radius", 0.0);
  S.geometry.plenum_wall_thickness =
      cfg.num("reservoir.wall_thickness", S.geometry.plenum_wall_thickness);
  const int ref_level = cfg.integer("mesh.reference_level", 3);
  const int max_level = cfg.integer("mesh.max_level", 4);
  S.geometry.mesh_level = ref_level;
  S.materials = DielectricMaterials::reference(lib);
  S.materials.emitter_dielectric = lib.get(cfg.str("emitter.material", "su8"));
  S.materials.extractor_carrier = lib.get(cfg.str("extractor.material", "su8"));
  if (cfg.has("reservoir.material")) S.materials.reservoir_body = lib.get(cfg.str(
      "reservoir.material", "su8"));
  S.metallisation = metallisation_from(cfg.str("extractor.metallisation", "front_and_aperture"));
  S.far_field = FarField::Asymptotic;
  if (!cfg.has("field.V_emitter") || !cfg.has("field.V_extractor"))
    throw std::runtime_error("field.V_emitter und field.V_extractor sind Pflichtangaben.");
  S.V_emitter = cfg.num("field.V_emitter", 0.0);
  S.V_extractor = cfg.num("field.V_extractor", 0.0);

  // --- the comparison ------------------------------------------------------
  //
  // Four plenum sizes, each about a factor 2.9 in liquid volume above the last;
  // the fourth is the additional enlargement step that checks whether the third
  // was already converged.  The truncated column runs first, as a diagnosis of
  // the superseded arrangement.  A partly filled plenum is included so that the
  // fill level is exercised rather than merely offered.
  std::vector<Variant> variants{
      {"saeule", "abgeschnittene Saeule (Diagnose)", ReservoirModel::TruncatedColumn, 0.0, 0.0,
       1.0, true},
      {"plenum_a", "Plenum r = 2.5 mm, Tiefe 1.0 mm", ReservoirModel::AxisymmetricPlenum,
       2.5e-3, 1.0e-3, 1.0, false},
      {"plenum_b", "Plenum r = 3.5 mm, Tiefe 1.5 mm", ReservoirModel::AxisymmetricPlenum,
       3.5e-3, 1.5e-3, 1.0, false},
      {"plenum_c", "Plenum r = 5.0 mm, Tiefe 2.5 mm", ReservoirModel::AxisymmetricPlenum,
       5.0e-3, 2.5e-3, 1.0, false},
      {"plenum_d", "Plenum r = 7.0 mm, Tiefe 4.0 mm (Vergroesserungsstufe)",
       ReservoirModel::AxisymmetricPlenum, 7.0e-3, 4.0e-3, 1.0, false},
      {"plenum_c_halb", "Plenum r = 5.0 mm, Tiefe 2.5 mm, halb gefuellt",
       ReservoirModel::AxisymmetricPlenum, 5.0e-3, 2.5e-3, 0.5, false},
  };

  auto setup_of = [&](const Variant& v, int level) {
    DielectricSetup s = S;
    s.geometry.mesh_level = level;
    s.geometry.reservoir = v.model;
    s.geometry.plenum_radius = v.plenum_radius;
    s.geometry.plenum_depth = v.plenum_depth;
    s.geometry.plenum_fill_fraction = v.fill;
    return s;
  };

  std::printf("P2c -- der Fluessigkeitsvorrat als modellierter Koerper\n\n");
  lib.print(stdout);
  std::printf("\n");
  S.materials.print(stdout);
  std::printf("\n");

  int exit_code = 0;
  auto fail = [&exit_code](const char* what) {
    std::fprintf(stderr, "PRUEFUNG FEHLGESCHLAGEN: %s\n", what);
    exit_code = 2;
  };

  std::vector<DielectricSolution> sol;
  for (const Variant& v : variants) {
    sol.push_back(solve_dielectric(setup_of(v, ref_level)));
    const DielectricSolution& s = sol.back();
    if (!s.mesh.validate().all_passed()) fail("Volumennetz einer Variante");
    if (!s.audit.ok()) fail("Randbedingungs-Audit einer Variante");
    if (s.audit.n_polymer_dirichlet != 0) fail("festgehaltene Knoten ohne Fluessigkeitskontakt");
    std::printf("  %-14s %6lld x %-6lld Knoten, Q_liquid = %.6e C, Q_extractor = %.6e C\n",
                v.tag.c_str(), static_cast<long long>(s.mesh.grid.nr),
                static_cast<long long>(s.mesh.grid.nz), s.Q_emitter, s.Q_extractor);
  }

  // --- the identity proof ---------------------------------------------------
  //
  // Bitwise on every node coordinate and exact on every cell material, over the
  // whole device in front of the base body -- taper, bore, base body, gap and
  // extractor.  If this is not exactly zero, no comparison below means anything.
  Real worst_dr = 0.0, worst_dz = 0.0;
  Index cell_mismatch = 0, rows_compared = 0, cells_compared = 0;
  bool shape_ok = true;
  {
    const DeviceVolumeMesh& a = sol[1].mesh;   // the first plenum is the reference
    std::FILE* f = std::fopen((outdir + "/front_identity.csv").c_str(), "w");
    std::fprintf(f, "# Nachweis, dass die Vorratsgroesse die vordere Geraetegeometrie nicht "
                    "beruehrt.\n");
    std::fprintf(f, "# Verglichen wird gegen plenum_a, bitweise auf den Knotenkoordinaten und "
                    "exakt auf\n# dem Material jeder Zelle, fuer alle Zeilen ab der "
                    "Rueckflaeche des Grundkoerpers\n# und alle Radien bis zum "
                    "Extraktoraussenrand.\n");
    std::fprintf(f, "variant,rows_compared,cells_compared,max_abs_dr_m,max_abs_dz_m,"
                    "cell_material_mismatches,bitwise_identical\n");
    for (std::size_t k = 0; k < variants.size(); ++k) {
      const DeviceVolumeMesh& b = sol[k].mesh;
      Real dr = 0.0, dz = 0.0;
      Index mism = 0, rows = 0, cells = 0;
      bool ok = (a.i_ext_outer == b.i_ext_outer) &&
                (a.grid.nz - a.j_base == b.grid.nz - b.j_base);
      if (ok) {
        for (Index j = a.j_base, jb = b.j_base; j < a.grid.nz; ++j, ++jb) {
          dz = std::max(dz, std::abs(a.grid.z_of_row(j) - b.grid.z_of_row(jb)));
          for (Index i = 0; i <= a.i_ext_outer; ++i)
            dr = std::max(dr, std::abs(a.grid.at(i, j).r - b.grid.at(i, jb).r));
          ++rows;
          if (j + 1 >= a.grid.nz) continue;
          for (Index i = 0; i + 1 <= a.i_ext_outer; ++i) {
            ++cells;
            if (a.cell_region[static_cast<std::size_t>(a.grid.cell(i, j))] !=
                b.cell_region[static_cast<std::size_t>(b.grid.cell(i, jb))])
              ++mism;
          }
        }
      } else {
        shape_ok = false;
      }
      worst_dr = std::max(worst_dr, dr);
      worst_dz = std::max(worst_dz, dz);
      cell_mismatch += mism;
      rows_compared = rows;
      cells_compared = cells;
      std::fprintf(f, "%s,%lld,%lld,%.9e,%.9e,%lld,%d\n", variants[k].tag.c_str(),
                   static_cast<long long>(rows), static_cast<long long>(cells), dr, dz,
                   static_cast<long long>(mism),
                   (ok && dr == 0.0 && dz == 0.0 && mism == 0) ? 1 : 0);
    }
    std::fclose(f);
    if (!shape_ok || worst_dr != 0.0 || worst_dz != 0.0 || cell_mismatch != 0)
      fail("die vordere Geraetegeometrie ist NICHT identisch geblieben");
    std::printf("\nFrontgeometrie: %lld Zeilen, %lld Zellen je Variante verglichen; "
                "max|dr| = %.1e m, max|dz| = %.1e m, %lld Materialabweichungen\n",
                static_cast<long long>(rows_compared), static_cast<long long>(cells_compared),
                worst_dr, worst_dz, static_cast<long long>(cell_mismatch));
  }

  // --- the compared quantities ---------------------------------------------
  const Index k_tip = probe_index(sol[0].probes, "axis_2_bore_radii");
  const Index k_mid = probe_index(sol[0].probes, "axis_gap_mid");
  auto Ez_surface_axis = [](const DielectricSolution& s) {
    return s.surface_Ez.empty() ? 0.0 : s.surface_Ez.front();
  };
  auto Emag_at_two_bore_radii = [&](const DielectricSolution& s) {
    return s.Emag_probe[static_cast<std::size_t>(k_tip)];
  };
  // |E| at the second, independent evaluation radius on the flat liquid
  // surface: half way out to the pinned edge, still outside the excluded zone.
  auto Ez_surface_half = [](const DielectricSolution& s) {
    if (s.surface_Ez.empty()) return 0.0;
    return s.surface_Ez[s.surface_Ez.size() / 2];
  };

  {
    std::FILE* f = std::fopen((outdir + "/variants.csv").c_str(), "w");
    std::fprintf(f, "# Vergleich der Vorratsmodelle bei identischer Frontgeometrie.\n");
    std::fprintf(f, "# LOKALE Feldgroessen und GLOBALE Ladungsgroessen sind getrennt zu "
                    "beurteilen:\n# die Gesamtladung eines endlichen Leiters im offenen Raum "
                    "ist eine Funktion seiner\n# Groesse und konvergiert nicht, das lokale "
                    "Extraktionsfeld soll es.\n");
    std::fprintf(f, "variant,label,reservoir_model,plenum_radius_m,plenum_depth_m,"
                    "fill_fraction,liquid_volume_m3,reservoir_volume_m3,nodes,"
                    "Ez_axis_surface_V_per_m,Ez_surface_half_V_per_m,"
                    "Emag_2_bore_radii_V_per_m,phi_gap_mid_V,Q_liquid_C,Q_extractor_C,Q_net_C");
    for (const Probe& p : sol[0].probes)
      std::fprintf(f, ",phi_%s_V,E_%s_V_per_m", p.name.c_str(), p.name.c_str());
    std::fprintf(f, "\n");
    for (std::size_t k = 0; k < variants.size(); ++k) {
      const DielectricSolution& s = sol[k];
      std::fprintf(f, "%s,\"%s\",%s,%.9e,%.9e,%.9g,%.9e,%.9e,%lld,%.9e,%.9e,%.9e,%.9e,"
                      "%.9e,%.9e,%.9e",
                   variants[k].tag.c_str(), variants[k].label.c_str(),
                   to_string(variants[k].model), variants[k].plenum_radius,
                   variants[k].plenum_depth, variants[k].fill,
                   s.mesh.analytic_volume_of(Region::Liquid),
                   s.mesh.analytic_volume_of(Region::ReservoirSolid),
                   static_cast<long long>(s.fem.n_nodes), Ez_surface_axis(s),
                   Ez_surface_half(s), Emag_at_two_bore_radii(s),
                   s.phi_probe[static_cast<std::size_t>(k_mid)], s.Q_emitter, s.Q_extractor,
                   s.Q_net);
      for (std::size_t p = 0; p < s.probes.size(); ++p)
        std::fprintf(f, ",%.9e,%.9e", s.phi_probe[p], s.Emag_probe[p]);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }

  // Growth steps a -> b -> c -> d.  Every quantity is judged on the change per
  // step, against the tolerances fixed in advance.
  Real worst_phi_step = 0.0, worst_E_step = 0.0;
  Real last_phi_step = 0.0, last_E_step = 0.0;
  Real last_meniscus_step = 0.0;
  std::string worst_E_probe, worst_phi_probe;
  {
    // The verdict is the WORST probe, and which probe that is has to be visible
    // -- a run that quotes only the meniscus quantities and calls it converged
    // would be choosing its own criterion after the measurement.  So the file
    // carries the per-probe change as well as the maximum.
    std::FILE* f = std::fopen((outdir + "/reservoir_convergence.csv").c_str(), "w");
    std::fprintf(f, "# Aenderung je Vergroesserungsschritt des Plenums, bei bitgleicher "
                    "Frontgeometrie.\n");
    std::fprintf(f, "# Die Grenzen sind vor der Messung festgelegt (es::reservoir_convergence) "
                    "und werden\n# nicht nachtraeglich gelockert.  d_*_max_* ist das Maximum "
                    "ueber ALLE Sondenpunkte;\n# die Einzelspalten stehen daneben, damit "
                    "sichtbar ist, welcher Punkt das Maximum\n# treibt.\n");
    std::fprintf(f, "step,from,to,liquid_volume_ratio,d_phi_max_over_span,worst_phi_probe,"
                    "d_E_max_rel,worst_E_probe,d_Ez_axis_rel,d_Emag_2rb_rel,d_Q_liquid_rel,"
                    "d_Q_extractor_rel");
    for (const Probe& p : sol[0].probes)
      std::fprintf(f, ",d_phi_%s_over_span,d_E_%s_rel", p.name.c_str(), p.name.c_str());
    std::fprintf(f, "\n");
    for (std::size_t k = 2; k <= 4; ++k) {
      const DielectricSolution& a = sol[k - 1];
      const DielectricSolution& b = sol[k];
      Real dphi = 0.0, dE = 0.0;
      std::string pmax, emax;
      for (std::size_t p = 0; p < a.probes.size(); ++p) {
        const Real q = std::abs(b.phi_probe[p] - a.phi_probe[p]) / std::abs(S.applied_span());
        const Real e = rel_change(a.Emag_probe[p], b.Emag_probe[p]);
        if (q > dphi) {
          dphi = q;
          pmax = a.probes[p].name;
        }
        if (e > dE) {
          dE = e;
          emax = a.probes[p].name;
        }
      }
      worst_phi_step = std::max(worst_phi_step, dphi);
      worst_E_step = std::max(worst_E_step, dE);
      last_phi_step = dphi;
      last_E_step = dE;
      worst_phi_probe = pmax;
      worst_E_probe = emax;
      last_meniscus_step = std::max(rel_change(Ez_surface_axis(a), Ez_surface_axis(b)),
                                    rel_change(Emag_at_two_bore_radii(a),
                                               Emag_at_two_bore_radii(b)));
      std::fprintf(f, "%zu,%s,%s,%.6g,%.9e,%s,%.9e,%s,%.9e,%.9e,%.9e,%.9e", k - 1,
                   variants[k - 1].tag.c_str(), variants[k].tag.c_str(),
                   b.mesh.analytic_volume_of(Region::Liquid) /
                       a.mesh.analytic_volume_of(Region::Liquid),
                   dphi, pmax.c_str(), dE, emax.c_str(),
                   rel_change(Ez_surface_axis(a), Ez_surface_axis(b)),
                   rel_change(Emag_at_two_bore_radii(a), Emag_at_two_bore_radii(b)),
                   rel_change(a.Q_emitter, b.Q_emitter),
                   rel_change(a.Q_extractor, b.Q_extractor));
      for (std::size_t p = 0; p < a.probes.size(); ++p)
        std::fprintf(f, ",%.9e,%.9e",
                     std::abs(b.phi_probe[p] - a.phi_probe[p]) / std::abs(S.applied_span()),
                     rel_change(a.Emag_probe[p], b.Emag_probe[p]));
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }
  const bool local_converged = last_phi_step < reservoir_convergence::kTolPhiOverSpan &&
                               last_E_step < reservoir_convergence::kTolFieldRelative;
  const bool meniscus_converged =
      last_meniscus_step < reservoir_convergence::kTolFieldRelative;
  std::printf("\nLetzte Vergroesserungsstufe, Maximum ueber alle Sonden: phi %.3e (%s), "
              "|E| %.3e (%s)\n  Grenzen %.1e -> %s\n",
              last_phi_step, worst_phi_probe.c_str(), last_E_step, worst_E_probe.c_str(),
              reservoir_convergence::kTolPhiOverSpan,
              local_converged ? "EINGEHALTEN" : "NICHT EINGEHALTEN");
  std::printf("  davon am Meniskus (E_z auf der Achse und |E| bei zwei Bohrungsradien): "
              "%.3e -> %s\n",
              last_meniscus_step, meniscus_converged ? "eingehalten" : "NICHT eingehalten");

  // --- mesh convergence at the reference reservoir --------------------------
  Real mesh_phi = 0.0, mesh_E = 0.0;
  {
    std::vector<DielectricSolution> levels;
    std::FILE* f = std::fopen((outdir + "/convergence_mesh.csv").c_str(), "w");
    std::fprintf(f, "# Netzkonvergenz bei festem Plenum (plenum_c).  Sie trennt den "
                    "Diskretisierungsfehler\n# von der Wirkung der Vorratsgroesse.\n");
    std::fprintf(f, "level,size_scale,nodes,nr,nz,residual_C,Q_liquid_C,Q_extractor_C,"
                    "Ez_axis_surface_V_per_m");
    for (const Probe& p : sol[0].probes)
      std::fprintf(f, ",phi_%s_V,E_%s_V_per_m", p.name.c_str(), p.name.c_str());
    std::fprintf(f, "\n");
    for (int L = 0; L <= max_level; ++L) {
      levels.push_back(solve_dielectric(setup_of(variants[3], L)));
      const DielectricSolution& s = levels.back();
      std::fprintf(f, "%d,%.9g,%lld,%lld,%lld,%.9e,%.9e,%.9e,%.9e", L, mesh_level_scale(L),
                   static_cast<long long>(s.fem.n_nodes),
                   static_cast<long long>(s.mesh.grid.nr),
                   static_cast<long long>(s.mesh.grid.nz), s.fem.residual_inf, s.Q_emitter,
                   s.Q_extractor, Ez_surface_axis(s));
      for (std::size_t p = 0; p < s.probes.size(); ++p)
        std::fprintf(f, ",%.9e,%.9e", s.phi_probe[p], s.Emag_probe[p]);
      std::fprintf(f, "\n");
      std::printf("  Netzstufe %d: %lld Knoten, Q_liquid = %.6e C\n", L,
                  static_cast<long long>(s.fem.n_nodes), s.Q_emitter);
    }
    std::fclose(f);
    const std::size_t a = levels.size() - 2, b = levels.size() - 1;
    for (std::size_t p = 0; p < levels[0].probes.size(); ++p) {
      mesh_phi = std::max(mesh_phi, std::abs(levels[b].phi_probe[p] - levels[a].phi_probe[p]) /
                                        std::abs(S.applied_span()));
      mesh_E = std::max(mesh_E, rel_change(levels[a].Emag_probe[p], levels[b].Emag_probe[p]));
    }
  }

  // --- artefacts of the reference variant ----------------------------------
  const std::size_t k_ref = 3;   // plenum_c
  sol[k_ref].mesh.write_csv(outdir);
  sol[k_ref].write_csv(outdir);
  lib.write_csv(outdir + "/materials_library.csv");
  S.materials.write_csv(outdir + "/materials.csv");

  {
    std::FILE* f = std::fopen((outdir + "/variant_outlines.csv").c_str(), "w");
    std::fprintf(f, "# massstaebliche Meridianumrisse aller verglichenen Varianten\n");
    std::fprintf(f, "variant,name,i,r_m,z_m\n");
    for (std::size_t k = 0; k < variants.size(); ++k)
      append_outline(f, variants[k].tag, sol[k].mesh);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/reference_surface_fields.csv").c_str(), "w");
    std::fprintf(f, "# einseitiges E_z unmittelbar ueber der ebenen Fluessigkeitsoberflaeche "
                    "(z = 0+),\n# je Variante.  Kantennahe Zellen sind ausgeschlossen: die "
                    "unverrundete\n# Austrittskante hat dort ein divergierendes Feld, das der "
                    "Elementgroesse folgt.\n");
    std::fprintf(f, "variant,r_m,Ez_V_per_m\n");
    for (std::size_t k = 0; k < variants.size(); ++k)
      for (std::size_t q = 0; q < sol[k].surface_r.size(); ++q)
        std::fprintf(f, "%s,%.9e,%.9e\n", variants[k].tag.c_str(), sol[k].surface_r[q],
                     sol[k].surface_Ez[q]);
    std::fclose(f);
  }

  // --- identical near-field windows for every variant -----------------------
  {
    const Real L = S.geometry.device.extraction_distance;
    const Window near{0.0, 1.6 * L, 1.15 * S.geometry.base_z(),
                      1.25 * (L + S.geometry.device.extractor_thickness), 240, 240};
    for (std::size_t k = 0; k < variants.size(); ++k)
      write_grid_csv(outdir + "/field_near_" + variants[k].tag + ".csv", sol[k], near);
    std::FILE* f = std::fopen((outdir + "/field_window.csv").c_str(), "w");
    std::fprintf(f, "name,r0_m,r1_m,z0_m,z1_m\n");
    std::fprintf(f, "near,%.9e,%.9e,%.9e,%.9e\n", near.r0, near.r1, near.z0, near.z1);
    std::fclose(f);
  }

  // --- report ---------------------------------------------------------------
  {
    const DielectricSolution& R = sol[k_ref];
    std::FILE* f = std::fopen((outdir + "/report.txt").c_str(), "w");
    std::fprintf(f, "P2c -- Entkopplung von Geraetegeometrie und Fluessigkeitsvorrat\n");
    std::fprintf(f, "==============================================================\n\n");

    std::fprintf(f, "WAS KORRIGIERT WURDE\n");
    std::fprintf(f, "  Bisher verschob ein einziger Wert (liquid_feed_z) gleichzeitig die "
                    "Laenge der\n  leitfaehigen Fluessigkeitssaeule, die Laenge des "
                    "dielektrischen Rueckteils und die\n  rueckwaertige Geraetegeometrie. "
                    "Die daraus berichtete \"Konvergenz gegen die Lage\n  der Zulaufgrenze\" "
                    "war deshalb irrefuehrend: variiert wurde die Geometrie der\n  "
                    "Hochspannungselektrode.  Die Lage eines vollstaendig eingetauchten "
                    "Kontakts ist im\n  Modell des ideal leitfaehigen Fluessigkeitskoerpers "
                    "dagegen irrelevant und wird\n  gar nicht geometrisch dargestellt.\n\n");

    std::fprintf(f, "GETRENNTE PARAMETER\n");
    std::fprintf(f, "  vorderer Emitter        : phi_1 = %.4g m, phi_2 = %.4g m, "
                    "phi_3 = %.4g m, H = %.4g m\n",
                 S.geometry.device.phi_1, S.geometry.device.phi_2, S.geometry.device.phi_3,
                 S.geometry.device.emitter_height);
    std::fprintf(f, "  dielektrischer Grundkoerper: Dicke %.4g m -> Rueckflaeche bei "
                    "z = %.4g m   VORLAEUFIG\n",
                 S.geometry.base_plate_thickness, S.geometry.base_z());
    std::fprintf(f, "  fester Zulaufkanal      : Radius %.4g m, Laenge %.4g m   VORLAEUFIG\n",
                 R.mesh.r_channel, S.geometry.feed_channel_length);
    std::fprintf(f, "  Vorratsgeometrie        : Plenumradius, Plenumtiefe, Wandstaerke "
                    "%.4g m, Fuellstand\n", S.geometry.plenum_wall_thickness);
    std::fprintf(f, "  Extraktor               : L = %.4g m, D_a = %.4g m, t = %.4g m, "
                    "R_aussen = %.4g m\n\n",
                 S.geometry.device.extraction_distance,
                 S.geometry.device.extractor_aperture_diameter,
                 S.geometry.device.extractor_thickness,
                 S.geometry.device.extractor_outer_radius);

    std::fprintf(f, "NACHWEIS DER FESTEN FRONTGEOMETRIE\n");
    std::fprintf(f, "  %lld Gitterzeilen und %lld Zellen je Variante verglichen, von der "
                    "Rueckflaeche des\n  Grundkoerpers bis zum oberen Domaenenrand und bis "
                    "zum Extraktoraussenradius:\n  groesste Knotenabweichung %.1e m radial, "
                    "%.1e m axial; %lld Materialabweichungen.\n  Das ist bitgenaue "
                    "Gleichheit, nicht Uebereinstimmung innerhalb einer Toleranz.\n\n",
                 static_cast<long long>(rows_compared), static_cast<long long>(cells_compared),
                 worst_dr, worst_dz, static_cast<long long>(cell_mismatch));

    std::fprintf(f, "GERECHNETE VARIANTEN (Netzstufe %d)\n", ref_level);
    std::fprintf(f, "  %-14s %-46s %12s %12s %12s\n", "Kennung", "Beschreibung",
                 "V_fluid [m^3]", "E_z(0,0+) [V/m]", "Q_fluid [C]");
    for (std::size_t k = 0; k < variants.size(); ++k)
      std::fprintf(f, "  %-14s %-46s %12.5e %12.6g %12.5e\n", variants[k].tag.c_str(),
                   variants[k].label.c_str(), sol[k].mesh.analytic_volume_of(Region::Liquid),
                   Ez_surface_axis(sol[k]), sol[k].Q_emitter);
    std::fprintf(f, "\n");

    std::fprintf(f, "LOKALES EXTRAKTIONSFELD GEGEN VORRATSGROESSE\n");
    std::fprintf(f, "  Grenzen, vor der Messung festgelegt: phi %.1e der Spannweite, "
                    "|E| %.1e relativ.\n",
                 reservoir_convergence::kTolPhiOverSpan,
                 reservoir_convergence::kTolFieldRelative);
    std::fprintf(f, "  Beurteilt wird das MAXIMUM ueber alle kantenfernen Sondenpunkte, "
                    "nicht eine nach\n  der Messung ausgewaehlte Groesse.\n");
    std::fprintf(f, "  Letzte Vergroesserungsstufe (%s -> %s):\n"
                    "    phi %.3e der Spannweite (schlechtester Punkt: %s)\n"
                    "    |E| %.3e relativ        (schlechtester Punkt: %s)\n"
                    "    -> %s\n",
                 variants[3].tag.c_str(), variants[4].tag.c_str(), last_phi_step,
                 worst_phi_probe.c_str(), last_E_step, worst_E_probe.c_str(),
                 local_converged ? "EINGEHALTEN" : "NICHT EINGEHALTEN");
    std::fprintf(f, "  Groesster Schritt ueber alle Stufen : phi %.3e, |E| %.3e\n",
                 worst_phi_step, worst_E_step);
    std::fprintf(f, "  Getrennt ausgewiesen, weil es die Groesse ist, die ein Emissionsmodell "
                    "braucht:\n  am Meniskus selbst -- E_z auf der Achse unmittelbar ueber "
                    "der ebenen Oberflaeche\n  und |E| bei zwei Bohrungsradien -- betraegt "
                    "die letzte Aenderung %.3e relativ,\n  also %s die Grenze von %.1e. "
                    "Das ersetzt das Urteil oben NICHT.\n",
                 last_meniscus_step, meniscus_converged ? "innerhalb" : "ausserhalb",
                 reservoir_convergence::kTolFieldRelative);
    std::fprintf(f, "  Diskretisierungsfehler zum Vergleich (Netzstufe %d -> %d, festes "
                    "Plenum):\n    phi %.3e der Spannweite, |E| %.3e relativ.\n",
                 max_level - 1, max_level, mesh_phi, mesh_E);
    std::fprintf(f, "  Die Netz- und die Vorratswirkung sind damit getrennt bezifferbar; "
                    "die Frontnetze\n  der verglichenen Varianten sind ohnehin bitgleich, "
                    "der Diskretisierungsfehler faellt\n  in der Differenz also weitgehend "
                    "heraus.\n");
    if (!local_converged)
      std::fprintf(f, "  BEFUND: die Grenze wird NICHT erreicht.  Sie wird nicht gelockert, "
                      "es wird keine\n  kuenstliche Elektrode eingefuehrt und keine "
                      "\"Referenzgroesse\" stillschweigend\n  ausgewaehlt.  Jede berichtete "
                      "Zahl gilt fuer die genannte Vorratsgeometrie.\n");
    std::fprintf(f, "\n");

    std::fprintf(f, "GLOBALE LADUNGSGROESSEN -- weiterhin vorratsabhaengig, und das ist "
                    "richtig so\n");
    std::fprintf(f, "  %-14s %14s %14s %14s\n", "Kennung", "Q_fluid [C]", "Q_extr [C]",
                 "Summe [C]");
    for (std::size_t k = 0; k < variants.size(); ++k)
      std::fprintf(f, "  %-14s %14.6e %14.6e %14.6e\n", variants[k].tag.c_str(),
                   sol[k].Q_emitter, sol[k].Q_extractor, sol[k].Q_net);
    std::fprintf(f, "  Ein endlicher Leiter im offenen Raum ohne Rueckfuehrelektrode traegt "
                    "eine Ladung,\n  die mit seiner Groesse waechst.  Auf diese Groessen "
                    "wird deshalb keine Toleranz\n  gelegt; sie werden als Funktion der "
                    "Vorratsgeometrie berichtet.\n\n");

    std::fprintf(f, "WAS NICHT EINGEFUEHRT WURDE\n");
    std::fprintf(f, "  Keine leitfaehige Halterung, keine rueckwaertige Metallscheibe, keine "
                    "Basisplatte\n  auf Emitterpotential.  Der Vorratskoerper ist ein "
                    "Dielektrikum; das Audit zaehlt\n  %lld festgehaltene Knoten ohne "
                    "Fluessigkeitskontakt (Sollwert 0) und %lld auf den\n  benannten "
                    "Polymerflaechen (Sollwert 0).  Die P2a-Abschlussscheibe bleibt "
                    "abgelehnt.\n\n",
                 static_cast<long long>(R.audit.n_polymer_dirichlet),
                 static_cast<long long>(R.audit.n_named_surface_dirichlet));

    std::fprintf(f, "WOVON DIE ABSOLUTEN ZAHLEN ABHAENGEN\n");
    std::fprintf(f, "  Der Vernetzer verlangt, dass das Plenum RADIAL AUSSERHALB des "
                    "modellierten\n  Extraktors liegt (%.4g m), denn nur dann kann seine "
                    "Groesse keinen einzigen\n  Nahfeldknoten verschieben.  Der "
                    "Extraktoraussenradius ist selbst ein\n  BEISPIELWERT.  Mit dieser "
                    "Wahl ragt der Vorrat radial ueber die Elektrode\n  hinaus und wird "
                    "von ihr nicht abgeschirmt; das hebt das Potential in der\n  "
                    "Extraktionsstrecke deutlich an (phi in der Mitte: %.4g V mit "
                    "abgeschnittener\n  Saeule, %.4g V mit Plenum) und senkt E_z an der "
                    "Oberflaeche entsprechend\n  (%.4g V/m gegen %.4g V/m).  Das ist kein "
                    "Rechenfehler, sondern die Wirkung eines\n  grossen Leiters auf "
                    "Emitterpotential hinter einer kleinen Elektrode -- und es\n  zeigt, "
                    "wie wenig das frueher berichtete Saeulenergebnis ueber das Geraet "
                    "sagte.\n  Ein belegter Extraktoraussenradius und belegte "
                    "Vorratsabmessungen sind die\n  naechsten fehlenden Eingaben.\n\n",
                 S.geometry.device.extractor_outer_radius,
                 sol[0].phi_probe[static_cast<std::size_t>(k_mid)],
                 sol[k_ref].phi_probe[static_cast<std::size_t>(k_mid)],
                 Ez_surface_axis(sol[0]), Ez_surface_axis(sol[k_ref]));

    std::fprintf(f, "ERSATZGEOMETRIE, NICHT REKONSTRUKTION\n");
    std::fprintf(f, "  Die Dissertation zeigt in Abb. A.5 ein verbessertes, NICHT "
                    "rotationssymmetrisches\n  Reservoir von etwa 25 mm Hoehe sowie 10 mm "
                    "Aussen- und 8 mm Innenbreite.  Diese\n  Angaben sind hier "
                    "ausschliesslich zur Wahl der Groessenordnung benutzt worden --\n  "
                    "Millimeter, nicht Mikrometer.  Das achsensymmetrische Plenum ist eine "
                    "Ersatz-\n  geometrie und keine Rekonstruktion.  base_plate_thickness, "
                    "feed_channel_length,\n  Wandstaerke und Fuellstand sind VORLAEUFIGE "
                    "Beispielwerte, keine belegten Masse.\n\n");

    std::fprintf(f, "WAS HIER NICHT MODELLIERT IST\n");
    std::fprintf(f, "  %s\n", liquid_model::why_ideal_conductor_is_admissible_in_p2b());
    std::fprintf(f, "  Keine Oberflaechenspannung, keine Meniskusberechnung, keine "
                    "Stroemung, keine\n  Emission, keine Raumladung, keine 3D-Modellierung. "
                    "Die Ebene bei z = 0 ist die\n  anfaengliche ebene Fluessigkeitsoberflaeche."
                    "  Die Schwerkraft und die Benetzung im\n  Vorratsraum sind nicht "
                    "modelliert; der Fuellstand ist eine ebene Flaeche und ein\n  Parameter "
                    "des Ersatzmodells.\n");
    std::fclose(f);
  }

  {
    const DielectricSolution& R = sol[k_ref];
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_reservoir (P2c)\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "V_emitter_V=%.9e\n", S.V_emitter);
    std::fprintf(f, "V_extractor_V=%.9e\n", S.V_extractor);
    std::fprintf(f, "reference_level=%d\n", ref_level);
    std::fprintf(f, "max_level=%d\n", max_level);
    std::fprintf(f, "reference_size_scale=%.9g\n", mesh_level_scale(ref_level));
    std::fprintf(f, "reference_variant=%s\n", variants[k_ref].tag.c_str());
    std::fprintf(f, "nodes=%lld\n", static_cast<long long>(R.fem.n_nodes));
    std::fprintf(f, "nr=%lld\n", static_cast<long long>(R.mesh.grid.nr));
    std::fprintf(f, "nz=%lld\n", static_cast<long long>(R.mesh.grid.nz));
    std::fprintf(f, "base_plate_thickness_m=%.9e\n", S.geometry.base_plate_thickness);
    std::fprintf(f, "base_z_m=%.9e\n", S.geometry.base_z());
    std::fprintf(f, "feed_channel_radius_m=%.9e\n", R.mesh.r_channel);
    std::fprintf(f, "feed_channel_length_m=%.9e\n", S.geometry.feed_channel_length);
    std::fprintf(f, "plenum_wall_thickness_m=%.9e\n", S.geometry.plenum_wall_thickness);
    std::fprintf(f, "domain_radius_m=%.9e\n", S.geometry.device.domain_radius);
    std::fprintf(f, "domain_z_min_m=%.9e\n", S.geometry.device.domain_z_min);
    std::fprintf(f, "domain_z_max_m=%.9e\n", S.geometry.device.domain_z_max);
    std::fprintf(f, "front_rows_compared=%lld\n", static_cast<long long>(rows_compared));
    std::fprintf(f, "front_cells_compared=%lld\n", static_cast<long long>(cells_compared));
    std::fprintf(f, "front_max_dr_m=%.9e\n", worst_dr);
    std::fprintf(f, "front_max_dz_m=%.9e\n", worst_dz);
    std::fprintf(f, "front_cell_mismatches=%lld\n", static_cast<long long>(cell_mismatch));
    std::fprintf(f, "tol_phi_over_span=%.9e\n", reservoir_convergence::kTolPhiOverSpan);
    std::fprintf(f, "tol_E_rel=%.9e\n", reservoir_convergence::kTolFieldRelative);
    std::fprintf(f, "last_step_phi_over_span=%.9e\n", last_phi_step);
    std::fprintf(f, "last_step_E_rel=%.9e\n", last_E_step);
    std::fprintf(f, "last_step_worst_phi_probe=%s\n", worst_phi_probe.c_str());
    std::fprintf(f, "last_step_worst_E_probe=%s\n", worst_E_probe.c_str());
    std::fprintf(f, "last_step_meniscus_E_rel=%.9e\n", last_meniscus_step);
    std::fprintf(f, "meniscus_field_converged=%s\n", meniscus_converged ? "yes" : "no");
    std::fprintf(f, "worst_step_phi_over_span=%.9e\n", worst_phi_step);
    std::fprintf(f, "worst_step_E_rel=%.9e\n", worst_E_step);
    std::fprintf(f, "local_field_converged=%s\n", local_converged ? "yes" : "no");
    std::fprintf(f, "mesh_change_phi_over_span=%.9e\n", mesh_phi);
    std::fprintf(f, "mesh_change_E_rel=%.9e\n", mesh_E);
    std::fprintf(f, "Q_liquid_C=%.9e\n", R.Q_emitter);
    std::fprintf(f, "Q_extractor_C=%.9e\n", R.Q_extractor);
    std::fprintf(f, "Q_net_C=%.9e\n", R.Q_net);
    std::fprintf(f, "polymer_dirichlet_nodes=%lld\n",
                 static_cast<long long>(R.audit.n_polymer_dirichlet));
    std::fprintf(f, "named_surface_dirichlet_nodes=%lld\n",
                 static_cast<long long>(R.audit.n_named_surface_dirichlet));
    std::fprintf(f, "emitter_material=%s\n", S.materials.emitter_dielectric.name.c_str());
    std::fprintf(f, "emitter_eps_r=%.9g\n", S.materials.emitter_dielectric.relative_permittivity);
    std::fprintf(f, "emitter_eps_r_status=%s\n",
                 to_string(S.materials.emitter_dielectric.status));
    std::fprintf(f, "reservoir_material=%s\n", S.materials.reservoir_body.name.c_str());
    std::fprintf(f, "reservoir_eps_r=%.9g\n", S.materials.reservoir_body.relative_permittivity);
    std::fprintf(f, "state=dielectric electrostatics, rho_f=0, liquid an ideal conductor, "
                    "no surface tension, no meniscus, no flow, no emission, no space charge\n");
    std::fclose(f);
  }

  cfg.warn_about_unused(stdout, {"meta.", "fluid.", "beam.", "output.", "bem.", "wetting.",
                                 "species.", "feed."});
  std::printf("\ngeschrieben nach %s\n", outdir.c_str());
  return exit_code;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\n%s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "\nFehler: %s\n", e.what());
  return 2;
}
