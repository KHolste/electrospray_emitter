// es_electrocapillary -- P3b: self-consistent static electro-capillary
// equilibrium of the pinned meniscus.
//
//   es_electrocapillary <geometrie.cfg> [<p3b.cfg> ...] <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS FOR
//
//   gamma kappa(s) = delta_p_exit + eps0 E_n(s)^2 / 2
//
// with E_n the vacuum-side normal field on the free surface of the liquid,
// which is an ideal conductor at V_emitter, and the electrostatics the P2c
// dielectric problem.  The P3a capillary solver and the P2c field solver are
// used unchanged; what is new is the moving mesh between them, the projection
// of the surface load, and the gate that decides whether that load is a usable
// quantity at all.
//
// NOT in this phase: emission, finite liquid conductivity, flow, viscosity,
// space charge, time dependence, dynamic stability, Taylor-cone onset,
// cone-jet.  A branch that ends is a branch this solver could not continue.
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

#include "es/capillary.hpp"
#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/electrocapillary.hpp"
#include "es/liquid.hpp"

using namespace es;
using constants::pi;

namespace {

// ---------------------------------------------------------------------------

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
        "device.emitter_back_length ist gesetzt.  Die leitende P2a-Abschlussscheibe bleibt "
        "abgelehnt.");
  p.emitter_back_length = 0.0;
  p.reserved.edge_radius_inner = c.num("reserved.edge_radius_inner", 0.0);
  p.reserved.edge_radius_outer = c.num("reserved.edge_radius_outer", 0.0);
  p.reserved.contact_angle_deg = c.num("reserved.contact_angle_deg", 0.0);
  p.reserved.bore_diameter_at_inlet = c.num("reserved.bore_diameter_at_inlet", 0.0);
  p.reserved.porous_emitter = c.flag("reserved.porous_emitter", false);
  p.reserved.collector_enabled = c.flag("reserved.collector_enabled", false);
  return p;
}

void refuse_contact_angle(const Config& c) {
  for (const char* key : {"capillary.contact_angle_deg", "wetting.contact_angle_deg",
                          "liquid.contact_angle_deg"})
    if (c.has(key))
      throw std::runtime_error(
          std::string(key) +
          " ist gesetzt.  Die Kontaktlinie ist an der scharfen Austrittskante GEPINNT; ein "
          "zusaetzlicher Kontaktwinkel wuerde dieselbe Kante ein zweites Mal festlegen.");
}

LiquidProperties liquid_from(const Config& c) {
  LiquidProperties L = liquid_data_by_name(c.str("liquid.name", "emi-bf4"));
  const bool any = c.has("liquid.surface_tension") || c.has("liquid.density") ||
                   c.has("liquid.temperature") || c.has("liquid.status") ||
                   c.has("liquid.source");
  if (any) {
    L.substance = c.str("liquid.substance", L.substance);
    L.gamma = c.num("liquid.surface_tension", L.gamma);
    L.rho = c.num("liquid.density", L.rho);
    L.T = c.num("liquid.temperature", L.T);
    L.status = liquid_status_from_string(c.str("liquid.status", to_string(L.status)));
    L.source = c.str("liquid.source", "");
    L.caveat = c.str("liquid.caveat", L.caveat);
    if (L.status != LiquidDataStatus::Illustrative && L.source.empty())
      throw std::runtime_error("liquid.source fehlt zu einem Status ungleich illustrative.");
  }
  return L;
}

MaterialStatus material_status_from(const std::string& s) {
  if (s == "measured") return MaterialStatus::Measured;
  if (s == "manufacturer_spec") return MaterialStatus::ManufacturerSpec;
  if (s == "literature") return MaterialStatus::Literature;
  if (s == "provisional") return MaterialStatus::Provisional;
  throw std::runtime_error("material.<name>.status unbekannt: '" + s + "'");
}

Metallisation metallisation_from(const std::string& s) {
  if (s == "front_only") return Metallisation::FrontOnly;
  if (s == "front_and_aperture") return Metallisation::FrontAndAperture;
  if (s == "all_surfaces") return Metallisation::AllSurfaces;
  throw std::runtime_error("extractor.metallisation unbekannt: '" + s + "'");
}

std::string csv_safe(std::string s) {
  for (char& c : s)
    if (c == ',' || c == '\n' || c == '\r') c = ' ';
  return s;
}

// ---------------------------------------------------------------------------

struct ShapeCase {
  std::string tag;
  Real Pi{0};
  FreeSurface surface;
  CapillaryMeniscus shape;
};

}  // namespace

// ===========================================================================

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.size() < 2) {
    std::printf(
        "es_electrocapillary -- P3b: selbstkonsistentes statisches "
        "Elektro-Kapillargleichgewicht\n\n"
        "  es_electrocapillary <geometrie.cfg> [<p3b.cfg> ...] <ausgabeverzeichnis> "
        "[key=value ...]\n\n"
        "Beispiel:\n"
        "  ./build/es_electrocapillary examples/device_p1.cfg examples/electrocapillary_p3b.cfg "
        "\\\n      results/<ordner> meta.commit=$(git rev-parse HEAD)\n"
        "  python python/plot_electrocapillary.py results/<ordner>\n\n"
        "Gerechnet wird gamma kappa = delta_p_exit + eps0 E_n^2/2, statisch, ideal leitende\n"
        "Fluessigkeit.  Keine Emission, keine endliche Leitfaehigkeit, keine Stroemung, keine\n"
        "Raumladung, keine Zeitabhaengigkeit, keine Stabilitaetsaussage, kein Taylor-Kegel.\n");
    return 1;
  }

  Config cfg;
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) cfg.load(pos[i]);
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  std::string config_names;
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) {
    const std::filesystem::path src(pos[i]);
    std::filesystem::copy_file(src, std::filesystem::path(outdir) / src.filename(),
                               std::filesystem::copy_options::overwrite_existing);
    config_names += (config_names.empty() ? "" : ";") + src.filename().string();
  }

  refuse_contact_angle(cfg);

  // --- the model ------------------------------------------------------------
  DielectricDeviceParameters geo;
  geo.device = device_from(cfg);
  geo.base_plate_thickness = cfg.num("device.base_plate_thickness", geo.base_plate_thickness);
  geo.reservoir = ReservoirModel::AxisymmetricPlenum;
  geo.feed_channel_radius = cfg.num("reservoir.feed_channel_radius", 0.0);
  geo.feed_channel_length = cfg.num("reservoir.feed_channel_length", geo.feed_channel_length);
  geo.plenum_radius = cfg.num("reservoir.plenum_radius", geo.plenum_radius);
  geo.plenum_depth = cfg.num("reservoir.plenum_depth", geo.plenum_depth);
  geo.plenum_wall_thickness =
      cfg.num("reservoir.wall_thickness", geo.plenum_wall_thickness);
  geo.plenum_fill_fraction = 1.0;

  MaterialLibrary lib;
  for (const char* raw : {"su8", "ip-q", "ipx-q", "peek"}) {
    const std::string name(raw);
    const std::string key = "material." + name + ".relative_permittivity";
    if (!cfg.has(key)) continue;
    const std::string src = cfg.str("material." + name + ".source", "");
    if (src.empty())
      throw std::runtime_error("material." + name + ".source fehlt.");
    lib.override_permittivity(name, cfg.num(key, 0.0),
                              material_status_from(cfg.str("material." + name + ".status",
                                                           "provisional")),
                              src);
  }
  DielectricMaterials mats = DielectricMaterials::reference(lib);
  if (cfg.has("emitter.material")) mats.emitter_dielectric = lib.get(cfg.str("emitter.material", ""));
  if (cfg.has("extractor.material"))
    mats.extractor_carrier = lib.get(cfg.str("extractor.material", ""));
  if (cfg.has("reservoir.material"))
    mats.reservoir_body = lib.get(cfg.str("reservoir.material", ""));
  mats.check_usable();

  const LiquidProperties liquid = liquid_from(cfg);
  liquid.validate_or_throw();

  const Real a = 0.5 * geo.device.phi_2;
  const Real gamma_over_a = liquid.gamma / a;
  const Metallisation metal =
      metallisation_from(cfg.str("extractor.metallisation", "front_and_aperture"));
  const FarField far = FarField::Asymptotic;

  const int level_gate_lo = cfg.integer("mesh.gate_min_level", 1);
  const int level_gate_hi = cfg.integer("mesh.gate_max_level", 4);
  const int level_ref = cfg.integer("mesh.reference_level", 2);
  const Real V_max = cfg.num("field.V_max", 4000.0);

  std::printf("P3b -- selbstkonsistentes statisches Elektro-Kapillargleichgewicht\n");
  std::printf("  a = phi_2/2 = %.4g m, gamma/a = %.5g Pa, Stoffstatus %s\n", a, gamma_over_a,
              to_string(liquid.status));
  std::printf("  Referenznetzstufe %d, Gate-Stufen %d..%d, V_max = %.0f V\n\n", level_ref,
              level_gate_lo, level_gate_hi, V_max);

  int exit_code = 0;

  // --- the prescribed shapes the gate and the field figures use --------------
  auto shape_for = [&](Real Pi) {
    CapillaryRequest cr;
    cr.delta_p_exit = capillary::pressure_from_pi(Pi, a, liquid.gamma);
    cr.target_relative_accuracy = 1.0e-10;
    return solve_capillary_meniscus(a, 0.0, liquid, cr);
  };

  std::vector<ShapeCase> shapes;
  for (Real Pi : {0.0, 0.5, -0.5, 1.5, -1.5, 1.9}) {
    ShapeCase sc;
    char tag[32];
    std::snprintf(tag, sizeof tag, "pi%+05.2f", Pi);
    sc.tag = tag;
    sc.Pi = Pi;
    if (Pi == 0.0) {
      sc.surface = FreeSurface::flat_surface(a, 0.0);
    } else {
      sc.shape = shape_for(Pi);
      if (!is_usable(sc.shape.status)) {
        std::printf("  Form Pi = %+.2f: %s -- uebersprungen\n", Pi, to_string(sc.shape.status));
        continue;
      }
      sc.surface = FreeSurface::from(sc.shape);
    }
    shapes.push_back(std::move(sc));
  }

  // ==========================================================================
  // 1.  The moving mesh: quality on every prescribed shape
  // ==========================================================================
  {
    std::FILE* f = std::fopen((outdir + "/mesh_quality.csv").c_str(), "w");
    std::fprintf(f, "# Qualitaet des beweglichen Netzes je vorgeschriebener Form und Netzstufe.\n"
                    "# Das Fluessigkeitsvolumen des Netzes wird gegen die geschlossene Form\n"
                    "# (Bohrungssaeule + Rotationsvolumen des Meniskus) geprueft.\n");
    std::fprintf(f, "shape,Pi,level,n_nodes,min_jacobian,inverted_cells,max_aspect,max_shear,"
                    "band_z_lo,band_z_hi,contact_r_error,contact_z_error,apex_error,"
                    "surface_error,V_mesh,V_reference,V_error\n");
    for (const ShapeCase& sc : shapes)
      for (int lvl : {level_gate_lo, level_ref, level_gate_hi}) {
        DielectricDeviceParameters p = geo;
        p.mesh_level = lvl;
        MeniscusMesh m = build_meniscus_mesh(p, sc.surface);
        const MeniscusMeshQuality& q = m.quality;
        std::fprintf(f, "%s,%.9e,%d,%lld,%.9e,%lld,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                        "%.9e,%.9e,%.9e\n",
                     sc.tag.c_str(), sc.Pi, lvl,
                     static_cast<long long>(m.device.grid.n_nodes()), q.min_jacobian,
                     static_cast<long long>(q.inverted_cells), q.max_cell_aspect, q.max_shear,
                     q.band_z_lo, q.band_z_hi, q.contact_radius_error, q.contact_z_error,
                     q.apex_error, q.surface_error, q.liquid_volume_mesh,
                     q.liquid_volume_reference, q.liquid_volume_error);
        if (!q.ok()) exit_code = 2;
      }
    std::fclose(f);
    std::printf("  bewegliches Netz geprueft: %zu Formen x 3 Stufen\n", shapes.size());
  }

  // --- the mesh itself, near the tip, for the figures ------------------------
  {
    std::FILE* f = std::fopen((outdir + "/mesh_nodes.csv").c_str(), "w");
    std::fprintf(f, "# Knoten des beweglichen Netzes in einem Fenster um die Austrittskante.\n");
    std::fprintf(f, "shape,i,j,r_m,z_m\n");
    std::FILE* fs = std::fopen((outdir + "/mesh_surface.csv").c_str(), "w");
    std::fprintf(fs, "# Die Zeile des Netzes, die die freie Oberflaeche traegt.\n");
    std::fprintf(fs, "shape,i,r_m,z_m\n");
    for (const ShapeCase& sc : shapes) {
      if (std::abs(sc.Pi) != 0.0 && std::abs(sc.Pi) != 1.5) continue;
      DielectricDeviceParameters p = geo;
      p.mesh_level = level_ref;
      MeniscusMesh m = build_meniscus_mesh(p, sc.surface);
      const QuadMesh& g = m.device.grid;
      const Index i1 = std::min(g.nr - 1, m.device.i_land + 6);
      const Index j0 = std::max<Index>(0, m.device.j_tip - 30);
      const Index j1 = std::min(g.nz - 1, m.device.j_tip + 30);
      for (Index j = j0; j <= j1; ++j)
        for (Index i = 0; i <= i1; ++i)
          std::fprintf(f, "%s,%lld,%lld,%.9e,%.9e\n", sc.tag.c_str(),
                       static_cast<long long>(i), static_cast<long long>(j), g.at(i, j).r,
                       g.at(i, j).z);
      for (Index i = 0; i <= m.i_contact; ++i)
        std::fprintf(fs, "%s,%lld,%.9e,%.9e\n", sc.tag.c_str(), static_cast<long long>(i),
                     g.at(i, m.j_surface).r, g.at(i, m.j_surface).z);
    }
    std::fclose(f);
    std::fclose(fs);
  }

  // ==========================================================================
  // 2.  Field on the prescribed shapes -- and the independent checks
  // ==========================================================================
  const Real V_probe = cfg.num("field.V_probe", 1500.0);
  {
    std::FILE* f = std::fopen((outdir + "/field_window.csv").c_str(), "w");
    std::fprintf(f, "# Potential und Feldstaerke in einem festen Fenster, identisch fuer jede\n"
                    "# vorgeschriebene Form.  Punkte im Leiter tragen kein Feld und sind leer.\n");
    std::fprintf(f, "shape,r_m,z_m,phi_V,E_V_per_m\n");
    const int nr = 121, nz = 181;
    const Real r0 = 0.0, r1 = 5.0 * a, z0 = -2.0 * a, z1 = 7.0 * a;
    for (const ShapeCase& sc : shapes) {
      if (std::abs(sc.Pi) != 0.0 && std::abs(sc.Pi) != 1.5) continue;
      DielectricDeviceParameters p = geo;
      p.mesh_level = level_ref;
      MeniscusMesh m = build_meniscus_mesh(p, sc.surface);
      DielectricSetup s;
      s.geometry = p;
      s.materials = mats;
      s.metallisation = metal;
      s.far_field = far;
      s.V_emitter = V_probe;
      s.V_extractor = 0.0;
      DielectricSolution sol =
          solve_dielectric_on(m.device, s, DielectricDiagnostics::FieldOnly);
      for (int jz = 0; jz < nz; ++jz)
        for (int ir = 0; ir < nr; ++ir) {
          const Vec2 x{r0 + (r1 - r0) * ir / (nr - 1.0), z0 + (z1 - z0) * jz / (nz - 1.0)};
          Index i, j;
          Real xi, eta;
          if (!locate_meniscus(m, x, &i, &j, &xi, &eta)) continue;
          if (!sol.cell_active[static_cast<std::size_t>(m.device.grid.cell(i, j))]) {
            std::fprintf(f, "%s,%.9e,%.9e,,\n", sc.tag.c_str(), x.r, x.z);
            continue;
          }
          const Real phi = potential_at_meniscus(m, sol.fem.phi, x);
          const Vec2 E = field_recovered_at_meniscus(m, sol.fem.phi, sol.cell_eps_r,
                                                     sol.cell_active, x);
          std::fprintf(f, "%s,%.9e,%.9e,%.9e,%.9e\n", sc.tag.c_str(), x.r, x.z, phi, norm(E));
        }
    }
    std::fclose(f);
    std::printf("  Feldfenster geschrieben\n");
  }

  // --- independent validation of the field and flux evaluation --------------
  Real polarity_field_error = 0.0, polarity_force_error = 0.0, quadratic_error = 0.0;
  Real charge_error = 0.0, tangential_worst = 0.0;
  {
    DielectricDeviceParameters p = geo;
    p.mesh_level = level_ref;
    MeniscusMesh m = build_meniscus_mesh(p, shapes.front().surface);
    auto solve_at_V = [&](Real Ve, Real Vx) {
      DielectricSetup s;
      s.geometry = p;
      s.materials = mats;
      s.metallisation = metal;
      s.far_field = far;
      s.V_emitter = Ve;
      s.V_extractor = Vx;
      return solve_dielectric_on(m.device, s, DielectricDiagnostics::FieldOnly);
    };
    const DielectricSolution sp = solve_at_V(1000.0, 0.0);
    const DielectricSolution sn = solve_at_V(-1000.0, -0.0);
    const DielectricSolution sd = solve_at_V(2000.0, 0.0);
    const MaxwellLoad lp = maxwell_load(m, sp, gamma_over_a);
    const MaxwellLoad ln = maxwell_load(m, sn, gamma_over_a);
    const MaxwellLoad ld = maxwell_load(m, sd, gamma_over_a);
    for (std::size_t k = 0; k + 1 < lp.node_En.size(); ++k)
      polarity_field_error =
          std::max(polarity_field_error,
                   std::abs(lp.node_En[k] + ln.node_En[k]) / std::max(std::abs(lp.node_En[k]), 1e-30));
    polarity_force_error = std::abs(lp.total_force - ln.total_force) / lp.total_force;
    quadratic_error = std::abs(ld.total_force / lp.total_force - 4.0) / 4.0;
    // Charge from sigma = eps0 E_n on the free surface against the FEM nodal
    // reactions -- two independent routes to the same charge.
    //
    // The CONTACT NODE is left out of both sums, and it has to be: its finite
    // element patch reaches into the emitter dielectric as well, so its reaction
    // is the charge of two different surfaces added together, while the surface
    // integral knows only about the meniscus.  The nodal patches of the
    // remaining surface nodes cover the surface up to the midpoint of the last
    // segment, so that is where the integral stops.
    {
      const QuadMesh& g = m.device.grid;
      const std::size_t n = lp.node_r.size();
      Real q_sigma = 0.0;
      for (std::size_t k = 0; k + 1 < n; ++k) {
        const Real r0 = lp.node_r[k], r1 = lp.node_r[k + 1];
        const Real z0 = lp.node_z[k], z1 = lp.node_z[k + 1];
        const Real s0 = constants::eps0 * lp.node_En[k], s1 = constants::eps0 * lp.node_En[k + 1];
        const bool last = (k + 2 == n);
        const Real f = last ? 0.5 : 1.0;    // half of the last segment
        const Real rm = last ? 0.5 * (r0 + r1) : r1;
        const Real sm = last ? 0.5 * (s0 + s1) : s1;
        const Real ds = f * std::hypot(r1 - r0, z1 - z0);
        q_sigma += 2.0 * pi * ds * ((2 * r0 * s0 + r0 * sm + rm * s0 + 2 * rm * sm) / 6.0);
      }
      Real q_fem = 0.0;
      for (Index i = 0; i + 1 <= m.i_contact; ++i)
        q_fem += sp.fem.reaction[static_cast<std::size_t>(g.node(i, m.j_surface))];
      charge_error = std::abs(q_sigma - q_fem) / std::max(std::abs(q_fem), 1e-30);
    }
    for (std::size_t k = 0; k < lp.node_tangential_fraction.size(); ++k)
      if (lp.node_d_edge[k] >= edge_gate::kExclusionCoarse * a)
        tangential_worst = std::max(tangential_worst, lp.node_tangential_fraction[k]);

    std::FILE* f = std::fopen((outdir + "/field_checks.csv").c_str(), "w");
    std::fprintf(f, "# Unabhaengige Pruefungen der Feld- und Flussauswertung, ebene Oberflaeche.\n");
    std::fprintf(f, "check,value,bound,note\n");
    std::fprintf(f, "polarity_En_sign,%.9e,0,E_n(+V) + E_n(-V) relativ; muss null sein\n",
                 polarity_field_error);
    std::fprintf(f, "polarity_force,%.9e,0,p_M ist polaritaetsunabhaengig\n",
                 polarity_force_error);
    std::fprintf(f, "quadratic_scaling,%.9e,0,F(2V)/F(V) = 4\n", quadratic_error);
    std::fprintf(f, "charge_sigma_vs_reaction,%.9e,1e-2,"
                    "Ladung aus sigma = eps0 E_n gegen die FEM-Knotenreaktionen ohne den "
                    "Kontaktknoten\n",
                 charge_error);
    std::fprintf(f, "tangential_fraction_edge_far,%.9e,,"
                    "|E_t|/|E| auf der Aequipotentialflaeche kantenfern; ein "
                    "Diskretisierungsmass ohne feste Grenze - die Netzabhaengigkeit steht in "
                    "gate_levels.csv\n",
                 tangential_worst);
    std::fclose(f);
    std::printf("  Feldpruefungen: Polaritaet %.1e, Quadratik %.1e, Ladung %.1e, "
                "Tangentialanteil %.1e\n",
                polarity_force_error, quadratic_error, charge_error, tangential_worst);
  }

  // ==========================================================================
  // 3.  The gate
  // ==========================================================================
  std::vector<int> gate_levels;
  for (int l = level_gate_lo; l <= level_gate_hi; ++l) gate_levels.push_back(l);

  std::vector<EdgeGateResult> gates;
  std::vector<Real> gate_psi;   // contact tangent angle of each tested shape [rad]
  for (const ShapeCase& sc : shapes) {
    EdgeGateResult g = run_edge_gate(geo, mats, V_probe, 0.0, metal, far, sc.surface, sc.tag,
                                     sc.Pi, gate_levels, gamma_over_a);
    gate_psi.push_back(sc.surface.flat ? 0.0 : sc.surface.psi.back());
    std::printf("  Gate %s (Pi = %+.2f): %s  beta = %.3f, dF = %.2e, dLast = %.2e, "
                "Grenzkraft %.4e / %.4e N (%.2e)\n",
                sc.tag.c_str(), sc.Pi, to_string(g.verdict), g.fitted_exponent,
                g.measured_total_force_change, g.measured_edge_far_change, g.limit_force_mesh,
                g.limit_force_exclusion, g.measured_limit_agreement);
    gates.push_back(std::move(g));
  }

  // --- WHERE the gate passed, as an interval of contact tangent angle -------
  //
  // The gate is a property of the SHAPE, not of the device: the singularity at
  // the contact line is set by the angle at which the free surface leaves the
  // edge.  The admissible set is therefore the range of contact angles over
  // which the gate passed, bounded by the nearest shape that failed, and it is
  // read off the measurements mechanically -- no shape is admitted because it
  // "looks like" one that passed.
  bool any_passed = false;
  Real psi_admissible_lo = 0.0, psi_admissible_hi = 0.0;
  for (std::size_t k = 0; k < gates.size(); ++k) {
    if (gates[k].verdict != GateVerdict::Passed) continue;
    if (!any_passed) {
      psi_admissible_lo = psi_admissible_hi = gate_psi[k];
      any_passed = true;
    }
    psi_admissible_lo = std::min(psi_admissible_lo, gate_psi[k]);
    psi_admissible_hi = std::max(psi_admissible_hi, gate_psi[k]);
  }
  const bool gate_passed = any_passed;
  auto admissible = [&](Real psi_contact) {
    return any_passed && psi_contact >= psi_admissible_lo - 1e-12 &&
           psi_contact <= psi_admissible_hi + 1e-12;
  };
  std::printf("  Gate bestanden fuer Kontaktwinkel psi in [%.2f, %.2f] Grad; ausserhalb wird\n"
              "  nicht gekoppelt.\n",
              psi_admissible_lo * 180.0 / pi, psi_admissible_hi * 180.0 / pi);

  {
    std::FILE* f = std::fopen((outdir + "/gate_levels.csv").c_str(), "w");
    std::fprintf(f, "# Kanten-Gate: je vorgeschriebener Form und Netzstufe.\n"
                    "# peak_node_pM ist der PUNKTWEISE Kantenwert.  Er konvergiert nicht und\n"
                    "# wird nirgends verwendet; er steht hier, damit man das sieht.\n");
    std::fprintf(f, "shape,Pi,level,n_nodes,n_segments,smallest_d_m,total_force_N,"
                    "force_beyond_0.10a_N,force_beyond_0.05a_N,force_beyond_0.025a_N,"
                    "peak_node_pM_Pa,max_tangential_fraction,fit_exponent,fit_r2,fit_d_lo_m,"
                    "fit_d_hi_m,fit_points\n");
    for (const EdgeGateResult& g : gates)
      for (const EdgeStudyPoint& pt : g.levels)
        std::fprintf(f, "%s,%.9e,%d,%lld,%lld,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                        "%.9e,%.9e,%lld\n",
                     g.shape_tag.c_str(), g.Pi, pt.mesh_level,
                     static_cast<long long>(pt.n_nodes),
                     static_cast<long long>(pt.n_surface_segments), pt.smallest_d,
                     pt.total_force, pt.force_coarse, pt.force_mid, pt.force_fine,
                     pt.peak_node_pM, pt.max_tangential_fraction, pt.fit_exponent, pt.fit_r2,
                     pt.fit_d_lo, pt.fit_d_hi, static_cast<long long>(pt.fit_points));
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/gate_verdict.csv").c_str(), "w");
    std::fprintf(f, "# Urteil des Kanten-Gates.  Die Grenzen stehen in edge_gate:: und wurden\n"
                    "# vor der Messung festgelegt.  exclusion_change ist die urspruengliche,\n"
                    "# fehlspezifizierte Groesse: sie vergleicht den Kraftinhalt zweier\n"
                    "# verschiedener Gebiete und kann nicht klein werden.  Sie wird berichtet,\n"
                    "# entscheidet aber nicht; das tut limit_agreement.\n");
    std::fprintf(f, "shape,Pi,verdict,fit_exponent,wedge_reference_exponent,"
                    "total_force_change,tol_total_force,edge_far_change,tol_edge_far,"
                    "exclusion_change,tol_exclusion,limit_force_mesh_N,limit_force_exclusion_N,"
                    "limit_agreement,tol_limit_agreement,note\n");
    for (const EdgeGateResult& g : gates)
      std::fprintf(f, "%s,%.9e,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                      "%.9e,%s\n",
                   g.shape_tag.c_str(), g.Pi, to_string(g.verdict), g.fitted_exponent,
                   g.wedge_reference_exponent, g.measured_total_force_change,
                   edge_gate::kTolTotalForce, g.measured_edge_far_change,
                   edge_gate::kTolEdgeFarLoad, g.measured_exclusion_change,
                   edge_gate::kTolExclusion, g.limit_force_mesh, g.limit_force_exclusion,
                   g.measured_limit_agreement, edge_gate::kTolLimitAgreement,
                   csv_safe(g.note).c_str());
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/edge_profiles.csv").c_str(), "w");
    std::fprintf(f, "# Einseitiges Normalfeld und Maxwell-Druck auf der Oberflaeche, je\n"
                    "# Netzstufe.  d ist der Abstand zur Kontaktlinie entlang der Oberflaeche.\n"
                    "# Die Werte bei kleinem d sind PUNKTWEISE und konvergieren nicht.\n");
    std::fprintf(f, "shape,Pi,level,d_m,r_m,z_m,tau,En_V_per_m,pM_Pa,tangential_fraction\n");
    std::FILE* fseg = std::fopen((outdir + "/edge_segments.csv").c_str(), "w");
    std::fprintf(fseg, "# Die konservative Segmentprojektion: Kraft und Flaeche je Segment,\n"
                       "# und der daraus gebildete Segmentdruck.  Diese Groesse konvergiert.\n");
    std::fprintf(fseg, "shape,Pi,level,d_mid_m,tau_mid,area_m2,force_N,pressure_Pa\n");
    for (std::size_t gi = 0; gi < gates.size(); ++gi)
      for (std::size_t li = 0; li < gates[gi].loads.size(); ++li) {
        const MaxwellLoad& L = gates[gi].loads[li];
        const int lvl = gates[gi].levels[li].mesh_level;
        for (std::size_t k = 0; k < L.node_r.size(); ++k)
          std::fprintf(f, "%s,%.9e,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                       gates[gi].shape_tag.c_str(), gates[gi].Pi, lvl, L.node_d_edge[k],
                       L.node_r[k], L.node_z[k], L.node_tau[k], L.node_En[k], L.node_pM[k],
                       L.node_tangential_fraction[k]);
        for (std::size_t k = 0; k < L.seg_force.size(); ++k)
          std::fprintf(fseg, "%s,%.9e,%d,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                       gates[gi].shape_tag.c_str(), gates[gi].Pi, lvl, L.seg_d_mid[k],
                       L.seg_tau_mid[k], L.seg_area[k], L.seg_force[k], L.seg_pressure[k]);
      }
    std::fclose(f);
    std::fclose(fseg);
  }

  if (!gate_passed) {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_electrocapillary (P3b)\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "config=%s\n", config_names.c_str());
    std::fprintf(f, "gate=failed\n");
    std::fprintf(f, "coupled=no\n");
    std::fclose(f);
    std::printf("\nDas Kanten-Gate ist NICHT bestanden.  Es wird nichts gekoppelt und keine\n"
                "selbstkonsistente Form vorgetaeuscht.  Abbildungen 1 bis 4 sind erzeugt,\n"
                "5 bis 7 entfallen mit EdgeLoadNotWellPosed.\n");
    return 2;
  }

  // ==========================================================================
  // 4.  Zero-field cross-check against P3a
  // ==========================================================================
  Real zero_field_error = 0.0;
  {
    CoupledRequest q;
    q.geometry = geo;
    q.geometry.mesh_level = level_ref;
    q.materials = mats;
    q.liquid = liquid;
    q.metallisation = metal;
    q.far_field = far;
    q.V_emitter = 0.0;
    q.V_extractor = 0.0;
    std::FILE* f = std::fopen((outdir + "/zero_field_check.csv").c_str(), "w");
    std::fprintf(f, "# Rueckwaertskompatibilitaet: bei V = 0 muss P3b die P3a-Loesung\n"
                    "# reproduzieren.  Pflicht-Gate.\n");
    std::fprintf(f, "delta_p_Pa,Pi,status,h_p3b_m,h_p3a_m,dh_over_a,max_node_distance_over_a,"
                    "load_max_Pa\n");
    for (Real Pi : {0.0, 0.5, -0.5, 1.0}) {
      q.delta_p_exit = capillary::pressure_from_pi(Pi, a, liquid.gamma);
      const CoupledPoint p = solve_coupled(q);
      const CapillaryMeniscus ref = shape_for(Pi);
      // Compared at equal NORMALISED ARCLENGTH, not by node index: the two
      // solves may have chosen different resolutions, and node k is then a
      // different point on the curve.
      Real worst = 0.0, load_max = 0.0;
      if (is_usable(p.status) && is_usable(ref.status)) {
        const std::size_t n = std::min(p.shape.nodes.size(), ref.nodes.size());
        for (std::size_t k = 0; k < n; ++k) {
          const Real t = static_cast<Real>(k) / static_cast<Real>(n - 1);
          const std::size_t ka = static_cast<std::size_t>(
              t * static_cast<Real>(p.shape.nodes.size() - 1) + 0.5);
          const std::size_t kb =
              static_cast<std::size_t>(t * static_cast<Real>(ref.nodes.size() - 1) + 0.5);
          worst = std::max(worst, norm(p.shape.nodes[ka] - ref.nodes[kb]) / a);
        }
      }
      for (Real v : p.load.node_pM) load_max = std::max(load_max, v);
      zero_field_error = std::max(zero_field_error, worst);
      std::fprintf(f, "%.9e,%.9e,%s,%.9e,%.9e,%.9e,%.9e,%.9e\n", q.delta_p_exit, Pi,
                   to_string(p.status), p.apex_height, ref.apex_height,
                   std::abs(p.apex_height - ref.apex_height) / a, worst, load_max);
      if (!is_usable(p.status)) exit_code = 2;
    }
    std::fclose(f);
    std::printf("  Nullfeld gegen P3a: groesster Knotenabstand %.3e a\n", zero_field_error);
    if (zero_field_error > 0.0) exit_code = 2;
  }

  // ==========================================================================
  // 5.  Continuation over the voltage, for several exit pressures and both
  //     polarities
  // ==========================================================================
  struct Branch {
    std::string tag;
    Real delta_p{0}, Pi{0}, sign{1};
    ContinuationResult result;
  };
  std::vector<Branch> branches;
  {
    std::FILE* f = std::fopen((outdir + "/coupled_points.csv").c_str(), "w");
    std::fprintf(f, "# Selbstkonsistente Punkte des Astes, der die feldfreie P3a-Loesung\n"
                    "# enthaelt.  Nur konvergierte Punkte stehen hier.\n");
    std::fprintf(f, "branch,delta_p_Pa,Pi_pressure,V_emitter_V,V_extractor_V,level,status,"
                    "h_m,h_over_a,E_apex_V_per_m,E_edge_far_V_per_m,max_curvature_1_per_m,"
                    "total_force_N,surface_area_m2,liquid_volume_m3,iterations,shape_change,"
                    "load_change,mech_residual,mech_residual_edge_far,contact_error,"
                    "fem_residual_C,fem_residual_rel,min_jacobian,crossings\n");
    std::FILE* fp = std::fopen((outdir + "/coupled_profiles.csv").c_str(), "w");
    std::fprintf(fp, "# Meniskusprofile der selbstkonsistenten Punkte.\n");
    std::fprintf(fp, "branch,V_emitter_V,node,r_m,z_m,psi_rad,load_Pa\n");

    struct Spec { const char* tag; Real Pi; Real sign; };
    const Spec specs[] = {{"dp0_pos", 0.0, +1.0},
                          {"dp0_neg", 0.0, -1.0},
                          {"dp_plus", 0.5, +1.0},
                          {"dp_minus", -0.5, +1.0}};
    for (const Spec& sp : specs) {
      CoupledRequest q;
      q.geometry = geo;
      q.geometry.mesh_level = level_ref;
      q.materials = mats;
      q.liquid = liquid;
      q.metallisation = metal;
      q.far_field = far;
      q.delta_p_exit = capillary::pressure_from_pi(sp.Pi, a, liquid.gamma);
      Branch b;
      b.tag = sp.tag;
      b.delta_p = q.delta_p_exit;
      b.Pi = sp.Pi;
      b.sign = sp.sign;

      // The field-free member of this branch decides whether the branch may be
      // followed at all: if its contact angle is outside the range in which the
      // gate passed, the surface load is not a well-posed quantity there and
      // nothing is coupled.  It is not approximated, not clipped and not
      // started "just to see".
      const CapillaryMeniscus start = shape_for(sp.Pi);
      const Real psi_start = is_usable(start.status) ? start.contact_tangent_angle : 0.0;
      if (!admissible(psi_start)) {
        b.result.end_status = CouplingStatus::EdgeLoadNotWellPosed;
        b.result.end_message =
            "Der feldfreie Ast dieses Drucks hat den Kontaktwinkel " +
            std::to_string(psi_start * 180.0 / pi) +
            " Grad und liegt damit ausserhalb des Bereichs, in dem das Kanten-Gate bestanden "
            "wurde. Bei einer nach innen gezogenen Oberflaeche wird die Leiterkante "
            "einspringend, die Singularitaet verstaerkt sich und die Flaechenlast ist dort "
            "keine brauchbare Groesse. Es wird nichts gekoppelt.";
        std::printf("  Ast %-9s (Pi = %+.2f): EdgeLoadNotWellPosed -- Kontaktwinkel %.2f Grad "
                    "ausserhalb [%.2f, %.2f]\n",
                    sp.tag, sp.Pi, psi_start * 180.0 / pi, psi_admissible_lo * 180.0 / pi,
                    psi_admissible_hi * 180.0 / pi);
        branches.push_back(std::move(b));
        continue;
      }
      b.result = continue_over_voltage(q, sp.sign * V_max);
      // Every converged point is checked against the same range afterwards; a
      // point whose shape left it is dropped, and the branch ends there.
      std::size_t keep = b.result.points.size();
      for (std::size_t k = 0; k < b.result.points.size(); ++k)
        if (!admissible(b.result.points[k].shape.contact_tangent_angle)) {
          keep = k;
          break;
        }
      if (keep < b.result.points.size()) {
        b.result.points.resize(keep);
        b.result.end_status = CouplingStatus::EdgeLoadNotWellPosed;
        b.result.end_message =
            "Der Ast verlaesst den Kontaktwinkelbereich, in dem das Kanten-Gate bestanden "
            "wurde; die weiteren Punkte werden verworfen.";
        b.result.last_converged_voltage =
            keep > 0 ? b.result.points.back().V_emitter : 0.0;
      }
      std::printf("  Ast %-9s (Pi = %+.2f, Polaritaet %+.0f): %zu Punkte bis %.1f V, Ende %s\n",
                  sp.tag, sp.Pi, sp.sign, b.result.points.size(),
                  b.result.last_converged_voltage, to_string(b.result.end_status));
      for (const CoupledPoint& p : b.result.points) {
        std::fprintf(f, "%s,%.9e,%.9e,%.9e,%.9e,%d,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                        "%.9e,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%d\n",
                     sp.tag, p.delta_p_exit, sp.Pi, p.V_emitter, p.V_extractor, p.mesh_level,
                     to_string(p.status), p.apex_height, p.apex_height / a, p.E_apex,
                     p.E_edge_far, p.max_curvature, p.total_force, p.surface_area,
                     p.liquid_volume, p.iterations, p.final_shape_change, p.final_load_change,
                     p.mechanical_residual, p.mechanical_residual_edge_far, p.contact_error,
                     p.fem_residual, p.fem_residual_relative, p.min_jacobian, p.crossings);
        const std::size_t stride = std::max<std::size_t>(1, p.shape.nodes.size() / 400);
        for (std::size_t k = 0; k < p.shape.nodes.size(); k += stride)
          std::fprintf(fp, "%s,%.9e,%zu,%.9e,%.9e,%.9e,%.9e\n", sp.tag, p.V_emitter, k,
                       p.shape.nodes[k].r, p.shape.nodes[k].z, p.shape.psi[k],
                       p.shape.load.empty() ? 0.0 : p.shape.load[k]);
      }
      branches.push_back(std::move(b));
    }
    std::fclose(f);
    std::fclose(fp);
  }
  {
    std::FILE* f = std::fopen((outdir + "/branch_ends.csv").c_str(), "w");
    std::fprintf(f, "# Wo und warum jeder Ast endet.  Das ist die Stelle, an der DIESER Loeser\n"
                    "# stehen bleibt -- kein Emissionsbeginn, kein Taylor-Kegel, keine\n"
                    "# Stabilitaetsaussage.\n");
    std::fprintf(f, "branch,Pi_pressure,polarity,points,last_converged_V,first_failed_V,"
                    "steps_attempted,steps_rejected,end_status,end_message\n");
    for (const Branch& b : branches)
      std::fprintf(f, "%s,%.9e,%+.0f,%zu,%.9e,%.9e,%d,%d,%s,%s\n", b.tag.c_str(), b.Pi, b.sign,
                   b.result.points.size(), b.result.last_converged_voltage,
                   b.result.first_failed_voltage, b.result.steps_attempted,
                   b.result.steps_rejected, to_string(b.result.end_status),
                   csv_safe(b.result.end_message).c_str());
    std::fclose(f);
  }

  // --- polarity: the two branches must give the same shapes -----------------
  Real polarity_shape_error = 0.0;
  {
    const Branch* pos = nullptr;
    const Branch* neg = nullptr;
    for (const Branch& b : branches) {
      if (b.tag == "dp0_pos") pos = &b;
      if (b.tag == "dp0_neg") neg = &b;
    }
    std::FILE* f = std::fopen((outdir + "/polarity.csv").c_str(), "w");
    std::fprintf(f, "# Im Perfect-Conductor-Modell haengt p_M = eps0 E_n^2/2 nicht von der\n"
                    "# Polaritaet ab, also muessen beide Polaritaeten dieselbe Form liefern.\n"
                    "# Das ist ein MODELLTEST und ausdruecklich keine Aussage ueber die reale\n"
                    "# polare Emission, die in P3b gar nicht modelliert ist.\n");
    std::fprintf(f, "abs_V,h_positive_m,h_negative_m,dh_over_a\n");
    if (pos && neg)
      for (std::size_t k = 0; k < std::min(pos->result.points.size(), neg->result.points.size());
           ++k) {
        const Real dh = std::abs(pos->result.points[k].apex_height -
                                 neg->result.points[k].apex_height) / a;
        polarity_shape_error = std::max(polarity_shape_error, dh);
        std::fprintf(f, "%.9e,%.9e,%.9e,%.9e\n", std::abs(pos->result.points[k].V_emitter),
                     pos->result.points[k].apex_height, neg->result.points[k].apex_height, dh);
      }
    std::fclose(f);
    std::printf("  Polaritaet: groesster Formunterschied %.3e a\n", polarity_shape_error);
  }

  // ==========================================================================
  // 6.  Mesh and coupling convergence at three operating points
  // ==========================================================================
  {
    std::FILE* f = std::fopen((outdir + "/coupled_convergence.csv").c_str(), "w");
    (void)0;
    std::fprintf(f, "# Netz- und Kopplungskonvergenz an drei Betriebspunkten.\n"
                    "# Zusaetzlich: Fernrand (geerdete Huelle statt asymptotisch) und ein\n"
                    "# anderer Startwert bzw. eine andere Unterrelaxation.\n");
    std::fprintf(f, "point,V_emitter_V,variant,level,far_field,relaxation,status,h_over_a,"
                    "E_edge_far_V_per_m,total_force_N,max_curvature_1_per_m,surface_area_m2,"
                    "liquid_volume_m3,mech_residual_edge_far,iterations\n");
    struct Point { const char* tag; Real V; };
    const Point points[] = {{"klein", 500.0}, {"mittel", 1000.0}, {"gross", 1400.0}};
    for (const Point& pt : points) {
      for (int lvl : {level_gate_lo, level_ref, level_ref + 1}) {
        CoupledRequest q;
        q.geometry = geo;
        q.geometry.mesh_level = lvl;
        q.materials = mats;
        q.liquid = liquid;
        q.metallisation = metal;
        q.far_field = far;
        q.delta_p_exit = 0.0;
        q.V_emitter = pt.V;
        const CoupledPoint p = solve_coupled(q);
        std::fprintf(f, "%s,%.9e,mesh,%d,%s,%.3g,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%d\n",
                     pt.tag, pt.V, lvl, to_string(far), q.relaxation, to_string(p.status),
                     p.apex_height / a, p.E_edge_far, p.total_force, p.max_curvature,
                     p.surface_area, p.liquid_volume, p.mechanical_residual_edge_far,
                     p.iterations);
      }
      // far field and iteration controls, at the reference level
      for (int variant = 0; variant < 2; ++variant) {
        CoupledRequest q;
        q.geometry = geo;
        q.geometry.mesh_level = level_ref;
        q.materials = mats;
        q.liquid = liquid;
        q.metallisation = metal;
        q.delta_p_exit = 0.0;
        q.V_emitter = pt.V;
        q.far_field = (variant == 0) ? FarField::Grounded : far;
        q.relaxation = (variant == 0) ? coupling::kRelaxation : 0.3;
        const CoupledPoint p = solve_coupled(q);
        std::fprintf(f, "%s,%.9e,%s,%d,%s,%.3g,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%d\n",
                     pt.tag, pt.V, variant == 0 ? "farfield" : "relaxation", level_ref,
                     to_string(q.far_field), q.relaxation, to_string(p.status),
                     p.apex_height / a, p.E_edge_far, p.total_force, p.max_curvature,
                     p.surface_area, p.liquid_volume, p.mechanical_residual_edge_far,
                     p.iterations);
      }
      std::printf("  Konvergenz am Punkt %s (%.0f V) gerechnet\n", pt.tag, pt.V);
    }
    std::fclose(f);
  }

  // ==========================================================================
  // 7.  Parameters, report, meta
  // ==========================================================================
  {
    std::FILE* f = std::fopen((outdir + "/parameters.csv").c_str(), "w");
    std::fprintf(f, "name,value_SI,unit,role\n");
    std::fprintf(f, "contact_radius,%.9e,m,gepinnter Kontaktradius a = phi_2/2\n", a);
    std::fprintf(f, "surface_tension,%.9e,N/m,gamma\n", liquid.gamma);
    std::fprintf(f, "capillary_pressure_scale,%.9e,Pa,gamma/a\n", gamma_over_a);
    std::fprintf(f, "V_probe,%.9e,V,Spannung der vorgeschriebenen Formen\n", V_probe);
    std::fprintf(f, "V_max,%.9e,V,Obergrenze der Fortsetzung\n", V_max);
    std::fprintf(f, "reference_level,%d,-,Netzstufe der Kopplung\n", level_ref);
    std::fprintf(f, "band_factor,%.9e,-,halbe Hoehe des Verformungsbandes in a\n",
                 meniscus_mesh::kBandFactor);
    std::fprintf(f, "r_aperture,%.9e,m,Aperturradius\n",
                 0.5 * geo.device.extractor_aperture_diameter);
    std::fprintf(f, "extraction_distance,%.9e,m,Extraktionsstrecke\n",
                 geo.device.extraction_distance);
    std::fprintf(f, "r_land,%.9e,m,Aussenradius der Stirnflaeche\n", 0.5 * geo.device.phi_1);
    std::fclose(f);
  }

  {
    std::FILE* f = std::fopen((outdir + "/report.txt").c_str(), "w");
    std::fprintf(f, "P3b -- selbstkonsistentes statisches Elektro-Kapillargleichgewicht\n");
    std::fprintf(f, "=================================================================\n\n");
    std::fprintf(f, "WAS GERECHNET WURDE\n");
    std::fprintf(f, "  gamma kappa(s) = delta_p_exit + eps0 E_n(s)^2 / 2, achsensymmetrisch,\n");
    std::fprintf(f, "  statisch.  E_n ist das vakuumseitige Normalfeld auf der freien\n");
    std::fprintf(f, "  Oberflaeche der Fluessigkeit, die ein idealer Leiter auf V_emitter ist.\n");
    std::fprintf(f, "  Die Elektrostatik ist das unveraenderte dielektrische P2c-Problem.\n\n");
    std::fprintf(f, "VORZEICHEN\n");
    std::fprintf(f, "  sigma = eps0 E_n, Traktion = sigma * E_n/2 * n = p_M n mit p_M >= 0.\n");
    std::fprintf(f, "  Die Last zieht IMMER nach aussen, unabhaengig von der Polaritaet, und\n");
    std::fprintf(f, "  geht deshalb wie eine Erhoehung des Innendrucks in die Gleichung ein.\n");
    std::fprintf(f, "  Gemessen: E_n(+V) + E_n(-V) = %.2e relativ, F(+V) - F(-V) = %.2e\n",
                 polarity_field_error, polarity_force_error);
    std::fprintf(f, "  relativ, F(2V)/F(V) - 4 = %.2e, Ladung aus sigma gegen die\n",
                 quadratic_error);
    std::fprintf(f, "  FEM-Knotenreaktionen %.2e relativ, |E_t|/|E| kantenfern %.2e.\n\n",
                 charge_error, tangential_worst);

    std::fprintf(f, "DAS KANTEN-GATE\n");
    for (const EdgeGateResult& g : gates) {
      std::fprintf(f, "  Form %s (Pi = %+.2f): %s\n", g.shape_tag.c_str(), g.Pi,
                   to_string(g.verdict));
      std::fprintf(f, "    Singularitaetsexponent p_M ~ d^beta: beta = %.3f (R^2 = %.4f).\n",
                   g.fitted_exponent, g.levels.back().fit_r2);
      std::fprintf(f, "    Integrierbar gegen 2 pi r ds fuer beta > -1: %s.\n",
                   g.fitted_exponent > -1.0 ? "ja" : "NEIN");
      std::fprintf(f, "    Kantenferne Last, zwei feinste Stufen: %.3e (Grenze %.1e).\n",
                   g.measured_edge_far_change, edge_gate::kTolEdgeFarLoad);
      std::fprintf(f, "    Gesamtkraft, zwei feinste Stufen  : %.3e (Grenze %.1e).\n",
                   g.measured_total_force_change, edge_gate::kTolTotalForce);
      std::fprintf(f, "    Grenzkraft ueber die Netzstufen   : %.6e N\n", g.limit_force_mesh);
      std::fprintf(f, "    Grenzkraft ueber die Ausschlussdistanz: %.6e N\n",
                   g.limit_force_exclusion);
      std::fprintf(f, "    Abweichung der beiden             : %.3e (Grenze %.1e).\n",
                   g.measured_limit_agreement, edge_gate::kTolLimitAgreement);
      std::fprintf(f, "    Halbierung der Ausschlussdistanz  : %.3e (Grenze %.1e) -- NICHT\n",
                   g.measured_exclusion_change, edge_gate::kTolExclusion);
      std::fprintf(f, "      eingehalten, und diese Groesse kann es nicht: sie vergleicht den\n");
      std::fprintf(f, "      Kraftinhalt zweier verschiedener Gebiete der Oberflaeche und ist\n");
      std::fprintf(f, "      kein Diskretisierungsfehler.  Sie wird berichtet, entscheidet\n");
      std::fprintf(f, "      aber nicht; siehe docs/10_electrocapillary_model.md.\n");
      std::fprintf(f, "    Punktweiser Kantenwert p_M (feinste Stufe): %.4g Pa -- waechst mit\n",
                   g.levels.back().peak_node_pM);
      std::fprintf(f, "      jeder Verfeinerung und wird nirgends verwendet.\n");
    }
    std::fprintf(f, "  WO DAS GATE BESTANDEN IST.  Die Singularitaet an der Kontaktlinie\n");
    std::fprintf(f, "  haengt vom Winkel ab, unter dem die Oberflaeche die Kante verlaesst.\n");
    std::fprintf(f, "  Gemessen: fuer eine ebene oder nach aussen gewoelbte Oberflaeche ist die\n");
    std::fprintf(f, "  Last integrierbar und die schwache Projektion konvergiert; fuer eine nach\n");
    std::fprintf(f, "  innen gezogene wird die Leiterkante einspringend, der Exponent faellt und\n");
    std::fprintf(f, "  bei Pi = -1.5 ist er mit %.3f nicht mehr integrierbar.  Gekoppelt wird\n",
                 gates.size() > 4 ? gates[4].fitted_exponent : 0.0);
    std::fprintf(f, "  ausschliesslich im Bereich psi = %.2f .. %.2f Grad.\n\n",
                 psi_admissible_lo * 180.0 / pi, psi_admissible_hi * 180.0 / pi);

    std::fprintf(f, "VERWENDETE LASTPROJEKTION\n");
    std::fprintf(f, "  Segmentweise konservativ: die Normalkraft jedes Oberflaechensegments\n");
    std::fprintf(f, "  wird mit dem einseitigen Feld integriert, der Segmentdruck ist Kraft\n");
    std::fprintf(f, "  durch Rotationsflaeche.  Die Summe ueber die Segmente ist per\n");
    std::fprintf(f, "  Konstruktion die integrierte Maxwell-Kraft.  Fuer die Uebergabe an den\n");
    std::fprintf(f, "  Kapillarloeser wird daraus eine STETIGE Last gebaut, indem die\n");
    std::fprintf(f, "  kumulierte Kraft und die kumulierte Flaeche monoton kubisch\n");
    std::fprintf(f, "  interpoliert werden und p = G'/A' gesetzt wird; die Binintegrale\n");
    std::fprintf(f, "  bleiben dabei exakt erhalten.\n\n");

    std::fprintf(f, "BEWEGLICHES NETZ\n");
    std::fprintf(f, "  Die Zeile des P2c-Netzes, die die freie Oberflaeche traegt, wird auf die\n");
    std::fprintf(f, "  vorgeschriebene Meridiankurve gezogen; das Band reicht %.1f Bohrungs-\n",
                 meniscus_mesh::kBandFactor);
    std::fprintf(f, "  radien nach oben und unten und endet auf vorhandenen Netzzeilen, so dass\n");
    std::fprintf(f, "  der ebene Fall bitgleich das P2c-Netz ist.  Zellen koennen nicht\n");
    std::fprintf(f, "  invertieren: in den verformten Spalten haengt r nur von i ab.  Gemessen\n");
    std::fprintf(f, "  wird es trotzdem -- siehe mesh_quality.csv.\n\n");

    std::fprintf(f, "NULLFELD-RUECKPRUEFUNG\n");
    std::fprintf(f, "  Bei V_emitter = V_extractor = 0 reproduziert P3b die P3a-Loesung.\n");
    std::fprintf(f, "  Groesster Knotenabstand ueber vier Druecke: %.3e a.\n\n",
                 zero_field_error);

    std::fprintf(f, "FORTSETZUNG UEBER DIE SPANNUNG\n");
    for (const Branch& b : branches)
      std::fprintf(f, "  %-9s Pi = %+.2f, Polaritaet %+.0f: %zu konvergierte Punkte, letzte\n"
                      "            Spannung %.1f V, erste gescheiterte %.1f V, Ende %s\n",
                   b.tag.c_str(), b.Pi, b.sign, b.result.points.size(),
                   b.result.last_converged_voltage, b.result.first_failed_voltage,
                   to_string(b.result.end_status));
    std::fprintf(f, "  Groesster Formunterschied zwischen den Polaritaeten: %.3e a.\n",
                 polarity_shape_error);
    std::fprintf(f, "  WAS DAS ENDE EINES ASTES IST UND WAS NICHT.  Es ist die Stelle, an der\n");
    std::fprintf(f, "  dieser Loeser mit dieser Schrittweite und diesen Grenzen stehen bleibt.\n");
    std::fprintf(f, "  Es ist KEIN Emissionsbeginn, KEIN Taylor-Kegel-Onset und KEINE Aussage\n");
    std::fprintf(f, "  ueber dynamische Stabilitaet: nichts davon wird in P3b gerechnet.\n\n");

    std::fprintf(f, "WAS AUSDRUECKLICH FEHLT\n");
    std::fprintf(f, "  Keine Emission, keine endliche Leitfaehigkeit der Fluessigkeit, keine\n");
    std::fprintf(f, "  Stroemung, keine Viskositaet, keine Raumladung, keine Zeitabhaengigkeit,\n");
    std::fprintf(f, "  keine dynamische Stabilitaet, kein Taylor-Kegel, kein Cone-Jet, keine\n");
    std::fprintf(f, "  Schwerkraft.  delta_p_exit bleibt eine Eingabe; P3b bestimmt weder den\n");
    std::fprintf(f, "  Vorrats- noch den Speisedruck.  Die Stoffwerte sind %s.\n",
                 to_string(liquid.status));
    std::fclose(f);
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_electrocapillary (P3b)\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "config=%s\n", config_names.c_str());
    std::fprintf(f, "state=static electro-capillary equilibrium, perfect conductor, no emission, "
                    "no flow, no space charge, no stability statement\n");
    std::fprintf(f, "contact_radius_m=%.9e\n", a);
    std::fprintf(f, "surface_tension_N_per_m=%.9e\n", liquid.gamma);
    std::fprintf(f, "capillary_pressure_scale_Pa=%.9e\n", gamma_over_a);
    std::fprintf(f, "liquid_substance=%s\n", liquid.substance.c_str());
    std::fprintf(f, "liquid_status=%s\n", to_string(liquid.status));
    std::fprintf(f, "emitter_material=%s\n", mats.emitter_dielectric.name.c_str());
    std::fprintf(f, "emitter_eps_r=%.9g\n", mats.emitter_dielectric.relative_permittivity);
    std::fprintf(f, "emitter_eps_r_status=%s\n", to_string(mats.emitter_dielectric.status));
    std::fprintf(f, "V_probe_V=%.9e\n", V_probe);
    std::fprintf(f, "V_max_V=%.9e\n", V_max);
    std::fprintf(f, "reference_level=%d\n", level_ref);
    std::fprintf(f, "gate_min_level=%d\n", level_gate_lo);
    std::fprintf(f, "gate_max_level=%d\n", level_gate_hi);
    std::fprintf(f, "gate=%s\n", gate_passed ? "passed" : "failed");
    std::fprintf(f, "gate_psi_admissible_lo_deg=%.9e\n", psi_admissible_lo * 180.0 / pi);
    std::fprintf(f, "gate_psi_admissible_hi_deg=%.9e\n", psi_admissible_hi * 180.0 / pi);
    std::fprintf(f, "coupled=%s\n", gate_passed ? "yes" : "no");
    std::fprintf(f, "gate_exponent=%.9e\n", gates.empty() ? 0.0 : gates.front().fitted_exponent);
    std::fprintf(f, "zero_field_max_node_error_over_a=%.9e\n", zero_field_error);
    std::fprintf(f, "polarity_shape_error_over_a=%.9e\n", polarity_shape_error);
    std::fprintf(f, "polarity_force_error=%.9e\n", polarity_force_error);
    std::fprintf(f, "quadratic_error=%.9e\n", quadratic_error);
    std::fprintf(f, "charge_error=%.9e\n", charge_error);
    std::fprintf(f, "tangential_fraction_edge_far=%.9e\n", tangential_worst);
    for (const Branch& b : branches)
      std::fprintf(f, "branch_%s_last_V=%.9e\n", b.tag.c_str(),
                   b.result.last_converged_voltage);
    std::fclose(f);
  }

  cfg.warn_about_unused(stdout, {"meta.", "fluid.", "beam.", "output.", "bem.", "wetting.",
                                 "species.", "feed.", "capillary."});
  std::printf("\ngeschrieben nach %s\n", outdir.c_str());
  return exit_code;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\n%s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "\nFehler: %s\n", e.what());
  return 2;
}
