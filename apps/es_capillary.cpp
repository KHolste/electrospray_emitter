// es_capillary -- P3a: the static shape of the free liquid surface at the exit
// edge, held by surface tension alone.
//
//   es_capillary <geometrie.cfg> [<p3a.cfg> ...] <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS FOR
//
// One thing only: the first checkable capillary ground floor.  The bore is
// assumed full up to the sharp exit edge, the contact line is pinned there, and
// the shape follows from
//
//     gamma * (dpsi/ds + sin psi / r) = delta_p_exit,
//     delta_p_exit = p_liquid - p_vacuum.
//
// NOT in this phase, and not to be read into any number below: no electric
// field, no Maxwell pressure, no coupling to the electrostatic solver, no
// emission, no space charge, no flow, no time dependence, no stability, no
// Taylor cone, no cone-jet.  Gravity is not coupled either -- the Bond number
// that justifies leaving it out is computed and printed.
//
// The P2c reservoir geometry is untouched: it is not built, not meshed and not
// solved here.  The upstream state enters through delta_p_exit and through
// nothing else.
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
#include "es/device_geometry.hpp"
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
    throw std::runtime_error("device.extractor_outer_radius fehlt (Pflichtangabe der Geometrie).");
  p.extractor_outer_radius = c.num("device.extractor_outer_radius", 0.0);
  p.emitter_back_length = c.num("device.emitter_back_length", 0.0);
  p.reserved.edge_radius_inner = c.num("reserved.edge_radius_inner", 0.0);
  p.reserved.edge_radius_outer = c.num("reserved.edge_radius_outer", 0.0);
  p.reserved.contact_angle_deg = c.num("reserved.contact_angle_deg", 0.0);
  p.reserved.bore_diameter_at_inlet = c.num("reserved.bore_diameter_at_inlet", 0.0);
  p.reserved.porous_emitter = c.flag("reserved.porous_emitter", false);
  p.reserved.collector_enabled = c.flag("reserved.collector_enabled", false);
  return p;
}

/// A contact angle must not be smuggled in beside the pinned contact line.
/// DeviceGeometry::build already refuses reserved.contact_angle_deg; this
/// catches the same intent under the names a user would reach for here.
void refuse_contact_angle(const Config& c) {
  for (const char* key : {"capillary.contact_angle_deg", "wetting.contact_angle_deg",
                          "liquid.contact_angle_deg"}) {
    if (c.has(key))
      throw std::runtime_error(
          std::string(key) +
          " ist gesetzt.  Die Kontaktlinie ist an der scharfen Austrittskante GEPINNT; ein "
          "zusaetzlicher Kontaktwinkel wuerde dieselbe Kante ein zweites Mal festlegen und das "
          "Problem ueberbestimmen.  Pinning und Young-Winkel sind zwei einander ausschliessende "
          "Beschreibungen, und es wird keine davon stillschweigend bevorzugt.");
  }
}

LiquidProperties liquid_from(const Config& c) {
  LiquidProperties L = liquid_data_by_name(c.str("liquid.name", "emi-bf4"));
  const bool any_override = c.has("liquid.surface_tension") || c.has("liquid.density") ||
                            c.has("liquid.temperature") || c.has("liquid.status") ||
                            c.has("liquid.source");
  if (any_override) {
    L.substance = c.str("liquid.substance", L.substance);
    L.gamma = c.num("liquid.surface_tension", L.gamma);
    L.rho = c.num("liquid.density", L.rho);
    L.T = c.num("liquid.temperature", L.T);
    L.status = liquid_status_from_string(c.str("liquid.status", to_string(L.status)));
    L.source = c.str("liquid.source", "");
    L.caveat = c.str("liquid.caveat", L.caveat);
    if (L.status != LiquidDataStatus::Illustrative && L.source.empty())
      throw std::runtime_error(
          "liquid.status ist '" + std::string(to_string(L.status)) +
          "', aber liquid.source fehlt.  Ein Stoffwert ohne Fundstelle ist keine Eingabe, "
          "sondern eine Behauptung.");
  }
  L.documented_only.mu = c.num("liquid.viscosity", L.documented_only.mu);
  L.documented_only.K = c.num("liquid.conductivity", L.documented_only.K);
  L.documented_only.eps_r = c.num("liquid.relative_permittivity", L.documented_only.eps_r);
  return L;
}

// ---------------------------------------------------------------------------

struct Case {
  std::string tag;
  Real Pi{0};
  CapillaryMeniscus m;
  SphericalCap cap;
  ResidualProfile res;
  Real err_normal{0}, err_z{0};
  Real kappa_min{0}, kappa_max{0}, kappa_mean{0};
};

/// Mean curvature measured on the polyline, recovered from the residual so that
/// the same independent evaluation backs both numbers.
void measure_curvature(Case& k) {
  const Real a = k.m.contact_radius, g = k.m.gamma;
  bool first = true;
  Real sum = 0.0;
  for (Real r : k.res.residual) {
    const Real kappa = r / a + k.m.delta_p_exit / g;
    if (first) { k.kappa_min = k.kappa_max = kappa; first = false; }
    k.kappa_min = std::min(k.kappa_min, kappa);
    k.kappa_max = std::max(k.kappa_max, kappa);
    sum += kappa;
  }
  k.kappa_mean = k.res.residual.empty() ? 0.0 : sum / static_cast<Real>(k.res.residual.size());
}

Case solve_case(const DeviceGeometry& g, const LiquidProperties& L, Real Pi, Real accuracy,
                const std::string& tag) {
  Case k;
  k.tag = tag;
  k.Pi = Pi;
  CapillaryRequest q;
  q.delta_p_exit = capillary::pressure_from_pi(Pi, g.contact_radius(), L.gamma);
  q.target_relative_accuracy = accuracy;
  k.m = solve_capillary_meniscus(g, L, q);
  if (is_usable(k.m.status)) {
    k.cap = spherical_cap(k.m.contact_radius, k.m.delta_p_exit, L.gamma);
    k.res = young_laplace_residual(k.m);
    k.err_normal = profile_error_against_cap(k.m);
    k.err_z = profile_z_error_against_cap(k.m);
    measure_curvature(k);
  }
  return k;
}

std::string fmt(Real x) {
  char b[32];
  std::snprintf(b, sizeof b, "%.9e", x);
  return b;
}

}  // namespace

// ===========================================================================

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.size() < 2) {
    std::printf(
        "es_capillary -- P3a: statischer Kapillarmeniskus ohne elektrisches Feld\n\n"
        "  es_capillary <geometrie.cfg> [<p3a.cfg> ...] <ausgabeverzeichnis> [key=value ...]\n\n"
        "Beispiel:\n"
        "  ./build/es_capillary examples/device_p1.cfg examples/capillary_p3a.cfg \\\n"
        "      results/<ordner> meta.commit=$(git rev-parse HEAD)\n"
        "  python python/plot_capillary.py results/<ordner>\n\n"
        "Gerechnet wird ausschliesslich das statische Kapillargleichgewicht:\n"
        "  gamma (dpsi/ds + sin psi / r) = delta_p_exit = p_fluessig - p_vakuum.\n"
        "Kein elektrisches Feld, keine Emission, keine Stroemung, keine Schwerkraft,\n"
        "keine Stabilitaetsaussage, kein Taylor-Kegel.\n");
    return 1;
  }

  Config cfg;
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) cfg.load(pos[i]);
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  // The figures have to name the configuration they were produced from, so the
  // files themselves travel with the results.
  std::string config_names;
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) {
    const std::filesystem::path src(pos[i]);
    std::filesystem::copy_file(src, std::filesystem::path(outdir) / src.filename(),
                               std::filesystem::copy_options::overwrite_existing);
    config_names += (config_names.empty() ? "" : ";") + src.filename().string();
  }

  refuse_contact_angle(cfg);
  const DeviceGeometry geom = DeviceGeometry::build(device_from(cfg));
  const LiquidProperties liquid = liquid_from(cfg);
  liquid.validate_or_throw();

  const Real a = geom.contact_radius();
  const Real p_scale = liquid.capillary_pressure_scale(a);   // gamma / a
  const Real dp_max = capillary::kMaxAbsPi * p_scale;        // hemispherical limit
  const Real Bo = liquid.bond_number(a);
  const Real accuracy = cfg.num("capillary.accuracy", 1.0e-11);
  if (!(accuracy > 0.0) || accuracy > 1.0e-4)
    throw std::runtime_error("capillary.accuracy muss in (0, 1e-4] liegen.");

  std::printf("P3a -- statischer Kapillarmeniskus, KEIN elektrisches Feld\n");
  std::printf("  Pinningradius a = phi_2/2 = %.6g m (aus DeviceParameters)\n", a);
  std::printf("  Kapillardruckskala gamma/a = %.6g Pa, Grenze |delta_p| <= %.6g Pa\n",
              p_scale, dp_max);
  std::printf("  Bond-Zahl rho g a^2 / gamma = %.4e\n", Bo);
  liquid.print(stdout);

  geom.write_csv(outdir);

  int exit_code = 0;

  // ------------------------------------------------------------------ profiles
  const std::vector<Real> pi_list{-1.98, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 1.98};
  std::vector<Case> cases;
  for (Real Pi : pi_list) {
    char tag[32];
    std::snprintf(tag, sizeof tag, "pi%+05.2f", Pi);
    cases.push_back(solve_case(geom, liquid, Pi, accuracy, tag));
    const Case& k = cases.back();
    std::printf("  %-9s Pi=%+6.3f  Status %-30s", k.tag.c_str(), k.Pi, to_string(k.m.status));
    if (is_usable(k.m.status))
      std::printf(" h/a=%+10.6f  n=%5d  Profilfehler %.2e\n", k.m.apex_height / a,
                  k.m.n_intervals, k.err_normal);
    else
      std::printf("\n");
    if (!is_usable(k.m.status)) exit_code = 2;
  }

  {
    std::FILE* f = std::fopen((outdir + "/profiles.csv").c_str(), "w");
    std::fprintf(f, "# P3a: numerische Meniskusprofile, Apex -> gepinnte Kontaktlinie.\n"
                    "# z_analytic ist die geschlossene Kugelkappe am selben Radius.\n"
                    "# dz_m und dn_m sind die Abweichungen davon, HIER in voller doppelter\n"
                    "# Genauigkeit gebildet: als Spaltendifferenz waeren sie von der\n"
                    "# Ausgabegenauigkeit der Koordinaten verdeckt.  dn_m ist der Abstand\n"
                    "# senkrecht zur Kugelflaeche, dz_m die senkrechte Differenz bei gleichem r.\n");
    std::fprintf(f, "variant,Pi,delta_p_Pa,node,r_m,z_m,psi_rad,z_analytic_m,dz_m,dn_m\n");
    for (const Case& k : cases) {
      if (!is_usable(k.m.status)) continue;
      const Real absR = std::abs(k.cap.sphere_radius);
      const Real sg = (k.cap.sphere_radius > 0.0) ? 1.0 : -1.0;
      const Real z_centre =
          k.m.contact_z - sg * std::sqrt(std::max(0.0, absR * absR - a * a));
      for (std::size_t i = 0; i < k.m.nodes.size(); ++i) {
        const Real r = std::min(k.m.nodes[i].r, a);
        const Real za = k.cap.z_at_radius(r);
        const Real dn = (k.m.delta_p_exit == 0.0)
                            ? k.m.nodes[i].z - k.m.contact_z
                            : std::hypot(k.m.nodes[i].r, k.m.nodes[i].z - z_centre) - absR;
        std::fprintf(f, "%s,%.9e,%.9e,%zu,%.9e,%.9e,%.9e,%.9e,%.17e,%.17e\n", k.tag.c_str(),
                     k.Pi, k.m.delta_p_exit, i, k.m.nodes[i].r, k.m.nodes[i].z, k.m.psi[i],
                     za, (k.m.nodes[i].z - k.m.contact_z) - za, dn);
      }
    }
    std::fclose(f);
  }

  {
    std::FILE* f = std::fopen((outdir + "/profile_summary.csv").c_str(), "w");
    std::fprintf(f, "# Numerisch gegen geschlossene Kugelkappe. Fehler sind auf a bzw. a^2, a^3 "
                    "bezogen.\n");
    std::fprintf(f, "variant,Pi,delta_p_Pa,status,n_intervals,est_rel_error,"
                    "h_num_m,h_ana_m,s_num_m,s_ana_m,A_num_m2,A_ana_m2,V_num_m3,V_ana_m3,"
                    "A_polyline_m2,V_polyline_m3,psi_contact_rad,psi_contact_ana_rad,"
                    "err_profile_normal,err_profile_z,err_h_rel,err_A_rel,err_V_rel,"
                    "residual_max,residual_rms,kappa_min,kappa_max,kappa_ana,"
                    "contact_radius_error_rel\n");
    for (const Case& k : cases) {
      if (!is_usable(k.m.status)) {
        std::fprintf(f, "%s,%.9e,%.9e,%s,,,,,,,,,,,,,,,,,,,,,,,,,\n", k.tag.c_str(), k.Pi,
                     capillary::pressure_from_pi(k.Pi, a, liquid.gamma),
                     to_string(k.m.status));
        continue;
      }
      const Real eh = std::abs(k.m.apex_height - k.cap.apex_height) /
                      std::max(std::abs(k.cap.apex_height), 1e-300);
      const Real eA = std::abs(k.m.revolved_area - k.cap.revolved_area) / k.cap.revolved_area;
      const Real eV = k.cap.revolved_volume != 0.0
                          ? std::abs(k.m.revolved_volume - k.cap.revolved_volume) /
                                std::abs(k.cap.revolved_volume)
                          : 0.0;
      std::fprintf(f,
                   "%s,%.9e,%.9e,%s,%d,%.9e,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"
                   "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                   k.tag.c_str(), k.Pi, k.m.delta_p_exit, to_string(k.m.status),
                   k.m.n_intervals, k.m.estimated_relative_error,
                   fmt(k.m.apex_height).c_str(), fmt(k.cap.apex_height).c_str(),
                   fmt(k.m.arclength).c_str(), fmt(k.cap.arclength).c_str(),
                   fmt(k.m.revolved_area).c_str(), fmt(k.cap.revolved_area).c_str(),
                   fmt(k.m.revolved_volume).c_str(), fmt(k.cap.revolved_volume).c_str(),
                   fmt(k.m.polyline_area).c_str(), fmt(k.m.polyline_volume).c_str(),
                   fmt(k.m.contact_tangent_angle).c_str(),
                   fmt(k.cap.contact_tangent_angle).c_str(),
                   k.err_normal, k.err_z, eh, eA, eV, k.res.max_abs, k.res.rms,
                   k.kappa_min, k.kappa_max, k.cap.curvature,
                   std::abs(k.m.contact().r - a) / a);
    }
    std::fclose(f);
  }

  {
    std::FILE* f = std::fopen((outdir + "/residuals.csv").c_str(), "w");
    std::fprintf(f, "# Young-Laplace-Residuum (gamma*kappa - delta_p) * a / gamma, ausgewertet\n"
                    "# ALLEIN aus den Knotenkoordinaten: Tangentenwinkel aus den Sehnen,\n"
                    "# Meridiankruemmung aus deren Drehung.  Teilt keinen Code mit dem Loeser.\n");
    std::fprintf(f, "variant,Pi,s_m,s_over_arclength,residual\n");
    for (const Case& k : cases) {
      if (!is_usable(k.m.status)) continue;
      for (std::size_t i = 0; i < k.res.s.size(); ++i)
        std::fprintf(f, "%s,%.9e,%.9e,%.9e,%.9e\n", k.tag.c_str(), k.Pi, k.res.s[i],
                     k.res.s[i] / k.m.arclength, k.res.residual[i]);
    }
    std::fclose(f);
  }

  // --------------------------------------------------------------- convergence
  const std::vector<int> levels{16, 32, 64, 128, 256};
  const std::vector<Real> conv_pi{-1.5, -0.5, 0.5, 1.5};
  {
    std::FILE* f = std::fopen((outdir + "/convergence.csv").c_str(), "w");
    std::fprintf(f, "# Netzkonvergenz bei vorgegebener Intervallzahl.  Die Intervallzahl ist\n"
                    "# KEINE Benutzereingabe: im Normalbetrieb waehlt der Loeser sie aus der\n"
                    "# geforderten Genauigkeit.  Hier wird sie fuer die Studie erzwungen.\n");
    std::fprintf(f, "Pi,n_intervals,err_profile_normal,err_profile_z,err_h_rel,err_A_rel,"
                    "err_V_rel,residual_max,residual_rms\n");
    for (Real Pi : conv_pi) {
      for (int n : levels) {
        CapillaryRequest q;
        q.delta_p_exit = capillary::pressure_from_pi(Pi, a, liquid.gamma);
        q.forced_intervals = n;
        const CapillaryMeniscus m = solve_capillary_meniscus(geom, liquid, q);
        if (!is_usable(m.status)) { exit_code = 2; continue; }
        const SphericalCap c = spherical_cap(a, m.delta_p_exit, liquid.gamma);
        const ResidualProfile R = young_laplace_residual(m);
        std::fprintf(f, "%.9e,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", Pi, n,
                     profile_error_against_cap(m), profile_z_error_against_cap(m),
                     std::abs(m.apex_height - c.apex_height) / std::abs(c.apex_height),
                     std::abs(m.revolved_area - c.revolved_area) / c.revolved_area,
                     std::abs(m.revolved_volume - c.revolved_volume) /
                         std::abs(c.revolved_volume),
                     R.max_abs, R.rms);
      }
    }
    std::fclose(f);
  }

  // -------------------------------------------------------------- pressure sweep
  Real last_solved_pi = 0.0;
  int n_refused = 0;
  {
    std::FILE* f = std::fopen((outdir + "/sweep.csv").c_str(), "w");
    std::fprintf(f, "# Apexhoehe und Kruemmung ueber dem dimensionslosen Druck Pi.\n"
                    "# Punkte ausserhalb |Pi| <= 2 werden NICHT gerechnet, sondern mit ihrem\n"
                    "# Status ausgewiesen -- keine Ersatzform, kein letzter Iterationsstand.\n");
    std::fprintf(f, "Pi,delta_p_Pa,status,h_over_a_num,h_over_a_ana,kappa_num_1_per_m,"
                    "kappa_ana_1_per_m,area_over_a2_num,volume_over_a3_num,psi_contact_deg\n");
    for (int i = -210; i <= 210; ++i) {
      const Real Pi = 0.01 * static_cast<Real>(i);
      CapillaryRequest q;
      q.delta_p_exit = capillary::pressure_from_pi(Pi, a, liquid.gamma);
      q.target_relative_accuracy = 1.0e-9;
      const CapillaryMeniscus m = solve_capillary_meniscus(geom, liquid, q);
      if (!is_usable(m.status)) {
        ++n_refused;
        std::fprintf(f, "%.9e,%.9e,%s,,,,%.9e,,,\n", Pi, q.delta_p_exit, to_string(m.status),
                     q.delta_p_exit / liquid.gamma);
        continue;
      }
      last_solved_pi = std::max(last_solved_pi, std::abs(Pi));
      const SphericalCap c = spherical_cap(a, m.delta_p_exit, liquid.gamma);
      Case k;
      k.m = m;
      k.res = young_laplace_residual(m);
      measure_curvature(k);
      std::fprintf(f, "%.9e,%.9e,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", Pi, m.delta_p_exit,
                   to_string(m.status), m.apex_height / a, c.apex_height / a, k.kappa_mean,
                   c.curvature, m.revolved_area / (a * a), m.revolved_volume / (a * a * a),
                   m.contact_tangent_angle * 180.0 / pi);
    }
    std::fclose(f);
  }

  // ------------------------------------------------- the material example (Kunze)
  const LiquidProperties kunze = emibf4_illustrative();
  const std::vector<Real> dp_list{-17000, -9000, -4500, -1000, 0, 1000, 4500, 9000, 17000};
  {
    std::FILE* f = std::fopen((outdir + "/liquid_example.csv").c_str(), "w");
    std::fprintf(f, "# STOFFBEISPIEL, kein belegter Betriebspunkt.  Die Stoffidentitaet\n"
                    "# EMI-BF4 ist in der Kunze-Dissertation belegt, die Zahlenwerte NICHT.\n"
                    "# Status: %s.\n", to_string(kunze.status));
    std::fprintf(f, "delta_p_Pa,Pi,status,h_m,h_over_a,kappa_1_per_m,area_m2,volume_m3,"
                    "psi_contact_deg\n");
    for (Real dp : dp_list) {
      CapillaryRequest q;
      q.delta_p_exit = dp;
      q.target_relative_accuracy = accuracy;
      const CapillaryMeniscus m = solve_capillary_meniscus(geom, kunze, q);
      if (!is_usable(m.status)) {
        std::fprintf(f, "%.9e,%.9e,%s,,,,,,\n", dp, dp * a / kunze.gamma, to_string(m.status));
        continue;
      }
      std::fprintf(f, "%.9e,%.9e,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", dp, m.Pi,
                   to_string(m.status), m.apex_height, m.apex_height / a, m.delta_p_exit / kunze.gamma,
                   m.revolved_area, m.revolved_volume, m.contact_tangent_angle * 180.0 / pi);
    }
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/liquid_example_profiles.csv").c_str(), "w");
    std::fprintf(f, "# Profile zum STOFFBEISPIEL EMI-BF4 (Status %s).\n", to_string(kunze.status));
    std::fprintf(f, "delta_p_Pa,Pi,node,r_m,z_m\n");
    for (Real dp : dp_list) {
      CapillaryRequest q;
      q.delta_p_exit = dp;
      q.target_relative_accuracy = accuracy;
      const CapillaryMeniscus m = solve_capillary_meniscus(geom, kunze, q);
      if (!is_usable(m.status)) continue;
      for (std::size_t i = 0; i < m.nodes.size(); ++i)
        std::fprintf(f, "%.9e,%.9e,%zu,%.9e,%.9e\n", dp, m.Pi, i, m.nodes[i].r, m.nodes[i].z);
    }
    std::fclose(f);
  }

  // ------------------------------------------------------- the refusal, on record
  {
    std::FILE* f = std::fopen((outdir + "/range.csv").c_str(), "w");
    std::fprintf(f, "# Verhalten an und jenseits der Grenze des darstellbaren Bereichs.\n"
                    "# Erwartet: |Pi| < 2 loesbar, |Pi| = 2 HemisphericalLimit, |Pi| > 2\n"
                    "# PressureOutsideCapillaryRange -- ohne Ersatzform.\n");
    std::fprintf(f, "Pi,delta_p_Pa,status,shape_returned,message\n");
    for (Real Pi : {1.9, 1.999, 1.999999999, 2.0, 2.000000001, 2.001, 2.05, 3.0, -2.05, -3.0}) {
      CapillaryRequest q;
      q.delta_p_exit = capillary::pressure_from_pi(Pi, a, liquid.gamma);
      q.target_relative_accuracy = 1.0e-9;
      const CapillaryMeniscus m = solve_capillary_meniscus(geom, liquid, q);
      std::string msg = m.message;
      for (char& ch : msg)
        if (ch == ',' || ch == '\n') ch = ' ';
      std::fprintf(f, "%.12e,%.9e,%s,%s,%s\n", Pi, q.delta_p_exit, to_string(m.status),
                   m.nodes.empty() ? "no" : "yes", msg.c_str());
    }
    std::fclose(f);
  }

  // -------------------------------------------------------------------- liquid
  {
    std::FILE* f = std::fopen((outdir + "/liquid.csv").c_str(), "w");
    std::fprintf(f, "# Stoffdatensatz mit Status und Fundstelle.  Ein Wert ohne Herkunft ist\n"
                    "# keine Eingabe, sondern eine Behauptung.\n");
    std::fprintf(f, "key,value\n");
    std::string src = liquid.source, cav = liquid.caveat;
    for (std::string* s : {&src, &cav})
      for (char& ch : *s)
        if (ch == ',' || ch == '\n') ch = ' ';
    std::fprintf(f, "substance,%s\n", liquid.substance.c_str());
    std::fprintf(f, "status,%s\n", to_string(liquid.status));
    std::fprintf(f, "status_meaning,%s\n", explain(liquid.status));
    std::fprintf(f, "temperature_K,%.9e\n", liquid.T);
    std::fprintf(f, "surface_tension_N_per_m,%.9e\n", liquid.gamma);
    std::fprintf(f, "density_kg_per_m3,%.9e\n", liquid.rho);
    std::fprintf(f, "source,%s\n", src.c_str());
    std::fprintf(f, "caveat,%s\n", cav.c_str());
    std::fprintf(f, "not_used_viscosity_Pa_s,%.9e\n", liquid.documented_only.mu);
    std::fprintf(f, "not_used_conductivity_S_per_m,%.9e\n", liquid.documented_only.K);
    std::fprintf(f, "not_used_relative_permittivity,%.9e\n", liquid.documented_only.eps_r);
    std::fclose(f);
  }

  // ---------------------------------------------------------------- parameters
  {
    std::FILE* f = std::fopen((outdir + "/capillary_parameters.csv").c_str(), "w");
    std::fprintf(f, "# Groessen des Kapillarproblems, SI.\n");
    std::fprintf(f, "name,value_SI,unit,role\n");
    std::fprintf(f, "contact_radius,%.9e,m,gepinnter Kontaktradius a = phi_2/2\n", a);
    std::fprintf(f, "contact_z,%.9e,m,Lage der gepinnten Kontaktlinie\n",
                 geom.feature(FeatureId::PinnedContactEdge).z);
    std::fprintf(f, "surface_tension,%.9e,N/m,gamma\n", liquid.gamma);
    std::fprintf(f, "capillary_pressure_scale,%.9e,Pa,gamma/a\n", p_scale);
    std::fprintf(f, "delta_p_max,%.9e,Pa,hemisphaerische Grenze 2 gamma/a\n", dp_max);
    std::fprintf(f, "bond_number,%.9e,-,rho g a^2 / gamma (Schwerkraft NICHT gekoppelt)\n", Bo);
    std::fprintf(f, "hydrostatic_over_capillary,%.9e,-,rho g a / (gamma/a) = Bond-Zahl\n", Bo);
    std::fprintf(f, "requested_accuracy,%.9e,-,geforderte relative Profilgenauigkeit\n",
                 accuracy);
    std::fclose(f);
  }

  // -------------------------------------------------------------------- report
  {
    std::FILE* f = std::fopen((outdir + "/report.txt").c_str(), "w");
    std::fprintf(f, "P3a -- statischer Kapillarmeniskus ohne elektrisches Feld\n");
    std::fprintf(f, "=========================================================\n\n");

    std::fprintf(f, "WAS GERECHNET WURDE\n");
    std::fprintf(f, "  gamma (dpsi/ds + sin psi / r) = delta_p_exit, achsensymmetrisch, in\n");
    std::fprintf(f, "  Bogenlaengenparametrisierung (r(s), z(s), psi(s)) vom Apex zur Kante.\n");
    std::fprintf(f, "  delta_p_exit = p_fluessig - p_vakuum an der Austrittsebene z = 0.\n");
    std::fprintf(f, "  Positives delta_p woelbt die Oberflaeche nach +z (zum Extraktor),\n");
    std::fprintf(f, "  negatives zieht sie in die Bohrung; delta_p = 0 ergibt exakt eben.\n");
    std::fprintf(f, "  Der Grenzwert sin(psi)/r -> dpsi/ds auf der Achse ist analytisch\n");
    std::fprintf(f, "  eingesetzt (dpsi/ds|Apex = kappa/2); es wird nirgends durch r = 0\n");
    std::fprintf(f, "  geteilt.\n\n");

    std::fprintf(f, "GEOMETRIE UND RANDBEDINGUNGEN\n");
    std::fprintf(f, "  Pinningradius a = phi_2/2 = %.6g m, direkt aus DeviceParameters.\n", a);
    std::fprintf(f, "  Kontaktlinie gepinnt bei (r, z) = (%.6g, %.6g) m.\n", a,
                 geom.feature(FeatureId::PinnedContactEdge).z);
    std::fprintf(f, "  Kein Kontaktwinkel: an einer gepinnten Kante waere er eine zweite,\n");
    std::fprintf(f, "  widerspruechliche Festlegung.  Beides gleichzeitig wird abgelehnt.\n");
    std::fprintf(f, "  Symmetrie bei r = 0: psi(0) = 0.\n\n");

    std::fprintf(f, "STOFFDATEN\n");
    std::fprintf(f, "  Stoff  : %s\n", liquid.substance.c_str());
    std::fprintf(f, "  Status : %s -- %s\n", to_string(liquid.status), explain(liquid.status));
    std::fprintf(f, "  T = %.4g K, gamma = %.6g N/m, rho = %.6g kg/m^3\n", liquid.T,
                 liquid.gamma, liquid.rho);
    std::fprintf(f, "  Fundstelle: %s\n", liquid.source.c_str());
    if (!liquid.caveat.empty()) std::fprintf(f, "  Vorbehalt : %s\n", liquid.caveat.c_str());
    std::fprintf(f, "  mu, K und eps_r sind vorgemerkt und gehen in P3a NICHT ein.\n\n");

    std::fprintf(f, "SCHWERKRAFT\n");
    std::fprintf(f, "  Bond-Zahl Bo = rho g a^2 / gamma = %.4e fuer a = %.4g m.\n", Bo, a);
    std::fprintf(f, "  Der hydrostatische Druck ueber eine Bohrungsradiushoehe betraegt damit\n");
    std::fprintf(f, "  %.4e der Kapillardruckskala gamma/a = %.5g Pa, also %.3g Pa.\n", Bo,
                 p_scale, Bo * p_scale);
    {
      Real h_max = 0.0;
      for (const Case& k : cases)
        if (is_usable(k.m.status)) h_max = std::max(h_max, std::abs(k.m.apex_height));
      std::fprintf(f, "  Ueber die groesste hier gerechnete Apexhoehe |h| = %.4g m ist der\n",
                   h_max);
      std::fprintf(f, "  hydrostatische Druck rho g |h| = %.3g Pa, also %.2e der Kapillar-\n",
                   liquid.rho * constants::g0 * h_max,
                   liquid.rho * constants::g0 * h_max / p_scale);
      std::fprintf(f, "  druckskala.  Die Vernachlaessigung ist damit fuer diese Geometrie\n");
      std::fprintf(f, "  quantitativ gerechtfertigt -- fuer eine millimetergrosse Geometrie\n");
      std::fprintf(f, "  waere sie es NICHT (Bo skaliert mit a^2).\n");
    }
    std::fprintf(f, "  Nicht behauptet wird damit der kapillare Aufstieg vom Vorrat bis zur\n");
    std::fprintf(f, "  Kante: der gefuellte Zulauf ist eine VORAUSSETZUNG dieses Laufs.\n\n");

    std::fprintf(f, "GUELTIGER DRUCKBEREICH\n");
    std::fprintf(f, "  |Pi| = |delta_p| a / gamma <= 2, also |delta_p| <= 2 gamma/a = %.6g Pa.\n",
                 dp_max);
    std::fprintf(f, "  Pi = 2 ist die Halbkugel (psi = 90 Grad an der Kante).  Darueber gibt es\n");
    std::fprintf(f, "  keine glatte, an a gepinnte Flaeche: der Kugelradius 2 gamma/delta_p\n");
    std::fprintf(f, "  waere kleiner als a.  Der Loeser stellt das an der integrierten Form\n");
    std::fprintf(f, "  fest -- die Meridiankurve wird senkrecht, bevor sie a erreicht -- und\n");
    std::fprintf(f, "  meldet PressureOutsideCapillaryRange.  Im Sweep wurden %d von %d\n",
                 n_refused, 421);
    std::fprintf(f, "  Punkten so abgelehnt; der groesste geloeste |Pi| war %.3f.\n\n",
                 last_solved_pi);

    std::fprintf(f, "NUMERIK GEGEN GESCHLOSSENE LOESUNG\n");
    std::fprintf(f, "  %-9s %8s %8s %12s %12s %12s %12s %10s\n", "Kennung", "Pi", "h/a",
                 "Profil/a", "|dz|/a", "dA/A", "dV/V", "Residuum");
    for (const Case& k : cases) {
      if (!is_usable(k.m.status)) {
        std::fprintf(f, "  %-9s %8.3f %8s %12s %12s %12s %12s %10s  %s\n", k.tag.c_str(), k.Pi,
                     "-", "-", "-", "-", "-", "-", to_string(k.m.status));
        continue;
      }
      const Real eA = std::abs(k.m.revolved_area - k.cap.revolved_area) / k.cap.revolved_area;
      const Real eV = k.cap.revolved_volume != 0.0
                          ? std::abs(k.m.revolved_volume - k.cap.revolved_volume) /
                                std::abs(k.cap.revolved_volume)
                          : 0.0;
      std::fprintf(f, "  %-9s %8.3f %8.5f %12.3e %12.3e %12.3e %12.3e %10.2e\n", k.tag.c_str(),
                   k.Pi, k.m.apex_height / a, k.err_normal, k.err_z, eA, eV, k.res.max_abs);
    }
    std::fprintf(f, "  Das Residuum wird ALLEIN aus den Knotenkoordinaten gebildet und teilt\n");
    std::fprintf(f, "  keinen Code mit dem Integrator.\n\n");

    std::fprintf(f, "KONSTANTE MITTLERE KRUEMMUNG\n");
    for (const Case& k : cases) {
      if (!is_usable(k.m.status)) continue;
      std::fprintf(f, "  Pi = %+6.3f: kappa gemessen in [%.6e, %.6e] 1/m, analytisch %.6e 1/m, "
                      "Spannweite/kappa = %.2e\n",
                   k.Pi, k.kappa_min, k.kappa_max, k.cap.curvature,
                   k.cap.curvature != 0.0
                       ? std::abs(k.kappa_max - k.kappa_min) / std::abs(k.cap.curvature)
                       : 0.0);
    }
    std::fprintf(f, "\n");

    std::fprintf(f, "DISKRETISIERUNG\n");
    std::fprintf(f, "  Es gibt keine Netzeingabe.  Eingabe ist die geforderte relative\n");
    std::fprintf(f, "  Profilgenauigkeit (%.1e); die Intervallzahl wird durch Verdopplung\n",
                 accuracy);
    std::fprintf(f, "  bestimmt, bis die Aenderung zur halben Aufloesung darunter liegt.\n");
    for (const Case& k : cases)
      if (is_usable(k.m.status))
        std::fprintf(f, "    Pi = %+6.3f -> %5d Intervalle, Schaetzfehler %.2e\n", k.Pi,
                     k.m.n_intervals, k.m.estimated_relative_error);
    std::fprintf(f, "\n");

    std::fprintf(f, "WAS AUSDRUECKLICH FEHLT\n");
    std::fprintf(f, "  Kein elektrisches Feld, kein Maxwell-Druck, keine Kopplung an den\n");
    std::fprintf(f, "  FEM-Elektrostatiksolver, keine Betriebsspannung.  Keine Emission,\n");
    std::fprintf(f, "  keine Raumladung.  Keine Stroemung, kein viskoser Druckabfall, keine\n");
    std::fprintf(f, "  Speiseimpedanz.  Keine Zeitabhaengigkeit und keine Stabilitaetsaussage:\n");
    std::fprintf(f, "  dass eine Loesung hier existiert, sagt nichts darueber, ob sie stabil\n");
    std::fprintf(f, "  ist.  Kein Taylor-Kegel, kein Cone-Jet -- die Formen hier sind\n");
    std::fprintf(f, "  Kugelkappen.  Keine Schwerkraft (siehe Bond-Zahl oben).  Die\n");
    std::fprintf(f, "  P2c-Vorratsgeometrie ist unveraendert und wird nicht volumetrisch\n");
    std::fprintf(f, "  geloest.\n");
    std::fclose(f);
  }

  // ---------------------------------------------------------------------- meta
  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_capillary (P3a)\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "config=%s\n", config_names.c_str());
    std::fprintf(f, "state=static capillary equilibrium, no electric field, no emission, "
                    "no flow, no gravity, no stability\n");
    std::fprintf(f, "contact_radius_m=%.9e\n", a);
    std::fprintf(f, "contact_z_m=%.9e\n", geom.feature(FeatureId::PinnedContactEdge).z);
    std::fprintf(f, "phi_2_m=%.9e\n", geom.parameters().phi_2);
    std::fprintf(f, "phi_1_m=%.9e\n", geom.parameters().phi_1);
    std::fprintf(f, "phi_3_m=%.9e\n", geom.parameters().phi_3);
    std::fprintf(f, "emitter_height_m=%.9e\n", geom.parameters().emitter_height);
    std::fprintf(f, "surface_tension_N_per_m=%.9e\n", liquid.gamma);
    std::fprintf(f, "density_kg_per_m3=%.9e\n", liquid.rho);
    std::fprintf(f, "temperature_K=%.9e\n", liquid.T);
    std::fprintf(f, "liquid_substance=%s\n", liquid.substance.c_str());
    std::fprintf(f, "liquid_status=%s\n", to_string(liquid.status));
    std::fprintf(f, "capillary_pressure_scale_Pa=%.9e\n", p_scale);
    std::fprintf(f, "delta_p_max_Pa=%.9e\n", dp_max);
    std::fprintf(f, "pi_max=%.9e\n", capillary::kMaxAbsPi);
    std::fprintf(f, "bond_number=%.9e\n", Bo);
    std::fprintf(f, "requested_accuracy=%.9e\n", accuracy);
    std::fprintf(f, "example_liquid_substance=%s\n", kunze.substance.c_str());
    std::fprintf(f, "example_liquid_status=%s\n", to_string(kunze.status));
    std::fprintf(f, "example_liquid_gamma=%.9e\n", kunze.gamma);
    std::fprintf(f, "example_liquid_rho=%.9e\n", kunze.rho);
    std::fprintf(f, "example_liquid_bond_number=%.9e\n", kunze.bond_number(a));
    std::fprintf(f, "example_liquid_delta_p_max_Pa=%.9e\n",
                 capillary::kMaxAbsPi * kunze.capillary_pressure_scale(a));
    std::fprintf(f, "sweep_points_refused=%d\n", n_refused);
    std::fprintf(f, "largest_solved_abs_pi=%.9e\n", last_solved_pi);
    std::fclose(f);
  }

  cfg.warn_about_unused(stdout, {"meta.", "fluid.", "beam.", "output.", "bem.", "wetting.",
                                 "species.", "feed.", "field.", "mesh.", "emitter.",
                                 "extractor.", "reservoir.", "material."});
  std::printf("\ngeschrieben nach %s\n", outdir.c_str());
  return exit_code;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\n%s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "\nFehler: %s\n", e.what());
  return 2;
}
