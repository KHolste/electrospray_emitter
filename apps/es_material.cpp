// es_material -- P2: the material-data audit for the working liquid.
//
//   es_material <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS
//
// It writes out the material-data record: every source, every measured point,
// the measurement method, the sample purity and the water content the source
// stated -- and, just as importantly, the fields it did NOT state.  It computes
// no physics.
//
// WHAT IT SETTLES
//
// Surface tension enters the P3a/P3b equilibrium linearly, so an unsourced
// gamma makes every pressure, every apex height and (through V^2 ~ gamma) every
// voltage of those phases uncertain by the same factor.  This run puts the
// value that was used next to the values that are documented.
//
// Nothing is averaged.  Sources that disagree are all written out and the
// spread is reported as what it is: the honest uncertainty of the number, which
// is far larger than any single source's own error bar.
//
// TWO THINGS THIS RUN IS CAREFUL ABOUT
//
//   * nu = mu/rho is DERIVED, not measured, and is written to its own file with
//     both parents, both parent values, the conditions that were checked and
//     the propagated uncertainty.  It never appears in a column that could be
//     read as a measurement of nu.
//
//   * impact.csv scales an already computed dimensionless equilibrium.  It is
//     not a new coupled simulation and says so on every row.  In particular the
//     Maxwell traction at fixed geometry, voltage and permittivity does NOT
//     scale with gamma -- exponent 0 -- because gamma appears nowhere in the
//     field problem that determines it.
//
// Exit code 2 means a declared check failed.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/config.hpp"
#include "es/liquid.hpp"
#include "es/material_data.hpp"

using namespace es;

namespace {

void put(std::FILE* f, Real v) {
  if (std::isfinite(v))
    std::fprintf(f, ",%.9e", v);
  else
    std::fprintf(f, ",nan");
}

/// A quoted CSV field.  Replacing the commas inside a citation would destroy
/// the citation, so the field is quoted instead and the reader unquotes it.
std::string csv_quote(const std::string& in) {
  std::string out = "\"";
  for (char c : in) {
    if (c == '"')
      out += "\"\"";
    else if (c == '\n' || c == '\r')
      out += ' ';
    else
      out += c;
  }
  out += '"';
  return out;
}

}  // namespace

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.empty()) {
    std::printf("es_material -- P2: Stoffdatenaudit\n\n"
                "  es_material <ausgabeverzeichnis> [key=value ...]\n\n"
                "Schreibt den Stoffdatensatz mit Provenienz: jede Quelle, jeder Messpunkt,\n"
                "Methode, Reinheit, Wassergehalt, Druck und Frequenz.  Es wird nicht\n"
                "gemittelt und nichts extrapoliert.\n");
    return 1;
  }
  Config cfg;
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  const MaterialDataset& d = emibf4_sourced();
  const Real T = cfg.num("material.temperature", 298.15);
  int exit_code = 0;

  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& s) {
    std::printf("%s\n", s.c_str());
    std::fprintf(log, "%s\n", s.c_str());
  };

  say("P2 -- Stoffdatenaudit fuer EMI-BF4");
  say("");
  d.print(stdout);
  d.print(log);

  // --- every source ---------------------------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/sources.csv").c_str(), "w");
    std::fprintf(f,
                 "# Jede Quelle je Stoffgroesse, mit Provenienz.  Ein LEERES Feld heisst,\n"
                 "# dass die Quelle die Angabe NICHT gemacht hat -- es heisst nie null.\n"
                 "# selected=yes markiert die eine Quelle, die eine Rechnung benutzt; die\n"
                 "# Auswahlregel steht in include/es/material_data.hpp und wird von\n"
                 "# tools/fetch_material_data.py mechanisch angewandt.\n"
                 "# n_non_ambient zaehlt Punkte abseits des Umgebungsdrucks: sie bleiben im\n"
                 "# Datensatz, tragen aber keinen Umgebungsdruckwert.\n");
    // TWO value columns, because they answer different questions.
    //   value_at_T  -- strict: interpolated between measured points, NaN outside
    //                  the measured range.  This is what a computation may use.
    //   value_near_T -- the nearest AMBIENT point within the same +-2 K window
    //                  the literature band uses.  A single measurement taken at
    //                  296.6 K is not a value at 298.15 K, but it IS part of the
    //                  scatter, and a figure that drew only the first column
    //                  would show a narrower band than the one reported.
    std::fprintf(f, "property,unit,database_id,selected,status,value_at_T,value_near_T,"
                    "T_near_K,uncertainty_at_T,"
                    "T_min_K,T_max_K,n_points,n_non_ambient,frequency_resolved,"
                    "states_method,states_purity,states_water,method,purity,water_content,"
                    "sample_source,constraints,reference,paper_title\n");
    for (std::size_t i = 0; i < d.n_properties; ++i) {
      const MaterialProperty& p = d.properties[i];
      for (std::size_t k = 0; k < p.n_sources; ++k) {
        const PropertySource& s = p.sources[k];
        Real unc = 0.0;
        Real near = std::numeric_limits<Real>::quiet_NaN();
        Real T_near = std::numeric_limits<Real>::quiet_NaN();
        Real best = 2.0;
        for (std::size_t j = 0; j < s.n_points; ++j) {
          if (!s.points[j].ambient()) continue;
          const Real dT = std::abs(s.points[j].T - T);
          if (dT > 2.0) continue;
          if (s.points[j].has_uncertainty()) unc = s.points[j].uncertainty;
          if (dT <= best) {
            best = dT;
            near = s.points[j].value;
            T_near = s.points[j].T;
          }
        }
        std::fprintf(f, "%s,%s,%s,%s,%s", to_string(p.kind), si_unit(p.kind), s.database_id,
                     (p.has_selection() && static_cast<std::size_t>(p.selected) == k) ? "yes"
                                                                                     : "no",
                     to_string(s.status));
        put(f, s.value_at(T));
        put(f, near);
        put(f, T_near);
        put(f, unc > 0.0 ? unc : std::numeric_limits<Real>::quiet_NaN());
        put(f, s.T_min());
        put(f, s.T_max());
        std::fprintf(f, ",%zu,%zu,%s,%s,%s,%s", s.n_points, s.n_non_ambient(),
                     s.is_frequency_resolved() ? "yes" : "no", s.states_method() ? "yes" : "no",
                     s.states_purity() ? "yes" : "no", s.states_water_content() ? "yes" : "no");
        std::fprintf(f, ",%s,%s,%s,%s,%s,%s,%s\n", csv_quote(s.method).c_str(),
                     csv_quote(s.purity).c_str(), csv_quote(s.water_content).c_str(),
                     csv_quote(s.sample_source).c_str(), csv_quote(s.constraints).c_str(),
                     csv_quote(s.reference).c_str(), csv_quote(s.paper_title).c_str());
      }
    }
    std::fclose(f);
    say("  sources.csv geschrieben");
  }

  // --- every measured point -------------------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/points.csv").c_str(), "w");
    std::fprintf(f, "# Jeder einzelne Messpunkt.  uncertainty = 0 heisst 'nicht angegeben',\n"
                    "# nicht 'null'.  pressure_Pa = 0 heisst 'nicht angegeben', was in\n"
                    "# diesen Zusammenstellungen Umgebungsdruck bedeutet.\n");
    std::fprintf(f, "property,database_id,selected,T_K,value,uncertainty,pressure_Pa,"
                    "frequency_Hz,ambient\n");
    for (std::size_t i = 0; i < d.n_properties; ++i) {
      const MaterialProperty& p = d.properties[i];
      for (std::size_t k = 0; k < p.n_sources; ++k) {
        const PropertySource& s = p.sources[k];
        const bool sel = p.has_selection() && static_cast<std::size_t>(p.selected) == k;
        for (std::size_t j = 0; j < s.n_points; ++j) {
          const PropertyPoint& q = s.points[j];
          std::fprintf(f, "%s,%s,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%s\n", to_string(p.kind),
                       s.database_id, sel ? "yes" : "no", q.T, q.value, q.uncertainty,
                       q.pressure_Pa, q.frequency_Hz, q.ambient() ? "yes" : "no");
        }
      }
    }
    std::fclose(f);
    say("  points.csv geschrieben");
  }

  // --- the summary that decides what a solver may use -----------------------
  {
    std::FILE* f = std::fopen((outdir + "/summary.csv").c_str(), "w");
    std::fprintf(f,
                 "# Was eine Rechnung benutzen darf, und wie unsicher es ist.\n"
                 "# illustrative_value ist der bisher benutzte Wert ohne Quelle aus\n"
                 "# src/liquid.cpp; er bleibt illustrative.\n"
                 "#\n"
                 "# value_kind unterscheidet, WIE der Wert zustande kommt:\n"
                 "#   direct_measurement -- die Groesse selbst ist gemessen; relative_spread\n"
                 "#                         ist dann die LITERATURSTREUUNG und damit die\n"
                 "#                         ehrliche Unsicherheit, groesser als jede\n"
                 "#                         einzelne Fehlerangabe.\n"
                 "#   derived            -- die Groesse selbst ist NICHT gemessen, sondern\n"
                 "#                         ueber eine genannte Identitaet aus anderen\n"
                 "#                         ausgewaehlten Groessen dieses Datensatzes\n"
                 "#                         berechnet; relative_spread ist dann die\n"
                 "#                         FORTGEPFLANZTE Unsicherheit der Eltern und keine\n"
                 "#                         Literaturstreuung.  Vollstaendiger Nachweis mit\n"
                 "#                         beiden Eltern in kinematic_viscosity.csv.\n"
                 "#   none               -- MissingMaterialData; es gibt keinen Wert.\n");
    std::fprintf(f, "property,unit,T_K,status,value_kind,selected_value,band_lo,band_hi,"
                    "relative_spread,n_sources_total,n_sources_at_T,illustrative_value,"
                    "illustrative_vs_selected,selected_reference\n");
    const LiquidProperties ill = emibf4_illustrative();
    for (std::size_t i = 0; i < d.n_properties; ++i) {
      const MaterialProperty& p = d.properties[i];
      // nu is not measured here, but it IS available: the contract derives it
      // from mu and rho under conditions it checks.  Reporting it as absent
      // would be as wrong as reporting it as measured, so it is reported as
      // what it is, in its own column.
      const MaterialValue v = (p.kind == PropertyKind::KinematicViscosity)
                                  ? derived_kinematic_viscosity(d, T)
                                  : material_value(d, p.kind, T);
      const char* value_kind = !v.usable() ? "none"
                               : (v.status == MaterialDataStatus::Derived) ? "derived"
                                                                           : "direct_measurement";
      Real illv = std::numeric_limits<Real>::quiet_NaN();
      switch (p.kind) {
        case PropertyKind::SurfaceTension: illv = ill.gamma; break;
        case PropertyKind::Density: illv = ill.rho; break;
        case PropertyKind::DynamicViscosity: illv = ill.documented_only.mu; break;
        case PropertyKind::ElectricalConductivity: illv = ill.documented_only.K; break;
        case PropertyKind::RelativePermittivity: illv = ill.documented_only.eps_r; break;
        default: break;
      }
      std::fprintf(f, "%s,%s,%.9e,%s,%s", to_string(p.kind), si_unit(p.kind), T,
                   to_string(v.status), value_kind);
      put(f, v.usable() ? v.value : std::numeric_limits<Real>::quiet_NaN());
      put(f, p.min_at(T));
      put(f, p.max_at(T));
      put(f, (v.status == MaterialDataStatus::Derived) ? v.relative_spread
                                                       : p.relative_spread_at(T));
      std::fprintf(f, ",%zu,%zu", p.n_sources, p.n_sources_at(T));
      put(f, illv);
      put(f, (v.usable() && std::isfinite(illv) && std::abs(v.value) > 0.0)
                 ? (illv - v.value) / v.value
                 : std::numeric_limits<Real>::quiet_NaN());
      if (v.status == MaterialDataStatus::Derived) {
        const KinematicViscosityDerivation r = derive_kinematic_viscosity(d, T);
        std::fprintf(f, ",%s\n",
                     csv_quote(std::string("ABGELEITET nu = mu/rho aus: ") +
                               r.mu_source->reference + " (mu) und " + r.rho_source->reference +
                               " (rho)")
                         .c_str());
      } else {
        std::fprintf(f, ",%s\n", v.source ? csv_quote(v.source->reference).c_str() : "");
      }
    }
    std::fclose(f);
    say("  summary.csv geschrieben");
  }

  // --- nu = mu/rho: derived, and only under stated conditions ---------------
  //
  // nu is not an independent property.  The contract does not select a directly
  // measured value -- none of the three sources states purity and water content
  // -- but mu and rho of THIS data set do have selections, and here they are
  // even the same publication and the same sample.  The derivation is therefore
  // performed, labelled `derived`, and written out with both parents, both
  // parent values, the conditions and the propagated uncertainty.  It is never
  // written into a column that a reader could mistake for a measurement of nu.
  {
    const KinematicViscosityDerivation r = derive_kinematic_viscosity(d, T);
    r.print(stdout);
    r.print(log);
    std::FILE* f = std::fopen((outdir + "/kinematic_viscosity.csv").c_str(), "w");
    std::fprintf(f,
                 "# nu = mu/rho -- ABGELEITET, nicht gemessen.\n"
                 "#\n"
                 "# Der Vertrag waehlt keine direkt gemessene kinematische Viskositaet aus:\n"
                 "# keine der drei Quellen nennt Reinheit und Wassergehalt.  Aus mu und rho\n"
                 "# darf ein Wert aber abgeleitet werden, sobald gezeigt ist, dass beide\n"
                 "# Werte dieselbe Fluessigkeit im selben Zustand beschreiben.  Die\n"
                 "# Bedingungen C1..C4 stehen in include/es/material_data.hpp; ist eine\n"
                 "# verletzt, wird kein Wert geliefert, sondern die Bedingung genannt.\n"
                 "#\n"
                 "# status=derived.  Diese Zeile darf nicht als Messung der kinematischen\n"
                 "# Viskositaet zitiert werden.  direct_band_lo/hi sind zum Vergleich die\n"
                 "# direkt gemessenen Werte, die die Auswahlregel NICHT erfuellen; sie gehen\n"
                 "# in keine Rechnung ein.\n");
    std::fprintf(f, "status,T_K,identity,nu,uncertainty_quadratic,uncertainty_linear,"
                    "relative_uncertainty,uncertainty_propagated,mu,mu_uncertainty,rho,"
                    "rho_uncertainty,same_publication,conditions,mu_reference,rho_reference,"
                    "mu_method,rho_method,direct_band_lo,direct_band_hi,"
                    "deviation_from_direct_band_hi,blocker\n");
    const MaterialProperty* kv = d.find(PropertyKind::KinematicViscosity);
    const Real dlo = kv ? kv->min_at(T) : std::numeric_limits<Real>::quiet_NaN();
    const Real dhi = kv ? kv->max_at(T) : std::numeric_limits<Real>::quiet_NaN();
    if (r.ok) {
      std::fprintf(f, "derived,%.9e,%s", T, r.identity);
      put(f, r.value);
      put(f, r.uncertainty);
      put(f, r.uncertainty_linear);
      put(f, r.relative_uncertainty);
      std::fprintf(f, ",%s", r.uncertainty_propagated ? "yes" : "no");
      put(f, r.mu);
      put(f, r.mu_uncertainty);
      put(f, r.rho);
      put(f, r.rho_uncertainty);
      std::fprintf(f, ",%s,%s,%s,%s,%s,%s", r.same_publication ? "yes" : "no",
                   csv_quote(r.conditions).c_str(), csv_quote(r.mu_source->reference).c_str(),
                   csv_quote(r.rho_source->reference).c_str(),
                   csv_quote(r.mu_source->method).c_str(),
                   csv_quote(r.rho_source->method).c_str());
      put(f, dlo);
      put(f, dhi);
      put(f, (std::isfinite(dhi) && dhi > 0.0) ? (r.value - dhi) / dhi
                                               : std::numeric_limits<Real>::quiet_NaN());
      std::fprintf(f, ",\n");
      char line[512];
      std::snprintf(line, sizeof line,
                    "  nu = %.6g m^2/s +- %.3g (%.2f %%, fortgepflanzt), abgeleitet aus mu "
                    "und rho derselben Publikation.  Status derived, NICHT measured.",
                    r.value, r.uncertainty, 100.0 * r.relative_uncertainty);
      say("");
      say(line);
      if (std::isfinite(dhi)) {
        std::snprintf(line, sizeof line,
                      "  Zum Vergleich das direkt gemessene Band (erfuellt die Auswahlregel "
                      "nicht): %.6g .. %.6g m^2/s; der abgeleitete Wert liegt %+.1f %% ueber "
                      "dessen oberer Kante.",
                      dlo, dhi, 100.0 * (r.value - dhi) / dhi);
        say(line);
      }
    } else {
      std::fprintf(f, "MissingMaterialData,%.9e,%s,nan,nan,nan,nan,no,nan,nan,nan,nan,no,,,,,",
                   T, r.identity);
      put(f, dlo);
      put(f, dhi);
      std::fprintf(f, ",nan,%s\n", csv_quote(r.blocker).c_str());
      say("  nu ist NICHT ableitbar: " + r.blocker);
      exit_code = 2;
    }
    std::fclose(f);
    say("  kinematic_viscosity.csv geschrieben");
  }

  // --- what this changes for P3a/P3b ---------------------------------------
  {
    const LiquidProperties ill = emibf4_illustrative();
    const MaterialProperty* g = d.find(PropertyKind::SurfaceTension);
    const MaterialValue gv = material_value(d, PropertyKind::SurfaceTension, T);
    std::FILE* f = std::fopen((outdir + "/impact.csv").c_str(), "w");
    std::fprintf(f,
                 "# Was der belegte gamma-Wert an den P3a/P3b-Zahlen aendert.\n"
                 "#\n"
                 "# JEDE Zeile ist eine SKALIERUNG einer bereits gerechneten dimensionslosen\n"
                 "# Loesung, KEINE neue gekoppelte Rechnung.  Die Spalte `recomputed` steht\n"
                 "# deshalb auf jeder Zeile auf no.\n"
                 "#\n"
                 "# `category` unterscheidet, worueber eine Zeile ueberhaupt spricht:\n"
                 "#   scaling_of_dimensionless_equilibrium\n"
                 "#       exakte Umskalierung der nichtdimensionalen P3a/P3b-Loesung; gamma\n"
                 "#       tritt dort nur ueber die Druckskala gamma/a auf.\n"
                 "#   dimensionless_group_at_fixed_field\n"
                 "#       eine dimensionslose Gruppe, in der gamma explizit steht, bei\n"
                 "#       festgehaltenem FELD ausgewertet.\n"
                 "#   invariant_at_fixed_geometry_voltage_permittivity\n"
                 "#       skaliert GAR NICHT mit gamma.  Die Maxwell-Traktion eps0 E^2/2 ist\n"
                 "#       bei festgehaltener Geometrie, Spannung und Permittivitaetsverteilung\n"
                 "#       ein Funktional allein des Feldes, und in der Laplace-/Poisson-\n"
                 "#       Gleichung kommt gamma nicht vor.  Exponent 0, Faktor exakt 1.\n"
                 "#       Eine fruehere Fassung dieser Datei fuehrte genau diese Zeile\n"
                 "#       faelschlich als 'linear in gamma'.\n"
                 "#\n"
                 "# Die Skalierungstabelle steht in src/material_data.cpp\n"
                 "# (gamma_scaling_rows) und wird in tests/test_material_data.cpp geprueft.\n");
    std::fprintf(f, "quantity,category,law,exponent,recomputed,held_fixed,illustrative,"
                    "selected,band_lo,band_hi,factor_selected,factor_lo,factor_hi,note\n");
    if (gv.usable()) {
      const Real lo = g->min_at(T), hi = g->max_at(T), sel = gv.value;
      std::size_t n_rows = 0;
      const GammaScalingRow* rows = gamma_scaling_rows(n_rows);
      for (std::size_t i = 0; i < n_rows; ++i) {
        const GammaScalingRow& r = rows[i];
        std::fprintf(f, "%s,%s,%s,%.1f,%s,%s", r.quantity, to_string(r.category), r.law,
                     r.exponent, r.recomputed ? "yes" : "no",
                     csv_quote(r.what_is_held_fixed).c_str());
        put(f, ill.gamma);
        put(f, sel);
        put(f, lo);
        put(f, hi);
        put(f, std::pow(sel / ill.gamma, r.exponent));
        put(f, std::pow(lo / ill.gamma, r.exponent));
        put(f, std::pow(hi / ill.gamma, r.exponent));
        std::fprintf(f, ",%s\n", csv_quote(r.note).c_str());
      }
      say("");
      say("  gamma illustrativ " + std::to_string(ill.gamma) + " N/m, belegt " +
          std::to_string(sel) + " N/m: Faktor " + std::to_string(sel / ill.gamma) +
          " auf jede Druckskala und " + std::to_string(std::sqrt(sel / ill.gamma)) +
          " auf jede Spannung.");
      say("  Die Maxwell-Traktion bei festgehaltener Geometrie, Spannung und "
          "Permittivitaetsverteilung skaliert dagegen NICHT mit gamma: Exponent 0, "
          "Faktor 1.  Keine Zeile dieser Tabelle ist eine neue gekoppelte Rechnung.");
      say("  Belegtes Band bei " + std::to_string(T) + " K: " + std::to_string(lo) + " .. " +
          std::to_string(hi) + " N/m ueber " + std::to_string(g->n_sources_at(T)) +
          " Quellen; Streuung " + std::to_string(100.0 * g->relative_spread_at(T)) + " %.");
      if (!(ill.gamma >= lo && ill.gamma <= hi)) {
        say("  BEFUND: der bisher benutzte gamma-Wert liegt AUSSERHALB des gesamten "
            "belegten Bandes.");
      }
    } else {
      std::fprintf(f, "# gamma ist MissingMaterialData; es gibt nichts zu skalieren.\n");
      exit_code = 2;
    }
    std::fclose(f);
    say("  impact.csv geschrieben");
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_material (P2)\nphase=P2\ncommit=%s\n",
                 cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "substance=%s\ncas=%s\nT_K=%.9e\n", d.name, d.cas, T);
    std::fprintf(f, "database=NIST ILThermo (Standard Reference Database #147) v2.0 "
                    "und IoLiTec-Datenblatt IL-0006\n");
    for (std::size_t i = 0; i < d.n_properties; ++i) {
      const MaterialProperty& p = d.properties[i];
      const MaterialValue v = material_value(d, p.kind, T);
      std::fprintf(f, "%s_status=%s\n", to_string(p.kind), to_string(v.status));
      if (v.usable()) {
        std::fprintf(f, "%s_value=%.9e\n", to_string(p.kind), v.value);
        std::fprintf(f, "%s_spread=%.9e\n", to_string(p.kind), v.relative_spread);
      }
      std::fprintf(f, "%s_n_sources=%zu\n", to_string(p.kind), p.n_sources);
    }
    {
      const KinematicViscosityDerivation r = derive_kinematic_viscosity(d, T);
      std::fprintf(f, "kinematic_viscosity_derived=%s\n", r.ok ? "yes" : "no");
      if (r.ok) {
        std::fprintf(f, "kinematic_viscosity_derived_value=%.9e\n", r.value);
        std::fprintf(f, "kinematic_viscosity_derived_rel_uncertainty=%.9e\n",
                     r.relative_uncertainty);
        std::fprintf(f, "kinematic_viscosity_derived_status=derived\n");
      } else {
        std::fprintf(f, "kinematic_viscosity_derived_blocker=%s\n", r.blocker.c_str());
      }
    }
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_material: %s\n", e.what());
  return 2;
}
