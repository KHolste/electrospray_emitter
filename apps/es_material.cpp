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
    std::fprintf(f, "# Was eine Rechnung benutzen darf, und wie unsicher es ist.\n"
                    "# illustrative_value ist der bisher benutzte Wert ohne Quelle aus\n"
                    "# src/liquid.cpp; er bleibt illustrative.  relative_spread ist die\n"
                    "# LITERATURSTREUUNG und damit die ehrliche Unsicherheit -- sie ist\n"
                    "# groesser als jede einzelne Fehlerangabe.\n");
    std::fprintf(f, "property,unit,T_K,status,selected_value,band_lo,band_hi,relative_spread,"
                    "n_sources_total,n_sources_at_T,illustrative_value,"
                    "illustrative_vs_selected,selected_reference\n");
    const LiquidProperties ill = emibf4_illustrative();
    for (std::size_t i = 0; i < d.n_properties; ++i) {
      const MaterialProperty& p = d.properties[i];
      const MaterialValue v = material_value(d, p.kind, T);
      Real illv = std::numeric_limits<Real>::quiet_NaN();
      switch (p.kind) {
        case PropertyKind::SurfaceTension: illv = ill.gamma; break;
        case PropertyKind::Density: illv = ill.rho; break;
        case PropertyKind::DynamicViscosity: illv = ill.documented_only.mu; break;
        case PropertyKind::ElectricalConductivity: illv = ill.documented_only.K; break;
        case PropertyKind::RelativePermittivity: illv = ill.documented_only.eps_r; break;
        default: break;
      }
      std::fprintf(f, "%s,%s,%.9e,%s", to_string(p.kind), si_unit(p.kind), T,
                   to_string(v.status));
      put(f, v.usable() ? v.value : std::numeric_limits<Real>::quiet_NaN());
      put(f, p.min_at(T));
      put(f, p.max_at(T));
      put(f, p.relative_spread_at(T));
      std::fprintf(f, ",%zu,%zu", p.n_sources, p.n_sources_at(T));
      put(f, illv);
      put(f, (v.usable() && std::isfinite(illv) && std::abs(v.value) > 0.0)
                 ? (illv - v.value) / v.value
                 : std::numeric_limits<Real>::quiet_NaN());
      std::fprintf(f, ",%s\n", v.source ? csv_quote(v.source->reference).c_str() : "");
    }
    std::fclose(f);
    say("  summary.csv geschrieben");
  }

  // --- what this changes for P3a/P3b ---------------------------------------
  {
    const LiquidProperties ill = emibf4_illustrative();
    const MaterialProperty* g = d.find(PropertyKind::SurfaceTension);
    const MaterialValue gv = material_value(d, PropertyKind::SurfaceTension, T);
    std::FILE* f = std::fopen((outdir + "/impact.csv").c_str(), "w");
    std::fprintf(f, "# Was der belegte gamma-Wert an den P3a/P3b-Zahlen aendert.  gamma geht\n"
                    "# LINEAR in die Gleichgewichtsgleichung ein, die Spannung skaliert mit\n"
                    "# sqrt(gamma).  Das sind Skalierungen der vorhandenen Rechnung, keine\n"
                    "# neuen Rechnungen.\n");
    std::fprintf(f, "quantity,scaling,illustrative,selected,band_lo,band_hi,factor_selected,"
                    "factor_lo,factor_hi\n");
    if (gv.usable()) {
      const Real lo = g->min_at(T), hi = g->max_at(T), sel = gv.value;
      auto row = [&](const char* q, const char* sc, Real e) {
        std::fprintf(f, "%s,%s", q, sc);
        put(f, ill.gamma);
        put(f, sel);
        put(f, lo);
        put(f, hi);
        put(f, std::pow(sel / ill.gamma, e));
        put(f, std::pow(lo / ill.gamma, e));
        put(f, std::pow(hi / ill.gamma, e));
        std::fprintf(f, "\n");
      };
      row("pressure_scale_gamma_over_a", "linear in gamma", 1.0);
      row("delta_p_exit_for_fixed_Pi", "linear in gamma", 1.0);
      row("voltage_for_fixed_shape", "sqrt(gamma)", 0.5);
      row("maxwell_force_for_fixed_shape", "linear in gamma", 1.0);
      say("");
      say("  gamma illustrativ " + std::to_string(ill.gamma) + " N/m, belegt " +
          std::to_string(sel) + " N/m: Faktor " + std::to_string(sel / ill.gamma) +
          " auf jede Druckskala und " + std::to_string(std::sqrt(sel / ill.gamma)) +
          " auf jede Spannung.");
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
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_material: %s\n", e.what());
  return 2;
}
