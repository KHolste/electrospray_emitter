// es_emission_audit -- P5: the ion-emission contract and its blocker.
//
//   es_emission_audit <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN DOES.  It writes out the emission contract in its shipped
// state -- DISABLED and BLOCKED -- together with the audit of why, and a
// strictly DIMENSIONLESS sensitivity study of the quoted rate form.
//
// WHAT IT DOES NOT DO.  It produces no current, for either polarity, at any
// field.  Every path that could is exercised and every one fails closed; the
// run records which status each returned, so the absence is a measured result
// and not an omission.
//
// The sensitivity figure is dimensionless on purpose: quoting an absolute
// current would require the activation barrier, and there is no sourced one.
//
// Exit code 2 means a declared check failed -- which here means a path
// PRODUCED a number when it should not have.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/emission_contract.hpp"

using namespace es;
using constants::e;
using constants::kB;

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.empty()) {
    std::printf("es_emission_audit -- P5: Emissionsvertrag (blockiert)\n\n"
                "  es_emission_audit <ausgabeverzeichnis> [key=value ...]\n");
    return 1;
  }
  Config cfg;
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);
  int exit_code = 0;

  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& s) {
    std::printf("%s\n", s.c_str());
    std::fprintf(log, "%s\n", s.c_str());
  };

  say("P5 -- Emissionsvertrag.  Status: BLOCKED, Modell abgeschaltet.");
  say("");
  say("WARUM.  Ein Emissionsmodell darf nur aus einer vollstaendig geprueften");
  say("Gleichung mit belegten Parametern rechnen.  Beides fehlt:");
  say("  * Die Ratengleichung von Iribarne und Thomson wurde in diesem Lauf an");
  say("    KEINER Primaerquelle gelesen -- Vorfaktor, genaue Definition der");
  say("    Flaechenladungsdichte und Gueltigkeitsbereich sind unbelegt.  Eine");
  say("    gelesene Sekundaerquelle nennt das Modell fuer KLEINE Ionen gut");
  say("    gestuetzt und macht ausdruecklich einen Vorbehalt fuer groessere --");
  say("    also genau fuer den Fall ionischer Fluessigkeiten.");
  say("  * Fuer EMI-BF4 gibt es keine belegte Aktivierungsbarriere.  Der Wert");
  say("    1,09 eV in src/fluid.cpp hat keine Quelle, und die genannte Spanne");
  say("    1,0 bis 1,4 eV ist bei 298 K ein Faktor 1e7 im Strom.");
  say("");

  // --- every path that could produce a number ------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/contract.csv").c_str(), "w");
    std::fprintf(f, "# Jeder Pfad, der eine Zahl liefern koennte, und was er stattdessen\n"
                    "# meldet.  current_density ist ueberall null, weil kein Pfad rechnet --\n"
                    "# nicht, weil eine Rechnung null ergaebe.\n");
    std::fprintf(f, "case,polarity,enabled,equation_validated,E_n_V_per_m,status,"
                    "current_density_A_per_m2\n");
    auto row = [&](const char* tag, const EmissionModel& m, Real En) {
      const EmissionResult r = emitted_current_density(m, En);
      std::fprintf(f, "%s,%s,%s,%s,%.9e,%s,%.9e\n", tag, to_string(m.polarity),
                   m.enabled ? "yes" : "no", m.equation_validated ? "yes" : "no", En,
                   to_string(r.status), r.current_density);
      if (r.current_density != 0.0) {
        say(std::string("  FEHLER: der Pfad '") + tag + "' hat eine Zahl geliefert.");
        exit_code = 2;
      }
      return r.status;
    };
    for (Polarity p : {Polarity::Positive, Polarity::Negative}) {
      const EmissionModel& m = emibf4_emission_blocked(p);
      m.print(stdout);
      m.print(log);
      row("ausgeliefert", m, (p == Polarity::Positive) ? 1.0e9 : -1.0e9);
      EmissionModel on = m;
      on.enabled = true;
      on.temperature = 298.15;
      row("eingeschaltet", on, (p == Polarity::Positive) ? 1.0e9 : -1.0e9);
      EmissionModel forced = on;
      forced.equation_validated = true;   // even this is not enough
      row("Gleichung_freigegeben_Barriere_fehlt", forced,
          (p == Polarity::Positive) ? 1.0e9 : -1.0e9);
    }
    std::fclose(f);
    say("  contract.csv geschrieben: kein Pfad liefert eine Zahl.");
  }

  // --- the dimensionless sensitivity ---------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/sensitivity.csv").c_str(), "w");
    std::fprintf(f, "# DIMENSIONSLOSE Sensitivitaet der zitierten Ratenform.  Es ist KEINE\n"
                    "# Betriebsprognose und keine Stoffaussage: x = E/E* mit\n"
                    "# E* = 4 pi eps0 dG^2 / e^3, und b = dG/(kT).  Die Kurve zeigt, wie\n"
                    "# steil die Form ist -- und damit, warum eine unbelegte Barriere keine\n"
                    "# Vorhersage traegt.\n");
    std::fprintf(f, "b,dG_eV_at_298K,x,j_over_j0\n");
    const Real T = 298.15;
    for (Real dG_eV : {0.8, 1.0, 1.09, 1.2, 1.4}) {
      const Real b = dG_eV * e / (kB * T);
      for (int k = 0; k <= 200; ++k) {
        const Real x = 0.2 + 1.0 * static_cast<Real>(k) / 200.0;
        std::fprintf(f, "%.9e,%.4f,%.9e,%.9e\n", b, dG_eV, x,
                     iribarne_thomson_dimensionless(x, b));
      }
    }
    std::fclose(f);

    std::FILE* g = std::fopen((outdir + "/barrier_leverage.csv").c_str(), "w");
    std::fprintf(g, "# Wie stark eine unbelegte Barriere durchschlaegt.  Bei FESTEM Feld\n"
                    "# wird die Rate zweier Barrieren verglichen; die Zahlen sind\n"
                    "# Verhaeltnisse und tragen deshalb keine unbelegte absolute Groesse.\n");
    std::fprintf(g, "dG_eV,E_star_V_per_m,rate_ratio_to_1.09eV_at_1e9_V_per_m\n");
    const Real ref = iribarne_thomson_rate(1.0e9, 1.09 * e, T);
    for (Real dG_eV : {0.8, 0.9, 1.0, 1.09, 1.2, 1.3, 1.4}) {
      std::fprintf(g, "%.4f,%.9e,%.9e\n", dG_eV, barrier_free_field(dG_eV * e),
                   iribarne_thomson_rate(1.0e9, dG_eV * e, T) / ref);
    }
    std::fclose(g);
    say("  sensitivity.csv und barrier_leverage.csv geschrieben (dimensionslos bzw. "
        "als Verhaeltnis).");
    say("");
    say("  Hebelwirkung der Barriere: von 1,0 auf 1,4 eV aendert sich die Rate bei");
    say("  festem Feld um den Faktor " +
        std::to_string(iribarne_thomson_rate(1.0e9, 1.0 * e, T) /
                       iribarne_thomson_rate(1.0e9, 1.4 * e, T)) + ".");
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_emission_audit (P5)\nphase=P5\nstatus=blocked\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "model_enabled=no\nequation_validated=no\n");
    std::fprintf(f, "barrier_sourced=no\nmass_sourced=yes\n");
    std::fprintf(f, "blocker=Ratengleichung nicht an einer Primaerquelle geprueft; "
                    "keine belegte Aktivierungsbarriere fuer EMI-BF4\n");
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  say("");
  say(exit_code == 0 ? "Kein Pfad hat eine Zahl geliefert -- der Vertrag haelt."
                     : "MINDESTENS EIN PFAD HAT EINE ZAHL GELIEFERT.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_emission_audit: %s\n", e.what());
  return 2;
}
