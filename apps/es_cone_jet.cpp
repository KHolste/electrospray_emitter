// es_cone_jet -- P8: the cone-jet contract and its status map.
//
//   es_cone_jet <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS.  A STATUS AND VALIDITY MAP, not a regime prediction.  It
// writes which sub-models a cone-jet computation would need and whether this
// project has them, and it evaluates the dimensionless numbers that are
// DEFINITIONS -- so they need no literature source and no unchecked prefactor.
//
// WHAT IT IS NOT.  There is no cone-jet current, no jet radius and no droplet
// diameter.  The mode is blocked for two independent reasons: the physics is
// not solved (no two-phase free-surface flow, no surface charge transport), and
// the empirical scaling is not adopted because its original equation, erratum,
// prefactor and validity range were not checked at the source.
//
// Exit code 2 means a declared check failed -- here that means the refusal did
// NOT hold.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "es/config.hpp"
#include "es/cone_jet_contract.hpp"
#include "es/constants.hpp"
#include "es/material_data.hpp"
#include "es/status.hpp"

using namespace es;

namespace {
void put(std::FILE* f, Real v) {
  if (std::isfinite(v))
    std::fprintf(f, ",%.9e", v);
  else
    std::fprintf(f, ",nan");
}
std::string q(const std::string& s) {
  std::string o = "\"";
  for (char c : s) o += (c == '"') ? std::string("\"\"") : std::string(1, c);
  return o + "\"";
}
}  // namespace

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.empty()) {
    std::printf("es_cone_jet -- P8: Cone-Jet-Vertrag (blockiert)\n\n"
                "  es_cone_jet <ausgabeverzeichnis> [key=value ...]\n");
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
  say("P8 -- Cone-Jet und Tropfenbetrieb.  Status: BLOCKED.");
  say("");
  say("STRIKT GETRENNT von der Pure-Ion-Emission (P5) und vom statischen");
  say("Meniskus (P3a/P3b).  Ein Cone-Jet HAT einen Jet; ein gepinnter statischer");
  say("Meniskus hat keinen.  Die beiden sind verschiedene Zustaende derselben");
  say("Fluessigkeit und lassen sich nicht auseinander berechnen.");
  say("");

  // --- the refusal, exercised ----------------------------------------------
  {
    bool threw = false;
    try {
      cone_jet_current();
    } catch (const NotImplementedInThisPhase& e) {
      threw = true;
      say(std::string("  ") + e.what());
    }
    if (!threw) {
      say("  FEHLER: die Verweigerung hat nicht gehalten.");
      exit_code = 2;
    }
  }

  // --- the requirement map --------------------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/requirements.csv").c_str(), "w");
    std::fprintf(f, "# Welches Teilmodell eine Cone-Jet-Rechnung braucht, und ob dieses\n"
                    "# Projekt es hat.  Die Verfuegbarkeit ist gegen die anderen Phasen\n"
                    "# gemessen, nicht behauptet.\n");
    std::fprintf(f, "requirement,available,provided_by,why\n");
    Index missing = 0;
    for (const ConeJetRequirement& r : cone_jet_requirements()) {
      std::fprintf(f, "%s,%s,%s,%s\n", q(r.name).c_str(), r.available ? "yes" : "no",
                   q(r.provided_by).c_str(), q(r.why).c_str());
      if (!r.available) ++missing;
    }
    std::fclose(f);
    say("");
    say("  requirements.csv: " + std::to_string(missing) + " Teilmodelle fehlen.");
  }

  // --- the dimensionless diagnosis -----------------------------------------
  const Real a = 0.5 * cfg.num("device.phi_2", 10.0e-6);
  const Real R = cfg.num("feed.channel_radius", 5.0e-6);
  const Real E = cfg.num("cone_jet.surface_field", 5.0e7);
  const Real T = cfg.num("material.temperature", 298.15);
  {
    std::FILE* f = std::fopen((outdir + "/diagnosis.csv").c_str(), "w");
    std::fprintf(f, "# DIMENSIONSLOSE DIAGNOSE ueber dem Volumenstrom.  Jede Kennzahl ist\n"
                    "# eine DEFINITION und braucht keine Literaturquelle; die Herleitung\n"
                    "# steht in include/es/cone_jet_contract.hpp.  Es ist eine Diagnose,\n"
                    "# welche Physik dominieren wuerde, und KEINE Regimevorhersage.\n"
                    "#\n"
                    "# eps_r ist als EINZELWERT weiterhin MissingMaterialData -- keine Quelle\n"
                    "# nennt Reinheit und Wassergehalt.  Was von eps_r NUR UEBER die\n"
                    "# Ladungsrelaxationszeit abhaengt, ist seit der P3-Korrektur trotzdem\n"
                    "# bestimmt: tau_e loest die implizite Gleichung\n"
                    "#     tau = eps0 eps_r(1/(2 pi tau)) / K\n"
                    "# auf der GEMESSENEN Dispersionskurve, und\n"
                    "#     r* = (gamma eps0^2 eps_r^2/(rho K^2))^(1/3) = (gamma tau^2/rho)^(1/3)\n"
                    "# haengt von eps_r ausschliesslich ueber tau ab.  Beide tragen deshalb\n"
                    "# einen Wert UND ein Band, und beides ist kein Ersatzwert.\n"
                    "# tau_e_self_consistent=yes markiert genau diese Herkunft.  Eine fruehere\n"
                    "# Fassung fuehrte beide Spalten als nan mit der Notiz 'eps_r fehlt'.\n");
    std::fprintf(f, "Q_m3_per_s,status,tau_e_s,t_capillary_s,t_viscous_s,t_residence_s,"
                    "r_star_m,Oh,Bo_E,Re,Ca,tau_e_lo_s,tau_e_hi_s,r_star_lo_m,r_star_hi_m,"
                    "tau_e_self_consistent,eps_r_at_f_star,f_star_Hz\n");
    for (int k = 0; k <= 60; ++k) {
      const Real Q = 1.0e-16 * std::pow(10.0, 4.0 * static_cast<Real>(k) / 60.0);
      const ConeJetDiagnosis c = diagnose_cone_jet(emibf4_sourced(), T, a, Q, R, E);
      std::fprintf(f, "%.9e,%s", Q, to_string(c.status));
      put(f, c.tau_e);
      put(f, c.t_capillary);
      put(f, c.t_viscous);
      put(f, c.t_residence);
      put(f, c.r_star);
      put(f, c.Oh);
      put(f, c.Bo_E);
      put(f, c.Re);
      put(f, c.Ca);
      put(f, c.tau_e_lo);
      put(f, c.tau_e_hi);
      put(f, c.r_star_lo);
      put(f, c.r_star_hi);
      std::fprintf(f, ",%s", c.tau_e_self_consistent ? "yes" : "no");
      put(f, c.eps_r_at_f_star > 0.0 ? c.eps_r_at_f_star
                                     : std::numeric_limits<Real>::quiet_NaN());
      put(f, c.f_star > 0.0 ? c.f_star : std::numeric_limits<Real>::quiet_NaN());
      std::fprintf(f, "\n");
    }
    std::fclose(f);

    // The same numbers against the SURFACE FIELD, which is the one quantity
    // P3b actually computes -- so this column is where the diagnosis touches
    // something this project has measured.
    std::FILE* g = std::fopen((outdir + "/electric_bond.csv").c_str(), "w");
    std::fprintf(g, "# Die elektrische Bondzahl ueber dem Oberflaechenfeld.  Bo_E ist genau\n"
                    "# das Verhaeltnis der beiden Terme, die P3b bilanziert:\n"
                    "#   Bo_E = (eps0 E^2 / 2) / (gamma / a).\n"
                    "# Bei Bo_E = 1 ist der Maxwell-Druck so gross wie die Kapillarskala.\n"
                    "# Das ist eine Ablesung, keine Stabilitaetsaussage.\n");
    std::fprintf(g, "E_V_per_m,Bo_E\n");
    for (int k = 0; k <= 100; ++k) {
      const Real Ee = 1.0e6 * std::pow(10.0, 2.0 * static_cast<Real>(k) / 100.0);
      std::fprintf(g, "%.9e,%.9e\n", Ee, electric_bond(Ee, a, 0.05401));
    }
    std::fclose(g);

    const ConeJetDiagnosis c = diagnose_cone_jet(emibf4_sourced(), T, a, 1.0e-13, R, E);
    c.print(stdout);
    c.print(log);
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_cone_jet (P8)\nphase=P8\nstatus=blocked\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "contact_radius_m=%.9e\nchannel_radius_m=%.9e\nsurface_field_V_per_m=%.9e\n",
                 a, R, E);
    std::fprintf(f, "correlation_adopted=no\n");
    std::fprintf(f, "blocker=keine Zweiphasenstroemung mit freier Oberflaeche, kein "
                    "Oberflaechenladungstransport, kein Zerfallsmodell; empirische "
                    "Skalierung nicht an der Quelle geprueft (Ganan-Calvo 1997 und Erratum "
                    "2000 nicht im Volltext erreichbar)\n");
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  say("");
  say(exit_code == 0 ? "Die Verweigerung haelt; es wurde kein Strom und kein "
                       "Tropfendurchmesser erzeugt."
                     : "MINDESTENS EINE DEKLARIERTE PRUEFUNG IST FEHLGESCHLAGEN.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_cone_jet: %s\n", e.what());
  return 2;
}
