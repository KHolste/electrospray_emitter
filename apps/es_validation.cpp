// es_validation -- P9: the axisymmetric/3D separation and the validation scheme.
//
//   es_validation <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN WRITES.  The validation matrix, the revolution reference (the
// 2 pi r weighting of every axisymmetric integral in this project checked
// against an explicit 3D quadrature), and the import contract exercised on a
// deliberately incomplete example set.
//
// WHAT IT IS NOT.  There is no 3D mesh, no 3D solver and no 3D result.  The
// revolution reference solves nothing new; it checks a weighting.
//
// Exit code 2 means a declared check failed.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/validation.hpp"

using namespace es;
using constants::pi;

namespace {
std::string q(const std::string& s) {
  std::string o = "\"";
  for (char c : s) o += (c == '"') ? std::string("\"\"") : std::string(1, c);
  return o + "\"";
}
}  // namespace

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.empty()) {
    std::printf("es_validation -- P9: 2D/3D-Trennung und Validierungsschema\n\n"
                "  es_validation <ausgabeverzeichnis> [key=value ...]\n");
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
  say("P9 -- achsensymmetrisch gegen 3D, Validierungsschema, Importvertrag");
  say("");
  say("ES GIBT KEIN 3D-NETZ, KEINEN 3D-LOESER UND KEIN 3D-ERGEBNIS.  Das Etikett");
  say("ThreeDimensional existiert, damit nachweislich nichts hier eines erzeugt.");

  // --- the label gate -------------------------------------------------------
  {
    const LabelledResult ax(0.343960, GeometryKind::Axisymmetric, "h/a bei 1400 V", "-");
    bool threw = false;
    try {
      (void)ax.value_as_three_dimensional();
    } catch (const std::exception& e) {
      threw = true;
      say("");
      say(std::string("  ") + e.what());
    }
    if (!threw) {
      say("  FEHLER: ein achsensymmetrisches Ergebnis kam als 3D-Ergebnis durch.");
      exit_code = 2;
    }
  }

  // --- the revolution reference --------------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/revolution.csv").c_str(), "w");
    std::fprintf(f, "# Rotationsreferenz.  Dieselbe Groesse einmal in der\n"
                    "# achsensymmetrischen Form int f 2 pi r ds und einmal als explizite\n"
                    "# 3D-Quadratur ueber den Azimut.  Der Unterschied ist die\n"
                    "# Summationsrundung, nicht ein Modellunterschied -- und genau das ist\n"
                    "# die Aussage: die 2 pi r-Wichtung IST das dreidimensionale Integral.\n");
    std::fprintf(f, "quantity,n_azimuth,axisymmetric,three_dimensional,relative_difference,"
                    "closed_form\n");
    const Real R = cfg.num("revolution.radius", 5.0e-6);
    std::vector<Vec2> meridian;
    const int n = cfg.integer("revolution.n_meridian", 2001);
    for (int k = 0; k < n; ++k) {
      const Real th = 0.5 * pi * static_cast<Real>(k) / static_cast<Real>(n - 1);
      meridian.push_back({R * std::sin(th), R * std::cos(th)});
    }
    Real worst = 0.0;
    for (int na : {3, 12, 60, 360, 997}) {
      const RevolutionCheck a =
          revolution_surface_integral(meridian, [](Vec2) { return 1.0; }, na);
      std::fprintf(f, "Flaeche,%d,%.12e,%.12e,%.3e,%.12e\n", na, a.axisymmetric,
                   a.three_dimensional, a.relative_difference, 2.0 * pi * R * R);
      const RevolutionCheck b =
          revolution_surface_integral(meridian, [R](Vec2 x) { return x.z / R; }, na);
      std::fprintf(f, "Integral z/R,%d,%.12e,%.12e,%.3e,%.12e\n", na, b.axisymmetric,
                   b.three_dimensional, b.relative_difference, pi * R * R);
      const RevolutionCheck v = revolution_volume(meridian, 0.0, na);
      std::fprintf(f, "Volumen,%d,%.12e,%.12e,%.3e,%.12e\n", na, v.axisymmetric,
                   v.three_dimensional, v.relative_difference, 2.0 / 3.0 * pi * R * R * R);
      worst = std::max({worst, a.relative_difference, b.relative_difference,
                        v.relative_difference});
    }
    std::fclose(f);
    char buf[200];
    std::snprintf(buf, sizeof buf,
                  "  Rotationsreferenz: groesster Unterschied zwischen der 2 pi r-Form und "
                  "der 3D-Quadratur %.3e", worst);
    say("");
    say(buf);
    if (!(worst < 1.0e-11)) exit_code = 2;
  }

  // --- the validation matrix ------------------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/validation_matrix.csv").c_str(), "w");
    std::fprintf(f, "# Welche Groessen sich zwischen einer achsensymmetrischen Rechnung und\n"
                    "# einem dreidimensionalen Geraet ueberhaupt vergleichen lassen, unter\n"
                    "# welcher Bedingung, und was dieses Projekt davon rechnet.\n");
    std::fprintf(f, "quantity,unit,comparability,condition,computed_by,status,measurable\n");
    Index n_direct = 0, n_reduced = 0, n_none = 0;
    for (const ValidationEntry& e : validation_matrix()) {
      std::fprintf(f, "%s,%s,%s,%s,%s,%s,%s\n", q(e.quantity).c_str(), q(e.unit).c_str(),
                   to_string(e.comparability), q(e.condition).c_str(),
                   q(e.computed_by).c_str(), q(e.status).c_str(),
                   e.measurable ? "yes" : "no");
      switch (e.comparability) {
        case Comparability::Direct: ++n_direct; break;
        case Comparability::AfterStatedReduction: ++n_reduced; break;
        case Comparability::NotComparable: ++n_none; break;
      }
    }
    std::fclose(f);
    say("  Validierungsmatrix: " + std::to_string(n_direct) + " direkt vergleichbar, " +
        std::to_string(n_reduced) + " nach ausgesprochener Reduktion, " +
        std::to_string(n_none) + " grundsaetzlich nicht.");
  }

  // --- the import contract, exercised --------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/import_contract.csv").c_str(), "w");
    std::fprintf(f, "# Der Importvertrag, an einem absichtlich unvollstaendigen Beispielsatz\n"
                    "# vorgefuehrt.  KEINER dieser Punkte ist eine Messung: sie zeigen, was\n"
                    "# ein Datensatz tragen muss und woran ein Import scheitert.\n");
    std::fprintf(f, "case,quantity,value,unit,uncertainty,uncertainty_type,geometry,"
                    "status\n");
    auto base = []() {
      MeasuredPoint p;
      p.quantity = "Extraktionsspannung";
      p.value = 1250.0;
      p.unit = "V";
      p.uncertainty = 15.0;
      p.uncertainty_type = UncertaintyType::TypeB;
      p.coverage_factor = 2.0;
      p.provenance = "BEISPIELEINTRAG -- keine Messung.  Er zeigt nur die Pflichtfelder.";
      p.geometry = GeometryKind::ThreeDimensional;
      p.geometry_stated = true;
      p.conditions = "T = 298 K, positive Polaritaet";
      return p;
    };
    struct C { const char* tag; std::function<void(MeasuredPoint&)> f; };
    const C cases[] = {
        {"vollstaendig", [](MeasuredPoint&) {}},
        {"ohne_Einheit", [](MeasuredPoint& p) { p.unit.clear(); }},
        {"ohne_Unsicherheit", [](MeasuredPoint& p) { p.uncertainty = 0.0; }},
        {"ohne_Unsicherheitstyp",
         [](MeasuredPoint& p) { p.uncertainty_type = UncertaintyType::NotStated; }},
        {"ohne_Fundstelle", [](MeasuredPoint& p) { p.provenance.clear(); }},
        {"ohne_Geometrieart", [](MeasuredPoint& p) { p.geometry_stated = false; }},
    };
    for (const C& c : cases) {
      MeasuredPoint p = base();
      c.f(p);
      const ImportResult r = import_measurements({p});
      std::fprintf(f, "%s,%s,%.9e,%s,%.9e,%s,%s,%s\n", c.tag, q(p.quantity).c_str(), p.value,
                   q(p.unit).c_str(), p.uncertainty, to_string(p.uncertainty_type),
                   p.geometry_stated ? to_string(p.geometry) : "NICHT ANGEGEBEN",
                   to_string(r.status));
      if (std::string(c.tag) == "vollstaendig") {
        p.print(stdout);
        p.print(log);
        if (!is_usable(r.status)) exit_code = 2;
      } else if (is_usable(r.status)) {
        say(std::string("  FEHLER: '") + c.tag + "' wurde importiert.");
        exit_code = 2;
      }
    }
    std::fclose(f);
    say("  import_contract.csv geschrieben: nur der vollstaendige Punkt wird importiert.");
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_validation (P9)\nphase=P9\nstatus=infrastructure_only\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "three_dimensional_solver=none\n");
    std::fprintf(f, "revolution_reference=yes\n");
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  say("");
  say(exit_code == 0 ? "Alle deklarierten Pruefungen dieses Laufs bestanden."
                     : "MINDESTENS EINE DEKLARIERTE PRUEFUNG IST FEHLGESCHLAGEN.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_validation: %s\n", e.what());
  return 2;
}
