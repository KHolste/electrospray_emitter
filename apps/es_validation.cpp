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
  //
  // SIX INDEPENDENT VERDICTS PER ROW.  The earlier one-dimensional version was
  // coloured by comparability alone, which made the total current -- directly
  // comparable, and not computable here at all -- look like a success.
  {
    std::FILE* f = std::fopen((outdir + "/validation_matrix.csv").c_str(), "w");
    std::fprintf(f,
                 "# Die Validierungsmatrix.  Sie beantwortet SECHS unabhaengige Fragen je\n"
                 "# Groesse und gerade nicht eine einzige:\n"
                 "#   comparable_geometry  -- laesst sich die Groesse zwischen einer\n"
                 "#                           achsensymmetrischen Rechnung und einem\n"
                 "#                           3D-Geraet ueberhaupt vergleichen?\n"
                 "#   implemented          -- rechnet dieses Projekt sie?\n"
                 "#   converged            -- ist das numerische Ergebnis nach einem VORAB\n"
                 "#                           festgelegten Kriterium konvergiert?\n"
                 "#   comparable_with_data -- liesse sie sich mit einer Messung vergleichen?\n"
                 "#   validated            -- ist sie TATSAECHLICH mit Messdaten verglichen\n"
                 "#                           worden und hat innerhalb der angegebenen\n"
                 "#                           Unsicherheiten uebereingestimmt?\n"
                 "#   blocked + reason     -- ist sie blockiert, und wodurch?\n"
                 "#\n"
                 "# INVARIANTE, im Code geprueft statt in Prosa versprochen: validated=yes\n"
                 "# verlangt implemented=yes UND converged=yes (oder n/a) UND\n"
                 "# comparable_with_data=yes UND nicht blockiert.  Eine Groesse kann NIE als\n"
                 "# validiert erscheinen, nur weil sie vergleichbar ist.\n"
                 "#\n"
                 "# Werte je Achse: yes / partial / no / n-a.\n");
    std::fprintf(f, "quantity,unit,comparability,comparable_geometry,implemented,converged,"
                    "comparable_with_data,validated,blocked,blocked_reason,computed_by,"
                    "phase_status,condition,convergence_note,validation_note,measurable\n");
    const std::vector<ValidationEntry> vm = validation_matrix();
    for (const ValidationEntry& e : vm) {
      const std::string bad = inconsistency(e);
      if (!bad.empty()) {
        say(std::string("  FEHLER in der Matrixzeile '") + e.quantity + "': " + bad);
        exit_code = 2;
      }
      std::fprintf(f, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                   q(e.quantity).c_str(), q(e.unit).c_str(), to_string(e.comparability),
                   to_string(as_assessment(e.comparability)), to_string(e.implemented),
                   to_string(e.converged), to_string(e.comparable_with_data),
                   to_string(e.validated), e.blocked ? "yes" : "no",
                   q(e.blocked_reason).c_str(), q(e.computed_by).c_str(),
                   q(e.phase_status).c_str(), q(e.condition).c_str(),
                   q(e.convergence_note).c_str(), q(e.validation_note).c_str(),
                   e.measurable ? "yes" : "no");
    }
    std::fclose(f);

    const ValidationTally t = tally(vm);
    say("  Validierungsmatrix, " + std::to_string(t.n_rows) + " Zeilen:");
    say("    vergleichbar (ganz oder nach Reduktion): " + std::to_string(t.n_comparable));
    say("    implementiert:                          " + std::to_string(t.n_implemented));
    say("    numerisch konvergiert:                  " + std::to_string(t.n_converged));
    say("    mit Messdaten vergleichbar:             " +
        std::to_string(t.n_comparable_with_data));
    say("    TATSAECHLICH VALIDIERT:                 " + std::to_string(t.n_validated));
    say("    blockiert:                              " + std::to_string(t.n_blocked));
    say("  Der Abstand zwischen der ersten und der vorletzten Zahl ist der Punkt dieser "
        "Tabelle.");
    // Nothing here is validated, and that is a checked statement rather than a
    // hope: no measured data have been imported at all.
    if (t.n_validated != 0) {
      say("  FEHLER: eine Zeile behauptet eine Validierung, obwohl keine Messdaten "
          "importiert sind.");
      exit_code = 2;
    }
    {
      std::FILE* g = std::fopen((outdir + "/validation_tally.csv").c_str(), "w");
      std::fprintf(g, "# Wie viele Zeilen jede Achse erreichen.  Die vorletzte Zahl ist\n"
                      "# null, und das ist der Befund von P9.\n");
      std::fprintf(g, "axis,n_rows,n_reached\n");
      std::fprintf(g, "comparable_geometry,%lld,%lld\n", (long long)t.n_rows,
                   (long long)t.n_comparable);
      std::fprintf(g, "implemented,%lld,%lld\n", (long long)t.n_rows,
                   (long long)t.n_implemented);
      std::fprintf(g, "converged,%lld,%lld\n", (long long)t.n_rows,
                   (long long)t.n_converged);
      std::fprintf(g, "comparable_with_data,%lld,%lld\n", (long long)t.n_rows,
                   (long long)t.n_comparable_with_data);
      std::fprintf(g, "validated,%lld,%lld\n", (long long)t.n_rows,
                   (long long)t.n_validated);
      std::fprintf(g, "blocked,%lld,%lld\n", (long long)t.n_rows, (long long)t.n_blocked);
      std::fclose(g);
    }
  }

  // --- the import contract, exercised --------------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/import_contract.csv").c_str(), "w");
    std::fprintf(f,
                 "# Der Importvertrag, an einem absichtlich unvollstaendigen Beispielsatz\n"
                 "# vorgefuehrt.  KEINER dieser Punkte ist eine Messung: sie zeigen, was ein\n"
                 "# Datensatz tragen muss und woran ein Import scheitert.\n"
                 "#\n"
                 "# HARTE FEHLER (der ganze Satz wird abgelehnt): fehlende Einheit, fehlende\n"
                 "# Fundstelle, fehlende Geometrieart, widerspruechliche Einheiten.  Jeder\n"
                 "# davon heisst, dass der Datensatz gar nicht zu deuten ist.\n"
                 "#\n"
                 "# KEIN harter Fehler: eine in der Publikation NICHT ANGEGEBENE\n"
                 "# Unsicherheit.  Das ist eine Tatsache ueber die Quelle, kein Defekt des\n"
                 "# Datensatzes.  Solche Punkte werden importiert, archiviert und duerfen mit\n"
                 "# sichtbarem Status qualitativ dargestellt werden -- aber sie tragen keine\n"
                 "# quantitative Validierung.  Eine fruehere Fassung dieses Vertrags lehnte\n"
                 "# sie hart ab und warf damit echte Messungen weg.\n");
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
        {"Unsicherheit_nicht_berichtet", [](MeasuredPoint& p) { p.uncertainty = 0.0; }},
        {"Unsicherheitstyp_nicht_berichtet",
         [](MeasuredPoint& p) { p.uncertainty_type = UncertaintyType::NotReported; }},
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
      const std::string tag = c.tag;
      const bool expected_quantitative = (tag == "vollstaendig");
      const bool expected_qualitative = (tag == "Unsicherheit_nicht_berichtet" ||
                                         tag == "Unsicherheitstyp_nicht_berichtet");
      if (expected_quantitative) {
        p.print(stdout);
        p.print(log);
        if (!usable_quantitatively(r.status)) {
          say("  FEHLER: der vollstaendige Punkt ist nicht quantitativ verwendbar.");
          exit_code = 2;
        }
      } else if (expected_qualitative) {
        // NOT a hard error any more, and that is the correction: a publication
        // that omits an error bar has produced an incomplete record, not a
        // broken one.  It is archived, may be drawn, and carries no number.
        p.print(stdout);
        p.print(log);
        if (!is_usable(r.status) || usable_quantitatively(r.status)) {
          say(std::string("  FEHLER: '") + c.tag +
              "' muss importierbar, aber nur qualitativ verwendbar sein.");
          exit_code = 2;
        }
      } else if (is_usable(r.status)) {
        say(std::string("  FEHLER: '") + c.tag + "' wurde importiert.");
        exit_code = 2;
      }
    }

    // A mixed set: the hard requirements hold for every point, one point has no
    // reported uncertainty.  The set is imported WHOLE, and the counts say how
    // much of it may carry a number.
    {
      MeasuredPoint a = base();
      MeasuredPoint b = base();
      b.quantity = "Strahlstrom";
      b.unit = "A";
      b.value = 2.1e-7;
      b.uncertainty = 0.0;
      b.uncertainty_type = UncertaintyType::NotReported;
      b.provenance = "BEISPIELEINTRAG -- keine Messung.  Eine Publikation ohne Fehlerbalken.";
      const ImportResult r = import_measurements({a, b});
      std::fprintf(f, "gemischter_Satz,%s,%.9e,%s,%.9e,%s,%s,%s\n",
                   q("2 Punkte, 1 ohne berichtete Unsicherheit").c_str(), b.value,
                   q(b.unit).c_str(), b.uncertainty, to_string(b.uncertainty_type),
                   to_string(b.geometry), to_string(r.status));
      say("  Gemischter Satz: " + std::to_string(r.n_quantitative) +
          " quantitativ verwendbar, " + std::to_string(r.n_qualitative_only) +
          " nur qualitativ -- und der Satz wird NICHT als Ganzes verworfen.");
      if (r.points.size() != 2 || r.n_quantitative != 1 || r.n_qualitative_only != 1) {
        say("  FEHLER: der gemischte Satz ist nicht wie erwartet importiert worden.");
        exit_code = 2;
      }
    }
    std::fclose(f);
    say("  import_contract.csv geschrieben.  Fehlende Einheit, Fundstelle oder "
        "Geometrieart sind harte Fehler; eine nicht berichtete Unsicherheit ist ein "
        "eigener Zustand und kein Grund, die Messung wegzuwerfen.");
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
