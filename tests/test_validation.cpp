// tests/test_validation.cpp -- P9: axisymmetric vs 3D, and the import contract
//
//   1. A result carries its geometry kind, and an axisymmetric one CANNOT be
//      read as a three-dimensional one.  Nothing in this project can produce a
//      ThreeDimensional label -- the test shows the path is closed.
//   2. The REVOLUTION REFERENCE: the 2 pi r weighting used in every
//      axisymmetric integral here is checked against an explicit 3D quadrature,
//      and the revolved field is checked for rotational invariance.
//   3. The import contract fails closed on a missing unit, uncertainty,
//      provenance or geometry kind -- and rejects the SET, not the point.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/validation.hpp"

using namespace es;
using constants::pi;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  [%s] %s\n", ok ? "ok" : "FEHLER", what.c_str());
  if (!ok) ++failures;
}

void check_rel(Real got, Real want, Real tol, const std::string& what) {
  const Real x =
      (std::abs(want) > 0.0) ? std::abs(got - want) / std::abs(want) : std::abs(got - want);
  std::printf("  [%s] %s: %.9e gegen %.9e, rel. %.3e (Grenze %.1e)\n",
              (x <= tol) ? "ok" : "FEHLER", what.c_str(), got, want, x, tol);
  if (!(x <= tol)) ++failures;
}

MeasuredPoint good_point() {
  MeasuredPoint p;
  p.quantity = "Extraktionsspannung";
  p.value = 1250.0;
  p.unit = "V";
  p.uncertainty = 15.0;
  p.uncertainty_type = UncertaintyType::TypeB;
  p.coverage_factor = 2.0;
  p.provenance = "Beispieleintrag, keine Messung -- er zeigt nur die Pflichtfelder";
  p.geometry = GeometryKind::ThreeDimensional;
  p.geometry_stated = true;
  p.conditions = "T = 298 K, positive Polaritaet, Q nicht angegeben";
  return p;
}

}  // namespace

int main() {
  std::printf("P9 -- achsensymmetrisch gegen 3D, und der Importvertrag\n\n");

  // =========================================================================
  std::printf("1. Ein Ergebnis traegt seine Geometrieart, und sie ist nicht aenderbar\n");
  {
    const LabelledResult ax(0.3439, GeometryKind::Axisymmetric, "h/a", "-");
    const LabelledResult rev(7.79e-7, GeometryKind::RevolvedAxisymmetric,
                             "integrierte Maxwell-Kraft", "N");
    ax.print(stdout);
    rev.print(stdout);
    for (const LabelledResult* r : {&ax, &rev}) {
      bool threw = false;
      std::string msg;
      try {
        (void)r->value_as_three_dimensional();
      } catch (const std::exception& e) {
        threw = true;
        msg = e.what();
      }
      check(threw, std::string(to_string(r->kind())) +
                       " kann nicht als dreidimensionales Ergebnis gelesen werden");
      check(msg.find(r->quantity()) != std::string::npos,
            "und die Meldung nennt die Groesse beim Namen");
    }
    // The only way to get a 3D value is to have computed one -- and nothing
    // here does.  The test constructs one by hand to show the gate opens for a
    // real 3D result and only for that.
    const LabelledResult td(1.0, GeometryKind::ThreeDimensional, "nur zum Test", "-");
    check(td.value_as_three_dimensional() == 1.0,
          "ein echtes 3D-Ergebnis kommt durch das Gatter -- und nur ein solches");
    check(!is_three_dimensional(GeometryKind::RevolvedAxisymmetric),
          "eine Rotationsreferenz ist ausdruecklich KEINE 3D-Rechnung");
  }

  // =========================================================================
  std::printf("\n2. Die Rotationsreferenz: 2 pi r IST das 3D-Integral\n");
  {
    // A hemisphere of radius R: the meridian is a quarter circle.
    const Real R = 5.0e-6;
    std::vector<Vec2> meridian;
    const int n = 2001;
    for (int k = 0; k < n; ++k) {
      const Real th = 0.5 * pi * static_cast<Real>(k) / static_cast<Real>(n - 1);
      meridian.push_back({R * std::sin(th), R * std::cos(th)});
    }

    // (a) area: f = 1.  The closed form is 2 pi R^2.
    const RevolutionCheck a =
        revolution_surface_integral(meridian, [](Vec2) { return 1.0; }, 720);
    std::printf("    Flaeche: achsensym. %.9e, 3D %.9e, Unterschied %.3e; exakt %.9e\n",
                a.axisymmetric, a.three_dimensional, a.relative_difference,
                2.0 * pi * R * R);
    // The bound is set by the ROUND-OFF of the summation, not by machine
    // epsilon: the 3D quadrature adds 720 x 2000 terms where the axisymmetric
    // form adds 2000, so a few times 1e-13 is what a double can hold.  A
    // tighter bound would be measuring the addition order.
    check(a.relative_difference < 1.0e-11,
          "die achsensymmetrische Form und die explizite 3D-Quadratur stimmen bis auf die "
          "Summationsrundung ueberein");
    check_rel(a.axisymmetric, 2.0 * pi * R * R, 1.0e-6,
              "und beide treffen die geschlossene Form 2 pi R^2");

    // (b) a non-constant integrand, so that the agreement is not an accident of
    //     f = 1: f = z / R, whose integral over the hemisphere is pi R^2.
    const RevolutionCheck b = revolution_surface_integral(
        meridian, [R](Vec2 x) { return x.z / R; }, 360);
    std::printf("    int (z/R) dA: achsensym. %.9e, 3D %.9e, Unterschied %.3e; exakt %.9e\n",
                b.axisymmetric, b.three_dimensional, b.relative_difference, pi * R * R);
    check(b.relative_difference < 1.0e-11, "auch mit nicht konstantem Integranden");
    check_rel(b.axisymmetric, pi * R * R, 1.0e-6, "und gegen die geschlossene Form pi R^2");

    // (c) the volume: 2/3 pi R^3.
    const RevolutionCheck v = revolution_volume(meridian, 0.0, 360);
    std::printf("    Volumen: achsensym. %.9e, 3D %.9e, Unterschied %.3e; exakt %.9e\n",
                v.axisymmetric, v.three_dimensional, v.relative_difference,
                2.0 / 3.0 * pi * R * R * R);
    check(v.relative_difference < 1.0e-11, "dasselbe fuer das Rotationsvolumen");
    check_rel(v.axisymmetric, 2.0 / 3.0 * pi * R * R * R, 1.0e-6,
              "und gegen (2/3) pi R^3");

    // (d) the azimuthal count must not matter -- if it did, the integrand would
    //     not be azimuthally invariant, i.e. the field would not be
    //     axisymmetric after all.
    const RevolutionCheck c3 =
        revolution_surface_integral(meridian, [](Vec2) { return 1.0; }, 3);
    const RevolutionCheck c9 =
        revolution_surface_integral(meridian, [](Vec2) { return 1.0; }, 997);
    check_rel(c3.three_dimensional, c9.three_dimensional, 1.0e-11,
              "3 und 997 Azimute geben dasselbe -- die Mittelpunktregel ist exakt fuer "
              "einen azimutunabhaengigen Integranden");
  }

  // =========================================================================
  std::printf("\n3. Der Importvertrag schlaegt geschlossen fehl\n");
  {
    const ImportResult ok = import_measurements({good_point()});
    check(is_usable(ok.status) && ok.points.size() == 1,
          "ein vollstaendiger Punkt wird importiert");

    struct Case {
      const char* what;
      ImportStatus want;
      std::function<void(MeasuredPoint&)> break_it;
    };
    // THE HARD REQUIREMENTS.  Each of them means the record cannot be
    // interpreted at all, and each rejects the whole set.
    const Case cases[] = {
        {"Einheit", ImportStatus::MissingUnit, [](MeasuredPoint& p) { p.unit.clear(); }},
        {"Fundstelle", ImportStatus::MissingProvenance,
         [](MeasuredPoint& p) { p.provenance.clear(); }},
        {"Geometrieart", ImportStatus::MissingGeometryKind,
         [](MeasuredPoint& p) { p.geometry_stated = false; }},
    };
    for (const Case& c : cases) {
      MeasuredPoint p = good_point();
      c.break_it(p);
      const ImportResult r = import_measurements({p});
      std::printf("    ohne %-16s -> %s\n", c.what, to_string(r.status));
      check(r.status == c.want, std::string("fehlende ") + c.what + " wird abgefangen");
      check(is_hard_error(r.status), "und zwar als HARTER Fehler");
      check(r.points.empty(), "und es wird NICHTS importiert");
    }

    // THE UNCERTAINTY IS NOT ONE OF THEM.  An earlier version of this contract
    // treated it the same way and threw real measurements away: a publication
    // that reports a current without an error bar has produced an INCOMPLETE
    // record, not a broken one.  It is imported, archived, may be drawn -- and
    // carries no quantitative claim.
    std::printf("    -- eine nicht berichtete Unsicherheit ist KEIN harter Fehler --\n");
    {
      const Case soft[] = {
          {"Unsicherheitswert", ImportStatus::OkUncertaintyNotReported,
           [](MeasuredPoint& p) { p.uncertainty = 0.0; }},
          {"Unsicherheitstyp", ImportStatus::OkUncertaintyNotReported,
           [](MeasuredPoint& p) { p.uncertainty_type = UncertaintyType::NotReported; }},
      };
      for (const Case& c : soft) {
        MeasuredPoint p = good_point();
        c.break_it(p);
        const ImportResult r = import_measurements({p});
        std::printf("    ohne %-16s -> %s\n", c.what, to_string(r.status));
        check(r.status == c.want,
              std::string("ohne ") + c.what + ": der Zustand heisst NotReported und nicht "
              "'fehlt'");
        check(is_usable(r.status),
              std::string("ohne ") + c.what + ": der Punkt wird IMPORTIERT und archiviert");
        check(!usable_quantitatively(r.status),
              std::string("ohne ") + c.what +
                  ": aber er darf keine quantitative Validierung tragen");
        check(!is_hard_error(r.status), "und der Satz wird nicht als Ganzes verworfen");
        check(r.points.size() == 1 && r.n_qualitative_only == 1 && r.n_quantitative == 0,
              "die Zaehlung sagt, wie viel des Satzes eine Zahl tragen darf");
      }
    }

    // A MIXED SET stays whole and stays honest about which half is which.
    {
      MeasuredPoint no_unc = good_point();
      no_unc.quantity = "Strahlstrom";
      no_unc.unit = "A";
      no_unc.value = 2.1e-7;
      no_unc.uncertainty = 0.0;
      no_unc.uncertainty_type = UncertaintyType::NotReported;
      const ImportResult r = import_measurements({good_point(), no_unc});
      std::printf("    gemischter Satz -> %s, %lld quantitativ, %lld nur qualitativ\n",
                  to_string(r.status), static_cast<long long>(r.n_quantitative),
                  static_cast<long long>(r.n_qualitative_only));
      check(r.points.size() == 2,
            "ein Satz aus einem vollstaendigen und einem unsicherheitslosen Punkt wird "
            "GANZ importiert");
      check(r.n_quantitative == 1 && r.n_qualitative_only == 1,
            "und beide Haelften werden getrennt gezaehlt");
      check(r.status == ImportStatus::OkUncertaintyNotReported,
            "der Satzstatus verschweigt nicht, dass ein Punkt keine Unsicherheit traegt");
      check(usable_quantitatively(r.points[0].check()) &&
                !usable_quantitatively(r.points[1].check()),
            "jeder Punkt traegt seinen eigenen Status -- der Satz ist nicht homogen");
    }

    // The SET is rejected on a HARD error, not the point.
    MeasuredPoint bad = good_point();
    bad.provenance.clear();
    const ImportResult set = import_measurements({good_point(), bad, good_point()});
    check(!is_usable(set.status) && set.points.empty(),
          "ein Satz mit einem schlechten Punkt wird als GANZES abgelehnt");
    check(set.first_bad == 1, "und der schlechte Punkt wird benannt");
    check(set.message.find("GANZES") != std::string::npos,
          "und die Meldung sagt, warum nicht nur der eine Punkt weggelassen wird");

    // A unit mismatch within one quantity.
    MeasuredPoint mv = good_point();
    mv.unit = "kV";
    const ImportResult mm = import_measurements({good_point(), mv});
    check(mm.status == ImportStatus::UnitMismatch,
          "zwei Einheiten fuer dieselbe Groesse werden abgefangen");
  }

  // =========================================================================
  // THE VALIDATION MATRIX -- six independent verdicts per row.
  //
  // The earlier version carried ONE comparability per row and the figure
  // coloured the row by it, so a quantity that is comparable in principle but
  // blocked appeared green.  The total current is the clearest case: directly
  // comparable, and not computable by this project at all.
  std::printf("\n4. Die Validierungsmatrix: sechs Fragen, nicht eine\n");
  {
    const std::vector<ValidationEntry> m = validation_matrix();
    Index direct = 0, reduced = 0, none = 0, computed = 0;
    std::printf("    %-40s %-4s %-4s %-4s %-4s %-4s %s\n", "Groesse", "geo", "impl", "konv",
                "vgl", "VAL", "blockiert");
    for (const ValidationEntry& e : m) {
      switch (e.comparability) {
        case Comparability::Direct: ++direct; break;
        case Comparability::AfterStatedReduction: ++reduced; break;
        case Comparability::NotComparable: ++none; break;
      }
      if (std::string(e.computed_by) != "-") ++computed;
      std::printf("    %-40s %-4s %-4s %-4s %-4s %-4s %s\n", e.quantity,
                  symbol(as_assessment(e.comparability)), symbol(e.implemented),
                  symbol(e.converged), symbol(e.comparable_with_data), symbol(e.validated),
                  e.blocked ? "JA" : "-");
    }
    check(m.size() >= 10, "die Matrix ist nicht bloss eine Geste");
    check(none >= 3,
          "es gibt Groessen, die achsensymmetrisch grundsaetzlich nicht darstellbar sind -- "
          "und sie stehen benannt in der Matrix");
    check(computed < static_cast<Index>(m.size()),
          "und Groessen, die dieses Projekt gar nicht rechnet");

    // --- every row is internally consistent -------------------------------
    for (const ValidationEntry& e : m) {
      const std::string bad = inconsistency(e);
      check(bad.empty(), std::string("'") + e.quantity + "': die Zeile ist in sich stimmig" +
                             (bad.empty() ? "" : " -- " + bad));
      check(e.blocked == (e.blocked_reason[0] != '\0'),
            std::string("'") + e.quantity +
                "': blockiert genau dann, wenn ein Grund genannt ist");
      check(e.convergence_note[0] != '\0' && e.validation_note[0] != '\0',
            std::string("'") + e.quantity +
                "': zu Konvergenz und Validierung steht je ein Satz da, statt einer Farbe");
    }

    // --- THE POINT: comparability is not validation -----------------------
    const ValidationTally t = tally(m);
    std::printf("    vergleichbar %lld, implementiert %lld, konvergiert %lld, "
                "mit Daten vergleichbar %lld, VALIDIERT %lld, blockiert %lld (von %lld)\n",
                static_cast<long long>(t.n_comparable),
                static_cast<long long>(t.n_implemented),
                static_cast<long long>(t.n_converged),
                static_cast<long long>(t.n_comparable_with_data),
                static_cast<long long>(t.n_validated), static_cast<long long>(t.n_blocked),
                static_cast<long long>(t.n_rows));
    check(t.n_validated == 0,
          "NICHTS ist validiert: es sind ueberhaupt keine Messdaten importiert, und die "
          "Matrix sagt das statt es durch eine Farbe zu verdecken");
    check(t.n_comparable > t.n_validated,
          "und es gibt Groessen, die vergleichbar und trotzdem nicht validiert sind -- "
          "genau die Unterscheidung, die die alte einachsige Karte verwischt hat");
    check(t.n_blocked > 0, "blockierte Groessen sind als solche gefuehrt");

    // The row that made the old figure wrong: directly comparable AND blocked.
    const ValidationEntry* current = nullptr;
    for (const ValidationEntry& e : m)
      if (std::string(e.quantity) == "Gesamtstrom") current = &e;
    check(current != nullptr, "der Gesamtstrom steht in der Matrix");
    if (current) {
      check(current->comparability == Comparability::Direct && current->blocked,
            "er ist DIREKT vergleichbar UND blockiert -- in der alten einachsigen Karte "
            "erschien er deshalb gruen");
      check(current->implemented == Assessment::No && current->validated == Assessment::No,
            "auf den eigenen Achsen steht jetzt, dass er weder gerechnet noch validiert ist");
      check(std::string(current->blocked_reason).find("P5") != std::string::npos,
            "und der Blockergrund benennt P5");
    }

    // The invariant is not vacuous: a hand-built row that claims validation
    // without the things it requires must be caught.
    {
      ValidationEntry fake = m.front();
      fake.validated = Assessment::Yes;
      fake.comparable_with_data = Assessment::No;
      check(!inconsistency(fake).empty(),
            "eine Zeile, die Validierung ohne Vergleichbarkeit mit Daten behauptet, wird "
            "abgefangen -- die Invariante ist nicht leer");
      ValidationEntry fake2 = m.front();
      fake2.validated = Assessment::Yes;
      fake2.implemented = Assessment::No;
      check(!inconsistency(fake2).empty(),
            "und ebenso eine, die Validierung ohne Implementierung behauptet");
    }
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
