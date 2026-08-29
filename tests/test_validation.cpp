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
    const Case cases[] = {
        {"Einheit", ImportStatus::MissingUnit, [](MeasuredPoint& p) { p.unit.clear(); }},
        {"Unsicherheit", ImportStatus::MissingUncertainty,
         [](MeasuredPoint& p) { p.uncertainty = 0.0; }},
        {"Unsicherheitstyp", ImportStatus::MissingUncertainty,
         [](MeasuredPoint& p) { p.uncertainty_type = UncertaintyType::NotStated; }},
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
      check(r.points.empty(), "und es wird NICHTS importiert");
    }

    // The SET is rejected, not the point.
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
  std::printf("\n4. Die Validierungsmatrix\n");
  {
    const std::vector<ValidationEntry> m = validation_matrix();
    Index direct = 0, reduced = 0, none = 0, computed = 0;
    for (const ValidationEntry& e : m) {
      switch (e.comparability) {
        case Comparability::Direct: ++direct; break;
        case Comparability::AfterStatedReduction: ++reduced; break;
        case Comparability::NotComparable: ++none; break;
      }
      if (std::string(e.computed_by) != "-") ++computed;
      std::printf("    %-38s %-22s %s\n", e.quantity, to_string(e.comparability), e.status);
    }
    std::printf("    %lld direkt, %lld nach Reduktion, %lld nicht vergleichbar; "
                "%lld werden hier gerechnet\n",
                static_cast<long long>(direct), static_cast<long long>(reduced),
                static_cast<long long>(none), static_cast<long long>(computed));
    check(m.size() >= 10, "die Matrix ist nicht bloss eine Geste");
    check(none >= 3,
          "es gibt Groessen, die achsensymmetrisch grundsaetzlich nicht darstellbar sind -- "
          "und sie stehen benannt in der Matrix");
    check(computed < static_cast<Index>(m.size()),
          "und Groessen, die dieses Projekt gar nicht rechnet");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
