#include "es/validation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::pi;

const char* to_string(GeometryKind k) {
  switch (k) {
    case GeometryKind::Axisymmetric: return "Axisymmetric";
    case GeometryKind::ThreeDimensional: return "ThreeDimensional";
    case GeometryKind::RevolvedAxisymmetric: return "RevolvedAxisymmetric";
  }
  return "?";
}

const char* explain(GeometryKind k) {
  switch (k) {
    case GeometryKind::Axisymmetric:
      return "Die Meridianhalbebene mit Rotationssymmetrie.  Alles, was dieses Projekt "
             "rechnet, ist von dieser Art.";
    case GeometryKind::ThreeDimensional:
      return "Eine echte dreidimensionale Geometrie.  NICHTS in diesem Projekt erzeugt eine; "
             "das Etikett existiert, damit der Weg dorthin nachweislich geschlossen ist.";
    case GeometryKind::RevolvedAxisymmetric:
      return "Ein achsensymmetrisches Ergebnis, das als 3D-Feld ausgewertet wurde.  Es ist "
             "immer noch achsensymmetrische Physik -- das Etikett verhindert, dass eine "
             "Rotationsreferenz als 3D-Rechnung ausgegeben wird.";
  }
  return "?";
}

Real LabelledResult::value_as_three_dimensional() const {
  if (!is_three_dimensional(kind_))
    throw std::runtime_error(
        "'" + quantity_ + "' ist ein Ergebnis der Art " + to_string(kind_) +
        " und darf nicht als dreidimensionales Ergebnis ausgegeben werden.  " +
        explain(kind_) +
        "  Ein achsensymmetrisches Modell sieht genau das nicht, worin sich eine reale "
        "Anordnung von ihm unterscheidet; die Zahl saehe richtig aus und beantwortete eine "
        "andere Frage.");
  return value_;
}

void LabelledResult::print(std::FILE* out) const {
  std::fprintf(out, "  %s = %.9g %s  [%s]\n", quantity_.c_str(), value_, unit_.c_str(),
               to_string(kind_));
}

// ---------------------------------------------------------------------------

Real RevolvedField::potential(Real x, Real y, Real z) const {
  const Real r = std::hypot(x, y);
  return potential_at(*mesh, *phi, Vec2{r, z});
}

void RevolvedField::field(Real x, Real y, Real z, Real* Ex, Real* Ey, Real* Ez) const {
  const Real r = std::hypot(x, y);
  Index i, j;
  Real xi, eta;
  if (!locate(*mesh, Vec2{r, z}, &i, &j, &xi, &eta)) {
    *Ex = *Ey = *Ez = 0.0;
    return;
  }
  const Vec2 e00 = field_recovered_at_node(*mesh, *phi, *eps_r, *active, i, j, 1.0);
  const Vec2 e10 = field_recovered_at_node(*mesh, *phi, *eps_r, *active, i + 1, j, 1.0);
  const Vec2 e11 = field_recovered_at_node(*mesh, *phi, *eps_r, *active, i + 1, j + 1, 1.0);
  const Vec2 e01 = field_recovered_at_node(*mesh, *phi, *eps_r, *active, i, j + 1, 1.0);
  const Real n0 = (1 - xi) * (1 - eta), n1 = xi * (1 - eta), n2 = xi * eta, n3 = (1 - xi) * eta;
  const Vec2 E = n0 * e00 + n1 * e10 + n2 * e11 + n3 * e01;
  // The radial component points along (x, y)/r.  On the axis it must vanish;
  // there is no direction for it to point in.
  if (r > 0.0) {
    *Ex = E.r * x / r;
    *Ey = E.r * y / r;
  } else {
    *Ex = 0.0;
    *Ey = 0.0;
  }
  *Ez = E.z;
}

Real RevolvedField::azimuthal_variation(Real r, Real z, int n_azimuth) const {
  if (n_azimuth < 2) return 0.0;
  Real lo = std::numeric_limits<Real>::max(), hi = -lo;
  for (int k = 0; k < n_azimuth; ++k) {
    const Real t = 2.0 * pi * static_cast<Real>(k) / static_cast<Real>(n_azimuth);
    const Real v = potential(r * std::cos(t), r * std::sin(t), z);
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  const Real scale = std::max(std::abs(lo), std::abs(hi));
  return (scale > 0.0) ? (hi - lo) / scale : (hi - lo);
}

// ---------------------------------------------------------------------------

RevolutionCheck revolution_surface_integral(const std::vector<Vec2>& meridian,
                                            const std::function<Real(Vec2)>& f,
                                            int n_azimuth) {
  RevolutionCheck c;
  c.n_azimuth = n_azimuth;
  c.n_meridian = static_cast<int>(meridian.size());
  if (meridian.size() < 2 || n_azimuth < 3) return c;

  // (a) the axisymmetric form: int f 2 pi r ds, trapezoidal along the meridian.
  for (std::size_t k = 0; k + 1 < meridian.size(); ++k) {
    const Vec2 a = meridian[k], b = meridian[k + 1];
    const Real ds = norm(b - a);
    const Real fa = f(a), fb = f(b);
    // int_a^b f 2 pi r ds with f and r both linear in the segment parameter.
    c.axisymmetric += 2.0 * pi * ds *
                      ((2.0 * a.r * fa + a.r * fb + b.r * fa + 2.0 * b.r * fb) / 6.0);
  }

  // (b) an explicit 3D quadrature: the same surface, swept over the azimuth,
  //     with the area element r dt ds.  The azimuthal integral is a midpoint
  //     rule, which is EXACT for an integrand independent of the azimuth -- and
  //     that independence is what is being checked on the other side.
  const Real dt = 2.0 * pi / static_cast<Real>(n_azimuth);
  for (int m = 0; m < n_azimuth; ++m) {
    for (std::size_t k = 0; k + 1 < meridian.size(); ++k) {
      const Vec2 a = meridian[k], b = meridian[k + 1];
      const Real ds = norm(b - a);
      const Real fa = f(a), fb = f(b);
      c.three_dimensional +=
          dt * ds * ((2.0 * a.r * fa + a.r * fb + b.r * fa + 2.0 * b.r * fb) / 6.0);
    }
  }
  c.relative_difference =
      (std::abs(c.axisymmetric) > 0.0)
          ? std::abs(c.three_dimensional - c.axisymmetric) / std::abs(c.axisymmetric)
          : std::abs(c.three_dimensional);
  return c;
}

RevolutionCheck revolution_volume(const std::vector<Vec2>& meridian, Real z_base,
                                  int n_azimuth) {
  RevolutionCheck c;
  c.n_azimuth = n_azimuth;
  c.n_meridian = static_cast<int>(meridian.size());
  if (meridian.size() < 2 || n_azimuth < 3) return c;

  // (a) the truncated-cone formula, which is the exact revolved volume of a
  //     polyline.
  for (std::size_t k = 0; k + 1 < meridian.size(); ++k) {
    const Real r0 = meridian[k].r, r1 = meridian[k + 1].r;
    const Real z0 = meridian[k].z - z_base, z1 = meridian[k + 1].z - z_base;
    c.axisymmetric += pi / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1) * (z0 - z1);
  }
  // (b) the same, integrated over the azimuth as (1/3)(r0^2+r0 r1+r1^2)(z0-z1)
  //     dt / 2 per wedge -- the 3D form with pi replaced by the summed dt/2.
  const Real dt = 2.0 * pi / static_cast<Real>(n_azimuth);
  for (int m = 0; m < n_azimuth; ++m)
    for (std::size_t k = 0; k + 1 < meridian.size(); ++k) {
      const Real r0 = meridian[k].r, r1 = meridian[k + 1].r;
      const Real z0 = meridian[k].z - z_base, z1 = meridian[k + 1].z - z_base;
      c.three_dimensional += dt / 6.0 * (r0 * r0 + r0 * r1 + r1 * r1) * (z0 - z1);
    }
  c.relative_difference =
      (std::abs(c.axisymmetric) > 0.0)
          ? std::abs(c.three_dimensional - c.axisymmetric) / std::abs(c.axisymmetric)
          : std::abs(c.three_dimensional);
  return c;
}

// ---------------------------------------------------------------------------

const char* to_string(UncertaintyType u) {
  switch (u) {
    case UncertaintyType::TypeA: return "TypeA";
    case UncertaintyType::TypeB: return "TypeB";
    case UncertaintyType::NotStated: return "NotStated";
  }
  return "?";
}

const char* to_string(ImportStatus s) {
  switch (s) {
    case ImportStatus::Ok: return "Ok";
    case ImportStatus::MissingUnit: return "MissingUnit";
    case ImportStatus::MissingUncertainty: return "MissingUncertainty";
    case ImportStatus::MissingProvenance: return "MissingProvenance";
    case ImportStatus::MissingGeometryKind: return "MissingGeometryKind";
    case ImportStatus::UnitMismatch: return "UnitMismatch";
  }
  return "?";
}

const char* explain(ImportStatus s) {
  switch (s) {
    case ImportStatus::Ok:
      return "Der Datensatz traegt Einheit, Unsicherheit mit Typ, Fundstelle und "
             "Geometrieart.";
    case ImportStatus::MissingUnit:
      return "Ohne Einheit ist der Wert eine Zahl.  Der Import schlaegt geschlossen fehl.";
    case ImportStatus::MissingUncertainty:
      return "Ohne Unsicherheit mit Typ (A oder B nach GUM) ist der Wert eine Anekdote.  "
             "Ein Vergleich mit einer Rechnung braucht beide Fehlerbalken.";
    case ImportStatus::MissingProvenance:
      return "Ohne Fundstelle ist der Wert ein Geruecht: er laesst sich nicht wiederfinden.";
    case ImportStatus::MissingGeometryKind:
      return "Ohne die Angabe, zu welcher Geometrie die Messung gehoert, laesst sie sich mit "
             "nichts vergleichen -- eine achsensymmetrische Rechnung und eine reale "
             "Anordnung sind nicht dasselbe.";
    case ImportStatus::UnitMismatch:
      return "Die Einheiten zweier Punkte derselben Groesse stimmen nicht ueberein.";
  }
  return "?";
}

ImportStatus MeasuredPoint::check() const {
  if (unit.empty()) return ImportStatus::MissingUnit;
  if (!(uncertainty > 0.0) || uncertainty_type == UncertaintyType::NotStated)
    return ImportStatus::MissingUncertainty;
  if (provenance.empty()) return ImportStatus::MissingProvenance;
  if (!geometry_stated) return ImportStatus::MissingGeometryKind;
  return ImportStatus::Ok;
}

void MeasuredPoint::print(std::FILE* out) const {
  std::fprintf(out, "  %-28s %.6g +- %.3g %s (%s", quantity.c_str(), value, uncertainty,
               unit.c_str(), to_string(uncertainty_type));
  if (coverage_factor > 0.0) std::fprintf(out, ", k = %.3g", coverage_factor);
  std::fprintf(out, ")  [%s]  %s\n", to_string(geometry), check() == ImportStatus::Ok
                                                              ? "importierbar"
                                                              : to_string(check()));
  if (!provenance.empty()) std::fprintf(out, "      %s\n", provenance.c_str());
  if (!conditions.empty()) std::fprintf(out, "      Bedingungen: %s\n", conditions.c_str());
}

ImportResult import_measurements(const std::vector<MeasuredPoint>& raw) {
  ImportResult r;
  for (std::size_t k = 0; k < raw.size(); ++k) {
    const ImportStatus s = raw[k].check();
    if (s != ImportStatus::Ok) {
      r.status = s;
      r.first_bad = static_cast<Index>(k);
      r.message = "Punkt " + std::to_string(k) + " ('" + raw[k].quantity +
                  "'): " + explain(s) +
                  "  Der Satz wird als GANZES abgelehnt: ein stillschweigend "
                  "weggelassener Punkt waere ein Vergleich mit einem anderen Datensatz.";
      return r;
    }
  }
  // A unit mismatch between two points of the same quantity.
  for (std::size_t a = 0; a < raw.size(); ++a)
    for (std::size_t b = a + 1; b < raw.size(); ++b)
      if (raw[a].quantity == raw[b].quantity && raw[a].unit != raw[b].unit) {
        r.status = ImportStatus::UnitMismatch;
        r.first_bad = static_cast<Index>(b);
        r.message = "'" + raw[a].quantity + "' erscheint mit den Einheiten '" + raw[a].unit +
                    "' und '" + raw[b].unit + "'.  " + explain(r.status);
        return r;
      }
  r.points = raw;
  r.message = explain(ImportStatus::Ok);
  return r;
}

// ---------------------------------------------------------------------------

const char* to_string(Comparability c) {
  switch (c) {
    case Comparability::Direct: return "Direct";
    case Comparability::AfterStatedReduction: return "AfterStatedReduction";
    case Comparability::NotComparable: return "NotComparable";
  }
  return "?";
}

std::vector<ValidationEntry> validation_matrix() {
  return {
      {"Extraktorspannung", "V", Comparability::Direct,
       "Eine angelegte Spannung ist eine Randbedingung und in beiden Geometrien dieselbe "
       "Zahl.", "P2a/P2b/P2c/P3b", "geprueft", true},
      {"Gesamtstrom", "A", Comparability::Direct,
       "Ein Gesamtstrom ist ein Integral ueber die ganze Anordnung; er hat in beiden "
       "Geometrien dieselbe Bedeutung.  Dieses Projekt sagt ihn NICHT voraus (P5 blockiert).",
       "-", "blocked (P5)", true},
      {"Transmission durch die Blende", "-", Comparability::AfterStatedReduction,
       "Vergleichbar, sobald gesagt ist, ueber welche Startverteilung gemittelt wurde.  In "
       "3D kann die Blende ausserdem exzentrisch sein, was achsensymmetrisch nicht "
       "darstellbar ist.", "P7", "validated_subset (Transportantwort)", true},
      {"integrierte Maxwell-Kraft", "N", Comparability::Direct,
       "Ein Flaechenintegral; die 2 pi r-Wichtung IST das 3D-Integral, und die "
       "Rotationsreferenz prueft genau das.", "P3b, P0", "DiscretizationNotConverged", false},
      {"Apexhoehe h/a", "-", Comparability::Direct,
       "Auf der Achse, wo beide Geometrien uebereinstimmen -- SOFERN die reale Anordnung "
       "achsensymmetrisch ist.", "P3a/P3b", "qualitativ (P0: 1-%-Ziel verfehlt)", false},
      {"Oberflaechenfeld am Apex", "V/m", Comparability::Direct,
       "Auf der Achse.  P6 hat allerdings gemessen, dass die Feldrekonstruktion dort nur "
       "erster Ordnung ist.", "P2b/P2c/P3b", "erster Ordnung achsennah (P6)", false},
      {"Strahldivergenz", "rad", Comparability::AfterStatedReduction,
       "Achsensymmetrisch ist sie ein einziger Winkel; in 3D ist sie eine Verteilung ueber "
       "den Azimut.  Vergleichbar nur nach einer ausgesprochenen Mittelung.", "P7",
       "validated_subset", true},
      {"Auftreffverteilung auf dem Extraktor", "A/m^2", Comparability::AfterStatedReduction,
       "Achsensymmetrisch ist sie eine Funktion von r allein.  Jede azimutale Struktur -- "
       "die in einer realen Anordnung genau die interessante ist -- fehlt.", "P7",
       "validated_subset", true},
      {"azimutale Asymmetrie des Strahls", "-", Comparability::NotComparable,
       "Sie IST die Groesse, die eine achsensymmetrische Rechnung nicht hat.  Sie kann nur "
       "gemessen und nur in 3D gerechnet werden.", "-", "nicht gerechnet", true},
      {"Versatz von Emitter und Blende", "m", Comparability::NotComparable,
       "Eine Exzentrizitaet bricht die Rotationssymmetrie.  Achsensymmetrisch ist sie nicht "
       "darstellbar, auch nicht naeherungsweise.", "-", "nicht gerechnet", true},
      {"Neigung des Emitters", "rad", Comparability::NotComparable,
       "Wie der Versatz: sie bricht die Symmetrie.", "-", "nicht gerechnet", true},
      {"Emitter-zu-Emitter-Uebersprechen im Array", "-", Comparability::NotComparable,
       "Ein Array ist nicht rotationssymmetrisch.  Ein einzelner achsensymmetrischer "
       "Emitter kann es grundsaetzlich nicht abbilden.", "-", "nicht gerechnet", true},
      {"Ladungsrelaxationszeit", "s", Comparability::Direct,
       "Eine Stoffgroesse, unabhaengig von der Geometrie.  Nicht berechenbar, weil eps_r "
       "fehlt.", "P3", "MissingMaterialData", true},
  };
}

}  // namespace es
