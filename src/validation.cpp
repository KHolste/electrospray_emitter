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
    case UncertaintyType::NotReported: return "NotReported";
  }
  return "?";
}

const char* explain(UncertaintyType u) {
  switch (u) {
    case UncertaintyType::TypeA:
      return "Aus der statistischen Auswertung wiederholter Beobachtungen gewonnen "
             "(GUM Typ A).";
    case UncertaintyType::TypeB:
      return "Auf andere Weise gewonnen: Geraetespezifikation, Kalibrierschein, "
             "Sachurteil (GUM Typ B).";
    case UncertaintyType::NotReported:
      return "Die Publikation gibt keine an.  Das ist eine Tatsache ueber die QUELLE und "
             "kein Fehler des Imports: der Punkt wird archiviert und darf qualitativ "
             "dargestellt werden, aber keine quantitative Validierung darf auf ihm ruhen.";
  }
  return "?";
}

const char* to_string(ImportStatus s) {
  switch (s) {
    case ImportStatus::Ok: return "Ok";
    case ImportStatus::OkUncertaintyNotReported: return "OkUncertaintyNotReported";
    case ImportStatus::MissingUnit: return "MissingUnit";
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
             "Geometrieart.  Nur mit diesem Status darf ein QUANTITATIVER Vergleich "
             "gerechnet werden.";
    case ImportStatus::OkUncertaintyNotReported:
      return "Einheit, Fundstelle und Geometrieart sind da; die Publikation gibt aber KEINE "
             "Unsicherheit an.  Der Punkt wird importiert und archiviert und darf mit "
             "sichtbarem Status qualitativ dargestellt werden.  Er traegt keine "
             "quantitative Validierung: usable_quantitatively() ist fuer ihn falsch, und "
             "jede Abweichung, jedes Chi-Quadrat und jedes Bestanden/Durchgefallen muss "
             "genau das abfragen.";
    case ImportStatus::MissingUnit:
      return "Ohne Einheit ist der Wert eine Zahl.  Der Import schlaegt geschlossen fehl.";
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
  // The HARD requirements first.  Each of them means the record cannot be
  // interpreted at all, and each rejects the whole set.
  if (unit.empty()) return ImportStatus::MissingUnit;
  if (provenance.empty()) return ImportStatus::MissingProvenance;
  if (!geometry_stated) return ImportStatus::MissingGeometryKind;
  // The uncertainty is NOT one of them.  A publication that omits it has
  // produced an incomplete record, not a broken one, and throwing the
  // measurement away would lose real data.  It gets its own state instead.
  if (!(uncertainty > 0.0) || uncertainty_type == UncertaintyType::NotReported)
    return ImportStatus::OkUncertaintyNotReported;
  return ImportStatus::Ok;
}

void MeasuredPoint::print(std::FILE* out) const {
  std::fprintf(out, "  %-28s %.6g +- %.3g %s (%s", quantity.c_str(), value, uncertainty,
               unit.c_str(), to_string(uncertainty_type));
  if (coverage_factor > 0.0) std::fprintf(out, ", k = %.3g", coverage_factor);
  std::fprintf(out, ")  [%s]  %s\n", to_string(geometry),
               usable_quantitatively(check())
                   ? "quantitativ verwendbar"
                   : (is_usable(check()) ? "nur qualitativ (Unsicherheit NICHT berichtet)"
                                         : to_string(check())));
  if (!provenance.empty()) std::fprintf(out, "      %s\n", provenance.c_str());
  if (!conditions.empty()) std::fprintf(out, "      Bedingungen: %s\n", conditions.c_str());
}

ImportResult import_measurements(const std::vector<MeasuredPoint>& raw) {
  ImportResult r;
  for (std::size_t k = 0; k < raw.size(); ++k) {
    const ImportStatus s = raw[k].check();
    if (is_hard_error(s)) {
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
  for (const MeasuredPoint& q : raw) {
    if (usable_quantitatively(q.check()))
      ++r.n_quantitative;
    else
      ++r.n_qualitative_only;
  }
  // The set's status is the worst one in it -- and it is NOT allowed to hide
  // the fact that some points carry no uncertainty.
  r.status = (r.n_qualitative_only > 0) ? ImportStatus::OkUncertaintyNotReported
                                        : ImportStatus::Ok;
  r.message = explain(r.status);
  if (r.n_qualitative_only > 0)
    r.message += "  " + std::to_string(r.n_qualitative_only) + " von " +
                 std::to_string(raw.size()) +
                 " Punkten tragen keine berichtete Unsicherheit und stehen deshalb nur "
                 "qualitativ zur Verfuegung; " + std::to_string(r.n_quantitative) +
                 " sind quantitativ verwendbar.";
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

const char* to_string(Assessment a) {
  switch (a) {
    case Assessment::No: return "no";
    case Assessment::Partial: return "partial";
    case Assessment::Yes: return "yes";
    case Assessment::NotApplicable: return "n/a";
  }
  return "?";
}

const char* symbol(Assessment a) {
  switch (a) {
    case Assessment::No: return "-";
    case Assessment::Partial: return "~";
    case Assessment::Yes: return "+";
    case Assessment::NotApplicable: return ".";
  }
  return "?";
}

Assessment as_assessment(Comparability c) {
  switch (c) {
    case Comparability::Direct: return Assessment::Yes;
    case Comparability::AfterStatedReduction: return Assessment::Partial;
    case Comparability::NotComparable: return Assessment::No;
  }
  return Assessment::No;
}

std::string inconsistency(const ValidationEntry& e) {
  // THE INVARIANT.  Validation is not implied by comparability, and a row may
  // never claim it without the three things it actually requires.
  if (e.validated == Assessment::Yes) {
    if (e.implemented != Assessment::Yes)
      return "als validiert gefuehrt, ohne implementiert zu sein";
    if (e.converged != Assessment::Yes && e.converged != Assessment::NotApplicable)
      return "als validiert gefuehrt, ohne numerisch konvergiert zu sein";
    if (e.comparable_with_data != Assessment::Yes)
      return "als validiert gefuehrt, ohne mit Messdaten vergleichbar zu sein";
    if (e.blocked) return "als validiert gefuehrt, obwohl blockiert";
  }
  if (e.blocked && e.blocked_reason[0] == '\0')
    return "blockiert, ohne dass ein Grund genannt ist";
  if (!e.blocked && e.blocked_reason[0] != '\0')
    return "nicht blockiert, traegt aber einen Blockergrund";
  if (e.implemented == Assessment::No && e.converged == Assessment::Yes)
    return "als konvergiert gefuehrt, ohne implementiert zu sein";
  if (e.condition[0] == '\0') return "ohne Bedingung zur Vergleichbarkeit";
  return "";
}

ValidationTally tally(const std::vector<ValidationEntry>& m) {
  ValidationTally t;
  t.n_rows = static_cast<Index>(m.size());
  for (const ValidationEntry& e : m) {
    const Assessment c = as_assessment(e.comparability);
    if (c == Assessment::Yes || c == Assessment::Partial) ++t.n_comparable;
    if (e.implemented == Assessment::Yes) ++t.n_implemented;
    if (e.converged == Assessment::Yes) ++t.n_converged;
    if (e.comparable_with_data == Assessment::Yes ||
        e.comparable_with_data == Assessment::Partial)
      ++t.n_comparable_with_data;
    if (e.validated == Assessment::Yes) ++t.n_validated;
    if (e.blocked) ++t.n_blocked;
  }
  return t;
}

// ---------------------------------------------------------------------------
// The matrix.
//
// SIX INDEPENDENT VERDICTS PER ROW.  The earlier one-dimensional version
// coloured a row by its comparability alone, which made the total current --
// directly comparable, and not computable by this project at all -- appear as a
// success.  Every axis is now answered on its own, and `validated` is answered
// last and separately.  Nothing in this project is validated: no measured data
// have been imported, and the test checks that the matrix says so.

std::vector<ValidationEntry> validation_matrix() {
  return {
      {"Extraktorspannung", "V", Comparability::Direct,
       "Eine angelegte Spannung ist eine Randbedingung und in beiden Geometrien dieselbe "
       "Zahl.",
       Assessment::Yes,            // implemented: it is an input every phase carries
       Assessment::NotApplicable,  // converged: an input does not converge
       Assessment::Yes, Assessment::No, false, "",
       "P2a/P2b/P2c/P3b", "geprueft",
       "Eine vorgegebene Randbedingung; es gibt nichts zu konvergieren.",
       "Es sind ueberhaupt keine Messdaten importiert -- P9 stellt nur den Vertrag dafuer "
       "auf.",
       true},

      {"Gesamtstrom", "A", Comparability::Direct,
       "Ein Gesamtstrom ist ein Integral ueber die ganze Anordnung; er hat in beiden "
       "Geometrien dieselbe Bedeutung.",
       Assessment::No, Assessment::No, Assessment::Yes, Assessment::No, true,
       "P5 ist blockiert: es gibt keine an einer Primaerquelle gepruefte Emissionsrate und "
       "kein belegtes Delta G fuer EMI-BF4.  Dieses Projekt sagt keinen Strom voraus.",
       "-", "blocked (P5)",
       "Nicht gerechnet, also nichts zu konvergieren.",
       "Blockiert; ausserdem sind keine Messdaten importiert.",
       true},

      {"Transmission durch die Blende", "-", Comparability::AfterStatedReduction,
       "Vergleichbar, sobald gesagt ist, ueber welche Startverteilung gemittelt wurde.  In "
       "3D kann die Blende ausserdem exzentrisch sein, was achsensymmetrisch nicht "
       "darstellbar ist.",
       Assessment::Yes, Assessment::No, Assessment::Partial, Assessment::No, false, "",
       "P7", "validated_subset (Transportantwort)",
       "Fuer den Integrator ist die Zeitordnung gemessen; fuer die Transmission SELBST gibt "
       "es keine Netz- oder Teilchenzahlstudie.  Nicht als konvergiert gefuehrt.",
       "Keine Messdaten importiert; ausserdem haengt die Zahl an einer Startverteilung, die "
       "ohne P5 nicht physikalisch ist.",
       true},

      {"integrierte Maxwell-Kraft", "N", Comparability::Direct,
       "Ein Flaechenintegral; die 2 pi r-Wichtung IST das 3D-Integral, und die "
       "Rotationsreferenz prueft genau das.",
       Assessment::Yes, Assessment::No, Assessment::No, Assessment::No, false, "",
       "P3b, P0", "DiscretizationNotConverged",
       "P0 hat das vorab gesetzte 1-%-Ziel gemessen und VERFEHLT: 4,6 bis 6,1 %.  Die "
       "Randsingularitaet setzt die Rate auf etwa 1+beta = 0,55.",
       "Nicht konvergiert und keine Messdaten; eine Kraft ist ausserdem auf einem Pruefstand "
       "nicht direkt messbar.",
       false},

      {"Apexhoehe h/a", "-", Comparability::Direct,
       "Auf der Achse, wo beide Geometrien uebereinstimmen -- SOFERN die reale Anordnung "
       "achsensymmetrisch ist.",
       Assessment::Yes, Assessment::No, Assessment::No, Assessment::No, false, "",
       "P3a/P3b", "qualitativ (P0: 1-%-Ziel verfehlt)",
       "Keine Apexhoehe traegt drei Stellen; bei 1000 V liegen h/a und E_n nicht einmal im "
       "asymptotischen Bereich, dort ist gar kein Fehler schaetzbar.",
       "Nicht konvergiert; und eine Apexhoehe im Betrieb zu messen ist selbst eine offene "
       "Aufgabe.",
       false},

      {"Oberflaechenfeld am Apex", "V/m", Comparability::Direct,
       "Auf der Achse.  P6 hat allerdings gemessen, dass die Feldrekonstruktion dort nur "
       "erster Ordnung ist.",
       Assessment::Yes, Assessment::No, Assessment::No, Assessment::No, false, "",
       "P2b/P2c/P3b", "erster Ordnung achsennah (P6)",
       "Genau auf der Achse ist die Rekonstruktion erster Ordnung -- die 2 pi r-Gewichtung "
       "der Zellvolumina macht sie dort unsymmetrisch.",
       "Nicht konvergiert; ein lokales Feld am Apex ist ausserdem nicht direkt messbar.",
       false},

      {"Strahldivergenz", "rad", Comparability::AfterStatedReduction,
       "Achsensymmetrisch ist sie ein einziger Winkel; in 3D ist sie eine Verteilung ueber "
       "den Azimut.  Vergleichbar nur nach einer ausgesprochenen Mittelung.",
       Assessment::Yes, Assessment::No, Assessment::Partial, Assessment::No, false, "",
       "P7", "validated_subset",
       "Wie die Transmission: der Integrator ist geprueft, die Groesse selbst hat keine "
       "Konvergenzstudie.",
       "Keine Messdaten importiert.",
       true},

      {"Auftreffverteilung auf dem Extraktor", "A/m^2", Comparability::AfterStatedReduction,
       "Achsensymmetrisch ist sie eine Funktion von r allein.  Jede azimutale Struktur -- "
       "die in einer realen Anordnung genau die interessante ist -- fehlt.",
       Assessment::Yes, Assessment::No, Assessment::Partial, Assessment::No, false, "",
       "P7", "validated_subset",
       "Keine Konvergenzstudie der Verteilung selbst.",
       "Keine Messdaten importiert; und ohne P5 traegt die Verteilung keinen Strom.",
       true},

      {"azimutale Asymmetrie des Strahls", "-", Comparability::NotComparable,
       "Sie IST die Groesse, die eine achsensymmetrische Rechnung nicht hat.  Sie kann nur "
       "gemessen und nur in 3D gerechnet werden.",
       Assessment::No, Assessment::NotApplicable, Assessment::No, Assessment::No, true,
       "Es gibt kein 3D-Netz und keinen 3D-Loeser.  ThreeDimensional ist ein Etikett, das "
       "nichts erzeugen kann -- und der Test haelt diesen Weg ausdruecklich geschlossen.",
       "-", "nicht gerechnet",
       "Nicht gerechnet.",
       "Grundsaetzlich nicht mit einer achsensymmetrischen Rechnung vergleichbar.",
       true},

      {"Versatz von Emitter und Blende", "m", Comparability::NotComparable,
       "Eine Exzentrizitaet bricht die Rotationssymmetrie.  Achsensymmetrisch ist sie nicht "
       "darstellbar, auch nicht naeherungsweise.",
       Assessment::No, Assessment::NotApplicable, Assessment::No, Assessment::No, true,
       "Kein 3D-Netz, kein 3D-Loeser.",
       "-", "nicht gerechnet", "Nicht gerechnet.",
       "Grundsaetzlich nicht vergleichbar.", true},

      {"Neigung des Emitters", "rad", Comparability::NotComparable,
       "Wie der Versatz: sie bricht die Symmetrie.",
       Assessment::No, Assessment::NotApplicable, Assessment::No, Assessment::No, true,
       "Kein 3D-Netz, kein 3D-Loeser.",
       "-", "nicht gerechnet", "Nicht gerechnet.",
       "Grundsaetzlich nicht vergleichbar.", true},

      {"Emitter-zu-Emitter-Uebersprechen im Array", "-", Comparability::NotComparable,
       "Ein Array ist nicht rotationssymmetrisch.  Ein einzelner achsensymmetrischer "
       "Emitter kann es grundsaetzlich nicht abbilden.",
       Assessment::No, Assessment::NotApplicable, Assessment::No, Assessment::No, true,
       "Kein 3D-Netz, kein 3D-Loeser, und ein einzelner Emitter waere auch mit einem nicht "
       "genug.",
       "-", "nicht gerechnet", "Nicht gerechnet.",
       "Grundsaetzlich nicht vergleichbar.", true},

      {"Ladungsrelaxationszeit", "s", Comparability::Direct,
       "Eine Stoffgroesse, unabhaengig von der Geometrie.",
       Assessment::Yes, Assessment::NotApplicable, Assessment::Yes, Assessment::No, false, "",
       "P3", "Band belegt, Einzelwert fehlt",
       "Geschlossen loesbar, kein Netz beteiligt: tau folgt aus der impliziten Gleichung "
       "tau = eps0 eps_r(1/(2 pi tau))/K auf der gemessenen Dispersionskurve, mit "
       "verschwindendem Selbstkonsistenzresiduum.",
       "Ein EINZELNER eps_r-Wert bleibt MissingMaterialData -- keine Quelle nennt Reinheit "
       "und Wassergehalt.  Belegt ist ein Band (tau = 4,2e-11 .. 1,3e-10 s); eine "
       "Validierung gegen eine Messung von tau selbst gibt es nicht.",
       true},
  };
}

}  // namespace es
