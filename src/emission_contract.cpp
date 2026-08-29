#include "es/emission_contract.hpp"

#include <cmath>
#include <limits>

#include "es/constants.hpp"

namespace es {

using constants::amu;
using constants::e;
using constants::eps0;
using constants::h_planck;
using constants::kB;
using constants::pi;

namespace {
constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();
}

const char* to_string(Polarity p) {
  return p == Polarity::Positive ? "positive" : "negative";
}

const char* to_string(EmissionStatus s) {
  switch (s) {
    case EmissionStatus::Disabled: return "Disabled";
    case EmissionStatus::MissingEmissionParameters: return "MissingEmissionParameters";
    case EmissionStatus::EquationNotValidated: return "EquationNotValidated";
    case EmissionStatus::PolarityMismatch: return "PolarityMismatch";
    case EmissionStatus::ParameterOutOfRange: return "ParameterOutOfRange";
    case EmissionStatus::Ok: return "Ok";
  }
  return "?";
}

const char* explain(EmissionStatus s) {
  switch (s) {
    case EmissionStatus::Disabled:
      return "Das Emissionsmodell ist abgeschaltet.  Das ist der Auslieferungszustand und "
             "keine Panne: eine Emission darf erst rechnen, wenn Gleichung und Parameter "
             "belegt sind.";
    case EmissionStatus::MissingEmissionParameters:
      return "Mindestens eine Pflichtangabe fehlt -- Ladungszahl, Masse oder "
             "Aktivierungsbarriere mit Fundstelle.  Es wird KEIN Ersatzwert gesetzt und "
             "nicht auf die quellenlosen Diagnosewerte des alten Pfades zurueckgegriffen.  "
             "Fuer EMI-BF4 fehlt die Barriere: es gibt keine belegte Zahl, und j haengt "
             "exponentiell von ihr ab.";
    case EmissionStatus::EquationNotValidated:
      return "Die Ratengleichung selbst wurde in diesem Projekt nicht an einer Primaerquelle "
             "geprueft.  Vorfaktor, genaue Definition der Flaechenladungsdichte und "
             "Gueltigkeitsbereich sind damit unbelegt; eine vollstaendige Parameterliste "
             "reicht dann nicht.";
    case EmissionStatus::PolarityMismatch:
      return "Das Vorzeichen der Spezies passt nicht zur angeforderten Polaritaet, oder das "
             "Feld zeigt in die Richtung, in der diese Spezies nicht emittiert wuerde.  Die "
             "beiden Polaritaeten werden nie aus denselben Speziesdaten gerechnet.";
    case EmissionStatus::ParameterOutOfRange:
      return "Ein Parameter liegt ausserhalb des Bereichs, fuer den die Gleichung angegeben "
             "ist.";
    case EmissionStatus::Ok:
      return "Alle Pruefungen sind bestanden.";
  }
  return "?";
}

// ---------------------------------------------------------------------------

EmissionStatus EmissionModel::status() const {
  if (!enabled) return EmissionStatus::Disabled;
  if (n_species == 0 || species == nullptr) return EmissionStatus::MissingEmissionParameters;
  if (!(temperature > 0.0)) return EmissionStatus::MissingEmissionParameters;
  bool any_matching = false;
  for (std::size_t k = 0; k < n_species; ++k) {
    const EmittedSpecies& s = species[k];
    if (s.polarity() != polarity) continue;
    any_matching = true;
    if (!s.complete()) return EmissionStatus::MissingEmissionParameters;
  }
  if (!any_matching) return EmissionStatus::PolarityMismatch;
  if (!equation_validated) return EmissionStatus::EquationNotValidated;
  return EmissionStatus::Ok;
}

void EmissionModel::print(std::FILE* out) const {
  std::fprintf(out, "  Emissionsmodell (%s): %s\n", to_string(polarity), to_string(status()));
  std::fprintf(out, "    aktiviert %s, Gleichung geprueft %s, T = %.4g K, %zu Spezies\n",
               enabled ? "ja" : "NEIN", equation_validated ? "ja" : "NEIN", temperature,
               n_species);
  for (std::size_t k = 0; k < n_species; ++k) {
    const EmittedSpecies& s = species[k];
    std::fprintf(out, "      %-10s z = %+d, m = %.6g kg [%s]\n", s.name, s.charge_number,
                 s.mass, s.mass_source[0] ? "belegt" : "OHNE QUELLE");
    if (s.has_barrier())
      std::fprintf(out, "                 dG = %.6g J = %.4f eV [%s]\n", s.activation_barrier,
                   s.activation_barrier / e, s.barrier_source);
    else
      std::fprintf(out, "                 dG FEHLT -- keine belegte Zahl\n");
    if (s.note[0]) std::fprintf(out, "                 %s\n", s.note);
  }
  std::fprintf(out, "    %s\n", explain(status()));
}

// ---------------------------------------------------------------------------
// The kernel
// ---------------------------------------------------------------------------

Real schottky_barrier_lowering(Real E) {
  if (!(E >= 0.0)) return kNaN;
  return std::sqrt(e * e * e * E / (4.0 * pi * eps0));
}

Real barrier_free_field(Real dG) {
  if (!(dG > 0.0)) return kNaN;
  return 4.0 * pi * eps0 * dG * dG / (e * e * e);
}

Real iribarne_thomson_rate(Real E, Real dG, Real T) {
  if (!(T > 0.0) || !(E >= 0.0) || !(dG > 0.0)) return kNaN;
  const Real kT = kB * T;
  const Real sigma_s = eps0 * E;
  const Real barrier = dG - schottky_barrier_lowering(E);
  return (kT / h_planck) * sigma_s * std::exp(-barrier / kT);
}

Real iribarne_thomson_dimensionless(Real x, Real b) {
  if (!(x >= 0.0) || !(b > 0.0)) return kNaN;
  // j / j0 = (E/E*) exp[-b (1 - sqrt(E/E*))], with the sigma_s = eps0 E factor
  // carried along so that the form is the same function, only scaled.
  return x * std::exp(-b * (1.0 - std::sqrt(x)));
}

// ---------------------------------------------------------------------------

EmissionResult emitted_current_density(const EmissionModel& m, Real E_n) {
  EmissionResult r;
  r.status = m.status();
  r.message = explain(r.status);
  if (!is_usable(r.status)) return r;   // no number, ever, on this path

  // Field direction.  A cation leaves where the outward field is positive; an
  // anion where it is negative.  The magnitude alone would emit both, which is
  // exactly the defect the prototype had.
  const bool right_way = (m.polarity == Polarity::Positive) ? (E_n > 0.0) : (E_n < 0.0);
  if (!right_way) {
    r.status = EmissionStatus::PolarityMismatch;
    r.message = explain(r.status);
    return r;
  }
  const Real E = std::abs(E_n);
  Real total = 0.0;
  for (std::size_t k = 0; k < m.n_species; ++k) {
    const EmittedSpecies& s = m.species[k];
    if (s.polarity() != m.polarity) continue;
    r.barrier_lowering = schottky_barrier_lowering(E);
    r.effective_barrier = s.activation_barrier - r.barrier_lowering;
    total += iribarne_thomson_rate(E, s.activation_barrier, m.temperature);
  }
  r.current_density = total;
  return r;
}

// ---------------------------------------------------------------------------
// The shipped data
// ---------------------------------------------------------------------------

namespace {

// Molar masses of the two ions of EMI-BF4.  The compound's molar mass,
// 197.97 g/mol, is stated by the IoLiTec data sheet IL-0006 (page 1/4) and by
// the ILThermo component record; the split into 111.17 and 86.81 g/mol is the
// formula composition C6H11N2+ and BF4-.  Only the MASSES are sourced this way;
// the barriers are not, and are left absent.
constexpr Real kMassEmi = 111.17e-3 / 6.02214076e23;
constexpr Real kMassBf4 = 86.81e-3 / 6.02214076e23;

const EmittedSpecies kEmiBf4Positive[] = {
    {"EMI+", +1, kMassEmi,
     "Formelmasse C6H11N2+; Verbindungsmolmasse 197,97 g/mol aus IoLiTec IL-0006 Seite 1/4 "
     "und dem ILThermo-Komponentendatensatz",
     0.0, "", 0.0,
     "Die Aktivierungsbarriere FEHLT.  Es gibt keine belegte Zahl fuer EMI-BF4, und j haengt "
     "exponentiell von ihr ab: die in der Literatur genannte Spanne 1,0 bis 1,4 eV ist bei "
     "298 K ein Faktor 1e7 im Strom."},
};

const EmittedSpecies kEmiBf4Negative[] = {
    {"BF4-", -1, kMassBf4,
     "Formelmasse BF4-; Verbindungsmolmasse 197,97 g/mol aus IoLiTec IL-0006 Seite 1/4 "
     "und dem ILThermo-Komponentendatensatz",
     0.0, "", 0.0,
     "Die Aktivierungsbarriere FEHLT, und sie ist NICHT die des Kations: Masse, "
     "Solvatationsenergie und die Zusammensetzung der emittierten Cluster unterscheiden sich "
     "zwischen den Polaritaeten."},
};

const EmissionModel kBlockedPositive = {
    false, false, 0.0, Polarity::Positive, kEmiBf4Positive, 1,
};
const EmissionModel kBlockedNegative = {
    false, false, 0.0, Polarity::Negative, kEmiBf4Negative, 1,
};

}  // namespace

const EmissionModel& emibf4_emission_blocked(Polarity p) {
  return (p == Polarity::Positive) ? kBlockedPositive : kBlockedNegative;
}

SyntheticEmissionModel synthetic_complete_model(Polarity p, Real dG, Real T) {
  SyntheticEmissionModel out;
  out.species = EmittedSpecies{"synthetisch",
                               (p == Polarity::Positive) ? +1 : -1,
                               1.0e-25,
                               "keine -- ausdruecklich synthetisch",
                               dG,
                               "keine -- ausdruecklich synthetisch",
                               0.0,
                               "Kein Stoff.  Diese Spezies existiert nur, damit ein Test den "
                               "Zweig EquationNotValidated erreichen kann."};
  // The source strings are deliberately non-empty so that has_mass() and
  // has_barrier() pass: the point of this model is to show that COMPLETE
  // parameters are still not enough while the equation itself is unverified.
  out.model.enabled = true;
  out.model.equation_validated = false;
  out.model.temperature = T;
  out.model.polarity = p;
  out.model.species = &out.species;
  out.model.n_species = 1;
  return out;
}

}  // namespace es
