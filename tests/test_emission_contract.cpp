// tests/test_emission_contract.cpp -- P5: the ion-emission contract
//
// The model is BLOCKED, so the tests fall into two groups that must not be
// confused:
//
//   A. THE CONTRACT.  Every path that could produce a number must fail closed,
//      and it must fail for the RIGHT reason.  Missing barrier, unvalidated
//      equation, wrong polarity, wrong field direction, model switched off --
//      each has its own status and none of them returns a current.
//
//   B. THE MATHEMATICAL KERNEL.  The quoted functional form is a pure function
//      of its arguments and can be tested without any material data:
//      dimensions, limits, monotonicity in E and in dG, and the sqrt(E) law of
//      the barrier lowering.  Testing the FORM is not validating the MODEL, and
//      the two are kept apart on purpose.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/emission_contract.hpp"

using namespace es;
using constants::e;
using constants::eps0;
using constants::h_planck;
using constants::kB;
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

}  // namespace

int main() {
  std::printf("P5 -- Emissionsvertrag: der Zustand ist BLOCKED\n\n");

  // =========================================================================
  std::printf("A. Der Vertrag schlaegt geschlossen fehl -- und aus dem richtigen Grund\n");
  {
    for (Polarity p : {Polarity::Positive, Polarity::Negative}) {
      const EmissionModel& m = emibf4_emission_blocked(p);
      m.print(stdout);
      check(m.status() == EmissionStatus::Disabled,
            std::string("EMI-BF4 ") + to_string(p) + ": ausgeliefert wird Disabled");
      const EmissionResult r = emitted_current_density(m, 1.0e9);
      check(!r.usable() && r.current_density == 0.0,
            "und es kommt kein Strom heraus");
    }

    // Switching it on is NOT enough.
    EmissionModel on = emibf4_emission_blocked(Polarity::Positive);
    on.enabled = true;
    on.temperature = 298.15;
    std::printf("    nach dem Einschalten: %s\n", to_string(on.status()));
    check(on.status() == EmissionStatus::MissingEmissionParameters,
          "eingeschaltet meldet es MissingEmissionParameters -- die Barriere fehlt");
    check(emitted_current_density(on, 1.0e9).current_density == 0.0,
          "und liefert weiterhin keinen Strom");

    // Complete parameters are STILL not enough while the equation is unverified.
    const SyntheticEmissionModel syn_owner =
        synthetic_complete_model(Polarity::Positive, 1.09 * e, 298.15);
    const EmissionModel& syn = *syn_owner;
    std::printf("    vollstaendige Parameter, ungepruefte Gleichung: %s\n",
                to_string(syn.status()));
    check(syn.status() == EmissionStatus::EquationNotValidated,
          "vollstaendige Parameter reichen nicht: EquationNotValidated");
    check(emitted_current_density(syn, 1.0e9).current_density == 0.0,
          "und auch dieser Pfad liefert keinen Strom");

    // Polarity is never taken from a magnitude.
    SyntheticEmissionModel wrong_owner =
        synthetic_complete_model(Polarity::Positive, 1.09 * e, 298.15);
    wrong_owner.model.polarity = Polarity::Negative;
    const EmissionModel& wrong = *wrong_owner;
    check(wrong.status() == EmissionStatus::PolarityMismatch,
          "eine Kationenspezies mit negativer Polaritaet ist ein PolarityMismatch");

    // Even a fully validated model must respect the field DIRECTION.
    SyntheticEmissionModel val_owner =
        synthetic_complete_model(Polarity::Positive, 1.09 * e, 298.15);
    val_owner.model.equation_validated = true;
    const EmissionModel& val = *val_owner;
    check(val.status() == EmissionStatus::Ok, "ein synthetisch freigegebenes Modell ist Ok");
    check(emitted_current_density(val, +1.0e9).usable(),
          "und emittiert Kationen bei nach aussen zeigendem Feld");
    const EmissionResult back = emitted_current_density(val, -1.0e9);
    check(back.status == EmissionStatus::PolarityMismatch && back.current_density == 0.0,
          "bei umgekehrtem Feld emittiert es NICHT -- nicht denselben Betrag");

    // And the two polarities are not the same computation.
    SyntheticEmissionModel neg_owner =
        synthetic_complete_model(Polarity::Negative, 1.09 * e, 298.15);
    neg_owner.model.equation_validated = true;
    const EmissionModel& neg = *neg_owner;
    const EmissionResult pos_r = emitted_current_density(val, +1.0e9);
    const EmissionResult neg_r = emitted_current_density(neg, -1.0e9);
    check(pos_r.usable() && neg_r.usable(),
          "beide Polaritaeten rechnen, jede mit IHREN eigenen Speziesdaten");
    check(pos_r.current_density > 0.0 && neg_r.current_density > 0.0,
          "und beide liefern einen Strom -- ein geteilter Puffer haette einen davon auf "
          "null gesetzt, und genau so wurde der Aliasingfehler gefunden");
    std::printf("    zum Vergleich: j+ = %.4e, j- = %.4e A/m^2 (beide synthetisch)\n",
                pos_r.current_density, neg_r.current_density);
  }

  // =========================================================================
  std::printf("\nB. Der mathematische Kern -- Dimensionen, Grenzen, Monotonie\n");
  {
    // Barrier lowering: zero at zero field, sqrt(E), and the closed form.
    check(schottky_barrier_lowering(0.0) == 0.0,
          "die Barrierensenkung ist bei E = 0 exakt null");
    const Real E1 = 1.0e9;
    check_rel(schottky_barrier_lowering(E1), std::sqrt(e * e * e * E1 / (4.0 * pi * eps0)),
              1.0e-15, "G(E) = sqrt(e^3 E / (4 pi eps0))");
    check_rel(schottky_barrier_lowering(4.0 * E1), 2.0 * schottky_barrier_lowering(E1),
              1.0e-15, "sie waechst wie sqrt(E): vierfaches Feld, doppelte Senkung");

    // The barrier-free field, and that it is what it says it is.
    const Real dG = 1.09 * e;
    const Real Estar = barrier_free_field(dG);
    std::printf("    dG = %.4f eV  ->  E* = %.4e V/m\n", dG / e, Estar);
    check_rel(schottky_barrier_lowering(Estar), dG, 1.0e-12,
              "bei E* ist die Senkung genau die Barriere");

    // Dimensions.  j = (kT/h) eps0 E exp(...) must come out in A/m^2, checked
    // by rebuilding the prefactor from its parts.
    const Real T = 298.15;
    const Real j = iribarne_thomson_rate(Estar, dG, T);
    const Real prefactor = (kB * T / h_planck) * eps0 * Estar;
    check_rel(j, prefactor, 1.0e-12,
              "bei E = E* ist die Rate genau der Vorfaktor (kT/h) eps0 E");
    std::printf("    Vorfaktor bei E*: %.4e A/m^2\n", prefactor);

    // Monotonicity in E, over ten decades of the current.
    {
      bool monotone = true;
      Real last = -1.0;
      for (int k = 0; k <= 60; ++k) {
        const Real E = 0.2 * Estar * std::pow(10.0, 0.02 * k);
        const Real v = iribarne_thomson_rate(E, dG, T);
        if (!(v > last)) monotone = false;
        last = v;
      }
      check(monotone, "die Rate waechst streng monoton mit dem Feld");
    }
    // Monotonicity in dG, the other way round, and it is EXPONENTIAL.
    {
      const Real E = 0.8 * Estar;
      const Real j1 = iribarne_thomson_rate(E, 1.0 * e, T);
      const Real j2 = iribarne_thomson_rate(E, 1.4 * e, T);
      check(j2 < j1, "eine hoehere Barriere gibt weniger Strom");
      std::printf("    dG = 1,0 eV -> %.3e,  dG = 1,4 eV -> %.3e A/m^2, Verhaeltnis %.3e\n",
                  j1, j2, j1 / j2);
      // THE POINT OF THE BLOCKER: the literature span of dG is a factor 1e7.
      check(j1 / j2 > 1.0e6,
            "die in der Literatur genannte Spanne 1,0 bis 1,4 eV ist ueber sechs "
            "Groessenordnungen im Strom -- deshalb ist eine unbelegte Barriere keine "
            "Vorhersage");
    }
    // Temperature enters twice, in the prefactor and in the exponent.
    {
      const Real E = 0.8 * Estar;
      check(iribarne_thomson_rate(E, dG, 350.0) > iribarne_thomson_rate(E, dG, 250.0),
            "hoehere Temperatur gibt mehr Strom");
    }
    // Fail-closed of the kernel itself.
    check(!std::isfinite(iribarne_thomson_rate(E1, dG, 0.0)), "T = 0 gibt nan, nicht null");
    check(!std::isfinite(iribarne_thomson_rate(E1, 0.0, T)), "dG = 0 gibt nan");
    check(iribarne_thomson_rate(0.0, dG, T) == 0.0,
          "bei E = 0 ist die Rate exakt null -- sigma_s = eps0 E ist dort null");

    // The dimensionless form is the same function, only scaled.
    {
      const Real b = dG / (kB * T);
      const Real x = 0.7;
      const Real jd = iribarne_thomson_dimensionless(x, b);
      const Real j0 = (kB * T / h_planck) * eps0 * Estar;
      check_rel(jd * j0, iribarne_thomson_rate(x * Estar, dG, T), 1.0e-12,
                "die dimensionslose Form ist dieselbe Funktion, nur skaliert");
      check_rel(iribarne_thomson_dimensionless(1.0, b), 1.0, 1.0e-15,
                "und ist bei x = 1 genau eins");
    }
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
