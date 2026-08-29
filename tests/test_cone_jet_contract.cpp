// tests/test_cone_jet_contract.cpp -- P8: the cone-jet contract
//
// The mode is BLOCKED, so the tests are of two kinds:
//   A. the refusal, and that it names both of its reasons;
//   B. the dimensionless numbers, which are DEFINITIONS and can therefore be
//      checked against their own derivations -- dimensions, limits, and the
//      balances they express -- without any literature value.
//
// Checking a definition is not validating a correlation.  No number here is a
// regime prediction.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/cone_jet_contract.hpp"
#include "es/material_data.hpp"
#include "es/status.hpp"

using namespace es;
using constants::eps0;
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
  std::printf("P8 -- Cone-Jet: eigener Vertrag, eigener Status (BLOCKED)\n\n");

  // =========================================================================
  std::printf("A. Die Verweigerung nennt BEIDE Gruende\n");
  {
    bool threw = false;
    std::string msg;
    try {
      cone_jet_current();
    } catch (const NotImplementedInThisPhase& x) {
      threw = true;
      msg = x.what();
    }
    check(threw, "cone_jet_current() wirft NotImplementedInThisPhase");
    check(msg.find("Zweiphasenstroemung") != std::string::npos,
          "und nennt die fehlende Zweiphasenstroemung");
    check(msg.find("Erratum") != std::string::npos,
          "und dass auch die empirische Skalierung nicht uebernommen wird -- samt Erratum");
    check(msg.find("Suchtreffer") != std::string::npos,
          "und dass Schnipselzahlen nicht als Quelle gelten");
    check(msg.find("ConeJetModel") != std::string::npos,
          "und dass der alte Block unangetastet und unbenutzt bleibt");

    const std::vector<ConeJetRequirement> req = cone_jet_requirements();
    std::printf("    %zu Teilmodelle im Vertrag:\n", req.size());
    int missing = 0;
    for (const ConeJetRequirement& r : req) {
      std::printf("      [%s] %-50s %s\n", r.available ? "da " : "FEHLT", r.name,
                  r.provided_by);
      if (!r.available) ++missing;
    }
    check(missing >= 4, "mindestens vier Teilmodelle fehlen -- daher blockiert");
    check(req.size() >= 6, "der Vertrag listet die noetigen Teilmodelle vollstaendig auf");
  }

  // =========================================================================
  std::printf("\nB. Die dimensionslosen Zahlen sind Definitionen\n");
  {
    const Real gamma = 0.05401, rho = 1280.9, mu = 0.03637, K = 1.5584, eps_r = 12.8;
    const Real a = 5.0e-6;

    // Each one against its own derivation, written out a second time.
    check_rel(charge_relaxation_time_cj(eps_r, K), eps0 * eps_r / K, 1.0e-15,
              "tau_e = eps0 eps_r / K");
    check_rel(capillary_inertial_time(rho, gamma, a), std::sqrt(rho * a * a * a / gamma),
              1.0e-15, "t_c = sqrt(rho a^3 / gamma)");
    check_rel(visco_capillary_time(mu, gamma, a), mu * a / gamma, 1.0e-15,
              "t_v = mu a / gamma");
    check_rel(ohnesorge(mu, rho, gamma, a),
              visco_capillary_time(mu, gamma, a) / capillary_inertial_time(rho, gamma, a),
              1.0e-14, "Oh ist das Verhaeltnis der beiden Zeiten -- das ist seine Bedeutung");

    // THE ELECTROHYDRODYNAMIC LENGTH, checked against the balance it is defined
    // by: at r = r*, the charge relaxation time equals the capillary-inertial
    // time.  That is the derivation, and it is what the test verifies.
    const Real rstar = electrohydrodynamic_length(gamma, eps_r, rho, K);
    std::printf("    r* = %.6e m\n", rstar);
    check_rel(capillary_inertial_time(rho, gamma, rstar),
              charge_relaxation_time_cj(eps_r, K), 1.0e-12,
              "bei r* sind Ladungsrelaxations- und Kapillarzeit gleich -- die Definition");
    check(rstar < a, "r* liegt weit unter dem Bohrungsradius");

    // The electric Bond number IS the ratio P3b balances.
    const Real E = 5.0e7;
    check_rel(electric_bond(E, a, gamma), (0.5 * eps0 * E * E) / (gamma / a), 1.0e-15,
              "Bo_E = (eps0 E^2/2) / (gamma/a) -- genau die beiden Terme von P3b");

    // Feed numbers against their definitions.
    const Real Q = 1.0e-13, R = 5.0e-6;
    check_rel(feed_reynolds(rho, Q, R, mu), rho * (Q / (pi * R * R)) * (2.0 * R) / mu, 1.0e-14,
              "Re = rho u d / mu mit u = Q/(pi R^2)");
    check_rel(feed_capillary(mu, Q, R, gamma), mu * (Q / (pi * R * R)) / gamma, 1.0e-14,
              "Ca = mu u / gamma");

    // Fail closed on every missing input.
    check(!std::isfinite(charge_relaxation_time_cj(0.0, K)), "ohne eps_r kein tau_e");
    check(!std::isfinite(electrohydrodynamic_length(gamma, 0.0, rho, K)), "ohne eps_r kein r*");
    check(!std::isfinite(ohnesorge(0.0, rho, gamma, a)), "ohne mu kein Oh");
    check(!std::isfinite(feed_reynolds(rho, Q, 0.0, mu)), "ohne Radius kein Re");
  }

  // =========================================================================
  std::printf("\nC. Die Diagnose mit den belegten Stoffdaten\n");
  {
    const ConeJetDiagnosis c =
        diagnose_cone_jet(emibf4_sourced(), 298.15, 5.0e-6, 1.0e-13, 5.0e-6, 5.0e7);
    c.print(stdout);
    check(c.status == ConeJetStatus::MissingMaterialData,
          "mit den belegten Daten fehlt eps_r, also MissingMaterialData");
    check(!std::isfinite(c.tau_e) && !std::isfinite(c.r_star),
          "tau_e und r* bleiben nan -- kein Ersatzwert");
    check(std::isfinite(c.Oh) && std::isfinite(c.Bo_E) && std::isfinite(c.Re),
          "die Kennzahlen, die eps_r NICHT brauchen, sind trotzdem auswertbar");
    check(c.K_status != MaterialDataStatus::MissingMaterialData,
          "die Leitfaehigkeit ist belegt -- es fehlt genau eine Groesse");
    std::printf("    Oh = %.4e, Bo_E = %.4e, Re = %.4e, Ca = %.4e\n", c.Oh, c.Bo_E, c.Re,
                c.Ca);
    check(c.Oh > 1.0,
          "Oh > 1: die Viskositaet dominiert die Traegheit auf dieser Laengenskala -- eine "
          "DIAGNOSE, keine Regimevorhersage");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
