// tests/test_feed.cpp -- P1: the pressure budget at the exit plane
//
// The budget has exactly one non-trivial term, the laminar pressure drop of a
// straight filled circular channel.  It is tested against the closed forms it
// claims to implement, and against an INDEPENDENT numerical integration of the
// parabolic profile -- so the resistance is checked against the velocity field
// it is derived from, not against a second copy of the same formula.
//
// The signs are tested one term at a time, because a budget with the wrong sign
// on one term still looks like a budget.

#include <cmath>
#include <cstdio>
#include <string>

#include "es/constants.hpp"
#include "es/feed.hpp"
#include "es/liquid.hpp"

using namespace es;
using constants::pi;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  [%s] %s\n", ok ? "ok" : "FEHLER", what.c_str());
  if (!ok) ++failures;
}

void check_rel(Real got, Real want, Real tol, const std::string& what) {
  const Real e =
      (std::abs(want) > 0.0) ? std::abs(got - want) / std::abs(want) : std::abs(got - want);
  std::printf("  [%s] %s: %.9e gegen %.9e, rel. %.3e (Grenze %.1e)\n",
              (e <= tol) ? "ok" : "FEHLER", what.c_str(), got, want, e, tol);
  if (!(e <= tol)) ++failures;
}

LiquidProperties test_liquid() {
  LiquidProperties L = emibf4_illustrative();
  return L;
}

}  // namespace

// ===========================================================================

int main() {
  std::printf("P1 -- Druckhaushalt am Austritt\n\n");
  const LiquidProperties liquid = test_liquid();
  const Real mu = liquid.documented_only.mu;
  const Real rho = liquid.rho;
  const Real a = 5.0e-6;

  // =========================================================================
  // 1.  Hagen-Poiseuille against its own velocity profile
  // =========================================================================
  std::printf("1. Hagen-Poiseuille gegen das Geschwindigkeitsprofil\n");
  {
    FeedChannel ch;
    ch.radius = 5.0e-6;
    ch.length = 300.0e-6;
    const Real Q = 0.5e-12;  // 0.5 pL/s

    const Real Rh = ch.hydraulic_resistance(mu);
    const Real dp = Rh * Q;
    std::printf("    R_h = %.6e Pa s/m^3, dp = %.6e Pa bei Q = %.3e m^3/s\n", Rh, dp, Q);

    // The closed form, written out here a second time from the definition
    // rather than called: 8 mu L / (pi R^4).
    check_rel(Rh, 8.0 * mu * ch.length / (pi * std::pow(ch.radius, 4.0)), 1.0e-15,
              "R_h = 8 mu L / (pi R^4)");

    // INDEPENDENT PATH.  Integrate the parabolic profile the resistance is
    // derived from, u(r) = (dp/(4 mu L)) (R^2 - r^2), over the cross section by
    // Simpson, and check that it returns Q.  This closes the loop between the
    // profile and the resistance without reusing either formula.
    {
      const int n = 20000;  // even
      const Real h = ch.radius / static_cast<Real>(n);
      auto u = [&](Real r) { return dp / (4.0 * mu * ch.length) * (ch.radius * ch.radius - r * r); };
      Real acc = 0.0;
      for (int k = 0; k <= n; ++k) {
        const Real r = static_cast<Real>(k) * h;
        const Real w = (k == 0 || k == n) ? 1.0 : ((k % 2) ? 4.0 : 2.0);
        acc += w * u(r) * 2.0 * pi * r;
      }
      acc *= h / 3.0;
      check_rel(acc, Q, 1.0e-9, "Integral des Profils ueber den Querschnitt ergibt Q");
    }

    check_rel(ch.mean_velocity(Q), Q / (pi * ch.radius * ch.radius), 1.0e-15,
              "mittlere Geschwindigkeit");
    check_rel(ch.centreline_velocity(Q), 2.0 * ch.mean_velocity(Q), 1.0e-15,
              "Mittengeschwindigkeit ist das Doppelte der mittleren");
    check_rel(ch.velocity_at(0.0, Q), ch.centreline_velocity(Q), 1.0e-15,
              "Profil auf der Achse");
    check(std::abs(ch.velocity_at(ch.radius, Q)) < 1.0e-30, "Haftbedingung an der Wand");

    // Wall shear stress from two independent forms: 4 mu u_mean / R and
    // dp R / (2 L), which are equal only if the resistance is right.
    check_rel(ch.wall_shear_stress(mu, Q), dp * ch.radius / (2.0 * ch.length), 1.0e-14,
              "Wandschubspannung: 4 mu u/R gegen dp R / (2 L)");

    // Scaling: R -> 2R must divide the resistance by 16, L -> 2L must double it.
    FeedChannel wide = ch;
    wide.radius *= 2.0;
    check_rel(wide.hydraulic_resistance(mu), Rh / 16.0, 1.0e-14,
              "doppelter Radius: ein Sechzehntel des Widerstands");
    FeedChannel longer = ch;
    longer.length *= 2.0;
    check_rel(longer.hydraulic_resistance(mu), 2.0 * Rh, 1.0e-14,
              "doppelte Laenge: doppelter Widerstand");
  }

  // =========================================================================
  // 2.  Hydrostatics, sign by sign
  // =========================================================================
  std::printf("\n2. Hydrostatik: Vorzeichen und Referenzhoehe\n");
  {
    FeedRequest q;
    q.mode = PressureMode::Budget;
    q.channel.radius = 5.0e-6;
    q.channel.length = 300.0e-6;
    q.p_reservoir = 1000.0;
    q.p_vacuum = 0.0;
    q.Q = 0.0;                    // no flow: the viscous term must vanish exactly
    q.gravity_axial = -constants::g0;   // emitter pointing up, on the bench
    q.z_exit = 0.0;
    q.z_reservoir = -0.010;       // reservoir 10 mm below the exit plane

    const PressureBudget b = solve_pressure_budget(q, liquid, a);
    b.print(stdout);
    check(is_usable(b.status), "der Haushalt ist auswertbar");
    check(b.viscous == 0.0, "ohne Stroemung ist der viskose Term exakt null");
    check_rel(b.hydrostatic, rho * constants::g0 * 0.010, 1.0e-14,
              "hydrostatischer Verlust rho g H beim Anheben um 10 mm");
    check_rel(b.delta_p_exit, 1000.0 - rho * constants::g0 * 0.010, 1.0e-14,
              "delta_p_exit = (p_res - p_vak) - rho g H");
    check(b.hydrostatic > 0.0,
          "Anheben kostet Druck: der hydrostatische Term ist ein Verlust, kein Gewinn");

    // Reservoir ABOVE the exit plane: the term must change sign.
    FeedRequest q2 = q;
    q2.z_reservoir = +0.010;
    const PressureBudget b2 = solve_pressure_budget(q2, liquid, a);
    check_rel(b2.hydrostatic, -rho * constants::g0 * 0.010, 1.0e-14,
              "Vorrat oberhalb des Austritts: der Term kehrt das Vorzeichen um");

    // No gravity -- the case the device is actually built for.
    FeedRequest q3 = q;
    q3.gravity_axial = 0.0;
    const PressureBudget b3 = solve_pressure_budget(q3, liquid, a);
    check(b3.hydrostatic == 0.0, "ohne Schwerkraft ist der hydrostatische Term exakt null");
    check_rel(b3.delta_p_exit, 1000.0, 1.0e-15, "und delta_p_exit ist der reine Antrieb");
  }

  // =========================================================================
  // 3.  The viscous term in the budget, and the flow direction
  // =========================================================================
  std::printf("\n3. Viskoser Term und Stroemungsrichtung\n");
  {
    FeedRequest q;
    q.mode = PressureMode::Budget;
    q.channel.radius = 5.0e-6;
    q.channel.length = 300.0e-6;
    q.p_reservoir = 5000.0;
    q.p_vacuum = 0.0;
    q.gravity_axial = 0.0;
    q.Q = 1.0e-15;
    const PressureBudget b = solve_pressure_budget(q, liquid, a);
    b.print(stdout);
    check(is_usable(b.status), "der Haushalt ist auswertbar");
    check_rel(b.viscous, q.channel.hydraulic_resistance(mu) * q.Q, 1.0e-15,
              "viskoser Term ist R_h Q");
    check(b.delta_p_exit < b.driving,
          "Foerdern zum Austritt kostet Druck: delta_p_exit liegt unter dem Antrieb");

    FeedRequest back = q;
    back.Q = -q.Q;
    const PressureBudget bb = solve_pressure_budget(back, liquid, a);
    check_rel(bb.viscous, -b.viscous, 1.0e-15,
              "umgekehrte Stroemungsrichtung kehrt den viskosen Term um");
    check(bb.delta_p_exit > b.driving, "und Rueckstroemen erhoeht den Austrittsdruck");

    // Linearity of the whole budget in Q -- it is an affine function and must
    // be exactly that.
    FeedRequest twice = q;
    twice.Q = 2.0 * q.Q;
    const PressureBudget bt = solve_pressure_budget(twice, liquid, a);
    check_rel(b.driving - bt.delta_p_exit, 2.0 * (b.driving - b.delta_p_exit), 1.0e-14,
              "der Haushalt ist exakt affin in Q");
  }

  // =========================================================================
  // 4.  Validity bounds, and that they fail closed
  // =========================================================================
  std::printf("\n4. Gueltigkeitsgrenzen\n");
  {
    FeedRequest q;
    q.mode = PressureMode::Budget;
    q.channel.radius = 5.0e-6;
    q.channel.length = 300.0e-6;
    q.p_reservoir = 0.0;
    q.gravity_axial = 0.0;

    // Push the Reynolds number over the laminar bound.  With a 5 um channel and
    // this viscosity that needs an absurd flow rate -- which is the point: the
    // check exists, and it is reachable.
    const Real Q_turb = feed::kReynoldsLaminar * mu * pi * q.channel.radius /
                        (2.0 * rho) * 1.2;
    q.Q = Q_turb;
    const PressureBudget bt = solve_pressure_budget(q, liquid, a);
    std::printf("    Re = %.4g bei Q = %.3e m^3/s -> %s\n", bt.reynolds, q.Q,
                to_string(bt.status));
    check(bt.status == FeedStatus::NotLaminar,
          "ueber der laminaren Schranke gibt es NotLaminar statt eines Ergebnisses");
    check(!is_usable(bt.status), "und der Status ist nicht brauchbar");

    // A channel with no geometry.
    FeedRequest bad;
    bad.mode = PressureMode::Budget;
    const PressureBudget bb = solve_pressure_budget(bad, liquid, a);
    check(bb.status == FeedStatus::ChannelGeometryInvalid,
          "ein Kanal ohne Radius oder Laenge schlaegt geschlossen fehl");
    check(!std::isfinite(bb.delta_p_exit),
          "und liefert nan statt einer null -- null waere ein physikalischer Wert");

    // A liquid without a viscosity, asked for a flow.
    LiquidProperties dry = liquid;
    dry.documented_only.mu = 0.0;
    FeedRequest flow;
    flow.mode = PressureMode::Budget;
    flow.channel.radius = 5.0e-6;
    flow.channel.length = 300.0e-6;
    flow.Q = 1.0e-15;
    const PressureBudget bd = solve_pressure_budget(flow, dry, a);
    check(bd.status == FeedStatus::MissingLiquidProperty,
          "ohne Viskositaet gibt es keinen viskosen Term, sondern MissingLiquidProperty");

    // A contact angle must be refused, not ignored.
    FeedRequest ang;
    ang.mode = PressureMode::Budget;
    ang.channel.radius = 5.0e-6;
    ang.channel.length = 300.0e-6;
    ang.contact_angle_requested = true;
    const PressureBudget ba = solve_pressure_budget(ang, liquid, a);
    check(ba.status == FeedStatus::MissingFeedInput,
          "ein Kontaktwinkel im gefuellten Kanal wird abgelehnt, nicht stillschweigend "
          "uebergangen");
  }

  // =========================================================================
  // 5.  The direct mode stays exactly what it was
  // =========================================================================
  std::printf("\n5. Direkter Modus\n");
  {
    FeedRequest q;
    q.mode = PressureMode::Direct;
    q.delta_p_exit_direct = -1234.5;
    const PressureBudget b = solve_pressure_budget(q, liquid, a);
    check(is_usable(b.status), "der direkte Modus ist brauchbar");
    check(b.delta_p_exit == -1234.5, "und liefert die Eingabe bitgenau zurueck");
    check(!std::isfinite(b.driving) && !std::isfinite(b.hydrostatic) &&
              !std::isfinite(b.viscous),
          "die Einzelterme sind dort nicht bekannt und stehen als nan, nicht als null");
  }

  // =========================================================================
  // 6.  Where the budget sits on the capillary scale
  // =========================================================================
  std::printf("\n6. Einordnung gegen gamma/a\n");
  {
    FeedRequest q;
    q.mode = PressureMode::Budget;
    q.channel.radius = 5.0e-6;
    q.channel.length = 300.0e-6;
    q.gravity_axial = 0.0;
    q.Q = 0.0;
    q.p_reservoir = liquid.gamma / a;   // exactly one capillary pressure scale
    const PressureBudget b = solve_pressure_budget(q, liquid, a);
    check_rel(b.Pi, 1.0, 1.0e-14, "Pi = delta_p_exit / (gamma/a)");
    check(b.within_capillary_range, "|Pi| <= 2 liegt im Kapillarbereich");

    q.p_reservoir = 3.0 * liquid.gamma / a;
    const PressureBudget b3 = solve_pressure_budget(q, liquid, a);
    check(!b3.within_capillary_range,
          "|Pi| = 3 liegt ausserhalb; das steht im Ergebnis, bevor der Meniskusloeser "
          "gefragt wird");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
