// tests/test_surface_kinematics.cpp -- P4: the kinematic boundary condition
//
// Only ONE thing is implemented in P4 and only that is tested: the kinematic
// condition dx/dt . n = u . n on a PRESCRIBED velocity field.  Two fields are
// used, and they separate two different failure modes:
//
//   dilation  div u = 3 alpha, exact map x -> x e^{alpha t}.  The volume must
//             follow a known NON-ZERO change; a solver that conserved volume
//             by accident would fail here.
//   squeeze   div u = 0 exactly, exact map (r,z) -> (r e^{-alpha t/2},
//             z e^{alpha t}).  The volume must be conserved while the shape
//             changes a lot; a solver that got the shape wrong but the volume
//             right would fail the shape check.
//
// And the three prohibitions are tested as ABSENCES: a zero field must move
// nothing exactly, the tangential redistribution must not move the surface, and
// the dynamic solver must fail closed.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/surface_kinematics.hpp"

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

/// Largest distance between the advected nodes and the exact Lagrangian map,
/// relative to the initial radius.
Real shape_error(const SurfacePolyline& got, const SurfacePolyline& start,
                 Vec2 (*exact)(Vec2, Real, Real), Real alpha, Real t, Real R) {
  Real worst = 0.0;
  for (std::size_t k = 0; k < got.nodes.size(); ++k)
    worst = std::max(worst, norm(got.nodes[k] - exact(start.nodes[k], alpha, t)) / R);
  return worst;
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  std::printf("P4 -- kinematische Randbedingung auf vorgeschriebenem Feld\n\n");
  const Real R = 5.0e-6;

  // =========================================================================
  std::printf("1. Die Vorbedingung: der dynamische Loeser schlaegt geschlossen fehl\n");
  {
    bool threw = false;
    std::string msg;
    try {
      solve_dynamic_meniscus();
    } catch (const NotImplementedInThisPhase& e) {
      threw = true;
      msg = e.what();
    }
    check(threw, "solve_dynamic_meniscus wirft NotImplementedInThisPhase");
    check(msg.find("q_s E_t") != std::string::npos,
          "und nennt die fehlende tangentiale Traktion beim Namen");
    check(msg.find("eps_r") != std::string::npos, "und die fehlende Permittivitaet");
    check(msg.find("Mobilitaet") != std::string::npos,
          "und sagt ausdruecklich, dass es keine Mobilitaet gibt");
    check(msg.find("du_z/dz") != std::string::npos,
          "und nennt den strukturellen Grund, warum die P3-Stroemung nicht traegt");
  }

  // =========================================================================
  std::printf("\n2. Ein Nullfeld bewegt exakt nichts\n");
  {
    const SurfacePolyline s0 = hemisphere(R, 41);
    const VelocityField zero = [](Vec2, Real) { return Vec2{0.0, 0.0}; };
    const AdvectionResult a =
        advect_surface(s0, zero, 1.0e-9, 100, KinematicMode::Lagrangian, ContactLine::Free, 0.0);
    check(a.status == StepStatus::Ok, "der Lauf ist zulaessig");
    Real worst = 0.0;
    for (std::size_t k = 0; k < s0.nodes.size(); ++k)
      worst = std::max(worst, norm(a.surface.nodes[k] - s0.nodes[k]));
    check(worst == 0.0, "kein Knoten hat sich bewegt -- exakt null, nicht klein");
    check(a.volume_change == 0.0, "und das Volumen ist bitgenau unveraendert");
  }

  // =========================================================================
  std::printf("\n3. Dilatation: eine BEKANNTE Volumenaenderung\n");
  {
    const Real alpha = 1.0e6, T = 3.0e-7;   // alpha*T = 0.3
    const SurfacePolyline s0 = hemisphere(R, 81);
    const Real V0 = s0.revolved_volume(0.0);
    check_rel(V0, 2.0 / 3.0 * pi * R * R * R, 1.0e-3,
              "das Anfangsvolumen der Halbkugel ist (2/3) pi R^3");

    // The step counts are chosen so that the CFL-like bound of the time-step
    // contract is satisfied throughout: with alpha = 1e6 and 81 nodes on a 5 um
    // hemisphere the shortest segment is about 1e-7 m, so a node may move at
    // most 2.5e-8 m per step.
    std::vector<Real> verr, serr;
    bool all_ok = true;
    for (int n : {80, 160, 320, 640}) {
      const AdvectionResult a = advect_surface(s0, dilation_field(alpha), T / n, n,
                                              KinematicMode::Lagrangian, ContactLine::Free, 0.0);
      if (a.status != StepStatus::Ok) {
        std::printf("    n = %d: %s\n", n, to_string(a.status));
        all_ok = false;
        break;
      }
      const Real want = V0 * std::exp(3.0 * alpha * T);
      verr.push_back(std::abs(a.volume_final - want) / want);
      serr.push_back(shape_error(a.surface, s0, dilation_exact, alpha, T, R));
    }
    check(all_ok && verr.size() == 4, "alle vier Laeufe sind zulaessig");
    if (!all_ok || verr.size() != 4) {
      std::printf("\n%s: %d Fehler\n", "FEHLGESCHLAGEN", ++failures);
      return 1;
    }
    std::printf("    Volumenfehler gegen V0 e^{3 alpha T}: ");
    for (Real e : verr) std::printf("%.3e ", e);
    std::printf("\n    Formfehler gegen die exakte Abbildung: ");
    for (Real e : serr) std::printf("%.3e ", e);
    std::printf("\n");
    check(verr.back() < 1.0e-12, "das Volumen folgt der bekannten Aenderung");
    check(serr.back() < 1.0e-12, "und jeder Knoten der exakten Lagrange-Abbildung");
    // RK4 on an exponential is fourth order, but only while the discretisation
    // error is ABOVE the round-off floor.  At 640 steps the shape error is a few
    // times 1e-15 of R, i.e. of order 1e-20 m -- that is arithmetic, not
    // integration, and an order read off there would be meaningless.  So the
    // order is measured on the two COARSEST levels and the floor is stated.
    {
      const Real order = std::log(serr[0] / serr[1]) / std::log(2.0);
      std::printf("    beobachtete Zeitordnung (80 -> 160 Schritte): %.2f\n", order);
      std::printf("    bei 640 Schritten liegt der Fehler bei %.1e von R und damit auf dem "
                  "Rundungsboden\n", serr[3]);
      check(order > 3.5, "die Zeitintegration ist vierter Ordnung, solange sie es sein kann");
    }
  }

  // =========================================================================
  std::printf("\n4. Squeeze: div u = 0, also EXAKT erhaltenes Volumen\n");
  {
    const Real alpha = 1.0e6, T = 5.0e-7;
    const SurfacePolyline s0 = hemisphere(R, 81);
    const Real V0 = s0.revolved_volume(0.0);
    // The squeeze shortens the segments near the equator as it goes, so the
    // CFL-like bound TIGHTENS during the run: 200 steps stops at 180 with
    // StepTooLarge, which is the contract working and not a failure.
    const AdvectionResult a = advect_surface(s0, squeeze_field(alpha), T / 600, 600,
                                            KinematicMode::Lagrangian, ContactLine::Free, 0.0);
    a.print(stdout);
    check(a.status == StepStatus::Ok, "der Lauf ist zulaessig");
    std::printf("    Formaenderung: r schrumpft um %.3f, z waechst um %.3f\n",
                std::exp(-0.5 * alpha * T), std::exp(alpha * T));
    check(std::abs(a.volume_change) < 1.0e-12,
          "das Volumen ist erhalten, obwohl sich die Form stark aendert");
    check_rel(a.volume_final, V0, 1.0e-12, "und zwar auf den Anfangswert");
    const Real se = shape_error(a.surface, s0, squeeze_exact, alpha, T, R);
    std::printf("    Formfehler gegen die exakte Abbildung: %.3e\n", se);
    check(se < 1.0e-12, "und jeder Knoten trifft die exakte Abbildung");
  }

  // =========================================================================
  std::printf("\n5. Die tangentiale Umverteilung ist Netzbewegung, keine Physik\n");
  {
    // Same field, same time, two modes.  The SURFACE must be the same; only the
    // parametrisation differs.  Compared by the volume and the area, which are
    // properties of the surface and not of the nodes.
    const Real alpha = 1.0e6, T = 4.0e-7;
    const SurfacePolyline s0 = hemisphere(R, 161);
    const AdvectionResult lag = advect_surface(s0, squeeze_field(alpha), T / 400, 400,
                                              KinematicMode::Lagrangian, ContactLine::Free, 0.0);
    const AdvectionResult nrm = advect_surface(s0, squeeze_field(alpha), T / 400, 400,
                                              KinematicMode::NormalOnly, ContactLine::Free, 0.0);
    check(lag.status == StepStatus::Ok && nrm.status == StepStatus::Ok,
          "beide Laeufe sind zulaessig");
    std::printf("    Volumen Lagrange %.9e, NormalOnly %.9e\n", lag.volume_final,
                nrm.volume_final);
    std::printf("    Flaeche  Lagrange %.9e, NormalOnly %.9e\n", lag.surface.revolved_area(),
                nrm.surface.revolved_area());
    std::printf("    Ungleichmaessigkeit der Knotenabstaende: Lagrange %.3e, NormalOnly %.3e\n",
                lag.spacing_nonuniformity, nrm.spacing_nonuniformity);
    check(nrm.spacing_nonuniformity < 0.01 * lag.spacing_nonuniformity,
          "die Umverteilung macht die Knotenabstaende gleichmaessiger -- das ist ihr "
          "einziger Zweck");

    // HOW EXACTLY is "the surface is the same"?  Re-interpolating onto a
    // POLYLINE cuts corners, so the redistribution preserves the surface only to
    // the order of the node spacing.  That is MEASURED rather than asserted: the
    // discrepancy must fall with the node count, and at the expected rate.
    std::vector<Real> dv;
    for (std::size_t nn : {41u, 81u, 161u, 321u}) {
      const SurfacePolyline sn = hemisphere(R, nn);
      const AdvectionResult l = advect_surface(sn, squeeze_field(alpha), T / 900, 900,
                                              KinematicMode::Lagrangian, ContactLine::Free, 0.0);
      const AdvectionResult q = advect_surface(sn, squeeze_field(alpha), T / 900, 900,
                                              KinematicMode::NormalOnly, ContactLine::Free, 0.0);
      if (l.status != StepStatus::Ok || q.status != StepStatus::Ok) {
        std::printf("    n = %zu: %s / %s\n", nn, to_string(l.status), to_string(q.status));
        dv.clear();
        break;
      }
      dv.push_back(std::abs(q.volume_final - l.volume_final) / l.volume_final);
    }
    check(dv.size() == 4, "alle vier Vergleichslaeufe sind zulaessig");
    if (dv.size() == 4) {
      std::printf("    Volumenunterschied der beiden Moden, 41/81/161/321 Knoten: ");
      for (Real e : dv) std::printf("%.3e ", e);
      const Real order = std::log(dv[0] / dv[3]) / std::log(8.0);
      std::printf("\n    beobachtete Ordnung in der Knotenzahl: %.2f\n", order);
      check(dv[3] < dv[0], "der Unterschied faellt mit der Knotenzahl");
      // MEASURED, not assumed: the observed order is about 1.2, not 2.  The
      // per-redistribution loss is O(h^2), but it accumulates over the steps and
      // the accumulation is not h-independent, so the net rate is lower.  That
      // is a real limitation of this mesh motion and it is written down rather
      // than tightened away: a dynamic solver that redistributes every step
      // must either use a curved reconstruction or account for the loss.
      check(order > 1.0,
            "die Umverteilung ist flaechenerhaltend nur bis auf die Diskretisierung, "
            "und die gemessene Ordnung ist rund 1,2 -- nicht 2");
      check(dv[3] < 1.0e-3,
            "bei 321 Knoten liegt der Unterschied unter 1e-3 des Volumens");
    }
  }

  // =========================================================================
  std::printf("\n6. Der Zeitschrittvertrag lehnt ab statt zu warnen\n");
  {
    const SurfacePolyline s0 = hemisphere(R, 41);
    // A step so large that a node moves further than a quarter of the shortest
    // segment.  It must be REJECTED, not accepted with a note.
    const AdvectionResult big = advect_surface(s0, dilation_field(1.0e6), 1.0e-6, 5,
                                              KinematicMode::Lagrangian, ContactLine::Free, 0.0);
    std::printf("    grosser Schritt: %s nach %d Schritten\n", to_string(big.status),
                big.steps);
    check(big.status == StepStatus::StepTooLarge, "ein zu grosser Schritt wird abgelehnt");
    check(big.steps < 5, "und die Integration hoert dort auf");

    // A field that drives the surface through the axis.
    const VelocityField inward = [](Vec2, Real) { return Vec2{-1.0, 0.0}; };
    const AdvectionResult ax = advect_surface(s0, inward, 1.0e-9, 100000,
                                             KinematicMode::Lagrangian, ContactLine::Free, 0.0);
    check(ax.status == StepStatus::NodeCrossedAxis || ax.status == StepStatus::StepTooLarge,
          "eine Bewegung durch die Achse wird abgefangen");
  }

  // =========================================================================
  std::printf("\n7. Die gepinnte Kontaktlinie bleibt gepinnt\n");
  {
    const SurfacePolyline s0 = hemisphere(R, 41);
    const Vec2 pinned = s0.nodes.back();
    const AdvectionResult a = advect_surface(s0, squeeze_field(1.0e5), 1.0e-9, 200,
                                            KinematicMode::Lagrangian, ContactLine::Pinned, 0.0);
    check(a.status == StepStatus::Ok, "der Lauf ist zulaessig");
    check(a.surface.nodes.back().r == pinned.r && a.surface.nodes.back().z == pinned.z,
          "der Kontaktknoten ist bitgenau unveraendert");
    check(a.surface.nodes[0].r == 0.0, "und der Apex sitzt weiterhin exakt auf der Achse");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
