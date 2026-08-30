// tests/test_transport.cpp -- P3: finite conductivity and pipe flow
//
// Both halves are checked against CLOSED FORMS, and the closed forms are
// written out here a second time from their definitions rather than called, so
// that the test cannot agree with the code by sharing a bug with it.
//
//   * pipe flow: the parabolic profile, the flow rate, the wall shear stress
//     and the mesh convergence of all three;
//   * charge transport: the relaxation time, the exponential decay, the steady
//     current in a cylinder, and -- the one that matters for the physics -- that
//     no current leaves a surface that carries the zero-flux condition.
//
// The last group is the honest limit statement: tau is computable only where
// BOTH eps_r and sigma are documented, and for EMI-BF4 eps_r is not.  The
// perfect-conductor limit of P3b is therefore not justified by this project's
// own material data -- it is assumed.  That is tested, so it cannot be lost.

#include <cmath>
#include <utility>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/material_data.hpp"
#include "es/transport.hpp"

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
  const Real e =
      (std::abs(want) > 0.0) ? std::abs(got - want) / std::abs(want) : std::abs(got - want);
  std::printf("  [%s] %s: %.9e gegen %.9e, rel. %.3e (Grenze %.1e)\n",
              (e <= tol) ? "ok" : "FEHLER", what.c_str(), got, want, e, tol);
  if (!(e <= tol)) ++failures;
}

}  // namespace

int main() {
  std::printf("P3 -- endliche Leitfaehigkeit und Zulaufstroemung\n\n");

  // =========================================================================
  std::printf("1. Voll ausgebildete Rohrstroemung gegen Hagen-Poiseuille\n");
  {
    const Real R = 5.0e-6, L = 3.0e-5, mu = 0.03637, dpdz = -1.0e9;
    const PipeFlowSolution f = solve_pipe_flow(R, L, mu, dpdz, 41, 21);
    std::printf("    %lld Knoten, FEM-Residuum %.3e\n", static_cast<long long>(f.n_nodes),
                f.fem_residual);
    // The closed forms, written out from their definitions.
    const Real u0 = -dpdz * R * R / (4.0 * mu);
    const Real Q = -dpdz * pi * R * R * R * R / (8.0 * mu);
    const Real tau_w = std::abs(dpdz) * R / 2.0;
    check_rel(f.centreline_velocity, u0, 2.0e-3, "Mittengeschwindigkeit");
    check_rel(f.flow_rate, Q, 2.0e-3, "Volumenstrom");
    check_rel(f.wall_shear_stress, tau_w, 5.0e-2, "Wandschubspannung (Einpunktgradient)");
    check(f.max_profile_error < 5.0e-3, "das Profil ist die Parabel");

    // Mesh convergence.  A Q1 element on a quadratic exact solution should be
    // second order in the flow rate; measured, not assumed.
    std::vector<Real> err;
    for (Index n : {11, 21, 41, 81}) {
      const PipeFlowSolution g = solve_pipe_flow(R, L, mu, dpdz, n, 9);
      err.push_back(std::abs(g.flow_rate - Q) / Q);
    }
    std::printf("    Fehler des Volumenstroms ueber die Netzstufen: ");
    for (Real e : err) std::printf("%.3e ", e);
    std::printf("\n");
    const Real order = std::log(err[0] / err[3]) / std::log(8.0);
    std::printf("    beobachtete Ordnung: %.3f\n", order);
    check(order > 1.7, "der Fehler faellt mindestens zweiter Ordnung");

    // The reduction is linear in dp/dz and in 1/mu -- both exactly.
    const PipeFlowSolution g = solve_pipe_flow(R, L, mu, 2.0 * dpdz, 41, 21);
    check_rel(g.flow_rate, 2.0 * f.flow_rate, 1.0e-12, "der Volumenstrom ist linear in dp/dz");
    const PipeFlowSolution h = solve_pipe_flow(R, L, 2.0 * mu, dpdz, 41, 21);
    check_rel(h.flow_rate, 0.5 * f.flow_rate, 1.0e-12, "und umgekehrt proportional zu mu");

    // No flow without a pressure gradient -- exactly zero, not small.
    const PipeFlowSolution z = solve_pipe_flow(R, L, mu, 0.0, 21, 9);
    check(z.flow_rate == 0.0 && z.centreline_velocity == 0.0,
          "ohne Druckgradient ist die Geschwindigkeit exakt null");
  }

  // =========================================================================
  std::printf("\n2. Der Widerstand von P1 gegen die geloeste Stroemung\n");
  {
    // P1 asserts R_h = 8 mu L / (pi R^4).  Here the field is solved and the
    // profile integrated.  The two share NO code.
    const Real R = 5.0e-6, L = 3.0e-4, mu = 0.03637;
    const Real dp = 1.0e4;                 // a pressure drop over the length
    const PipeFlowSolution f = solve_pipe_flow(R, L, mu, -dp / L, 81, 9);
    const Real Rh_solved = dp / f.flow_rate;
    const Real Rh_p1 = 8.0 * mu * L / (pi * std::pow(R, 4.0));
    check_rel(Rh_solved, Rh_p1, 2.0e-3,
              "der geloeste hydraulische Widerstand trifft die geschlossene Form von P1");
  }

  // =========================================================================
  std::printf("\n3. Ladungsrelaxation\n");
  {
    const Real eps_r = 12.8, sigma = 1.5584;
    const Real tau = charge_relaxation_time(eps_r, sigma);
    check_rel(tau, eps0 * eps_r / sigma, 1.0e-15, "tau = eps0 eps_r / sigma");
    std::printf("    tau = %.4e s bei eps_r = %.2f und sigma = %.4f S/m\n", tau, eps_r, sigma);

    // The closed-form decay, and that it IS the solution of the ODE it claims.
    const Real rho0 = 1.0e-3;
    check_rel(relaxed_charge_density(rho0, 0.0, tau), rho0, 1.0e-15, "rho(0) = rho0");
    check_rel(relaxed_charge_density(rho0, tau, tau), rho0 / std::exp(1.0), 1.0e-14,
              "rho(tau) = rho0 / e");
    {
      // d rho/dt + rho/tau = 0, checked by a central difference of the returned
      // function -- the ODE, not a second copy of the exponential.
      const Real t = 0.7 * tau, h = 1.0e-6 * tau;
      const Real d = (relaxed_charge_density(rho0, t + h, tau) -
                      relaxed_charge_density(rho0, t - h, tau)) / (2.0 * h);
      check_rel(d, -relaxed_charge_density(rho0, t, tau) / tau, 1.0e-8,
                "die Zerfallsfunktion loest d rho/dt = -rho/tau");
    }
    check(!std::isfinite(charge_relaxation_time(0.0, sigma)),
          "ohne Permittivitaet gibt es kein tau, sondern nan");
    check(!std::isfinite(charge_relaxation_time(eps_r, 0.0)),
          "ohne Leitfaehigkeit ebenso");
  }

  // =========================================================================
  std::printf("\n4. Der Perfect-Conductor-Grenzfall ist ein VERHAELTNIS\n");
  {
    const Real eps_r = 12.8, sigma = 1.5584;
    const Real tau = charge_relaxation_time(eps_r, sigma);
    const RelaxationVerdict fast =
        judge_conductor_limit_explicit(eps_r, sigma, 1.0e4 * tau);
    fast.print(stdout);
    check(fast.limit == ConductorLimit::PerfectConductorJustified,
          "bei einer Prozesszeit weit ueber tau ist der Grenzfall gerechtfertigt");
    const RelaxationVerdict slow = judge_conductor_limit_explicit(eps_r, sigma, tau);
    check(slow.limit == ConductorLimit::FiniteConductivityRequired,
          "bei einer Prozesszeit von der Groesse tau nicht mehr");
    check(slow.tau == fast.tau, "tau haengt nicht von der Prozesszeit ab -- das Urteil schon");

    // AND THE POINT: for this project's own material data it is not computable.
    const RelaxationVerdict real =
        judge_conductor_limit(emibf4_sourced(), 298.15, 1.0e-6);
    real.print(stdout);
    check(real.limit == ConductorLimit::MissingMaterialData,
          "mit dem belegten Stoffdatensatz ist tau NICHT berechenbar: eps_r fehlt");
    check(!std::isfinite(real.tau), "und tau ist nan, nicht ein plausibler Ersatzwert");
    check(real.sigma_status != MaterialDataStatus::MissingMaterialData,
          "die Leitfaehigkeit ist dagegen belegt -- es fehlt genau eine der beiden Zahlen");
  }

  // =========================================================================
  // WHICH permittivity belongs in tau_q.  The earlier version of this project
  // demanded one "DC permittivity" and threw the 1-18 GHz measurement away as
  // "not DC".  Both halves of that were wrong: the free charge decays on the
  // time scale tau itself, so the permittivity that belongs in tau is the one
  // at f* = 1/(2 pi tau) -- which for this liquid lies in the low GHz range,
  // i.e. exactly where the discarded data were measured.
  std::printf("\n4b. Welche Permittivitaet in tau_q gehoert\n");
  {
    const MaterialDataset& d = emibf4_sourced();

    // The band: no selection, no average, and nothing below the electrode
    // polarisation floor -- concept (1) is a property of the cell.
    const PermittivityBand band = permittivity_band(d, 298.15);
    band.print(stdout);
    check(band.ok, "es gibt ein begruendetes eps_r-Band");
    check(band.n_frequency_points > 1,
          "der Datensatz enthaelt eine frequenzaufgeloeste Messung -- sie ist nicht "
          "geloescht worden");
    check(band.n_static_points > 0,
          "und daneben als statisch berichtete Werte aus Mikrowellenspektren");
    check(!band.any_below_polarisation_floor,
          "kein einziger Punkt liegt unter der Elektrodenpolarisationsschwelle: der "
          "Datensatz enthaelt gar keine kHz-Scheinpermittivitaet");
    check(band.lo > 0.0 && band.hi > band.lo, "das Band ist nicht entartet");

    // The single-value query STILL fails closed.  Nothing here relaxes that.
    const MaterialValue single = material_value(d, PropertyKind::RelativePermittivity, 298.15);
    check(single.status == MaterialDataStatus::MissingMaterialData,
          "ein EINZELNER eps_r-Wert bleibt MissingMaterialData: keine Quelle nennt "
          "Reinheit und Wassergehalt");

    // The implicit equation, solved on the measured curve.
    const SelfConsistentRelaxation sc = self_consistent_relaxation(d, 298.15);
    sc.print(stdout);
    check(sc.ok, "die implizite Gleichung tau = eps0 eps_r(1/(2 pi tau))/K ist loesbar");
    if (sc.ok) {
      check(sc.residual < 1e-10,
            "die Fixpunktiteration ist konvergiert: eps_r(f*) trifft das benutzte eps_r");
      // The two defining identities hold EXACTLY, checked rather than trusted.
      check(sc.tau == constants::eps0 * sc.eps_r / sc.sigma,
            "das Ergebnis erfuellt tau = eps0 eps_r / K exakt, mit genau dem eps_r bei f*");
      check(sc.f_star == 1.0 / (2.0 * constants::pi * sc.tau),
            "und f* = 1/(2 pi tau) exakt, mit genau diesem tau");
      check(sc.f_star_inside_measured,
            "f* liegt INNERHALB des gemessenen Frequenzbereichs 1-18 GHz -- die Kurve "
            "wird nicht extrapoliert und keine Dispersionsfunktion erfunden");
      check(sc.eps_r < sc.eps_static,
            "eps_r bei f* ist kleiner als der statische Wert, weil die Dispersion bei "
            "GHz bereits abgefallen ist");
      check(sc.tau < sc.tau_static,
            "das selbstkonsistente tau ist deshalb kuerzer als das mit eps_s gerechnete");
      // The whole difference between the two is a few tens of per cent -- worth
      // reporting, and far too small to change the verdict below.
      const Real rel = (sc.tau_static - sc.tau) / sc.tau;
      std::printf("    Unterschied zwischen tau(eps_s) und tau(eps(f*)): %+.1f %%\n",
                  100.0 * rel);
      check(rel > 0.0 && rel < 1.0,
            "der Unterschied ist merklich, aber keine Groessenordnung");
    }

    // THE VERDICT, taken at the corner of the band that is WORST for the
    // approximation.  This is the point of the whole section: the
    // equipotential treatment of P3b can be justified from documented data
    // over the entire justified range, without any single unsourced eps_r.
    const Real a = 5.0e-6, rho = 1280.9, gamma = 0.05401, mu = 0.03637;
    const Real t_cap = std::sqrt(rho * a * a * a / gamma);
    const Real t_visc = mu * a / gamma;
    for (auto tp : {std::pair<const char*, Real>{"t_kap", t_cap},
                    std::pair<const char*, Real>{"t_vis", t_visc}}) {
      std::printf("    -- Prozesszeit %s = %.4e s\n", tp.first, tp.second);
      const BandedRelaxationVerdict v =
          judge_conductor_limit_over_band(d, 298.15, tp.second);
      v.print(stdout);
      check(v.limit == ConductorLimit::PerfectConductorJustified,
            std::string("ueber das GANZE begruendete Band ist der Perfect-Conductor-Fall "
                        "fuer ") + tp.first + " gerechtfertigt");
      check(v.ratio_min <= v.ratio_max, "die Verhaeltnisse sind richtig herum sortiert");
      check(v.tau_max >= v.tau_min, "und die Zeiten ebenso");
      check(v.tau_self_consistent >= v.tau_min && v.tau_self_consistent <= v.tau_max,
            "das selbstkonsistente tau liegt im Band der vier Ecken");
      // The worst corner is the largest eps_r with the smallest K.  Checked,
      // not asserted: a sign error in that reasoning would be invisible.
      check(std::abs(v.tau_max - constants::eps0 * v.eps_hi / v.sigma_lo) <=
                1e-14 * v.tau_max,
            "die ungueenstigste Ecke ist tatsaechlich (eps_hi, K_lo)");
      check(v.ratio_min > 100.0 * transport::kPerfectConductorMargin,
            std::string("und sie verfehlt die Schranke nicht knapp, sondern um mehr als "
                        "zwei Groessenordnungen (") + tp.first + ")");
    }

    // A counter-example, so that the verdict is not simply always positive:
    // for a process time of the order of tau itself it must flip.
    const BandedRelaxationVerdict fastp =
        judge_conductor_limit_over_band(d, 298.15, 1.0e-10);
    check(fastp.limit == ConductorLimit::FiniteConductivityRequired,
          "bei einer Prozesszeit von 0,1 ns kippt das Urteil -- die Schranke wirkt");
  }

  // =========================================================================
  std::printf("\n5. Stationaerer Leitungsstrom im Zylinder\n");
  {
    const Real R = 5.0e-6, L = 3.0e-4, sigma = 1.5584, V = 1.0;
    const ConductionSolution c = solve_cylinder_conduction(R, L, sigma, V, 21, 41);
    std::printf("    %lld Knoten, FEM-Residuum %.3e\n", static_cast<long long>(c.n_nodes),
                c.fem_residual);
    // Closed forms, written from their definitions.
    const Real I = sigma * pi * R * R * V / L;
    const Real Rohm = L / (sigma * pi * R * R);
    check_rel(c.current, I, 1.0e-9, "Strom I = sigma A V / L, aus der Knotenreaktion");
    check_rel(c.resistance, Rohm, 1.0e-9, "Widerstand R = L / (sigma A)");
    std::printf("    groesster Knotenfehler des Potentials: %.3e von V\n",
                c.max_potential_error);
    // The Q1 solution of a linear exact solution IS that solution, so what is
    // left is the round-off of the direct solve -- not a discretisation error.
    check(c.max_potential_error < 1.0e-10,
          "das Potential ist die exakte lineare Loesung bis auf die Rundung des Loesers");

    // NO CURRENT LEAVES A SURFACE THAT MAY NOT CARRY ONE.  This is the physics
    // statement, not a numerical nicety: without emission there is nowhere for
    // charge crossing the free surface to go.
    std::printf("    seitliche Leckstromdichte / axiale: %.3e\n", c.lateral_leakage);
    check(c.lateral_leakage < 1.0e-12,
          "durch die zero-flux-Flaeche flieszt kein Strom -- ohne Emission darf keiner");

    // Ohm's law, term by term.
    const ConductionSolution c2 = solve_cylinder_conduction(R, L, sigma, 2.0 * V, 21, 41);
    check_rel(c2.current, 2.0 * c.current, 1.0e-12, "der Strom ist linear in der Spannung");
    const ConductionSolution cs = solve_cylinder_conduction(R, L, 2.0 * sigma, V, 21, 41);
    check_rel(cs.current, 2.0 * c.current, 1.0e-9, "und linear in sigma");
    const ConductionSolution cl = solve_cylinder_conduction(R, 2.0 * L, sigma, V, 21, 41);
    check_rel(cl.current, 0.5 * c.current, 1.0e-9,
              "und umgekehrt proportional zur Laenge");
    const ConductionSolution cr = solve_cylinder_conduction(2.0 * R, L, sigma, V, 21, 41);
    check_rel(cr.current, 4.0 * c.current, 1.0e-9, "und proportional zur Querschnittsflaeche");

    // Polarity.
    const ConductionSolution cn = solve_cylinder_conduction(R, L, sigma, -V, 21, 41);
    check_rel(cn.current, c.current, 1.0e-12,
              "die Stromstaerke haengt nicht vom Vorzeichen der Spannung ab");
    check(cn.phi[cn.phi.size() - 1] < 0.0, "das Potential dagegen schon");

    // Zero drive, zero current -- exactly.
    const ConductionSolution c0 = solve_cylinder_conduction(R, L, sigma, 0.0, 11, 11);
    check(c0.current == 0.0, "ohne Spannung flieszt exakt kein Strom");
  }

  // =========================================================================
  std::printf("\n6. Die Zeitskalen nebeneinander (Rechenbeispiel)\n");
  {
    // Everything here is a RATIO of quantities computed above; no new physics.
    // tau comes from 4b: the self-consistent solution on the MEASURED dispersion
    // curve.  No unsourced eps_r is used here any more.
    const SelfConsistentRelaxation sc6 = self_consistent_relaxation(emibf4_sourced(), 298.15);
    check(sc6.ok, "die selbstkonsistente Loesung aus 4b liegt vor");
    const Real eps_r = sc6.eps_r;
    const Real tau = sc6.tau;
    const Real a = 5.0e-6, rho = 1280.9, gamma = 0.05401, mu = 0.03637;
    const Real t_cap = std::sqrt(rho * a * a * a / gamma);   // capillary-inertial
    const Real t_visc = mu * a / gamma;                       // visco-capillary
    std::printf("    tau (Ladung, selbstkonsistent bei f* = %.3e Hz, eps_r = %.4f)"
                " = %.4e s\n", sc6.f_star, eps_r, tau);
    std::printf("    t_kap = sqrt(rho a^3/gamma) = %.4e s\n", t_cap);
    std::printf("    t_vis = mu a / gamma        = %.4e s\n", t_visc);
    std::printf("    t_kap/tau = %.3e, t_vis/tau = %.3e\n", t_cap / tau, t_visc / tau);
    check(t_cap / tau > transport::kPerfectConductorMargin,
          "gegen die kapillare Zeitskala ist die Ladungsrelaxation sehr schnell");
    check(t_visc / tau > transport::kPerfectConductorMargin,
          "gegen die viskokapillare ebenso");
    std::printf("    Das rechtfertigt den Aequipotentialansatz fuer die STATISCHE Form --\n"
                "    und nach 4b sogar ueber das ganze begruendete eps_r/K-Band, ohne\n"
                "    dass ein einzelner unbelegter eps_r-Wert benutzt wird.\n"
                "    Es sagt NICHTS ueber einen emittierenden Betrieb: dort ist die\n"
                "    Prozesszeit die Transitzeit durch die Emissionszone, und die ist\n"
                "    hier nicht gerechnet.\n");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
