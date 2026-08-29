// tests/test_load_projection.cpp -- P0: the load projection of P3b
//
// WHAT IS TESTED HERE, AND WHY IT IS TESTED WITH MANUFACTURED LOADS
//
// P3b hands the capillary solver a surface load and claims three things about
// it: the handed load is CONTINUOUS, it carries the INTEGRATED MAXWELL FORCE,
// and nothing is clipped or excluded near the contact edge.  With a solved
// field none of the three can be checked, because the true answer is unknown.
// With a PRESCRIBED load on a prescribed surface all three can, because the
// integral is known in closed form.
//
// Two loads, as required:
//   * a smooth analytic one, p = p0 (1 + (r/a)^2) on the flat disc, whose force
//     is (3/2) pi p0 a^2;
//   * an integrable singularity p = C d^beta with -1 < beta < 0, whose force is
//     2 pi C a^(2+beta) / ((1+beta)(2+beta)).
//
// The last groups settle the criterion P3b reported as missed: the closed form
// of edge_gate::kTolExclusion's measure is derived and pinned against the
// measured one, and the measure is shown to converge to a NON-ZERO constant
// under refinement while the discretisation error of the very same load falls.
// That is why the criterion was replaced, and it is a derivation, not a swap.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/load_projection.hpp"

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

constexpr Real kA = 7.5e-6;          ///< pinning radius of the P1 device [m]
constexpr Real kGammaOverA = 5.0e3;  ///< only a reporting scale here

FreeSurface flat() { return FreeSurface::flat_surface(kA, 0.0); }

}  // namespace

// ===========================================================================

int main() {
  std::printf("P0 -- Lastprojektion: Stetigkeit, Krafterhaltung, Kante\n\n");

  // =========================================================================
  // 1.  A smooth analytic load: the whole chain against a closed form
  // =========================================================================
  std::printf("1. Glatte analytische Last p = p0 (1 + (r/a)^2)\n");
  {
    const Real p0 = 1.0e4;
    const FreeSurface fs = flat();
    auto p_of_d = [&](Real d) {
      const Real r = kA - d;  // flat surface: arclength IS the radius
      return p0 * (1.0 + (r / kA) * (r / kA));
    };
    const Real F_exact = flat_disc_smooth_force(p0, kA);

    // Convergence of the segment projection: trapezoidal on a smooth integrand,
    // so second order.  Measured, not assumed.
    std::vector<Real> err;
    for (Index n : {16, 32, 64, 128}) {
      const MaxwellLoad L =
          manufactured_load(fs, uniform_radius_nodes(kA, n), p_of_d, kGammaOverA);
      err.push_back(std::abs(L.total_force - F_exact) / F_exact);
    }
    std::printf("    Fehler der Segmentprojektion: ");
    for (Real e : err) std::printf("%.3e ", e);
    std::printf("\n");
    const Real order = std::log(err[0] / err[3]) / std::log(8.0);
    std::printf("    beobachtete Ordnung ueber 16 -> 128 Segmente: %.3f\n", order);
    check(order > 1.8 && order < 2.2, "Segmentprojektion ist zweiter Ordnung");

    const MaxwellLoad L =
        manufactured_load(fs, uniform_radius_nodes(kA, 64), p_of_d, kGammaOverA);
    const LoadProjectionAudit A = audit_projection(L, fs, "glatt", F_exact);
    A.print(stdout);

    check_rel(A.segment_force, F_exact, 1.0e-3, "Segmentkraft gegen die geschlossene Form");
    check_rel(A.bin_force, A.segment_force, 1.0e-13, "Binkraft gegen die Segmentkraft");
    check_rel(A.handed_force_reconstructed, A.segment_force, 1.0e-13,
              "uebergebene Last traegt die Kraft (gegen A')");
    check(A.error_handed_true < 1.0e-4,
          "uebergebene Last traegt die Kraft auch gegen 2 pi r ds");

    // CONTINUITY.  The DECAY of the jump with the probe offset is the
    // statement, not a threshold: a continuous load with a kink gives about
    // 0.1 when the offset is divided by ten, a staircase gives about 1.
    check(A.handed_jump_decay < 0.2,
          "der Sprung der uebergebenen Last faellt mit dem Probenabstand -- sie ist stetig");
    check(A.staircase_jump_decay > 0.9,
          "der Sprung der Treppenfunktion faellt NICHT -- der Test oben hat Zaehne");
    check(A.handed_jump_ratio < 1.0e-4 * A.staircase_jump_ratio,
          "und er ist vier Groessenordnungen kleiner als der der Treppenfunktion");

    // NOT A STAIRCASE.  A piecewise-constant reconstruction would be flat
    // between the bin centres.  The handed load must follow the prescribed one.
    {
      Real worst = 0.0;
      const ProjectedLoad pl = ProjectedLoad::from(L);
      for (int k = 1; k < 120; ++k) {
        const Real tau = static_cast<Real>(k) / 128.0 + 0.5 / 128.0;
        const Real want = p_of_d(kA - tau * kA);
        worst = std::max(worst, std::abs(pl.at(tau) - want) / want);
      }
      std::printf("    groesster Abstand der uebergebenen Last zur vorgeschriebenen: %.3e\n",
                  worst);
      check(worst < 1.0e-2, "uebergebene Last folgt der vorgeschriebenen Last");
    }
  }

  // =========================================================================
  // 2.  An integrable singularity p = C d^beta, -1 < beta < 0
  // =========================================================================
  std::printf("\n2. Integrable Singularitaet p = C d^beta\n");
  for (Real beta : {-0.25, -0.44, -0.75}) {
    const Real C = 1.0e4 * std::pow(kA, -beta);  // so that p(a) ~ 1e4 Pa
    const FreeSurface fs = flat();
    auto p_of_d = [&](Real d) { return (d > 0.0) ? C * std::pow(d, beta) : 0.0; };
    const Real F_exact = flat_disc_power_law_force(C, beta, kA);

    std::printf("  beta = %+.2f, F_exakt = %.9e N\n", beta, F_exact);
    std::vector<Real> err;
    for (Index n : {32, 64, 128, 256}) {
      const MaxwellLoad L =
          manufactured_load(fs, uniform_radius_nodes(kA, n), p_of_d, kGammaOverA);
      err.push_back(std::abs(L.total_force - F_exact) / F_exact);
    }
    std::printf("    Fehler der Segmentprojektion: ");
    for (Real e : err) std::printf("%.3e ", e);
    std::printf("\n");
    // The last segment carries the singularity, and a trapezoid on it is wrong
    // by O(h^(1+beta)) of the total.  That is the rate the error must show --
    // an honest statement of what the projection can do, not a claim that it is
    // second order.
    const Real order = std::log(err[0] / err[3]) / std::log(8.0);
    std::printf("    beobachtete Ordnung: %.3f, erwartet 1+beta = %.3f\n", order, 1.0 + beta);
    check(std::abs(order - (1.0 + beta)) < 0.15,
          "der Fehler faellt wie h^(1+beta) -- die Singularitaet ist integriert, nicht gekappt");

    const MaxwellLoad L =
        manufactured_load(fs, uniform_radius_nodes(kA, 128), p_of_d, kGammaOverA);
    const LoadProjectionAudit A = audit_projection(L, fs, "singulaer", F_exact);

    check_rel(A.bin_force, A.segment_force, 1.0e-13, "Binkraft gegen die Segmentkraft");
    check_rel(A.handed_force_reconstructed, A.segment_force, 1.0e-13,
              "uebergebene Last traegt die Kraft (gegen A')");
    check(A.handed_jump_decay < 0.2, "uebergebene Last ist stetig");
    check(A.staircase_jump_decay > 0.9, "die Treppenfunktion springt");

    // NO HIDDEN EXCLUSION ZONE.  The bin nearest the edge must carry force, the
    // tau range must reach 1, and no bin may be empty.
    check(A.empty_bins == 0, "kein leeres Bin: die Projektion deckt tau in [0,1] ab");
    check(A.tau_last == 1.0, "die Segmente reichen bis tau = 1, also bis an die Kontaktlinie");
    check(A.last_bin_force_fraction > 0.0,
          "das kantennaechste Bin traegt Kraft -- keine Ausschlusszone");
    check(A.area_gap < 1.0e-13, "die Binflaechen summieren sich zur Segmentflaeche");

    // NO CLIPPING.  The handed load must reach the bin maximum; a cap anywhere
    // would pin it below.
    check(A.max_handed_pressure >= 0.98 * A.max_bin_pressure,
          "die uebergebene Last erreicht das Binmaximum -- keine Kappung");
  }

  // =========================================================================
  // 3.  The contact edge: what the projection does with the last node
  // =========================================================================
  std::printf("\n3. Behandlung der Kontaktkante\n");
  {
    const Real beta = -0.44;
    const Real C = 1.0e4 * std::pow(kA, -beta);
    const FreeSurface fs = flat();
    auto p_of_d = [&](Real d) { return (d > 0.0) ? C * std::pow(d, beta) : 0.0; };

    // The pointwise nodal value at the edge GROWS without bound under
    // refinement -- that is the singularity, and it is what P3b refuses to use.
    // The FORCE the last segment carries must fall to zero like h^(1+beta).
    // Both are measured.
    std::vector<Real> peak, last_force;
    for (Index n : {32, 64, 128, 256}) {
      const MaxwellLoad L =
          manufactured_load(fs, uniform_radius_nodes(kA, n), p_of_d, kGammaOverA);
      peak.push_back(L.node_pM.back());
      last_force.push_back(L.seg_force.back() / L.total_force);
    }
    std::printf("    punktweiser Kantenwert: ");
    for (Real v : peak) std::printf("%.4e ", v);
    std::printf("Pa\n    Kraftanteil des letzten Segments: ");
    for (Real v : last_force) std::printf("%.4e ", v);
    std::printf("\n");
    check(peak[3] > peak[0], "der punktweise Kantenwert waechst mit der Verfeinerung");
    check(last_force[3] < last_force[0],
          "der Kraftanteil des letzten Segments faellt -- die Last ist integrabel");
    const Real order = std::log(last_force[0] / last_force[3]) / std::log(8.0);
    std::printf("    Ordnung des Kraftanteils: %.3f, erwartet 1+beta = %.3f\n", order,
                1.0 + beta);
    check(std::abs(order - (1.0 + beta)) < 0.15,
          "der Kraftanteil des letzten Segments faellt wie h^(1+beta)");
  }

  // =========================================================================
  // 4.  Why edge_gate::kTolExclusion cannot be met -- the closed form
  // =========================================================================
  std::printf("\n4. Das Kriterium kTolExclusion: geschlossene Form gegen Messung\n");
  {
    const Real d0 = edge_gate::kExclusionMid;  // 0.05 a
    for (Real beta : {-0.25, -0.44, -0.75}) {
      const Real C = 1.0e4 * std::pow(kA, -beta);
      const FreeSurface fs = flat();
      auto p_of_d = [&](Real d) { return (d > 0.0) ? C * std::pow(d, beta) : 0.0; };

      const Real limit = exclusion_halving_limit(beta, d0);

      // Measured exactly as run_edge_gate() measures it, on four refinements.
      std::vector<Real> measured;
      for (Index n : {64, 128, 256, 512}) {
        const MaxwellLoad L =
            manufactured_load(fs, uniform_radius_nodes(kA, n), p_of_d, kGammaOverA);
        const Real f_mid = L.force_beyond(edge_gate::kExclusionMid * kA);
        const Real f_fine = L.force_beyond(edge_gate::kExclusionFine * kA);
        measured.push_back(std::abs(f_fine - f_mid) / std::abs(f_fine));
      }
      std::printf("  beta = %+.2f: Grenzwert %.4f, gemessen ", beta, limit);
      for (Real v : measured) std::printf("%.4f ", v);
      std::printf("\n");

      // (a) The measure converges to the closed-form limit.
      check_rel(measured[3], limit, 2.0e-2,
                "gemessene Halbierungsaenderung gegen die geschlossene Form");
      // (b) It does NOT go to zero: it settles well above the bound.
      check(measured[3] > 1.5 * edge_gate::kTolExclusion,
            "die Halbierungsaenderung bleibt ueber kTolExclusion -- "
            "sie ist kein Diskretisierungsfehler");
      // (c) It stops moving under refinement, while the discretisation error of
      //     the SAME load keeps falling at the rate the singularity allows.
      //     That contrast is the whole argument.
      const Real drift = std::abs(measured[3] - measured[2]) / measured[3];
      const Real F_exact = flat_disc_power_law_force(C, beta, kA);
      const MaxwellLoad L512 =
          manufactured_load(fs, uniform_radius_nodes(kA, 512), p_of_d, kGammaOverA);
      const MaxwellLoad L64 =
          manufactured_load(fs, uniform_radius_nodes(kA, 64), p_of_d, kGammaOverA);
      const Real e512 = std::abs(L512.total_force - F_exact) / F_exact;
      const Real e64 = std::abs(L64.total_force - F_exact) / F_exact;
      const Real expected = std::pow(8.0, -(1.0 + beta));
      std::printf("    Drift der Messgroesse 256 -> 512: %.3e; Diskretisierungsfehler "
                  "%.3e -> %.3e, Verhaeltnis %.4f (erwartet %.4f)\n",
                  drift, e64, e512, e512 / e64, expected);
      check(drift < 0.05, "die Halbierungsaenderung steht unter Verfeinerung still");
      check(e512 / e64 < 1.25 * expected,
            "der Diskretisierungsfehler derselben Last faellt mit der erwarteten Rate weiter");
    }
    std::printf("  Damit ist kTolExclusion aufgeloest: die Groesse, die es prueft, hat einen\n"
                "  endlichen Grenzwert ungleich null und kann die Schranke nicht einhalten.\n"
                "  Die Schranke bleibt deklariert und wird als nicht eingehalten berichtet.\n");
  }

  // =========================================================================
  // 5.  The replacement criterion, against a KNOWN limit
  // =========================================================================
  std::printf("\n5. Das ersetzende Kriterium gegen einen bekannten Grenzwert\n");
  {
    // The gate's second extrapolation fits F(d0) = F_inf - K d0^(1+beta) over
    // the three exclusion distances.  On a manufactured load F_inf is known, so
    // the criterion can be checked against ground truth instead of against
    // itself.  That is what a replacement criterion has to survive.
    const Real beta = -0.44;
    const Real C = 1.0e4 * std::pow(kA, -beta);
    const FreeSurface fs = flat();
    auto p_of_d = [&](Real d) { return (d > 0.0) ? C * std::pow(d, beta) : 0.0; };
    const Real F_exact = flat_disc_power_law_force(C, beta, kA);
    const MaxwellLoad L =
        manufactured_load(fs, uniform_radius_nodes(kA, 256), p_of_d, kGammaOverA);

    const Real d[3] = {edge_gate::kExclusionCoarse * kA, edge_gate::kExclusionMid * kA,
                       edge_gate::kExclusionFine * kA};
    const Real F[3] = {L.force_beyond(d[0]), L.force_beyond(d[1]), L.force_beyond(d[2])};
    Real sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int k = 0; k < 3; ++k) {
      const Real x = std::pow(d[k], 1.0 + beta);
      sx += x;
      sy += F[k];
      sxx += x * x;
      sxy += x * F[k];
    }
    const Real n = 3.0;
    const Real slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    const Real F_inf = (sy - slope * sx) / n;
    std::printf("    Extrapolation ueber die Ausschlussdistanz: %.9e N\n", F_inf);
    std::printf("    exakt                                    : %.9e N\n", F_exact);
    check_rel(F_inf, F_exact, edge_gate::kTolLimitAgreement,
              "die Extrapolation ueber die Ausschlussdistanz trifft den bekannten Grenzwert");
  }

  // =========================================================================
  // 6.  Richardson and the discretisation verdict
  // =========================================================================
  std::printf("\n6. Richardson-Extrapolation\n");
  {
    // A manufactured second-order sequence on the sqrt(2) refinement ladder.
    const Real f_inf = 3.0, K = 0.5;
    std::vector<Real> v;
    for (int l = 0; l < 4; ++l) {
      const Real h = std::pow(kMeshLevelRatio, -l);
      v.push_back(f_inf - K * h * h);
    }
    const RichardsonEstimate e = richardson(v, kMeshLevelRatio);
    std::printf("    Ordnung %.6f, extrapoliert %.9f, Fehler der feinsten Stufe %.3e\n",
                e.observed_order, e.extrapolated, e.relative_error_finest);
    check(e.usable, "die Extrapolation ist auswertbar");
    check_rel(e.observed_order, 2.0, 1.0e-9, "beobachtete Ordnung");
    check_rel(e.extrapolated, f_inf, 1.0e-12, "extrapolierter Grenzwert");

    // A non-monotone sequence must NOT produce an error estimate.
    const RichardsonEstimate bad = richardson({1.0, 2.0, 1.5}, kMeshLevelRatio);
    check(!bad.usable, "eine nicht monotone Folge liefert keine Fehlerschaetzung");
    check(verdict_of(bad) == DiscretizationVerdict::NotInAsymptoticRange,
          "und sie bekommt den Status NotInAsymptoticRange");
    check(std::isfinite(bad.last_relative_change),
          "die blosse Aenderung bleibt bekannt und wird berichtet");

    // The verdict must key on the pre-declared one-per-cent target.
    std::vector<Real> slow;
    for (int l = 0; l < 4; ++l) slow.push_back(1.0 - 0.30 * std::pow(kMeshLevelRatio, -l));
    const RichardsonEstimate se = richardson(slow, kMeshLevelRatio);
    std::printf("    langsame Folge: Fehler der feinsten Stufe %.3e\n",
                se.relative_error_finest);
    check(verdict_of(se) == DiscretizationVerdict::DiscretizationNotConverged,
          "eine Folge ueber 1 % bekommt DiscretizationNotConverged");
    check(verdict_of(richardson({1.0, 2.0}, kMeshLevelRatio)) ==
              DiscretizationVerdict::InsufficientLevels,
          "zwei Stufen reichen nicht und sagen das");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
