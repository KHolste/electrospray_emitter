// tests/test_space_charge.cpp -- P6: Poisson with free space charge
//
// The six mandatory checks, in order:
//   1. rho = 0 reproduces the Laplace solution -- BITWISE, not approximately;
//   2. a manufactured Poisson solution with known phi and rho;
//   3. the deposited total charge equals the sum of the macroparticle charges;
//   4. no divergence on approaching a single macroparticle;
//   5. mesh convergence of potential and field;
//   6. the polarity signs.
//
// Plus the one that separates the two contributions: the self field and the
// electrode field are solved separately and their superposition is exact,
// because the problem is linear and nothing in the discretisation breaks that.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/space_charge.hpp"

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

constexpr Real kR = 2.0e-5, kL = 4.0e-5;

struct Box {
  QuadMesh mesh;
  std::vector<Real> eps_r;
  std::vector<char> active, fixed;
  std::vector<Real> fixed_value;
};

/// A grounded cylinder, optionally with the top end held at `V`.
Box grounded_box(Index nr, Index nz, Real V = 0.0) {
  Box b;
  b.mesh = cylinder_mesh_symmetric(kR, kL, nr, nz);
  b.eps_r.assign(static_cast<std::size_t>(b.mesh.n_cells()), 1.0);
  b.active.assign(static_cast<std::size_t>(b.mesh.n_cells()), 1);
  b.fixed.assign(static_cast<std::size_t>(b.mesh.n_nodes()), 0);
  b.fixed_value.assign(static_cast<std::size_t>(b.mesh.n_nodes()), 0.0);
  for (Index j = 0; j < b.mesh.nz; ++j) {
    const Index n = b.mesh.node(b.mesh.nr - 1, j);
    b.fixed[static_cast<std::size_t>(n)] = 1;
  }
  for (Index i = 0; i < b.mesh.nr; ++i) {
    const Index a = b.mesh.node(i, 0), c = b.mesh.node(i, b.mesh.nz - 1);
    b.fixed[static_cast<std::size_t>(a)] = 1;
    b.fixed[static_cast<std::size_t>(c)] = 1;
    b.fixed_value[static_cast<std::size_t>(c)] = V;
  }
  return b;
}

}  // namespace

int main() {
  std::printf("P6 -- Poisson mit freier Raumladung\n\n");

  // =========================================================================
  std::printf("1. rho = 0 reproduziert die Laplace-Loesung\n");
  {
    const Box b = grounded_box(31, 61, 100.0);
    const std::vector<Real> zero_q(static_cast<std::size_t>(b.mesh.n_nodes()), 0.0);
    const SpaceChargeSolution s = solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed,
                                                         b.fixed_value, zero_q, {});
    Real worst = 0.0;
    for (std::size_t k = 0; k < s.phi.size(); ++k)
      worst = std::max(worst, std::abs(s.phi[k] - s.phi_no_charge[k]));
    std::printf("    groesster Knotenunterschied: %.3e V\n", worst);
    check(worst == 0.0,
          "mit einer identisch verschwindenden Ladung ist die Loesung BITGLEICH die "
          "Laplace-Loesung");
    check(s.max_field_shift == 0.0, "und das Feld ebenfalls bitgleich");
  }

  // =========================================================================
  std::printf("\n2. Hergestellte Poisson-Loesung mit bekanntem phi und rho\n");
  {
    const Real phi0 = 250.0;
    std::vector<Real> err_phi, err_E, err_E_axis, nodes;
    for (Index n : {21, 41, 81, 161}) {
      Box b = grounded_box(n, 2 * n - 1);
      // The manufactured solution vanishes on r = R and z = +-L, which is
      // exactly where this box is grounded -- so the boundary data are the
      // exact ones and the only error left is the discretisation.
      std::vector<Real> rho(static_cast<std::size_t>(b.mesh.n_nodes()), 0.0);
      for (Index j = 0; j < b.mesh.nz; ++j)
        for (Index i = 0; i < b.mesh.nr; ++i)
          rho[static_cast<std::size_t>(b.mesh.node(i, j))] =
              manufactured_charge_density(b.mesh.at(i, j), kR, kL, phi0);
      const SpaceChargeSolution s =
          solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed, b.fixed_value, {}, rho);
      // TWO field errors, because the nodal recovery is not the same thing at
      // the axis as in the interior.  In the interior it averages cells on all
      // sides and is superconvergent; on the axis (i = 0) and on the outer
      // column it is ONE-SIDED, and a one-sided average of a superconvergent
      // centre value is only first order.  Reporting a single number would hide
      // which of the two is being measured.
      Real wp = 0.0, we = 0.0, we_axis = 0.0;
      for (Index j = 1; j + 1 < b.mesh.nz; ++j)
        for (Index i = 0; i + 1 < b.mesh.nr; ++i) {
          const Vec2 x = b.mesh.at(i, j);
          wp = std::max(wp, std::abs(s.phi[static_cast<std::size_t>(b.mesh.node(i, j))] -
                                     manufactured_potential(x, kR, kL, phi0)));
          // dphi/dr and dphi/dz of the closed form, for the field comparison.
          const Vec2 Eexact{2.0 * phi0 * x.r * (kL * kL - x.z * x.z) / (kR * kR * kL * kL),
                            2.0 * phi0 * (kR * kR - x.r * x.r) * x.z / (kR * kR * kL * kL)};
          const Vec2 Egot =
              field_recovered_at_node(b.mesh, s.phi, b.eps_r, b.active, i, j, 1.0);
          const Real d = norm(Egot - Eexact);
          if (x.r < 0.25 * kR)
            we_axis = std::max(we_axis, d);   // the axis NEIGHBOURHOOD
          else
            we = std::max(we, d);             // away from the axis
        }
      err_phi.push_back(wp / phi0);
      err_E.push_back(we / (phi0 / kR));
      err_E_axis.push_back(we_axis / (phi0 / kR));
      nodes.push_back(static_cast<Real>(n));
    }
    std::printf("    Potentialfehler     : ");
    for (Real e : err_phi) std::printf("%.3e ", e);
    std::printf("\n    Feldfehler r > R/4  : ");
    for (Real e : err_E) std::printf("%.3e ", e);
    std::printf("\n    Feldfehler r < R/4  : ");
    for (Real e : err_E_axis) std::printf("%.3e ", e);
    std::printf("\n");
    const Real op = std::log(err_phi[0] / err_phi[3]) / std::log(8.0);
    const Real oe = std::log(err_E[0] / err_E[3]) / std::log(8.0);
    const Real oa = std::log(err_E_axis[0] / err_E_axis[3]) / std::log(8.0);
    std::printf("    Ordnung: Potential %.2f, Feld r > R/4 %.2f, Feld r < R/4 %.2f\n",
                op, oe, oa);
    check(err_phi.back() < 1.0e-4, "das Potential trifft die hergestellte Loesung");
    check(op > 1.8, "und konvergiert zweiter Ordnung");
    // A MEASURED PROPERTY OF THE EXISTING RECOVERY, and the measurement locates
    // it.  field_recovered_at_node() averages the superconvergent cell-centre
    // gradients with the CELL VOLUME as weight.  In axisymmetry that weight
    // carries a factor 2 pi r, so the two cells flanking a node at radius r get
    // weights in the ratio (r + h/2)/(r - h/2): an asymmetry of order h/r.
    //
    // Measured: away from the axis (r > R/4) the recovery is SECOND order, and
    // within a quarter radius of the axis it drops to FIRST order.  That is
    // exactly the h/r prediction, and it was found by splitting the error by
    // radius rather than by quoting one number for the whole domain.
    //
    // This matters beyond P6: the same recovery produced every surface field of
    // P2b, P2c and P3b.  It does not invalidate them -- their quantities are
    // integrals or far from the axis -- but a first-order field is what a PIC
    // loop would sample, and that has to be known before it is used.
    check(oe > 1.8,
          "abseits der Achse (r > R/4) ist die Rekonstruktion ZWEITER Ordnung");
    check(oa > 0.9 && oa < 1.4,
          "in Achsennaehe (r < R/4) faellt sie auf ERSTE Ordnung -- genau wie die "
          "Unsymmetrie h/r der Volumengewichtung es verlangt");
    check(err_E_axis.back() > err_E.back(),
          "und der Fehler ist in Achsennaehe groesser -- die Unsymmetrie waechst wie h/r");
  }

  // =========================================================================
  std::printf("\n3. Ladungserhaltung der Deposition\n");
  {
    const Box b = grounded_box(41, 81);
    std::vector<Macroparticle> parts;
    // A spread of positions: on the axis, off the axis, inside cells, exactly
    // on nodes and exactly on cell edges.
    for (int k = 0; k < 50; ++k) {
      const Real t = static_cast<Real>(k) / 49.0;
      parts.push_back({{0.9 * kR * t, -0.8 * kL + 1.6 * kL * t}, 1.6e-19 * (k + 1)});
    }
    parts.push_back({{0.0, 0.0}, 3.2e-19});                       // on the axis
    parts.push_back({b.mesh.at(7, 11), -5.0e-19});                // exactly on a node
    const DepositionResult d = deposit(b.mesh, parts);
    std::printf("    %zu Teilchen, Summe %.9e C, deponiert %.9e C\n", parts.size(),
                d.total_particles, d.total_deposited);
    check_rel(d.total_deposited, d.total_particles, 1.0e-14,
              "die deponierte Gesamtladung ist die Summe der Makropartikelladungen");
    check(d.partition_of_unity_error < 1.0e-15,
          "die Formfunktionen bilden eine Zerlegung der Eins -- das ist der Grund");
    check(d.n_outside == 0, "kein Teilchen fiel aus dem Netz");
    check(d.n_on_axis >= 1, "mindestens eines sass auf der Achse und wurde deponiert");

    // A particle OUTSIDE is counted and not deposited.
    std::vector<Macroparticle> out = {{{2.0 * kR, 0.0}, 1.0e-18}};
    const DepositionResult od = deposit(b.mesh, out);
    check(od.n_outside == 1 && od.total_deposited == 0.0,
          "ein Teilchen ausserhalb wird GEZAEHLT und nicht deponiert");
    check(od.conservation_error == 1.0,
          "und die Erhaltungspruefung schlaegt dann aus, statt zu schweigen");
  }

  // =========================================================================
  std::printf("\n4. Keine Divergenz beim Annaehern an ein Makropartikel\n");
  {
    const Box b = grounded_box(41, 81);
    const Vec2 xp{0.35 * kR, 0.0};
    const std::vector<Macroparticle> one = {{xp, 1.0e-15}};
    const DepositionResult d = deposit(b.mesh, one);
    const SpaceChargeSolution s = solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed,
                                                         b.fixed_value, d.node_charge, {});
    // Sample the potential and the field along a line through the particle at
    // ever smaller distance.  A ring singularity would grow without bound.
    std::printf("    d/R          phi [V]        |E| [V/m]\n");
    Real last_phi = 0.0, last_E = 0.0;
    bool bounded = true;
    for (int k = 0; k < 8; ++k) {
      const Real dr = 0.3 * kR * std::pow(0.5, k);
      const Vec2 x{xp.r + dr, xp.z};
      const Real ph = potential_at(b.mesh, s.phi, x);
      const Real Em = norm(interpolated_field(b.mesh, s.phi, b.eps_r, b.active, x));
      std::printf("    %-12.3e %-14.6e %-14.6e\n", dr / kR, ph, Em);
      if (k > 0 && (std::abs(ph) > 10.0 * std::abs(last_phi) + 1.0 ||
                    Em > 10.0 * last_E + 1.0))
        bounded = false;
      last_phi = ph;
      last_E = Em;
    }
    check(bounded,
          "Potential und Feld bleiben beim Annaehern beschraenkt -- die FEM-Loesung eines "
          "Knotenlastvektors ist stueckweise bilinear und hat gar keine Singularitaet");
    // And the peak is a MESH property, not a distance property: it grows under
    // refinement.  Stated, because it is the price of the method.
    std::vector<Real> peak;
    for (Index n : {21, 41, 81}) {
      Box c = grounded_box(n, 2 * n - 1);
      const DepositionResult dd = deposit(c.mesh, one);
      const SpaceChargeSolution ss = solve_with_space_charge(
          c.mesh, c.eps_r, c.active, c.fixed, c.fixed_value, dd.node_charge, {});
      Real m = 0.0;
      for (Real v : ss.phi) m = std::max(m, std::abs(v));
      peak.push_back(m);
    }
    std::printf("    Spitzenpotential ueber die Netzstufen: ");
    for (Real v : peak) std::printf("%.4e ", v);
    std::printf("V\n");
    check(peak[2] > peak[0],
          "das Selbstfeld eines Makropartikels waechst unter Verfeinerung -- es ist eine "
          "Netzgroesse, und eine PIC-Schleife darf ein Teilchen sein eigenes Feld nicht "
          "ungefiltert fuehlen lassen");
  }

  // =========================================================================
  std::printf("\n5. Selbstfeld und Elektrodenfeld getrennt, und ihre Ueberlagerung\n");
  {
    const Box b0 = grounded_box(41, 81, 0.0);      // charge only
    const Box bV = grounded_box(41, 81, 250.0);    // electrodes only
    const std::vector<Macroparticle> parts = {{{0.3 * kR, 0.1 * kL}, 2.0e-15},
                                              {{0.0, -0.2 * kL}, -1.0e-15}};
    const DepositionResult d = deposit(b0.mesh, parts);

    const SpaceChargeSolution self = solve_with_space_charge(
        b0.mesh, b0.eps_r, b0.active, b0.fixed, b0.fixed_value, d.node_charge, {});
    const SpaceChargeSolution elec = solve_with_space_charge(
        bV.mesh, bV.eps_r, bV.active, bV.fixed, bV.fixed_value, {}, {});
    const SpaceChargeSolution both = solve_with_space_charge(
        bV.mesh, bV.eps_r, bV.active, bV.fixed, bV.fixed_value, d.node_charge, {});

    Real worst = 0.0, scale = 0.0;
    for (std::size_t k = 0; k < both.phi.size(); ++k) {
      worst = std::max(worst, std::abs(both.phi[k] - (self.phi[k] + elec.phi[k])));
      scale = std::max(scale, std::abs(both.phi[k]));
    }
    std::printf("    |phi_beide - (phi_selbst + phi_elektroden)| = %.3e V von %.3e V\n", worst,
                scale);
    check(worst / scale < 1.0e-12,
          "das Problem ist linear: Selbstfeld und Elektrodenfeld ueberlagern sich exakt");
    check(elec.max_potential_shift == 0.0,
          "ohne Ladung ist die Verschiebung gegen die Laplace-Loesung exakt null");
    check(self.max_potential_shift > 0.0 && both.max_potential_shift > 0.0,
          "mit Ladung ist sie es nicht");
    check_rel(both.max_potential_shift, self.max_potential_shift, 1.0e-12,
              "und die Verschiebung ist in beiden Faellen dieselbe -- sie haengt nur an der "
              "Ladung, nicht an den Elektroden");
  }

  // =========================================================================
  std::printf("\n6. Polaritaet\n");
  {
    const Box b = grounded_box(41, 81);
    const Vec2 xp{0.3 * kR, 0.0};
    for (Real q : {+1.0e-15, -1.0e-15}) {
      const DepositionResult d = deposit(b.mesh, {{xp, q}});
      const SpaceChargeSolution s = solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed,
                                                           b.fixed_value, d.node_charge, {});
      const Real ph = potential_at(b.mesh, s.phi, xp);
      std::printf("    q = %+.1e C  ->  phi(Teilchen) = %+.6e V\n", q, ph);
      check((q > 0.0) == (ph > 0.0),
            "eine positive Ladung hebt das Potential, eine negative senkt es");
    }
    // And the solution is odd in the charge, exactly.
    const DepositionResult dp = deposit(b.mesh, {{xp, +1.0e-15}});
    const DepositionResult dm = deposit(b.mesh, {{xp, -1.0e-15}});
    const SpaceChargeSolution sp = solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed,
                                                          b.fixed_value, dp.node_charge, {});
    const SpaceChargeSolution sm = solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed,
                                                          b.fixed_value, dm.node_charge, {});
    Real worst = 0.0, scale = 0.0;
    for (std::size_t k = 0; k < sp.phi.size(); ++k) {
      worst = std::max(worst, std::abs(sp.phi[k] + sm.phi[k]));
      scale = std::max(scale, std::abs(sp.phi[k]));
    }
    check(worst / scale < 1.0e-12, "die Loesung ist exakt ungerade in der Ladung");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
