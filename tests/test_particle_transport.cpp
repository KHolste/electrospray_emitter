// tests/test_particle_transport.cpp -- P7: electrostatic particle transport
//
// The mandatory checks:
//   * a uniform field, against the closed-form parabola;
//   * the energy gain q dV;
//   * polarity reversal;
//   * time-step and mesh convergence;
//   * hit detection on all four surfaces;
//   * extracted + intercepted + still flying = launched, exactly.
//
// Everything is a TRANSPORT RESPONSE to a prescribed launch distribution.  P5
// is blocked, so there is no physical particle source and no current here is a
// prediction.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/particle_transport.hpp"
#include "es/space_charge.hpp"

using namespace es;
using constants::amu;
using constants::e;

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

constexpr Real kR = 2.0e-5, kZ = 1.0e-4;

struct Field {
  QuadMesh mesh;
  std::vector<Real> eps_r, phi;
  std::vector<char> active;
  TransportDomain domain;
};

/// A parallel-plate field: phi = V z / Z exactly, so E_z = -V/Z is uniform.
/// It is IMPOSED at every node rather than solved, so the trajectory test
/// measures the integrator and not the field solver.
Field uniform_field(Real V, Index nr = 41, Index nz = 81, Real aperture = 0.5 * kR) {
  Field f;
  f.mesh.nr = nr;
  f.mesh.nz = nz;
  f.mesh.nodes.resize(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nz));
  for (Index j = 0; j < nz; ++j)
    for (Index i = 0; i < nr; ++i)
      f.mesh.nodes[static_cast<std::size_t>(j * nr + i)] =
          Vec2{kR * static_cast<Real>(i) / static_cast<Real>(nr - 1),
               kZ * static_cast<Real>(j) / static_cast<Real>(nz - 1)};
  f.eps_r.assign(static_cast<std::size_t>(f.mesh.n_cells()), 1.0);
  f.active.assign(static_cast<std::size_t>(f.mesh.n_cells()), 1);
  f.phi.resize(static_cast<std::size_t>(f.mesh.n_nodes()));
  for (Index j = 0; j < nz; ++j)
    for (Index i = 0; i < nr; ++i)
      f.phi[static_cast<std::size_t>(f.mesh.node(i, j))] = V * f.mesh.at(i, j).z / kZ;
  f.domain = {kR, kZ, aperture};
  return f;
}

const TransportSpecies kCation{"EMI+ (Testspezies)", +e, 111.17e-3 / 6.02214076e23};
const TransportSpecies kAnion{"BF4- (Testspezies)", -e, 86.81e-3 / 6.02214076e23};

}  // namespace

int main() {
  std::printf("P7 -- elektrostatischer Teilchentransport\n\n");

  // =========================================================================
  std::printf("1. Gleichfoermiges Feld gegen die geschlossene Parabel\n");
  {
    // A cation is pulled towards LOW potential, so the extractor must be
    // negative for it to accelerate away from z = 0.
    const Real V = -1000.0;
    const Field f = uniform_field(V);
    const Real Ez = -V / kZ;             // E = -grad phi, uniform
    const Real a = kCation.charge * Ez / kCation.mass;
    std::printf("    E_z = %.4e V/m, a = %.4e m/s^2\n", Ez, a);

    std::vector<TracedParticle> launch;
    TracedParticle p;
    p.x = {0.3 * kR, 1.0e-6};
    p.v = {0.0, 0.0};
    p.current = 1.0e-9;
    launch.push_back(p);

    // The exact time to reach z = Z from rest: Z - z0 = a t^2 / 2.
    const Real t_exact = std::sqrt(2.0 * (kZ - p.x.z) / a);
    std::vector<Real> err;
    for (int n : {200, 400, 800, 1600}) {
      const Real dt = t_exact / n;
      const TransportResult r = transport_particles(f.mesh, f.phi, f.eps_r, f.active, f.domain,
                                                    kCation, launch, dt, 4 * n);
      const TracedParticle& q = r.particles[0];
      // Compare the SPEED at the exit against the closed form v = sqrt(2 a s).
      const Real v_exact = std::sqrt(2.0 * a * (q.x_end.z - p.x.z));
      err.push_back(std::abs(norm(q.v) - v_exact) / v_exact);
    }
    std::printf("    Geschwindigkeitsfehler ueber die Zeitschritte: ");
    for (Real x : err) std::printf("%.3e ", x);
    std::printf("\n");
    // IN A UNIFORM FIELD VELOCITY VERLET IS EXACT: the acceleration is
    // constant, so the update x + dt v + dt^2 a/2 IS the closed-form parabola.
    // The measured error is therefore at round-off (1e-15 to 5e-14) and NO time
    // order can be read from it -- trying to would measure the accumulation of
    // round-off, which is what a first version of this test did and why it
    // reported a negative order.  The order is measured below, in a field that
    // is not uniform.
    check(err.back() < 1.0e-12,
          "im gleichfoermigen Feld ist der Integrator exakt: der Fehler liegt auf dem "
          "Rundungsboden");

    // The radial coordinate must not move at all: E_r = 0 exactly.
    const TransportResult r = transport_particles(f.mesh, f.phi, f.eps_r, f.active, f.domain,
                                                  kCation, launch, t_exact / 800, 4000);
    check(std::abs(r.particles[0].x_end.r - r.particles[0].x_start.r) < 1.0e-12 * kR,
          "im gleichfoermigen Feld bewegt sich die radiale Koordinate nicht");
  }

  // =========================================================================
  std::printf("\n2. Energiegewinn q dV\n");
  {
    const Real V = -1000.0;
    const Field f = uniform_field(V);
    std::vector<TracedParticle> launch;
    for (Real r0 : {0.05, 0.2, 0.4}) {
      TracedParticle p;
      p.x = {r0 * kR, 1.0e-7};
      p.v = {0.0, 0.0};
      p.current = 1.0e-9;
      launch.push_back(p);
    }
    const TransportResult r = transport_particles(f.mesh, f.phi, f.eps_r, f.active, f.domain,
                                                  kCation, launch, 1.0e-12, 200000);
    for (const TracedParticle& q : r.particles) {
      const Real dV = q.potential_start - q.potential_end;
      std::printf("    r0 = %.3e: dV = %+.4f V, Energiegewinn %.6e J, q dV = %.6e J\n",
                  q.x_start.r, dV, q.energy_gain, kCation.charge * dV);
      check_rel(q.energy_gain, kCation.charge * dV, 1.0e-5,
                "der Energiegewinn ist q dV");
    }
    std::printf("    groesster Energiefehler entlang der Bahnen: %.3e\n",
                r.balance.max_energy_error);
    check(r.balance.max_energy_error < 1.0e-5,
          "die Gesamtenergie bleibt entlang der Bahn erhalten");
  }

  // =========================================================================
  std::printf("\n3. Polaritaetsumkehr\n");
  {
    // Reversing BOTH the species charge and the applied potential must leave
    // the trajectory invariant: the force q E is unchanged.  Reversing only one
    // of them must reverse the motion.
    std::vector<TracedParticle> launch;
    TracedParticle p;
    p.x = {0.3 * kR, 1.0e-6};
    p.v = {0.0, 0.0};
    p.current = 1.0e-9;
    launch.push_back(p);

    // Same mass for both so that only the SIGN differs -- otherwise the test
    // would measure the mass ratio instead of the polarity.
    const TransportSpecies pos{"positiv", +e, kCation.mass};
    const TransportSpecies neg{"negativ", -e, kCation.mass};

    const Field fm = uniform_field(-1000.0);
    const Field fp = uniform_field(+1000.0);
    const TransportResult a = transport_particles(fm.mesh, fm.phi, fm.eps_r, fm.active,
                                                  fm.domain, pos, launch, 1.0e-12, 200000);
    const TransportResult b = transport_particles(fp.mesh, fp.phi, fp.eps_r, fp.active,
                                                  fp.domain, neg, launch, 1.0e-12, 200000);
    std::printf("    Kation im negativen Feld : %s bei z = %.4e nach %.4e s\n",
                to_string(a.particles[0].fate), a.particles[0].x_end.z, a.particles[0].time);
    std::printf("    Anion  im positiven Feld : %s bei z = %.4e nach %.4e s\n",
                to_string(b.particles[0].fate), b.particles[0].x_end.z, b.particles[0].time);
    check(a.particles[0].fate == b.particles[0].fate,
          "die Umkehr BEIDER Vorzeichen laesst die Bahn unveraendert");
    check_rel(b.particles[0].time, a.particles[0].time, 1.0e-12, "gleiche Flugzeit");
    check(std::abs(b.particles[0].x_end.z - a.particles[0].x_end.z) < 1.0e-12 * kZ,
          "und gleicher Auftreffort");

    // Only the species reversed: the particle must run into the emitter.
    const TransportResult c = transport_particles(fm.mesh, fm.phi, fm.eps_r, fm.active,
                                                  fm.domain, neg, launch, 1.0e-12, 200000);
    std::printf("    Anion im negativen Feld  : %s\n", to_string(c.particles[0].fate));
    check(c.particles[0].fate == Fate::HitEmitter,
          "nur die Spezies umgekehrt: das Teilchen laeuft zurueck auf den Emitter");
  }

  // =========================================================================
  std::printf("\n4. Treffererkennung auf allen vier Flaechen\n");
  {
    const Field f = uniform_field(-1000.0, 41, 81, 0.3 * kR);
    std::vector<TracedParticle> launch;
    auto add = [&](Vec2 x, Vec2 v) {
      TracedParticle p;
      p.x = x;
      p.v = v;
      p.current = 1.0e-9;
      launch.push_back(p);
    };
    add({0.1 * kR, 1.0e-6}, {0.0, 0.0});                 // straight through the aperture
    add({0.8 * kR, 1.0e-6}, {0.0, 0.0});                 // onto the extractor
    add({0.5 * kR, 0.5 * kZ}, {5.0e4, 0.0});             // sideways onto the polymer
    add({0.2 * kR, 1.0e-7}, {0.0, -5.0e3});              // back onto the emitter
    const TransportResult r = transport_particles(f.mesh, f.phi, f.eps_r, f.active, f.domain,
                                                  kCation, launch, 1.0e-12, 400000);
    for (const TracedParticle& q : r.particles)
      std::printf("    Start (%.3e, %.3e) -> %s bei (%.3e, %.3e)\n", q.x_start.r,
                  q.x_start.z, to_string(q.fate), q.x_end.r, q.x_end.z);
    check(r.particles[0].fate == Fate::LeftThroughAperture, "durch die Blende");
    check(r.particles[1].fate == Fate::HitExtractor, "auf den Extraktor");
    check(r.particles[2].fate == Fate::HitPolymer, "auf das Polymer");
    check(r.particles[3].fate == Fate::HitEmitter, "zurueck auf den Emitter");
  }

  // =========================================================================
  std::printf("\n5. Die Bilanz schliesst\n");
  {
    const Field f = uniform_field(-1000.0, 41, 81, 0.3 * kR);
    std::vector<TracedParticle> launch;
    for (int k = 0; k < 60; ++k) {
      TracedParticle p;
      const Real t = (static_cast<Real>(k) + 0.5) / 60.0;
      p.x = {0.95 * kR * t, 1.0e-6};
      p.v = {0.0, 0.0};
      p.current = 1.0e-9 * (1.0 + 0.5 * t);   // deliberately not all equal
      launch.push_back(p);
    }
    // A short budget on purpose, so that some particles are STILL FLYING and
    // the balance has to account for them separately instead of folding them in.
    const TransportResult r = transport_particles(f.mesh, f.phi, f.eps_r, f.active, f.domain,
                                                  kCation, launch, 1.0e-12, 3000);
    r.balance.print(stdout);
    check(r.balance.closure_error < 1.0e-14,
          "extrahiert + abgefangen + noch fliegend = gestartet, exakt");
    check(r.balance.n_launched ==
              r.balance.n_extracted + r.balance.n_intercepted + r.balance.n_flying,
          "und die Teilchenzahlen ebenso");
    check(r.balance.n_flying > 0,
          "es gibt hier tatsaechlich noch fliegende Teilchen -- sie werden getrennt "
          "gezaehlt und nicht als verloren verbucht");

    const TransportResult full = transport_particles(f.mesh, f.phi, f.eps_r, f.active,
                                                     f.domain, kCation, launch, 1.0e-12,
                                                     400000);
    check(full.balance.n_flying == 0, "mit genug Schritten fliegt keines mehr");
    check(full.balance.closure_error < 1.0e-14, "und die Bilanz schliesst weiterhin");
    check_rel(full.balance.extracted + full.balance.intercepted, full.balance.launched,
              1.0e-14, "extrahiert + abgefangen = gestartet");
  }

  // =========================================================================
  std::printf("\n6. Netzkonvergenz und die unvollstaendige Spezies\n");
  {
    // The trajectory in the IMPOSED uniform field must be mesh independent,
    // because the field is exact at every node and bilinear in between.
    std::vector<TracedParticle> launch;
    TracedParticle p;
    p.x = {0.3 * kR, 1.0e-6};
    p.v = {0.0, 0.0};
    p.current = 1.0e-9;
    launch.push_back(p);
    std::vector<Real> zend;
    for (Index n : {21, 41, 81}) {
      const Field f = uniform_field(-1000.0, n, 2 * n - 1);
      const TransportResult r = transport_particles(f.mesh, f.phi, f.eps_r, f.active, f.domain,
                                                    kCation, launch, 1.0e-12, 200000);
      zend.push_back(norm(r.particles[0].v));
    }
    std::printf("    Endgeschwindigkeit ueber die Netzstufen: ");
    for (Real v : zend) std::printf("%.9e ", v);
    std::printf("\n");
    check(std::abs(zend[2] - zend[0]) / zend[0] < 1.0e-6,
          "im exakt dargestellten Feld ist die Bahn netzunabhaengig");

    // THE TIME ORDER, measured where it can be: a field that is NOT uniform.
    // The manufactured Poisson potential of P6 is quadratic in r and z, so the
    // force varies along the path and the integrator is no longer exact.  There
    // is no closed-form trajectory, so the order is measured by
    // SELF-CONVERGENCE against the finest step.
    {
      const Real phi0 = 4000.0, Rm = kR, Lm = 0.5 * kZ;
      Field g = uniform_field(0.0, 81, 161);   // reuse the mesh, replace the field
      for (Index j = 0; j < g.mesh.nz; ++j)
        for (Index i = 0; i < g.mesh.nr; ++i) {
          const Vec2 x = g.mesh.at(i, j);
          const Vec2 xs{x.r, x.z - 0.5 * kZ};
          g.phi[static_cast<std::size_t>(g.mesh.node(i, j))] =
              manufactured_potential(xs, Rm, Lm, phi0);
        }
      std::vector<TracedParticle> lg;
      TracedParticle q;
      q.x = {0.4 * kR, 0.35 * kZ};
      q.v = {0.0, 0.0};
      q.current = 1.0e-9;
      lg.push_back(q);
      const Real T = 2.0e-9;
      std::vector<Vec2> ends;
      for (int nsteps : {500, 1000, 2000, 4000, 8000}) {
        const TransportResult rr = transport_particles(g.mesh, g.phi, g.eps_r, g.active,
                                                       g.domain, kCation, lg, T / nsteps,
                                                       nsteps);
        ends.push_back(rr.particles[0].x);
      }
      std::vector<Real> de;
      for (int k = 0; k < 4; ++k) de.push_back(norm(ends[k] - ends[4]) / kR);
      std::printf("    Selbstkonvergenz im NICHT gleichfoermigen Feld: ");
      for (Real x : de) std::printf("%.3e ", x);
      const Real ord = std::log(de[0] / de[3]) / std::log(8.0);
      std::printf("\n    beobachtete Zeitordnung: %.2f\n", ord);
      check(ord > 1.7 && ord < 2.4,
            "Velocity-Verlet ist zweiter Ordnung -- gemessen dort, wo der Integrator nicht "
            "ohnehin exakt ist");
    }

    // Fail closed on an incomplete species.
    bool threw = false;
    try {
      const Field f = uniform_field(-1000.0);
      TransportSpecies bad{"ohne Ladung", 0.0, 1.0e-25};
      (void)transport_particles(f.mesh, f.phi, f.eps_r, f.active, f.domain, bad, launch,
                                1.0e-12, 10);
    } catch (const std::exception&) {
      threw = true;
    }
    check(threw, "eine Spezies ohne Ladung wird abgelehnt, nicht mit einer Vorgabe ersetzt");
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
