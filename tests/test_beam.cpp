#include <cmath>
#include <cstdio>

#include "es/beam.hpp"
#include "es/constants.hpp"
#include "es/fluid.hpp"
#include "es/status.hpp"
#include <string>

using namespace es;
using constants::eps0;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-46s got=%-15.7g want=%-15.7g relerr=%-9.2g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-46s %s\n", what, ok ? "OK" : "FAIL");
}

// ---------------------------------------------------------------------------
static void test_ring_kernel() {
  std::printf("\n=== ring macroparticle kernel ===\n");
  const Real Q = 1e-15;
  const Vec2 src{2.0e-4, 1.0e-4};

  // Far field: a ring looks like a point charge.
  for (Real d : {0.05, 0.2}) {
    const Vec2 x{0.0, src.z + d};
    char buf[96];
    std::snprintf(buf, sizeof buf, "V at %.2f m looks like a point charge", d);
    check(buf, ring_potential(Q, x, src), Q / (4.0 * pi * eps0 * d), 1e-5);
  }

  // Field must be minus the gradient of the potential.
  for (const Vec2& x : {Vec2{3.0e-4, 3.0e-4}, Vec2{5.0e-5, -2.0e-4}, Vec2{1.0e-9, 5.0e-4}}) {
    const Vec2 E = ring_field(Q, x, src);
    const Real hr = 1e-7 * 2.0e-4, hz = 1e-7 * 2.0e-4;
    const Real Er = -(ring_potential(Q, {x.r + hr, x.z}, src) -
                      ring_potential(Q, {x.r - hr, x.z}, src)) / (2 * hr);
    const Real Ez = -(ring_potential(Q, {x.r, x.z + hz}, src) -
                      ring_potential(Q, {x.r, x.z - hz}, src)) / (2 * hz);
    const Real scale = std::max(std::abs(Er), std::abs(Ez));
    char buf[96];
    std::snprintf(buf, sizeof buf, "E = -grad V at (%.1e,%.1e)", x.r, x.z);
    const Real err = std::max(std::abs(E.r - Er), std::abs(E.z - Ez)) / scale;
    if (err > 1e-4) ++failures;
    std::printf("  %-46s relerr=%-9.2g %s\n", buf, err, err <= 1e-4 ? "OK" : "FAIL");
  }

  // A ring degenerating onto the axis must match the point-charge limit.
  const Vec2 axis{0.0, 0.0};
  const Vec2 tiny{1.0e-9, 0.0};
  const Vec2 probe{0.0, 1.0e-3};
  check("axis ring matches the r -> 0 limit", ring_potential(Q, probe, tiny),
        ring_potential(Q, probe, axis), 1e-8);
}

// ---------------------------------------------------------------------------
// Energy conservation is the sharpest test of the whole chain: launch a
// particle from rest and its kinetic energy per charge on arrival must equal
// the potential it fell through.  This exercises the field evaluation, the
// adaptive stepping and the Verlet integrator simultaneously.
// ---------------------------------------------------------------------------
static BemSolver build_needle_extractor(Real U, Real gap, Real aperture) {
  NeedleParams np;
  np.tip_radius = 2.0e-6;
  np.half_angle = 20.0 * pi / 180.0;
  np.shank_radius = 1.0e-4;
  np.length = 8.0e-4;
  np.z_tip = 0.0;
  ExtractorParams ep;
  ep.aperture_radius = aperture;
  ep.outer_radius = 2.0e-3;
  ep.thickness = 1.0e-4;
  ep.z_plate = gap;
  ep.h_edge = aperture / 25.0;
  BemSolver bem(merge({make_needle(np), make_extractor(ep)}));
  bem.solve({U, 0.0, 0.0});
  return bem;
}

static void test_energy_conservation() {
  std::printf("\n=== energy conservation along the trajectory ===\n");
  const Real U = 1500.0;
  BemSolver bem = build_needle_extractor(U, 3.0e-4, 1.5e-4);
  std::printf("  N = %d elements, peak |E_n| on the needle = %.4g V/m\n",
              static_cast<int>(bem.size()), bem.peak_emitter_field());

  // Launch from the apex region of the needle with a uniform weight.
  std::vector<Real> w(static_cast<std::size_t>(bem.size()), 0.0);
  int launched = 0;
  for (Index i = 0; i < bem.size(); ++i) {
    const Element& e = bem.mesh().elems[static_cast<std::size_t>(i)];
    if (e.tag == Tag::Emitter && e.mid.z > -6.0e-6 && e.mid.r < 4.0e-6) {
      w[static_cast<std::size_t>(i)] = 1e-9;
      ++launched;
    }
  }
  std::printf("  launching from %d elements near the apex\n", launched);
  expect("found emitting elements near the apex", launched >= 3);

  BeamParams bp;
  bp.z_end = 1.5e-3;
  bp.r_max = 4.0e-3;
  bp.cfl = 0.03;
  bp.max_steps = 60000;
  BeamResult res = trace_beam_with_weights(bem, w, {{"ion", 5.0e5, 1.0}}, bp);
  res.print(stdout);

  Real worst = 0.0;
  int checked = 0;
  for (const Ray& r : res.rays) {
    if (r.status != RayStatus::Transmitted) continue;
    // Energy gained = potential drop between launch point and arrival point.
    const Real drop = bem.potential_at(r.x0) - bem.potential_at(r.x);
    if (!(drop > 1.0)) continue;
    worst = std::max(worst, std::abs(r.energy_eV - drop) / drop);
    ++checked;
  }
  std::printf("  checked %d transmitted rays; worst energy error = %.3g\n", checked, worst);
  expect("at least one ray was transmitted", checked > 0);
  if (checked > 0) {
    if (worst > 2e-3) ++failures;
    std::printf("  %-46s %s\n", "kinetic energy matches the potential drop",
                worst <= 2e-3 ? "OK" : "FAIL");
  }
  // Arrival energy must be close to the full applied voltage once the beam is
  // out in the field-free region beyond the extractor.
  expect("mean beam energy is within 5% of the applied voltage",
         std::abs(res.mean_energy_eV - U) < 0.05 * U);
}

// ---------------------------------------------------------------------------
static void test_interception() {
  std::printf("\n=== extractor interception vs aperture size ===\n");
  Real prev = -1.0;
  for (Real ap : {6.0e-5, 1.2e-4, 2.5e-4}) {
    BemSolver bem = build_needle_extractor(1500.0, 3.0e-4, ap);
    std::vector<Real> w(static_cast<std::size_t>(bem.size()), 0.0);
    for (Index i = 0; i < bem.size(); ++i) {
      const Element& e = bem.mesh().elems[static_cast<std::size_t>(i)];
      if (e.tag == Tag::Emitter && e.mid.z > -1.0e-5 && e.mid.r < 8.0e-6)
        w[static_cast<std::size_t>(i)] = 1e-9;
    }
    BeamParams bp;
    bp.z_end = 1.5e-3;
    bp.r_max = 4.0e-3;
    bp.cfl = 0.04;
    bp.max_steps = 40000;
    const BeamResult r = trace_beam_with_weights(bem, w, {{"ion", 5.0e5, 1.0}}, bp);
    std::printf("  aperture %6.1f um: interception %6.2f %%, half-angle(95%%) %5.2f deg\n",
                ap * 1e6, 100.0 * r.interception_fraction,
                r.half_angle_95 * 180.0 / pi);
    if (prev >= 0.0)
      expect("interception falls as the aperture opens up", r.interception_fraction <= prev + 1e-9);
    prev = r.interception_fraction;
  }
}

// ---------------------------------------------------------------------------
// Space charge is disabled: the ring model is not well posed (see beam.hpp).
// What is tested here is that it fails closed rather than producing numbers.
// ---------------------------------------------------------------------------
static void test_space_charge_is_closed() {
  std::printf("\n=== space charge: disabled, must fail closed ===\n");
  BemSolver bem = build_needle_extractor(1500.0, 3.0e-4, 2.0e-4);
  std::vector<Real> w(static_cast<std::size_t>(bem.size()), 0.0);
  for (Index i = 0; i < bem.size(); ++i) {
    const Element& e = bem.mesh().elems[static_cast<std::size_t>(i)];
    if (e.tag == Tag::Emitter && e.mid.z > -1.0e-5 && e.mid.r < 8.0e-6)
      w[static_cast<std::size_t>(i)] = 1e-9;
  }
  BeamParams bp;
  bp.z_end = 1.0e-3;
  bp.max_steps = 2000;

  // Without space charge the trace still works.
  const BeamResult a = trace_beam_with_weights(bem, w, {{"ion", 5.0e5, 1.0}}, bp);
  expect("Laplace-only trace still runs", a.current_launched > 0.0);

  // With space charge requested it must refuse, and say why.
  bp.space_charge_iters = 3;
  bool threw = false;
  std::string msg;
  try {
    (void)trace_beam_with_weights(bem, w, {{"ion", 5.0e5, 1.0}}, bp);
  } catch (const NotImplementedInThisPhase& e) {
    threw = true;
    msg = e.what();
  }
  expect("space_charge_iters > 0 throws NotImplementedInThisPhase", threw);
  expect("message names phase P4", msg.find("P4") != std::string::npos);

  // The ring kernels themselves stay, and stay singular -- which is the reason
  // the model is closed rather than regularised by hand.
  const Vec2 xp{5.0e-5, 1.0e-4};
  const Vec2 E0 = ring_field(1e-15, xp, xp);
  expect("ring self-field is still singular (the reason for the lock)",
         !std::isfinite(E0.r) || !std::isfinite(E0.z));
}

// ---------------------------------------------------------------------------
static void test_droplets_are_closed() {
  std::printf("\n=== droplet species: locked until the cone-jet coupling ===\n");
  BemSolver bem = build_needle_extractor(1500.0, 3.0e-4, 2.0e-4);
  std::vector<Real> w(static_cast<std::size_t>(bem.size()), 0.0);
  for (Index i = 0; i < bem.size(); ++i)
    if (bem.mesh().elems[static_cast<std::size_t>(i)].tag == Tag::Emitter)
      w[static_cast<std::size_t>(i)] = 1e-9;
  BeamParams bp;
  bp.z_end = 1.0e-3;
  bp.max_steps = 2000;

  bool threw = false;
  std::string msg;
  try {
    (void)trace_beam_with_weights(bem, w, {{"droplet", 1e4, 1.0, SpeciesKind::Droplet}}, bp);
  } catch (const NotImplementedInThisPhase& e) { threw = true; msg = e.what(); }
  expect("droplet species throws NotImplementedInThisPhase", threw);
  expect("message names phase P6", msg.find("P6") != std::string::npos);
}

int main() {
  test_ring_kernel();
  test_energy_conservation();
  test_interception();
  test_space_charge_is_closed();
  test_droplets_are_closed();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
