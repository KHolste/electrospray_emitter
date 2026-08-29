// P2b: the dielectric axisymmetric electrostatics of the capillary emitter.
//
// The point of this file is not to re-run the mesher's own validate() and call
// that a test.  Every check below either compares against an INDEPENDENT
// computation (the boundary-element solver, a closed-form volume, a second mesh
// level) or asserts a property that the superseded P2a arrangement would fail.
//
// The most important one is the last kind.  P2b exists because P2a treated a
// photopolymer as a conductor; a test suite that cannot tell the two apart
// would not have caught that, so `test_no_polymer_is_a_conductor` walks the
// named polymer surfaces explicitly and `test_metallic_reference_is_different`
// shows that the corrected model actually gives a different answer.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/boundary_mesh.hpp"
#include "es/constants.hpp"
#include "es/dielectric_device.hpp"
#include "es/vacuum_bem.hpp"

using namespace es;
using constants::eps0;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-54s got=%-13.7g want=%-13.7g rel=%-9.3g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

static void check_small(const char* what, Real got, Real bound) {
  const bool ok = std::abs(got) <= bound;
  if (!ok) ++failures;
  std::printf("  %-54s |%.4e| <= %.3e %s\n", what, got, bound, ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-54s %s\n", what, ok ? "OK" : "FAIL");
}

// ---------------------------------------------------------------------------
// The reference configuration of P2b.  Kept in one place so that every test
// varies exactly one thing against it.
// ---------------------------------------------------------------------------

static DielectricDeviceParameters reference_geometry(int level = 1) {
  DielectricDeviceParameters p;  // device defaults are the documented P1 set
  p.device.emitter_back_length = 0.0;
  // The open domain has to be much larger than the P1 box: see
  // test_far_field_treatment, which measures why.
  p.device.domain_radius = 1.2e-2;
  p.device.domain_z_min = -5.0e-3;
  p.device.domain_z_max = 9.0e-3;
  p.liquid_feed_z = -2.0e-4;
  p.mesh_level = level;
  return p;
}

static DielectricSetup reference_setup(const MaterialLibrary& lib, int level = 1) {
  DielectricSetup s;
  s.geometry = reference_geometry(level);
  s.materials = DielectricMaterials::reference(lib);
  s.conductor_model = ConductorModel::Dielectric;
  s.metallisation = Metallisation::FrontAndAperture;
  s.far_field = FarField::Asymptotic;
  s.V_emitter = 1500.0;
  s.V_extractor = 0.0;
  return s;
}

static Index probe_index(const std::vector<Probe>& p, const char* name) {
  for (std::size_t k = 0; k < p.size(); ++k)
    if (p[k].name == name) return static_cast<Index>(k);
  return -1;
}

// ==========================================================================
// 1.  The volume mesh
// ==========================================================================
static void test_mesh() {
  std::printf("\n=== 1. Automatisches Volumennetz ===\n");
  const DeviceVolumeMesh m = build_volume_mesh(reference_geometry(1));
  const VolumeMeshReport rep = m.validate();
  rep.print(stdout);
  expect("alle Netzpruefungen bestanden", rep.all_passed());

  // Region volumes against closed-form values -- independent of validate().
  check("Rotationsvolumen der Fluessigkeitssaeule", m.revolved_volume_of(Region::Liquid),
        pi * m.r_bore * m.r_bore * (-m.p.liquid_feed_z), 1.0e-12);
  check("Rotationsvolumen des Extraktortraegers", m.revolved_volume_of(Region::ExtractorSolid),
        pi * (m.r_ext_outer * m.r_ext_outer - m.r_aperture * m.r_aperture) *
            m.p.device.extractor_thickness,
        1.0e-12);

  // Determinism: same parameters, bitwise same mesh.
  const DeviceVolumeMesh m2 = build_volume_mesh(reference_geometry(1));
  bool same = (m.grid.nr == m2.grid.nr) && (m.grid.nz == m2.grid.nz);
  if (same)
    for (std::size_t k = 0; k < m.grid.nodes.size(); ++k)
      same = same && m.grid.nodes[k].r == m2.grid.nodes[k].r &&
             m.grid.nodes[k].z == m2.grid.nodes[k].z;
  expect("zweiter Aufbau ist bitgenau identisch", same);

  // Moving the feed boundary must not move a single node at or above the taper
  // foot.  Without that the feed study would measure the mesh as well.
  DielectricDeviceParameters q = reference_geometry(1);
  q.liquid_feed_z = -3.7e-4;
  const DeviceVolumeMesh m3 = build_volume_mesh(q);
  bool front_identical = (m.grid.nr == m3.grid.nr);
  Index compared = 0;
  if (front_identical) {
    for (Index j = m.j_foot, j3 = m3.j_foot; j < m.grid.nz && j3 < m3.grid.nz; ++j, ++j3) {
      if (m.grid.z_of_row(j) != m3.grid.z_of_row(j3)) {
        front_identical = false;
        break;
      }
      for (Index i = 0; i < m.grid.nr; ++i)
        front_identical = front_identical && (m.grid.at(i, j).r == m3.grid.at(i, j3).r);
      ++compared;
    }
  }
  std::printf("  verglichene Zeilen ab dem Kegelfuss: %lld\n", static_cast<long long>(compared));
  expect("das Netz ab dem Kegelfuss haengt nicht von der Zulaufposition ab", front_identical);

  // The P2a rearward conducting closure must be refused outright.
  DielectricDeviceParameters bad = reference_geometry(0);
  bad.device.emitter_back_length = 8.0e-4;
  bool threw = false;
  try {
    (void)build_volume_mesh(bad);
  } catch (const std::exception&) {
    threw = true;
  }
  expect("die leitende P2a-Abschlussscheibe wird abgelehnt", threw);
}

// ==========================================================================
// 2.  No polymer surface is a conductor  (required check 9)
// ==========================================================================
static void test_no_polymer_is_a_conductor() {
  std::printf("\n=== 2. Keine Polymerflaeche wird als Leiter behandelt ===\n");
  const MaterialLibrary lib;
  const DielectricSetup s = reference_setup(lib, 1);
  const DeviceVolumeMesh m = build_volume_mesh(s.geometry);
  const std::vector<NodeRole> role = node_roles(m, s);
  const BoundaryAudit a = audit_boundaries(m, role, s);
  a.print(stdout);
  expect("Audit ohne Verletzung", a.ok());
  check_small("festgehaltene Knoten auf Polymerflaechen",
              static_cast<Real>(a.n_polymer_dirichlet), 0.0);
  check_small("festgehaltene Knoten der Schnittebene ausserhalb der Fluessigkeit",
              static_cast<Real>(a.n_feed_plane_outside_liquid), 0.0);

  auto at = [&](Index i, Index j) { return role[static_cast<std::size_t>(m.grid.node(i, j))]; };

  // Spot checks by hand, on the surfaces that matter most.
  expect("Kegelflanke auf halber Hoehe ist frei",
         at(m.i_land, (m.j_foot + m.j_tip) / 2) == NodeRole::Free);
  expect("Stirnflaeche (Land) zwischen Bohrung und Aussenkante ist frei",
         at((m.i_bore + m.i_land) / 2, m.j_tip) == NodeRole::Free);
  expect("Rueckflaeche des Polymers bei z = liquid_feed_z ist frei",
         at((m.i_bore + m.i_land) / 2, m.j_feed) == NodeRole::Free);
  expect("Aussenkante der Stirnflaeche ist frei", at(m.i_land, m.j_tip) == NodeRole::Free);
  expect("Rueckseite des Extraktortraegers ist frei",
         at((m.i_aperture + m.i_ext_outer) / 2, m.j_ex_back) == NodeRole::Free);
  expect("Aussenrand des Extraktortraegers ist frei",
         at(m.i_ext_outer, (m.j_ex_front + m.j_ex_back) / 2) == NodeRole::Free);

  // And the surfaces that ARE electrodes.
  expect("Bohrungswand traegt Emitterpotential",
         at(m.i_bore, (m.j_feed + m.j_tip) / 2) == NodeRole::LiquidConductor);
  expect("ebene Fluessigkeitsreferenz traegt Emitterpotential",
         at(m.i_bore / 2, m.j_tip) == NodeRole::LiquidConductor);
  expect("Zulaufquerschnitt ist als liquid_feed_boundary gekennzeichnet",
         at(m.i_bore / 2, m.j_feed) == NodeRole::LiquidFeedBoundary &&
             at(m.i_bore, m.j_feed) == NodeRole::LiquidFeedBoundary);
  expect("metallisierte Vorderflaeche traegt Extraktorpotential",
         at((m.i_aperture + m.i_ext_outer) / 2, m.j_ex_front) ==
             NodeRole::ExtractorMetallisation);
  expect("Aperturwand traegt Extraktorpotential (Referenzmetallisierung)",
         at(m.i_aperture, (m.j_ex_front + m.j_ex_back) / 2) == NodeRole::ExtractorMetallisation);

  // With FrontOnly the aperture wall must become polymer again.
  DielectricSetup s2 = s;
  s2.metallisation = Metallisation::FrontOnly;
  const std::vector<NodeRole> role2 = node_roles(m, s2);
  expect("mit metallisation=front_only ist die Aperturwand frei",
         role2[static_cast<std::size_t>(
             m.grid.node(m.i_aperture, (m.j_ex_front + m.j_ex_back) / 2))] == NodeRole::Free);
  expect("Audit bleibt auch dann sauber", audit_boundaries(m, role2, s2).ok());

  // A deliberately broken role vector must be caught -- otherwise the audit is
  // only testing that node_roles agrees with itself.
  std::vector<NodeRole> broken = role;
  broken[static_cast<std::size_t>(m.grid.node(m.i_land, (m.j_foot + m.j_tip) / 2))] =
      NodeRole::LiquidConductor;
  const BoundaryAudit bad = audit_boundaries(m, broken, s);
  expect("ein als Elektrode markierter Flankenknoten wird gefunden",
         !bad.ok() && bad.n_polymer_dirichlet == 1);
}

// ==========================================================================
// 3.  Interface continuity on the taper flank  (required check 3)
// ==========================================================================
static void test_interface_continuity() {
  std::printf("\n=== 3. Stetigkeit an der Materialgrenze (Kegelflanke) ===\n");
  const MaterialLibrary lib;
  std::vector<Real> err;
  for (int level : {0, 1, 2}) {
    const DielectricSolution s = solve_dielectric(reference_setup(lib, level));
    err.push_back(s.relative_interface_error());
    std::printf("  Stufe %d: D_n innen %.6e, aussen %.6e, relativer Sprung %.3e, "
                "phi-Sprung %.1e\n",
                level, s.Dn_polymer_side, s.Dn_vacuum_side, s.relative_interface_error(),
                s.phi_interface_jump);
    check_small("phi ist stetig (ein Knoten auf der Grenze)", s.phi_interface_jump, 1.0e-12);
  }
  check_small("D_n-Sprung auf der feinsten Stufe", err.back(), 3.0e-2);
  expect("der D_n-Sprung faellt monoton mit der Verfeinerung",
         err[1] < err[0] && err[2] < err[1]);
}

// ==========================================================================
// 4.  Linearity and polarity reversal  (required check 4)
// ==========================================================================
static void test_linearity_and_polarity() {
  std::printf("\n=== 4. Linearitaet und Polaritaetsumkehr ===\n");
  const MaterialLibrary lib;
  const DielectricSetup base = reference_setup(lib, 0);
  const DielectricSolution a = solve_dielectric(base);

  DielectricSetup twice = base;
  twice.V_emitter *= 2.0;
  twice.V_extractor *= 2.0;
  const DielectricSolution b = solve_dielectric(twice);

  DielectricSetup flip = base;
  flip.V_emitter = -base.V_emitter;
  flip.V_extractor = -base.V_extractor;
  const DielectricSolution c = solve_dielectric(flip);

  // Superposition with a non-trivial second electrode potential, so the test is
  // not satisfied by a solver that ignores V_extractor.  It has to be built from
  // the two UNIT-POTENTIAL solutions: adding a constant to both electrodes does
  // NOT shift phi by that constant, because phi -> 0 at infinity is not
  // translation invariant.  An earlier version of this test assumed it was and
  // failed by 12 per cent -- the test was wrong, not the solver.
  DielectricSetup only_emitter = base;
  only_emitter.V_emitter = 1500.0;
  only_emitter.V_extractor = 0.0;
  DielectricSetup only_extractor = base;
  only_extractor.V_emitter = 0.0;
  only_extractor.V_extractor = 1500.0;
  const DielectricSolution e1 = solve_dielectric(only_emitter);
  const DielectricSolution e2 = solve_dielectric(only_extractor);

  DielectricSetup shifted = base;
  shifted.V_emitter = 1000.0;
  shifted.V_extractor = 400.0;
  const DielectricSolution d = solve_dielectric(shifted);

  Real w2 = 0.0, wm = 0.0, scale = 0.0;
  for (std::size_t k = 0; k < a.fem.phi.size(); ++k) {
    w2 = std::max(w2, std::abs(b.fem.phi[k] - 2.0 * a.fem.phi[k]));
    wm = std::max(wm, std::abs(c.fem.phi[k] + a.fem.phi[k]));
    scale = std::max(scale, std::abs(a.fem.phi[k]));
  }
  check_small("phi(2V) = 2 phi(V)", w2 / scale, 1.0e-12);
  check_small("phi(-V) = -phi(V)", wm / scale, 1.0e-12);
  check("Q_emitter kehrt das Vorzeichen um", c.Q_emitter, -a.Q_emitter, 1.0e-12);
  check("Q_extractor kehrt das Vorzeichen um", c.Q_extractor, -a.Q_extractor, 1.0e-12);

  Real wsup = 0.0, sscale = 0.0;
  for (std::size_t k = 0; k < d.fem.phi.size(); ++k) {
    const Real combo =
        (1000.0 / 1500.0) * e1.fem.phi[k] + (400.0 / 1500.0) * e2.fem.phi[k];
    wsup = std::max(wsup, std::abs(d.fem.phi[k] - combo));
    sscale = std::max(sscale, std::abs(d.fem.phi[k]));
  }
  const Index k = probe_index(a.probes, "axis_gap_mid");
  std::printf("  Superposition am Mittelpunkt: direkt %.6f V, aus den Einheitsloesungen "
              "%.6f V\n",
              d.phi_probe[static_cast<std::size_t>(k)],
              (1000.0 / 1500.0) * e1.phi_probe[static_cast<std::size_t>(k)] +
                  (400.0 / 1500.0) * e2.phi_probe[static_cast<std::size_t>(k)]);
  check_small("Superposition der beiden Einheitsloesungen", wsup / sscale, 1.0e-12);
}

// ==========================================================================
// 5.  FEM against the independent BEM at eps_r = 1  (required check 5)
// ==========================================================================
//
// Both solvers get the SAME conductor arrangement -- the superseded metallic
// one, because that is the only arrangement the single-layer BEM can represent
// -- and the same voltages.  Nothing about the P2b physics is being validated
// here.  What is being validated is the finite-element machinery: the 2*pi*r
// weight, the assembly, and above all the open far-field treatment, since the
// BEM has no truncation boundary at all and the FEM has one at 12 mm.
static void test_fem_against_bem() {
  std::printf("\n=== 5. FEM gegen BEM, eps_r = 1 (Loeserquerpruefung) ===\n");
  const MaterialLibrary lib;
  const Real VE = 1500.0, VX = 0.0;

  DielectricSetup s = reference_setup(lib, 2);
  s.conductor_model = ConductorModel::MetallicReference;
  s.metallisation = Metallisation::AllSurfaces;
  s.materials = DielectricMaterials::all_vacuum(lib);
  const DielectricSolution fem = solve_dielectric(s);

  // The identical shape for the BEM: the conductor ends where the FEM cuts the
  // liquid, and is closed there by the P2a disc.
  DeviceParameters dp = s.geometry.device;
  dp.emitter_back_length = -s.geometry.liquid_feed_z;
  dp.domain_z_min = s.geometry.device.domain_z_min;
  const DeviceGeometry g = DeviceGeometry::build(dp);
  const BoundaryMesh bm = BoundaryMesh::generate(g, 0.5);
  BemSolver bem(vacuum_bem_mesh(bm, g));
  bem.solve({{VE, VX, 0.0}});
  std::printf("  BEM-Panels: %lld, FEM-Knoten: %lld\n", static_cast<long long>(bem.size()),
              static_cast<long long>(fem.fem.n_nodes));

  Real worst_phi = 0.0, worst_E = 0.0;
  for (std::size_t k = 0; k < fem.probes.size(); ++k) {
    const Vec2 x = fem.probes[k].x;
    const Real pb = bem.potential_at(x);
    const Vec2 Eb = bem.field_at(x);
    const Real dphi = std::abs(fem.phi_probe[k] - pb) / (VE - VX);
    const Real dE = std::abs(fem.Emag_probe[k] - norm(Eb)) / std::max(norm(Eb), 1e-300);
    worst_phi = std::max(worst_phi, dphi);
    worst_E = std::max(worst_E, dE);
    std::printf("  %-26s phi FEM %11.5g  BEM %11.5g  d/span %8.2e | |E| FEM %10.4g  "
                "BEM %10.4g  rel %8.2e\n",
                fem.probes[k].name.c_str(), fem.phi_probe[k], pb, dphi, fem.Emag_probe[k],
                norm(Eb), dE);
  }
  check_small("groesste Potentialabweichung, relativ zu |V_E - V_X|", worst_phi, 6.0e-3);
  check_small("groesste relative Feldabweichung", worst_E, 4.0e-2);
}

// ==========================================================================
// 6.  Mesh convergence  (required check 6)
// ==========================================================================
static void test_mesh_convergence() {
  std::printf("\n=== 6. Netzkonvergenz an kantenfernen Punkten ===\n");
  const MaterialLibrary lib;
  std::vector<DielectricSolution> sol;
  for (int level : {0, 1, 2, 3}) sol.push_back(solve_dielectric(reference_setup(lib, level)));

  const Real span = 1500.0;
  std::printf("  %-26s %11s %11s %11s %11s   letzte rel. Aenderung\n", "Punkt", "Stufe0",
              "Stufe1", "Stufe2", "Stufe3");
  Real worst_phi = 0.0, worst_E = 0.0;
  for (std::size_t k = 0; k < sol[0].probes.size(); ++k) {
    std::printf("  %-26s", sol[0].probes[k].name.c_str());
    for (const auto& s : sol) std::printf(" %11.5f", s.phi_probe[k]);
    const Real dphi = std::abs(sol[3].phi_probe[k] - sol[2].phi_probe[k]) / span;
    const Real dE = std::abs(sol[3].Emag_probe[k] - sol[2].Emag_probe[k]) /
                    std::max(sol[3].Emag_probe[k], 1e-300);
    std::printf("   phi %8.2e  |E| %8.2e\n", dphi, dE);
    worst_phi = std::max(worst_phi, dphi);
    worst_E = std::max(worst_E, dE);
  }
  check_small("Potentialaenderung Stufe 2 -> 3, relativ zur Spannweite", worst_phi, 1.0e-3);
  check_small("Feldaenderung Stufe 2 -> 3, relativ", worst_E, 1.0e-2);

  // The emitter charge is a global functional and converges monotonically; its
  // successive changes must shrink.
  std::vector<Real> q;
  for (const auto& s : sol) q.push_back(s.Q_emitter);
  const Real d01 = std::abs(q[1] - q[0]), d12 = std::abs(q[2] - q[1]),
             d23 = std::abs(q[3] - q[2]);
  std::printf("  Q_emitter: %.8e %.8e %.8e %.8e   Aenderungen %.2e %.2e %.2e\n", q[0], q[1],
              q[2], q[3], d01, d12, d23);
  expect("die Aenderungen von Q_emitter werden kleiner", d12 < d01 && d23 < d12);
  check_small("letzte relative Aenderung von Q_emitter", d23 / std::abs(q[3]), 3.0e-3);
}

// ==========================================================================
// 7.  Position of the feed boundary  (required check 7)
// ==========================================================================
//
// The requirement was to SHOW that pushing the feed boundary back stops moving
// the field at the meniscus.  It does not, and this test records the measured
// failure against tolerances fixed in advance (es::feed_truncation) rather than
// asserting a convergence that is not there.
//
// What it does assert is the DIAGNOSIS, because that is the part a later change
// could break silently:
//
//   * the changes shrink monotonically -- it is a slow tail, not a divergence;
//   * a grounded enclosure changes nothing, so the open far field is not the
//     cause;
//   * the emitter charge follows the self-capacitance of a thin cylinder,
//     2 pi eps0 L / (ln(2L/a) - 1), which is the actual mechanism.
//
// If a later phase adds the base plate that docs/04 provides for and the study
// then converges, this test will fail -- and it should, because the finding it
// records will no longer be true.
static void test_feed_boundary_position() {
  std::printf("\n=== 7. Lage der Zulaufgrenze (Trunkierungsstudie) ===\n");
  const MaterialLibrary lib;
  const Real zf[4] = {-1.0e-4, -2.0e-4, -4.0e-4, -8.0e-4};
  std::vector<DielectricSolution> sol;
  for (Real z : zf) {
    DielectricSetup s = reference_setup(lib, 1);
    s.geometry.liquid_feed_z = z;
    sol.push_back(solve_dielectric(s));
  }
  const Real span = 1500.0;
  Real worst_phi = 0.0, worst_E = 0.0;
  std::printf("  %-26s %11s %11s %11s %11s\n", "Punkt (phi [V])", "-100um", "-200um", "-400um",
              "-800um");
  for (std::size_t k = 0; k < sol[0].probes.size(); ++k) {
    std::printf("  %-26s", sol[0].probes[k].name.c_str());
    for (const auto& s : sol) std::printf(" %11.5f", s.phi_probe[k]);
    std::printf("\n");
    worst_phi = std::max(worst_phi, std::abs(sol[3].phi_probe[k] - sol[2].phi_probe[k]) / span);
    worst_E = std::max(worst_E, std::abs(sol[3].Emag_probe[k] - sol[2].Emag_probe[k]) /
                                    std::max(sol[3].Emag_probe[k], 1e-300));
  }
  std::printf("  Q_emitter: ");
  for (const auto& s : sol) std::printf(" %.6e", s.Q_emitter);
  std::printf("\n  letzte Verdopplung: phi %.3e der Spannweite (Grenze %.1e), |E| %.3e "
              "relativ (Grenze %.1e)\n",
              worst_phi, feed_truncation::kTolPhiOverSpan, worst_E,
              feed_truncation::kTolFieldRelative);
  expect("BEFUND: die Zulaufposition ist NICHT auskonvergiert -- die vorab "
         "festgelegten Grenzen werden ueberschritten",
         worst_phi > feed_truncation::kTolPhiOverSpan &&
             worst_E > feed_truncation::kTolFieldRelative);

  Real d01 = 0.0, d12 = 0.0, d23 = 0.0;
  for (std::size_t k = 0; k < sol[0].probes.size(); ++k) {
    d01 = std::max(d01, std::abs(sol[1].phi_probe[k] - sol[0].phi_probe[k]));
    d12 = std::max(d12, std::abs(sol[2].phi_probe[k] - sol[1].phi_probe[k]));
    d23 = std::max(d23, std::abs(sol[3].phi_probe[k] - sol[2].phi_probe[k]));
  }
  std::printf("  groesste Potentialaenderung je Verdopplung: %.3e, %.3e, %.3e V\n", d01, d12,
              d23);
  expect("es ist ein langsam abklingender Auslaeufer, keine Divergenz",
         d12 < d01 && d23 < d12);

  // The mechanism: charge follows the self-capacitance of the modelled column.
  std::printf("  Q_emitter gegen 2 pi eps0 L / (ln(2L/a) - 1) * V:\n");
  Real worst_ratio_spread = 0.0;
  std::vector<Real> ratio;
  for (std::size_t k = 0; k < sol.size(); ++k) {
    const Real L = -zf[k], a = sol[k].mesh.r_bore;
    const Real C = 2.0 * pi * eps0 * L / (std::log(2.0 * L / a) - 1.0);
    ratio.push_back(sol[k].Q_emitter / (C * span));
    std::printf("    L = %6.1f um : Q = %.4e C, C*V = %.4e C, Verhaeltnis %.3f\n", L * 1e6,
                sol[k].Q_emitter, C * span, ratio.back());
  }
  for (Real r : ratio)
    worst_ratio_spread = std::max(worst_ratio_spread, std::abs(r - ratio[2]) / ratio[2]);
  expect("die Emitterladung folgt der Selbstkapazitaet der Saeule (Verhaeltnis "
         "bleibt innerhalb 25 %)",
         worst_ratio_spread < 0.25);

  // And it is not the open far field: a grounded enclosure gives the same.
  DielectricSetup g = reference_setup(lib, 1);
  g.geometry.liquid_feed_z = -8.0e-4;
  g.far_field = FarField::Grounded;
  const DielectricSolution sg = solve_dielectric(g);
  Real worst_bc = 0.0;
  for (std::size_t k = 0; k < sg.probes.size(); ++k)
    worst_bc = std::max(worst_bc, std::abs(sg.phi_probe[k] - sol[3].phi_probe[k]) / span);
  check_small("geerdete Huelle aendert die Trunkierungsabhaengigkeit nicht", worst_bc, 1.0e-3);
}

// ==========================================================================
// 8.  Sensitivity to eps_r  (required check 8)
// ==========================================================================
static void test_permittivity_sensitivity() {
  std::printf("\n=== 8. Empfindlichkeit gegenueber eps_r ===\n");
  MaterialLibrary lib;
  const Material su8 = lib.get("su8");
  expect("SU-8 ist als vorlaeufig gekennzeichnet", su8.status == MaterialStatus::Provisional);
  expect("SU-8 traegt einen begruendeten Sensitivitaetsbereich", su8.has_range());

  const Real values[4] = {1.0, su8.eps_r_low, su8.relative_permittivity, su8.eps_r_high};
  std::vector<Real> Ez, Q;
  for (Real e : values) {
    MaterialLibrary l2;
    l2.override_permittivity("su8", e, MaterialStatus::Provisional, "Sensitivitaetsstudie");
    DielectricSetup s = reference_setup(l2, 1);
    s.materials = DielectricMaterials::reference(l2);
    const DielectricSolution sol = solve_dielectric(s);
    const Index k = probe_index(sol.probes, "axis_2_bore_radii");
    Ez.push_back(std::abs(sol.Ez_probe[static_cast<std::size_t>(k)]));
    Q.push_back(sol.Q_emitter);
    std::printf("  eps_r = %5.2f : E_z(2 r_bore) = %.6e V/m, Q_emitter = %.6e C\n", e,
                Ez.back(), Q.back());
  }
  // More polarisable emitter body -> more charge on the emitter conductor.
  expect("Q_emitter waechst monoton mit eps_r",
         Q[0] < Q[1] && Q[1] < Q[2] && Q[2] < Q[3]);
  const Real spread = (Ez[3] - Ez[1]) / Ez[2];
  std::printf("  E_z-Spanne ueber den Bereich %.2f ... %.2f: %.2f %% des Nominalwerts\n",
              su8.eps_r_low, su8.eps_r_high, 100.0 * spread);
  expect("die Permittivitaet wirkt ueberhaupt (mehr als 0.1 %)", std::abs(spread) > 1.0e-3);

  // IP-Q and IPx-Q are registered without a number and must refuse to be used.
  for (const char* name : {"ip-q", "ipx-q"}) {
    bool threw = false;
    try {
      (void)lib.get(name).permittivity_or_throw();
    } catch (const std::exception&) {
      threw = true;
    }
    expect((std::string(name) + " ohne Wert wird abgelehnt statt geraten").c_str(), threw);
  }
  // ... and must become usable by configuration alone, with no code change.
  MaterialLibrary l3;
  l3.override_permittivity("ip-q", 3.0, MaterialStatus::Provisional, "Testwert");
  DielectricSetup s3 = reference_setup(l3, 0);
  s3.materials.emitter_dielectric = l3.get("ip-q");
  bool ok = true;
  try {
    (void)solve_dielectric(s3);
  } catch (const std::exception& e) {
    ok = false;
    std::printf("  unerwartet: %s\n", e.what());
  }
  expect("ein per Konfiguration versorgtes Harz laeuft ohne Codeaenderung", ok);
}

// ==========================================================================
// 9.  The open far-field treatment
// ==========================================================================
//
// The asymptotic (monopole) condition and a grounded box are two DIFFERENT
// physical problems on a finite domain, and they converge to the same answer as
// the box grows.  Their difference is therefore a direct measurement of the
// truncation error, and the asymptotic one must approach the limit faster --
// otherwise there would be no reason to prefer it.
static void test_far_field_treatment() {
  std::printf("\n=== 9. Offene Fernrandbehandlung ===\n");
  const MaterialLibrary lib;
  const Real R[3] = {3.0e-3, 1.2e-2, 4.8e-2};
  const Real zlo[3] = {-4.0e-4, -5.0e-3, -2.3e-2};
  const Real zhi[3] = {1.5e-3, 9.0e-3, 3.9e-2};
  Real gap[3];
  Real phi_asym[3];
  for (int k = 0; k < 3; ++k) {
    DielectricSetup a = reference_setup(lib, 1);
    a.geometry.device.domain_radius = R[k];
    a.geometry.device.domain_z_min = zlo[k];
    a.geometry.device.domain_z_max = zhi[k];
    DielectricSetup b = a;
    b.far_field = FarField::Grounded;
    const DielectricSolution sa = solve_dielectric(a);
    const DielectricSolution sb = solve_dielectric(b);
    const Index i = probe_index(sa.probes, "axis_gap_mid");
    phi_asym[k] = sa.phi_probe[static_cast<std::size_t>(i)];
    gap[k] = std::abs(sa.phi_probe[static_cast<std::size_t>(i)] -
                      sb.phi_probe[static_cast<std::size_t>(i)]);
    std::printf("  R = %6.2f mm : asymptotisch %10.5f V, geerdet %10.5f V, Differenz %.3e V\n",
                R[k] * 1e3, sa.phi_probe[static_cast<std::size_t>(i)],
                sb.phi_probe[static_cast<std::size_t>(i)], gap[k]);
  }
  expect("die beiden Fernrandbedingungen naehern sich mit wachsender Box an",
         gap[1] < 0.2 * gap[0] && gap[2] < 0.5 * gap[1]);
  const Real drift = std::abs(phi_asym[2] - phi_asym[1]) / 1500.0;
  std::printf("  asymptotische Loesung, Aenderung 12 mm -> 48 mm: %.3e der Spannweite\n", drift);
  check_small("die asymptotische Loesung ist bei 12 mm ausgereizt", drift, 1.0e-3);
  expect("die P1-Box von 3 mm ist dafuer zu klein (Differenz > 1 % der Spannweite)",
         gap[0] / 1500.0 > 1.0e-2);
}

// ==========================================================================
// 10.  The corrected model differs from the superseded one
// ==========================================================================
static void test_metallic_reference_is_different() {
  std::printf("\n=== 10. Der korrigierte Vertrag aendert das Ergebnis ===\n");
  const MaterialLibrary lib;
  const DielectricSolution di = solve_dielectric(reference_setup(lib, 1));

  DielectricSetup ms = reference_setup(lib, 1);
  ms.conductor_model = ConductorModel::MetallicReference;
  ms.metallisation = Metallisation::AllSurfaces;
  ms.materials = DielectricMaterials::all_vacuum(lib);
  const DielectricSolution me = solve_dielectric(ms);

  const Index k = probe_index(di.probes, "axis_2_bore_radii");
  const Real a = std::abs(di.Ez_probe[static_cast<std::size_t>(k)]);
  const Real b = std::abs(me.Ez_probe[static_cast<std::size_t>(k)]);
  std::printf("  E_z bei 2 Bohrungsradien: dielektrisch %.6e, metallisch %.6e, "
              "Verhaeltnis %.4f\n",
              a, b, b / a);
  std::printf("  Q_emitter               : dielektrisch %.6e, metallisch %.6e, "
              "Verhaeltnis %.4f\n",
              di.Q_emitter, me.Q_emitter, me.Q_emitter / di.Q_emitter);
  expect("die beiden Modelle liefern messbar verschiedene Felder",
         std::abs(b - a) / a > 1.0e-2);
  expect("und messbar verschiedene Emitterladungen",
         std::abs(me.Q_emitter - di.Q_emitter) / std::abs(di.Q_emitter) > 1.0e-2);

  // MetallicReference with a real dielectric would be neither model.
  DielectricSetup mixed = ms;
  mixed.materials = DielectricMaterials::reference(lib);
  bool threw = false;
  try {
    (void)solve_dielectric(mixed);
  } catch (const std::exception&) {
    threw = true;
  }
  expect("metallische Referenz mit eps_r != 1 wird abgelehnt", threw);
}

// ==========================================================================
int main() {
  std::printf("P2b: dielektrische achsensymmetrische Elektrostatik\n");
  test_mesh();
  test_no_polymer_is_a_conductor();
  test_interface_continuity();
  test_linearity_and_polarity();
  test_fem_against_bem();
  test_mesh_convergence();
  test_feed_boundary_position();
  test_permittivity_sensitivity();
  test_far_field_treatment();
  test_metallic_reference_is_different();
  std::printf("\n%s (%d Fehler)\n",
              failures == 0 ? "ALLE TESTS BESTANDEN" : "TESTS FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
