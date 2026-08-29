// P2a: adapter from the P1 boundary mesh to the vacuum BEM, and the numerical
// properties the resulting solve must have.
//
// The point is NOT to re-run the mesher's own validate() and call that a test.
// Everything here checks a property from the outside: against a closed-form
// two-conductor solution, against a deliberately wrong input, against the same
// problem at a different voltage, polarity or mesh density.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/boundary_mesh.hpp"
#include "es/constants.hpp"
#include "es/device_geometry.hpp"
#include "es/vacuum_bem.hpp"

using namespace es;
using constants::eps0;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-60s got=%-14.7g want=%-14.7g %s\n", what, got, want, ok ? "OK" : "FAIL");
}

static void below(const char* what, Real got, Real limit) {
  const bool ok = got <= limit;
  if (!ok) ++failures;
  std::printf("  %-60s got=%-14.7g <= %-12.4g %s\n", what, got, limit, ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-60s %s\n", what, ok ? "OK" : "FAIL");
}

static DeviceParameters example() {
  DeviceParameters p;  // the documented P1 example set
  return p;
}

// ==========================================================================
// 1.  The mandatory extractor outer radius
// ==========================================================================
static bool rejects(DeviceParameters p) {
  try {
    (void)DeviceGeometry::build(p);
    return false;
  } catch (const std::exception&) {
    return true;
  }
}

static void test_extractor_outer_radius_is_mandatory() {
  std::printf("\n=== Aussenradius der Extraktionselektrode ist Pflicht ===\n");
  {
    DeviceParameters p = example();
    p.extractor_outer_radius = 0.0;
    expect("0 (frueher 'bis zum Domaenenrand') -> abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example();
    p.extractor_outer_radius = p.domain_radius;
    expect("gleich domain_radius -> abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example();
    p.extractor_outer_radius = 1.01 * p.domain_radius;
    expect("groesser als domain_radius -> abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example();
    p.extractor_outer_radius = 0.4 * p.extractor_aperture_diameter;
    expect("kleiner als der Aperturradius -> abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example();
    const DeviceGeometry g = DeviceGeometry::build(p);
    expect("gueltiger Wert -> akzeptiert", g.extractor_outer_radius() == p.extractor_outer_radius);
    // With a finite outer radius the electrode is a CLOSED body of revolution:
    // aperture wall, both faces and the outer rim, and nothing of it lies on
    // the domain edge.
    bool has_rim = false, extractor_on_domain_edge = false;
    for (const BoundaryCurve& c : g.boundaries()) {
      if (c.name == "extractor_surface.rim") has_rim = true;
      if (c.id == BoundaryId::OpenBoundary &&
          (c.side_a == Region::ExtractorSolid || c.side_b == Region::ExtractorSolid))
        extractor_on_domain_edge = true;
    }
    expect("die Elektrode hat eine eigene Mantelflaeche (rim)", has_rim);
    expect("kein Stueck der Elektrode liegt auf dem offenen Domaenenrand",
           !extractor_on_domain_edge);
  }
}

// ==========================================================================
// 2.  The adapter selects by tag, and only conductor/vacuum interfaces
// ==========================================================================
static void test_selection() {
  std::printf("\n=== Randauswahl fuer die Vakuum-BEM ===\n");
  const DeviceGeometry g = DeviceGeometry::build(example());
  const BoundaryMesh bm = BoundaryMesh::generate(g);
  VacuumSelectionReport rep;
  const Mesh m = vacuum_bem_mesh(bm, g, &rep);

  std::set<std::string> accepted, rejected;
  for (const VacuumSelectionReport::CurveDecision& c : rep.curves)
    (c.accepted ? accepted : rejected).insert(c.curve);

  for (const char* n : {"free_surface_reference", "emitter_tip_land", "emitter_outer_surface",
                        "extractor_surface.aperture", "extractor_surface.front",
                        "extractor_surface.back", "extractor_surface.rim"}) {
    char buf[96];
    std::snprintf(buf, sizeof buf, "uebernommen: %s", n);
    expect(buf, accepted.count(n) == 1);
  }
  for (const char* n : {"symmetry_axis.liquid", "symmetry_axis.vacuum", "bore_wall",
                        "liquid_inlet", "open_boundary.z_min.emitter",
                        "open_boundary.z_min.vacuum", "open_boundary.r_max",
                        "open_boundary.z_max"}) {
    char buf[96];
    std::snprintf(buf, sizeof buf, "NICHT uebernommen: %s", n);
    expect(buf, rejected.count(n) == 1);
  }

  // No element may reach the solver without a conductor tag: Tag::Other and
  // Tag::Collector both map onto Electrode::Collector in bem.cpp, so an
  // untagged panel would silently become a third, grounded electrode.
  bool only_known = true, no_axis = true, positive_area = true;
  for (const Element& e : m.elems) {
    only_known = only_known && (e.tag == Tag::Emitter || e.tag == Tag::FreeSurface ||
                                e.tag == Tag::Extractor);
    no_axis = no_axis && e.mid.r > 0.0;
    positive_area = positive_area && e.area > 0.0 && e.len > 0.0;
  }
  expect("nur Emitter-, Freiflaechen- und Extraktorkennungen", only_known);
  expect("kein Achsenelement, alle Mittelpunkte bei r > 0", no_axis);
  expect("alle Rotationsflaechen und Laengen positiv", positive_area);

  // Areas must be exactly what the mesher computed, summed per boundary id.
  Real a_free = 0.0, a_emit = 0.0, a_ext = 0.0;
  for (BoundaryId id : {BoundaryId::EmitterTipLand, BoundaryId::EmitterOuterSurface})
    a_emit += bm.stats_of(id).total_revolved_area;
  a_free = bm.stats_of(BoundaryId::FreeSurfaceReference).total_revolved_area;
  a_ext = bm.stats_of(BoundaryId::ExtractorSurface).total_revolved_area;
  check("Rotationsflaeche Emitter aus dem Netz uebernommen", rep.revolved_area_emitter, a_emit,
        1e-15);
  check("Rotationsflaeche ebene Fluessigkeitsoberflaeche", rep.revolved_area_free_surface, a_free,
        1e-15);
  check("Rotationsflaeche Extraktor", rep.revolved_area_extractor, a_ext, 1e-15);

  // The flat liquid surface must be an annulus of the bore radius, exactly.
  check("Flaeche der ebenen Fluessigkeitsoberflaeche = pi r_bore^2",
        rep.revolved_area_free_surface, pi * g.contact_radius() * g.contact_radius(), 1e-13);

  // Outward normals must point into the vacuum.  Checked against the region
  // the geometry says lies on either side, not against a hard-coded direction.
  bool normals_ok = true;
  for (std::size_t k = 0; k < rep.panels.size(); ++k) {
    const VacuumPanel& p = rep.panels[k];
    const BoundaryElement& be = bm.elements()[static_cast<std::size_t>(p.mesh_element)];
    const Element& el = m.elems[static_cast<std::size_t>(p.bem_element)];
    const Vec2 into_vacuum = (be.side_b == Region::Vacuum) ? be.normal : -1.0 * be.normal;
    normals_ok = normals_ok && dot(el.normal, into_vacuum) > 0.999999;
    normals_ok = normals_ok && std::abs(norm(el.normal) - 1.0) < 1e-12;
    // The BEM convention: outward normal = perp(tangent).
    normals_ok = normals_ok && norm(el.normal - perp(el.tangent)) < 1e-12;
  }
  expect("Normale zeigt ins Vakuum und ist perp(tangent)", normals_ok);

  // Changing the electrode radius must not change WHICH boundaries are taken.
  DeviceParameters q = example();
  q.extractor_outer_radius = 1.0e-3;
  const DeviceGeometry g2 = DeviceGeometry::build(q);
  VacuumSelectionReport rep2;
  (void)vacuum_bem_mesh(BoundaryMesh::generate(g2), g2, &rep2);
  std::set<std::string> accepted2;
  for (const VacuumSelectionReport::CurveDecision& c : rep2.curves)
    if (c.accepted) accepted2.insert(c.curve);
  expect("Auswahl haengt an den Kennungen, nicht an den Abmessungen", accepted == accepted2);

  // And it must not depend on the mesh density either.
  VacuumSelectionReport rep3;
  (void)vacuum_bem_mesh(BoundaryMesh::generate(g, 0.3), g, &rep3);
  std::set<std::string> accepted3;
  for (const VacuumSelectionReport::CurveDecision& c : rep3.curves)
    if (c.accepted) accepted3.insert(c.curve);
  expect("Auswahl haengt nicht von der Netzstufe ab", accepted == accepted3);
  expect("feineres Netz liefert mehr Panels", rep3.n_selected > rep.n_selected);
}

// ==========================================================================
// 3.  The far-field boundary condition really is V -> 0
// ==========================================================================
static void test_boundary_condition_at_infinity() {
  std::printf("\n=== Randbedingung im Unendlichen ===\n");
  const Real R = 1.0e-3, V = 1.0;
  Mesh sphere = make_sphere(R, V, 240);
  BemSolver bem(sphere);
  bem.solve({{V, 0.0, 0.0}});
  check("isolierte Kugel: C = 4 pi eps0 R", bem.total_charge() / V, 4.0 * pi * eps0 * R, 2e-3);
  // A free-space kernel means the potential falls off as Q/(4 pi eps0 d) with
  // no additive constant.  A solver with a grounded box at finite distance
  // would show a deficit here.
  const Real Q = bem.total_charge();
  for (Real f : {50.0, 500.0, 5000.0}) {
    const Real d = f * R;
    char buf[96];
    std::snprintf(buf, sizeof buf, "V(%.0f R) = Q/(4 pi eps0 d)", f);
    check(buf, bem.potential_at({0.0, d}), Q / (4.0 * pi * eps0 * d), 3e-3);
  }
}

// ==========================================================================
// 4.  Maxwell capacitance matrix against a closed-form two-conductor case
// ==========================================================================
static void test_maxwell_matrix_two_spheres() {
  std::printf("\n=== Maxwell-Kapazitaetsmatrix, zwei weit getrennte Kugeln ===\n");
  const Real R1 = 1.0e-3, R2 = 2.0e-3, d = 8.0e-2;  // d >> R1, R2
  Mesh a = make_sphere(R1, 0.0, 200);
  Mesh b = make_sphere(R2, 0.0, 200);
  for (Element& e : b.elems) { e.a.z += d; e.b.z += d; e.mid.z += d; e.tag = Tag::Extractor; }
  for (Element& e : b.elems) e.body = 1;
  Mesh both = a;
  both.elems.insert(both.elems.end(), b.elems.begin(), b.elems.end());

  BemSolver bem(both);
  bem.solve({{1.0, 0.0, 0.0}});
  const CapacitanceMatrix c = maxwell_capacitance(bem);

  // Leading order in R/d:  c_ii = 4 pi eps0 R_i, c_ij = -4 pi eps0 R_i R_j / d.
  check("c_EE = 4 pi eps0 R1", c.c_EE, 4.0 * pi * eps0 * R1, 3e-3);
  check("c_XX = 4 pi eps0 R2", c.c_XX, 4.0 * pi * eps0 * R2, 3e-3);
  check("c_EX = -4 pi eps0 R1 R2 / d", c.c_EX, -4.0 * pi * eps0 * R1 * R2 / d, 5e-3);
  below("Reziprozitaet c_EX = c_XE", c.reciprocity_error(), 1e-9);
  expect("Selbstkoeffizienten positiv, Influenzkoeffizient negativ",
         c.c_EE > 0 && c.c_XX > 0 && c.c_EX < 0);
  expect("C_m = -c_EX > 0", c.mutual() > 0.0);
}

// ==========================================================================
// 5.  The device solve: linearity, polarity, residual, screening
// ==========================================================================
struct Device {
  DeviceGeometry g;
  BoundaryMesh bm;
  BemSolver bem;
  VacuumSelectionReport sel;
  std::vector<EdgeZone> zones;
};

static void build(Device& d, Real scale, const std::array<Real, 3>& V) {
  d.g = DeviceGeometry::build(example());
  d.bm = BoundaryMesh::generate(d.g, scale);
  Mesh panels = vacuum_bem_mesh(d.bm, d.g, &d.sel);
  d.zones = edge_zones(d.g, d.bm, panels);
  d.bem.set_mesh(std::move(panels));
  d.bem.solve(V);
}

static void test_linearity_and_polarity() {
  std::printf("\n=== Linearitaet und Polaritaetsumkehr ===\n");
  const Real VE = 1500.0, VX = 0.0;
  Device d;
  build(d, 0.5, {{VE, VX, 0.0}});
  std::printf("  (%d Panels)\n", static_cast<int>(d.bem.size()));

  const std::vector<Real> s1 = d.bem.sigma_for({{VE, VX, 0.0}});
  Real smax = 0.0;
  for (Real v : s1) smax = std::max(smax, std::abs(v));

  for (Real f : {0.5, 2.0, -1.0, -3.7}) {
    const std::vector<Real> s = d.bem.sigma_for({{f * VE, f * VX, 0.0}});
    Real dev = 0.0;
    for (std::size_t i = 0; i < s.size(); ++i)
      dev = std::max(dev, std::abs(s[i] - f * s1[i]) / (std::abs(f) * smax));
    char buf[96];
    std::snprintf(buf, sizeof buf, "sigma(%.2f V) = %.2f sigma(V)", f, f);
    below(buf, dev, 1e-12);
  }

  // A second, independent voltage direction: superposition must hold when the
  // extractor is not grounded either.
  {
    const std::vector<Real> sa = d.bem.sigma_for({{1000.0, 0.0, 0.0}});
    const std::vector<Real> sb = d.bem.sigma_for({{0.0, -400.0, 0.0}});
    const std::vector<Real> sc = d.bem.sigma_for({{1000.0, -400.0, 0.0}});
    Real dev = 0.0, sc_max = 0.0;
    for (Real v : sc) sc_max = std::max(sc_max, std::abs(v));
    for (std::size_t i = 0; i < sc.size(); ++i)
      dev = std::max(dev, std::abs(sc[i] - sa[i] - sb[i]) / sc_max);
    below("Superposition beider Elektroden", dev, 1e-12);
  }

  // Full polarity reversal.
  BemSolver plus = d.bem, minus = d.bem;
  plus.solve({{VE, VX, 0.0}});
  minus.solve({{-VE, -VX, 0.0}});
  Real dsig = 0.0, dabs = 0.0;
  for (Index i = 0; i < plus.size(); ++i) {
    const Real a = plus.sigma()[static_cast<std::size_t>(i)];
    const Real b = minus.sigma()[static_cast<std::size_t>(i)];
    dsig = std::max(dsig, std::abs(a + b) / smax);
    dabs = std::max(dabs, std::abs(std::abs(a) - std::abs(b)) / smax);
  }
  below("sigma wechselt exakt das Vorzeichen", dsig, 1e-12);
  below("|sigma| bleibt gleich", dabs, 1e-12);

  const Vec2 x{0.0, 0.1 * d.g.contact_radius()};
  const Vec2 ep = plus.field_at(x), em = minus.field_at(x);
  below("E wechselt das Vorzeichen", norm(ep + em) / norm(ep), 1e-12);
  below("|E| bleibt gleich", std::abs(norm(ep) - norm(em)) / norm(ep), 1e-12);
  below("V wechselt das Vorzeichen",
        std::abs(plus.potential_at(x) + minus.potential_at(x)) / std::abs(plus.potential_at(x)),
        1e-12);
  check("Ladungen kehren sich um", minus.charge_on(Tag::Emitter), -plus.charge_on(Tag::Emitter),
        1e-12);
}

static void test_residual_and_screening() {
  std::printf("\n=== Dirichlet-Residuum und Abschirmung des offenen Emitterbogens ===\n");
  const Real VE = 1500.0;
  Device d;
  build(d, 0.5, {{VE, 0.0, 0.0}});

  const PotentialResidual r = potential_residual(d.bem, d.zones);
  std::printf("  ausserhalb der Kantenzonen: Emitter %.3e V, Extraktor %.3e V\n",
              r.max_emitter_clear, r.max_extractor_clear);
  std::printf("  einschliesslich:            %.3e V / %.3e V, groesstes bei r=%.3e z=%.3e\n",
              r.max_emitter, r.max_extractor, r.worst_position.r, r.worst_position.z);
  below("relatives Residuum, Kantenzonen ausgeschlossen", r.relative(), 1e-3);
  below("relatives Residuum einschliesslich Kantenzonen", r.relative_including_edges(), 2e-2);

  // The edge-dominated residual must not be a plateau: it has to fall with
  // refinement, or the singularity is not being resolved at all.
  Device fine;
  build(fine, 0.125, {{VE, 0.0, 0.0}});
  const PotentialResidual rf = potential_residual(fine.bem, fine.zones);
  std::printf("  Verfeinerung 0.5 -> 0.125: %.3e -> %.3e V (mit Kantenzonen)\n",
              std::max(r.max_emitter, r.max_extractor),
              std::max(rf.max_emitter, rf.max_extractor));
  expect("Residuum an den Singularitaeten faellt mit der Verfeinerung",
         rf.relative_including_edges() < r.relative_including_edges());
  expect("Residuum ausserhalb der Kantenzonen faellt ebenfalls",
         rf.relative() < r.relative());

  // The emitter arc is open where the domain floor cuts the shank.  If the
  // cavity behind it were not screened, sigma/eps0 would be a two-sided sum
  // instead of the vacuum-side field, and the whole tip result with it.
  Real worst = 0.0;
  const Real H = d.g.parameters().emitter_height;
  for (int k = 1; k <= 10; ++k) {
    const Real z = -H * (0.05 + 0.9 * k / 10.0);
    for (Real rr : {0.0, 0.4 * d.g.contact_radius()})
      worst = std::max(worst, std::abs(d.bem.potential_at({rr, z}) - VE));
  }
  std::printf("  max |V - V_emitter| im Schafthohlraum: %.3e V von %.4g V\n", worst, VE);
  below("Hohlraum abgeschirmt (relativ)", worst / VE, 1e-3);
}

// ==========================================================================
// 6.  Mesh convergence
// ==========================================================================
static void test_mesh_convergence() {
  std::printf("\n=== Netzkonvergenz ueber drei automatisch erzeugte Stufen ===\n");
  const Real VE = 1500.0;
  const std::array<Real, 3> V{{VE, 0.0, 0.0}};
  const DeviceGeometry g = DeviceGeometry::build(example());
  const Vec2 p_ref{0.0, 0.1 * g.contact_radius()};

  std::vector<Real> cEE, Cm, Eref;
  std::vector<int> n;
  for (Real s : {1.0, 0.5, 0.25}) {
    const BoundaryMesh bm = BoundaryMesh::generate(g, s);
    BemSolver bem(vacuum_bem_mesh(bm, g));
    bem.solve(V);
    const CapacitanceMatrix c = maxwell_capacitance(bem);
    cEE.push_back(c.c_EE);
    Cm.push_back(c.mutual());
    Eref.push_back(bem.field_at(p_ref).z);
    n.push_back(static_cast<int>(bem.size()));
    std::printf("  scale %.3f : %4d Panels, c_EE = %.6e F, C_m = %.6e F, E_z(ref) = %.6e V/m\n", s,
                n.back(), cEE.back(), Cm.back(), Eref.back());
  }
  expect("jede Stufe hat mehr Panels als die vorige", n[1] > n[0] && n[2] > n[1]);

  auto shrinks = [](const std::vector<Real>& v, const char* what, Real tol) {
    const Real d1 = std::abs(v[1] - v[0]) / std::abs(v[2]);
    const Real d2 = std::abs(v[2] - v[1]) / std::abs(v[2]);
    std::printf("  %-34s |d1| = %.3e -> |d2| = %.3e\n", what, d1, d2);
    if (!(d2 < d1)) { ++failures; std::printf("    FAIL: Aenderung nimmt nicht ab\n"); }
    if (!(d2 <= tol)) { ++failures; std::printf("    FAIL: letzte Aenderung > %.3g\n", tol); }
    return d2;
  };
  shrinks(cEE, "c_EE", 5e-3);
  shrinks(Cm, "C_m", 5e-3);
  shrinks(Eref, "E_z am Referenzpunkt", 5e-3);
}

// ==========================================================================
// 7.  Edge zones are marked, and the peak field outside them is the one used
// ==========================================================================
static void test_edge_zones() {
  std::printf("\n=== Kantenzonen ===\n");
  const DeviceGeometry g = DeviceGeometry::build(example());
  const BoundaryMesh bm = BoundaryMesh::generate(g, 0.5);
  const Mesh panels = vacuum_bem_mesh(bm, g);
  const std::vector<EdgeZone> z = edge_zones(g, bm, panels);
  int sharp = 0, trunc = 0;
  for (const EdgeZone& e : z) (e.kind == EdgeKind::SharpFeature ? sharp : trunc)++;
  expect("eine Zone je benanntem Merkmal", sharp == static_cast<int>(g.features().size()));
  // The emitter arc is cut exactly once, at the domain floor.  The electrode is
  // closed, and the axis end of the emitter arc is a closure, not an edge.
  expect("genau ein offenes Bogenende markiert", trunc == 1);
  bool at_floor = true;
  for (const EdgeZone& e : z)
    if (e.kind == EdgeKind::TruncationEnd)
      at_floor = at_floor && e.position.z == g.parameters().domain_z_min &&
                 e.position.r == 0.5 * g.parameters().phi_3;
  expect("das offene Ende liegt am Schnitt durch die Domaenensohle", at_floor);
  bool geometric = true;
  for (const EdgeZone& e : z)
    geometric = geometric && e.radius > 0.0 && e.radius < 0.5 * e.local_feature_size;
  expect("Radius rein geometrisch und positiv", geometric);

  // The marking must not move when the mesh is refined -- otherwise the
  // "converged" peak field would be reported over a shrinking region.
  const BoundaryMesh bmf = BoundaryMesh::generate(g, 0.125);
  const std::vector<EdgeZone> z2 = edge_zones(g, bmf, vacuum_bem_mesh(bmf, g));
  bool same = z.size() == z2.size();
  for (std::size_t i = 0; same && i < z.size(); ++i)
    same = same && z[i].radius == z2[i].radius && norm(z[i].position - z2[i].position) == 0.0;
  expect("Kantenzonen sind netzunabhaengig", same);

  Device d;
  build(d, 0.5, {{1500.0, 0.0, 0.0}});
  Index which = -1;
  const Real outside = peak_field_outside_edges(d.bem, Electrode::Emitter, d.zones, &which);
  const Real everywhere = d.bem.peak_emitter_field();
  expect("das Maximum liegt an einer scharfen Kante", everywhere > outside);
  expect("der berichtete Wert stammt von ausserhalb der Zonen",
         which >= 0 &&
             !in_edge_zone(d.zones,
                           d.bem.mesh().elems[static_cast<std::size_t>(which)].mid));

  // And the edge value must indeed fail to converge, while the reported one does.
  Device fine;
  build(fine, 0.125, {{1500.0, 0.0, 0.0}});
  const Real edge_c = d.bem.peak_emitter_field();
  const Real edge_f = fine.bem.peak_emitter_field();
  const Real out_c = outside;
  const Real out_f = peak_field_outside_edges(fine.bem, Electrode::Emitter, fine.zones);
  std::printf("  an der Kante  : %.4e -> %.4e V/m  (Aenderung %.1f %%)\n", edge_c, edge_f,
              100.0 * std::abs(edge_f - edge_c) / edge_c);
  std::printf("  ausserhalb    : %.4e -> %.4e V/m  (Aenderung %.1f %%)\n", out_c, out_f,
              100.0 * std::abs(out_f - out_c) / out_f);
  expect("Kantenwert waechst mit der Verfeinerung, der berichtete nicht",
         std::abs(edge_f - edge_c) / edge_c > 4.0 * std::abs(out_f - out_c) / out_f);
}

// ==========================================================================
int main() {
  std::printf("=== P2a: Vakuum-BEM auf der achsensymmetrischen P1-Geometrie ===\n");
  test_extractor_outer_radius_is_mandatory();
  test_selection();
  test_boundary_condition_at_infinity();
  test_maxwell_matrix_two_spheres();
  test_linearity_and_polarity();
  test_residual_and_screening();
  test_mesh_convergence();
  test_edge_zones();
  std::printf("\n%s (%d Fehler)\n", failures ? "FEHLGESCHLAGEN" : "bestanden", failures);
  return failures ? 1 : 0;
}
