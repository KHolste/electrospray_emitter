// Automatic axisymmetric boundary mesher: size function, topology, the
// axisymmetric contract, and determinism.
//
// The mesher owns its own validate(), and the point of these tests is NOT to
// re-run it and call that a test.  validate() is exercised once, and everything
// else here checks a property from the outside -- against a closed-form value,
// against a deliberately broken input, or against a second run.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/boundary_mesh.hpp"
#include "es/constants.hpp"
#include "es/device_geometry.hpp"

using namespace es;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-58s got=%-14.7g want=%-14.7g %s\n", what, got, want, ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-58s %s\n", what, ok ? "OK" : "FAIL");
}

static DeviceParameters example() {
  DeviceParameters p;  // the defaults are the documented P1 example set
  return p;
}

// ==========================================================================
// 1.  Size function
// ==========================================================================
static void test_size_field() {
  std::printf("\n=== Groessenfunktion ===\n");
  const DeviceParameters p = example();
  const DeviceGeometry g = DeviceGeometry::build(p);
  const SizeField h = SizeField::from_geometry(g);

  // Every source size is positive and every evaluation is bounded by h_max.
  bool positive = true;
  for (const auto& s : h.sources()) positive = positive && s.h > 0.0 && s.local_feature_size > 0.0;
  expect("alle Quellgroessen positiv", positive);

  const Real diag = std::sqrt(p.domain_radius * p.domain_radius +
                              (p.domain_z_max - p.domain_z_min) * (p.domain_z_max - p.domain_z_min));
  check("h_max = Diagonale / 40", h.h_max(), diag / mesher::kDomainDivisions, 1e-14);

  // The pinned exit edge is the finest place in the device.  Its local feature
  // size is min(bore radius, land width) = 5 um, refined by kFeatureDivisions.
  const Vec2 tip = g.feature(FeatureId::PinnedContactEdge);
  const Real lfs_tip = std::min(0.5 * p.phi_2, 0.5 * (p.phi_1 - p.phi_2));
  check("h an der gepinnten Austrittskante", h(tip), lfs_tip / mesher::kFeatureDivisions, 1e-12);

  // The aperture edges are refined by their own local scale, the electrode
  // thickness, not by the aperture diameter.
  const Vec2 ap = g.feature(FeatureId::ExtractorApertureEdgeFront);
  check("h an der vorderen Aperturkante", h(ap),
        p.extractor_thickness / mesher::kFeatureDivisions, 1e-12);

  // Strong refinement means: markedly finer than the far field.
  expect("Austrittskante mindestens 100x feiner als das Fernfeld",
         h(tip) * 100.0 < h.h_max());
  expect("Aperturkante deutlich feiner als das Fernfeld", h(ap) * 10.0 < h.h_max());

  // Far away the ceiling is actually reached, so the far field is coarse.
  check("h am fernen Domaenenrand ist h_max",
        h({0.5 * p.domain_radius, 0.9 * p.domain_z_max}), h.h_max(), 1e-14);

  // The field is G-Lipschitz by construction; verify it on random-free,
  // deterministic sample pairs spanning the domain.
  Real worst = 0.0;
  for (int i = 0; i <= 40; ++i)
    for (int j = 0; j <= 40; ++j) {
      const Vec2 x{p.domain_radius * i / 40.0,
                   p.domain_z_min + (p.domain_z_max - p.domain_z_min) * j / 40.0};
      const Vec2 y{p.domain_radius * j / 40.0,
                   p.domain_z_min + (p.domain_z_max - p.domain_z_min) * i / 40.0};
      const Real d = norm(x - y);
      if (d > 0.0) worst = std::max(worst, std::abs(h(x) - h(y)) / d);
    }
  expect("Groessenfeld ist G-Lipschitz (Wachstum <= G)", worst <= mesher::kGradation * 1.000001);
  std::printf("  %-58s %.6g <= %.6g\n", "  groesste beobachtete Steigung", worst,
              mesher::kGradation);

  // No user knob exists: the field is a pure function of the geometry.  Halving
  // the bore must refine the tip, without any parameter being passed.
  DeviceParameters q = example();
  q.phi_2 = 0.5 * p.phi_2;
  const DeviceGeometry g2 = DeviceGeometry::build(q);
  const SizeField h2 = SizeField::from_geometry(g2);
  expect("kleinere Bohrung verfeinert die Spitze von selbst",
         h2(g2.feature(FeatureId::PinnedContactEdge)) < h(tip));
}

// ==========================================================================
// 2.  The axisymmetric contract
// ==========================================================================
static void test_axisymmetric_contract() {
  std::printf("\n=== Achsensymmetrischer Vertrag ===\n");
  const DeviceParameters p = example();
  const DeviceGeometry g = DeviceGeometry::build(p);
  const BoundaryMesh m = BoundaryMesh::generate(g);

  expect("Netz ist nicht leer", m.elements().size() > 50);

  // Axis elements are symmetry elements, not rings.
  int axis = 0, ring = 0;
  bool axis_clean = true, ring_clean = true;
  for (const BoundaryElement& e : m.elements()) {
    if (e.is_axis()) {
      ++axis;
      axis_clean = axis_clean && e.a.r == 0.0 && e.b.r == 0.0 && e.revolved_area == 0.0 &&
                   e.id == BoundaryId::SymmetryAxis && e.side_a == e.side_b;
    } else {
      ++ring;
      ring_clean = ring_clean && e.mid_radius() > 0.0 && e.revolved_area > 0.0 &&
                   e.meridian_length > 0.0;
    }
  }
  expect("Achsenelemente: r = 0, Rotationsflaeche exakt 0, kein Interface", axis_clean);
  expect("Ringelemente: mittlerer Radius > 0 und Flaeche > 0", ring_clean);
  std::printf("  %-58s %d Symmetrie- / %d Ringelemente\n", "  Aufteilung", axis, ring);

  // Every element on the axis really belongs to the axis boundary, and vice
  // versa -- a ring element with two nodes at r = 0 would be a silent zero.
  bool consistent = true;
  for (const BoundaryElement& e : m.elements())
    consistent = consistent &&
                 ((e.a.r == 0.0 && e.b.r == 0.0) == (e.kind == ElementKind::AxisSymmetry));
  expect("Kennzeichnung Symmetrie/Ring folgt exakt aus r = 0", consistent);

  // Surfaces of revolution are preserved exactly, boundary by boundary.  These
  // are closed-form values, not the mesher repeating itself.
  const Real r1 = 0.5 * p.phi_1, r2 = 0.5 * p.phi_2, r3 = 0.5 * p.phi_3;
  auto area_of = [&m](BoundaryId id) { return m.stats_of(id).total_revolved_area; };

  check("Bohrungswand: Sum(2 pi r ds) = 2 pi r2 L", area_of(BoundaryId::BoreWall),
        2 * pi * r2 * (0.0 - p.domain_z_min), 1e-13);
  check("Stirnflaeche: Sum = pi (r1^2 - r2^2)", area_of(BoundaryId::EmitterTipLand),
        pi * (r1 * r1 - r2 * r2), 1e-13);
  check("anfaengliche ebene Fluessigkeitsoberflaeche: Sum = pi r2^2",
        area_of(BoundaryId::FreeSurfaceReference), pi * r2 * r2, 1e-13);
  check("Zulaufschnitt: Sum = pi r2^2", area_of(BoundaryId::LiquidInlet), pi * r2 * r2, 1e-13);
  check("Symmetrieachse: Rotationsflaeche exakt 0",
        area_of(BoundaryId::SymmetryAxis) + 1.0, 1.0, 0.0);

  // Emitter flank: truncated cone plus cylinder, both closed form.
  const Real slant = std::sqrt((r3 - r1) * (r3 - r1) + p.emitter_height * p.emitter_height);
  check("Emitteraussenflaeche: Kegelstumpf + Zylinder",
        area_of(BoundaryId::EmitterOuterSurface),
        pi * (r1 + r3) * slant + 2 * pi * r3 * (-p.domain_z_min - p.emitter_height), 1e-13);

  // Meridian lengths add up per curve.
  Real worst_len = 0.0;
  for (std::size_t c = 0; c < g.boundaries().size(); ++c) {
    const Real got = m.stats_of_curve(static_cast<int>(c)).total_meridian_length;
    const Real want = g.boundaries()[c].meridian_length();
    worst_len = std::max(worst_len, std::abs(got - want) / want);
  }
  check("Meridianlaengen je Kurve erhalten", worst_len + 1.0, 1.0, 1e-12);

  // The normal points from side_a into side_b and equals perp(tangent).
  bool normals_ok = true;
  for (const BoundaryElement& e : m.elements()) {
    if (e.is_axis()) continue;
    const Vec2 want = perp(e.tangent);
    normals_ok = normals_ok && std::abs(e.normal.r - want.r) < 1e-14 &&
                 std::abs(e.normal.z - want.z) < 1e-14 &&
                 std::abs(norm(e.normal) - 1.0) < 1e-14 && std::abs(dot(e.normal, e.tangent)) < 1e-14;
  }
  expect("Normale = perp(Tangente), normiert und orthogonal", normals_ok);

  // Spot-check the direction against physics-free geometry: the initial flat
  // liquid surface separates liquid (below) from vacuum (above), so its normal
  // must point towards +z.
  bool fs_ok = true;
  int fs_n = 0;
  for (const BoundaryElement* e : m.elements_with(BoundaryId::FreeSurfaceReference)) {
    ++fs_n;
    fs_ok = fs_ok && e->side_a == Region::Liquid && e->side_b == Region::Vacuum &&
            e->normal.z > 0.999 && e->a.z == 0.0 && e->b.z == 0.0;
  }
  expect("ebene Fluessigkeitsoberflaeche: Normale zeigt in +z ins Vakuum", fs_ok && fs_n > 0);

  // And the bore wall separates liquid from emitter solid, normal towards +r.
  bool bore_ok = true;
  for (const BoundaryElement* e : m.elements_with(BoundaryId::BoreWall))
    bore_ok = bore_ok && e->side_a == Region::Liquid && e->side_b == Region::EmitterSolid &&
              e->normal.r > 0.999;
  expect("Bohrungswand: Normale zeigt von der Fluessigkeit in den Festkoerper", bore_ok);
}

// ==========================================================================
// 3.  Topology: corners, connectivity, no merging across identifiers
// ==========================================================================
static void test_topology() {
  std::printf("\n=== Topologie ===\n");
  const DeviceParameters p = example();
  const DeviceGeometry g = DeviceGeometry::build(p);
  const BoundaryMesh m = BoundaryMesh::generate(g);

  // Every geometric corner and every named feature is a node, bit for bit.
  auto exact_node = [&m](Vec2 x) {
    for (const MeshNode& n : m.nodes())
      if (n.p.r == x.r && n.p.z == x.z) return true;
    return false;
  };
  bool corners = true;
  int n_corner = 0;
  for (const BoundaryCurve& c : g.boundaries())
    for (const Vec2& v : c.points) {
      ++n_corner;
      corners = corners && exact_node(v);
    }
  expect("jede Geometrieecke ist bitgenau ein Netzknoten", corners);
  std::printf("  %-58s %d\n", "  gepruefte Eckpunkte", n_corner);

  bool feats = true;
  for (const NamedFeature& f : g.features()) feats = feats && exact_node(f.position);
  expect("jedes benannte Merkmal ist bitgenau ein Netzknoten", feats);

  // A node marked as carrying a feature must sit exactly on it.
  bool feat_marks = true;
  int marked = 0;
  for (const MeshNode& n : m.nodes())
    if (n.feature >= 0) {
      ++marked;
      const Vec2 f = g.features()[n.feature].position;
      feat_marks = feat_marks && n.p.r == f.r && n.p.z == f.z;
    }
  expect("Merkmalsmarkierungen sitzen auf den Merkmalen", feat_marks);
  check("Anzahl markierter Merkmalsknoten", marked, static_cast<Real>(g.features().size()), 0.0);

  // No element may be shared between two boundary identifiers or two curves,
  // and no element may straddle a corner: both endpoints of every element lie
  // on one straight geometric segment of its own curve.
  bool one_segment = true;
  for (const BoundaryElement& e : m.elements()) {
    const auto& pts = g.boundaries()[e.curve].points;
    const Vec2 s0 = pts[e.segment], s1 = pts[e.segment + 1];
    const Real L = norm(s1 - s0);
    const Vec2 u = (s1 - s0) / L;
    for (Vec2 q : {e.a, e.b}) {
      const Real t = dot(q - s0, u);
      const Real off = norm(q - (s0 + t * u));
      one_segment = one_segment && off <= 1e-12 * L && t >= -1e-12 * L && t <= L * (1 + 1e-12);
    }
    one_segment = one_segment && e.id == g.boundaries()[e.curve].id;
  }
  expect("kein Element ueberspannt eine Ecke oder Materialgrenze", one_segment);

  // Every node is used by at least two elements: no dangling ends, no gaps.
  std::vector<int> degree(m.nodes().size(), 0);
  for (const BoundaryElement& e : m.elements()) {
    ++degree[e.node_a];
    ++degree[e.node_b];
  }
  int min_deg = 1 << 30, junctions = 0;
  for (int d : degree) {
    min_deg = std::min(min_deg, d);
    if (d > 2) ++junctions;
  }
  check("kleinster Knotengrad", min_deg, 2.0, 0.0);
  expect("Verzweigungsknoten existieren (Dreifachpunkte an Kanten)", junctions > 0);

  // Elements of one curve form an unbroken chain.
  bool chained = true;
  for (std::size_t c = 0; c < g.boundaries().size(); ++c) {
    const auto es = m.elements_of_curve(static_cast<int>(c));
    chained = chained && !es.empty();
    for (std::size_t i = 1; i < es.size(); ++i) chained = chained && es[i]->node_a == es[i - 1]->node_b;
  }
  expect("Elemente einer Randkurve bilden eine ununterbrochene Kette", chained);

  // Node coordinates are unique -- two nodes closer than the snap tolerance
  // would mean the same point was inserted twice.
  bool unique = true;
  for (std::size_t i = 0; i < m.nodes().size(); ++i)
    for (std::size_t j = i + 1; j < m.nodes().size(); ++j)
      if (norm(m.nodes()[i].p - m.nodes()[j].p) <= 0.0) unique = false;
  expect("keine doppelten Knoten", unique);

  // Size ratio between elements sharing a node stays bounded.
  const Real ratio = m.max_neighbour_ratio();
  expect("Elementgroessenverhaeltnis unter der Schranke", ratio <= mesher::kMaxNeighbourRatio);
  std::printf("  %-58s %.4g <= %.2f\n", "  groesstes Verhaeltnis", ratio,
              mesher::kMaxNeighbourRatio);

  // Minimum element count per straight segment is respected.
  int worst_seg = 1 << 30;
  for (std::size_t c = 0; c < g.boundaries().size(); ++c) {
    const auto& pts = g.boundaries()[c].points;
    for (std::size_t s = 0; s + 1 < pts.size(); ++s) {
      int n = 0;
      for (const BoundaryElement& e : m.elements())
        if (e.curve == static_cast<int>(c) && e.segment == static_cast<int>(s)) ++n;
      worst_seg = std::min(worst_seg, n);
    }
  }
  expect("jedes Geometriesegment traegt mindestens vier Elemente",
         worst_seg >= mesher::kMinElementsPerSegment);
  std::printf("  %-58s %d\n", "  kleinste Elementzahl je Segment", worst_seg);
}

// ==========================================================================
// 4.  Closed, correctly oriented region boundaries
// ==========================================================================
static void test_region_loops() {
  std::printf("\n=== Geschlossene Gebietsraender ===\n");
  const DeviceParameters p = example();
  const DeviceGeometry g = DeviceGeometry::build(p);
  const BoundaryMesh m = BoundaryMesh::generate(g);

  // Re-derive each region's meridian area and revolved volume from the mesh
  // alone, by walking its adjacent elements with the region on the left.  This
  // is a closure test, an orientation test and an area test at once: a gap, a
  // flipped element or a missing boundary all break it.
  for (Region r : {Region::Liquid, Region::EmitterSolid, Region::ExtractorSolid,
                   Region::Vacuum}) {
    Real area = 0.0, vol = 0.0;
    std::vector<int> outd(m.nodes().size(), 0), ind(m.nodes().size(), 0);
    for (const BoundaryElement& e : m.elements()) {
      bool flip;
      if (e.is_axis()) {
        if (e.side_a != r) continue;
        flip = false;
      } else if (e.side_a == r) {
        flip = false;
      } else if (e.side_b == r) {
        flip = true;
      } else {
        continue;
      }
      const Vec2 a = flip ? e.b : e.a, b = flip ? e.a : e.b;
      ++outd[flip ? e.node_b : e.node_a];
      ++ind[flip ? e.node_a : e.node_b];
      area += 0.5 * (a.r * b.z - b.r * a.z);
      vol += pi * (a.r * a.r + a.r * b.r + b.r * b.r) / 3.0 * (b.z - a.z);
    }
    bool closed = true;
    for (std::size_t k = 0; k < m.nodes().size(); ++k)
      closed = closed && outd[k] == ind[k] && outd[k] <= 1;

    char buf[128];
    std::snprintf(buf, sizeof buf, "%s: Rand geschlossen und einfach", to_string(r));
    expect(buf, closed);
    std::snprintf(buf, sizeof buf, "%s: Meridianflaeche aus dem Netz", to_string(r));
    check(buf, area, g.region(r).meridian_area(), 1e-13);
    std::snprintf(buf, sizeof buf, "%s: Rotationsvolumen aus dem Netz", to_string(r));
    check(buf, vol, g.region(r).revolved_volume(), 1e-13);
    std::snprintf(buf, sizeof buf, "%s: Umlauf ist mathematisch positiv (CCW)", to_string(r));
    expect(buf, area > 0.0);
  }

  // The four region volumes still fill the domain -- the same conservation the
  // geometry tests pin, now recomputed from mesh elements.
  Real sum = 0.0;
  for (Region r : {Region::Liquid, Region::EmitterSolid, Region::ExtractorSolid, Region::Vacuum})
    sum += g.region(r).revolved_volume();
  check("Summe der Gebietsvolumina = Domaenenvolumen", sum, g.domain_revolved_volume(), 1e-12);
}

// ==========================================================================
// 5.  Determinism, and behaviour on other parameter sets
// ==========================================================================
static void test_determinism_and_variants() {
  std::printf("\n=== Determinismus und Parametervarianten ===\n");
  const DeviceGeometry g = DeviceGeometry::build(example());
  const BoundaryMesh a = BoundaryMesh::generate(g);
  const BoundaryMesh b = BoundaryMesh::generate(g);

  bool identical = a.nodes().size() == b.nodes().size() &&
                   a.elements().size() == b.elements().size();
  if (identical)
    for (std::size_t i = 0; i < a.nodes().size(); ++i)
      identical = identical && a.nodes()[i].p.r == b.nodes()[i].p.r &&
                  a.nodes()[i].p.z == b.nodes()[i].p.z;
  if (identical)
    for (std::size_t i = 0; i < a.elements().size(); ++i)
      identical = identical && a.elements()[i].node_a == b.elements()[i].node_a &&
                  a.elements()[i].node_b == b.elements()[i].node_b &&
                  a.elements()[i].meridian_length == b.elements()[i].meridian_length;
  expect("zwei Laeufe liefern ein bitweise identisches Netz", identical);

  // A second, independently built geometry object must give the same mesh --
  // nothing may depend on object identity or on allocation order.
  const DeviceGeometry g2 = DeviceGeometry::build(example());
  const BoundaryMesh c = BoundaryMesh::generate(g2);
  bool same_geom = c.nodes().size() == a.nodes().size();
  if (same_geom)
    for (std::size_t i = 0; i < a.nodes().size(); ++i)
      same_geom = same_geom && a.nodes()[i].p.r == c.nodes()[i].p.r &&
                  a.nodes()[i].p.z == c.nodes()[i].p.z;
  expect("frisch gebaute Geometrie liefert dasselbe Netz", same_geom);

  // Variants: a cylindrical emitter (phi_1 == phi_3), a finite extractor ring
  // (extractor_outer_radius < domain_radius) and a thin electrode must all mesh
  // and pass every check.  These change the boundary topology, not just sizes.
  struct Variant {
    const char* name;
    DeviceParameters p;
  };
  std::vector<Variant> variants;
  {
    DeviceParameters q = example();
    q.phi_1 = q.phi_3;
    variants.push_back({"zylindrischer Emitter (phi_1 = phi_3)", q});
  }
  {
    DeviceParameters q = example();
    q.extractor_outer_radius = 1.0e-3;
    variants.push_back({"endlicher Extraktorring (Rand + Fase)", q});
  }
  {
    DeviceParameters q = example();
    q.extractor_thickness = 1.0e-5;
    q.extractor_aperture_diameter = 1.0e-4;
    variants.push_back({"duenne Elektrode, kleine Apertur", q});
  }
  {
    DeviceParameters q = example();
    q.phi_2 = 1.0e-6;
    q.phi_1 = 4.0e-6;
    variants.push_back({"sehr feine Bohrung (1 um)", q});
  }
  for (const Variant& v : variants) {
    const DeviceGeometry gv = DeviceGeometry::build(v.p);
    const BoundaryMesh mv = BoundaryMesh::generate(gv);
    const MeshReport rep = mv.validate(gv);
    char buf[160];
    std::snprintf(buf, sizeof buf, "%s: alle Pruefungen bestanden", v.name);
    expect(buf, rep.all_passed());
    if (!rep.all_passed()) rep.print(stdout);
    std::snprintf(buf, sizeof buf, "  %s", v.name);
    std::printf("  %-58s %d Elemente, %d Knoten\n", buf,
                static_cast<int>(mv.elements().size()), static_cast<int>(mv.nodes().size()));
  }
}

// ==========================================================================
// 6.  The mesher's own report, run once on the example set
// ==========================================================================
static void test_report() {
  std::printf("\n=== Selbstpruefung des Vernetzers (Beispielsatz) ===\n");
  const DeviceGeometry g = DeviceGeometry::build(example());
  const BoundaryMesh m = BoundaryMesh::generate(g);
  const MeshReport rep = m.validate(g);
  rep.print(stdout);
  expect("alle Pruefungen des Randvernetzers bestanden", rep.all_passed());

  // The report must actually be able to fail: a report with no checks, or one
  // that reports success unconditionally, would be worthless.
  expect("der Bericht enthaelt Pruefungen", rep.checks.size() >= 10);
  MeshReport broken = rep;
  broken.checks.push_back({"kuenstlicher Fehlschlag", false, ""});
  expect("ein fehlgeschlagener Eintrag laesst den Bericht fehlschlagen", !broken.all_passed());
  check("Fehlerzaehlung", broken.failures(), 1.0, 0.0);
}

// ==========================================================================
// 7.  Naming: the surface at z = 0 must not be sold as a meniscus
// ==========================================================================
static void test_naming() {
  std::printf("\n=== Benennung der Fluessigkeitsoberflaeche ===\n");
  const std::string s = boundary_long_name(BoundaryId::FreeSurfaceReference);
  expect("Langname nennt sie 'anfaenglich' und 'eben'",
         s.find("anfaengliche") != std::string::npos && s.find("eben") != std::string::npos);
  expect("Langname sagt ausdruecklich: noch kein berechneter Meniskus",
         s.find("noch kein berechneter Meniskus") != std::string::npos);
  std::printf("  Langname: %s\n", s.c_str());

  // Nothing else may claim to be a meniscus either.
  bool clean = true;
  for (BoundaryId id : {BoundaryId::SymmetryAxis, BoundaryId::EmitterOuterSurface,
                        BoundaryId::EmitterTipLand, BoundaryId::BoreWall,
                        BoundaryId::LiquidInlet, BoundaryId::ExtractorSurface,
                        BoundaryId::OpenBoundary}) {
    const std::string t = boundary_long_name(id);
    clean = clean && t.find("Meniskus") == std::string::npos;
  }
  expect("kein anderer Rand nennt sich Meniskus", clean);
}

int main() {
  test_size_field();
  test_axisymmetric_contract();
  test_topology();
  test_region_loops();
  test_determinism_and_variants();
  test_report();
  test_naming();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
