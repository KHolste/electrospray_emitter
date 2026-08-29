#include <algorithm>
// P1 device geometry: measures of revolution, boundary topology, validation.
//
// The point of the analytic tests is that a planar-2D mistake cannot pass them:
// every expected value carries the 2*pi*r Jacobian, and a planar formula is off
// by exactly that factor.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/device_geometry.hpp"
#include "es/status.hpp"

using namespace es;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-52s got=%-14.7g want=%-14.7g %s\n", what, got, want, ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-52s %s\n", what, ok ? "OK" : "FAIL");
}

// ==========================================================================
// 1.  Measures of revolution against closed-form solids
// ==========================================================================
static void test_revolved_measures() {
  std::printf("\n=== Rotationsmasse gegen analytische Koerper ===\n");
  const Real R = 7.0e-5, H = 2.3e-4;

  // --- cylinder -----------------------------------------------------------
  {
    const std::vector<Vec2> loop{{0, 0}, {R, 0}, {R, H}, {0, H}};
    check("Zylinder: Volumen = pi R^2 H", revolved_volume(loop), pi * R * R * H, 1e-14);
    const std::vector<Vec2> lateral{{R, 0}, {R, H}};
    check("Zylinder: Mantelflaeche = 2 pi R H", revolved_area(lateral), 2 * pi * R * H, 1e-14);
    const std::vector<Vec2> cap{{0, H}, {R, H}};
    check("Zylinder: Deckflaeche = pi R^2", revolved_area(cap), pi * R * R, 1e-14);

    // A planar reading would give the meridian area R*H and the meridian length
    // H.  Both must be clearly different, or the axisymmetric factor is missing.
    const Real planar_area = R * H;
    expect("Meridianflaeche != Rotationsvolumen (2 pi r fehlt nicht)",
           std::abs(planar_area - revolved_volume(loop)) > 0.5 * planar_area);
    check("Meridianlaenge der Mantellinie", meridian_length(lateral), H, 1e-14);
    check("Verhaeltnis Mantelflaeche / Meridianlaenge = 2 pi R",
          revolved_area(lateral) / meridian_length(lateral), 2 * pi * R, 1e-13);
  }

  // --- truncated cone -----------------------------------------------------
  {
    const Real R1 = 9.0e-5, R2 = 3.0e-5, Hc = 1.7e-4;
    const std::vector<Vec2> loop{{0, 0}, {R1, 0}, {R2, Hc}, {0, Hc}};
    check("Kegelstumpf: V = pi H (R1^2 + R1 R2 + R2^2)/3", revolved_volume(loop),
          pi * Hc * (R1 * R1 + R1 * R2 + R2 * R2) / 3.0, 1e-14);
    const std::vector<Vec2> flank{{R1, 0}, {R2, Hc}};
    const Real slant = std::sqrt((R1 - R2) * (R1 - R2) + Hc * Hc);
    check("Kegelstumpf: Mantel = pi (R1 + R2) s", revolved_area(flank),
          pi * (R1 + R2) * slant, 1e-14);
  }

  // --- full cone as the degenerate case -----------------------------------
  {
    const Real Rc = 5.0e-5, Hc = 1.1e-4;
    const std::vector<Vec2> loop{{0, 0}, {Rc, 0}, {0, Hc}};
    check("Kegel: V = pi R^2 H / 3", revolved_volume(loop), pi * Rc * Rc * Hc / 3.0, 1e-14);
    const std::vector<Vec2> flank{{Rc, 0}, {0, Hc}};
    check("Kegel: Mantel = pi R s", revolved_area(flank),
          pi * Rc * std::sqrt(Rc * Rc + Hc * Hc), 1e-14);
  }

  // --- orientation --------------------------------------------------------
  {
    std::vector<Vec2> ccw{{0, 0}, {R, 0}, {R, H}, {0, H}};
    std::vector<Vec2> cw(ccw.rbegin(), ccw.rend());
    expect("Meridianflaeche ist fuer CCW positiv", meridian_signed_area(ccw) > 0);
    expect("und fuer CW negativ", meridian_signed_area(cw) < 0);
    check("Betrag des Volumens ist orientierungsunabhaengig",
          std::abs(revolved_volume(cw)), std::abs(revolved_volume(ccw)), 1e-14);
  }
}

// ==========================================================================
// 2.  The device: regions, volumes, closure
// ==========================================================================
static DeviceParameters example() {
  DeviceParameters p;
  p.phi_3 = 4.0e-5;
  p.phi_1 = 2.0e-5;
  p.phi_2 = 1.0e-5;
  p.emitter_height = 6.0e-5;
  p.extraction_distance = 5.0e-4;
  p.extractor_aperture_diameter = 4.0e-4;
  p.extractor_thickness = 1.0e-4;
  p.domain_radius = 3.0e-3;
  p.domain_z_min = -4.0e-4;
  p.domain_z_max = 1.5e-3;
  return p;
}

static void test_device_regions() {
  std::printf("\n=== Geraetegeometrie: Gebiete und Volumenbilanz ===\n");
  const DeviceParameters p = example();
  const DeviceGeometry g = DeviceGeometry::build(p);

  const Real r1 = 0.5 * p.phi_1, r2 = 0.5 * p.phi_2, r3 = 0.5 * p.phi_3;
  const Real H = p.emitter_height, zmin = p.domain_z_min;

  // Liquid: a cylinder of radius r2 from z_min to 0.
  check("Fluessigkeitsvolumen = pi r2^2 (0 - z_min)",
        g.region(Region::Liquid).revolved_volume(), pi * r2 * r2 * (0.0 - zmin), 1e-13);

  // Emitter: shank ring plus frustum ring, both hollow with bore radius r2.
  const Real v_shank = pi * (r3 * r3 - r2 * r2) * (-H - zmin);
  const Real v_frustum = pi * H * (r3 * r3 + r3 * r1 + r1 * r1) / 3.0 - pi * r2 * r2 * H;
  check("Emittervolumen = Schaft + Kegelstumpf, jeweils abzueglich Bohrung",
        g.region(Region::EmitterSolid).revolved_volume(), v_shank + v_frustum, 1e-12);

  // Extractor: annulus.
  const Real rext = g.extractor_outer_radius();
  check("Extraktorvolumen = pi (r_ext^2 - r_a^2) t",
        g.region(Region::ExtractorSolid).revolved_volume(),
        pi * (rext * rext - 0.25 * p.extractor_aperture_diameter *
                                 p.extractor_aperture_diameter) *
            p.extractor_thickness, 1e-13);

  // Closure: the four regions must tile the domain exactly.
  Real sum = 0.0;
  for (const RegionBody& b : g.regions()) sum += b.revolved_volume();
  check("Summe aller Gebietsvolumina = Domaenenvolumen", sum, g.domain_revolved_volume(), 1e-12);

  expect("alle Volumina sind positiv",
         g.region(Region::Vacuum).revolved_volume() > 0 &&
             g.region(Region::Liquid).revolved_volume() > 0 &&
             g.region(Region::EmitterSolid).revolved_volume() > 0 &&
             g.region(Region::ExtractorSolid).revolved_volume() > 0);

  // Derived quantities.
  check("Kegelhalbwinkel = atan((r3 - r1)/H)", g.cone_half_angle(),
        std::atan2(r3 - r1, H), 1e-14);
  check("Stirnflaechenbreite = r1 - r2", g.land_width(), r1 - r2, 1e-14);
  check("Kontaktradius = phi_2/2", g.contact_radius(), r2, 1e-14);
}

// ==========================================================================
// 3.  Boundary identifiers and topology
// ==========================================================================
static void test_boundaries() {
  std::printf("\n=== Randkennungen und Topologie ===\n");
  const DeviceParameters p = example();
  const DeviceGeometry g = DeviceGeometry::build(p);

  for (BoundaryId id : {BoundaryId::SymmetryAxis, BoundaryId::EmitterOuterSurface,
                        BoundaryId::EmitterTipLand, BoundaryId::BoreWall,
                        BoundaryId::FreeSurfaceReference, BoundaryId::LiquidInlet,
                        BoundaryId::ExtractorSurface, BoundaryId::OpenBoundary}) {
    char buf[96];
    std::snprintf(buf, sizeof buf, "Rand '%s' ist vorhanden", to_string(id));
    expect(buf, !g.boundaries_with(id).empty());
  }

  // Every boundary name must be unique -- otherwise a mesher cannot address them.
  {
    std::vector<std::string> names;
    for (const BoundaryCurve& b : g.boundaries()) names.push_back(b.name);
    std::sort(names.begin(), names.end());
    expect("alle Randnamen sind eindeutig",
           std::adjacent_find(names.begin(), names.end()) == names.end());
  }

  // The axis really is at r = 0 and covers the whole domain height.
  {
    Real axis_len = 0.0;
    bool on_axis = true;
    for (const BoundaryCurve* b : g.boundaries_with(BoundaryId::SymmetryAxis)) {
      axis_len += b->meridian_length();
      for (const Vec2& q : b->points) on_axis = on_axis && (q.r == 0.0);
      // A curve on the axis has zero surface of revolution.
      check("Rotationsflaeche der Achse ist null", b->revolved_area() + 1.0, 1.0, 1e-12);
    }
    expect("Achse liegt exakt bei r = 0", on_axis);
    check("Achse ueberspannt die Domaenenhoehe", axis_len,
          p.domain_z_max - p.domain_z_min, 1e-13);
  }

  // Interfaces separate the regions they claim to separate.
  {
    const BoundaryCurve* bore = g.boundaries_with(BoundaryId::BoreWall).front();
    expect("Bohrungswand trennt Fluessigkeit und Emitter",
           (bore->side_a == Region::Liquid && bore->side_b == Region::EmitterSolid));
    check("Bohrungswand: Rotationsflaeche = 2 pi r2 L", bore->revolved_area(),
          2 * pi * 0.5 * p.phi_2 * (0.0 - p.domain_z_min), 1e-13);

    const BoundaryCurve* land = g.boundaries_with(BoundaryId::EmitterTipLand).front();
    check("Stirnflaeche = pi (r1^2 - r2^2)", land->revolved_area(),
          pi * (0.25 * p.phi_1 * p.phi_1 - 0.25 * p.phi_2 * p.phi_2), 1e-13);

    const BoundaryCurve* fs = g.boundaries_with(BoundaryId::FreeSurfaceReference).front();
    check("Referenzflaeche des Flu.-Austritts = pi r2^2", fs->revolved_area(),
          pi * 0.25 * p.phi_2 * p.phi_2, 1e-13);
    expect("Referenzflaeche trennt Fluessigkeit und Vakuum",
           fs->side_a == Region::Liquid && fs->side_b == Region::Vacuum);
  }

  // Features sit where the drawing says they do.
  check("gepinnte Austrittskante bei r = phi_2/2",
        g.feature(FeatureId::PinnedContactEdge).r, 0.5 * p.phi_2, 1e-14);
  check("und bei z = 0", g.feature(FeatureId::PinnedContactEdge).z + 1.0, 1.0, 1e-14);
  check("aeussere Stirnkante bei r = phi_1/2", g.feature(FeatureId::EmitterOuterEdge).r,
        0.5 * p.phi_1, 1e-14);
  check("Aperturkante vorne bei z = extraction_distance",
        g.feature(FeatureId::ExtractorApertureEdgeFront).z, p.extraction_distance, 1e-14);
  check("Aperturkante hinten bei z = extraction_distance + t",
        g.feature(FeatureId::ExtractorApertureEdgeBack).z,
        p.extraction_distance + p.extractor_thickness, 1e-14);
}

// ==========================================================================
// 4.  Validation and fail-closed behaviour
// ==========================================================================
static bool rejects(DeviceParameters p) {
  try {
    (void)DeviceGeometry::build(p);
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

static void test_validation() {
  std::printf("\n=== Validierung der Eingaben ===\n");
  expect("gueltiger Satz wird akzeptiert", !rejects(example()));

  {
    DeviceParameters p = example(); p.phi_2 = 0.0;
    expect("phi_2 = 0 wird abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example(); p.phi_2 = p.phi_1;
    expect("phi_2 = phi_1 wird abgelehnt (Stirnflaeche ohne Breite)", rejects(p));
  }
  {
    DeviceParameters p = example(); p.phi_1 = p.phi_3 * 1.5;
    expect("phi_1 > phi_3 wird abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example(); p.phi_1 = p.phi_3;
    expect("phi_1 = phi_3 ist zulaessig (zylindrischer Emitter)", !rejects(p));
  }
  {
    DeviceParameters p = example(); p.emitter_height = -1e-6;
    expect("negative emitter_height wird abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example(); p.extraction_distance = 0.0;
    expect("extraction_distance = 0 wird abgelehnt (Kollision)", rejects(p));
  }
  {
    DeviceParameters p = example(); p.domain_z_min = -0.5 * p.emitter_height;
    expect("Domaene endet ueber dem Emitterfuss -> abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example(); p.domain_z_max = p.extraction_distance;
    expect("Domaene endet vor dem Extraktor -> abgelehnt", rejects(p));
  }
  {
    DeviceParameters p = example(); p.domain_radius = 0.4 * p.phi_3;
    expect("Domaenenradius kleiner als der Emitterfuss -> abgelehnt", rejects(p));
  }

  // Reserved parameters must fail closed, not silently do nothing.
  const char* names[] = {"edge_radius_inner", "edge_radius_outer", "contact_angle_deg",
                         "bore_diameter_at_inlet", "porous_emitter", "collector_enabled"};
  for (int k = 0; k < 6; ++k) {
    DeviceParameters p = example();
    switch (k) {
      case 0: p.reserved.edge_radius_inner = 1e-7; break;
      case 1: p.reserved.edge_radius_outer = 1e-7; break;
      case 2: p.reserved.contact_angle_deg = 30.0; break;
      case 3: p.reserved.bore_diameter_at_inlet = 2e-5; break;
      case 4: p.reserved.porous_emitter = true; break;
      default: p.reserved.collector_enabled = true; break;
    }
    bool threw = false;
    std::string msg;
    try { (void)DeviceGeometry::build(p); }
    catch (const NotImplementedInThisPhase& e) { threw = true; msg = e.what(); }
    catch (const std::exception&) {}
    char buf[96];
    std::snprintf(buf, sizeof buf, "reserviert '%s' schlaegt geschlossen fehl", names[k]);
    expect(buf, threw);
    if (k == 0)
      expect("Meldung nennt die spaetere Phase", msg.find("P3") != std::string::npos);
  }
}

int main() {
  test_revolved_measures();
  test_device_regions();
  test_boundaries();
  test_validation();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
