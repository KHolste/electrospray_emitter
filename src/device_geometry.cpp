#include "es/device_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::pi;

// ---------------------------------------------------------------------------
// Names
// ---------------------------------------------------------------------------

const char* to_string(Region r) {
  switch (r) {
    case Region::Vacuum: return "vacuum";
    case Region::Liquid: return "liquid";
    case Region::EmitterSolid: return "emitter_solid";
    case Region::ExtractorSolid: return "extractor_solid";
    default: return "outside";
  }
}

const char* to_string(BoundaryId b) {
  switch (b) {
    case BoundaryId::SymmetryAxis: return "symmetry_axis";
    case BoundaryId::EmitterOuterSurface: return "emitter_outer_surface";
    case BoundaryId::EmitterTipLand: return "emitter_tip_land";
    case BoundaryId::BoreWall: return "bore_wall";
    case BoundaryId::FreeSurfaceReference: return "free_surface_reference";
    case BoundaryId::LiquidInlet: return "liquid_inlet";
    case BoundaryId::ExtractorSurface: return "extractor_surface";
    case BoundaryId::NumericalEmitterBackClosure: return "numerical_emitter_back_closure";
    case BoundaryId::OpenBoundary: return "open_boundary";
  }
  return "open_boundary";
}

const char* to_string(FeatureId f) {
  switch (f) {
    case FeatureId::PinnedContactEdge: return "pinned_contact_edge";
    case FeatureId::EmitterOuterEdge: return "emitter_outer_edge";
    case FeatureId::ExtractorApertureEdgeFront: return "extractor_aperture_edge_front";
    case FeatureId::ExtractorApertureEdgeBack: return "extractor_aperture_edge_back";
    case FeatureId::NumericalBackClosureEdge: return "numerical_back_closure_edge";
  }
  return "numerical_back_closure_edge";
}

// ---------------------------------------------------------------------------
// Measures of revolution
// ---------------------------------------------------------------------------

Real meridian_length(const std::vector<Vec2>& poly) {
  Real s = 0.0;
  for (std::size_t i = 1; i < poly.size(); ++i) s += norm(poly[i] - poly[i - 1]);
  return s;
}

Real revolved_area(const std::vector<Vec2>& poly) {
  // 2 pi \int r ds.  Over a straight segment r varies linearly, so the exact
  // contribution is 2 pi * mean radius * segment length -- the truncated-cone
  // lateral area pi (r_a + r_b) * slant.
  Real a = 0.0;
  for (std::size_t i = 1; i < poly.size(); ++i) {
    const Real ds = norm(poly[i] - poly[i - 1]);
    a += pi * (poly[i - 1].r + poly[i].r) * ds;
  }
  return a;
}

Real meridian_signed_area(const std::vector<Vec2>& loop) {
  if (loop.size() < 3) return 0.0;
  Real a = 0.0;
  for (std::size_t i = 0; i < loop.size(); ++i) {
    const Vec2& p = loop[i];
    const Vec2& q = loop[(i + 1) % loop.size()];
    a += p.r * q.z - q.r * p.z;
  }
  return 0.5 * a;
}

Real revolved_volume(const std::vector<Vec2>& loop) {
  // Green's theorem with Q = pi r^2:  oint pi r^2 dz = int_A 2 pi r dr dz.
  // Over a straight segment  \int_0^1 r(t)^2 dt = (r_a^2 + r_a r_b + r_b^2)/3,
  // so the result is exact for a polygonal meridian contour.
  if (loop.size() < 3) return 0.0;
  Real v = 0.0;
  for (std::size_t i = 0; i < loop.size(); ++i) {
    const Vec2& p = loop[i];
    const Vec2& q = loop[(i + 1) % loop.size()];
    v += pi * (p.r * p.r + p.r * q.r + q.r * q.r) / 3.0 * (q.z - p.z);
  }
  return v;
}

Real BoundaryCurve::meridian_length() const { return es::meridian_length(points); }
Real BoundaryCurve::revolved_area() const { return es::revolved_area(points); }

Real RegionBody::meridian_area() const {
  Real a = std::abs(meridian_signed_area(outer_loop));
  for (const auto& h : holes) a -= std::abs(meridian_signed_area(h));
  return a;
}

Real RegionBody::revolved_volume() const {
  Real v = std::abs(es::revolved_volume(outer_loop));
  for (const auto& h : holes) v -= std::abs(es::revolved_volume(h));
  return v;
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

namespace {

void require(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error("DeviceGeometry: " + what);
}

/// Orient a closed meridian loop counter-clockwise, so that revolved_volume()
/// comes out positive and every region is described the same way.
void make_ccw(std::vector<Vec2>& loop) {
  if (meridian_signed_area(loop) < 0.0) std::reverse(loop.begin(), loop.end());
}

}  // namespace

DeviceGeometry DeviceGeometry::build(const DeviceParameters& p) {
  // --- reserved parameters must not be pretended ---------------------------
  const auto& rv = p.reserved;
  if (rv.edge_radius_inner != 0.0 || rv.edge_radius_outer != 0.0 ||
      rv.contact_angle_deg != 0.0 || rv.bore_diameter_at_inlet != 0.0 ||
      rv.porous_emitter || rv.collector_enabled) {
    throw NotImplementedInThisPhase(
        "Reservierte Geometrieparameter (Kantenradius, Kontaktwinkel, verjuengte "
        "Bohrung, poroeser Emitter, Kollektor)",
        "spaetere Phasen -- Kantenradius und Kontaktwinkel in P3, verjuengte Bohrung "
        "und poroeser Emitter danach, Kollektor mit dem Strahlmodell",
        "Die Felder existieren, damit der Parametersatz spaeter nicht umgebaut werden "
        "muss. Umgesetzt ist keiner davon. Ein Wert ungleich dem Standardwert wuerde "
        "eine Geometrie vortaeuschen, die nicht gebaut wird.");
  }

  // --- validation -----------------------------------------------------------
  require(p.phi_2 > 0.0, "phi_2 must be positive");
  require(p.phi_1 > p.phi_2, "need phi_2 < phi_1 (the tip land must have width)");
  require(p.phi_3 >= p.phi_1, "need phi_1 <= phi_3 (the taper must not widen towards the tip)");
  require(p.emitter_height > 0.0, "emitter_height must be positive");
  require(p.extraction_distance > 0.0, "extraction_distance must be positive");
  require(p.extractor_thickness > 0.0, "extractor_thickness must be positive");
  require(p.extractor_aperture_diameter > 0.0, "extractor_aperture_diameter must be positive");
  require(p.domain_radius > 0.5 * p.phi_3,
          "domain_radius must exceed the emitter foot radius");
  require(p.domain_radius > 0.5 * p.extractor_aperture_diameter,
          "domain_radius must exceed the extractor aperture radius");
  require(p.domain_z_min < -p.emitter_height,
          "domain_z_min must lie below the emitter foot");
  // The numerical rearward continuation.  Either it is off, or it is long
  // enough to have a cylindrical shank behind the taper and short enough to
  // stay strictly inside the plotted domain.
  require(p.emitter_back_length >= 0.0,
          "emitter_back_length must not be negative");
  if (p.emitter_back_length > 0.0) {
    require(p.emitter_back_length > p.emitter_height,
            "emitter_back_length must exceed emitter_height, otherwise the closing "
            "disc would cut into the taper instead of closing a cylindrical shank");
    require(p.domain_z_min < -p.emitter_back_length,
            "domain_z_min must lie STRICTLY below the numerical back closure: the closure is a "
            "conductor, the open domain edge is not, and the two must not coincide");
  }
  require(p.domain_z_max > p.extraction_distance + p.extractor_thickness,
          "domain_z_max must lie above the extractor");
  // The electrode's outer radius is a physical dimension of the device and is
  // mandatory.  Equating it with the edge of the computational box would hide a
  // missing dimension behind something that is not a conductor at all.
  require(p.extractor_outer_radius > 0.0,
          "extractor_outer_radius is mandatory and must be positive; 0 no longer means "
          "'out to the domain boundary'");
  require(p.extractor_outer_radius > 0.5 * p.extractor_aperture_diameter,
          "extractor_outer_radius must exceed the aperture radius");
  require(p.domain_radius > p.extractor_outer_radius,
          "domain_radius must be STRICTLY greater than extractor_outer_radius: the electrode "
          "is a conductor, the open domain edge is not, and the two must not coincide");

  DeviceGeometry g;
  g.p_ = p;

  const Real r1 = 0.5 * p.phi_1;   // tip land outer radius
  const Real r2 = 0.5 * p.phi_2;   // bore radius = contact radius
  const Real r3 = 0.5 * p.phi_3;   // foot radius
  const Real H = p.emitter_height;
  const Real zmin = p.domain_z_min;
  const Real zmax = p.domain_z_max;
  const Real R = p.domain_radius;
  const Real ra = 0.5 * p.extractor_aperture_diameter;
  const Real ze = p.extraction_distance;
  const Real zt = ze + p.extractor_thickness;
  const Real rext = p.extractor_outer_radius;  // mandatory, strictly inside the domain

  // The emitter occupies z <= 0 and the extractor z >= extraction_distance > 0,
  // so the two cannot collide.  Assert it anyway rather than assume it.
  require(ze > 0.0, "extractor would intersect the emitter");

  // Rear end of the emitter conductor.  Without the numerical closure it is the
  // domain floor and the conductor is cut open there; with it the conductor
  // ends at z_back, strictly above the floor, and is capped by a disc.
  const bool closed = p.emitter_back_length > 0.0;
  const Real z_back = closed ? -p.emitter_back_length : zmin;

  // --- regions --------------------------------------------------------------
  std::vector<Vec2> emitter{{r2, z_back}, {r2, 0.0}, {r1, 0.0}, {r3, -H}, {r3, z_back}};
  std::vector<Vec2> liquid{{0.0, z_back}, {r2, z_back}, {r2, 0.0}, {0.0, 0.0}};
  std::vector<Vec2> extractor{{ra, ze}, {rext, ze}, {rext, zt}, {ra, zt}};
  std::vector<Vec2> domain{{0.0, zmin}, {R, zmin}, {R, zmax}, {0.0, zmax}};
  make_ccw(emitter);
  make_ccw(liquid);
  make_ccw(extractor);
  make_ccw(domain);

  g.regions_.push_back({Region::EmitterSolid, emitter, {}});
  g.regions_.push_back({Region::Liquid, liquid, {}});
  g.regions_.push_back({Region::ExtractorSolid, extractor, {}});
  g.regions_.push_back({Region::Vacuum, domain, {emitter, liquid, extractor}});

  // --- boundaries -----------------------------------------------------------
  auto add = [&g](BoundaryId id, std::string name, std::vector<Vec2> pts, Region a, Region b) {
    g.boundaries_.push_back({id, std::move(name), std::move(pts), a, b});
  };

  add(BoundaryId::SymmetryAxis, "symmetry_axis.liquid", {{0.0, z_back}, {0.0, 0.0}},
      Region::Liquid, Region::Liquid);
  add(BoundaryId::SymmetryAxis, "symmetry_axis.vacuum", {{0.0, 0.0}, {0.0, zmax}},
      Region::Vacuum, Region::Vacuum);
  if (closed)
    // Behind the closed conductor the axis runs through vacuum again.  It is a
    // symmetry line there like everywhere else, not an interface.
    add(BoundaryId::SymmetryAxis, "symmetry_axis.vacuum_behind_closure",
        {{0.0, zmin}, {0.0, z_back}}, Region::Vacuum, Region::Vacuum);

  add(BoundaryId::FreeSurfaceReference, "free_surface_reference", {{0.0, 0.0}, {r2, 0.0}},
      Region::Liquid, Region::Vacuum);
  add(BoundaryId::BoreWall, "bore_wall", {{r2, 0.0}, {r2, z_back}}, Region::Liquid,
      Region::EmitterSolid);
  if (!closed)
    add(BoundaryId::LiquidInlet, "liquid_inlet", {{0.0, zmin}, {r2, zmin}}, Region::Liquid,
        Region::Outside);

  add(BoundaryId::EmitterTipLand, "emitter_tip_land", {{r2, 0.0}, {r1, 0.0}},
      Region::EmitterSolid, Region::Vacuum);
  add(BoundaryId::EmitterOuterSurface, "emitter_outer_surface",
      {{r1, 0.0}, {r3, -H}, {r3, z_back}}, Region::EmitterSolid, Region::Vacuum);

  if (closed) {
    // The closing disc.  It spans the WHOLE cross section, r = 0 out to the
    // foot radius, because the conductor it closes is the union of the emitter
    // metal and the liquid column: both are at V_emitter in P2a.  It is split
    // at the bore radius only because a boundary curve carries one pair of
    // adjacent regions, not because the two halves differ physically.
    //
    // In the closed configuration there is no liquid inlet: the liquid column
    // is terminated by this cap.  A hydraulic feed boundary belongs to the flow
    // model and returns with it; pretending to have one here, on a surface that
    // is simultaneously a conductor facing vacuum, would be a contradiction.
    add(BoundaryId::NumericalEmitterBackClosure, "numerical_emitter_back_closure.liquid",
        {{0.0, z_back}, {r2, z_back}}, Region::Liquid, Region::Vacuum);
    add(BoundaryId::NumericalEmitterBackClosure, "numerical_emitter_back_closure.solid",
        {{r2, z_back}, {r3, z_back}}, Region::EmitterSolid, Region::Vacuum);
  }

  add(BoundaryId::ExtractorSurface, "extractor_surface.aperture", {{ra, ze}, {ra, zt}},
      Region::ExtractorSolid, Region::Vacuum);
  add(BoundaryId::ExtractorSurface, "extractor_surface.front", {{ra, ze}, {rext, ze}},
      Region::ExtractorSolid, Region::Vacuum);
  add(BoundaryId::ExtractorSurface, "extractor_surface.back", {{ra, zt}, {rext, zt}},
      Region::ExtractorSolid, Region::Vacuum);
  // The electrode is a closed body of revolution: aperture wall, both faces and
  // the outer rim.  The rim exists unconditionally now that rext < R is enforced.
  add(BoundaryId::ExtractorSurface, "extractor_surface.rim", {{rext, ze}, {rext, zt}},
      Region::ExtractorSolid, Region::Vacuum);

  // Outer edges of the open domain.  The pieces are split where a solid or the
  // liquid crosses them, so every piece knows what it touches.
  if (closed) {
    // Nothing reaches the floor any more: the conductor stops at z_back.
    add(BoundaryId::OpenBoundary, "open_boundary.z_min", {{0.0, zmin}, {R, zmin}},
        Region::Vacuum, Region::Outside);
  } else {
    add(BoundaryId::OpenBoundary, "open_boundary.z_min.emitter", {{r2, zmin}, {r3, zmin}},
        Region::EmitterSolid, Region::Outside);
    add(BoundaryId::OpenBoundary, "open_boundary.z_min.vacuum", {{r3, zmin}, {R, zmin}},
        Region::Vacuum, Region::Outside);
  }
  // r = R is now vacuum along its whole length: no solid reaches the box edge.
  add(BoundaryId::OpenBoundary, "open_boundary.r_max", {{R, zmin}, {R, zmax}},
      Region::Vacuum, Region::Outside);
  add(BoundaryId::OpenBoundary, "open_boundary.z_max", {{0.0, zmax}, {R, zmax}},
      Region::Vacuum, Region::Outside);

  // --- features -------------------------------------------------------------
  g.features_.push_back({FeatureId::PinnedContactEdge, {r2, 0.0}});
  g.features_.push_back({FeatureId::EmitterOuterEdge, {r1, 0.0}});
  g.features_.push_back({FeatureId::ExtractorApertureEdgeFront, {ra, ze}});
  g.features_.push_back({FeatureId::ExtractorApertureEdgeBack, {ra, zt}});
  if (closed)
    // Not a device edge.  It is where the numerical continuation is capped, and
    // it is marked so that no field value can be read off it by accident.
    g.features_.push_back({FeatureId::NumericalBackClosureEdge, {r3, z_back}});

  return g;
}

// ---------------------------------------------------------------------------

const RegionBody& DeviceGeometry::region(Region r) const {
  for (const RegionBody& b : regions_)
    if (b.region == r) return b;
  throw std::runtime_error(std::string("DeviceGeometry: no region ") + to_string(r));
}

Vec2 DeviceGeometry::feature(FeatureId f) const {
  for (const NamedFeature& n : features_)
    if (n.id == f) return n.position;
  throw std::runtime_error(std::string("DeviceGeometry: no feature ") + to_string(f));
}

std::vector<const BoundaryCurve*> DeviceGeometry::boundaries_with(BoundaryId id) const {
  std::vector<const BoundaryCurve*> out;
  for (const BoundaryCurve& b : boundaries_)
    if (b.id == id) out.push_back(&b);
  return out;
}

Real DeviceGeometry::cone_half_angle() const {
  return std::atan2(0.5 * (p_.phi_3 - p_.phi_1), p_.emitter_height);
}

Real DeviceGeometry::extractor_outer_radius() const { return p_.extractor_outer_radius; }

Real DeviceGeometry::back_closure_z() const {
  if (!has_back_closure())
    throw std::runtime_error(
        "DeviceGeometry: there is no numerical back closure (emitter_back_length = 0)");
  return -p_.emitter_back_length;
}

Real DeviceGeometry::back_closure_clearance() const {
  return has_back_closure()
             ? evaluation_z_min() - back_closure_z()
             : evaluation_z_min() - p_.domain_z_min;
}

Real DeviceGeometry::domain_revolved_volume() const {
  return pi * p_.domain_radius * p_.domain_radius * (p_.domain_z_max - p_.domain_z_min);
}

void DeviceGeometry::print(std::FILE* out) const {
  std::fprintf(out, "parametrische P1-Geometrie (achsensymmetrisch, r-z, SI)\n");
  std::fprintf(out, "  phi_3 (Fuss)              : %10.4g m\n", p_.phi_3);
  std::fprintf(out, "  phi_1 (Stirnflaeche)      : %10.4g m\n", p_.phi_1);
  std::fprintf(out, "  phi_2 (Bohrung)           : %10.4g m\n", p_.phi_2);
  std::fprintf(out, "  emitter_height            : %10.4g m\n", p_.emitter_height);
  std::fprintf(out, "  extraction_distance       : %10.4g m\n", p_.extraction_distance);
  std::fprintf(out, "  extractor_aperture_diam.  : %10.4g m\n", p_.extractor_aperture_diameter);
  std::fprintf(out, "  extractor_thickness       : %10.4g m\n", p_.extractor_thickness);
  std::fprintf(out, "  extractor_outer_radius    : %10.4g m  (Pflichtangabe, < domain_radius)\n",
               p_.extractor_outer_radius);
  if (has_back_closure())
    std::fprintf(out,
                 "  emitter_back_length: %8.4g m  NUMERISCHE Rueckverlaengerung,\n"
                 "      Abschlussscheibe bei z = %.4g m, Abstand zur ausgewerteten Region\n"
                 "      (z >= %.4g m) = %.4g m.  KEINE physikalische Emitterlaenge.\n",
                 p_.emitter_back_length, back_closure_z(), evaluation_z_min(),
                 back_closure_clearance());
  else
    std::fprintf(out, "  emitter_back_length: 0        -- keine Rueckschliessung, der "
                      "Leiter endet offen am Modellschnitt\n");
  std::fprintf(out, "  Domaene r x z             : %10.4g m x [%.4g, %.4g] m\n", p_.domain_radius,
               p_.domain_z_min, p_.domain_z_max);
  std::fprintf(out, "  abgeleitet: Kegelhalbwinkel %.3f deg, Stirnflaechenbreite %.4g m,\n"
                    "              Kontaktradius %.4g m\n",
               cone_half_angle() * 180.0 / pi, land_width(), contact_radius());

  std::fprintf(out, "\nGebiete (Rotationsvolumen)\n");
  Real sum = 0.0;
  for (const RegionBody& b : regions_) {
    std::fprintf(out, "  %-16s meridian_area = %11.5g m^2   revolved_volume = %11.5g m^3\n",
                 to_string(b.region), b.meridian_area(), b.revolved_volume());
    sum += b.revolved_volume();
  }
  std::fprintf(out, "  %-16s %48.5g m^3\n", "Summe", sum);
  std::fprintf(out, "  %-16s %48.5g m^3   (Abweichung %.3g)\n", "Domaene",
               domain_revolved_volume(),
               std::abs(sum - domain_revolved_volume()) / domain_revolved_volume());

  std::fprintf(out, "\nRaender (Rotationsflaechen)\n");
  for (const BoundaryCurve& b : boundaries_)
    std::fprintf(out, "  %-34s %-22s meridian_length = %10.4g m  revolved_area = %10.4g m^2\n",
                 b.name.c_str(), to_string(b.id), b.meridian_length(), b.revolved_area());

  std::fprintf(out, "\nMerkmale (Kreise)\n");
  for (const NamedFeature& n : features_)
    std::fprintf(out, "  %-30s r = %10.4g m, z = %10.4g m\n", to_string(n.id), n.position.r,
                 n.position.z);
}

void DeviceGeometry::write_csv(const std::string& dir) const {
  const std::string d = dir + "/";

  {
    std::FILE* f = std::fopen((d + "regions.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open regions.csv");
    std::fprintf(f, "# closed meridian loops; loop 0 is the outer contour, >0 are holes\n");
    std::fprintf(f, "region,loop,i,r_m,z_m\n");
    for (const RegionBody& b : regions_) {
      for (std::size_t i = 0; i < b.outer_loop.size(); ++i)
        std::fprintf(f, "%s,0,%zu,%.9e,%.9e\n", to_string(b.region), i, b.outer_loop[i].r,
                     b.outer_loop[i].z);
      for (std::size_t k = 0; k < b.holes.size(); ++k)
        for (std::size_t i = 0; i < b.holes[k].size(); ++i)
          std::fprintf(f, "%s,%zu,%zu,%.9e,%.9e\n", to_string(b.region), k + 1, i,
                       b.holes[k][i].r, b.holes[k][i].z);
    }
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "boundaries.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open boundaries.csv");
    std::fprintf(f, "# labelled boundary polylines in the meridian half-plane\n");
    std::fprintf(f, "name,id,side_a,side_b,i,r_m,z_m\n");
    for (const BoundaryCurve& b : boundaries_)
      for (std::size_t i = 0; i < b.points.size(); ++i)
        std::fprintf(f, "%s,%s,%s,%s,%zu,%.9e,%.9e\n", b.name.c_str(), to_string(b.id),
                     to_string(b.side_a), to_string(b.side_b), i, b.points[i].r, b.points[i].z);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "features.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open features.csv");
    std::fprintf(f, "# zero-dimensional features (circles in 3D) a mesher must resolve\n");
    std::fprintf(f, "name,r_m,z_m\n");
    for (const NamedFeature& n : features_)
      std::fprintf(f, "%s,%.9e,%.9e\n", to_string(n.id), n.position.r, n.position.z);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "parameters.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open parameters.csv");
    std::fprintf(f, "# self-describing parameter set, SI units\n");
    std::fprintf(f, "name,value_SI,unit,role\n");
    std::fprintf(f, "phi_3,%.9e,m,outer diameter at the emitter foot\n", p_.phi_3);
    std::fprintf(f, "phi_1,%.9e,m,outer diameter of the tip land\n", p_.phi_1);
    std::fprintf(f, "phi_2,%.9e,m,diameter of the exit bore\n", p_.phi_2);
    std::fprintf(f, "emitter_height,%.9e,m,axial length of the outer taper\n", p_.emitter_height);
    std::fprintf(f, "extraction_distance,%.9e,m,tip plane to extractor face\n",
                 p_.extraction_distance);
    std::fprintf(f, "extractor_aperture_diameter,%.9e,m,\n", p_.extractor_aperture_diameter);
    std::fprintf(f, "extractor_thickness,%.9e,m,\n", p_.extractor_thickness);
    std::fprintf(f, "extractor_outer_radius,%.9e,m,mandatory; strictly inside domain_radius\n",
                 p_.extractor_outer_radius);
    std::fprintf(f, "emitter_back_length,%.9e,m,NUMERICAL rearward continuation of the "
                    "emitter conductor; 0 = none. Not a physical emitter length\n",
                 p_.emitter_back_length);
    if (has_back_closure()) {
      std::fprintf(f, "back_closure_z,%.9e,m,z of the conducting end cap; derived\n",
                   back_closure_z());
      std::fprintf(f, "evaluation_z_min,%.9e,m,rear limit of the physically evaluated region "
                      "(taper foot); derived\n", evaluation_z_min());
      std::fprintf(f, "back_closure_clearance,%.9e,m,distance from the end cap to that "
                      "region; derived\n", back_closure_clearance());
    }
    std::fprintf(f, "domain_radius,%.9e,m,open computational domain\n", p_.domain_radius);
    std::fprintf(f, "domain_z_min,%.9e,m,open computational domain\n", p_.domain_z_min);
    std::fprintf(f, "domain_z_max,%.9e,m,open computational domain\n", p_.domain_z_max);
    std::fprintf(f, "cone_half_angle,%.9e,rad,derived\n", cone_half_angle());
    std::fprintf(f, "land_width,%.9e,m,derived\n", land_width());
    std::fprintf(f, "contact_radius,%.9e,m,derived\n", contact_radius());
    std::fclose(f);
  }
}

}  // namespace es
