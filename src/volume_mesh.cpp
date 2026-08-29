#include "es/volume_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::pi;
using namespace volume_mesher;

namespace {

void require(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error("build_volume_mesh: " + what);
}

// ---------------------------------------------------------------------------
// One-dimensional size field and node placement
// ---------------------------------------------------------------------------
//
//   h(x) = min( h_max, min_s [ h_s + G |x - x_s| ] )
//
// A minimum of G-Lipschitz functions is G-Lipschitz, so the field bounds its
// own growth rate and no smoothing pass is needed -- the same argument the
// boundary mesher rests on (docs/07_mesher_decision.md, 7.2).
// The mesh level multiplies the FINISHED field, gradation cone included.
// Scaling only the sources and h_max, which is the obvious thing to do, does
// not refine anything: at a distance d from the nearest source the size is
// h_s + G d, and away from the features the G d term dominates and does not
// move with the level at all.  The mid-gap element size then stays at about
// a quarter of the distance to the tip however far the level is pushed, and a
// mesh study measures interpolation jitter instead of convergence.  Scaling the
// whole field also tightens the neighbour ratio bound, to (1 + G s/2)/(1 - G s/2).
struct SizeField1D {
  std::vector<Real> x, h;   ///< UNSCALED source sizes
  Real h_max{0.0};          ///< UNSCALED cap
  Real scale{1.0};

  Real operator()(Real p) const {
    Real v = h_max;
    for (std::size_t k = 0; k < x.size(); ++k)
      v = std::min(v, h[k] + kGradation * std::abs(p - x[k]));
    return scale * v;
  }
};

/// Quadrature breakpoints for \int dx/h over [a, b], refined until every piece
/// is short against the local size.
///
/// A UNIFORM sampling does not work here.  The size field spans four decades
/// inside a single key interval -- 0.09 um at the exit edge against 90 um in the
/// far field -- so a fixed number of samples resolves 1/h nowhere near the
/// refined end, and the first elements of the interval come out several times
/// too large.  That is not a cosmetic defect: it puts a jump of a factor two
/// into the element size exactly at the tip plane, which is where the answer is
/// read off.  The refinement below is local, so the cost follows the number of
/// elements rather than the ratio of the scales.
void refine_breakpoints(Real a, Real b, const SizeField1D& h, std::vector<Real>& pts,
                        int depth) {
  const Real mid = 0.5 * (a + b);
  const Real hm = std::min(std::min(h(a), h(mid)), h(b));
  if (depth < 40 && (b - a) > 0.25 * hm) {
    refine_breakpoints(a, mid, h, pts, depth + 1);
    refine_breakpoints(mid, b, h, pts, depth + 1);
  } else {
    pts.push_back(b);
  }
}

/// Nodes inside [a, b], endpoints included, placed at equal increments of
/// \int dx/h.  The endpoints are always nodes, so every key level survives
/// exactly.  Returns the interior nodes plus b (a is assumed already present).
std::vector<Real> subdivide(Real a, Real b, const SizeField1D& h, Index min_cells) {
  std::vector<Real> x{a};
  refine_breakpoints(a, b, h, x, 0);
  const std::size_t M = x.size() - 1;
  std::vector<Real> s(M + 1, 0.0);
  for (std::size_t k = 1; k <= M; ++k) {
    const Real dx = x[k] - x[k - 1];
    s[k] = s[k - 1] + dx / h(0.5 * (x[k] + x[k - 1]));
  }
  const Real total = s[M];
  Index n = std::max(static_cast<Index>(std::llround(total)), min_cells);

  std::vector<Real> out;
  out.reserve(static_cast<std::size_t>(n));
  std::size_t k = 0;
  for (Index q = 1; q < n; ++q) {
    const Real target = total * static_cast<Real>(q) / static_cast<Real>(n);
    while (k < M && s[k + 1] < target) ++k;
    const Real t = (s[k + 1] > s[k]) ? (target - s[k]) / (s[k + 1] - s[k]) : 0.0;
    out.push_back(x[k] + t * (x[k + 1] - x[k]));
  }
  out.push_back(b);
  return out;
}

/// Local feature size of every entry of a sorted key-level list: the distance
/// to the nearer neighbour.  Purely geometric, no user input.
std::vector<Real> local_feature_sizes(const std::vector<Real>& key) {
  const std::size_t n = key.size();
  std::vector<Real> lfs(n, 0.0);
  for (std::size_t k = 0; k < n; ++k) {
    Real d = 0.0;
    if (k > 0) d = key[k] - key[k - 1];
    if (k + 1 < n) d = (k > 0) ? std::min(d, key[k + 1] - key[k]) : key[k + 1] - key[k];
    lfs[k] = d;
  }
  return lfs;
}

Index index_of(const std::vector<Real>& v, Real x) {
  for (std::size_t k = 0; k < v.size(); ++k)
    if (v[k] == x) return static_cast<Index>(k);
  return -1;
}

}  // namespace

Real mesh_level_scale(int level) { return std::pow(2.0, -0.5 * static_cast<Real>(level)); }

// ---------------------------------------------------------------------------

bool VolumeMeshReport::all_passed() const {
  for (const VolumeMeshCheck& c : checks)
    if (!c.passed) return false;
  return true;
}

void VolumeMeshReport::print(std::FILE* out) const {
  std::fprintf(out, "Volumennetz -- Pruefungen\n");
  for (const VolumeMeshCheck& c : checks)
    std::fprintf(out, "  [%s] %-42s gemessen %11.4g  Grenze %11.4g  %s\n",
                 c.passed ? "ok" : "FEHLER", c.name.c_str(), c.measured, c.bound,
                 c.note.c_str());
}

void VolumeMeshReport::write_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# Pruefliste des Volumenvernetzers\n");
  std::fprintf(f, "check,passed,measured,bound,note\n");
  for (const VolumeMeshCheck& c : checks)
    std::fprintf(f, "%s,%d,%.9e,%.9e,\"%s\"\n", c.name.c_str(), c.passed ? 1 : 0, c.measured,
                 c.bound, c.note.c_str());
  std::fclose(f);
}

// ---------------------------------------------------------------------------

Real DeviceVolumeMesh::emitter_outer_radius_at(Real z) const {
  const Real H = p.device.emitter_height;
  if (z > 0.0 || z < p.liquid_feed_z) return 0.0;
  if (z <= -H) return r_foot;
  return r_land + (r_foot - r_land) * (-z / H);
}

Real DeviceVolumeMesh::warp(Real rr, Real z) const {
  // Identity above the tip plane; below it the reference radius r_land is
  // carried to the local outer radius of the emitter body.  Piecewise linear,
  // continuous, strictly increasing, pinned at r_bore and r_anchor.
  const Real H = p.device.emitter_height;
  const Real zc = std::min(0.0, std::max(-H, z));
  const Real rho = (z >= 0.0) ? r_land : (r_land + (r_foot - r_land) * (-zc / H));
  if (rr <= r_bore || rr >= r_anchor) return rr;
  if (rr <= r_land) return r_bore + (rr - r_bore) * (rho - r_bore) / (r_land - r_bore);
  return rho + (rr - r_land) * (r_anchor - rho) / (r_anchor - r_land);
}

Region DeviceVolumeMesh::region_at(Vec2 x) const {
  const Real ze = p.device.extraction_distance;
  const Real zt = ze + p.device.extractor_thickness;
  if (x.z >= p.liquid_feed_z && x.z <= 0.0) {
    if (x.r <= r_bore) return Region::Liquid;
    if (x.r <= emitter_outer_radius_at(x.z)) return Region::EmitterSolid;
  }
  if (x.z >= ze && x.z <= zt && x.r >= r_aperture && x.r <= r_ext_outer)
    return Region::ExtractorSolid;
  return Region::Vacuum;
}

Index DeviceVolumeMesh::n_cells_of(Region r) const {
  Index n = 0;
  for (Region c : cell_region)
    if (c == r) ++n;
  return n;
}

Real DeviceVolumeMesh::revolved_volume_of(Region r) const {
  Real v = 0.0;
  for (Index j = 0; j + 1 < grid.nz; ++j)
    for (Index i = 0; i + 1 < grid.nr; ++i)
      if (cell_region[static_cast<std::size_t>(grid.cell(i, j))] == r)
        v += grid.cell_revolved_volume(i, j);
  return v;
}

Real DeviceVolumeMesh::analytic_volume_of(Region r) const {
  const Real H = p.device.emitter_height;
  const Real zf = p.liquid_feed_z;
  const Real t = p.device.extractor_thickness;
  switch (r) {
    case Region::Liquid:
      return pi * r_bore * r_bore * (0.0 - zf);
    case Region::EmitterSolid: {
      // Cylindrical shank plus truncated cone, both less the bore.
      const Real shank = pi * (r_foot * r_foot - r_bore * r_bore) * (-H - zf);
      const Real cone = pi * H *
                        ((r_land * r_land + r_land * r_foot + r_foot * r_foot) / 3.0 -
                         r_bore * r_bore);
      return shank + cone;
    }
    case Region::ExtractorSolid:
      return pi * (r_ext_outer * r_ext_outer - r_aperture * r_aperture) * t;
    case Region::Vacuum: {
      const Real total = pi * p.device.domain_radius * p.device.domain_radius *
                         (p.device.domain_z_max - p.device.domain_z_min);
      return total - analytic_volume_of(Region::Liquid) -
             analytic_volume_of(Region::EmitterSolid) -
             analytic_volume_of(Region::ExtractorSolid);
    }
    default:
      return 0.0;
  }
}

// ---------------------------------------------------------------------------

DeviceVolumeMesh build_volume_mesh(const DielectricDeviceParameters& pin) {
  const DeviceParameters& d = pin.device;

  // A dielectric emitter has no conducting rear closure.  Refusing the P2a
  // parameter here is what keeps the superseded model from creeping back in.
  if (d.emitter_back_length != 0.0)
    throw std::runtime_error(
        "build_volume_mesh: device.emitter_back_length ist gesetzt.  Die leitende "
        "Abschlussscheibe gehoert zum metallischen P2a-Emitter und ist fuer den "
        "dielektrischen Emitter physikalisch falsch.  Die Fluessigkeitssaeule wird "
        "stattdessen an liquid_feed_z abgeschnitten; dort gilt die Dirichlet-Bedingung "
        "AUSSCHLIESSLICH auf dem Fluessigkeitsquerschnitt.");
  const DeviceParameters::Reserved& rv = d.reserved;
  if (rv.edge_radius_inner != 0.0 || rv.edge_radius_outer != 0.0 ||
      rv.contact_angle_deg != 0.0 || rv.bore_diameter_at_inlet != 0.0 || rv.porous_emitter ||
      rv.collector_enabled)
    throw NotImplementedInThisPhase(
        "Reservierte Geometrieparameter im dielektrischen Modell",
        "spaetere Phasen (Kantenradius und Kontaktwinkel in P3)",
        "Der Volumenvernetzer bildet keinen davon ab. Ein Wert ungleich dem Standardwert "
        "wuerde eine Geometrie vortaeuschen, die nicht gebaut wird.");

  DeviceVolumeMesh m;
  m.p = pin;
  m.r_bore = 0.5 * d.phi_2;
  m.r_land = 0.5 * d.phi_1;
  m.r_foot = 0.5 * d.phi_3;
  m.r_aperture = 0.5 * d.extractor_aperture_diameter;
  m.r_ext_outer = d.extractor_outer_radius;

  const Real H = d.emitter_height;
  const Real zf = pin.liquid_feed_z;
  const Real ze = d.extraction_distance;
  const Real zt = ze + d.extractor_thickness;
  const Real R = d.domain_radius;
  const Real zmin = d.domain_z_min, zmax = d.domain_z_max;

  require(m.r_bore > 0.0 && m.r_land > m.r_bore && m.r_foot >= m.r_land,
          "need 0 < r_bore < r_land <= r_foot");
  require(H > 0.0 && ze > 0.0 && d.extractor_thickness > 0.0, "H, L, t must be positive");
  require(m.r_aperture > m.r_foot,
          "the extractor aperture radius must exceed the emitter foot radius; the radial warp "
          "pins the grid at the first key radius outside the emitter");
  require(m.r_ext_outer > m.r_aperture, "extractor_outer_radius must exceed the aperture radius");
  require(R > m.r_ext_outer, "domain_radius must strictly exceed extractor_outer_radius");
  require(zf < -H, "liquid_feed_z must lie below the taper foot, so that the modelled body "
                   "ends on a cylindrical shank");
  require(zmin < zf, "domain_z_min must lie strictly below the liquid feed boundary: the feed "
                     "boundary is a device condition, the open domain edge is not");
  require(zmax > zt, "domain_z_max must lie above the extractor");
  require(pin.mesh_level >= 0 && pin.mesh_level <= 8, "mesh_level must be in 0..8");

  // ---- size field ---------------------------------------------------------
  const Real scale = mesh_level_scale(pin.mesh_level);
  const Real diag = std::hypot(R, zmax - zmin);
  const Real h_max = diag / static_cast<Real>(kDomainDivisions);

  // Radial key levels, and axial ones.  liquid_feed_z is deliberately NOT among
  // the levels that FEED THE SIZE FIELD: if it were, moving the feed boundary
  // would change the element size at the tip as well, and the feed study would
  // be measuring two things at once.  It is added to the list that places
  // NODES, further down, so that it is still an exact mesh level.
  std::vector<Real> rkey{0.0, m.r_bore, m.r_land, m.r_foot, m.r_aperture, m.r_ext_outer, R};
  std::vector<Real> zkey{zmin, -H, 0.0, ze, zt, zmax};
  std::sort(rkey.begin(), rkey.end());
  std::sort(zkey.begin(), zkey.end());
  rkey.erase(std::unique(rkey.begin(), rkey.end()), rkey.end());
  zkey.erase(std::unique(zkey.begin(), zkey.end()), zkey.end());
  const std::vector<Real> rlfs = local_feature_sizes(rkey);
  const std::vector<Real> zlfs = local_feature_sizes(zkey);

  auto lfs_at = [](const std::vector<Real>& key, const std::vector<Real>& lfs, Real x) {
    for (std::size_t k = 0; k < key.size(); ++k)
      if (key[k] == x) return lfs[k];
    return 0.0;
  };

  SizeField1D hr, hz;
  hr.h_max = hz.h_max = h_max;
  hr.scale = hz.scale = scale;
  for (std::size_t k = 0; k < rkey.size(); ++k) {
    hr.x.push_back(rkey[k]);
    hr.h.push_back(rlfs[k] / static_cast<Real>(kKeyDivisions));
  }
  for (std::size_t k = 0; k < zkey.size(); ++k) {
    hz.x.push_back(zkey[k]);
    hz.h.push_back(zlfs[k] / static_cast<Real>(kKeyDivisions));
  }
  // Named zero-dimensional features: the unrounded edges the field concentrates
  // at.  Their target size is the 2-D local feature size over kFeatureDivisions,
  // applied to BOTH directions, so the cells there stay close to isotropic.
  const Vec2 features[4] = {{m.r_bore, 0.0},       // pinned contact edge
                            {m.r_land, 0.0},       // outer edge of the tip land
                            {m.r_aperture, ze},    // aperture edge, front
                            {m.r_aperture, zt}};   // aperture edge, back
  for (const Vec2& f : features) {
    const Real lfs2 = std::min(lfs_at(rkey, rlfs, f.r), lfs_at(zkey, zlfs, f.z));
    const Real h = lfs2 / static_cast<Real>(kFeatureDivisions);
    hr.x.push_back(f.r);
    hr.h.push_back(h);
    hz.x.push_back(f.z);
    hz.h.push_back(h);
  }

  // ---- node lists ----------------------------------------------------------
  m.r_ref.push_back(rkey.front());
  for (std::size_t k = 0; k + 1 < rkey.size(); ++k) {
    const std::vector<Real> seg = subdivide(rkey[k], rkey[k + 1], hr, kMinCellsPerInterval);
    m.r_ref.insert(m.r_ref.end(), seg.begin(), seg.end());
  }
  // The node-placement levels are the size-field levels PLUS the feed boundary.
  // Each interval is subdivided independently by the very same size field, so
  // every node at z >= -H comes out bit-identical no matter where the feed
  // boundary sits -- which is what makes the feed study a study of the boundary
  // and not of the mesh -- while the grading across the feed plane is the same
  // graded field as everywhere else.
  std::vector<Real> zplace = zkey;
  zplace.push_back(zf);
  std::sort(zplace.begin(), zplace.end());
  require(std::adjacent_find(zplace.begin(), zplace.end()) == zplace.end(),
          "liquid_feed_z coincides with a geometric level; it must lie strictly between them");
  m.z_ref.push_back(zplace.front());
  for (std::size_t k = 0; k + 1 < zplace.size(); ++k) {
    const std::vector<Real> seg = subdivide(zplace[k], zplace[k + 1], hz, kMinCellsPerInterval);
    m.z_ref.insert(m.z_ref.end(), seg.begin(), seg.end());
  }

  // ---- landmarks -----------------------------------------------------------
  m.r_anchor = m.r_aperture;
  m.i_axis = index_of(m.r_ref, 0.0);
  m.i_bore = index_of(m.r_ref, m.r_bore);
  m.i_land = index_of(m.r_ref, m.r_land);
  m.i_foot = index_of(m.r_ref, m.r_foot);
  m.i_aperture = index_of(m.r_ref, m.r_aperture);
  m.i_ext_outer = index_of(m.r_ref, m.r_ext_outer);
  m.i_far = index_of(m.r_ref, R);
  m.j_min = index_of(m.z_ref, zmin);
  m.j_feed = index_of(m.z_ref, zf);
  m.j_foot = index_of(m.z_ref, -H);
  m.j_tip = index_of(m.z_ref, 0.0);
  m.j_ex_front = index_of(m.z_ref, ze);
  m.j_ex_back = index_of(m.z_ref, zt);
  m.j_max = index_of(m.z_ref, zmax);
  const Index landmarks[] = {m.i_axis,     m.i_bore,     m.i_land,   m.i_foot,
                             m.i_aperture, m.i_ext_outer, m.i_far,   m.j_min,
                             m.j_feed,     m.j_foot,     m.j_tip,    m.j_ex_front,
                             m.j_ex_back,  m.j_max};
  for (Index k : landmarks)
    require(k >= 0, "a key level did not end up as an exact mesh node");

  // ---- nodes ---------------------------------------------------------------
  m.grid.nr = static_cast<Index>(m.r_ref.size());
  m.grid.nz = static_cast<Index>(m.z_ref.size());
  m.grid.nodes.resize(static_cast<std::size_t>(m.grid.n_nodes()));
  for (Index j = 0; j < m.grid.nz; ++j)
    for (Index i = 0; i < m.grid.nr; ++i)
      m.grid.nodes[static_cast<std::size_t>(m.grid.node(i, j))] = {
          m.warp(m.r_ref[static_cast<std::size_t>(i)], m.z_ref[static_cast<std::size_t>(j)]),
          m.z_ref[static_cast<std::size_t>(j)]};

  // ---- cell regions, from indices alone ------------------------------------
  m.cell_region.assign(static_cast<std::size_t>(m.grid.n_cells()), Region::Vacuum);
  for (Index j = 0; j + 1 < m.grid.nz; ++j)
    for (Index i = 0; i + 1 < m.grid.nr; ++i) {
      Region rg = Region::Vacuum;
      if (j >= m.j_feed && j < m.j_tip) {
        if (i < m.i_bore)
          rg = Region::Liquid;
        else if (i < m.i_land)
          rg = Region::EmitterSolid;
      }
      if (j >= m.j_ex_front && j < m.j_ex_back && i >= m.i_aperture && i < m.i_ext_outer)
        rg = Region::ExtractorSolid;
      m.cell_region[static_cast<std::size_t>(m.grid.cell(i, j))] = rg;
    }
  return m;
}

// ---------------------------------------------------------------------------

VolumeMeshReport DeviceVolumeMesh::validate() const {
  VolumeMeshReport rep;
  auto add = [&rep](const char* name, bool ok, Real measured, Real bound, const char* note) {
    rep.checks.push_back({name, ok, measured, bound, note});
  };

  // 1 -- rows lie on constant z (the premise of exact point location)
  add("zeilen_auf_konstantem_z", grid.validate_level_rows(), 0.0, 0.0,
      "Voraussetzung der exakten Punktlokalisierung");

  // 2 -- every cell is convex and positively oriented at every Gauss point
  Real jmin = 0.0;
  bool jac_ok = true;
  try {
    jmin = grid.min_jacobian();
    jac_ok = jmin > 0.0;
  } catch (const std::exception&) {
    jac_ok = false;
  }
  add("jacobi_determinante_positiv", jac_ok, jmin, 0.0, "kein verdrehtes Element");

  // 3 -- strict monotonicity of the node coordinates
  bool mono = true;
  for (Index j = 0; j + 1 < grid.nz; ++j) mono &= (grid.z_of_row(j) < grid.z_of_row(j + 1));
  for (Index j = 0; j < grid.nz; ++j)
    for (Index i = 0; i + 1 < grid.nr; ++i) mono &= (grid.at(i, j).r < grid.at(i + 1, j).r);
  add("knoten_streng_monoton", mono, 0.0, 0.0, "kein ueberschlagenes Gitter");

  // 4 -- the mesh fills the domain exactly
  const Real vol = grid.total_revolved_volume();
  const Real vol_exact = pi * p.device.domain_radius * p.device.domain_radius *
                         (p.device.domain_z_max - p.device.domain_z_min);
  const Real dv = std::abs(vol - vol_exact) / vol_exact;
  add("rotationsvolumen_der_domaene", dv < 1e-12, dv, 1e-12, "Summe der Zellen = Zylinder");

  // 5 -- every region volume is reproduced exactly, which is the real test that
  //      the material interfaces lie ON grid lines and not near them
  for (Region rg : {Region::Liquid, Region::EmitterSolid, Region::ExtractorSolid}) {
    const Real a = analytic_volume_of(rg), b = revolved_volume_of(rg);
    const Real e = (a > 0.0) ? std::abs(a - b) / a : std::abs(a - b);
    add((std::string("gebietsvolumen_") + to_string(rg)).c_str(), e < 1e-11, e, 1e-11,
        "geschlossene Formel gegen Netzsumme");
  }

  // 6 -- index-based region assignment agrees with the closed-form predicate at
  //      the cell centroid.  The two are derived independently on purpose.
  Index mismatched = 0;
  for (Index j = 0; j + 1 < grid.nz; ++j)
    for (Index i = 0; i + 1 < grid.nr; ++i) {
      const auto c = grid.cell_nodes(i, j);
      Vec2 mid{0.0, 0.0};
      for (int k = 0; k < 4; ++k) mid += 0.25 * grid.nodes[static_cast<std::size_t>(c[k])];
      if (region_at(mid) != cell_region[static_cast<std::size_t>(grid.cell(i, j))]) ++mismatched;
    }
  add("gebietszuordnung_geometrisch_geprueft", mismatched == 0,
      static_cast<Real>(mismatched), 0.0, "Indexlogik gegen geschlossene Geometrie");

  // 7 -- neighbouring element sizes stay close
  Real worst = 1.0;
  for (Index i = 1; i + 2 < grid.nr; ++i) {
    const Real a = r_ref[static_cast<std::size_t>(i)] - r_ref[static_cast<std::size_t>(i - 1)];
    const Real b = r_ref[static_cast<std::size_t>(i + 1)] - r_ref[static_cast<std::size_t>(i)];
    worst = std::max(worst, std::max(a / b, b / a));
  }
  for (Index j = 1; j + 2 < grid.nz; ++j) {
    const Real a = z_ref[static_cast<std::size_t>(j)] - z_ref[static_cast<std::size_t>(j - 1)];
    const Real b = z_ref[static_cast<std::size_t>(j + 1)] - z_ref[static_cast<std::size_t>(j)];
    worst = std::max(worst, std::max(a / b, b / a));
  }
  add("elementgroessenverhaeltnis", worst <= kMaxNeighbourRatio, worst, kMaxNeighbourRatio,
      "benachbarte Elementgroessen");

  // 8 -- no cell is empty
  Index degenerate = 0;
  for (Index j = 0; j + 1 < grid.nz; ++j)
    for (Index i = 0; i + 1 < grid.nr; ++i)
      if (!(grid.cell_revolved_volume(i, j) > 0.0)) ++degenerate;
  add("keine_entarteten_zellen", degenerate == 0, static_cast<Real>(degenerate), 0.0,
      "positives Rotationsvolumen ueberall");

  return rep;
}

void DeviceVolumeMesh::print(std::FILE* out) const {
  std::fprintf(out, "Volumennetz (achsensymmetrisch, strukturiert, aus den Geraeteparametern)\n");
  std::fprintf(out, "  Netzstufe                : %d  (Groessenfeldfaktor %.6g)\n", p.mesh_level,
               mesh_level_scale(p.mesh_level));
  std::fprintf(out, "  Knoten                   : %lld x %lld = %lld\n",
               static_cast<long long>(grid.nr), static_cast<long long>(grid.nz),
               static_cast<long long>(grid.n_nodes()));
  std::fprintf(out, "  Zellen                   : %lld\n", static_cast<long long>(grid.n_cells()));
  std::fprintf(out, "  Zulaufgrenze liquid_feed_z: %.6g m\n", p.liquid_feed_z);
  for (Region rg : {Region::Vacuum, Region::Liquid, Region::EmitterSolid, Region::ExtractorSolid})
    std::fprintf(out, "  %-16s %8lld Zellen, Rotationsvolumen %.6g m^3 (Formel %.6g)\n",
                 to_string(rg), static_cast<long long>(n_cells_of(rg)), revolved_volume_of(rg),
                 analytic_volume_of(rg));
}

void DeviceVolumeMesh::write_csv(const std::string& dir) const {
  const std::string d = dir + "/";
  {
    std::FILE* f = std::fopen((d + "mesh_nodes.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open mesh_nodes.csv");
    std::fprintf(f, "# strukturiertes Volumennetz: Knoten (i,j) -> (r,z)\n");
    std::fprintf(f, "# nr=%lld nz=%lld\n", static_cast<long long>(grid.nr),
                 static_cast<long long>(grid.nz));
    std::fprintf(f, "i,j,r_m,z_m\n");
    for (Index j = 0; j < grid.nz; ++j)
      for (Index i = 0; i < grid.nr; ++i)
        std::fprintf(f, "%lld,%lld,%.9e,%.9e\n", static_cast<long long>(i),
                     static_cast<long long>(j), grid.at(i, j).r, grid.at(i, j).z);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "mesh_cells.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open mesh_cells.csv");
    std::fprintf(f, "# Materialgebiet je Zelle\n");
    std::fprintf(f, "i,j,region\n");
    for (Index j = 0; j + 1 < grid.nz; ++j)
      for (Index i = 0; i + 1 < grid.nr; ++i)
        std::fprintf(f, "%lld,%lld,%s\n", static_cast<long long>(i), static_cast<long long>(j),
                     to_string(cell_region[static_cast<std::size_t>(grid.cell(i, j))]));
    std::fclose(f);
  }
  {
    // Closed-form outlines of the P2b body, and the named boundary curves with
    // the electrical role each one plays.  Written here rather than taken from
    // DeviceGeometry, because the P1 geometry runs the emitter down to the
    // domain floor while P2b ends it at the feed boundary -- drawing the P1
    // outline over a P2b field would show a body that was never solved.
    const Real H = p.device.emitter_height;
    const Real zf = p.liquid_feed_z;
    const Real ze = p.device.extraction_distance;
    const Real zt = ze + p.device.extractor_thickness;
    std::FILE* f = std::fopen((d + "device_outline.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open device_outline.csv");
    std::fprintf(f, "# geschlossene Meridianumrisse der P2b-Gebiete und die benannten "
                    "Randkurven\n");
    std::fprintf(f, "kind,name,i,r_m,z_m\n");
    auto poly = [&](const char* kind, const char* name, std::vector<Vec2> pts) {
      for (std::size_t k = 0; k < pts.size(); ++k)
        std::fprintf(f, "%s,%s,%zu,%.9e,%.9e\n", kind, name, k, pts[k].r, pts[k].z);
    };
    poly("region", "liquid", {{0.0, zf}, {r_bore, zf}, {r_bore, 0.0}, {0.0, 0.0}});
    poly("region", "emitter_dielectric",
         {{r_bore, zf}, {r_foot, zf}, {r_foot, -H}, {r_land, 0.0}, {r_bore, 0.0}});
    poly("region", "extractor_carrier",
         {{r_aperture, ze}, {r_ext_outer, ze}, {r_ext_outer, zt}, {r_aperture, zt}});
    poly("boundary", "liquid_feed_boundary", {{0.0, zf}, {r_bore, zf}});
    poly("boundary", "free_surface_reference", {{0.0, 0.0}, {r_bore, 0.0}});
    poly("boundary", "bore_wall", {{r_bore, zf}, {r_bore, 0.0}});
    poly("boundary", "emitter_tip_land", {{r_bore, 0.0}, {r_land, 0.0}});
    poly("boundary", "emitter_outer_surface", {{r_land, 0.0}, {r_foot, -H}, {r_foot, zf}});
    poly("boundary", "emitter_rear_face", {{r_bore, zf}, {r_foot, zf}});
    poly("boundary", "extractor_front", {{r_aperture, ze}, {r_ext_outer, ze}});
    poly("boundary", "extractor_aperture_wall", {{r_aperture, ze}, {r_aperture, zt}});
    poly("boundary", "extractor_back", {{r_aperture, zt}, {r_ext_outer, zt}});
    poly("boundary", "extractor_rim", {{r_ext_outer, ze}, {r_ext_outer, zt}});
    poly("feature", "pinned_contact_edge", {{r_bore, 0.0}});
    poly("feature", "emitter_outer_edge", {{r_land, 0.0}});
    poly("feature", "extractor_aperture_edge_front", {{r_aperture, ze}});
    poly("feature", "extractor_aperture_edge_back", {{r_aperture, zt}});
    poly("feature", "liquid_feed_rim", {{r_bore, zf}});
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "parameters.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open parameters.csv");
    std::fprintf(f, "# selbstbeschreibender Parametersatz des P2b-Modells, SI\n");
    std::fprintf(f, "name,value_SI,unit,role\n");
    std::fprintf(f, "phi_3,%.9e,m,Aussendurchmesser am Emitterfuss\n", p.device.phi_3);
    std::fprintf(f, "phi_1,%.9e,m,Aussendurchmesser der Stirnflaeche\n", p.device.phi_1);
    std::fprintf(f, "phi_2,%.9e,m,Bohrungsdurchmesser\n", p.device.phi_2);
    std::fprintf(f, "emitter_height,%.9e,m,axiale Laenge der Verjuengung\n",
                 p.device.emitter_height);
    std::fprintf(f, "extraction_distance,%.9e,m,Stirnebene bis Extraktor\n",
                 p.device.extraction_distance);
    std::fprintf(f, "extractor_aperture_diameter,%.9e,m,\n",
                 p.device.extractor_aperture_diameter);
    std::fprintf(f, "extractor_thickness,%.9e,m,\n", p.device.extractor_thickness);
    std::fprintf(f, "extractor_outer_radius,%.9e,m,Pflichtangabe\n",
                 p.device.extractor_outer_radius);
    std::fprintf(f, "liquid_feed_z,%.9e,m,Lage der Zulaufgrenze; BEISPIELWERT, keine "
                    "gemessene Abmessung\n", p.liquid_feed_z);
    std::fprintf(f, "domain_radius,%.9e,m,offene Rechendomaene\n", p.device.domain_radius);
    std::fprintf(f, "domain_z_min,%.9e,m,offene Rechendomaene\n", p.device.domain_z_min);
    std::fprintf(f, "domain_z_max,%.9e,m,offene Rechendomaene\n", p.device.domain_z_max);
    std::fprintf(f, "r_bore,%.9e,m,abgeleitet\n", r_bore);
    std::fprintf(f, "r_land,%.9e,m,abgeleitet\n", r_land);
    std::fprintf(f, "r_foot,%.9e,m,abgeleitet\n", r_foot);
    std::fprintf(f, "r_aperture,%.9e,m,abgeleitet\n", r_aperture);
    std::fprintf(f, "mesh_level,%d,-,Netzstufe\n", p.mesh_level);
    std::fprintf(f, "mesh_size_scale,%.9e,-,Groessenfeldfaktor 2^(-Stufe/2)\n",
                 mesh_level_scale(p.mesh_level));
    std::fclose(f);
  }
  validate().write_csv(d + "mesh_checks.csv");
}

}  // namespace es
