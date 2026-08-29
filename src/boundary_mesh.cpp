#include "es/boundary_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "es/constants.hpp"

namespace es {

using constants::pi;

// ---------------------------------------------------------------------------
// Names
// ---------------------------------------------------------------------------

const char* to_string(ElementKind k) {
  return k == ElementKind::AxisSymmetry ? "axis_symmetry" : "ring";
}

const char* boundary_long_name(BoundaryId b) {
  switch (b) {
    case BoundaryId::SymmetryAxis:
      return "Symmetrieachse r = 0 (kein Interface)";
    case BoundaryId::EmitterOuterSurface:
      return "Emitteraussenflaeche (Kegelstumpf und Schaft)";
    case BoundaryId::EmitterTipLand:
      return "Stirnflaeche des Emitters";
    case BoundaryId::BoreWall:
      return "Bohrungswand (Fluessigkeit gegen Emitterfestkoerper)";
    case BoundaryId::FreeSurfaceReference:
      return "anfaengliche ebene Fluessigkeitsoberflaeche "
             "- noch kein berechneter Meniskus";
    case BoundaryId::LiquidInlet:
      return "Zulaufschnitt durch die Fluessigkeitssaeule";
    case BoundaryId::ExtractorSurface:
      return "Extraktorflaechen und Aperturwand";
    default:
      return "offener Domaenenrand";
  }
}

// ---------------------------------------------------------------------------
// Small geometric helpers
// ---------------------------------------------------------------------------

namespace {

Real cross2(Vec2 u, Vec2 v) { return u.r * v.z - u.z * v.r; }

/// Distance from x to the segment [p,q].  Handles the degenerate p == q.
Real dist_point_segment(Vec2 x, Vec2 p, Vec2 q) {
  const Vec2 d = q - p;
  const Real dd = norm2(d);
  if (dd <= 0.0) return norm(x - p);
  Real t = dot(x - p, d) / dd;
  t = std::min(std::max(t, 0.0), 1.0);
  return norm(x - (p + t * d));
}

/// Even-odd crossing test for a closed meridian polygon.  Points exactly on the
/// boundary are undefined by construction; every caller offsets away from the
/// boundary before asking.
bool point_in_loop(const std::vector<Vec2>& loop, Vec2 x) {
  bool inside = false;
  const std::size_t n = loop.size();
  for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
    const Vec2& a = loop[i];
    const Vec2& b = loop[j];
    if ((a.z > x.z) != (b.z > x.z)) {
      const Real rr = a.r + (x.z - a.z) / (b.z - a.z) * (b.r - a.r);
      if (x.r < rr) inside = !inside;
    }
  }
  return inside;
}

bool point_in_region_body(const RegionBody& body, Vec2 x) {
  if (!point_in_loop(body.outer_loop, x)) return false;
  for (const auto& h : body.holes)
    if (point_in_loop(h, x)) return false;
  return true;
}

/// Which material occupies x.  Region::Outside if x lies outside the domain.
/// The solid and liquid bodies are disjoint and are tested first; the vacuum
/// body is the domain rectangle with those three as holes, so a point inside a
/// solid can never be reported as vacuum.
Region region_at(const DeviceGeometry& g, Vec2 x) {
  for (Region r : {Region::Liquid, Region::EmitterSolid, Region::ExtractorSolid})
    if (point_in_region_body(g.region(r), x)) return r;
  if (point_in_region_body(g.region(Region::Vacuum), x)) return Region::Vacuum;
  return Region::Outside;
}

Real domain_diagonal(const DeviceParameters& p) {
  return std::sqrt(sqr(p.domain_radius) + sqr(p.domain_z_max - p.domain_z_min));
}

/// All straight geometric segments of every boundary curve.
struct Seg {
  Vec2 p, q;
  int curve;
  int index;  ///< segment index inside its curve
  Real len;
};

std::vector<Seg> all_segments(const DeviceGeometry& g) {
  std::vector<Seg> out;
  const auto& curves = g.boundaries();
  for (std::size_t c = 0; c < curves.size(); ++c) {
    const auto& pts = curves[c].points;
    for (std::size_t i = 1; i < pts.size(); ++i)
      out.push_back({pts[i - 1], pts[i], static_cast<int>(c), static_cast<int>(i - 1),
                     norm(pts[i] - pts[i - 1])});
  }
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Size field
// ---------------------------------------------------------------------------

SizeField SizeField::from_geometry(const DeviceGeometry& g) {
  const DeviceParameters& p = g.parameters();
  const Real diag = domain_diagonal(p);
  const Real snap = mesher::kNodeSnapRelative * diag;
  const std::vector<Seg> segs = all_segments(g);

  auto same = [snap](Vec2 u, Vec2 v) { return norm(u - v) <= snap; };

  // Local feature size at a polyline vertex: the smaller of the shortest
  // segment meeting there and the distance to the nearest segment that does
  // not.  Both are properties of the geometry alone.
  auto local_feature_size = [&](Vec2 x) {
    Real incident = std::numeric_limits<Real>::infinity();
    Real separate = std::numeric_limits<Real>::infinity();
    for (const Seg& s : segs) {
      if (same(s.p, x) || same(s.q, x))
        incident = std::min(incident, s.len);
      else
        separate = std::min(separate, dist_point_segment(x, s.p, s.q));
    }
    return std::min(incident, separate);
  };

  // Collect the distinct vertices, in the order the geometry lists them, so
  // the source list is deterministic.
  std::vector<Vec2> vertices;
  for (const BoundaryCurve& c : g.boundaries())
    for (const Vec2& v : c.points) {
      bool seen = false;
      for (const Vec2& w : vertices) seen = seen || same(v, w);
      if (!seen) vertices.push_back(v);
    }

  SizeField f;
  for (const Vec2& v : vertices) {
    const Real lfs = local_feature_size(v);
    if (!(lfs > 0.0) || !std::isfinite(lfs))
      throw std::runtime_error("BoundaryMesh: degenerate local feature size at a vertex");

    int feature = -1;
    for (std::size_t k = 0; k < g.features().size(); ++k)
      if (same(g.features()[k].position, v)) feature = static_cast<int>(k);

    Source s;
    s.x = v;
    s.local_feature_size = lfs;
    s.is_named_feature = (feature >= 0);
    if (feature >= 0) {
      s.h = lfs / mesher::kFeatureDivisions;
      s.origin = std::string("Merkmal ") + to_string(g.features()[feature].id);
    } else {
      s.h = lfs / mesher::kCornerDivisions;
      s.origin = "Geometrieecke";
    }
    f.sources_.push_back(std::move(s));
  }

  f.h_max_ = diag / mesher::kDomainDivisions;
  f.h_min_ = f.h_max_;
  for (const Source& s : f.sources_) f.h_min_ = std::min(f.h_min_, s.h);
  if (!(f.h_min_ > 0.0))
    throw std::runtime_error("BoundaryMesh: size field has a non-positive floor");
  return f;
}

Real SizeField::operator()(Vec2 x) const {
  Real h = h_max_;
  for (const Source& s : sources_) h = std::min(h, s.h + mesher::kGradation * norm(x - s.x));
  return h;
}

void SizeField::print(std::FILE* out) const {
  std::fprintf(out,
               "Groessenfunktion  h(x) = min( min_s [ h_s + G |x - x_s| ], h_max )\n"
               "  G (Gradation)                 = %.4g\n"
               "  h_max = Diagonale / %.0f       = %.5g m\n"
               "  kleinste Quellgroesse h_min   = %.5g m\n"
               "  Merkmal:  h = lfs / %.0f     Ecke: h = lfs / %.0f\n"
               "  Quellen (%zu):\n",
               mesher::kGradation, mesher::kDomainDivisions, h_max_, h_min_,
               mesher::kFeatureDivisions, mesher::kCornerDivisions, sources_.size());
  for (const Source& s : sources_)
    std::fprintf(out, "    r=%11.5g z=%11.5g  lfs=%10.4g  h=%10.4g  %s%s\n", s.x.r, s.x.z,
                 s.local_feature_size, s.h, s.origin.c_str(),
                 s.is_named_feature ? "  [verfeinert]" : "");
}

// ---------------------------------------------------------------------------
// Node registry
// ---------------------------------------------------------------------------

namespace {

class NodeRegistry {
 public:
  NodeRegistry(Real snap, std::vector<MeshNode>* nodes) : snap_(snap), nodes_(nodes) {}

  int insert(Vec2 p) {
    const long long kr = key(p.r), kz = key(p.z);
    for (long long dr = -1; dr <= 1; ++dr)
      for (long long dz = -1; dz <= 1; ++dz) {
        auto it = buckets_.find({kr + dr, kz + dz});
        if (it == buckets_.end()) continue;
        for (int idx : it->second)
          if (norm((*nodes_)[idx].p - p) <= snap_) return idx;
      }
    MeshNode n;
    n.p = p;
    n.on_axis = (p.r == 0.0);
    const int idx = static_cast<int>(nodes_->size());
    nodes_->push_back(n);
    buckets_[{kr, kz}].push_back(idx);
    return idx;
  }

 private:
  long long key(Real x) const { return static_cast<long long>(std::floor(x / snap_)); }
  Real snap_;
  std::vector<MeshNode>* nodes_;
  std::map<std::pair<long long, long long>, std::vector<int>> buckets_;
};

/// Node positions along a straight segment, graded by the size field.
///
/// The element count follows from the mesh-density integral N = int ds / h(s);
/// the nodes are then placed at equal increments of that integral, which is the
/// exact statement of "every element is one h long".  The segment endpoints are
/// always the first and last node, so a geometric corner survives exactly.
std::vector<Vec2> subdivide(Vec2 p, Vec2 q, const SizeField& h) {
  const Real L = norm(q - p);
  if (!(L > 0.0)) throw std::runtime_error("BoundaryMesh: zero-length geometric segment");
  const Vec2 u = (q - p) / L;

  // Sub-sampling resolution: fine enough that the smallest h on the segment is
  // covered by several samples.  Purely a function of the geometry, hence
  // deterministic.
  Real h_lo = h(p);
  const int kScan = 200;
  for (int i = 0; i <= kScan; ++i) h_lo = std::min(h_lo, h(p + (L * i / kScan) * u));
  long long m = static_cast<long long>(std::ceil(8.0 * L / h_lo));
  m = std::min<long long>(std::max<long long>(m, 512), 400000);
  const int M = static_cast<int>(m);

  // Cumulative density integral by the trapezoidal rule.
  std::vector<Real> phi(M + 1, 0.0);
  Real prev = 1.0 / h(p);
  for (int i = 1; i <= M; ++i) {
    const Real cur = 1.0 / h(p + (L * i / M) * u);
    phi[i] = phi[i - 1] + 0.5 * (prev + cur) * (L / M);
    prev = cur;
  }
  const Real total = phi[M];

  int n = static_cast<int>(std::ceil(total - 1e-9));
  n = std::max(n, mesher::kMinElementsPerSegment);

  std::vector<Vec2> out;
  out.reserve(n + 1);
  out.push_back(p);
  int i = 0;
  for (int k = 1; k < n; ++k) {
    const Real target = total * k / n;
    while (i < M && phi[i + 1] < target) ++i;
    const Real denom = phi[i + 1] - phi[i];
    const Real frac = denom > 0.0 ? (target - phi[i]) / denom : 0.0;
    const Real s = L * (i + frac) / M;
    out.push_back(p + s * u);
  }
  out.push_back(q);

  for (std::size_t k = 1; k < out.size(); ++k)
    if (!(norm(out[k] - out[k - 1]) > 0.0))
      throw std::runtime_error("BoundaryMesh: subdivision produced a zero-length element");
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Generation
// ---------------------------------------------------------------------------

BoundaryMesh BoundaryMesh::generate(const DeviceGeometry& g) {
  BoundaryMesh mesh;
  mesh.size_ = SizeField::from_geometry(g);

  const DeviceParameters& p = g.parameters();
  const Real diag = domain_diagonal(p);
  const Real snap = mesher::kNodeSnapRelative * diag;
  NodeRegistry reg(snap, &mesh.nodes_);

  const auto& curves = g.boundaries();
  for (std::size_t c = 0; c < curves.size(); ++c) {
    const BoundaryCurve& curve = curves[c];
    std::vector<Vec2> pts = curve.points;
    if (pts.size() < 2)
      throw std::runtime_error("BoundaryMesh: boundary curve with fewer than two points");

    const bool axis = std::all_of(pts.begin(), pts.end(), [](Vec2 v) { return v.r == 0.0; });
    bool reversed = false;

    // --- orientation ------------------------------------------------------
    // The order of the input polyline is not a contract, so the direction is
    // decided geometrically and then verified for every element.
    if (axis) {
      // Symmetry elements: traverse towards decreasing z, so that
      // perp(tangent) = (-1, 0) points out of the meridian half-plane.
      if (pts.front().z < pts.back().z) {
        std::reverse(pts.begin(), pts.end());
        reversed = true;
      }
      if (curve.side_a != curve.side_b)
        throw std::runtime_error("BoundaryMesh: axis curve '" + curve.name +
                                 "' claims two different regions");
    } else {
      const Vec2 t0 = normalized(pts[1] - pts[0]);
      const Vec2 n0 = perp(t0);
      const Vec2 m0 = 0.5 * (pts[0] + pts[1]);
      const Real d0 = mesher::kSideProbeRelative * norm(pts[1] - pts[0]);
      const Region plus = region_at(g, m0 + d0 * n0);
      const Region minus = region_at(g, m0 - d0 * n0);
      if (plus == curve.side_b && minus == curve.side_a) {
        // keep
      } else if (plus == curve.side_a && minus == curve.side_b) {
        std::reverse(pts.begin(), pts.end());
        reversed = true;
      } else {
        throw std::runtime_error("BoundaryMesh: material assignment of boundary '" +
                                 curve.name + "' is not reproduced by the geometry");
      }
    }

    // --- discretise segment by segment ------------------------------------
    const int n_seg = static_cast<int>(pts.size()) - 1;
    for (std::size_t s = 1; s < pts.size(); ++s) {
      // `segment` always indexes the ORIGINAL curve polyline, even where the
      // traversal was reversed above, so that a consumer can map an element
      // back onto the geometry it came from.
      const int seg = reversed ? n_seg - static_cast<int>(s) : static_cast<int>(s - 1);
      const std::vector<Vec2> nodes = subdivide(pts[s - 1], pts[s], mesh.size_);
      for (std::size_t k = 1; k < nodes.size(); ++k) {
        BoundaryElement e;
        e.id = curve.id;
        e.curve = static_cast<int>(c);
        e.segment = seg;
        e.a = nodes[k - 1];
        e.b = nodes[k];
        e.node_a = reg.insert(e.a);
        e.node_b = reg.insert(e.b);
        e.a = mesh.nodes_[e.node_a].p;
        e.b = mesh.nodes_[e.node_b].p;
        e.meridian_length = norm(e.b - e.a);
        if (!(e.meridian_length > 0.0))
          throw std::runtime_error("BoundaryMesh: zero-length element on '" + curve.name + "'");
        if (e.a.r < 0.0 || e.b.r < 0.0)
          throw std::runtime_error("BoundaryMesh: negative radius on '" + curve.name + "'");
        e.tangent = (e.b - e.a) / e.meridian_length;
        e.side_a = curve.side_a;
        e.side_b = curve.side_b;

        if (axis) {
          e.kind = ElementKind::AxisSymmetry;
          e.normal = Vec2{-1.0, 0.0};
          e.revolved_area = 0.0;  // exactly, not "small"
        } else {
          e.kind = ElementKind::Ring;
          e.normal = perp(e.tangent);
          e.revolved_area = pi * (e.a.r + e.b.r) * e.meridian_length;
          if (!(e.mid_radius() > 0.0))
            throw std::runtime_error("BoundaryMesh: ring element with zero mean radius on '" +
                                     curve.name + "'");
          const Real d = mesher::kSideProbeRelative * e.meridian_length;
          const Vec2 m = e.midpoint();
          if (region_at(g, m + d * e.normal) != e.side_b ||
              region_at(g, m - d * e.normal) != e.side_a)
            throw std::runtime_error("BoundaryMesh: element on '" + curve.name +
                                     "' does not separate the regions it claims to");
        }
        mesh.elements_.push_back(e);
      }
    }
  }

  // --- annotate nodes -----------------------------------------------------
  auto mark = [&](Vec2 x, bool corner, int feature) {
    for (MeshNode& n : mesh.nodes_)
      if (norm(n.p - x) <= snap) {
        n.is_corner = n.is_corner || corner;
        if (feature >= 0) n.feature = feature;
        return true;
      }
    return false;
  };
  for (const BoundaryCurve& c : curves)
    for (const Vec2& v : c.points)
      if (!mark(v, true, -1))
        throw std::runtime_error("BoundaryMesh: geometric corner lost during meshing");
  for (std::size_t k = 0; k < g.features().size(); ++k)
    if (!mark(g.features()[k].position, true, static_cast<int>(k)))
      throw std::runtime_error("BoundaryMesh: named feature lost during meshing");

  return mesh;
}

// ---------------------------------------------------------------------------
// Queries and statistics
// ---------------------------------------------------------------------------

std::vector<const BoundaryElement*> BoundaryMesh::elements_with(BoundaryId id) const {
  std::vector<const BoundaryElement*> out;
  for (const BoundaryElement& e : elements_)
    if (e.id == id) out.push_back(&e);
  return out;
}

std::vector<const BoundaryElement*> BoundaryMesh::elements_of_curve(int curve) const {
  std::vector<const BoundaryElement*> out;
  for (const BoundaryElement& e : elements_)
    if (e.curve == curve) out.push_back(&e);
  return out;
}

namespace {

LengthStats stats_from(const std::vector<const BoundaryElement*>& es) {
  LengthStats s;
  if (es.empty()) return s;
  std::vector<Real> len;
  len.reserve(es.size());
  for (const BoundaryElement* e : es) {
    len.push_back(e->meridian_length);
    s.total_meridian_length += e->meridian_length;
    s.total_revolved_area += e->revolved_area;
  }
  std::sort(len.begin(), len.end());
  s.count = static_cast<int>(len.size());
  s.min = len.front();
  s.max = len.back();
  const std::size_t n = len.size();
  s.median = (n % 2) ? len[n / 2] : 0.5 * (len[n / 2 - 1] + len[n / 2]);
  return s;
}

}  // namespace

LengthStats BoundaryMesh::stats_of(BoundaryId id) const { return stats_from(elements_with(id)); }

LengthStats BoundaryMesh::stats_of_curve(int curve) const {
  return stats_from(elements_of_curve(curve));
}

LengthStats BoundaryMesh::stats_total() const {
  std::vector<const BoundaryElement*> all;
  for (const BoundaryElement& e : elements_) all.push_back(&e);
  return stats_from(all);
}

Real BoundaryMesh::max_neighbour_ratio() const {
  std::vector<std::vector<int>> at(nodes_.size());
  for (std::size_t i = 0; i < elements_.size(); ++i) {
    at[elements_[i].node_a].push_back(static_cast<int>(i));
    at[elements_[i].node_b].push_back(static_cast<int>(i));
  }
  Real worst = 1.0;
  for (const auto& group : at)
    for (std::size_t i = 0; i < group.size(); ++i)
      for (std::size_t j = i + 1; j < group.size(); ++j) {
        const Real a = elements_[group[i]].meridian_length;
        const Real b = elements_[group[j]].meridian_length;
        worst = std::max(worst, std::max(a / b, b / a));
      }
  return worst;
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

bool MeshReport::all_passed() const {
  for (const MeshCheck& c : checks)
    if (!c.passed) return false;
  return true;
}

int MeshReport::failures() const {
  int n = 0;
  for (const MeshCheck& c : checks) n += c.passed ? 0 : 1;
  return n;
}

void MeshReport::print(std::FILE* out) const {
  for (const MeshCheck& c : checks)
    std::fprintf(out, "  [%s] %-52s %s\n", c.passed ? "ok  " : "FAIL", c.name.c_str(),
                 c.detail.c_str());
  std::fprintf(out, "  %d von %zu Pruefungen fehlgeschlagen\n", failures(), checks.size());
}

namespace {

std::string fmt(const char* f, Real a, Real b = 0.0, Real c = 0.0) {
  char buf[256];
  std::snprintf(buf, sizeof buf, f, a, b, c);
  return buf;
}

/// Do the interiors of [a1,b1] and [a2,b2] meet?  Endpoint contact is allowed;
/// a proper crossing or a shared stretch of positive length is not.
bool overlaps(Vec2 a1, Vec2 b1, Vec2 a2, Vec2 b2, Real tol) {
  const Real L1 = norm(b1 - a1), L2 = norm(b2 - a2);
  if (L1 <= 0.0 || L2 <= 0.0) return false;
  const Vec2 u1 = (b1 - a1) / L1, u2 = (b2 - a2) / L2;

  // Signed perpendicular distances, in metres.
  const Real d1 = cross2(u1, a2 - a1), d2 = cross2(u1, b2 - a1);
  const Real d3 = cross2(u2, a1 - a2), d4 = cross2(u2, b1 - a2);

  if (std::abs(d1) <= tol && std::abs(d2) <= tol) {  // collinear
    const Real t1 = dot(u1, a2 - a1), t2 = dot(u1, b2 - a1);
    const Real lo = std::max(0.0, std::min(t1, t2));
    const Real hi = std::min(L1, std::max(t1, t2));
    return (hi - lo) > tol;
  }
  const bool cross_a = (d1 > tol && d2 < -tol) || (d1 < -tol && d2 > tol);
  const bool cross_b = (d3 > tol && d4 < -tol) || (d3 < -tol && d4 > tol);
  return cross_a && cross_b;
}

}  // namespace

MeshReport BoundaryMesh::validate(const DeviceGeometry& g) const {
  MeshReport rep;
  auto add = [&rep](const char* name, bool ok, std::string detail) {
    rep.checks.push_back({name, ok, std::move(detail)});
  };

  const Real diag = domain_diagonal(g.parameters());
  const Real snap = mesher::kNodeSnapRelative * diag;

  // --- 1. no degenerate elements -----------------------------------------
  {
    int bad_len = 0, bad_radius = 0, bad_axis = 0;
    Real shortest = std::numeric_limits<Real>::infinity();
    for (const BoundaryElement& e : elements_) {
      if (!(e.meridian_length > 0.0)) ++bad_len;
      shortest = std::min(shortest, e.meridian_length);
      if (e.is_axis()) {
        if (e.a.r != 0.0 || e.b.r != 0.0 || e.revolved_area != 0.0) ++bad_axis;
      } else {
        if (!(e.mid_radius() > 0.0) || e.a.r < 0.0 || e.b.r < 0.0) ++bad_radius;
      }
    }
    add("keine Nulllaengen-Elemente", bad_len == 0,
        fmt("kuerzestes Element %.5g m", shortest));
    add("keine Nullradien-Ringelemente", bad_radius == 0,
        fmt("%.0f Ringelemente", static_cast<Real>(elements_.size() -
                                                   elements_with(BoundaryId::SymmetryAxis).size())));
    add("Achsenelemente exakt bei r = 0 mit Flaeche 0", bad_axis == 0,
        fmt("%.0f Symmetrieelemente",
            static_cast<Real>(elements_with(BoundaryId::SymmetryAxis).size())));
  }

  // --- 2. no duplicate elements ------------------------------------------
  {
    std::set<std::pair<int, int>> pairs;
    int dup = 0;
    for (const BoundaryElement& e : elements_) {
      auto key = std::minmax(e.node_a, e.node_b);
      if (!pairs.insert({key.first, key.second}).second) ++dup;
    }
    std::vector<std::pair<Real, Real>> mids;
    for (const BoundaryElement& e : elements_) mids.push_back({e.midpoint().r, e.midpoint().z});
    std::sort(mids.begin(), mids.end());
    int dup_mid = 0;
    for (std::size_t i = 1; i < mids.size(); ++i)
      if (std::abs(mids[i].first - mids[i - 1].first) <= snap &&
          std::abs(mids[i].second - mids[i - 1].second) <= snap)
        ++dup_mid;
    add("keine doppelten Elemente", dup == 0 && dup_mid == 0,
        fmt("%.0f Elemente, %.0f Knoten", static_cast<Real>(elements_.size()),
            static_cast<Real>(nodes_.size())));
  }

  // --- 3. no gaps: connectivity ------------------------------------------
  {
    std::vector<int> degree(nodes_.size(), 0);
    for (const BoundaryElement& e : elements_) {
      ++degree[e.node_a];
      ++degree[e.node_b];
    }
    int lonely = 0, min_deg = 1 << 30;
    for (int d : degree) {
      if (d < 2) ++lonely;
      min_deg = std::min(min_deg, d);
    }
    add("kein Knoten mit weniger als zwei Elementen (keine Luecken)", lonely == 0,
        fmt("kleinster Knotengrad %.0f", static_cast<Real>(min_deg)));

    // Each curve is one unbroken chain from its first to its last vertex.
    int broken = 0, endpoints_lost = 0;
    for (std::size_t c = 0; c < g.boundaries().size(); ++c) {
      const auto es = elements_of_curve(static_cast<int>(c));
      if (es.empty()) {
        ++broken;
        continue;
      }
      for (std::size_t i = 1; i < es.size(); ++i)
        if (es[i]->node_a != es[i - 1]->node_b) ++broken;
      const auto& pts = g.boundaries()[c].points;
      const Vec2 first = es.front()->a, last = es.back()->b;
      const bool forward = norm(first - pts.front()) <= snap && norm(last - pts.back()) <= snap;
      const bool reverse = norm(first - pts.back()) <= snap && norm(last - pts.front()) <= snap;
      if (!forward && !reverse) ++endpoints_lost;
    }
    add("jede Randkurve ist eine ununterbrochene Kette", broken == 0,
        fmt("%.0f Randkurven", static_cast<Real>(g.boundaries().size())));
    add("Kettenenden liegen auf den Kurvenenden", endpoints_lost == 0, "");
  }

  // --- 4. corners and features survive exactly ---------------------------
  {
    int lost_corner = 0, lost_feature = 0;
    auto node_at = [&](Vec2 x) {
      for (const MeshNode& n : nodes_)
        if (norm(n.p - x) <= snap) return true;
      return false;
    };
    for (const BoundaryCurve& c : g.boundaries())
      for (const Vec2& v : c.points)
        if (!node_at(v)) ++lost_corner;
    for (const NamedFeature& f : g.features())
      if (!node_at(f.position)) ++lost_feature;
    add("alle Geometrieecken sind exakt Netzknoten", lost_corner == 0, "");
    add("alle benannten Merkmale sind exakt Netzknoten", lost_feature == 0,
        fmt("%.0f Merkmale", static_cast<Real>(g.features().size())));

    // No element may straddle a corner: its two endpoints must lie on one
    // straight segment.  Guaranteed by construction; checked because merging
    // across a corner would silently smear a material boundary.
    int straddling = 0;
    for (const BoundaryElement& e : elements_) {
      const auto& pts = g.boundaries()[e.curve].points;
      bool on_one = false;
      for (std::size_t i = 1; i < pts.size(); ++i) {
        const Real da = dist_point_segment(e.a, pts[i - 1], pts[i]);
        const Real db = dist_point_segment(e.b, pts[i - 1], pts[i]);
        if (da <= 1e-9 * diag && db <= 1e-9 * diag) on_one = true;
      }
      if (!on_one) ++straddling;
    }
    add("kein Element ueberspannt eine Geometrieecke", straddling == 0, "");
  }

  // --- 5. no unintended overlaps -----------------------------------------
  {
    int hits = 0;
    const Real tol = 1e-9 * diag;
    for (std::size_t i = 0; i < elements_.size(); ++i)
      for (std::size_t j = i + 1; j < elements_.size(); ++j) {
        const BoundaryElement& u = elements_[i];
        const BoundaryElement& v = elements_[j];
        if (std::max(u.a.r, u.b.r) < std::min(v.a.r, v.b.r) - tol) continue;
        if (std::min(u.a.r, u.b.r) > std::max(v.a.r, v.b.r) + tol) continue;
        if (std::max(u.a.z, u.b.z) < std::min(v.a.z, v.b.z) - tol) continue;
        if (std::min(u.a.z, u.b.z) > std::max(v.a.z, v.b.z) + tol) continue;
        if (overlaps(u.a, u.b, v.a, v.b, tol)) ++hits;
      }
    add("keine ungewollten Ueberlappungen oder Kreuzungen", hits == 0,
        fmt("%.0f Elementpaare geprueft",
            static_cast<Real>(elements_.size() * (elements_.size() - 1) / 2)));
  }

  // --- 6. material assignment is unambiguous -----------------------------
  {
    int wrong = 0, ambiguous = 0;
    for (const BoundaryElement& e : elements_) {
      if (e.is_axis()) {
        if (e.side_a != e.side_b) ++ambiguous;
        continue;
      }
      if (e.side_a == e.side_b) ++ambiguous;
      const Real d = mesher::kSideProbeRelative * e.meridian_length;
      const Vec2 m = e.midpoint();
      if (region_at(g, m + d * e.normal) != e.side_b) ++wrong;
      if (region_at(g, m - d * e.normal) != e.side_a) ++wrong;
    }
    add("jedes Element trennt genau die angegebenen Gebiete", wrong == 0 && ambiguous == 0,
        "Normale zeigt von side_a nach side_b");

    // Every element carries exactly one boundary identifier, and every curve
    // name is unique -- otherwise an element could not be addressed.
    std::vector<std::string> names;
    for (const BoundaryCurve& c : g.boundaries()) names.push_back(c.name);
    std::sort(names.begin(), names.end());
    add("Randkennungen sind eindeutig zugeordnet",
        std::adjacent_find(names.begin(), names.end()) == names.end(),
        fmt("%.0f benannte Randkurven", static_cast<Real>(names.size())));
  }

  // --- 7. closed, correctly oriented region boundaries -------------------
  {
    bool all_ok = true;
    std::string detail;
    for (Region r : {Region::Liquid, Region::EmitterSolid, Region::ExtractorSolid,
                     Region::Vacuum}) {
      // Orient every adjacent element so the region lies on the left, i.e. so
      // that perp(tangent) points away from it.  For an axis element the region
      // is the r > 0 side, and the stored orientation already does that.
      std::vector<int> out_deg(nodes_.size(), 0), in_deg(nodes_.size(), 0);
      Real area = 0.0, vol = 0.0;
      int n_el = 0;
      for (const BoundaryElement& e : elements_) {
        bool flip = false;
        if (e.is_axis()) {
          if (e.side_a != r) continue;  // stored direction already keeps r on the left
        } else if (e.side_a == r) {
          flip = false;  // normal points away from r, so r is already on the left
        } else if (e.side_b == r) {
          flip = true;
        } else {
          continue;
        }
        const Vec2 a = flip ? e.b : e.a;
        const Vec2 b = flip ? e.a : e.b;
        const int ia = flip ? e.node_b : e.node_a;
        const int ib = flip ? e.node_a : e.node_b;
        ++n_el;
        ++out_deg[ia];
        ++in_deg[ib];
        area += 0.5 * (a.r * b.z - b.r * a.z);
        vol += pi * (a.r * a.r + a.r * b.r + b.r * b.r) / 3.0 * (b.z - a.z);
      }
      // Closed and manifold: every node used by this region is entered exactly
      // as often as it is left, and at most once each way.
      bool closed = true;
      for (std::size_t k = 0; k < nodes_.size(); ++k)
        if (out_deg[k] != in_deg[k] || out_deg[k] > 1) closed = false;

      const RegionBody& body = g.region(r);
      const Real want_a = body.meridian_area();
      const Real want_v = body.revolved_volume();
      const bool ok = closed && area > 0.0 &&
                      std::abs(area - want_a) <= 1e-12 * want_a &&
                      std::abs(vol - want_v) <= 1e-12 * want_v;
      all_ok = all_ok && ok;
      char buf[220];
      std::snprintf(buf, sizeof buf, "%s: %d Elemente, dA/A=%.2g, dV/V=%.2g%s; ", to_string(r),
                    n_el, std::abs(area - want_a) / want_a, std::abs(vol - want_v) / want_v,
                    closed ? "" : ", NICHT GESCHLOSSEN");
      detail += buf;
    }
    add("Gebietsraender geschlossen und korrekt orientiert (CCW)", all_ok, detail);
  }

  // --- 8. analytic surfaces of revolution preserved ----------------------
  {
    Real worst_area = 0.0, worst_len = 0.0;
    for (std::size_t c = 0; c < g.boundaries().size(); ++c) {
      const LengthStats s = stats_of_curve(static_cast<int>(c));
      const BoundaryCurve& curve = g.boundaries()[c];
      const Real la = curve.revolved_area(), ll = curve.meridian_length();
      if (la > 0.0) worst_area = std::max(worst_area, std::abs(s.total_revolved_area - la) / la);
      if (ll > 0.0) worst_len = std::max(worst_len, std::abs(s.total_meridian_length - ll) / ll);
    }
    add("Rotationsflaechen erhalten (Summe = analytisch)", worst_area <= 1e-12,
        fmt("groesster relativer Fehler %.3g", worst_area));
    add("Meridianlaengen erhalten", worst_len <= 1e-12,
        fmt("groesster relativer Fehler %.3g", worst_len));
  }

  // --- 9. bounded size ratio ---------------------------------------------
  {
    const Real ratio = max_neighbour_ratio();
    add("kontrolliertes Elementgroessenverhaeltnis", ratio <= mesher::kMaxNeighbourRatio,
        fmt("groesstes Verhaeltnis benachbarter Elemente %.4g (Grenze %.2f)", ratio,
            mesher::kMaxNeighbourRatio));
  }

  // --- 10. reproducibility -----------------------------------------------
  {
    const BoundaryMesh again = BoundaryMesh::generate(g);
    bool same = again.nodes_.size() == nodes_.size() &&
                again.elements_.size() == elements_.size();
    if (same)
      for (std::size_t i = 0; i < nodes_.size(); ++i)
        same = same && again.nodes_[i].p.r == nodes_[i].p.r &&
               again.nodes_[i].p.z == nodes_[i].p.z;
    if (same)
      for (std::size_t i = 0; i < elements_.size(); ++i)
        same = same && again.elements_[i].node_a == elements_[i].node_a &&
               again.elements_[i].node_b == elements_[i].node_b;
    add("reproduzierbare Knotenanordnung (bitweise)", same,
        "zweiter Aufbau aus denselben Parametern");
  }

  return rep;
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

void BoundaryMesh::print(std::FILE* out, const DeviceGeometry& g) const {
  std::fprintf(out, "Automatischer achsensymmetrischer Randvernetzer\n");
  std::fprintf(out, "  Knoten   : %zu\n", nodes_.size());
  std::fprintf(out, "  Elemente : %zu\n\n", elements_.size());
  size_.print(out);

  std::fprintf(out, "\nElementstatistik je Randkurve (Laengen in m)\n");
  std::fprintf(out,
               "  %-34s %-24s %5s %11s %11s %11s %12s %12s\n", "Randkurve", "Kennung", "n",
               "min", "median", "max", "Sum_Laenge", "Sum_Flaeche");
  for (std::size_t c = 0; c < g.boundaries().size(); ++c) {
    const LengthStats s = stats_of_curve(static_cast<int>(c));
    std::fprintf(out, "  %-34s %-24s %5d %11.4g %11.4g %11.4g %12.5g %12.5g\n",
                 g.boundaries()[c].name.c_str(), to_string(g.boundaries()[c].id), s.count, s.min,
                 s.median, s.max, s.total_meridian_length, s.total_revolved_area);
  }

  std::fprintf(out, "\nElementstatistik je Randkennung\n");
  for (BoundaryId id : {BoundaryId::SymmetryAxis, BoundaryId::EmitterOuterSurface,
                        BoundaryId::EmitterTipLand, BoundaryId::BoreWall,
                        BoundaryId::FreeSurfaceReference, BoundaryId::LiquidInlet,
                        BoundaryId::ExtractorSurface, BoundaryId::OpenBoundary}) {
    const LengthStats s = stats_of(id);
    std::fprintf(out, "  %-24s n=%-5d min=%10.4g median=%10.4g max=%10.4g\n", to_string(id),
                 s.count, s.min, s.median, s.max);
    std::fprintf(out, "  %-24s %s\n", "", boundary_long_name(id));
  }
  const LengthStats t = stats_total();
  std::fprintf(out, "\n  gesamt: n=%d, min=%.4g m, median=%.4g m, max=%.4g m\n", t.count, t.min,
               t.median, t.max);
  std::fprintf(out, "  groesstes Verhaeltnis benachbarter Elemente: %.4g\n",
               max_neighbour_ratio());
}

void BoundaryMesh::write_csv(const std::string& dir, const DeviceGeometry& g) const {
  const std::string d = dir + "/";
  {
    std::FILE* f = std::fopen((d + "mesh_nodes.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open mesh_nodes.csv");
    std::fprintf(f, "# boundary mesh nodes in the meridian half-plane (r, z), SI\n");
    std::fprintf(f, "index,r_m,z_m,on_axis,is_corner,feature\n");
    for (std::size_t i = 0; i < nodes_.size(); ++i)
      std::fprintf(f, "%zu,%.17e,%.17e,%d,%d,%s\n", i, nodes_[i].p.r, nodes_[i].p.z,
                   nodes_[i].on_axis ? 1 : 0, nodes_[i].is_corner ? 1 : 0,
                   nodes_[i].feature >= 0 ? to_string(g.features()[nodes_[i].feature].id) : "");
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "mesh_elements.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open mesh_elements.csv");
    std::fprintf(f, "# boundary elements; revolved_area = pi*(r_a+r_b)*L, exactly 0 on the axis\n");
    std::fprintf(f, "index,curve_name,boundary_id,kind,side_a,side_b,node_a,node_b,"
                    "r_a_m,z_a_m,r_b_m,z_b_m,r_mid_m,z_mid_m,n_r,n_z,"
                    "meridian_length_m,revolved_area_m2\n");
    for (std::size_t i = 0; i < elements_.size(); ++i) {
      const BoundaryElement& e = elements_[i];
      std::fprintf(f,
                   "%zu,%s,%s,%s,%s,%s,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,"
                   "%.9e,%.9e,%.17e,%.17e\n",
                   i, g.boundaries()[e.curve].name.c_str(), to_string(e.id), to_string(e.kind),
                   to_string(e.side_a), to_string(e.side_b), e.node_a, e.node_b, e.a.r, e.a.z,
                   e.b.r, e.b.z, e.midpoint().r, e.midpoint().z, e.normal.r, e.normal.z,
                   e.meridian_length, e.revolved_area);
    }
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "mesh_boundaries.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open mesh_boundaries.csv");
    std::fprintf(f, "# per-curve and per-identifier element statistics, SI units\n");
    std::fprintf(f, "scope,name,boundary_id,long_name,n_elements,len_min_m,len_median_m,"
                    "len_max_m,total_meridian_length_m,total_revolved_area_m2\n");
    for (std::size_t c = 0; c < g.boundaries().size(); ++c) {
      const LengthStats s = stats_of_curve(static_cast<int>(c));
      std::fprintf(f, "curve,%s,%s,\"%s\",%d,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                   g.boundaries()[c].name.c_str(), to_string(g.boundaries()[c].id),
                   boundary_long_name(g.boundaries()[c].id), s.count, s.min, s.median, s.max,
                   s.total_meridian_length, s.total_revolved_area);
    }
    for (BoundaryId id : {BoundaryId::SymmetryAxis, BoundaryId::EmitterOuterSurface,
                          BoundaryId::EmitterTipLand, BoundaryId::BoreWall,
                          BoundaryId::FreeSurfaceReference, BoundaryId::LiquidInlet,
                          BoundaryId::ExtractorSurface, BoundaryId::OpenBoundary}) {
      const LengthStats s = stats_of(id);
      std::fprintf(f, "boundary_id,%s,%s,\"%s\",%d,%.9e,%.9e,%.9e,%.9e,%.9e\n", to_string(id),
                   to_string(id), boundary_long_name(id), s.count, s.min, s.median, s.max,
                   s.total_meridian_length, s.total_revolved_area);
    }
    const LengthStats t = stats_total();
    std::fprintf(f, "total,all,all,\"gesamtes Randnetz\",%d,%.9e,%.9e,%.9e,%.9e,%.9e\n", t.count,
                 t.min, t.median, t.max, t.total_meridian_length, t.total_revolved_area);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "mesh_size_field.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open mesh_size_field.csv");
    std::fprintf(f, "# refinement sources of the size function h(x) = min_s (h_s + G|x-x_s|)\n");
    std::fprintf(f, "# G=%.4g, h_max=%.9e m\n", mesher::kGradation, size_.h_max());
    std::fprintf(f, "r_m,z_m,local_feature_size_m,h_m,named_feature,origin\n");
    for (const SizeField::Source& s : size_.sources())
      std::fprintf(f, "%.17e,%.17e,%.9e,%.9e,%d,\"%s\"\n", s.x.r, s.x.z, s.local_feature_size,
                   s.h, s.is_named_feature ? 1 : 0, s.origin.c_str());
    std::fclose(f);
  }
}

}  // namespace es
