#include "es/geometry.hpp"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

const char* tag_name(Tag t) {
  switch (t) {
    case Tag::Emitter: return "emitter";
    case Tag::FreeSurface: return "meniscus";
    case Tag::Extractor: return "extractor";
    case Tag::Collector: return "collector";
    default: return "other";
  }
}

// ---------------------------------------------------------------------------
// Grading
// ---------------------------------------------------------------------------
// Element size varies linearly along the run, h(s) = h0 + (h1-h0) s/L.  The
// number of elements is the integral of ds/h(s); nodes are placed so that each
// element carries the same amount of that integral.  This keeps the transition
// smooth instead of producing a jump between a fine and a coarse block.
std::vector<Real> graded_parameters(Real length, Real h0, Real h1) {
  if (!(length > 0.0)) return {0.0, 1.0};
  h0 = std::max(h0, 1e-14 * length);
  h1 = std::max(h1, 1e-14 * length);

  const Real dh = h1 - h0;
  Real ntot;
  if (std::abs(dh) < 1e-12 * h0) {
    ntot = length / h0;
  } else {
    ntot = length * std::log(h1 / h0) / dh;
  }
  const Index n = std::max<Index>(1, static_cast<Index>(std::ceil(ntot - 1e-9)));

  std::vector<Real> t(static_cast<std::size_t>(n) + 1);
  for (Index i = 0; i <= n; ++i) {
    const Real u = static_cast<Real>(i) * ntot / static_cast<Real>(n);
    Real s;
    if (std::abs(dh) < 1e-12 * h0) {
      s = u * h0;
    } else {
      s = (length * h0 / dh) * (std::exp(u * dh / length) - 1.0);
    }
    t[static_cast<std::size_t>(i)] = std::clamp(s / length, 0.0, 1.0);
  }
  t.front() = 0.0;
  t.back() = 1.0;
  return t;
}

// ---------------------------------------------------------------------------
// Mesh assembly
// ---------------------------------------------------------------------------

void Mesh::begin_body(Tag tag, Real potential) {
  nodes_pending_.clear();
  tag_pending_ = tag;
  pot_pending_ = potential;
}

void Mesh::add_node(Vec2 p) { nodes_pending_.push_back(p); }

void Mesh::line_to(Vec2 to, Real h_start, Real h_end) {
  if (nodes_pending_.empty()) { nodes_pending_.push_back(to); return; }
  const Vec2 from = nodes_pending_.back();
  const Vec2 d = to - from;
  const Real L = norm(d);
  if (L <= 0.0) return;
  const std::vector<Real> t = graded_parameters(L, h_start, h_end);
  for (std::size_t i = 1; i < t.size(); ++i) nodes_pending_.push_back(from + t[i] * d);
}

void Mesh::arc_to(Vec2 center, Real radius, Real a_start, Real a_end, Real h) {
  const Vec2 p0{center.r + radius * std::cos(a_start), center.z + radius * std::sin(a_start)};
  if (nodes_pending_.empty()) {
    nodes_pending_.push_back(p0);
  } else if (norm(nodes_pending_.back() - p0) > 1e-9 * std::max(radius, 1e-30)) {
    // Small gap between the previous run and the arc start: snap a node in
    // rather than silently skewing the first arc element.
    nodes_pending_.push_back(p0);
  }
  const Real sweep = a_end - a_start;
  const Real arclen = std::abs(sweep) * radius;
  const Index n = std::max<Index>(1, static_cast<Index>(std::ceil(arclen / std::max(h, 1e-30))));
  for (Index i = 1; i <= n; ++i) {
    const Real ang = a_start + sweep * static_cast<Real>(i) / static_cast<Real>(n);
    nodes_pending_.push_back({center.r + radius * std::cos(ang), center.z + radius * std::sin(ang)});
  }
}

void Mesh::end_body(bool close_loop) {
  if (nodes_pending_.size() < 2) { nodes_pending_.clear(); return; }
  if (close_loop) {
    Real scale = 0.0;
    for (const Vec2& p : nodes_pending_) scale = std::max(scale, std::abs(p.r) + std::abs(p.z));
    if (norm(nodes_pending_.front() - nodes_pending_.back()) > 1e-10 * std::max(scale, 1e-30))
      nodes_pending_.push_back(nodes_pending_.front());
  }
  const int body = body_counter_++;
  for (std::size_t i = 0; i + 1 < nodes_pending_.size(); ++i) {
    Element e;
    e.a = nodes_pending_[i];
    e.b = nodes_pending_[i + 1];
    e.tag = tag_pending_;
    e.potential = pot_pending_;
    e.body = body;
    elems.push_back(e);
  }
  nodes_pending_.clear();
  finalize();
}

void Mesh::finalize() {
  for (Element& e : elems) {
    const Vec2 d = e.b - e.a;
    e.len = norm(d);
    e.tangent = e.len > 0 ? d / e.len : Vec2{1, 0};
    e.normal = perp(e.tangent);
    e.mid = 0.5 * (e.a + e.b);
    e.area = 2.0 * constants::pi * std::max(e.mid.r, 0.0) * e.len;
  }
}

Real Mesh::total_area() const {
  Real s = 0.0;
  for (const Element& e : elems) s += e.area;
  return s;
}

void Mesh::write_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "i,r_a,z_a,r_b,z_b,r_mid,z_mid,len,area,nr,nz,tag,body,potential\n");
  for (std::size_t i = 0; i < elems.size(); ++i) {
    const Element& e = elems[i];
    std::fprintf(f, "%zu,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.6f,%.6f,%s,%d,%.9e\n", i,
                 e.a.r, e.a.z, e.b.r, e.b.z, e.mid.r, e.mid.z, e.len, e.area,
                 e.normal.r, e.normal.z, tag_name(e.tag), e.body, e.potential);
  }
  std::fclose(f);
}

// ---------------------------------------------------------------------------
// Verification geometries
// ---------------------------------------------------------------------------

Mesh make_sphere(Real R, Real V, int n_elem) {
  Mesh m;
  m.begin_body(Tag::Emitter, V);
  const int n = std::max(4, n_elem);
  for (int i = 0; i <= n; ++i) {
    const Real ang = -0.5 * constants::pi + constants::pi * static_cast<Real>(i) / n;
    m.add_node({R * std::cos(ang), R * std::sin(ang)});
  }
  m.end_body(false);  // closed through the axis
  return m;
}

Mesh make_prolate_spheroid(Real a, Real b, Real V, int n_elem) {
  Mesh m;
  m.begin_body(Tag::Emitter, V);
  const int n = std::max(8, n_elem);
  for (int i = 0; i <= n; ++i) {
    // Cluster nodes toward the poles, where the curvature is highest.
    const Real x = static_cast<Real>(i) / n;
    const Real f = 0.5 * (1.0 - std::cos(constants::pi * x));
    const Real th = -0.5 * constants::pi + constants::pi * f;
    m.add_node({b * std::cos(th), a * std::sin(th)});
  }
  m.end_body(false);
  return m;
}

Real spheroid_capacitance(Real a, Real b) {
  if (a <= b * (1.0 + 1e-9)) return 4.0 * constants::pi * constants::eps0 * a;
  const Real c = std::sqrt(a * a - b * b);
  return 8.0 * constants::pi * constants::eps0 * c / std::log((a + c) / (a - c));
}

Real spheroid_tip_field(Real a, Real b, Real V) {
  if (a <= b * (1.0 + 1e-9)) return V / a;
  const Real c = std::sqrt(a * a - b * b);
  return 2.0 * V * c / (b * b * std::log((a + c) / (a - c)));
}

// ---------------------------------------------------------------------------
// Emitter geometries
// ---------------------------------------------------------------------------

Mesh make_capillary(const CapillaryParams& p) {
  if (!(p.r_outer > p.r_inner && p.r_inner > 0.0))
    throw std::runtime_error("make_capillary: need 0 < r_inner < r_outer");

  const Real rr = (p.rim_radius > 0.0) ? std::min(p.rim_radius, 0.5 * (p.r_outer - p.r_inner))
                                       : 0.5 * (p.r_outer - p.r_inner);
  const Real h_tip = (p.h_tip > 0.0) ? p.h_tip : p.r_inner / 12.0;
  const Real h_far = (p.h_far > 0.0) ? p.h_far : p.shank_length / 25.0;
  const Real z_top = p.z_tip;
  const Real z_bot = p.z_tip - p.shank_length;

  Mesh m;
  m.begin_body(Tag::Emitter, p.potential);
  m.add_node({p.r_inner, z_top - rr});
  m.line_to({p.r_inner, z_bot}, h_tip, h_far);                 // bore wall, downward
  m.line_to({p.r_outer, z_bot}, h_far);                        // bottom closure
  m.line_to({p.r_outer, z_top - rr}, h_far, h_tip);            // outer wall, upward
  m.arc_to({p.r_outer - rr, z_top - rr}, rr, 0.0, 0.5 * constants::pi, h_tip);
  if (p.r_inner + rr < p.r_outer - rr)
    m.line_to({p.r_inner + rr, z_top}, h_tip);                 // flat annular rim face
  m.arc_to({p.r_inner + rr, z_top - rr}, rr, 0.5 * constants::pi, constants::pi, h_tip);
  m.end_body(true);
  return m;
}

Mesh make_capillary_open(const OpenCapillaryParams& p) {
  if (!(p.r_outer > p.r_bore && p.r_bore > 0.0))
    throw std::runtime_error("make_capillary_open: need 0 < r_bore < r_outer");
  const Real h_rim = (p.h_rim > 0.0) ? p.h_rim : p.r_bore / 12.0;
  const Real h_far = (p.h_far > 0.0) ? p.h_far : p.shank_length / 25.0;
  const Real z_bot = p.z_rim - p.shank_length;

  Mesh m;
  m.begin_body(Tag::Emitter, p.potential);
  m.add_node({0.0, z_bot});
  m.line_to({p.r_outer, z_bot}, h_far);            // fictitious bottom cap
  m.line_to({p.r_outer, p.z_rim}, h_far, h_rim);   // outer wall
  m.line_to({p.r_bore, p.z_rim}, h_rim);           // annular rim face
  m.end_body(false);                               // open: the meniscus continues
  return m;
}

Mesh make_needle(const NeedleParams& p) {
  if (!(p.tip_radius > 0.0)) throw std::runtime_error("make_needle: tip_radius must be > 0");
  if (!(p.half_angle > 1e-4 && p.half_angle < 0.5 * constants::pi - 1e-4))
    throw std::runtime_error("make_needle: half_angle out of range");
  if (!(p.shank_radius > p.tip_radius))
    throw std::runtime_error("make_needle: shank_radius must exceed tip_radius");

  const Real Rt = p.tip_radius;
  const Real th = p.half_angle;
  const Real h_tip = (p.h_tip > 0.0) ? p.h_tip : Rt / 8.0;
  const Real h_far = (p.h_far > 0.0) ? p.h_far : p.length / 40.0;

  const Vec2 cen{0.0, p.z_tip - Rt};
  // The cone is tangent to the apex sphere where their outward normals agree,
  // i.e. at polar angle th about the sphere centre, measured from the +r axis.
  const Vec2 ptan{Rt * std::cos(th), cen.z + Rt * std::sin(th)};
  // Walk down the cone generatrix until it reaches the shank radius.
  const Real s = (p.shank_radius - ptan.r) / std::sin(th);
  const Real z_cone_base = ptan.z - s * std::cos(th);
  const Real z_bot = p.z_tip - p.length;
  if (z_cone_base <= z_bot)
    throw std::runtime_error("make_needle: length too short for the requested cone");

  Mesh m;
  m.begin_body(Tag::Emitter, p.potential);
  m.add_node({0.0, z_bot});
  m.line_to({p.shank_radius, z_bot}, h_far);                    // bottom disc
  m.line_to({p.shank_radius, z_cone_base}, h_far);              // cylindrical shank
  m.line_to(ptan, h_far, h_tip);                                // cone flank
  m.arc_to(cen, Rt, th, 0.5 * constants::pi, h_tip);            // apex cap
  m.end_body(false);  // both ends sit on the axis
  return m;
}

Mesh make_extractor(const ExtractorParams& p) {
  if (!(p.outer_radius > p.aperture_radius && p.aperture_radius > 0.0))
    throw std::runtime_error("make_extractor: need 0 < aperture_radius < outer_radius");

  const Real er = (p.edge_radius > 0.0) ? std::min(p.edge_radius, 0.5 * p.thickness)
                                        : 0.5 * p.thickness;
  const Real h_e = (p.h_edge > 0.0) ? p.h_edge : p.aperture_radius / 20.0;
  const Real h_f = (p.h_far > 0.0) ? p.h_far : p.outer_radius / 20.0;
  const Real z1 = p.z_plate;                 // emitter-facing face
  const Real z2 = p.z_plate + p.thickness;   // downstream face
  const Real ra = p.aperture_radius;

  Mesh m;
  m.begin_body(Tag::Extractor, p.potential);
  m.add_node({ra, z2 - er});
  if (z2 - er > z1 + er) m.line_to({ra, z1 + er}, h_e);                      // aperture wall
  m.arc_to({ra + er, z1 + er}, er, constants::pi, 1.5 * constants::pi, h_e); // lower lip
  m.line_to({p.outer_radius, z1}, h_e, h_f);                                 // upstream face
  m.line_to({p.outer_radius, z2}, h_f);                                      // outer rim
  m.line_to({ra + er, z2}, h_f, h_e);                                        // downstream face
  m.arc_to({ra + er, z2 - er}, er, 0.5 * constants::pi, constants::pi, h_e); // upper lip
  m.end_body(true);
  return m;
}

Mesh merge(const std::vector<Mesh>& parts) {
  Mesh out;
  int body_offset = 0;
  for (const Mesh& p : parts) {
    int max_body = -1;
    for (const Element& e : p.elems) {
      Element c = e;
      c.body += body_offset;
      out.elems.push_back(c);
      max_body = std::max(max_body, e.body);
    }
    body_offset += max_body + 1;
  }
  out.finalize();
  return out;
}

}  // namespace es
