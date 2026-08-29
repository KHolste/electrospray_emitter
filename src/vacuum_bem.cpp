#include "es/vacuum_bem.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {
namespace {

/// Exactly one side must be vacuum for an element to be a vacuum-facing
/// conductor panel.  Region tags decide this, nothing else.
bool separates_vacuum(const BoundaryElement& e, Region* material) {
  const bool a_vac = e.side_a == Region::Vacuum;
  const bool b_vac = e.side_b == Region::Vacuum;
  if (a_vac == b_vac) return false;  // vacuum|vacuum (axis) or neither
  *material = a_vac ? e.side_b : e.side_a;
  return true;
}

/// Conductors of the P2a model, and which electrode tag they carry.  The liquid
/// column is a perfect conductor here: the charge relaxation time of an ionic
/// liquid is orders of magnitude shorter than any other time scale in the
/// problem, so its surface is an equipotential with the metal it wets.  Finite
/// conductivity is a later phase.
///
/// The tag follows the BOUNDARY, not only the material.  Tag::FreeSurface is
/// reserved for the one surface it names -- the flat liquid surface at z = 0.
/// The liquid half of the numerical back closure is liquid too, but calling it a
/// free surface would put a numerical cap into every number reported for the
/// free surface, so it is tagged as emitter metal is: both are simply the
/// emitter electrode.
bool is_conductor(Region r, BoundaryId id, Tag* tag) {
  switch (r) {
    case Region::EmitterSolid:
      *tag = Tag::Emitter;
      return true;
    case Region::Liquid:
      *tag = (id == BoundaryId::FreeSurfaceReference) ? Tag::FreeSurface : Tag::Emitter;
      return true;
    case Region::ExtractorSolid:
      *tag = Tag::Extractor;
      return true;
    default:
      return false;
  }
}

/// The boundary identifier must agree with the material pairing.  If it does
/// not, the geometry grew an interface this phase has no rule for, and guessing
/// a potential for it would be an invented boundary condition.
void require_consistent(BoundaryId id, Region material, const std::string& curve) {
  const bool ok =
      (material == Region::Liquid &&
       (id == BoundaryId::FreeSurfaceReference ||
        id == BoundaryId::NumericalEmitterBackClosure)) ||
      (material == Region::EmitterSolid &&
       (id == BoundaryId::EmitterTipLand || id == BoundaryId::EmitterOuterSurface ||
        id == BoundaryId::NumericalEmitterBackClosure)) ||
      (material == Region::ExtractorSolid && id == BoundaryId::ExtractorSurface);
  if (!ok)
    throw std::runtime_error(
        "vacuum_bem_mesh: boundary '" + curve + "' (" + to_string(id) + ") faces vacuum against " +
        to_string(material) +
        ", which this phase has no potential for.  Refusing to invent a boundary condition.");
}

/// What the electrode is called in P2a output.  es::tag_name() calls
/// Tag::FreeSurface a "meniscus", which is exactly what the plane at z = 0 is
/// not in this phase, so it is never used here.
const char* p2a_label(Tag t) {
  switch (t) {
    case Tag::Emitter: return "emitter";
    case Tag::FreeSurface: return "flat_liquid_surface_reference";
    case Tag::Extractor: return "extractor";
    default: return "-";
  }
}

const char* reject_reason(const BoundaryElement& e) {
  if (e.is_axis()) return "Symmetrieachse: kein Interface, keine Rotationsflaeche";
  if (e.side_a == Region::Outside || e.side_b == Region::Outside)
    return "offener Domaenenrand: kein Leiter, keine Randbedingung";
  if (e.side_a == e.side_b) return "kein Materialwechsel";
  return "innere Materialgrenze, kein Vakuumkontakt";
}

}  // namespace

const char* electrode_label(Tag t) { return p2a_label(t); }

// ---------------------------------------------------------------------------

void VacuumSelectionReport::print(std::FILE* out) const {
  std::fprintf(out, "Randauswahl fuer die Vakuum-BEM (nach Rand- und Gebietskennung)\n");
  std::fprintf(out, "  %-32s %-24s %-28s %6s  %s\n", "Randkurve", "Kennung", "Gebiete", "n",
               "Entscheidung");
  for (const CurveDecision& c : curves) {
    char sides[64];
    std::snprintf(sides, sizeof sides, "%s | %s", to_string(c.side_a), to_string(c.side_b));
    std::fprintf(out, "  %-32s %-24s %-28s %6d  %s%s%s\n", c.curve.c_str(),
                 to_string(c.boundary), sides, c.n_elements,
                 c.accepted ? "JA  -> " : "nein  ", c.accepted ? p2a_label(c.tag) : "",
                 c.accepted ? "" : c.reason.c_str());
  }
  std::fprintf(out,
               "\n  uebernommen: %d von %d Randelementen"
               "   (Emitter %d, ebene Fluessigkeitsoberflaeche %d, Extraktor %d)\n",
               n_selected, n_mesh_elements, n_emitter, n_free_surface, n_extractor);
  std::fprintf(out, "  Rotationsflaechen: Emitter %.6g m^2, Fluessigkeitsoberflaeche %.6g m^2, "
                    "Extraktor %.6g m^2\n",
               revolved_area_emitter, revolved_area_free_surface, revolved_area_extractor);
  std::fprintf(out,
               "  davon NUMERISCHE Rueckschliessung: %d Panels, %.6g m^2 -- mitgeloest, damit\n"
               "  der Leiter geschlossen ist; kein Bauteil und keine Quelle eines Ergebnisses\n",
               n_numerical_closure, revolved_area_numerical_closure);
  int n_eval = 0;
  for (const VacuumPanel& pn : panels) n_eval += pn.evaluable ? 1 : 0;
  std::fprintf(out, "  physikalisch auswertbare Panels: %d von %d\n", n_eval,
               n_selected);
}

void VacuumSelectionReport::write_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# selection of vacuum-facing conductor panels; decided from region and\n");
  std::fprintf(f, "# boundary tags only, never from coordinates or element order\n");
  std::fprintf(f, "# 'numerical' marks the rearward closure: solved, never reported\n");
  std::fprintf(f, "curve,boundary_id,side_a,side_b,n_elements,accepted,bem_tag,numerical,"
                  "reason\n");
  for (const CurveDecision& c : curves) {
    const bool num = c.accepted && c.boundary == BoundaryId::NumericalEmitterBackClosure;
    std::fprintf(f, "%s,%s,%s,%s,%d,%s,%s,%s,%s\n", c.curve.c_str(),
                 to_string(c.boundary), to_string(c.side_a), to_string(c.side_b), c.n_elements,
                 c.accepted ? "yes" : "no", c.accepted ? p2a_label(c.tag) : "-",
                 num ? "yes" : "no",
                 c.accepted ? (num ? "numerical closure of the conductor, not a device surface"
                                   : "vacuum-facing conductor surface")
                            : c.reason.c_str());
  }
  std::fclose(f);
}

// ---------------------------------------------------------------------------

std::vector<Vec2> open_arc_ends(const Mesh& bem_mesh) {
  if (bem_mesh.elems.empty()) return {};
  // The node tolerance comes from the extent of the mesh itself, so the answer
  // cannot be tuned by choosing a parameter until the check passes.
  Real rmin = bem_mesh.elems.front().a.r, rmax = rmin;
  Real zmin = bem_mesh.elems.front().a.z, zmax = zmin;
  for (const Element& e : bem_mesh.elems)
    for (Vec2 x : {e.a, e.b}) {
      rmin = std::min(rmin, x.r); rmax = std::max(rmax, x.r);
      zmin = std::min(zmin, x.z); zmax = std::max(zmax, x.z);
    }
  const Real snap = 1e-9 * std::max(rmax - rmin, zmax - zmin);

  std::vector<std::pair<Vec2, int>> ends;
  auto bump = [&](Vec2 x) {
    for (auto& e : ends)
      if (norm(e.first - x) <= snap) { ++e.second; return; }
    ends.push_back({x, 1});
  };
  for (const Element& e : bem_mesh.elems) { bump(e.a); bump(e.b); }

  std::vector<Vec2> out;
  for (const auto& e : ends)
    if (e.second == 1 && e.first.r > snap) out.push_back(e.first);
  return out;
}

// ---------------------------------------------------------------------------

Mesh vacuum_bem_mesh(const BoundaryMesh& bm, const DeviceGeometry& g,
                     VacuumSelectionReport* report) {
  VacuumSelectionReport rep;
  rep.n_mesh_elements = static_cast<int>(bm.elements().size());

  // One decision per boundary curve, so the report reads like the geometry
  // does.  The per-element loop below must agree with it exactly.
  const std::vector<BoundaryCurve>& curves = g.boundaries();
  std::vector<VacuumSelectionReport::CurveDecision> decisions(curves.size());
  for (std::size_t k = 0; k < curves.size(); ++k) {
    decisions[k].curve = curves[k].name;
    decisions[k].boundary = curves[k].id;
    decisions[k].side_a = curves[k].side_a;
    decisions[k].side_b = curves[k].side_b;
  }

  Mesh mesh;
  const std::vector<BoundaryElement>& src = bm.elements();
  for (std::size_t i = 0; i < src.size(); ++i) {
    const BoundaryElement& be = src[i];
    VacuumSelectionReport::CurveDecision* dec =
        (be.curve >= 0 && be.curve < static_cast<int>(decisions.size()))
            ? &decisions[static_cast<std::size_t>(be.curve)]
            : nullptr;
    if (dec) ++dec->n_elements;

    Region material = Region::Vacuum;
    Tag tag = Tag::Other;
    const bool vac = !be.is_axis() && separates_vacuum(be, &material);
    if (!vac || !is_conductor(material, be.id, &tag)) {
      if (dec && dec->reason.empty()) dec->reason = reject_reason(be);
      continue;
    }
    require_consistent(be.id, material, dec ? dec->curve : std::string("<unnamed>"));
    if (dec) { dec->accepted = true; dec->tag = tag; }

    // Orientation: the BEM's convention is outward normal = perp(tangent), and
    // "outward" means into the vacuum.  BoundaryElement::normal points from
    // side_a into side_b, so it must be flipped when side_a is the vacuum.  The
    // endpoints are then ordered to match, instead of trusting the mesher's
    // traversal direction.
    const Vec2 n_out = (be.side_b == Region::Vacuum) ? be.normal : -1.0 * be.normal;
    Vec2 a = be.a, b = be.b;
    if (dot(perp(normalized(b - a)), n_out) < 0.0) std::swap(a, b);

    Element el;
    el.a = a;
    el.b = b;
    el.mid = 0.5 * (a + b);
    el.tangent = normalized(b - a);
    el.normal = perp(el.tangent);
    el.len = be.meridian_length;
    el.area = be.revolved_area;  // taken from the mesher, not recomputed
    el.tag = tag;
    el.potential = 0.0;
    // One body per conductor: the emitter metal and the liquid surface it wets
    // form a single closed contour in the meridian plane, the electrode a
    // second one.  io.cpp's inside-conductor test counts crossings per body.
    el.body = (tag == Tag::Extractor) ? 1 : 0;

    VacuumPanel panel;
    panel.mesh_element = static_cast<int>(i);
    panel.bem_element = static_cast<int>(mesh.elems.size());
    panel.boundary = be.id;
    panel.material = material;
    panel.tag = tag;
    panel.numerical = (be.id == BoundaryId::NumericalEmitterBackClosure);
    panel.evaluable = false;  // decided by mark_evaluable_panels(), never guessed
    rep.panels.push_back(panel);

    if (panel.numerical) {
      ++rep.n_numerical_closure;
      rep.revolved_area_numerical_closure += el.area;
    }
    switch (tag) {
      case Tag::Emitter: ++rep.n_emitter; rep.revolved_area_emitter += el.area; break;
      case Tag::FreeSurface:
        ++rep.n_free_surface;
        rep.revolved_area_free_surface += el.area;
        break;
      default: ++rep.n_extractor; rep.revolved_area_extractor += el.area; break;
    }
    mesh.elems.push_back(el);
  }

  rep.n_selected = static_cast<int>(mesh.elems.size());
  if (rep.n_selected == 0)
    throw std::runtime_error("vacuum_bem_mesh: no vacuum-facing conductor panel selected");
  for (VacuumSelectionReport::CurveDecision& d : decisions)
    if (!d.accepted && d.reason.empty()) d.reason = "keine Elemente";
  rep.curves.assign(decisions.begin(), decisions.end());

  // The conductor must be closed.  Every arc endpoint is shared by two panels,
  // or it lies on r = 0, where the surface of revolution closes on itself.
  // Anything else is a free sheet edge, and a free sheet edge is not something
  // to note in a report -- it is a different boundary-value problem.
  for (const Vec2& x : open_arc_ends(mesh)) {
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "vacuum_bem_mesh: the conductor surface is OPEN at r = %.6g m, z = %.6g m.  "
                  "A single layer on an open sheet carries the SUM of both faces, so "
                  "sigma/eps0 is not a one-sided vacuum field, and the free edge adds a "
                  "1/sqrt(d) density singularity that is an artefact of the cut.  Close the "
                  "emitter conductor by setting device.emitter_back_length > 0 "
                  "(a numerical rearward continuation with a conducting end cap); do not "
                  "close it with the open domain edge, which is not a conductor.",
                  x.r, x.z);
    throw std::runtime_error(buf);
  }

  if (report) *report = std::move(rep);
  return mesh;
}

// ---------------------------------------------------------------------------

Real CapacitanceMatrix::reciprocity_error() const {
  const Real s = std::max(std::abs(c_EX), std::abs(c_XE));
  return s > 0.0 ? std::abs(c_EX - c_XE) / s : 0.0;
}

namespace {
/// Charge carried by an electrode for a given density vector.
Real charge_of(const Mesh& m, const std::vector<Real>& sigma, Electrode e) {
  Real q = 0.0;
  for (std::size_t i = 0; i < m.elems.size(); ++i)
    if (electrode_of(m.elems[i].tag) == e) q += sigma[i] * m.elems[i].area;
  return q;
}
}  // namespace

CapacitanceMatrix maxwell_capacitance(BemSolver& bem) {
  const std::vector<Real> sE = bem.sigma_for({{1.0, 0.0, 0.0}});
  const std::vector<Real> sX = bem.sigma_for({{0.0, 1.0, 0.0}});
  CapacitanceMatrix c;
  c.c_EE = charge_of(bem.mesh(), sE, Electrode::Emitter);
  c.c_XE = charge_of(bem.mesh(), sE, Electrode::Extractor);
  c.c_EX = charge_of(bem.mesh(), sX, Electrode::Emitter);
  c.c_XX = charge_of(bem.mesh(), sX, Electrode::Extractor);
  return c;
}

Real electrode_charge(const BemSolver& bem, Electrode e) {
  return charge_of(bem.mesh(), bem.sigma(), e);
}

// ---------------------------------------------------------------------------

const char* to_string(EdgeKind k) {
  switch (k) {
    case EdgeKind::SharpFeature: return "sharp_feature";
    case EdgeKind::TruncationEnd: return "truncation_end";
    case EdgeKind::NumericalClosure: return "numerical_closure";
  }
  return "sharp_feature";
}

PotentialResidual potential_residual(const BemSolver& bem, const std::vector<EdgeZone>& zones,
                                     const std::vector<char>& evaluable) {
  PotentialResidual r;
  const Mesh& m = bem.mesh();
  const std::array<Real, 3>& V = bem.applied();
  r.reference_span = std::abs(V[0] - V[1]);
  if (r.reference_span == 0.0)
    r.reference_span = std::max(std::abs(V[0]), std::abs(V[1]));

  Real se = 0.0, sx = 0.0, sp = 0.0;
  int ne = 0, nx = 0, np = 0;
  const Real ts[2] = {0.25, 0.75};
  for (std::size_t i = 0; i < m.elems.size(); ++i) {
    const Element& el = m.elems[i];
    const Electrode e = electrode_of(el.tag);
    const bool phys = i < evaluable.size() && evaluable[i];
    for (Real t : ts) {
      const Vec2 x = el.a + t * (el.b - el.a);
      const Real d = std::abs(bem.potential_at(x) - el.potential);
      const bool clear = !in_edge_zone(zones, x);
      if (phys && clear) {
        r.max_physical = std::max(r.max_physical, d);
        sp += d * d;
        ++np;
      }
      if (d > std::max(r.max_emitter, r.max_extractor)) r.worst_position = x;
      if (e == Electrode::Emitter) {
        r.max_emitter = std::max(r.max_emitter, d);
        if (clear) r.max_emitter_clear = std::max(r.max_emitter_clear, d);
        se += d * d;
        ++ne;
      } else {
        r.max_extractor = std::max(r.max_extractor, d);
        if (clear) r.max_extractor_clear = std::max(r.max_extractor_clear, d);
        sx += d * d;
        ++nx;
      }
    }
  }
  r.rms_emitter = ne ? std::sqrt(se / ne) : 0.0;
  r.rms_extractor = nx ? std::sqrt(sx / nx) : 0.0;
  r.rms_physical = np ? std::sqrt(sp / np) : 0.0;
  return r;
}

// ---------------------------------------------------------------------------

namespace {
/// Local feature size the mesher assigned to the geometric vertex nearest x.
/// A property of the geometry, so a radius derived from it is the same at every
/// refinement level.
Real local_feature_size_near(const BoundaryMesh& bm, Vec2 x) {
  Real lfs = 0.0, best = std::numeric_limits<Real>::infinity();
  for (const SizeField::Source& s : bm.size_field().sources()) {
    const Real d = norm(s.x - x);
    if (d < best) { best = d; lfs = s.local_feature_size; }
  }
  return lfs;
}
}  // namespace

std::vector<EdgeZone> edge_zones(const DeviceGeometry& g, const BoundaryMesh& bm,
                                 const Mesh& bem_mesh) {
  std::vector<EdgeZone> out;
  for (const NamedFeature& nf : g.features()) {
    EdgeZone z;
    // The rim of the numerical closure is a sharp edge too, but it is excluded
    // for a stronger reason than sharpness: the surface it sits on does not
    // exist on the device.  Keeping the two kinds apart keeps the figures
    // honest about which exclusion is a numerical caveat and which is a
    // statement that the surface is not real.
    z.kind = (nf.id == FeatureId::NumericalBackClosureEdge) ? EdgeKind::NumericalClosure
                                                            : EdgeKind::SharpFeature;
    z.name = to_string(nf.id);
    z.position = nf.position;
    z.local_feature_size = local_feature_size_near(bm, nf.position);
    z.radius = 0.25 * z.local_feature_size;
    out.push_back(z);
  }

  // Free edges of an open conductor arc.  vacuum_bem_mesh() refuses a mesh that
  // has any, so a solved problem produces none; the marking stays because the
  // detector is what proves that, and because edge_zones() may be handed a mesh
  // that never went through the adapter.
  for (const Vec2& x : open_arc_ends(bem_mesh)) {
    EdgeZone z;
    z.kind = EdgeKind::TruncationEnd;
    z.name = "truncation_end";
    z.position = x;
    z.local_feature_size = local_feature_size_near(bm, x);
    z.radius = 0.25 * z.local_feature_size;
    out.push_back(z);
  }
  return out;
}

// ---------------------------------------------------------------------------

std::vector<char> mark_evaluable_panels(const Mesh& bem_mesh, VacuumSelectionReport& report,
                                        const DeviceGeometry& g,
                                        const std::vector<EdgeZone>& zones) {
  std::vector<char> ok(static_cast<std::size_t>(bem_mesh.size()), 0);
  const Real z_eval = g.evaluation_z_min();
  for (VacuumPanel& pn : report.panels) {
    const std::size_t k = static_cast<std::size_t>(pn.bem_element);
    if (pn.bem_element < 0 || k >= ok.size()) continue;
    const Element& el = bem_mesh.elems[k];
    bool good = true;
    // (1) the numerical closure is not a device surface.
    if (pn.numerical) good = false;
    // (2) behind the taper foot the emitter is a shank whose length was chosen
    //     numerically.  Its surface field answers a question about the closure,
    //     not about the device, so it is not offered as a device field either.
    //     The extractor is unaffected: it lies entirely above z = 0.
    if (good && pn.tag != Tag::Extractor && el.mid.z < z_eval) good = false;
    // (3) inside a marked edge zone nothing converges.
    if (good && in_edge_zone(zones, el.mid)) good = false;
    pn.evaluable = good;
    ok[k] = good ? 1 : 0;
  }
  return ok;
}

Real peak_field_evaluable(const BemSolver& bem, Electrode e, const std::vector<char>& ok,
                          Index* which) {
  Real best = 0.0;
  Index arg = -1;
  for (Index i = 0; i < bem.size(); ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    if (k >= ok.size() || !ok[k]) continue;
    if (electrode_of(bem.mesh().elems[k].tag) != e) continue;
    const Real v = std::abs(bem.En(i));
    if (v > best) { best = v; arg = i; }
  }
  if (which) *which = arg;
  return best;
}

bool in_edge_zone(const std::vector<EdgeZone>& zones, Vec2 x) {
  for (const EdgeZone& z : zones)
    if (norm(x - z.position) < z.radius) return true;
  return false;
}

Real peak_field_outside_edges(const BemSolver& bem, Electrode e,
                              const std::vector<EdgeZone>& zones, Index* which) {
  Real best = 0.0;
  Index arg = -1;
  for (Index i = 0; i < bem.size(); ++i) {
    const Element& el = bem.mesh().elems[static_cast<std::size_t>(i)];
    if (electrode_of(el.tag) != e) continue;
    if (in_edge_zone(zones, el.mid)) continue;
    const Real v = std::abs(bem.En(i));
    if (v > best) { best = v; arg = i; }
  }
  if (which) *which = arg;
  return best;
}

}  // namespace es
