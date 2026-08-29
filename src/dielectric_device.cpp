#include "es/dielectric_device.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::eps0;

const char* to_string(ConductorModel c) {
  return c == ConductorModel::Dielectric ? "dielectric" : "metallic_reference";
}

const char* to_string(Metallisation m) {
  switch (m) {
    case Metallisation::FrontOnly: return "front_only";
    case Metallisation::FrontAndAperture: return "front_and_aperture";
    case Metallisation::AllSurfaces: return "all_surfaces";
  }
  return "front_and_aperture";
}

const char* to_string(NodeRole r) {
  switch (r) {
    case NodeRole::Free: return "free";
    case NodeRole::LiquidConductor: return "liquid_conductor";
    case NodeRole::LiquidFeedBoundary: return "liquid_feed_boundary";
    case NodeRole::ExtractorMetallisation: return "extractor_metallisation";
    case NodeRole::FarFieldDirichlet: return "far_field_dirichlet";
    case NodeRole::EmitterMetalReference: return "emitter_metal_reference";
  }
  return "free";
}

namespace {

bool is_fixed(NodeRole r) { return r != NodeRole::Free; }

/// Does the metallisation option coat this face?  The emitter-facing face is
/// coated in every option -- an extractor whose front is bare is not an
/// extraction electrode.
bool coated_front(Metallisation) { return true; }
bool coated_aperture(Metallisation m) { return m != Metallisation::FrontOnly; }
bool coated_back(Metallisation m) { return m == Metallisation::AllSurfaces; }
bool coated_rim(Metallisation m) { return m == Metallisation::AllSurfaces; }

}  // namespace

// ---------------------------------------------------------------------------

std::vector<NodeRole> node_roles(const DeviceVolumeMesh& m, const DielectricSetup& s) {
  const QuadMesh& g = m.grid;
  std::vector<NodeRole> role(static_cast<std::size_t>(g.n_nodes()), NodeRole::Free);
  auto set = [&](Index i, Index j, NodeRole r) {
    role[static_cast<std::size_t>(g.node(i, j))] = r;
  };

  // --- the ionic liquid: ONE equipotential body ------------------------------
  //
  // Its closure is the Dirichlet set, and the closure is taken from the CELL
  // REGIONS, not from a list of surfaces: a node belongs to the liquid
  // conductor exactly when one of the cells around it is liquid.  That way the
  // bore, the feed channel and the plenum are one body by construction, no
  // surface can be forgotten, and adding another liquid-filled part later needs
  // no change here.
  for (Index j = 0; j + 1 < g.nz; ++j)
    for (Index i = 0; i + 1 < g.nr; ++i) {
      if (m.cell_region[static_cast<std::size_t>(g.cell(i, j))] != Region::Liquid) continue;
      set(i, j, NodeRole::LiquidConductor);
      set(i + 1, j, NodeRole::LiquidConductor);
      set(i, j + 1, NodeRole::LiquidConductor);
      set(i + 1, j + 1, NodeRole::LiquidConductor);
    }
  // The rear cut of the superseded truncated column is labelled separately, so
  // that the figures and the audit can show that the emitter potential stops at
  // the liquid cross section and does not continue across the polymer.
  if (m.p.reservoir == ReservoirModel::TruncatedColumn)
    for (Index i = 0; i <= m.i_bore; ++i) set(i, m.j_base, NodeRole::LiquidFeedBoundary);

  // --- the metallisation on the extractor carrier ---------------------------
  const Metallisation mt =
      (s.conductor_model == ConductorModel::MetallicReference) ? Metallisation::AllSurfaces
                                                               : s.metallisation;
  if (coated_front(mt))
    for (Index i = m.i_aperture; i <= m.i_ext_outer; ++i)
      set(i, m.j_ex_front, NodeRole::ExtractorMetallisation);
  if (coated_aperture(mt))
    for (Index j = m.j_ex_front; j <= m.j_ex_back; ++j)
      set(m.i_aperture, j, NodeRole::ExtractorMetallisation);
  if (coated_back(mt))
    for (Index i = m.i_aperture; i <= m.i_ext_outer; ++i)
      set(i, m.j_ex_back, NodeRole::ExtractorMetallisation);
  if (coated_rim(mt))
    for (Index j = m.j_ex_front; j <= m.j_ex_back; ++j)
      set(m.i_ext_outer, j, NodeRole::ExtractorMetallisation);
  // A coating on every face does NOT turn the carrier into metal: the interior
  // stays a dielectric whose nodes are ordinary unknowns, and they come out at
  // V_extractor because they are enclosed.  Fixing them by hand would make the
  // polymer audit below vacuous.  The one case in which the whole block must be
  // fixed is the metallic reference arrangement, where those cells really are
  // conductor interior and carry no equation at all.
  if (s.conductor_model == ConductorModel::MetallicReference)
    for (Index j = m.j_ex_front; j <= m.j_ex_back; ++j)
      for (Index i = m.i_aperture; i <= m.i_ext_outer; ++i)
        set(i, j, NodeRole::ExtractorMetallisation);

  // --- the superseded metallic arrangement, for the BEM cross-check only -----
  if (s.conductor_model == ConductorModel::MetallicReference)
    for (Index j = m.j_base; j <= m.j_tip; ++j)
      for (Index i = 0; i <= m.i_land; ++i)
        if (!is_fixed(role[static_cast<std::size_t>(g.node(i, j))]) ||
            role[static_cast<std::size_t>(g.node(i, j))] == NodeRole::LiquidConductor)
          set(i, j, NodeRole::EmitterMetalReference);

  // --- a grounded enclosure, when that is what was asked for ----------------
  if (s.far_field == FarField::Grounded) {
    for (Index i = 0; i < g.nr; ++i) {
      set(i, 0, NodeRole::FarFieldDirichlet);
      set(i, g.nz - 1, NodeRole::FarFieldDirichlet);
    }
    for (Index j = 0; j < g.nz; ++j) set(g.nr - 1, j, NodeRole::FarFieldDirichlet);
  }
  return role;
}

// ---------------------------------------------------------------------------

BoundaryAudit audit_boundaries(const DeviceVolumeMesh& m, const std::vector<NodeRole>& role,
                               const DielectricSetup& s) {
  const QuadMesh& g = m.grid;
  BoundaryAudit a;
  a.n_nodes = g.n_nodes();
  for (NodeRole r : role) {
    if (is_fixed(r)) ++a.n_dirichlet;
    switch (r) {
      case NodeRole::LiquidConductor: ++a.n_liquid; break;
      case NodeRole::LiquidFeedBoundary: ++a.n_feed; break;
      case NodeRole::ExtractorMetallisation: ++a.n_metal; break;
      case NodeRole::FarFieldDirichlet: ++a.n_far; break;
      case NodeRole::EmitterMetalReference: ++a.n_emitter_metal_reference; break;
      default: break;
    }
  }

  auto at = [&](Index i, Index j) { return role[static_cast<std::size_t>(g.node(i, j))]; };
  auto complain = [&a](const std::string& text) { a.violations.push_back(text); };

  // Which nodes touch liquid.  Everything below is expressed through this, so
  // no surface has to be enumerated for the check that decides.
  std::vector<char> touches_liquid(static_cast<std::size_t>(g.n_nodes()), 0);
  for (Index j = 0; j + 1 < g.nz; ++j)
    for (Index i = 0; i + 1 < g.nr; ++i) {
      if (m.cell_region[static_cast<std::size_t>(g.cell(i, j))] != Region::Liquid) continue;
      for (Index dj = 0; dj < 2; ++dj)
        for (Index di = 0; di < 2; ++di)
          touches_liquid[static_cast<std::size_t>(g.node(i + di, j + dj))] = 1;
    }

  // The rear cut plane of the superseded truncated column is an electrode ONLY
  // on the liquid cross section.  Every other node of that plane must be free.
  if (m.p.reservoir == ReservoirModel::TruncatedColumn) {
    for (Index i = m.i_bore + 1; i < g.nr; ++i)
      if (is_fixed(at(i, m.j_base)) && at(i, m.j_base) != NodeRole::FarFieldDirichlet &&
          at(i, m.j_base) != NodeRole::EmitterMetalReference) {
        ++a.n_feed_plane_outside_liquid;
      }
    if (a.n_feed_plane_outside_liquid > 0)
      complain("Die Schnittebene am Ende der abgeschnittenen Saeule traegt " +
               std::to_string(static_cast<long long>(a.n_feed_plane_outside_liquid)) +
               " festgehaltene Knoten ausserhalb des Fluessigkeitsquerschnitts. Die Ebene "
               "ist keine Elektrode.");
  }

  if (s.conductor_model != ConductorModel::Dielectric) {
    // In the metallic reference arrangement the emitter body IS a conductor on
    // purpose.  The polymer audit does not apply, and saying so is safer than
    // running a check that silently passes for the wrong reason.
    a.n_polymer_dirichlet = -1;
    a.n_named_surface_dirichlet = -1;
    return a;
  }

  // --- the structural check -------------------------------------------------
  for (Index n = 0; n < g.n_nodes(); ++n) {
    const NodeRole r = role[static_cast<std::size_t>(n)];
    if (!is_fixed(r)) continue;
    if (r == NodeRole::ExtractorMetallisation || r == NodeRole::FarFieldDirichlet) continue;
    if (!touches_liquid[static_cast<std::size_t>(n)]) ++a.n_polymer_dirichlet;
  }
  if (a.n_polymer_dirichlet > 0)
    complain("Es tragen " + std::to_string(static_cast<long long>(a.n_polymer_dirichlet)) +
             " festgehaltene Knoten das Emitterpotential, ohne an eine Fluessigkeitszelle zu "
             "grenzen. Ein Dielektrikum ist keine Elektrode, und ein leitender Halter oder "
             "eine rueckwaertige Metallscheibe ist in diesem Modell ausgeschlossen.");

  // Book-keeping, not a check: how much of the liquid surface belongs to the
  // reservoir.  It has to be non-zero with a plenum, or the liquid would have
  // been split into disconnected pieces.
  if (m.has_plenum())
    for (Index n = 0; n < g.n_nodes(); ++n) {
      if (role[static_cast<std::size_t>(n)] != NodeRole::LiquidConductor) continue;
      const Index j = n / g.nr;
      if (j < m.j_base) ++a.n_reservoir_liquid_surface;
    }

  // --- the named surfaces ---------------------------------------------------
  struct Surface {
    const char* name;
    Index i0, i1, j0, j1;
  };
  const Metallisation mt = s.metallisation;
  std::vector<Surface> polymer{
      {"emitter_outer_surface", m.i_land, m.i_land, m.j_base, m.j_tip},
      {"emitter_tip_land", m.i_bore + 1, m.i_land, m.j_tip, m.j_tip},
      {"emitter_rear_face", m.i_bore + 1, m.i_land, m.j_base, m.j_base},
      {"emitter_body_interior", m.i_bore + 1, m.i_land - 1, m.j_base + 1, m.j_tip - 1},
  };
  if (m.has_plenum()) {
    polymer.push_back({"reservoir_top_face", m.i_land, m.i_plenum_outer, m.j_base, m.j_base});
    polymer.push_back({"reservoir_outer_rim", m.i_plenum_outer, m.i_plenum_outer,
                       m.j_block_bottom, m.j_base});
    polymer.push_back({"reservoir_underside", m.i_axis, m.i_plenum_outer, m.j_block_bottom,
                       m.j_block_bottom});
    polymer.push_back({"reservoir_top_wall_interior", m.i_channel + 1, m.i_plenum_outer - 1,
                       m.j_roof + 1, m.j_base - 1});
    polymer.push_back({"reservoir_side_wall_interior", m.i_plenum + 1, m.i_plenum_outer - 1,
                       m.j_cav_bottom + 1, m.j_roof - 1});
    polymer.push_back({"reservoir_floor_interior", m.i_axis, m.i_plenum_outer - 1,
                       m.j_block_bottom + 1, m.j_cav_bottom - 1});
    // With a partly filled plenum the cavity floor faces vacuum, not liquid, so
    // it is a polymer surface and must be free.  With a full plenum it is the
    // underside of the liquid and legitimately fixed, so it is not listed.
    if (m.p.plenum_fill_fraction < 1.0)
      polymer.push_back({"plenum_floor_dry", m.i_axis, m.i_plenum, m.j_cav_bottom,
                         m.j_cav_bottom});
  }
  if (!coated_back(mt))
    polymer.push_back({"extractor_back_face", m.i_aperture, m.i_ext_outer, m.j_ex_back,
                       m.j_ex_back});
  if (!coated_rim(mt))
    polymer.push_back({"extractor_rim", m.i_ext_outer, m.i_ext_outer, m.j_ex_front,
                       m.j_ex_back});
  if (!coated_aperture(mt))
    polymer.push_back({"extractor_aperture_wall", m.i_aperture, m.i_aperture, m.j_ex_front,
                       m.j_ex_back});
  polymer.push_back({"extractor_carrier_interior", m.i_aperture + 1, m.i_ext_outer - 1,
                     m.j_ex_front + 1, m.j_ex_back - 1});

  // Two adjacent faces share their corner node, and where a coated face meets an
  // uncoated one that node IS on the metal -- it is the rim of the film.  Those
  // nodes are enumerated first and excluded, so that the audit tests what it
  // means to test (a polymer FACE held at a potential) instead of failing on the
  // one node the two faces have in common.
  std::vector<char> on_metal(static_cast<std::size_t>(g.n_nodes()), 0);
  for (Index n = 0; n < g.n_nodes(); ++n)
    if (role[static_cast<std::size_t>(n)] == NodeRole::ExtractorMetallisation)
      on_metal[static_cast<std::size_t>(n)] = 1;
  auto is_film_rim = [&](Index i, Index j) {
    // A node counts as film rim only if it is a corner or edge of the electrode
    // block that a coated face reaches; interior polymer nodes never are.
    const bool front = (j == m.j_ex_front) && coated_front(mt);
    const bool aperture = (i == m.i_aperture) && coated_aperture(mt);
    const bool back = (j == m.j_ex_back) && coated_back(mt);
    const bool rim = (i == m.i_ext_outer) && coated_rim(mt);
    return front || aperture || back || rim;
  };

  Index shared = 0;
  for (const Surface& sf : polymer) {
    Index bad = 0;
    for (Index j = sf.j0; j <= sf.j1; ++j)
      for (Index i = sf.i0; i <= sf.i1; ++i) {
        if (i < 0 || j < 0 || i >= g.nr || j >= g.nz) continue;
        const NodeRole r = at(i, j);
        if (!is_fixed(r)) continue;
        if (on_metal[static_cast<std::size_t>(g.node(i, j))] && is_film_rim(i, j)) {
          ++shared;
          continue;
        }
        // A node of a named polymer surface that nevertheless touches liquid is
        // a corner shared with the liquid closure -- the rim of the bore at the
        // tip land, the mouth of the feed channel.  It is legitimately fixed and
        // is excluded here for the same reason as the film rim; the structural
        // check above is what would catch a genuinely conducting polymer face.
        if (touches_liquid[static_cast<std::size_t>(g.node(i, j))]) continue;
        // A far-field Dirichlet node on a polymer surface would mean the
        // computational box touches the device; the mesher forbids that, but
        // the check is cheap and the failure would be silent.
        ++bad;
      }
    a.n_named_surface_dirichlet += bad;
    if (bad > 0)
      complain(std::string("Polymerflaeche '") + sf.name + "' traegt " +
               std::to_string(static_cast<long long>(bad)) +
               " festgehaltene Knoten. Ein Dielektrikum ist keine Elektrode.");
  }
  a.n_film_rim_shared = shared;
  return a;
}

void BoundaryAudit::print(std::FILE* out) const {
  std::fprintf(out, "Randbedingungs-Audit\n");
  std::fprintf(out, "  Knoten gesamt                          : %lld\n",
               static_cast<long long>(n_nodes));
  std::fprintf(out, "  davon Dirichlet                        : %lld\n",
               static_cast<long long>(n_dirichlet));
  std::fprintf(out, "    Fluessigkeit (idealer Leiter)        : %lld\n",
               static_cast<long long>(n_liquid));
  std::fprintf(out, "    Zulaufgrenze (Fluessigkeitsquerschn.): %lld\n",
               static_cast<long long>(n_feed));
  std::fprintf(out, "    Extraktormetallisierung              : %lld\n",
               static_cast<long long>(n_metal));
  std::fprintf(out, "    geerdete Huelle (nur far_field=grounded): %lld\n",
               static_cast<long long>(n_far));
  std::fprintf(out, "    metallischer Referenzemitter (KEIN P2b-Modell): %lld\n",
               static_cast<long long>(n_emitter_metal_reference));
  if (n_polymer_dirichlet < 0) {
    std::fprintf(out, "  Polymerflaechen: Pruefung entfaellt (metallische Referenz)\n");
  } else {
    std::fprintf(out, "  festgehaltene Knoten ohne Fluessigkeitskontakt (strukturell): "
                      "%lld  (muss 0 sein)\n",
                 static_cast<long long>(n_polymer_dirichlet));
    std::fprintf(out, "  davon auf benannten Polymerflaechen                       : "
                      "%lld  (muss 0 sein)\n",
                 static_cast<long long>(n_named_surface_dirichlet));
  }
  std::fprintf(out, "  festgehaltene Knoten der Schnittebene ausserhalb der Fluessigkeit: "
                    "%lld  (muss 0 sein)\n",
               static_cast<long long>(n_feed_plane_outside_liquid));
  std::fprintf(out, "  Fluessigkeitsknoten hinter der Grundplatte (Vorrat)            : %lld\n",
               static_cast<long long>(n_reservoir_liquid_surface));
  for (const std::string& v : violations) std::fprintf(out, "  VERLETZUNG: %s\n", v.c_str());
}

void BoundaryAudit::write_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# Audit der Dirichlet-Belegung; polymer_dirichlet muss 0 sein\n");
  std::fprintf(f, "quantity,value\n");
  std::fprintf(f, "nodes,%lld\n", static_cast<long long>(n_nodes));
  std::fprintf(f, "dirichlet,%lld\n", static_cast<long long>(n_dirichlet));
  std::fprintf(f, "liquid_conductor,%lld\n", static_cast<long long>(n_liquid));
  std::fprintf(f, "liquid_feed_boundary,%lld\n", static_cast<long long>(n_feed));
  std::fprintf(f, "extractor_metallisation,%lld\n", static_cast<long long>(n_metal));
  std::fprintf(f, "far_field_dirichlet,%lld\n", static_cast<long long>(n_far));
  std::fprintf(f, "emitter_metal_reference,%lld\n",
               static_cast<long long>(n_emitter_metal_reference));
  std::fprintf(f, "polymer_dirichlet,%lld\n", static_cast<long long>(n_polymer_dirichlet));
  std::fprintf(f, "named_surface_dirichlet,%lld\n",
               static_cast<long long>(n_named_surface_dirichlet));
  std::fprintf(f, "reservoir_liquid_nodes,%lld\n",
               static_cast<long long>(n_reservoir_liquid_surface));
  std::fprintf(f, "feed_plane_outside_liquid,%lld\n",
               static_cast<long long>(n_feed_plane_outside_liquid));
  std::fprintf(f, "violations,%lld\n", static_cast<long long>(violations.size()));
  std::fclose(f);
  const std::string txt = path + ".violations.txt";
  std::FILE* g2 = std::fopen(txt.c_str(), "w");
  if (g2) {
    for (const std::string& v : violations) std::fprintf(g2, "%s\n", v.c_str());
    std::fclose(g2);
  }
}

// ---------------------------------------------------------------------------

std::vector<Probe> dielectric_probes(const DeviceVolumeMesh& m) {
  const DeviceParameters& d = m.p.device;
  const Real L = d.extraction_distance;
  const Real H = d.emitter_height;
  const Real w = m.r_land - m.r_bore;
  std::vector<Probe> out;
  auto add = [&](const char* name, Vec2 x, const char* note) {
    // Distance to the nearest unrounded edge of the device.
    const Vec2 edges[4] = {{m.r_bore, 0.0},
                           {m.r_land, 0.0},
                           {m.r_aperture, L},
                           {m.r_aperture, L + d.extractor_thickness}};
    Real cl = norm(x - edges[0]);
    for (int k = 1; k < 4; ++k) cl = std::min(cl, norm(x - edges[k]));
    out.push_back({name, x, cl, note});
  };
  add("axis_2_bore_radii", {0.0, 2.0 * m.r_bore}, "auf der Achse ueber der Referenzflaeche");
  add("axis_10_bore_radii", {0.0, 10.0 * m.r_bore}, "auf der Achse");
  add("axis_gap_quarter", {0.0, 0.25 * L}, "Achse, ein Viertel der Extraktionsstrecke");
  add("axis_gap_mid", {0.0, 0.50 * L}, "Achse, halbe Extraktionsstrecke");
  add("axis_gap_three_quarter", {0.0, 0.75 * L}, "Achse, drei Viertel");
  add("flank_beside_taper", {m.emitter_outer_radius_at(-0.5 * H) + 2.0 * w, -0.5 * H},
      "im Vakuum neben der Kegelflanke, auf halber Hoehe");
  add("off_axis_half_aperture", {0.5 * m.r_aperture, 0.5 * L}, "ausserhalb der Achse");
  return out;
}

// ---------------------------------------------------------------------------

DielectricSolution solve_dielectric(const DielectricSetup& s) {
  DielectricSolution out;
  out.mesh = build_volume_mesh(s.geometry);
  const DeviceVolumeMesh& m = out.mesh;
  const QuadMesh& g = m.grid;

  s.materials.check_usable();
  const Real eps_emitter = s.materials.emitter_dielectric.permittivity_or_throw();
  const Real eps_carrier = s.materials.extractor_carrier.permittivity_or_throw();
  const Real eps_reservoir = s.materials.reservoir_body.permittivity_or_throw();
  if (s.conductor_model == ConductorModel::MetallicReference &&
      (eps_emitter != 1.0 || eps_carrier != 1.0 || eps_reservoir != 1.0))
    throw std::runtime_error(
        "MetallicReference ist ausschliesslich der Sonderfall eps_r = 1, in dem die "
        "unabhaengige BEM dasselbe Problem loest. Mit eps_r != 1 waere es weder das eine "
        "noch das andere Modell.");

  out.role = node_roles(m, s);
  out.audit = audit_boundaries(m, out.role, s);
  if (!out.audit.ok()) {
    std::string msg = "solve_dielectric: das Randbedingungs-Audit ist fehlgeschlagen.\n";
    for (const std::string& v : out.audit.violations) msg += "  " + v + "\n";
    throw std::runtime_error(msg);
  }

  // --- cells ---------------------------------------------------------------
  AxisymProblem prob;
  prob.mesh = &g;
  prob.eps_r.assign(static_cast<std::size_t>(g.n_cells()), 1.0);
  prob.active.assign(static_cast<std::size_t>(g.n_cells()), 1);
  for (Index c = 0; c < g.n_cells(); ++c) {
    const Region rg = m.cell_region[static_cast<std::size_t>(c)];
    Real e = 1.0;
    bool act = true;
    switch (rg) {
      case Region::Liquid:
        act = false;  // interior of an ideal conductor: not a field region
        break;
      case Region::EmitterSolid:
        if (s.conductor_model == ConductorModel::MetallicReference)
          act = false;
        else
          e = eps_emitter;
        break;
      case Region::ExtractorSolid:
        if (s.conductor_model == ConductorModel::MetallicReference)
          act = false;
        else
          e = eps_carrier;
        break;
      case Region::ReservoirSolid:
        // A DIELECTRIC, in every conductor model.  The metallic reference
        // arrangement exists to reproduce the P2a vacuum problem the BEM can
        // solve, and it forces every eps_r to 1 anyway; turning the reservoir
        // body into metal there would invent a conducting holder, which is the
        // one thing this phase must not do.
        e = eps_reservoir;
        break;
      default:
        e = s.materials.vacuum.relative_permittivity;
        break;
    }
    prob.eps_r[static_cast<std::size_t>(c)] = e;
    prob.active[static_cast<std::size_t>(c)] = act ? 1 : 0;
  }

  // --- nodes ---------------------------------------------------------------
  prob.fixed.assign(static_cast<std::size_t>(g.n_nodes()), 0);
  prob.fixed_value.assign(static_cast<std::size_t>(g.n_nodes()), 0.0);
  out.emitter_mask.assign(static_cast<std::size_t>(g.n_nodes()), 0);
  out.extractor_mask.assign(static_cast<std::size_t>(g.n_nodes()), 0);
  for (Index n = 0; n < g.n_nodes(); ++n) {
    const NodeRole r = out.role[static_cast<std::size_t>(n)];
    if (!is_fixed(r)) continue;
    prob.fixed[static_cast<std::size_t>(n)] = 1;
    switch (r) {
      case NodeRole::LiquidConductor:
      case NodeRole::LiquidFeedBoundary:
      case NodeRole::EmitterMetalReference:
        prob.fixed_value[static_cast<std::size_t>(n)] = s.V_emitter;
        out.emitter_mask[static_cast<std::size_t>(n)] = 1;
        break;
      case NodeRole::ExtractorMetallisation:
        prob.fixed_value[static_cast<std::size_t>(n)] = s.V_extractor;
        out.extractor_mask[static_cast<std::size_t>(n)] = 1;
        break;
      case NodeRole::FarFieldDirichlet:
        prob.fixed_value[static_cast<std::size_t>(n)] = 0.0;
        break;
      default:
        break;
    }
  }

  // --- far field ------------------------------------------------------------
  prob.far_field = s.far_field;
  prob.far_field_origin = {0.0, 0.5 * s.geometry.device.extraction_distance};
  if (s.far_field == FarField::Asymptotic) {
    for (Index i = 0; i + 1 < g.nr; ++i) {
      prob.far_edges.push_back({g.node(i, 0), g.node(i + 1, 0)});
      prob.far_edges.push_back({g.node(i, g.nz - 1), g.node(i + 1, g.nz - 1)});
    }
    for (Index j = 0; j + 1 < g.nz; ++j)
      prob.far_edges.push_back({g.node(g.nr - 1, j), g.node(g.nr - 1, j + 1)});
  }

  out.fem = solve_axisym(prob);
  out.cell_eps_r = prob.eps_r;
  out.cell_active = prob.active;

  // --- charges --------------------------------------------------------------
  out.Q_emitter = charge_of(out.fem, out.emitter_mask);
  out.Q_extractor = charge_of(out.fem, out.extractor_mask);
  out.Q_net = out.Q_emitter + out.Q_extractor;

  // --- probes ---------------------------------------------------------------
  out.probes = dielectric_probes(m);
  for (const Probe& p : out.probes) {
    out.phi_probe.push_back(potential_at(g, out.fem.phi, p.x));
    const Vec2 E = field_recovered_at(g, out.fem.phi, out.cell_eps_r, out.cell_active, p.x);
    out.Ez_probe.push_back(E.z);
    out.Emag_probe.push_back(norm(E));
  }

  // --- one-sided normal field on the flat liquid reference plane -------------
  //
  // Evaluated in the VACUUM cell above z = 0 at eta = 0, i.e. the limit from the
  // vacuum side.  The last cells before the pinned edge sit inside the wedge
  // singularity of the unrounded exit edge; they are counted and excluded, not
  // quoted.  Nothing here is a converged peak field, and nothing claims to be.
  {
    const Index j = m.j_tip;
    const Real edge_zone = 0.25 * (m.r_land - m.r_bore);
    const Real eps_vac = s.materials.vacuum.relative_permittivity;
    for (Index i = 0; i <= m.i_bore; ++i) {
      const Real rn = g.at(i, j).r;
      if (m.r_bore - rn < edge_zone) {
        ++out.surface_edge_cells;
        continue;
      }
      // Nodal recovery restricted to VACUUM cells: the one-sided limit from
      // above, second order, and free of the jitter a raw cell gradient has.
      const Vec2 E = field_recovered_at_node(g, out.fem.phi, out.cell_eps_r, out.cell_active, i,
                                             j, eps_vac);
      out.surface_r.push_back(rn);
      out.surface_Ez.push_back(E.z);
    }
  }

  // --- continuity across the polymer/vacuum interface on the taper flank -----
  //
  // phi is continuous by construction (one nodal unknown on the interface).
  // The normal flux density is NOT: it is continuous only if the solution is
  // right, so it is the check with content.  Sampled at half the taper height,
  // one cell either side of the interface grid line.
  {
    const Index j = (m.j_foot + m.j_tip) / 2;
    const Index i = m.i_land;
    if (i >= 1 && i + 1 < g.nr && j + 1 < g.nz) {
      // Outward normal of the flank, from the meridian tangent of the interface.
      const Vec2 a{m.warp(m.r_land, g.z_of_row(j)), g.z_of_row(j)};
      const Vec2 b{m.warp(m.r_land, g.z_of_row(j + 1)), g.z_of_row(j + 1)};
      Vec2 n = normalized(perp(b - a));
      if (n.r < 0.0) n = -1.0 * n;  // point into the vacuum
      const Vec2 E_in = field_in_cell(g, out.fem.phi, i - 1, j, 1.0, 0.5);
      const Vec2 E_out = field_in_cell(g, out.fem.phi, i, j, 0.0, 0.5);
      out.Dn_polymer_side = eps0 * eps_emitter * dot(E_in, n);
      out.Dn_vacuum_side = eps0 * s.materials.vacuum.relative_permittivity * dot(E_out, n);
      out.phi_interface_jump =
          std::abs(potential_in_cell(g, out.fem.phi, i - 1, j, 1.0, 0.5) -
                   potential_in_cell(g, out.fem.phi, i, j, 0.0, 0.5));
    }
  }
  return out;
}

Real DielectricSolution::relative_interface_error() const {
  const Real s = std::max(std::abs(Dn_polymer_side), std::abs(Dn_vacuum_side));
  return s > 0.0 ? std::abs(Dn_polymer_side - Dn_vacuum_side) / s : 0.0;
}

void DielectricSolution::print(std::FILE* out) const {
  mesh.print(out);
  audit.print(out);
  std::fprintf(out, "Loesung\n");
  std::fprintf(out, "  freie Unbekannte             : %lld von %lld\n",
               static_cast<long long>(fem.n_free), static_cast<long long>(fem.n_nodes));
  std::fprintf(out, "  Halbbandbreite               : %lld, Faktor %.2f MiB\n",
               static_cast<long long>(fem.half_bandwidth), fem.factor_bytes / 1048576.0);
  std::fprintf(out, "  Residuum der freien Gleichungen: %.3e C\n", fem.residual_inf);
  std::fprintf(out, "  Q_emitter = %.6e C, Q_extractor = %.6e C, Summe = %.6e C\n", Q_emitter,
               Q_extractor, Q_net);
  std::fprintf(out, "  D_n am Dielektrikumsuebergang: innen %.6e, aussen %.6e, relativ %.2e\n",
               Dn_polymer_side, Dn_vacuum_side, relative_interface_error());
  std::fprintf(out, "  %-28s %11s %11s %11s %13s %13s\n", "Punkt", "r [m]", "z [m]",
               "Abstand [m]", "phi [V]", "|E| [V/m]");
  for (std::size_t k = 0; k < probes.size(); ++k)
    std::fprintf(out, "  %-28s %11.4g %11.4g %11.4g %13.6g %13.6g\n", probes[k].name.c_str(),
                 probes[k].x.r, probes[k].x.z, probes[k].clearance, phi_probe[k], Emag_probe[k]);
}

void DielectricSolution::write_csv(const std::string& dir) const {
  const std::string d = dir + "/";
  {
    std::FILE* f = std::fopen((d + "probes.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open probes.csv");
    std::fprintf(f, "# feste Auswertepunkte; clearance = Abstand zur naechsten unverrundeten "
                    "Kante\n");
    std::fprintf(f, "name,r_m,z_m,clearance_m,phi_V,Ez_V_per_m,Emag_V_per_m,note\n");
    for (std::size_t k = 0; k < probes.size(); ++k)
      std::fprintf(f, "%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,\"%s\"\n", probes[k].name.c_str(),
                   probes[k].x.r, probes[k].x.z, probes[k].clearance, phi_probe[k], Ez_probe[k],
                   Emag_probe[k], probes[k].note.c_str());
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "reference_surface_field.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open reference_surface_field.csv");
    std::fprintf(f, "# einseitiges E_z unmittelbar ueber der ebenen Fluessigkeitsreferenz "
                    "(z = 0+)\n");
    std::fprintf(f, "# %lld kantennahe Zellen sind ausgeschlossen: die unverrundete "
                    "Austrittskante\n# hat dort ein divergierendes Feld, das der "
                    "Elementgroesse folgt.\n",
                 static_cast<long long>(surface_edge_cells));
    std::fprintf(f, "r_m,Ez_V_per_m\n");
    for (std::size_t k = 0; k < surface_r.size(); ++k)
      std::fprintf(f, "%.9e,%.9e\n", surface_r[k], surface_Ez[k]);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "node_roles.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open node_roles.csv");
    std::fprintf(f, "# nur die festgehaltenen Knoten; alle uebrigen sind frei\n");
    std::fprintf(f, "r_m,z_m,role\n");
    for (Index j = 0; j < mesh.grid.nz; ++j)
      for (Index i = 0; i < mesh.grid.nr; ++i) {
        const NodeRole r = role[static_cast<std::size_t>(mesh.grid.node(i, j))];
        if (r == NodeRole::Free) continue;
        std::fprintf(f, "%.9e,%.9e,%s\n", mesh.grid.at(i, j).r, mesh.grid.at(i, j).z,
                     to_string(r));
      }
    std::fclose(f);
  }
  audit.write_csv(d + "boundary_audit.csv");
}

}  // namespace es
