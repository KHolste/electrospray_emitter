// P2c: the axisymmetric plenum -- the replacement for the truncated liquid
// column, and the separation of the parameters that used to be tangled.
//
// The defect this file exists to prevent is specific.  Until now a single
// number, the position of the "liquid feed boundary", moved the length of the
// CONDUCTING liquid column, the rearward extent of the DIELECTRIC body and the
// whole rear geometry of the device at the same time.  A study that varied it
// was therefore not moving a boundary condition; it was building a different
// high-voltage electrode each time, and the large field change it produced was
// reported as a failed convergence "against the position of the feed boundary".
//
// So the checks below are not about the plenum being pretty.  They are:
//
//   1  resizing the reservoir moves NO node of the front device, bitwise;
//   2  and changes NO cell's material in front of the base body;
//   3  every liquid region -- bore, feed channel, plenum -- sits at exactly
//      V_emitter, and they are ONE connected conductor;
//   4  no polymer or PEEK node carries a Dirichlet electrode condition, checked
//      structurally (a fixed node must touch liquid) and not by a surface list;
//   5  the P2a conducting back disc stays refused, and the plenum does not
//      smuggle one back in;
//   6  the mesh still represents every region volume exactly, from closed-form
//      formulae that know nothing about the mesh.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/dielectric_device.hpp"

using namespace es;
using constants::pi;

static int failures = 0;

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-62s %s\n", what, ok ? "OK" : "FAIL");
}

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-46s got=%-12.6g want=%-12.6g rel=%-9.3g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

// ---------------------------------------------------------------------------
// The fixed front geometry.  Every test below varies the reservoir against it
// and nothing else.
// ---------------------------------------------------------------------------

static DielectricDeviceParameters plenum_geometry(Real r_plenum, Real depth, int level = 1) {
  DielectricDeviceParameters p;
  p.device.emitter_back_length = 0.0;
  p.device.domain_radius = 1.2e-2;
  p.device.domain_z_min = -1.2e-2;
  p.device.domain_z_max = 9.0e-3;
  p.base_plate_thickness = 1.4e-4;    // rear face of the printed body at -200 um
  p.reservoir = ReservoirModel::AxisymmetricPlenum;
  p.feed_channel_radius = 0.0;        // = bore radius
  p.feed_channel_length = 3.0e-4;
  p.plenum_radius = r_plenum;
  p.plenum_depth = depth;
  p.plenum_wall_thickness = 5.0e-4;
  p.plenum_fill_fraction = 1.0;
  p.mesh_level = level;
  return p;
}

static DielectricSetup plenum_setup(const MaterialLibrary& lib, Real r_plenum, Real depth,
                                    int level = 1) {
  DielectricSetup s;
  s.geometry = plenum_geometry(r_plenum, depth, level);
  s.materials = DielectricMaterials::reference(lib);
  s.conductor_model = ConductorModel::Dielectric;
  s.metallisation = Metallisation::FrontAndAperture;
  s.far_field = FarField::Asymptotic;
  s.V_emitter = 1500.0;
  s.V_extractor = 0.0;
  return s;
}

// ==========================================================================
// 1.  The mesh represents the replacement geometry exactly
// ==========================================================================
static void test_plenum_mesh() {
  std::printf("\n=== 1. Netz des achsensymmetrischen Plenums ===\n");
  const DeviceVolumeMesh m = build_volume_mesh(plenum_geometry(2.5e-3, 1.0e-3, 1));
  m.print(stdout);
  const VolumeMeshReport rep = m.validate();
  rep.print(stdout);
  expect("alle Netzpruefungen bestanden", rep.all_passed());

  // Closed-form volumes, written out here independently of analytic_volume_of.
  const Real a = m.r_bore, ac = m.r_channel;
  const Real v_liquid = pi * a * a * (-m.z_base) + pi * ac * ac * (m.z_base - m.z_roof) +
                        pi * m.r_plenum * m.r_plenum * (m.z_roof - m.z_cav_bottom);
  check("Rotationsvolumen der gesamten Fluessigkeit",
        m.revolved_volume_of(Region::Liquid), v_liquid, 1.0e-11);
  const Real Ro = m.r_plenum_outer;
  const Real v_res = pi * Ro * Ro * (m.z_base - m.z_block_bottom) -
                     pi * m.r_plenum * m.r_plenum * (m.z_roof - m.z_cav_bottom) -
                     pi * ac * ac * (m.z_base - m.z_roof);
  check("Rotationsvolumen des Vorratskoerpers", m.revolved_volume_of(Region::ReservoirSolid),
        v_res, 1.0e-11);

  // Determinism: same parameters, bitwise same mesh.
  const DeviceVolumeMesh m2 = build_volume_mesh(plenum_geometry(2.5e-3, 1.0e-3, 1));
  bool same = (m.grid.nr == m2.grid.nr) && (m.grid.nz == m2.grid.nz);
  if (same)
    for (std::size_t k = 0; k < m.grid.nodes.size(); ++k)
      same = same && m.grid.nodes[k].r == m2.grid.nodes[k].r &&
             m.grid.nodes[k].z == m2.grid.nodes[k].z;
  expect("zweiter Aufbau ist bitgenau identisch", same);

  // A partly filled plenum has less liquid and a dry cavity floor, and the mesh
  // has to represent the liquid level exactly too.
  DielectricDeviceParameters half = plenum_geometry(2.5e-3, 1.0e-3, 1);
  half.plenum_fill_fraction = 0.4;
  const DeviceVolumeMesh mh = build_volume_mesh(half);
  expect("teilgefuelltes Plenum: alle Netzpruefungen bestanden", mh.validate().all_passed());
  check("Fluessigkeitsvolumen bei Fuellgrad 0.4",
        mh.revolved_volume_of(Region::Liquid),
        pi * a * a * (-mh.z_base) + pi * ac * ac * (mh.z_base - mh.z_roof) +
            pi * mh.r_plenum * mh.r_plenum * 0.4 * mh.p.plenum_depth,
        1.0e-11);
}

// ==========================================================================
// 2 and 3.  Resizing the reservoir changes nothing in front of the base body
// ==========================================================================
//
// This is the check the whole phase turns on.  The comparison is BITWISE on the
// node coordinates and EXACT on the material of every cell, over every row from
// the rear face of the base body forward -- taper, bore, base body, the whole
// gap and the extractor.  Anything weaker would leave room for the mesh to move
// underneath a "convergence" study, which is the mistake being repaired.
static void test_front_geometry_is_untouched() {
  std::printf("\n=== 2. Die Plenumgroesse laesst die vordere Geometrie unberuehrt ===\n");
  const std::vector<std::pair<Real, Real>> sizes{
      {2.5e-3, 1.0e-3}, {3.5e-3, 1.5e-3}, {5.0e-3, 2.5e-3}, {7.0e-3, 4.0e-3}};
  std::vector<DeviceVolumeMesh> m;
  for (const auto& s : sizes) m.push_back(build_volume_mesh(plenum_geometry(s.first, s.second, 1)));

  const DeviceVolumeMesh& ref = m.front();
  bool nodes_identical = true, regions_identical = true, radii_identical = true;
  Index rows = 0, cells = 0;
  for (std::size_t k = 1; k < m.size(); ++k) {
    const DeviceVolumeMesh& q = m[k];
    // The radial node list is a tensor factor: every radius out to the
    // extractor rim must be the same number, not merely the same count.
    if (q.i_ext_outer != ref.i_ext_outer) {
      radii_identical = false;
      continue;
    }
    for (Index i = 0; i <= ref.i_ext_outer; ++i)
      radii_identical = radii_identical && (q.r_ref[static_cast<std::size_t>(i)] ==
                                            ref.r_ref[static_cast<std::size_t>(i)]);
    // Rows from the rear face of the base body forward.
    if (ref.grid.nz - ref.j_base != q.grid.nz - q.j_base) {
      nodes_identical = false;
      continue;
    }
    for (Index j = ref.j_base, jq = q.j_base; j < ref.grid.nz; ++j, ++jq) {
      nodes_identical = nodes_identical && (ref.grid.z_of_row(j) == q.grid.z_of_row(jq));
      for (Index i = 0; i <= ref.i_ext_outer; ++i)
        nodes_identical = nodes_identical && (ref.grid.at(i, j).r == q.grid.at(i, jq).r);
      if (k == 1) ++rows;
      if (j + 1 >= ref.grid.nz) continue;
      for (Index i = 0; i + 1 <= ref.i_ext_outer; ++i) {
        const Region a = ref.cell_region[static_cast<std::size_t>(ref.grid.cell(i, j))];
        const Region b = q.cell_region[static_cast<std::size_t>(q.grid.cell(i, jq))];
        regions_identical = regions_identical && (a == b);
        if (k == 1) ++cells;
      }
    }
  }
  std::printf("  verglichen: %lld Zeilen ab der Grundplatte, %lld Zellen je Variante, "
              "%zu Varianten\n",
              static_cast<long long>(rows), static_cast<long long>(cells), sizes.size() - 1);
  expect("die Radien bis zum Extraktoraussenrand sind bitgenau gleich", radii_identical);
  expect("jeder Nahfeldknoten vor der Grundplatte ist bitgenau gleich", nodes_identical);
  expect("jede Materialzuordnung vor der Grundplatte ist gleich", regions_identical);

  // ... and the reservoir really did change, or the test above proves nothing.
  bool differs = false;
  for (std::size_t k = 1; k < m.size(); ++k)
    differs = differs || (m[k].revolved_volume_of(Region::Liquid) !=
                          ref.revolved_volume_of(Region::Liquid));
  expect("und der Fluessigkeitsvorrat unterscheidet sich tatsaechlich", differs);
}

// ==========================================================================
// 4.  Every liquid region sits at exactly V_emitter, as ONE conductor
// ==========================================================================
static void test_all_liquid_is_one_equipotential() {
  std::printf("\n=== 3. Saemtliche Fluessigkeitsbereiche liegen auf V_emitter ===\n");
  const MaterialLibrary lib;
  const DielectricSetup s = plenum_setup(lib, 2.5e-3, 1.0e-3, 1);
  const DielectricSolution sol = solve_dielectric(s);
  const DeviceVolumeMesh& m = sol.mesh;
  const QuadMesh& g = m.grid;

  // Every node of every liquid cell -- interior and closure alike.
  Index liquid_nodes = 0, wrong_role = 0;
  Real worst = 0.0;
  std::vector<char> seen(static_cast<std::size_t>(g.n_nodes()), 0);
  for (Index j = 0; j + 1 < g.nz; ++j)
    for (Index i = 0; i + 1 < g.nr; ++i) {
      if (m.cell_region[static_cast<std::size_t>(g.cell(i, j))] != Region::Liquid) continue;
      for (Index dj = 0; dj < 2; ++dj)
        for (Index di = 0; di < 2; ++di) {
          const Index n = g.node(i + di, j + dj);
          if (seen[static_cast<std::size_t>(n)]) continue;
          seen[static_cast<std::size_t>(n)] = 1;
          ++liquid_nodes;
          if (sol.role[static_cast<std::size_t>(n)] != NodeRole::LiquidConductor) ++wrong_role;
          worst = std::max(worst, std::abs(sol.fem.phi[static_cast<std::size_t>(n)] -
                                           s.V_emitter));
        }
    }
  std::printf("  %lld Fluessigkeitsknoten, groesste Abweichung von V_emitter: %.3e V\n",
              static_cast<long long>(liquid_nodes), worst);
  expect("jeder Fluessigkeitsknoten traegt die Rolle liquid_conductor", wrong_role == 0);
  expect("und liegt exakt auf V_emitter", worst == 0.0);

  // Spot checks by name, at three places that are far apart in the geometry.
  auto at = [&](Index i, Index j) { return sol.role[static_cast<std::size_t>(g.node(i, j))]; };
  expect("Bohrungswand", at(m.i_bore, (m.j_base + m.j_tip) / 2) == NodeRole::LiquidConductor);
  expect("Wand des Zulaufkanals",
         at(m.i_channel, (m.j_roof + m.j_base) / 2) == NodeRole::LiquidConductor);
  expect("Plenumwand",
         at(m.i_plenum, (m.j_cav_bottom + m.j_roof) / 2) == NodeRole::LiquidConductor);
  expect("Plenumdecke", at((m.i_channel + m.i_plenum) / 2, m.j_roof) ==
                            NodeRole::LiquidConductor);
  expect("Plenumboden", at(m.i_plenum / 2, m.j_cav_bottom) == NodeRole::LiquidConductor);
  expect("es gibt keinen Zulauf-Schnittquerschnitt mehr", sol.audit.n_feed == 0);
  expect("ein wesentlicher Teil der Fluessigkeit liegt hinter der Grundplatte",
         sol.audit.n_reservoir_liquid_surface > 0);
}

// ==========================================================================
// 5.  No polymer and no PEEK node is an electrode
// ==========================================================================
static void test_no_dielectric_is_an_electrode() {
  std::printf("\n=== 4. Kein Polymer- oder PEEK-Knoten ist eine Elektrode ===\n");
  const MaterialLibrary lib;
  const DielectricSetup s = plenum_setup(lib, 3.5e-3, 1.5e-3, 1);
  const DeviceVolumeMesh m = build_volume_mesh(s.geometry);
  const std::vector<NodeRole> role = node_roles(m, s);
  const BoundaryAudit a = audit_boundaries(m, role, s);
  a.print(stdout);
  expect("Audit ohne Verletzung", a.ok());
  expect("strukturell: kein festgehaltener Knoten ohne Fluessigkeitskontakt",
         a.n_polymer_dirichlet == 0);
  expect("benannte Polymerflaechen tragen keine Dirichlet-Bedingung",
         a.n_named_surface_dirichlet == 0);

  auto at = [&](Index i, Index j) { return role[static_cast<std::size_t>(m.grid.node(i, j))]; };
  expect("Oberseite des PEEK-Koerpers ist frei",
         at((m.i_land + m.i_plenum_outer) / 2, m.j_base) == NodeRole::Free);
  expect("Aussenmantel des PEEK-Koerpers ist frei",
         at(m.i_plenum_outer, (m.j_block_bottom + m.j_roof) / 2) == NodeRole::Free);
  expect("Unterseite des PEEK-Koerpers ist frei",
         at(m.i_plenum / 2, m.j_block_bottom) == NodeRole::Free);
  expect("Inneres der PEEK-Wand ist frei",
         at((m.i_plenum + m.i_plenum_outer) / 2, (m.j_cav_bottom + m.j_roof) / 2) ==
             NodeRole::Free);
  expect("Rueckflaeche des gedruckten Emitters ist frei",
         at((m.i_bore + m.i_land) / 2, m.j_base) == NodeRole::Free);

  // A single node of the PEEK outer skin declared an electrode must be caught.
  // Without this the audit would only be testing that node_roles agrees with
  // itself -- and a conducting holder is precisely what must not creep in.
  std::vector<NodeRole> broken = role;
  broken[static_cast<std::size_t>(
      m.grid.node(m.i_plenum_outer, (m.j_block_bottom + m.j_roof) / 2))] =
      NodeRole::LiquidConductor;
  const BoundaryAudit bad = audit_boundaries(m, broken, s);
  expect("ein als Elektrode markierter PEEK-Knoten wird gefunden",
         !bad.ok() && bad.n_polymer_dirichlet == 1);

  // The same for the underside, which no named surface of the emitter covers.
  std::vector<NodeRole> broken2 = role;
  broken2[static_cast<std::size_t>(m.grid.node(m.i_plenum / 2, m.j_block_bottom))] =
      NodeRole::LiquidConductor;
  expect("eine als Elektrode markierte PEEK-Unterseite wird gefunden",
         audit_boundaries(m, broken2, s).n_polymer_dirichlet == 1);

  // A partly filled plenum exposes the cavity floor to vacuum.  It is then a
  // dielectric surface and must lose its Dirichlet condition, not keep it.
  DielectricSetup h = s;
  h.geometry.plenum_fill_fraction = 0.4;
  const DeviceVolumeMesh mh = build_volume_mesh(h.geometry);
  const std::vector<NodeRole> rh = node_roles(mh, h);
  const BoundaryAudit ah = audit_boundaries(mh, rh, h);
  expect("teilgefuellt: Audit ohne Verletzung", ah.ok());
  expect("teilgefuellt: der trockene Plenumboden ist frei",
         rh[static_cast<std::size_t>(mh.grid.node(mh.i_plenum / 2, mh.j_cav_bottom))] ==
             NodeRole::Free);
  expect("teilgefuellt: der Fluessigkeitsspiegel traegt Emitterpotential",
         rh[static_cast<std::size_t>(mh.grid.node(mh.i_plenum / 2, mh.j_fill))] ==
             NodeRole::LiquidConductor);
}

// ==========================================================================
// 6.  The conducting back disc stays refused
// ==========================================================================
static void test_conducting_closure_stays_refused() {
  std::printf("\n=== 5. Die leitfaehige Abschlussscheibe bleibt abgelehnt ===\n");
  DielectricDeviceParameters bad = plenum_geometry(2.5e-3, 1.0e-3, 0);
  bad.device.emitter_back_length = 8.0e-4;
  bool threw = false;
  std::string what;
  try {
    (void)build_volume_mesh(bad);
  } catch (const std::exception& e) {
    threw = true;
    what = e.what();
  }
  expect("emitter_back_length wird auch mit Plenum abgelehnt", threw);
  expect("und die Ablehnung nennt den Grund",
         what.find("Abschlussscheibe") != std::string::npos);

  // The plenum must not lie inside the near field: a radial level there would
  // move near-field nodes when it is resized, which is what test 2 forbids.
  DielectricDeviceParameters inside = plenum_geometry(1.0e-3, 1.0e-3, 0);
  bool threw2 = false;
  try {
    (void)build_volume_mesh(inside);
  } catch (const std::exception&) {
    threw2 = true;
  }
  expect("ein Plenum innerhalb des Extraktoraussenradius wird abgelehnt", threw2);

  // And there is no way to ask for a conducting reservoir: the region carries a
  // permittivity, so a material without a number aborts instead of guessing.
  MaterialLibrary lib;
  bool threw3 = false;
  try {
    (void)lib.get("peek").permittivity_or_throw();
  } catch (const std::exception&) {
    threw3 = true;
  }
  expect("PEEK ohne belegten Wert wird abgelehnt statt geraten", threw3);
}

// ==========================================================================
int main() {
  std::printf("P2c: achsensymmetrisches Plenum als Vorrat-Ersatzmodell\n");
  test_plenum_mesh();
  test_front_geometry_is_untouched();
  test_all_liquid_is_one_equipotential();
  test_no_dielectric_is_an_electrode();
  test_conducting_closure_stays_refused();
  std::printf("\n%s (%d Fehler)\n",
              failures == 0 ? "ALLE TESTS BESTANDEN" : "TESTS FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
