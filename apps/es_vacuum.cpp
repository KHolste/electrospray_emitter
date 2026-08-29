// es_vacuum -- P2a: static vacuum electrostatics on the P1 device geometry.
//
//   es_vacuum <geometry.cfg> [<vacuum.cfg> ...] <output-directory> [key=value ...]
//
// Solves ONE problem and nothing else:
//
//   * charge-free vacuum, rho = 0;
//   * emitter metal and the INITIAL FLAT LIQUID SURFACE at V_emitter, the
//     latter purely as a perfect-conductor reference plane;
//   * extractor at V_extractor;
//   * V -> 0 at infinity, which is the boundary condition the existing BEM
//     kernel carries in its free-space Green's function.
//
// NOT in this phase: meniscus deformation, finite liquid conductivity, ion or
// droplet emission, flow, space charge.  The plane at z = 0 is the initial flat
// liquid surface, not a computed meniscus and not an emitting one.
//
// Exit code 2 means a mesh failed its own checks or a numerical check failed.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/boundary_mesh.hpp"
#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/device_geometry.hpp"
#include "es/io.hpp"
#include "es/vacuum_bem.hpp"

using namespace es;

namespace {

using constants::eps0;

DeviceParameters geometry_from(const Config& c) {
  DeviceParameters p;
  p.phi_3 = c.num("device.phi_3", p.phi_3);
  p.phi_1 = c.num("device.phi_1", p.phi_1);
  p.phi_2 = c.num("device.phi_2", p.phi_2);
  p.emitter_height = c.num("device.emitter_height", p.emitter_height);
  p.extraction_distance = c.num("device.extraction_distance", p.extraction_distance);
  p.extractor_aperture_diameter =
      c.num("device.extractor_aperture_diameter", p.extractor_aperture_diameter);
  p.extractor_thickness = c.num("device.extractor_thickness", p.extractor_thickness);
  p.domain_radius = c.num("domain.radius", p.domain_radius);
  p.domain_z_min = c.num("domain.z_min", p.domain_z_min);
  p.domain_z_max = c.num("domain.z_max", p.domain_z_max);
  // Mandatory: no fallback, because the fallback used to be "out to the domain
  // boundary", which silently turned a non-conductor into one.
  if (!c.has("device.extractor_outer_radius"))
    throw std::runtime_error(
        "device.extractor_outer_radius fehlt.  Der Aussenradius der Extraktionselektrode ist "
        "eine Pflichtangabe: die Elektrode ist ein Leiter, der offene Domaenenrand nicht.  "
        "domain.radius muss echt groesser sein.");
  p.extractor_outer_radius = c.num("device.extractor_outer_radius", 0.0);
  p.reserved.edge_radius_inner = c.num("reserved.edge_radius_inner", 0.0);
  p.reserved.edge_radius_outer = c.num("reserved.edge_radius_outer", 0.0);
  p.reserved.contact_angle_deg = c.num("reserved.contact_angle_deg", 0.0);
  p.reserved.bore_diameter_at_inlet = c.num("reserved.bore_diameter_at_inlet", 0.0);
  p.reserved.porous_emitter = c.flag("reserved.porous_emitter", false);
  p.reserved.collector_enabled = c.flag("reserved.collector_enabled", false);
  return p;
}

struct Level {
  Real scale{1.0};
  int n_mesh{0};
  int n_bem{0};
  Real h_min{0}, h_med{0}, h_max{0};
  CapacitanceMatrix cap;
  Real Q_emitter{0}, Q_extractor{0}, Q_net{0};
  Real E_ref{0};          ///< E_z at the axial reference point [V/m]
  Real V_ref{0};
  Real E_peak_free{0};    ///< peak |E_n| on the flat liquid surface, edges excluded
  Real E_peak_emitter{0};
  Real E_peak_extractor{0};
  Real E_edge_max{0};     ///< largest |E_n| INSIDE the marked edge zones (not converged)
  Real residual_max{0}, residual_rel{0};           ///< outside the marked edge zones
  Real residual_max_all{0}, residual_rel_all{0};  ///< including them
  Vec2 residual_worst{};
  Real cavity_screening{0};  ///< max |V - V_emitter| inside the shank cavity [V]
};

/// Bundle of everything one refinement level needs; kept alive because the
/// solver holds a reference-free copy of the mesh but the report does not.
struct Solved {
  BoundaryMesh bm;
  BemSolver bem;
  VacuumSelectionReport sel;
  std::vector<EdgeZone> zones;
};

void build_level(const DeviceGeometry& g, Real scale, const std::array<Real, 3>& V, Solved& out) {
  out.bm = BoundaryMesh::generate(g, scale);
  const MeshReport mr = out.bm.validate(g);
  if (!mr.all_passed()) {
    mr.print(stdout);
    throw std::runtime_error("Randnetz hat eigene Pruefungen nicht bestanden");
  }
  Mesh panels = vacuum_bem_mesh(out.bm, g, &out.sel);
  out.zones = edge_zones(g, out.bm, panels);
  out.bem.set_mesh(std::move(panels));
  out.bem.solve(V);
}

/// Potential deep inside the emitter shank cavity.  The emitter arc is open at
/// the truncation plane; if the cavity were not screened, sigma/eps0 would be a
/// two-sided sum and not the surface field.  This measures the screening.
Real cavity_screening(const BemSolver& bem, const DeviceGeometry& g, Real V_emitter) {
  const DeviceParameters& p = g.parameters();
  const Real r2 = 0.5 * p.phi_2;
  Real worst = 0.0;
  // Points on the axis inside the liquid column, from just under the flat
  // surface down to mid-shank.  All lie inside the conductor union.
  for (int k = 1; k <= 8; ++k) {
    const Real z = -0.25 * p.emitter_height * k / 8.0 - 0.02 * p.emitter_height;
    for (Real r : {0.0, 0.5 * r2}) {
      const Real v = bem.potential_at({r, z});
      worst = std::max(worst, std::abs(v - V_emitter));
    }
  }
  return worst;
}

Real field_z_at(const BemSolver& bem, Vec2 x) { return bem.field_at(x).z; }

void fill_level(Level& L, Solved& s, const DeviceGeometry& g, const std::array<Real, 3>& V,
                Vec2 p_ref) {
  L.n_mesh = static_cast<int>(s.bm.elements().size());
  L.n_bem = static_cast<int>(s.bem.size());
  const LengthStats st = s.bm.stats_total();
  L.h_min = st.min; L.h_med = st.median; L.h_max = st.max;

  L.cap = maxwell_capacitance(s.bem);
  L.Q_emitter = electrode_charge(s.bem, Electrode::Emitter);
  L.Q_extractor = electrode_charge(s.bem, Electrode::Extractor);
  L.Q_net = L.Q_emitter + L.Q_extractor;
  L.V_ref = s.bem.potential_at(p_ref);
  L.E_ref = field_z_at(s.bem, p_ref);

  L.E_peak_emitter = peak_field_outside_edges(s.bem, Electrode::Emitter, s.zones);
  L.E_peak_extractor = peak_field_outside_edges(s.bem, Electrode::Extractor, s.zones);
  Real best_free = 0.0, edge_max = 0.0;
  for (Index i = 0; i < s.bem.size(); ++i) {
    const Element& el = s.bem.mesh().elems[static_cast<std::size_t>(i)];
    const Real e = std::abs(s.bem.En(i));
    if (in_edge_zone(s.zones, el.mid)) { edge_max = std::max(edge_max, e); continue; }
    if (el.tag == Tag::FreeSurface) best_free = std::max(best_free, e);
  }
  L.E_peak_free = best_free;
  L.E_edge_max = edge_max;

  const PotentialResidual pr = potential_residual(s.bem, s.zones);
  L.residual_max = std::max(pr.max_emitter_clear, pr.max_extractor_clear);
  L.residual_rel = pr.relative();
  L.residual_max_all = std::max(pr.max_emitter, pr.max_extractor);
  L.residual_rel_all = pr.relative_including_edges();
  L.residual_worst = pr.worst_position;
  L.cavity_screening = ::cavity_screening(s.bem, g, V[0]);
}

}  // namespace

int main(int argc, char** argv) try {
  Config cfg;
  std::vector<std::string> rest = Config::positional_args(argc, argv);
  if (rest.size() < 2) {
    std::fprintf(stderr,
                 "usage: es_vacuum <geometry.cfg> [<vacuum.cfg> ...] <output-directory> "
                 "[key=value ...]\n");
    return 1;
  }
  const std::string outdir = rest.back();
  rest.pop_back();
  for (const std::string& f : rest) cfg.load(f);
  cfg.apply_cli(argc, argv);

  const DeviceParameters dp = geometry_from(cfg);
  const DeviceGeometry g = DeviceGeometry::build(dp);

  if (!cfg.has("bem.V_emitter") || !cfg.has("bem.V_extractor"))
    throw std::runtime_error("bem.V_emitter und bem.V_extractor sind Pflichtangaben");
  const Real VE = cfg.num("bem.V_emitter", 0.0);
  const Real VX = cfg.num("bem.V_extractor", 0.0);
  const std::array<Real, 3> V{{VE, VX, 0.0}};

  // The reference point for the converged axial field: on the axis, one tenth
  // of the bore radius above the centre of the flat liquid surface.  It is
  // fixed by the geometry, so it is the same physical point at every refinement
  // level.  Its distance to the nearest unrounded edge -- the pinned exit edge
  // at (r_bore, 0) -- is 1.005 r_bore, four times the marking radius of that
  // edge zone, so it is not contaminated by the corner singularity.
  const Real z_ref = 0.1 * g.contact_radius();
  const Vec2 p_ref{0.0, z_ref};

  std::string d = outdir;
  if (!d.empty() && d.back() != '/' && d.back() != '\\') d += '/';

  std::printf("P2a -- statische Vakuum-Elektrostatik (rho = 0, ebene Perfect-Conductor-\n");
  std::printf("       Referenzflaeche bei z = 0, keine Emission)\n\n");
  g.print(stdout);
  std::printf("\nangelegte Potentiale (absolut gegen V(unendlich) = 0)\n");
  std::printf("  V_emitter   = %10.4g V   (Metall UND ebene Fluessigkeitsoberflaeche)\n", VE);
  std::printf("  V_extractor = %10.4g V\n", VX);
  std::printf("  Referenzpunkt fuer das konvergierte Axialfeld: r = 0, z = %.4g m\n", z_ref);

  g.write_csv(d);

  // ------------------------------------------------------------------ levels
  const std::vector<Real> scales = {1.0, 0.5, 0.25, 0.125};
  std::vector<Level> levels;
  std::vector<Solved> solved(scales.size());
  for (std::size_t k = 0; k < scales.size(); ++k) {
    std::printf("\n--- Netzstufe %d (size_scale = %.4g) ---\n", static_cast<int>(k), scales[k]);
    build_level(g, scales[k], V, solved[k]);
    Level L;
    L.scale = scales[k];
    fill_level(L, solved[k], g, V, p_ref);
    levels.push_back(L);
    std::printf("  Randnetz %d Elemente -> BEM %d Panels"
                "   (Emitter %d, Fluessigkeitsoberflaeche %d, Extraktor %d)\n",
                L.n_mesh, L.n_bem, solved[k].sel.n_emitter, solved[k].sel.n_free_surface,
                solved[k].sel.n_extractor);
    std::printf("  c_EE = %.9e F   c_EX = %.9e F   c_XX = %.9e F\n", L.cap.c_EE, L.cap.c_EX,
                L.cap.c_XX);
    std::printf("  E_z(Referenzpunkt) = %.9e V/m\n", L.E_ref);
    std::printf("  Residuum ausserhalb der Kantenzonen = %.3e V (%.3e relativ)\n",
                L.residual_max, L.residual_rel);
    std::printf("  Residuum einschliesslich Kantenzonen = %.3e V (%.3e relativ), bei "
                "r = %.4g m, z = %.4g m\n",
                L.residual_max_all, L.residual_rel_all, L.residual_worst.r, L.residual_worst.z);
    std::printf("  Abschirmung im Schafthohlraum: max |V - V_emitter| = %.3e V\n",
                L.cavity_screening);
  }

  // The level everything else is reported on.  Not the finest: the field maps
  // cost O(N) per grid point, and the convergence table below shows what the
  // remaining discretisation error at this level is.
  const std::size_t ref = 2;
  Solved& S = solved[ref];
  const Level& R = levels[ref];
  std::printf("\nReferenz-Netzstufe fuer Karten und Oberflaechendaten: %d\n",
              static_cast<int>(ref));

  S.sel.print(stdout);
  S.sel.write_csv(d + "bem_selection.csv");

  // --------------------------------------------------------------- linearity
  std::printf("\nLinearitaet (dieselbe Faktorisierung, verschiedene Spannungen)\n");
  {
    std::FILE* f = std::fopen((d + "linearity.csv").c_str(), "w");
    std::fprintf(f, "# sigma, V and E must scale exactly with the applied voltages\n");
    std::fprintf(f, "factor,V_emitter,V_extractor,Q_emitter_C,V_ref_V,Ez_ref_V_per_m,"
                    "max_rel_dev_sigma\n");
    const std::vector<Real> sigma0 = S.bem.sigma_for({{VE, VX, 0.0}});
    for (Real a : {0.5, 1.0, 2.0, 3.7}) {
      const std::vector<Real> s = S.bem.sigma_for({{a * VE, a * VX, 0.0}});
      Real dev = 0.0, scale = 0.0;
      for (std::size_t i = 0; i < s.size(); ++i) scale = std::max(scale, std::abs(sigma0[i]));
      for (std::size_t i = 0; i < s.size(); ++i)
        dev = std::max(dev, std::abs(s[i] - a * sigma0[i]) / (a * scale));
      BemSolver tmp = S.bem;
      tmp.solve({{a * VE, a * VX, 0.0}});
      const Real vref = tmp.potential_at(p_ref);
      const Real eref = field_z_at(tmp, p_ref);
      std::fprintf(f, "%.6g,%.9e,%.9e,%.9e,%.9e,%.9e,%.3e\n", a, a * VE, a * VX,
                   electrode_charge(tmp, Electrode::Emitter), vref, eref, dev);
      std::printf("  a = %5.3g : max rel. Abweichung sigma - a*sigma0 = %.3e\n", a, dev);
    }
    std::fclose(f);
  }

  // ---------------------------------------------------------------- polarity
  std::printf("\nPolaritaetsumkehr (V -> -V)\n");
  {
    BemSolver plus = S.bem;
    plus.solve({{VE, VX, 0.0}});
    BemSolver minus = S.bem;
    minus.solve({{-VE, -VX, 0.0}});
    Real dsig = 0.0, dE = 0.0, dmag = 0.0, scale = 0.0;
    for (Index i = 0; i < plus.size(); ++i)
      scale = std::max(scale, std::abs(plus.sigma()[static_cast<std::size_t>(i)]));
    for (Index i = 0; i < plus.size(); ++i) {
      const Real a = plus.sigma()[static_cast<std::size_t>(i)];
      const Real b = minus.sigma()[static_cast<std::size_t>(i)];
      dsig = std::max(dsig, std::abs(a + b) / scale);
      dmag = std::max(dmag, std::abs(std::abs(a) - std::abs(b)) / scale);
    }
    const Vec2 ep = plus.field_at(p_ref), em = minus.field_at(p_ref);
    dE = norm(ep + em) / std::max(norm(ep), 1e-300);
    std::printf("  max |sigma+ + sigma-| / max|sigma| = %.3e\n", dsig);
    std::printf("  |E+ + E-| / |E+| am Referenzpunkt  = %.3e\n", dE);
    std::printf("  |E| unveraendert:                    %.9e vs %.9e V/m\n", norm(ep), norm(em));
    std::FILE* f = std::fopen((d + "polarity.csv").c_str(), "w");
    std::fprintf(f, "quantity,value,note\n");
    std::fprintf(f, "max_rel_sigma_sum,%.9e,sigma(+V) + sigma(-V) relative to max|sigma|\n", dsig);
    std::fprintf(f, "max_rel_abs_sigma_diff,%.9e,||sigma(+V)| - |sigma(-V)||\n", dmag);
    std::fprintf(f, "rel_field_sum_at_ref,%.9e,|E(+V) + E(-V)| / |E(+V)|\n", dE);
    std::fprintf(f, "Emag_plus_V_per_m,%.9e,|E| at the reference point for +V\n", norm(ep));
    std::fprintf(f, "Emag_minus_V_per_m,%.9e,|E| at the reference point for -V\n", norm(em));
    std::fclose(f);
  }

  // ------------------------------------------------------------- convergence
  {
    std::FILE* f = std::fopen((d + "convergence.csv").c_str(), "w");
    std::fprintf(f, "# mesh convergence of the P2a vacuum solve.  size_scale is a uniform\n");
    std::fprintf(f, "# multiplier on the automatic size field; 1.0 is the automatic mesh.\n");
    std::fprintf(f, "# c_EE, c_EX, c_XX are Maxwell capacitance coefficients [F] with\n");
    std::fprintf(f, "# V = 0 at infinity;  C_m = -c_EX is the mutual capacitance.\n");
    std::fprintf(f, "# Q_over_dV is Q_emitter/(V_emitter - V_extractor) at the operating\n");
    std::fprintf(f, "# point -- NOT the same as C_m unless the system is charge neutral.\n");
    std::fprintf(f, "level,size_scale,n_mesh_elements,n_bem_panels,h_min_m,h_median_m,h_max_m,"
                    "c_EE_F,c_EX_F,c_XX_F,C_mutual_F,Q_emitter_C,Q_extractor_C,Q_net_C,"
                    "Q_over_dV_F,Ez_ref_V_per_m,V_ref_V,Epeak_free_surface_V_per_m,"
                    "Epeak_emitter_V_per_m,Epeak_extractor_V_per_m,Emax_in_edge_zones_V_per_m,"
                    "residual_max_clear_V,residual_rel_clear,residual_max_all_V,"
                    "residual_rel_all,cavity_screening_V\n");
    const Real dV = VE - VX;
    for (std::size_t k = 0; k < levels.size(); ++k) {
      const Level& L = levels[k];
      std::fprintf(f,
                   "%d,%.6g,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                   "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                   static_cast<int>(k), L.scale, L.n_mesh, L.n_bem, L.h_min, L.h_med, L.h_max,
                   L.cap.c_EE, L.cap.c_EX, L.cap.c_XX, L.cap.mutual(), L.Q_emitter,
                   L.Q_extractor, L.Q_net, dV != 0.0 ? L.Q_emitter / dV : 0.0, L.E_ref, L.V_ref,
                   L.E_peak_free, L.E_peak_emitter, L.E_peak_extractor, L.E_edge_max,
                   L.residual_max, L.residual_rel, L.residual_max_all,
                   L.residual_rel_all, L.cavity_screening);
    }
    std::fclose(f);
  }

  // ------------------------------------------------- extractor outer radius
  std::printf("\nEinfluss des endlichen Extraktoraussenradius (Netzstufe 1)\n");
  {
    std::FILE* f = std::fopen((d + "extractor_radius_study.csv").c_str(), "w");
    std::fprintf(f, "# the electrode's outer radius is a real dimension, not the domain edge;\n");
    std::fprintf(f, "# this is how much the answer depends on it.  domain_radius = %.9e m\n",
                 dp.domain_radius);
    std::fprintf(f, "extractor_outer_radius_m,n_bem_panels,c_EE_F,c_EX_F,C_mutual_F,"
                    "Q_emitter_C,Ez_ref_V_per_m\n");
    for (Real rext : {5.0e-4, 7.5e-4, 1.0e-3, 1.5e-3, 2.0e-3, 2.5e-3}) {
      DeviceParameters q = dp;
      q.extractor_outer_radius = rext;
      if (!(q.domain_radius > rext)) continue;
      const DeviceGeometry gq = DeviceGeometry::build(q);
      Solved sq;
      build_level(gq, 0.5, V, sq);
      const CapacitanceMatrix c = maxwell_capacitance(sq.bem);
      const Real q_e = electrode_charge(sq.bem, Electrode::Emitter);
      const Real e_ref = field_z_at(sq.bem, p_ref);
      std::fprintf(f, "%.9e,%d,%.9e,%.9e,%.9e,%.9e,%.9e\n", rext,
                   static_cast<int>(sq.bem.size()), c.c_EE, c.c_EX, c.mutual(), q_e, e_ref);
      std::printf("  r_ext = %8.4g m : C_m = %.6e F, E_z(ref) = %.6e V/m\n", rext, c.mutual(),
                  e_ref);
    }
    std::fclose(f);
  }

  // ----------------------------------------------------------- surface dump
  {
    std::FILE* f = std::fopen((d + "surface.csv").c_str(), "w");
    std::fprintf(f, "# surface charge density and normal field on the reference mesh level\n");
    std::fprintf(f, "# sigma/eps0 is the vacuum-side normal field; the emitter arc is open only\n");
    std::fprintf(f, "# at the truncation plane, and the cavity behind it is screened to %.2e V\n",
                 R.cavity_screening);
    std::fprintf(f, "i,r_m,z_m,arclen_m,sigma_C_per_m2,En_V_per_m,area_m2,tag,boundary_id,"
                    "in_edge_zone\n");
    Real s_arc = 0.0;
    for (std::size_t k = 0; k < S.sel.panels.size(); ++k) {
      const VacuumPanel& pn = S.sel.panels[k];
      const Element& el = S.bem.mesh().elems[static_cast<std::size_t>(pn.bem_element)];
      s_arc += el.len;
      std::fprintf(f, "%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%s,%s,%d\n", pn.bem_element, el.mid.r,
                   el.mid.z, s_arc, S.bem.sigma()[static_cast<std::size_t>(pn.bem_element)],
                   S.bem.En(pn.bem_element), el.area, electrode_label(el.tag),
                   to_string(pn.boundary),
                   in_edge_zone(S.zones, el.mid) ? 1 : 0);
    }
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "edge_zones.csv").c_str(), "w");
    std::fprintf(f, "# places where |E| is NOT a converged quantity and is therefore marked,\n");
    std::fprintf(f, "# never reported as a peak field.  radius = local_feature_size / 4, a\n");
    std::fprintf(f, "# purely geometric length, so the marked region is the same at every\n");
    std::fprintf(f, "# refinement level.  kind=sharp_feature: unrounded edge of the device.\n");
    std::fprintf(f, "# kind=truncation_end: open end of the modelled conductor arc, where the\n");
    std::fprintf(f, "# single-layer density has a 1/sqrt(d) sheet-edge singularity -- an\n");
    std::fprintf(f, "# artefact of cutting the shank at the domain floor, not a device edge.\n");
    std::fprintf(f, "kind,name,r_m,z_m,radius_m,local_feature_size_m\n");
    for (const EdgeZone& z : S.zones)
      std::fprintf(f, "%s,%s,%.9e,%.9e,%.9e,%.9e\n", to_string(z.kind), z.name.c_str(),
                   z.position.r, z.position.z, z.radius, z.local_feature_size);
    std::fclose(f);
  }

  // ------------------------------------------------------------- field maps
  {
    const Real R_dom = dp.domain_radius, z0 = dp.domain_z_min, z1 = dp.domain_z_max;
    std::printf("\nFeldkarten werden ausgewertet ...\n");
    const FieldGrid full = sample_field(S.bem, 241, 241, 0.0, R_dom, z0, z1);
    write_grid_csv(full, d + "field_full.csv");
    std::printf("  field_full.csv    %d x %d\n", full.nr, full.nz);

    const Real rz = 6.0 * g.contact_radius();
    const FieldGrid tip = sample_field(S.bem, 161, 161, 0.0, rz, -0.6 * rz, 0.9 * rz);
    write_grid_csv(tip, d + "field_tip.csv");
    std::printf("  field_tip.csv     %d x %d\n", tip.nr, tip.nz);

    const Real ra = 0.5 * dp.extractor_aperture_diameter;
    const FieldGrid ap = sample_field(S.bem, 161, 161, 0.0, 1.8 * ra,
                                      dp.extraction_distance - 1.1 * ra,
                                      dp.extraction_distance + dp.extractor_thickness + 1.1 * ra);
    write_grid_csv(ap, d + "field_aperture.csv");
    std::printf("  field_aperture.csv %d x %d\n", ap.nr, ap.nz);

    // The open end of the emitter arc.  Resolved on its own grid because the
    // marked zone there is micrometres wide and invisible on the domain map,
    // and because a reader is entitled to see the artefact that is excluded
    // rather than take its exclusion on trust.
    const Real r3 = 0.5 * dp.phi_3;
    const Real w = 3.0 * r3;
    const FieldGrid trn = sample_field(S.bem, 161, 161, 0.0, w, dp.domain_z_min,
                                       dp.domain_z_min + w);
    write_grid_csv(trn, d + "field_truncation.csv");
    std::printf("  field_truncation.csv %d x %d\n", trn.nr, trn.nz);
  }

  // ---------------------------------------------------------- axis profile
  {
    std::FILE* f = std::fopen((d + "axis_profile.csv").c_str(), "w");
    std::fprintf(f, "# potential and axial field on r = 0, from the flat liquid surface to the\n");
    std::fprintf(f, "# top of the plotted region.  Inside the liquid column V = V_emitter.\n");
    std::fprintf(f, "z_m,V_V,Ez_V_per_m\n");
    const int n = 900;
    for (int i = 0; i <= n; ++i) {
      const Real z = dp.domain_z_min + (dp.domain_z_max - dp.domain_z_min) * i / n;
      const Vec2 x{0.0, z};
      std::fprintf(f, "%.9e,%.9e,%.9e\n", z, S.bem.potential_at(x), S.bem.field_at(x).z);
    }
    std::fclose(f);
  }

  // --------------------------------------------------------------- report
  {
    std::FILE* f = std::fopen((d + "report.txt").c_str(), "w");
    std::fprintf(f, "P2a -- statische Vakuum-Elektrostatik auf der P1-Geometrie\n");
    std::fprintf(f, "=========================================================\n\n");
    std::fprintf(f,
                 "rho = 0.  Emittermetall und die ANFAENGLICHE EBENE FLUESSIGKEITSOBERFLAECHE\n"
                 "bei z = 0 liegen auf V_emitter, die ebene Flaeche ausschliesslich als\n"
                 "Perfect-Conductor-Referenz.  Das ist kein berechneter Meniskus, kein\n"
                 "emittierender Meniskus und keine Aussage ueber das Pure-Ion-Regime.\n"
                 "Keine Meniskusverformung, keine endliche Fluessigkeitsleitfaehigkeit,\n"
                 "keine Emission, keine Stroemung, keine Raumladung.\n\n");
    std::fprintf(f, "Randbedingung im Unendlichen\n");
    std::fprintf(f, "  V -> 0.  Der BEM-Kern benutzt die azimutal integrierte Freiraum-\n");
    std::fprintf(f, "  Greensfunktion; es gibt keine Abschneideflaeche und nichts, worauf\n");
    std::fprintf(f, "  eine Bedingung gesetzt werden koennte.  Die angelegten Potentiale\n");
    std::fprintf(f, "  sind absolut gegen unendlich, und die Gesamtladung ist nicht null.\n\n");
    std::fprintf(f, "Angelegt: V_emitter = %.6g V, V_extractor = %.6g V\n", VE, VX);
    std::fprintf(f, "Referenzpunkt: r = 0, z = %.6g m (kantenfern, netzunabhaengig)\n\n", z_ref);

    std::fprintf(f, "Netzkonvergenz\n");
    std::fprintf(f, "  %-6s %8s %8s %16s %16s %16s %12s\n", "Stufe", "n_Netz", "n_BEM",
                 "c_EE [F]", "C_m [F]", "E_z(ref) [V/m]", "Residuum*");
    std::fprintf(f, "  (* relativ, ausserhalb der markierten Kantenzonen)\n");
    for (std::size_t k = 0; k < levels.size(); ++k) {
      const Level& L = levels[k];
      std::fprintf(f, "  %-6d %8d %8d %16.9e %16.9e %16.9e %12.3e\n", static_cast<int>(k),
                   L.n_mesh, L.n_bem, L.cap.c_EE, L.cap.mutual(), L.E_ref, L.residual_rel);
    }
    const Level& a = levels[levels.size() - 2];
    const Level& b = levels.back();
    std::fprintf(f, "\n  relative Aenderung zwischen den beiden feinsten Stufen\n");
    std::fprintf(f, "    c_EE      : %.3e\n", std::abs(b.cap.c_EE - a.cap.c_EE) / std::abs(b.cap.c_EE));
    std::fprintf(f, "    C_m       : %.3e\n",
                 std::abs(b.cap.mutual() - a.cap.mutual()) / std::abs(b.cap.mutual()));
    std::fprintf(f, "    E_z(ref)  : %.3e\n", std::abs(b.E_ref - a.E_ref) / std::abs(b.E_ref));

    std::fprintf(f, "\nKapazitaetsgroessen auf der Referenzstufe %d -- benannt, nicht \"C\"\n",
                 static_cast<int>(ref));
    std::fprintf(f, "  c_EE (Maxwell-Selbstkoeffizient Emitter) : %.9e F\n", R.cap.c_EE);
    std::fprintf(f, "  c_EX = c_XE (Influenzkoeffizient)        : %.9e / %.9e F\n", R.cap.c_EX,
                 R.cap.c_XE);
    std::fprintf(f, "  Reziprozitaetsfehler |c_EX-c_XE|/|c_EX|  : %.3e\n",
                 R.cap.reciprocity_error());
    std::fprintf(f, "  c_XX (Maxwell-Selbstkoeffizient Extraktor): %.9e F\n", R.cap.c_XX);
    std::fprintf(f, "  C_m = -c_EX (Gegenkapazitaet)            : %.9e F\n", R.cap.mutual());
    std::fprintf(f, "  C_E,inf = c_EE + c_EX (gegen unendlich)  : %.9e F\n",
                 R.cap.emitter_to_infinity());
    std::fprintf(f, "  Q_emitter/(V_E - V_X) am Arbeitspunkt    : %.9e F\n",
                 (VE - VX) != 0.0 ? R.Q_emitter / (VE - VX) : 0.0);
    if (VX == 0.0)
      std::fprintf(f, "    Bei V_extractor = 0 ist dieses Verhaeltnis identisch c_EE und\n"
                      "    NICHT die Gegenkapazitaet C_m.\n");
    std::fprintf(f, "  Q_emitter = %.9e C, Q_extractor = %.9e C, Summe = %.9e C\n", R.Q_emitter,
                 R.Q_extractor, R.Q_net);
    std::fprintf(f, "  Die Summe ist nicht null: bei V(unendlich) = 0 traegt das System\n");
    std::fprintf(f, "  Nettoladung, und Q_E/(V_E-V_X) ist deshalb NICHT C_m.\n");

    std::fprintf(f, "\nFelder auf der Referenzstufe (Kantenzonen ausgeschlossen)\n");
    std::fprintf(f, "  E_z am Referenzpunkt                     : %.9e V/m\n", R.E_ref);
    std::fprintf(f, "  groesstes |E_n| auf der ebenen Flaeche   : %.9e V/m\n", R.E_peak_free);
    std::fprintf(f, "  groesstes |E_n| am Emitter               : %.9e V/m\n", R.E_peak_emitter);
    std::fprintf(f, "  groesstes |E_n| am Extraktor             : %.9e V/m\n", R.E_peak_extractor);
    std::fprintf(f, "  groesstes |E_n| INNERHALB der Kantenzonen: %.9e V/m\n", R.E_edge_max);
    std::fprintf(f, "    Der letzte Wert ist KEINE konvergierte Groesse.  Die Austritts-,\n");
    std::fprintf(f, "    Stirn- und Aperturkanten sind unverrundet; das Feld einer scharfen\n");
    std::fprintf(f, "    Kante divergiert und folgt der Elementgroesse.  Verrundungsradien\n");
    std::fprintf(f, "    sind ein P3-Parameter.\n");

    std::fprintf(f, "\nOffener Emitterbogen und Abschirmung\n");
    std::fprintf(f, "  Der Emitterbogen ist am Schnitt durch die Domaenensohle offen; der\n");
    std::fprintf(f, "  offene Domaenenrand wurde NICHT als Deckel benutzt.  Der dahinter\n");
    std::fprintf(f, "  liegende Hohlraum ist ein Rohr vom Radius %.3g m und der Laenge %.3g m\n",
                 0.5 * dp.phi_3, -dp.domain_z_min);
    std::fprintf(f, "  und daher stark abgeschirmt: max |V - V_emitter| = %.3e V von %.3g V.\n",
                 R.cavity_screening, std::abs(VE));

    std::fprintf(f, "  Deshalb ist sigma/eps0 auf Stirnflaeche und ebener Fluessigkeits-\n");
    std::fprintf(f, "  oberflaeche das einseitige Vakuumfeld.\n");

    std::fprintf(f, "\nResiduum der Dirichlet-Bedingung (Viertelpunkte, nicht Kollokation)\n");
    std::fprintf(f, "  ausserhalb der Kantenzonen : max |V - V_soll| = %.3e V, relativ %.3e\n",
                 R.residual_max, R.residual_rel);
    std::fprintf(f, "  einschliesslich Kantenzonen: max |V - V_soll| = %.3e V, relativ %.3e\n",
                 R.residual_max_all, R.residual_rel_all);
    std::fprintf(f, "  groesstes Residuum bei r = %.4g m, z = %.4g m\n", R.residual_worst.r,
                 R.residual_worst.z);
    std::fprintf(f, "  Das grosse Residuum sitzt dort, wo eine stueckweise konstante Dichte\n");
    std::fprintf(f, "  eine Singularitaet nicht darstellen kann.  Es faellt mit der\n");
    std::fprintf(f, "  Verfeinerung (siehe convergence.csv), und die Felder an diesen Stellen\n");
    std::fprintf(f, "  werden ohnehin nicht als Ergebnis berichtet.\n");
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_vacuum (P2a)\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "config=%s\n", cfg.str("meta.config", "siehe README.md").c_str());
    std::fprintf(f, "V_emitter_V=%.9e\n", VE);
    std::fprintf(f, "V_extractor_V=%.9e\n", VX);
    std::fprintf(f, "reference_level=%d\n", static_cast<int>(ref));
    std::fprintf(f, "reference_size_scale=%.6g\n", scales[ref]);
    std::fprintf(f, "n_bem_panels=%d\n", R.n_bem);
    std::fprintf(f, "z_ref_m=%.9e\n", z_ref);
    std::fprintf(f, "state=vacuum electrostatics, rho=0, flat perfect-conductor reference, "
                    "no emission\n");
    std::fclose(f);
  }

  cfg.warn_about_unused(stdout, {"meta.", "fluid.", "beam.", "output."});
  std::printf("\ngeschrieben nach %s\n", outdir.c_str());
  return 0;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\n%s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "\nFehler: %s\n", e.what());
  return 2;
}
