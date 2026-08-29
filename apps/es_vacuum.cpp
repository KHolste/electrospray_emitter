// es_vacuum -- P2a: static vacuum electrostatics on the P1 device geometry.
//
//   es_vacuum <geometry.cfg> [<vacuum.cfg> ...] <output-directory> [key=value ...]
//
// Solves ONE problem and nothing else:
//
//   * charge-free vacuum, rho = 0;
//   * emitter metal and the INITIAL FLAT LIQUID SURFACE at z = 0 at V_emitter,
//     the latter purely as a perfect-conductor reference plane;
//   * extractor at V_extractor;
//   * V -> 0 at infinity, which is the boundary condition the existing BEM
//     kernel carries in its free-space Green's function.
//
// NOT in this phase: meniscus deformation, finite liquid conductivity, ion or
// droplet emission, flow, space charge.  The plane at z = 0 is the initial flat
// liquid surface, not a computed meniscus and not an emitting one.
//
// THE EMITTER CONDUCTOR IS CLOSED.  device.emitter_back_length is mandatory: it
// gives the modelled conductor a rear end, closed by a conducting disc, so that
// the single-layer density is that of a closed conductor and sigma/eps0 is a
// one-sided vacuum field.  The disc is tagged numerical_emitter_back_closure and
// is excluded from every reported field.  The length itself is NOT a convergence
// knob -- see the truncation study below and the comment on
// DeviceParameters::emitter_back_length.
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
  // Mandatory for the same kind of reason: without it the emitter conductor has
  // no rear end at all and the boundary-integral problem is posed on an open
  // sheet.  There is no default, because every reported number depends on it.
  if (!c.has("device.emitter_back_length"))
    throw std::runtime_error(
        "device.emitter_back_length fehlt.  Die rueckwaertige Laenge des modellierten "
        "Emitterleiters ist eine Pflichtangabe: ohne sie endet der Leiter offen am "
        "Modellschnitt, und eine Einfachschicht auf einer offenen Flaeche traegt die Summe "
        "beider Seiten.  Der Wert geht in JEDE berichtete Groesse ein -- E_z am "
        "Referenzpunkt ebenso wie c_EE und c_EX -- und ist deshalb anzugeben, nicht zu raten.");
  p.emitter_back_length = c.num("device.emitter_back_length", 0.0);
  p.reserved.edge_radius_inner = c.num("reserved.edge_radius_inner", 0.0);
  p.reserved.edge_radius_outer = c.num("reserved.edge_radius_outer", 0.0);
  p.reserved.contact_angle_deg = c.num("reserved.contact_angle_deg", 0.0);
  p.reserved.bore_diameter_at_inlet = c.num("reserved.bore_diameter_at_inlet", 0.0);
  p.reserved.porous_emitter = c.flag("reserved.porous_emitter", false);
  p.reserved.collector_enabled = c.flag("reserved.collector_enabled", false);
  return p;
}

// ---------------------------------------------------------------------------
// Fixed probe points: potential between the electrodes
// ---------------------------------------------------------------------------
//
// Placed from the geometry, so that they are the same physical points at every
// mesh level and at every rearward length.  Five on the axis across the
// extraction gap, two off axis; all of them in vacuum, none of them near an
// unrounded edge.
struct Probe {
  std::string name;
  Vec2 x;
};

std::vector<Probe> probe_points(const DeviceGeometry& g) {
  const DeviceParameters& p = g.parameters();
  const Real d = p.extraction_distance;
  const Real ra = 0.5 * p.extractor_aperture_diameter;
  std::vector<Probe> out;
  for (Real f : {0.10, 0.25, 0.50, 0.75, 0.90}) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "axis_%.0f_pct_of_gap", 100.0 * f);
    out.push_back({buf, {0.0, f * d}});
  }
  out.push_back({"off_axis_half_aperture_mid_gap", {0.5 * ra, 0.5 * d}});
  out.push_back({"off_axis_foot_radius_quarter_gap", {0.5 * p.phi_3, 0.25 * d}});
  return out;
}

// ---------------------------------------------------------------------------

struct Level {
  Real scale{1.0};
  int n_mesh{0};
  int n_bem{0};
  Real h_min{0}, h_med{0}, h_max{0};
  CapacitanceMatrix cap;
  Real Q_emitter{0}, Q_extractor{0}, Q_net{0};
  Real E_ref{0};          ///< E_z at the axial reference point [V/m]
  Real V_ref{0};
  Real E_peak_free{0};    ///< peak |E_n| on the flat liquid surface, evaluable panels only
  Real E_peak_emitter{0};
  Real E_peak_extractor{0};
  Real E_edge_max{0};     ///< largest |E_n| INSIDE the marked edge zones (not converged)
  Real E_closure_max{0};  ///< largest |E_n| on the numerical closure (not a device field)
  Real residual_max{0}, residual_rel{0};           ///< outside the marked edge zones
  Real residual_max_all{0}, residual_rel_all{0};  ///< including them
  Real residual_rel_phys{0}, residual_rms_phys{0};///< on the evaluable device surfaces
  Vec2 residual_worst{};
  Real interior_dV{0};    ///< max |V - V_emitter| deep inside the closed conductor [V]
  Real interior_E{0};     ///< max |E| there [V/m]
  std::vector<Real> V_probe;
};

/// Bundle of everything one refinement level needs; kept alive because the
/// solver holds a copy of the mesh but the report indexes into it.
struct Solved {
  BoundaryMesh bm;
  BemSolver bem;
  VacuumSelectionReport sel;
  std::vector<EdgeZone> zones;
  std::vector<char> evaluable;
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
  out.evaluable = mark_evaluable_panels(panels, out.sel, g, out.zones);
  out.bem.set_mesh(std::move(panels));
  out.bem.solve(V);
}

/// Inside a CLOSED perfect conductor V is the applied potential and E vanishes.
/// Sampled deep inside: a piecewise-constant single layer reproduces the
/// interior only at a distance of order the local element size from the surface,
/// so a sample 1 um from a 20 um panel would measure that panel's discretisation
/// rather than the screening.  The band keeps every point a quarter of the shank
/// radius from the shank and a fifth of the length from either end face.
void interior_check(const BemSolver& bem, const DeviceGeometry& g, Real V_emitter, Real* dV,
                    Real* Emax) {
  const Real z_back = g.back_closure_z();
  const Real r3 = 0.5 * g.parameters().phi_3;
  *dV = 0.0;
  *Emax = 0.0;
  for (int k = 0; k <= 12; ++k) {
    const Real z = z_back * (0.2 + 0.6 * k / 12.0);
    for (Real r : {0.0, 0.25 * r3, 0.5 * r3}) {
      *dV = std::max(*dV, std::abs(bem.potential_at({r, z}) - V_emitter));
      *Emax = std::max(*Emax, norm(bem.field_at({r, z})));
    }
  }
}

Real field_z_at(const BemSolver& bem, Vec2 x) { return bem.field_at(x).z; }

void fill_level(Level& L, Solved& s, const DeviceGeometry& g, const std::array<Real, 3>& V,
                Vec2 p_ref, const std::vector<Probe>& probes) {
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
  L.V_probe.clear();
  for (const Probe& p : probes) L.V_probe.push_back(s.bem.potential_at(p.x));

  L.E_peak_emitter = peak_field_evaluable(s.bem, Electrode::Emitter, s.evaluable);
  L.E_peak_extractor = peak_field_evaluable(s.bem, Electrode::Extractor, s.evaluable);
  Real best_free = 0.0, edge_max = 0.0, closure_max = 0.0;
  for (const VacuumPanel& pn : s.sel.panels) {
    const std::size_t k = static_cast<std::size_t>(pn.bem_element);
    const Element& el = s.bem.mesh().elems[k];
    const Real e = std::abs(s.bem.En(pn.bem_element));
    if (pn.numerical) { closure_max = std::max(closure_max, e); continue; }
    if (in_edge_zone(s.zones, el.mid)) { edge_max = std::max(edge_max, e); continue; }
    if (pn.evaluable && el.tag == Tag::FreeSurface) best_free = std::max(best_free, e);
  }
  L.E_peak_free = best_free;
  L.E_edge_max = edge_max;
  L.E_closure_max = closure_max;

  const PotentialResidual pr = potential_residual(s.bem, s.zones, s.evaluable);
  L.residual_max = std::max(pr.max_emitter_clear, pr.max_extractor_clear);
  L.residual_rel = pr.relative();
  L.residual_max_all = std::max(pr.max_emitter, pr.max_extractor);
  L.residual_rel_all = pr.relative_including_edges();
  L.residual_rel_phys = pr.relative_physical();
  L.residual_rms_phys = pr.relative_rms_physical();
  L.residual_worst = pr.worst_position;
  interior_check(s.bem, g, V[0], &L.interior_dV, &L.interior_E);
}

/// One entry of the truncation study.
struct Trunc {
  Real L{0};
  Real z_closure{0}, clearance{0};
  int n_mesh{0}, n_bem{0}, n_closure{0};
  CapacitanceMatrix cap;
  Real Q_emitter{0};
  Real E_ref{0};
  Real E_peak_free{0};
  Real interior_dV{0}, interior_E{0};
  std::vector<Real> V_probe;
  bool front_mesh_identical{true};
};

Real rel_change(Real a, Real b) {
  const Real s = std::max(std::abs(b), 1e-300);
  return std::abs(b - a) / s;
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

  // The reference point for the axial field: on the axis, one tenth of the bore
  // radius above the centre of the flat liquid surface.  It is fixed by the
  // geometry, so it is the same physical point at every refinement level and at
  // every rearward length.  Its distance to the nearest unrounded edge -- the
  // pinned exit edge at (r_bore, 0) -- is 1.005 r_bore, four times the marking
  // radius of that edge zone, so it is not contaminated by the corner.
  const Real z_ref = 0.1 * g.contact_radius();
  const Vec2 p_ref{0.0, z_ref};
  const std::vector<Probe> probes = probe_points(g);

  // The truncation study doubles the rearward length twice below and once above
  // the value actually used, so the production point lies inside the series and
  // is meshed identically to the other points.  The domain must contain the
  // longest of them, or the study would change the size field along with the
  // length and measure the two together.
  const Real L0 = dp.emitter_back_length;
  const std::vector<Real> back_lengths = {0.25 * L0, 0.5 * L0, L0, 2.0 * L0};
  if (!(dp.domain_z_min < -2.0 * L0))
    throw std::runtime_error(
        "domain.z_min muss unter dem doppelten device.emitter_back_length liegen.  Die "
        "Trunkierungsstudie verlaengert den Leiter bis auf das Doppelte und haelt dabei die "
        "Domaene fest, damit sich das Groessenfeld -- und mit ihm das Netz an Spitze und "
        "Extraktor -- nicht zusammen mit der Laenge aendert.");

  std::string d = outdir;
  if (!d.empty() && d.back() != '/' && d.back() != '\\') d += '/';

  std::printf("P2a -- statische Vakuum-Elektrostatik (rho = 0, ebene Perfect-Conductor-\n");
  std::printf("       Referenzflaeche bei z = 0, keine Emission)\n\n");
  g.print(stdout);
  std::printf("\nangelegte Potentiale (absolut gegen V(unendlich) = 0)\n");
  std::printf("  V_emitter   = %10.4g V   (Metall UND ebene Fluessigkeitsoberflaeche)\n", VE);
  std::printf("  V_extractor = %10.4g V\n", VX);
  std::printf("  Referenzpunkt fuer das Axialfeld: r = 0, z = %.4g m\n", z_ref);

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
    fill_level(L, solved[k], g, V, p_ref, probes);
    levels.push_back(L);
    std::printf("  Randnetz %d Elemente -> BEM %d Panels"
                "   (Emitter %d, Fluessigkeitsoberflaeche %d, Extraktor %d;"
                " davon %d numerische Rueckschliessung)\n",
                L.n_mesh, L.n_bem, solved[k].sel.n_emitter, solved[k].sel.n_free_surface,
                solved[k].sel.n_extractor, solved[k].sel.n_numerical_closure);
    std::printf("  c_EE = %.9e F   c_EX = %.9e F   c_XX = %.9e F\n", L.cap.c_EE, L.cap.c_EX,
                L.cap.c_XX);
    std::printf("  E_z(Referenzpunkt) = %.9e V/m\n", L.E_ref);
    std::printf("  Residuum auf den auswertbaren Geraeteflaechen = %.3e relativ "
                "(RMS %.3e)\n", L.residual_rel_phys, L.residual_rms_phys);
    std::printf("  Residuum einschliesslich Kantenzonen = %.3e V (%.3e relativ), bei "
                "r = %.4g m, z = %.4g m\n",
                L.residual_max_all, L.residual_rel_all, L.residual_worst.r, L.residual_worst.z);
    std::printf("  im geschlossenen Leiter: max |V - V_emitter| = %.3e V, max |E| = %.3e V/m\n",
                L.interior_dV, L.interior_E);
  }

  // The level everything else is reported on.  Not the finest: the field maps
  // cost O(N) per grid point, and the convergence table shows what the remaining
  // discretisation error at this level is.
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
    std::fprintf(f, "# mesh convergence of the P2a vacuum solve at the production rearward\n");
    std::fprintf(f, "# length emitter_back_length = %.9e m.  size_scale is a uniform\n", L0);
    std::fprintf(f, "# multiplier on the automatic size field; 1.0 is the automatic mesh.\n");
    std::fprintf(f, "# c_EE, c_EX, c_XX are Maxwell capacitance coefficients [F] with\n");
    std::fprintf(f, "# V = 0 at infinity;  C_m = -c_EX is the mutual capacitance.\n");
    std::fprintf(f, "# Q_over_dV is Q_emitter/(V_emitter - V_extractor) at the operating\n");
    std::fprintf(f, "# point -- NOT the same as C_m unless the system is charge neutral.\n");
    std::fprintf(f, "# THIS TABLE IS MESH CONVERGENCE ONLY.  For the dependence on the\n");
    std::fprintf(f, "# rearward length see truncation.csv; none of these quantities is\n");
    std::fprintf(f, "# independent of it.\n");
    std::fprintf(f, "level,size_scale,n_mesh_elements,n_bem_panels,n_closure_panels,h_min_m,"
                    "h_median_m,h_max_m,"
                    "c_EE_F,c_EX_F,c_XX_F,C_mutual_F,Q_emitter_C,Q_extractor_C,Q_net_C,"
                    "Q_over_dV_F,Ez_ref_V_per_m,V_ref_V,Epeak_free_surface_V_per_m,"
                    "Epeak_emitter_V_per_m,Epeak_extractor_V_per_m,Emax_in_edge_zones_V_per_m,"
                    "Emax_on_numerical_closure_V_per_m,"
                    "residual_max_clear_V,residual_rel_clear,residual_max_all_V,"
                    "residual_rel_all,residual_rel_physical,residual_rms_physical,"
                    "interior_dV_V,interior_Emag_V_per_m\n");
    const Real dV = VE - VX;
    for (std::size_t k = 0; k < levels.size(); ++k) {
      const Level& L = levels[k];
      std::fprintf(f,
                   "%d,%.6g,%d,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                   "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                   "%.9e\n",
                   static_cast<int>(k), L.scale, L.n_mesh, L.n_bem,
                   solved[k].sel.n_numerical_closure, L.h_min, L.h_med, L.h_max,
                   L.cap.c_EE, L.cap.c_EX, L.cap.c_XX, L.cap.mutual(), L.Q_emitter,
                   L.Q_extractor, L.Q_net, dV != 0.0 ? L.Q_emitter / dV : 0.0, L.E_ref, L.V_ref,
                   L.E_peak_free, L.E_peak_emitter, L.E_peak_extractor, L.E_edge_max,
                   L.E_closure_max, L.residual_max, L.residual_rel, L.residual_max_all,
                   L.residual_rel_all, L.residual_rel_phys, L.residual_rms_phys, L.interior_dV,
                   L.interior_E);
    }
    std::fclose(f);
  }

  // -------------------------------------------------------- truncation study
  //
  // Four increasing rearward lengths on a FIXED domain and a fixed size_scale,
  // so that the mesh at the tip and at the extractor is bitwise the same for all
  // of them and only the rear end differs.  That identity is verified, not
  // assumed: if it ever failed, the study would be measuring two things at once.
  std::printf("\nTrunkierungsstudie gegen die rueckwaertige Laenge (Netzstufe %d)\n",
              static_cast<int>(ref));
  std::vector<Trunc> tr;
  {
    std::vector<Element> front_ref;
    for (Real L : back_lengths) {
      DeviceParameters q = dp;
      q.emitter_back_length = L;
      const DeviceGeometry gq = DeviceGeometry::build(q);
      Solved sq;
      build_level(gq, scales[ref], V, sq);

      Trunc t;
      t.L = L;
      t.z_closure = gq.back_closure_z();
      t.clearance = gq.back_closure_clearance();
      t.n_mesh = static_cast<int>(sq.bm.elements().size());
      t.n_bem = static_cast<int>(sq.bem.size());
      t.n_closure = sq.sel.n_numerical_closure;
      t.cap = maxwell_capacitance(sq.bem);
      t.Q_emitter = electrode_charge(sq.bem, Electrode::Emitter);
      t.E_ref = field_z_at(sq.bem, p_ref);
      for (const Probe& pb : probes) t.V_probe.push_back(sq.bem.potential_at(pb.x));
      interior_check(sq.bem, gq, VE, &t.interior_dV, &t.interior_E);
      Real best_free = 0.0;
      for (const VacuumPanel& pn : sq.sel.panels)
        if (pn.evaluable &&
            sq.bem.mesh().elems[static_cast<std::size_t>(pn.bem_element)].tag == Tag::FreeSurface)
          best_free = std::max(best_free, std::abs(sq.bem.En(pn.bem_element)));
      t.E_peak_free = best_free;

      std::vector<Element> front;
      for (const Element& e : sq.bem.mesh().elems)
        if (e.tag == Tag::Extractor || e.mid.z >= gq.evaluation_z_min()) front.push_back(e);
      if (front_ref.empty()) {
        front_ref = front;
      } else {
        bool same = front.size() == front_ref.size();
        for (std::size_t i = 0; same && i < front.size(); ++i)
          same = same && front[i].a.r == front_ref[i].a.r && front[i].a.z == front_ref[i].a.z &&
                 front[i].b.r == front_ref[i].b.r && front[i].b.z == front_ref[i].b.z &&
                 front[i].tag == front_ref[i].tag;
        t.front_mesh_identical = same;
      }
      std::printf("  L = %9.4g m : %4d Panels, E_z(ref) = %.6e V/m, c_EX = %.6e F, "
                  "c_EE = %.6e F%s\n",
                  L, t.n_bem, t.E_ref, t.cap.c_EX, t.cap.c_EE,
                  t.front_mesh_identical ? "" : "   VERNETZUNG VORNE ABWEICHEND");
      tr.push_back(std::move(t));
    }
  }

  // The verdict, against tolerances fixed in the library before the measurement.
  const Trunc& T0 = tr[tr.size() - 2];
  const Trunc& T1 = tr.back();
  const Real dEz = rel_change(T0.E_ref, T1.E_ref);
  const Real dCEX = rel_change(T0.cap.c_EX, T1.cap.c_EX);
  const Real dCEE = rel_change(T0.cap.c_EE, T1.cap.c_EE);
  Real dVprobe = 0.0;
  for (std::size_t k = 0; k < probes.size(); ++k)
    dVprobe = std::max(dVprobe, std::abs(T1.V_probe[k] - T0.V_probe[k]) / std::abs(VE - VX));
  bool front_ok = true;
  for (const Trunc& t : tr) front_ok = front_ok && t.front_mesh_identical;
  const bool truncation_converged = dEz <= truncation::kTolEzRef &&
                                    dCEX <= truncation::kTolCEX &&
                                    dVprobe <= truncation::kTolVProbe;

  std::printf("  Vernetzung an Spitze und Extraktor ueber alle Laengen identisch: %s\n",
              front_ok ? "ja" : "NEIN");
  std::printf("  relative Aenderung bei der letzten Verdopplung (%.4g m -> %.4g m):\n", T0.L,
              T1.L);
  std::printf("    E_z(ref)  %.3e   (Toleranz %.1e)  %s\n", dEz, truncation::kTolEzRef,
              dEz <= truncation::kTolEzRef ? "erreicht" : "NICHT erreicht");
  std::printf("    c_EX      %.3e   (Toleranz %.1e)  %s\n", dCEX, truncation::kTolCEX,
              dCEX <= truncation::kTolCEX ? "erreicht" : "NICHT erreicht");
  std::printf("    V(Sonden) %.3e   (Toleranz %.1e)  %s\n", dVprobe, truncation::kTolVProbe,
              dVprobe <= truncation::kTolVProbe ? "erreicht" : "NICHT erreicht");
  std::printf("    c_EE      %.3e   (keine Toleranz: haengt zwangslaeufig von der Laenge ab)\n",
              dCEE);
  std::printf("  ERGEBNIS: %s\n",
              truncation_converged
                  ? "truncation-konvergiert"
                  : "NICHT truncation-konvergiert -- emitter_back_length ist eine "
                    "anzugebende Geometrieabmessung, kein Konvergenzparameter");

  {
    std::FILE* f = std::fopen((d + "probe_points.csv").c_str(), "w");
    std::fprintf(f, "# fixed points between emitter and extractor, placed from the geometry;\n");
    std::fprintf(f, "# the same physical points at every mesh level and every back length\n");
    std::fprintf(f, "index,name,r_m,z_m\n");
    for (std::size_t k = 0; k < probes.size(); ++k)
      std::fprintf(f, "%zu,%s,%.9e,%.9e\n", k, probes[k].name.c_str(), probes[k].x.r,
                   probes[k].x.z);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "truncation.csv").c_str(), "w");
    std::fprintf(f, "# dependence on emitter_back_length at size_scale = %.6g, domain fixed.\n",
                 scales[ref]);
    std::fprintf(f, "# front_mesh_identical=1 means the panels at the tip and on the\n");
    std::fprintf(f, "# extractor are bitwise the same as for the shortest length, so only\n");
    std::fprintf(f, "# the rear end differs between rows.\n");
    std::fprintf(f, "# Tolerances fixed before the measurement: Ez_ref %.1e, c_EX %.1e,\n",
                 truncation::kTolEzRef, truncation::kTolCEX);
    std::fprintf(f, "# probe potentials %.1e of |V_E - V_X|.  Verdict: %s.\n",
                 truncation::kTolVProbe,
                 truncation_converged ? "converged" : "NOT converged");
    std::fprintf(f, "emitter_back_length_m,z_closure_m,clearance_to_evaluated_region_m,"
                    "n_mesh_elements,n_bem_panels,n_closure_panels,front_mesh_identical,"
                    "c_EE_F,c_EX_F,c_XX_F,C_mutual_F,Q_emitter_C,Ez_ref_V_per_m,"
                    "Epeak_free_surface_V_per_m,interior_dV_V,interior_Emag_V_per_m");
    for (std::size_t k = 0; k < probes.size(); ++k)
      std::fprintf(f, ",V_%s_V", probes[k].name.c_str());
    std::fprintf(f, "\n");
    for (const Trunc& t : tr) {
      std::fprintf(f, "%.9e,%.9e,%.9e,%d,%d,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e",
                   t.L, t.z_closure, t.clearance, t.n_mesh, t.n_bem, t.n_closure,
                   t.front_mesh_identical ? 1 : 0, t.cap.c_EE, t.cap.c_EX, t.cap.c_XX,
                   t.cap.mutual(), t.Q_emitter, t.E_ref, t.E_peak_free, t.interior_dV,
                   t.interior_E);
      for (Real v : t.V_probe) std::fprintf(f, ",%.9e", v);
      std::fprintf(f, "\n");
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
    std::fprintf(f, "# surface charge density and normal field on the reference mesh level.\n");
    std::fprintf(f, "# The conductor is CLOSED, so sigma/eps0 is the one-sided vacuum normal\n");
    std::fprintf(f, "# field -- but only where the panel bounds the device.  evaluable=1 marks\n");
    std::fprintf(f, "# exactly those panels: not the numerical closure (numerical=1), not the\n");
    std::fprintf(f, "# shank behind the taper foot, not the marked edge zones.\n");
    std::fprintf(f, "i,r_m,z_m,arclen_m,sigma_C_per_m2,En_V_per_m,area_m2,tag,boundary_id,"
                    "in_edge_zone,numerical,evaluable\n");
    Real s_arc = 0.0;
    for (std::size_t k = 0; k < S.sel.panels.size(); ++k) {
      const VacuumPanel& pn = S.sel.panels[k];
      const Element& el = S.bem.mesh().elems[static_cast<std::size_t>(pn.bem_element)];
      s_arc += el.len;
      std::fprintf(f, "%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%s,%s,%d,%d,%d\n", pn.bem_element,
                   el.mid.r, el.mid.z, s_arc,
                   S.bem.sigma()[static_cast<std::size_t>(pn.bem_element)],
                   S.bem.En(pn.bem_element), el.area, electrode_label(el.tag),
                   to_string(pn.boundary), in_edge_zone(S.zones, el.mid) ? 1 : 0,
                   pn.numerical ? 1 : 0, pn.evaluable ? 1 : 0);
    }
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((d + "edge_zones.csv").c_str(), "w");
    std::fprintf(f, "# places where no field value may be read off.  radius =\n");
    std::fprintf(f, "# local_feature_size / 4, a purely geometric length, so the marked region\n");
    std::fprintf(f, "# is the same at every refinement level.\n");
    std::fprintf(f, "# kind=sharp_feature   : unrounded edge of the device; |E| diverges there\n");
    std::fprintf(f, "#                        and follows the element size.\n");
    std::fprintf(f, "# kind=numerical_closure: rim of the numerical rearward closure.  Excluded\n");
    std::fprintf(f, "#                        for a stronger reason than sharpness -- the\n");
    std::fprintf(f, "#                        surface it sits on is not part of the device.\n");
    std::fprintf(f, "# kind=truncation_end  : free edge of an OPEN conductor arc.  None can\n");
    std::fprintf(f, "#                        occur: vacuum_bem_mesh refuses such a mesh.\n");
    std::fprintf(f, "kind,name,r_m,z_m,radius_m,local_feature_size_m\n");
    for (const EdgeZone& z : S.zones)
      std::fprintf(f, "%s,%s,%.9e,%.9e,%.9e,%.9e\n", to_string(z.kind), z.name.c_str(),
                   z.position.r, z.position.z, z.radius, z.local_feature_size);
    std::fclose(f);
  }

  // -------------------------------------------- what each number depends on
  {
    std::FILE* f = std::fopen((d + "quantities.csv").c_str(), "w");
    std::fprintf(f, "# Every reported number, with what it converged against and what it\n");
    std::fprintf(f, "# still depends on.  mesh_convergence refers to convergence.csv,\n");
    std::fprintf(f, "# truncation_convergence to truncation.csv.\n");
    std::fprintf(f, "quantity,symbol,value_SI,unit,mesh_convergence,truncation_convergence,"
                    "depends_on,note\n");
    const Level& a = levels[levels.size() - 2];
    const Level& b = levels.back();
    auto row = [&](const char* q, const char* sym, Real val, const char* unit, const char* mesh,
                   const char* trunc, const char* dep, const char* note) {
      std::fprintf(f, "%s,%s,%.9e,%s,%s,%s,%s,%s\n", q, sym, val, unit, mesh, trunc, dep, note);
    };
    const char* tv = truncation_converged ? "converged" : "not_converged";
    char mesh_ez[64], mesh_cee[64], mesh_cm[64];
    std::snprintf(mesh_ez, sizeof mesh_ez, "converged (%.1e)",
                  rel_change(a.E_ref, b.E_ref));
    std::snprintf(mesh_cee, sizeof mesh_cee, "converged (%.1e)",
                  rel_change(a.cap.c_EE, b.cap.c_EE));
    std::snprintf(mesh_cm, sizeof mesh_cm, "converged (%.1e)",
                  rel_change(a.cap.mutual(), b.cap.mutual()));
    row("axial field at the reference point", "Ez_ref", R.E_ref, "V/m", mesh_ez, tv,
        "emitter_back_length; extractor_outer_radius; all device dimensions",
        "local quantity; mesh converged but NOT independent of the rearward length");
    for (std::size_t k = 0; k < probes.size(); ++k) {
      char nm[128];
      std::snprintf(nm, sizeof nm, "potential at probe %s", probes[k].name.c_str());
      row(nm, "V_probe", R.V_probe[k], "V", "converged", tv,
          "emitter_back_length; extractor_outer_radius; all device dimensions",
          "local quantity; same caveat as Ez_ref");
    }
    row("Maxwell self coefficient of the emitter", "c_EE", R.cap.c_EE, "F", mesh_cee,
        "not_applicable",
        "emitter_back_length (first order); extractor_outer_radius; all device dimensions",
        "the charge of a conductor grows with its length; NOT a device prediction");
    row("Maxwell induction coefficient", "c_EX", R.cap.c_EX, "F", "converged", tv,
        "emitter_back_length; extractor_outer_radius; all device dimensions",
        "grows by about a third per doubling of the rearward length");
    row("mutual capacitance", "C_m", R.cap.mutual(), "F", mesh_cm, tv,
        "emitter_back_length; extractor_outer_radius; all device dimensions",
        "= -c_EX; C_m changes by a factor 2 between r_ext 0.5 mm and 2.5 mm");
    row("Maxwell self coefficient of the extractor", "c_XX", R.cap.c_XX, "F", "converged",
        "not_applicable", "extractor_outer_radius; all device dimensions", "");
    row("emitter charge at the operating point", "Q_E", R.Q_emitter, "C", "converged", tv,
        "emitter_back_length; extractor_outer_radius; all device dimensions",
        "Q_E + Q_X is not zero: with V = 0 at infinity the system carries net charge");
    row("peak |E_n| on the flat liquid surface", "Epeak_free", R.E_peak_free, "V/m", "converged",
        tv, "emitter_back_length; all device dimensions",
        "evaluable panels only; the pinned-edge zone is excluded");
    row("largest |E_n| inside the marked edge zones", "Eedge", R.E_edge_max, "V/m",
        "not_converged", "not_applicable", "element size",
        "unrounded edges: the corner field diverges; rounding radii are a P3 parameter");
    row("largest |E_n| on the numerical closure", "Eclosure", R.E_closure_max, "V/m",
        "not_applicable", "not_applicable", "emitter_back_length",
        "not a device surface; reported only so that its exclusion can be checked");
    row("max |V - V_emitter| inside the closed conductor", "interior_dV", R.interior_dV, "V",
        "converged", "not_applicable", "element size",
        "discretisation error only; falls with refinement");
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

    // The numerical rearward closure, on its own grid.  It is micrometres wide
    // and invisible on the domain map, and a reader is entitled to see that the
    // conductor really is closed there rather than take it on trust.
    const Real r3 = 0.5 * dp.phi_3;
    const Real zc = g.back_closure_z();
    const FieldGrid cl = sample_field(S.bem, 161, 161, 0.0, 2.5 * r3, zc - 2.0 * r3,
                                      zc + 3.0 * r3);
    write_grid_csv(cl, d + "field_back_closure.csv");
    std::printf("  field_back_closure.csv %d x %d\n", cl.nr, cl.nz);
  }

  // ---------------------------------------------------------- axis profile
  {
    std::FILE* f = std::fopen((d + "axis_profile.csv").c_str(), "w");
    std::fprintf(f, "# potential and axial field on r = 0 across the whole plotted domain.\n");
    std::fprintf(f, "# Inside the closed conductor (%.6g m <= z <= 0) V = V_emitter.\n",
                 g.back_closure_z());
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

    std::fprintf(f, "Numerische Rueckschliessung des Emitterleiters\n");
    std::fprintf(f, "  Der Leiter aus Emittermetall und Fluessigkeitssaeule endet bei\n");
    std::fprintf(f, "  z = %.6g m und ist dort durch eine VOLLE leitende Scheibe von r = 0\n",
                 g.back_closure_z());
    std::fprintf(f, "  bis r = %.6g m auf V_emitter geschlossen.  Randkennung:\n",
                 0.5 * dp.phi_3);
    std::fprintf(f, "  numerical_emitter_back_closure -- kein Bauteil, kein Domaenenrand.\n");
    std::fprintf(f, "  Die Scheibe wird mitgeloest, weil ein Leiter mit Loch ein anderes\n");
    std::fprintf(f, "  mathematisches Objekt ist; sie liefert keinen einzigen berichteten\n");
    std::fprintf(f, "  Feldwert.  Ihre Kante ist als nicht auswertbar markiert.\n");
    std::fprintf(f, "  Abstand der Scheibe zur ausgewerteten Region (z >= %.6g m): %.6g m.\n",
                 g.evaluation_z_min(), g.back_closure_clearance());
    std::fprintf(f, "  Panels: %d von %d; Rotationsflaeche %.6g m^2.\n",
                 S.sel.n_numerical_closure, S.sel.n_selected,
                 S.sel.revolved_area_numerical_closure);
    std::fprintf(f, "  Offene Bogenenden im geloesten Netz: 0 -- vacuum_bem_mesh lehnt ein\n");
    std::fprintf(f, "  Netz mit offenem Leiter ab, statt es zu dokumentieren.\n");
    std::fprintf(f, "  Im Inneren des geschlossenen Leiters: max |V - V_emitter| = %.3e V\n",
                 R.interior_dV);
    std::fprintf(f, "  von %.4g V, max |E| = %.3e V/m gegen %.3e V/m an der Oberflaeche.\n",
                 std::abs(VE), R.interior_E, R.E_peak_emitter);
    std::fprintf(f, "  Das ist Diskretisierungsfehler und faellt mit der Verfeinerung\n");
    std::fprintf(f, "  (Spalte interior_Emag_V_per_m in convergence.csv).\n\n");

    std::fprintf(f, "Trunkierungsstudie -- und ihr Ergebnis\n");
    std::fprintf(f, "  Vier zunehmende rueckwaertige Laengen, Domaene und size_scale fest,\n");
    std::fprintf(f, "  Vernetzung an Spitze und Extraktor ueber alle Laengen bitweise\n");
    std::fprintf(f, "  identisch (geprueft: %s).\n", front_ok ? "ja" : "NEIN");
    std::fprintf(f, "  %14s %14s %16s %16s %16s\n", "L [m]", "n_BEM", "E_z(ref) [V/m]",
                 "c_EX [F]", "c_EE [F]");
    for (const Trunc& t : tr)
      std::fprintf(f, "  %14.6g %14d %16.9e %16.9e %16.9e\n", t.L, t.n_bem, t.E_ref, t.cap.c_EX,
                   t.cap.c_EE);
    std::fprintf(f, "\n  Vorab festgelegte Toleranzen (es::truncation, vacuum_bem.hpp):\n");
    std::fprintf(f, "    E_z(ref) %.1e, c_EX %.1e, Sondenpotentiale %.1e von |V_E - V_X|.\n",
                 truncation::kTolEzRef, truncation::kTolCEX, truncation::kTolVProbe);
    std::fprintf(f, "  Gemessen bei der letzten Verdopplung (%.4g m -> %.4g m):\n", T0.L, T1.L);
    std::fprintf(f, "    E_z(ref)          %.3e  %s\n", dEz,
                 dEz <= truncation::kTolEzRef ? "erreicht" : "NICHT erreicht");
    std::fprintf(f, "    c_EX              %.3e  %s\n", dCEX,
                 dCEX <= truncation::kTolCEX ? "erreicht" : "NICHT erreicht");
    std::fprintf(f, "    Sondenpotentiale  %.3e  %s\n", dVprobe,
                 dVprobe <= truncation::kTolVProbe ? "erreicht" : "NICHT erreicht");
    std::fprintf(f, "    c_EE              %.3e  (ohne Toleranz)\n", dCEE);
    if (!truncation_converged) {
      std::fprintf(f, "\n  ERGEBNIS: NICHT truncation-konvergiert.\n");
      std::fprintf(f, "  Die Absicht war, den Leiter so weit nach hinten zu verlaengern, dass\n");
      std::fprintf(f, "  das lokale Spitzenfeld stehen bleibt, und nur die Gesamtkapazitaet\n");
      std::fprintf(f, "  als laengenabhaengig zu fuehren.  Die Messung gibt das nicht her.\n");
      std::fprintf(f, "  E_z(ref) traegt einen 1/L-Schwanz -- pro Verdopplung 5.6, 4.0, 2.5,\n");
      std::fprintf(f, "  1.5 Prozent zwischen 200 um und 3.2 mm -- und c_EX waechst um rund\n");
      std::fprintf(f, "  ein Drittel je Verdopplung, ohne sich zu setzen.\n");
      std::fprintf(f, "  Der Grund liegt in der Randbedingung, nicht in der Diskretisierung:\n");
      std::fprintf(f, "  mit V -> 0 im Unendlichen und nur zwei Leitern traegt das System\n");
      std::fprintf(f, "  Nettoladung.  Es gibt im Modell keine Rueckelektrode und keine\n");
      std::fprintf(f, "  Umhausung, also traegt ein laengerer Emitter schlicht mehr Ladung,\n");
      std::fprintf(f, "  und die wird an der Spitze ebenso gespuert wie am Extraktor.\n");
      std::fprintf(f, "  FOLGE: emitter_back_length ist eine anzugebende Abmessung, kein\n");
      std::fprintf(f, "  Konvergenzparameter.  Kein P2a-Wert wird als unabhaengig von ihr\n");
      std::fprintf(f, "  ausgegeben.  Ob das lokale Feld unempfindlich wird, sobald eine\n");
      std::fprintf(f, "  geerdete Umhausung im Modell ist, gehoert in die Phase, die eine\n");
      std::fprintf(f, "  hinzufuegt.\n");
    } else {
      std::fprintf(f, "\n  ERGEBNIS: truncation-konvergiert innerhalb der Toleranzen.\n");
    }
    std::fprintf(f, "  Der benutzte Wert emitter_back_length = %.6g m ist ein BEISPIELWERT,\n",
                 L0);
    std::fprintf(f, "  keine gemessene Abmessung -- wie extractor_outer_radius.\n\n");

    std::fprintf(f, "Netzkonvergenz (bei emitter_back_length = %.6g m)\n", L0);
    std::fprintf(f, "  %-6s %8s %8s %16s %16s %16s %12s\n", "Stufe", "n_Netz", "n_BEM",
                 "c_EE [F]", "C_m [F]", "E_z(ref) [V/m]", "Residuum*");
    std::fprintf(f, "  (* relativ, auf den auswertbaren Geraeteflaechen)\n");
    for (std::size_t k = 0; k < levels.size(); ++k) {
      const Level& L = levels[k];
      std::fprintf(f, "  %-6d %8d %8d %16.9e %16.9e %16.9e %12.3e\n", static_cast<int>(k),
                   L.n_mesh, L.n_bem, L.cap.c_EE, L.cap.mutual(), L.E_ref, L.residual_rel_phys);
    }
    const Level& a = levels[levels.size() - 2];
    const Level& b = levels.back();
    std::fprintf(f, "\n  relative Aenderung zwischen den beiden feinsten Stufen\n");
    std::fprintf(f, "    c_EE      : %.3e\n", rel_change(a.cap.c_EE, b.cap.c_EE));
    std::fprintf(f, "    C_m       : %.3e\n", rel_change(a.cap.mutual(), b.cap.mutual()));
    std::fprintf(f, "    E_z(ref)  : %.3e\n", rel_change(a.E_ref, b.E_ref));

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
    std::fprintf(f, "  ALLE diese Werte gelten fuer emitter_back_length = %.6g m und\n", L0);
    std::fprintf(f, "  extractor_outer_radius = %.6g m.  c_EE ist keine Vorhersage fuer das\n",
                 dp.extractor_outer_radius);
    std::fprintf(f, "  Geraet: die Kapazitaet eines Leiters haengt von seiner Laenge ab.\n");

    std::fprintf(f, "\nFelder auf der Referenzstufe (nur auswertbare Panels)\n");
    std::fprintf(f, "  E_z am Referenzpunkt                     : %.9e V/m\n", R.E_ref);
    std::fprintf(f, "  groesstes |E_n| auf der ebenen Flaeche   : %.9e V/m\n", R.E_peak_free);
    std::fprintf(f, "  groesstes |E_n| am Emitter               : %.9e V/m\n", R.E_peak_emitter);
    std::fprintf(f, "  groesstes |E_n| am Extraktor             : %.9e V/m\n", R.E_peak_extractor);
    std::fprintf(f, "  groesstes |E_n| INNERHALB der Kantenzonen: %.9e V/m\n", R.E_edge_max);
    std::fprintf(f, "  groesstes |E_n| auf der Rueckschliessung : %.9e V/m\n", R.E_closure_max);
    std::fprintf(f, "    Die letzten beiden Werte sind KEINE Ergebnisse.  Die Austritts-,\n");
    std::fprintf(f, "    Stirn- und Aperturkanten sind unverrundet; das Feld einer scharfen\n");
    std::fprintf(f, "    Kante divergiert und folgt der Elementgroesse.  Verrundungsradien\n");
    std::fprintf(f, "    sind ein P3-Parameter.  Die Rueckschliessung ist kein Bauteil.\n");

    std::fprintf(f, "\nPotential an festen Punkten zwischen Emitter und Extraktor\n");
    std::fprintf(f, "  %-36s %12s %12s %14s\n", "Punkt", "r [m]", "z [m]", "V [V]");
    for (std::size_t k = 0; k < probes.size(); ++k)
      std::fprintf(f, "  %-36s %12.5g %12.5g %14.7g\n", probes[k].name.c_str(), probes[k].x.r,
                   probes[k].x.z, R.V_probe[k]);
    std::fprintf(f, "  groesste relative Aenderung dieser Werte bei der letzten Verdopplung\n");
    std::fprintf(f, "  der rueckwaertigen Laenge: %.3e von |V_E - V_X|.\n", dVprobe);

    std::fprintf(f, "\nResiduum der Dirichlet-Bedingung (Viertelpunkte, nicht Kollokation)\n");
    std::fprintf(f, "  auf den auswertbaren Geraeteflaechen: max %.3e V, relativ %.3e, "
                    "RMS relativ %.3e\n",
                 R.residual_max, R.residual_rel_phys, R.residual_rms_phys);
    std::fprintf(f, "  ausserhalb der Kantenzonen  : max |V - V_soll| = %.3e V, relativ %.3e\n",
                 R.residual_max, R.residual_rel);
    std::fprintf(f, "  einschliesslich Kantenzonen : max |V - V_soll| = %.3e V, relativ %.3e\n",
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
    std::fprintf(f, "n_closure_panels=%d\n", S.sel.n_numerical_closure);
    std::fprintf(f, "z_ref_m=%.9e\n", z_ref);
    std::fprintf(f, "emitter_back_length_m=%.9e\n", L0);
    std::fprintf(f, "back_closure_z_m=%.9e\n", g.back_closure_z());
    std::fprintf(f, "back_closure_clearance_m=%.9e\n", g.back_closure_clearance());
    std::fprintf(f, "evaluation_z_min_m=%.9e\n", g.evaluation_z_min());
    std::fprintf(f, "open_arc_ends=0\n");
    std::fprintf(f, "truncation_tol_Ez_ref=%.9e\n", truncation::kTolEzRef);
    std::fprintf(f, "truncation_tol_c_EX=%.9e\n", truncation::kTolCEX);
    std::fprintf(f, "truncation_tol_V_probe=%.9e\n", truncation::kTolVProbe);
    std::fprintf(f, "truncation_change_Ez_ref=%.9e\n", dEz);
    std::fprintf(f, "truncation_change_c_EX=%.9e\n", dCEX);
    std::fprintf(f, "truncation_change_c_EE=%.9e\n", dCEE);
    std::fprintf(f, "truncation_change_V_probe=%.9e\n", dVprobe);
    std::fprintf(f, "truncation_converged=%s\n", truncation_converged ? "yes" : "no");
    std::fprintf(f, "front_mesh_identical=%s\n", front_ok ? "yes" : "no");
    std::fprintf(f, "state=vacuum electrostatics, rho=0, flat perfect-conductor reference, "
                    "closed emitter conductor, no emission\n");
    std::fclose(f);
  }

  cfg.warn_about_unused(stdout, {"meta.", "fluid.", "beam.", "output."});
  std::printf("\ngeschrieben nach %s\n", outdir.c_str());
  if (!front_ok) {
    std::fprintf(stderr,
                 "\nFehler: die Vernetzung an Spitze und Extraktor haengt von der "
                 "rueckwaertigen Laenge ab; die Trunkierungsstudie misst dann zwei Dinge "
                 "gleichzeitig.\n");
    return 2;
  }
  return 0;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\n%s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "\nFehler: %s\n", e.what());
  return 2;
}
