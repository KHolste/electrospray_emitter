#include "es/transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::eps0;
using constants::pi;

namespace {

constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();

/// A uniform cylinder mesh: nr columns in r, nz rows in z, level rows.
QuadMesh cylinder_mesh(Real R, Real L, Index nr, Index nz) {
  QuadMesh m;
  m.nr = nr;
  m.nz = nz;
  m.nodes.resize(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nz));
  for (Index j = 0; j < nz; ++j)
    for (Index i = 0; i < nr; ++i)
      m.nodes[static_cast<std::size_t>(j * nr + i)] =
          Vec2{R * static_cast<Real>(i) / static_cast<Real>(nr - 1),
               L * static_cast<Real>(j) / static_cast<Real>(nz - 1)};
  return m;
}

}  // namespace

// ===========================================================================
// 1.  Fully developed pipe flow
// ===========================================================================

PipeFlowSolution solve_pipe_flow(Real R, Real L, Real mu, Real dpdz, Index nr, Index nz) {
  if (!(R > 0.0) || !(L > 0.0) || !(mu > 0.0))
    throw std::runtime_error("solve_pipe_flow: R, L und mu muessen positiv sein");
  if (nr < 3 || nz < 3) throw std::runtime_error("solve_pipe_flow: Netz zu grob");

  PipeFlowSolution out;
  out.radius = R;
  out.length = L;
  out.mu = mu;
  out.dpdz = dpdz;

  const QuadMesh m = cylinder_mesh(R, L, nr, nz);
  out.n_nodes = m.n_nodes();
  out.n_cells = m.n_cells();

  AxisymProblem p;
  p.mesh = &m;
  // THE COEFFICIENT IS THE VISCOSITY and the scale is one: see the header of
  // AxisymProblem.  Nothing about the assembly changes; only what the numbers
  // mean does.
  p.coefficient_scale = 1.0;
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), mu);
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  // -div(mu grad u) = -dp/dz
  p.cell_source.assign(static_cast<std::size_t>(m.n_cells()), -dpdz);
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  // No slip on the wall.  The two ends and the axis carry the NATURAL condition
  // du/dn = 0, which is exactly right here: the profile is z-independent, so
  // du/dz = 0 on the ends, and regularity gives du/dr = 0 on the axis.  Nothing
  // is imposed that the exact solution does not satisfy.
  for (Index j = 0; j < m.nz; ++j) {
    const Index n = m.node(m.nr - 1, j);
    p.fixed[static_cast<std::size_t>(n)] = 1;
    p.fixed_value[static_cast<std::size_t>(n)] = 0.0;
  }
  p.far_field = FarField::Insulated;
  p.check();

  const AxisymSolution sol = solve_axisym(p);
  out.u = sol.phi;
  out.fem_residual = sol.residual_inf;

  // --- the closed forms this is judged against ------------------------------
  const Real u0 = -dpdz * R * R / (4.0 * mu);
  out.centreline_closed_form = u0;
  out.flow_rate_closed_form = -dpdz * pi * std::pow(R, 4.0) / (8.0 * mu);

  out.centreline_velocity = out.u[static_cast<std::size_t>(m.node(0, m.nz / 2))];

  // Flow rate by integrating the nodal profile over one cross section with the
  // trapezoidal rule against 2 pi r dr.  It is a different operation from the
  // closed form and from the FEM assembly.
  {
    const Index j = m.nz / 2;
    Real q = 0.0;
    for (Index i = 0; i + 1 < m.nr; ++i) {
      const Real r0 = m.at(i, j).r, r1 = m.at(i + 1, j).r;
      const Real u0n = out.u[static_cast<std::size_t>(m.node(i, j))];
      const Real u1n = out.u[static_cast<std::size_t>(m.node(i + 1, j))];
      // int_{r0}^{r1} u(r) 2 pi r dr with u linear in r.
      q += 2.0 * pi * (r1 - r0) * ((2.0 * r0 * u0n + r0 * u1n + r1 * u0n + 2.0 * r1 * u1n) / 6.0);
    }
    out.flow_rate = q;
  }
  out.mean_velocity = out.flow_rate / (pi * R * R);

  // Wall shear stress from the last cell's gradient: tau_w = mu |du/dr| at r=R.
  {
    const Index j = m.nz / 2;
    const Real r0 = m.at(m.nr - 2, j).r, r1 = m.at(m.nr - 1, j).r;
    const Real ua = out.u[static_cast<std::size_t>(m.node(m.nr - 2, j))];
    const Real ub = out.u[static_cast<std::size_t>(m.node(m.nr - 1, j))];
    out.wall_shear_stress = mu * std::abs((ub - ua) / (r1 - r0));
  }

  // Profile error against the parabola, on the mid plane.
  {
    const Index j = m.nz / 2;
    Real worst = 0.0;
    for (Index i = 0; i < m.nr; ++i) {
      const Real r = m.at(i, j).r;
      const Real want = u0 * (1.0 - (r / R) * (r / R));
      const Real got = out.u[static_cast<std::size_t>(m.node(i, j))];
      if (std::abs(want) > 1e-12 * std::abs(u0))
        worst = std::max(worst, std::abs(got - want) / std::abs(u0));
    }
    out.max_profile_error = worst;
  }
  out.flow_rate_error = (std::abs(out.flow_rate_closed_form) > 0.0)
                            ? std::abs(out.flow_rate - out.flow_rate_closed_form) /
                                  std::abs(out.flow_rate_closed_form)
                            : kNaN;
  return out;
}

// ===========================================================================
// 2.  Charge transport
// ===========================================================================

Real charge_relaxation_time(Real eps_r, Real sigma) {
  if (!(eps_r > 0.0) || !(sigma > 0.0)) return kNaN;
  return eps0 * eps_r / sigma;
}

Real relaxed_charge_density(Real rho0, Real t, Real tau) {
  if (!(tau > 0.0)) return kNaN;
  return rho0 * std::exp(-t / tau);
}

const char* to_string(ConductorLimit v) {
  switch (v) {
    case ConductorLimit::NotAttempted: return "NotAttempted";
    case ConductorLimit::PerfectConductorJustified: return "PerfectConductorJustified";
    case ConductorLimit::FiniteConductivityRequired: return "FiniteConductivityRequired";
    case ConductorLimit::MissingMaterialData: return "MissingMaterialData";
  }
  return "?";
}

const char* explain(ConductorLimit v) {
  switch (v) {
    case ConductorLimit::NotAttempted:
      return "Der Grenzfall wurde nicht beurteilt.";
    case ConductorLimit::PerfectConductorJustified:
      return "Die Ladungsrelaxationszeit ist um mehr als den festgelegten Faktor kuerzer als "
             "die angegebene Prozesszeit.  Die Behandlung der Fluessigkeit als "
             "Aequipotentialflaeche ist FUER DIESE Prozesszeit gerechtfertigt -- und fuer "
             "keine andere.";
    case ConductorLimit::FiniteConductivityRequired:
      return "Die Ladungsrelaxationszeit ist vergleichbar mit der Prozesszeit oder laenger.  "
             "Die Behandlung als idealer Leiter ist dann eine Annahme und kein Grenzfall; "
             "das Problem braucht endliche Leitfaehigkeit.";
    case ConductorLimit::MissingMaterialData:
      return "Die Permittivitaet oder die Leitfaehigkeit fehlt.  tau ist nicht berechenbar, "
             "und es wird KEIN Ersatzwert gesetzt: ohne beide Zahlen ist der Grenzfall "
             "weder gerechtfertigt noch widerlegt, sondern unbekannt.";
  }
  return "?";
}

void RelaxationVerdict::print(std::FILE* out) const {
  std::fprintf(out, "  Ladungsrelaxation: %s\n", to_string(limit));
  std::fprintf(out, "    eps_r = %.6g [%s], sigma = %.6g S/m [%s]\n", eps_r,
               to_string(eps_status), sigma, to_string(sigma_status));
  if (std::isfinite(tau)) {
    std::fprintf(out, "    tau = eps0 eps_r / sigma = %.6e s\n", tau);
    std::fprintf(out, "    Prozesszeit %.6e s, Verhaeltnis t/tau = %.4g (verlangt > %.0f)\n",
                 process_time, ratio, transport::kPerfectConductorMargin);
  } else {
    std::fprintf(out, "    tau ist nicht berechenbar.\n");
  }
  std::fprintf(out, "    %s\n", message.c_str());
}

RelaxationVerdict judge_conductor_limit_explicit(Real eps_r, Real sigma, Real process_time) {
  RelaxationVerdict v;
  v.eps_r = eps_r;
  v.sigma = sigma;
  v.process_time = process_time;
  v.eps_status = (eps_r > 0.0) ? MaterialDataStatus::Literature
                               : MaterialDataStatus::MissingMaterialData;
  v.sigma_status = (sigma > 0.0) ? MaterialDataStatus::Literature
                                 : MaterialDataStatus::MissingMaterialData;
  v.tau = charge_relaxation_time(eps_r, sigma);
  if (!std::isfinite(v.tau)) {
    v.limit = ConductorLimit::MissingMaterialData;
    v.message = explain(v.limit);
    v.ratio = kNaN;
    return v;
  }
  v.ratio = process_time / v.tau;
  v.limit = (v.ratio > transport::kPerfectConductorMargin)
                ? ConductorLimit::PerfectConductorJustified
                : ConductorLimit::FiniteConductivityRequired;
  v.message = explain(v.limit);
  return v;
}

RelaxationVerdict judge_conductor_limit(const MaterialDataset& d, Real T, Real process_time) {
  const MaterialValue e = material_value(d, PropertyKind::RelativePermittivity, T);
  const MaterialValue s = material_value(d, PropertyKind::ElectricalConductivity, T);
  RelaxationVerdict v;
  v.process_time = process_time;
  v.eps_status = e.status;
  v.sigma_status = s.status;
  v.eps_r = e.usable() ? e.value : 0.0;
  v.sigma = s.usable() ? s.value : 0.0;
  if (!e.usable() || !s.usable()) {
    v.limit = ConductorLimit::MissingMaterialData;
    v.tau = kNaN;
    v.ratio = kNaN;
    v.message = std::string(explain(v.limit)) + "  ";
    if (!e.usable()) v.message += "eps_r: " + e.message + "  ";
    if (!s.usable()) v.message += "sigma: " + s.message;
    return v;
  }
  return judge_conductor_limit_explicit(v.eps_r, v.sigma, process_time);
}

// ===========================================================================
// 2b.  Which permittivity belongs in tau_q
// ===========================================================================

const char* to_string(PermittivityConcept c) {
  switch (c) {
    case PermittivityConcept::StaticApparentLowFrequency: return "static_apparent_low_frequency";
    case PermittivityConcept::IntrinsicStatic: return "intrinsic_static";
    case PermittivityConcept::FrequencyResolved: return "frequency_resolved_complex";
    case PermittivityConcept::ElectrodePolarisation: return "electrode_polarisation";
    case PermittivityConcept::DcConductivity: return "dc_conductivity";
  }
  return "?";
}

const char* explain(PermittivityConcept c) {
  switch (c) {
    case PermittivityConcept::StaticApparentLowFrequency:
      return "Was eine Kapazitaetsmessbruecke bei kHz anzeigt.  Bei K in der Groessenordnung "
             "1 S/m wird das vom Leitungsbeitrag und von der Elektrodenpolarisation "
             "beherrscht und erreicht Werte von 10^4 und mehr.  Eine Eigenschaft der "
             "MESSZELLE, nicht der Fluessigkeit; kein Modell darf sie benutzen.";
    case PermittivityConcept::IntrinsicStatic:
      return "eps_s: der Grenzwert der DIELEKTRISCHEN Dispersion fuer f -> 0, nachdem der "
             "Leitungsbeitrag abgetrennt wurde.  Bei ionischen Fluessigkeiten wird er nicht "
             "bei null Hertz gemessen, sondern durch Anpassung eines Relaxationsmodells an "
             "Mikrowellenspektren gewonnen und extrapoliert.";
    case PermittivityConcept::FrequencyResolved:
      return "eps*(f) = eps'(f) - i eps''(f), an einzelnen Frequenzen gemessen.  Traegt keinen "
             "Gleichstromwert, ist aber genau dann die richtige Groesse, wenn der Vorgang "
             "selbst bei dieser Frequenz stattfindet.";
    case PermittivityConcept::ElectrodePolarisation:
      return "Messartefakt: Ionen sammeln sich an den Elektroden und schirmen das angelegte "
             "Feld ab, was die scheinbare Kapazitaet aufblaeht.  Der Grund, warum fuer diese "
             "Fluessigkeiten ueberhaupt Mikrowellenspektroskopie benutzt wird.";
    case PermittivityConcept::DcConductivity:
      return "K, eine eigene Groesse mit eigener Quelle.  Sie als imaginaere Permittivitaet "
             "K/(eps0 omega) zu schreiben ist Buchhaltung, keine zweite Dielektrizitaetszahl.";
  }
  return "?";
}

namespace {

/// One admissible dielectric datum at temperature T.
struct DielectricPoint {
  Real f{0};        ///< [Hz]; 0 = reported as static (concept 2)
  Real value{0};
  Real uncertainty{0};
  const PropertySource* src{nullptr};
};

/// Collect the permittivity points at T that may back a DIELECTRIC statement.
/// Refused: points below the electrode-polarisation floor (concept 1), points
/// away from ambient pressure, points at another temperature.
std::vector<DielectricPoint> dielectric_points(const MaterialDataset& d, Real T, Real tol,
                                               bool* any_below_floor) {
  std::vector<DielectricPoint> out;
  if (any_below_floor) *any_below_floor = false;
  const MaterialProperty* p = d.find(PropertyKind::RelativePermittivity);
  if (p == nullptr) return out;
  for (std::size_t k = 0; k < p->n_sources; ++k) {
    const PropertySource& src = p->sources[k];
    for (std::size_t j = 0; j < src.n_points; ++j) {
      const PropertyPoint& q = src.points[j];
      if (!q.ambient()) continue;
      if (std::abs(q.T - T) > tol) continue;
      if (q.frequency_resolved() && q.frequency_Hz < transport::kElectrodePolarisationFloor) {
        if (any_below_floor) *any_below_floor = true;
        continue;   // concept (1): electrode polarisation, not a dielectric datum
      }
      out.push_back({q.frequency_Hz, q.value, q.uncertainty, &src});
    }
  }
  return out;
}

/// eps_r'(f) from the measured dispersion, log-linear between the two
/// neighbouring measured frequencies.  Outside the measured range it returns
/// the nearest measured value and sets `extrapolated` -- it does NOT invent a
/// dispersion model, and the caller reports the fact.
Real eps_at_frequency(const std::vector<DielectricPoint>& pts, Real f, bool* extrapolated) {
  std::vector<const DielectricPoint*> fr;
  for (const DielectricPoint& q : pts)
    if (q.f > 0.0) fr.push_back(&q);
  if (fr.empty()) return kNaN;
  std::sort(fr.begin(), fr.end(),
            [](const DielectricPoint* a, const DielectricPoint* b) { return a->f < b->f; });
  if (extrapolated) *extrapolated = false;
  if (f <= fr.front()->f) {
    if (extrapolated) *extrapolated = (f < fr.front()->f);
    return fr.front()->value;
  }
  if (f >= fr.back()->f) {
    if (extrapolated) *extrapolated = (f > fr.back()->f);
    return fr.back()->value;
  }
  for (std::size_t i = 1; i < fr.size(); ++i) {
    if (f <= fr[i]->f) {
      const Real w = (std::log(f) - std::log(fr[i - 1]->f)) /
                     (std::log(fr[i]->f) - std::log(fr[i - 1]->f));
      return (1.0 - w) * fr[i - 1]->value + w * fr[i]->value;
    }
  }
  return kNaN;
}

}  // namespace

void PermittivityBand::print(std::FILE* out) const {
  if (!ok) {
    std::fprintf(out, "  Permittivitaetsband: NICHT bestimmbar -- %s\n", blocker.c_str());
    return;
  }
  std::fprintf(out, "  eps_r bei %.2f K: begruendetes Band %.3g .. %.3g ueber %zu Quellen\n", T,
               lo, hi, n_sources);
  std::fprintf(out, "    %zu als statisch berichtete Punkte, %zu frequenzaufgeloeste "
                    "(%.3g .. %.3g Hz)\n",
               n_static_points, n_frequency_points, f_min, f_max);
  std::fprintf(out, "    Punkte unterhalb der Elektrodenpolarisationsschwelle (%.0e Hz): %s\n",
               transport::kElectrodePolarisationFloor,
               any_below_polarisation_floor ? "JA -- ausgeschlossen" : "keine im Datensatz");
  std::fprintf(out, "    %s\n", basis.c_str());
  std::fprintf(out, "    KEIN Einzelwert ist ausgewaehlt: keine Quelle nennt Reinheit und "
                    "Wassergehalt.\n");
}

PermittivityBand permittivity_band(const MaterialDataset& d, Real T) {
  PermittivityBand b;
  b.T = T;
  bool below = false;
  const std::vector<DielectricPoint> pts = dielectric_points(d, T, 2.0, &below);
  b.any_below_polarisation_floor = below;
  if (pts.empty()) {
    b.blocker = "kein zulaessiger dielektrischer Messpunkt bei dieser Temperatur.";
    return b;
  }
  std::vector<const PropertySource*> srcs;
  b.lo = pts[0].value;
  b.hi = pts[0].value;
  b.f_min = std::numeric_limits<Real>::max();
  b.f_max = 0.0;
  for (const DielectricPoint& q : pts) {
    b.lo = std::min(b.lo, q.value);
    b.hi = std::max(b.hi, q.value);
    if (q.f > 0.0) {
      ++b.n_frequency_points;
      b.f_min = std::min(b.f_min, q.f);
      b.f_max = std::max(b.f_max, q.f);
    } else {
      ++b.n_static_points;
    }
    if (std::find(srcs.begin(), srcs.end(), q.src) == srcs.end()) srcs.push_back(q.src);
  }
  if (b.n_frequency_points == 0) { b.f_min = 0.0; b.f_max = 0.0; }
  b.n_sources = srcs.size();
  b.ok = true;
  b.basis = "Band ueber alle als dielektrisch zulaessigen Punkte: die als statisch "
            "berichteten Werte (aus Mikrowellenspektren extrapoliert) und die "
            "frequenzaufgeloesten Messpunkte.  Es wird nicht gemittelt und keine "
            "Dispersionsfunktion angepasst.";
  return b;
}

Real dielectric_permittivity_at(const MaterialDataset& d, Real T, Real frequency_Hz,
                                bool* extrapolated) {
  const std::vector<DielectricPoint> pts = dielectric_points(d, T, 2.0, nullptr);
  return eps_at_frequency(pts, frequency_Hz, extrapolated);
}

// ---------------------------------------------------------------------------

void SelfConsistentRelaxation::print(std::FILE* out) const {
  if (!ok) {
    std::fprintf(out, "  selbstkonsistentes tau: NICHT bestimmbar -- %s\n", blocker.c_str());
    return;
  }
  std::fprintf(out, "  tau = eps0 eps_r(f*) / K, implizit mit f* = 1/(2 pi tau):\n");
  std::fprintf(out, "    tau   = %.6e s\n", tau);
  std::fprintf(out, "    f*    = %.6e Hz  (%s des Messbereichs %.3g .. %.3g Hz)\n", f_star,
               f_star_inside_measured ? "innerhalb" : "AUSSERHALB", f_measured_min,
               f_measured_max);
  std::fprintf(out, "    eps_r = %.4f bei f*\n", eps_r);
  std::fprintf(out, "    K     = %.6g S/m (ausgewaehlte Quelle)\n", sigma);
  std::fprintf(out, "    %d Iterationen; Selbstkonsistenz-Residuum "
                    "|eps_r(f*) - eps_r| / eps_r = %.3e\n",
               iterations, residual);
  std::fprintf(out, "    zum Vergleich mit dem intrinsisch STATISCHEN eps_s = %.4f: "
                    "tau = %.6e s (%+.1f %%)\n",
               eps_static, tau_static, 100.0 * (tau_static - tau) / tau);
}

SelfConsistentRelaxation self_consistent_relaxation(const MaterialDataset& d, Real T) {
  SelfConsistentRelaxation r;
  r.T = T;
  const MaterialValue sv = material_value(d, PropertyKind::ElectricalConductivity, T);
  if (!sv.usable()) {
    r.blocker = "die elektrische Leitfaehigkeit ist nicht belegt: " + sv.message;
    return r;
  }
  r.sigma = sv.value;

  bool below = false;
  const std::vector<DielectricPoint> pts = dielectric_points(d, T, 2.0, &below);
  Real fmin = std::numeric_limits<Real>::max(), fmax = 0.0;
  std::size_t n_static = 0, n_fr = 0;
  for (const DielectricPoint& q : pts) {
    if (q.f > 0.0) { ++n_fr; fmin = std::min(fmin, q.f); fmax = std::max(fmax, q.f); }
    else { ++n_static; }
  }
  if (n_fr == 0) {
    r.blocker = "es gibt keine frequenzaufgeloeste Permittivitaetsmessung; die implizite "
                "Gleichung fuer tau laesst sich dann nicht auf einer gemessenen Kurve "
                "loesen.";
    return r;
  }
  r.f_measured_min = fmin;
  r.f_measured_max = fmax;

  // The comparison value: the LARGEST of the intrinsic static values, i.e. the
  // one that makes tau longest and the perfect-conductor case weakest.  Not an
  // average -- this file does not average sources.
  Real eps_s = 0.0;
  for (const DielectricPoint& q : pts)
    if (q.f == 0.0) eps_s = std::max(eps_s, q.value);
  if (!(eps_s > 0.0) && n_static == 0) eps_s = eps_at_frequency(pts, fmin, nullptr);
  r.eps_static = eps_s;
  r.tau_static = eps0 * eps_s / r.sigma;

  // Fixed-point iteration.  Started from the static value, which is the
  // largest permittivity in play and therefore the slowest tau: the iteration
  // then walks DOWN onto the dispersion curve.
  Real tau = eps0 * eps_s / r.sigma;
  Real eps = eps_s;
  bool extrap = false;
  int it = 0;
  for (; it < 400; ++it) {
    const Real f_star = 1.0 / (2.0 * pi * tau);
    const Real e = eps_at_frequency(pts, f_star, &extrap);
    if (!std::isfinite(e) || !(e > 0.0)) {
      r.blocker = "die gemessene Dispersionskurve liefert bei der gesuchten Frequenz keinen "
                  "Wert.";
      return r;
    }
    const Real tau_new = eps0 * e / r.sigma;
    const Real step = std::abs(tau_new - tau) / tau_new;
    tau = tau_new;
    eps = e;
    if (step == 0.0) { ++it; break; }
  }
  r.iterations = it;

  // The three quantities cannot all three be exact at once in floating point,
  // so the two that are DEFINITIONS are made exact and the residual is reported
  // where the approximation actually sits: in the self-consistency of eps_r.
  //
  //   tau    = eps0 eps_r / K            -- exact, this is the formula
  //   f*     = 1 / (2 pi tau)            -- exact, this is the definition of f*
  //   eps_r  = eps_r(f*)                 -- satisfied to `residual`
  r.eps_r = eps;
  r.tau = eps0 * r.eps_r / r.sigma;
  r.f_star = 1.0 / (2.0 * pi * r.tau);
  const Real e_check = eps_at_frequency(pts, r.f_star, &extrap);
  r.residual = std::abs(e_check - r.eps_r) / r.eps_r;
  r.f_star_inside_measured = !extrap;
  r.ok = true;
  return r;
}

// ---------------------------------------------------------------------------

void BandedRelaxationVerdict::print(std::FILE* out) const {
  std::fprintf(out, "  Perfect-Conductor-Urteil ueber das ganze begruendete Band:\n");
  std::fprintf(out, "    eps_r %.3g .. %.3g,  K %.4g .. %.4g S/m\n", eps_lo, eps_hi, sigma_lo,
               sigma_hi);
  std::fprintf(out, "    tau   %.4e .. %.4e s  (Ecken des Bandes)\n", tau_min, tau_max);
  std::fprintf(out, "    selbstkonsistent: tau = %.4e s\n", tau_self_consistent);
  std::fprintf(out, "    t_process = %.4e s -> Verhaeltnis %.4g .. %.4g, Schranke %.0f\n",
               process_time, ratio_min, ratio_max, transport::kPerfectConductorMargin);
  std::fprintf(out, "    %s: %s\n", to_string(limit), message.c_str());
}

BandedRelaxationVerdict judge_conductor_limit_over_band(const MaterialDataset& d, Real T,
                                                        Real process_time) {
  BandedRelaxationVerdict v;
  v.T = T;
  v.process_time = process_time;

  const PermittivityBand band = permittivity_band(d, T);
  const MaterialProperty* pk = d.find(PropertyKind::ElectricalConductivity);
  if (!band.ok || pk == nullptr) {
    v.limit = ConductorLimit::MissingMaterialData;
    v.message = band.ok ? "keine Leitfaehigkeitsdaten." : band.blocker;
    v.tau_min = v.tau_max = v.ratio_min = v.ratio_max = kNaN;
    return v;
  }
  v.eps_lo = band.lo;
  v.eps_hi = band.hi;
  v.sigma_lo = pk->min_at(T);
  v.sigma_hi = pk->max_at(T);
  if (!(std::isfinite(v.sigma_lo) && v.sigma_lo > 0.0 && std::isfinite(v.sigma_hi))) {
    v.limit = ConductorLimit::MissingMaterialData;
    v.message = "das Leitfaehigkeitsband bei dieser Temperatur ist leer.";
    v.tau_min = v.tau_max = v.ratio_min = v.ratio_max = kNaN;
    return v;
  }

  // Four corners.  tau grows with eps_r and falls with K, so the worst corner
  // for the approximation is (eps_hi, sigma_lo) -- but it is computed, not
  // asserted, because a sign error in that reasoning would be invisible.
  const Real taus[4] = {eps0 * v.eps_lo / v.sigma_lo, eps0 * v.eps_lo / v.sigma_hi,
                        eps0 * v.eps_hi / v.sigma_lo, eps0 * v.eps_hi / v.sigma_hi};
  v.tau_min = taus[0];
  v.tau_max = taus[0];
  for (int i = 1; i < 4; ++i) {
    v.tau_min = std::min(v.tau_min, taus[i]);
    v.tau_max = std::max(v.tau_max, taus[i]);
  }
  v.ratio_min = process_time / v.tau_max;   // the worst case
  v.ratio_max = process_time / v.tau_min;

  const SelfConsistentRelaxation sc = self_consistent_relaxation(d, T);
  v.tau_self_consistent = sc.ok ? sc.tau : kNaN;

  if (v.ratio_min > transport::kPerfectConductorMargin) {
    v.limit = ConductorLimit::PerfectConductorJustified;
    v.message = "auch an der fuer die Naeherung UNGUENSTIGSTEN Ecke des Bandes ist tau um "
                "mehr als den geforderten Faktor kuerzer als die Prozesszeit.  Die "
                "Aequipotentialbehandlung ist damit belegt und nicht mehr nur angenommen -- "
                "und zwar ohne dass ein einzelner unbelegter eps_r-Wert benutzt wird.";
  } else {
    v.limit = ConductorLimit::FiniteConductivityRequired;
    v.message = "an mindestens einer Ecke des Bandes ist tau nicht klein genug.  Die "
                "Aequipotentialbehandlung ist dort nicht gerechtfertigt.";
  }
  return v;
}

// ===========================================================================
// 3.  Steady conduction in a cylinder
// ===========================================================================

ConductionSolution solve_cylinder_conduction(Real R, Real L, Real sigma, Real V, Index nr,
                                             Index nz) {
  if (!(R > 0.0) || !(L > 0.0) || !(sigma > 0.0))
    throw std::runtime_error("solve_cylinder_conduction: R, L und sigma muessen positiv sein");
  if (nr < 3 || nz < 3) throw std::runtime_error("solve_cylinder_conduction: Netz zu grob");

  ConductionSolution out;
  out.radius = R;
  out.length = L;
  out.sigma = sigma;
  out.voltage = V;

  const QuadMesh m = cylinder_mesh(R, L, nr, nz);
  out.n_nodes = m.n_nodes();

  AxisymProblem p;
  p.mesh = &m;
  // THE COEFFICIENT IS THE CONDUCTIVITY.  div(sigma grad phi) = 0 is the same
  // equation as div(eps grad phi) = 0; the nodal reaction is then a current in
  // ampere rather than a charge in coulomb.
  p.coefficient_scale = 1.0;
  p.eps_r.assign(static_cast<std::size_t>(m.n_cells()), sigma);
  p.active.assign(static_cast<std::size_t>(m.n_cells()), 1);
  p.fixed.assign(static_cast<std::size_t>(m.n_nodes()), 0);
  p.fixed_value.assign(static_cast<std::size_t>(m.n_nodes()), 0.0);
  for (Index i = 0; i < m.nr; ++i) {
    const Index a = m.node(i, 0), b = m.node(i, m.nz - 1);
    p.fixed[static_cast<std::size_t>(a)] = 1;
    p.fixed_value[static_cast<std::size_t>(a)] = 0.0;
    p.fixed[static_cast<std::size_t>(b)] = 1;
    p.fixed_value[static_cast<std::size_t>(b)] = V;
  }
  // The LATERAL surface carries the natural condition j . n = 0.  That is the
  // condition a non-emitting free surface must also carry: charge that crossed
  // it would have nowhere to go.
  p.far_field = FarField::Insulated;
  p.check();

  const AxisymSolution sol = solve_axisym(p);
  out.phi = sol.phi;
  out.fem_residual = sol.residual_inf;

  // The current through the driven end, from the OPERATOR reaction -- the same
  // route by which the electrostatic problem gets a conductor charge, and far
  // more accurate than differentiating phi and integrating the result.
  {
    std::vector<char> mask(static_cast<std::size_t>(m.n_nodes()), 0);
    for (Index i = 0; i < m.nr; ++i) mask[static_cast<std::size_t>(m.node(i, m.nz - 1))] = 1;
    out.current = std::abs(charge_of(sol, mask));
  }
  out.current_closed_form = sigma * pi * R * R * std::abs(V) / L;
  out.current_error = (out.current_closed_form > 0.0)
                          ? std::abs(out.current - out.current_closed_form) /
                                out.current_closed_form
                          : kNaN;
  out.resistance = (out.current > 0.0) ? std::abs(V) / out.current : kNaN;
  out.resistance_closed_form = L / (sigma * pi * R * R);

  // The exact solution is linear in z; check it node by node.
  {
    Real worst = 0.0;
    for (Index j = 0; j < m.nz; ++j)
      for (Index i = 0; i < m.nr; ++i) {
        const Real want = V * m.at(i, j).z / L;
        const Real got = out.phi[static_cast<std::size_t>(m.node(i, j))];
        worst = std::max(worst, std::abs(got - want));
      }
    out.max_potential_error = (std::abs(V) > 0.0) ? worst / std::abs(V) : worst;
  }

  // Lateral leakage: the radial current density on the outer column, against
  // the mean axial one.  In the exact solution it is identically zero.
  {
    const Real j_axial = sigma * std::abs(V) / L;
    Real worst = 0.0;
    for (Index j = 1; j + 1 < m.nz; ++j) {
      const Real r0 = m.at(m.nr - 2, j).r, r1 = m.at(m.nr - 1, j).r;
      const Real p0 = out.phi[static_cast<std::size_t>(m.node(m.nr - 2, j))];
      const Real p1 = out.phi[static_cast<std::size_t>(m.node(m.nr - 1, j))];
      worst = std::max(worst, sigma * std::abs((p1 - p0) / (r1 - r0)));
    }
    out.lateral_leakage = (j_axial > 0.0) ? worst / j_axial : kNaN;
  }
  return out;
}

}  // namespace es
