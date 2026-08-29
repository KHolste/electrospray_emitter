#include "es/load_projection.hpp"

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
}

// ===========================================================================
// ProjectedLoad -- the load handed to the capillary solver
// ===========================================================================

ProjectedLoad::ProjectedLoad()
    : p(static_cast<std::size_t>(kBins), 0.0), area(static_cast<std::size_t>(kBins), 0.0) {}

ProjectedLoad ProjectedLoad::from(const MaxwellLoad& L) {
  ProjectedLoad t;
  std::vector<Real> force(static_cast<std::size_t>(kBins), 0.0);
  const Real db = 1.0 / static_cast<Real>(kBins);
  for (std::size_t k = 0; k < L.seg_force.size(); ++k) {
    const Real t0 = L.seg_tau0[k], t1 = L.seg_tau1[k];
    const Real span = std::max(t1 - t0, 1e-300);
    const int b0 = std::max(0, std::min(kBins - 1, static_cast<int>(t0 / db)));
    const int b1 = std::max(0, std::min(kBins - 1, static_cast<int>(t1 / db)));
    for (int b = b0; b <= b1; ++b) {
      const Real lo = std::max(t0, static_cast<Real>(b) * db);
      const Real hi = std::min(t1, static_cast<Real>(b + 1) * db);
      const Real f = std::max(0.0, hi - lo) / span;
      force[static_cast<std::size_t>(b)] += f * L.seg_force[k];
      t.area[static_cast<std::size_t>(b)] += f * L.seg_area[k];
    }
  }
  for (int b = 0; b < kBins; ++b) {
    const std::size_t i = static_cast<std::size_t>(b);
    t.p[i] = t.area[i] > 0.0 ? force[i] / t.area[i] : 0.0;
  }
  return t;
}

/// Monotone cubic (Fritsch-Carlson) slopes for a set of cumulative values on a
/// uniform grid.  Monotone because the cumulative quantities are.
static std::vector<Real> monotone_slopes(const std::vector<Real>& y, Real h) {
  const std::size_t n = y.size();
  std::vector<Real> d(n, 0.0), delta(n - 1, 0.0);
  for (std::size_t k = 0; k + 1 < n; ++k) delta[k] = (y[k + 1] - y[k]) / h;
  d[0] = delta[0];
  d[n - 1] = delta[n - 2];
  for (std::size_t k = 1; k + 1 < n; ++k) {
    if (delta[k - 1] * delta[k] <= 0.0) {
      d[k] = 0.0;
    } else {
      const Real w1 = 2.0 * h + h, w2 = h + 2.0 * h;
      d[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
    }
  }
  return d;
}

void ProjectedLoad::build_cumulative() const {
  const Real h = 1.0 / static_cast<Real>(kBins);
  cum_force_.assign(static_cast<std::size_t>(kBins) + 1, 0.0);
  cum_area_.assign(static_cast<std::size_t>(kBins) + 1, 0.0);
  for (int b = 0; b < kBins; ++b) {
    const std::size_t i = static_cast<std::size_t>(b);
    cum_force_[i + 1] = cum_force_[i] + p[i] * area[i];
    cum_area_[i + 1] = cum_area_[i] + area[i];
  }
  slope_force_ = monotone_slopes(cum_force_, h);
  slope_area_ = monotone_slopes(cum_area_, h);
}

/// Derivative of the cubic Hermite through (y[b], y[b+1]) with slopes d, at the
/// local coordinate u in [0,1].  At u = 0 it returns d[b] and at u = 1 it
/// returns d[b+1]; neighbouring intervals share those slopes, so the derivative
/// is continuous across every bin boundary.  That is the continuity claim, and
/// it is a property of this formula, not an assertion about it.
Real ProjectedLoad::slope_of(const std::vector<Real>& y, const std::vector<Real>& d, int b,
                             Real u) const {
  const Real h = 1.0 / static_cast<Real>(kBins);
  const std::size_t i = static_cast<std::size_t>(b);
  const Real y0 = y[i], y1 = y[i + 1], d0 = d[i] * h, d1 = d[i + 1] * h;
  const Real g00 = 6 * u * u - 6 * u, g10 = 3 * u * u - 4 * u + 1;
  const Real g01 = -6 * u * u + 6 * u, g11 = 3 * u * u - 2 * u;
  return (g00 * y0 + g10 * d0 + g01 * y1 + g11 * d1) / h;
}

Real ProjectedLoad::at(Real tau) const {
  const Real t = std::min(1.0, std::max(0.0, tau));
  const Real h = 1.0 / static_cast<Real>(kBins);
  if (cum_force_.empty()) build_cumulative();
  int b = static_cast<int>(t / h);
  b = std::max(0, std::min(kBins - 1, b));
  const Real u = (t - static_cast<Real>(b) * h) / h;
  const Real dA = slope_of(cum_area_, slope_area_, b, u);
  if (!(dA > 0.0)) return p[static_cast<std::size_t>(b)];
  return slope_of(cum_force_, slope_force_, b, u) / dA;
}

Real ProjectedLoad::bin_pressure_at(Real tau) const {
  const Real t = std::min(1.0, std::max(0.0, tau));
  const Real h = 1.0 / static_cast<Real>(kBins);
  int b = static_cast<int>(t / h);
  b = std::max(0, std::min(kBins - 1, b));
  return p[static_cast<std::size_t>(b)];
}

Real ProjectedLoad::integrated_force() const {
  Real f = 0.0;
  for (std::size_t k = 0; k < p.size(); ++k) f += p[k] * area[k];
  return f;
}

Real ProjectedLoad::total_area() const {
  Real s = 0.0;
  for (Real v : area) s += v;
  return s;
}

Index ProjectedLoad::empty_bins() const {
  Index n = 0;
  for (Real v : area)
    if (!(v > 0.0)) ++n;
  return n;
}

Real ProjectedLoad::difference(const ProjectedLoad& a, const ProjectedLoad& b) {
  Real d = 0.0;
  for (std::size_t k = 0; k < a.p.size(); ++k) d = std::max(d, std::abs(a.p[k] - b.p[k]));
  return d;
}

ProjectedLoad ProjectedLoad::blend(const ProjectedLoad& old_load, const ProjectedLoad& fresh,
                                   Real w) {
  ProjectedLoad t;
  for (std::size_t k = 0; k < t.p.size(); ++k) {
    t.p[k] = (1.0 - w) * old_load.p[k] + w * fresh.p[k];
    t.area[k] = fresh.area[k];
  }
  return t;
}

bool ProjectedLoad::is_zero() const {
  for (Real v : p)
    if (v != 0.0) return false;
  return true;
}

const std::vector<Real>& ProjectedLoad::cumulative_force() const {
  if (cum_force_.empty()) build_cumulative();
  return cum_force_;
}

const std::vector<Real>& ProjectedLoad::cumulative_area() const {
  if (cum_area_.empty()) build_cumulative();
  return cum_area_;
}

// ===========================================================================
// Manufactured loads
// ===========================================================================

std::vector<Real> uniform_radius_nodes(Real contact_radius, Index n_segments) {
  if (n_segments < 1) throw std::runtime_error("uniform_radius_nodes: n_segments < 1");
  std::vector<Real> r(static_cast<std::size_t>(n_segments) + 1, 0.0);
  for (Index k = 0; k <= n_segments; ++k)
    r[static_cast<std::size_t>(k)] =
        contact_radius * static_cast<Real>(k) / static_cast<Real>(n_segments);
  return r;
}

MaxwellLoad manufactured_load(const FreeSurface& fs, const std::vector<Real>& node_r,
                              const PrescribedPressure& p_of_d, Real gamma_over_a) {
  if (node_r.size() < 2) throw std::runtime_error("manufactured_load: zu wenige Knoten");
  MaxwellLoad out;
  out.gamma_over_a = gamma_over_a;
  const Real L = fs.arclength;
  for (std::size_t k = 0; k < node_r.size(); ++k) {
    const Real r = node_r[k];
    const Real s = fs.s_at_radius(r);
    const Real d = std::max(0.0, L - s);
    out.node_r.push_back(r);
    out.node_z.push_back(fs.z_at_radius(r));
    out.node_s.push_back(s);
    out.node_tau.push_back(L > 0.0 ? s / L : 0.0);
    out.node_d_edge.push_back(d);
  }
  // The edge node gets the local mean over the last half segment instead of a
  // point value.  See the header: it is a stated property of the manufactured
  // INPUT and mirrors what a recovered FEM field does there.  Every other node
  // gets the exact prescribed value.
  const std::size_t n = out.node_r.size();
  const Real h_edge = out.node_d_edge[n - 2] - out.node_d_edge[n - 1];
  for (std::size_t k = 0; k < n; ++k) {
    Real value;
    if (k + 1 == n) {
      const int m = 4096;
      Real acc = 0.0;
      for (int q = 0; q < m; ++q) {
        const Real d = 0.5 * h_edge * (static_cast<Real>(q) + 0.5) / static_cast<Real>(m);
        acc += p_of_d(d);
      }
      value = acc / static_cast<Real>(m);
    } else {
      value = p_of_d(out.node_d_edge[k]);
    }
    out.node_pM.push_back(value);
    out.node_En.push_back(std::sqrt(std::max(0.0, 2.0 * value / eps0)));
    out.node_tangential_fraction.push_back(0.0);
  }
  assemble_load_segments(out, fs);
  return out;
}

Real flat_disc_power_law_force(Real C, Real beta, Real a) {
  if (!(beta > -1.0)) return kNaN;
  return 2.0 * pi * C * std::pow(a, 2.0 + beta) / ((1.0 + beta) * (2.0 + beta));
}

Real flat_disc_smooth_force(Real p0, Real a) { return 1.5 * pi * p0 * a * a; }

Real flat_disc_power_law_force_beyond(Real C, Real beta, Real a, Real d0) {
  // F(d0) = int_0^{a-d0} C (a-r)^beta 2 pi r dr, substituting u = a - r:
  //       = 2 pi C int_{d0}^{a} u^beta (a - u) du
  //       = 2 pi C [ a (a^{1+b} - d0^{1+b})/(1+b) - (a^{2+b} - d0^{2+b})/(2+b) ].
  if (!(beta > -1.0)) return kNaN;
  const Real b1 = 1.0 + beta, b2 = 2.0 + beta;
  return 2.0 * pi * C *
         (a * (std::pow(a, b1) - std::pow(d0, b1)) / b1 -
          (std::pow(a, b2) - std::pow(d0, b2)) / b2);
}

Real exclusion_halving_limit(Real beta, Real d0_over_a) {
  // Independent of C and of a: form the ratio with a = 1, C = 1.
  const Real F_total = flat_disc_power_law_force(1.0, beta, 1.0);
  const Real F0 = flat_disc_power_law_force_beyond(1.0, beta, 1.0, d0_over_a);
  const Real F1 = flat_disc_power_law_force_beyond(1.0, beta, 1.0, 0.5 * d0_over_a);
  if (!(std::abs(F1) > 0.0)) return kNaN;
  // The gate forms |F(d0/2) - F(d0)| / F(d0/2), the same normalisation
  // run_edge_gate() uses for measured_exclusion_change.
  (void)F_total;
  return std::abs(F1 - F0) / std::abs(F1);
}

// ===========================================================================
// The audit
// ===========================================================================

LoadProjectionAudit audit_projection(const MaxwellLoad& raw, const FreeSurface& fs,
                                     const std::string& tag, Real analytic_force) {
  LoadProjectionAudit A;
  A.tag = tag;
  A.n_nodes = static_cast<Index>(raw.node_r.size());
  A.n_segments = static_cast<Index>(raw.seg_force.size());
  A.analytic_force = analytic_force;
  A.segment_force = raw.total_force;

  const ProjectedLoad pl = ProjectedLoad::from(raw);
  A.bin_force = pl.integrated_force();
  A.bin_area = pl.total_area();
  A.empty_bins = pl.empty_bins();

  A.segment_area = 0.0;
  for (Real v : raw.seg_area) A.segment_area += v;
  A.surface_area = fs.revolved_area();

  // --- force against the reconstructed measure -------------------------------
  //
  // int p A' dtau = int G' dtau = G(1) - G(0), exactly, by the fundamental
  // theorem applied to the cubic Hermite.  Evaluating the cumulative
  // interpolant at the ends is therefore the integral itself and no quadrature
  // is needed -- the number below IS the reconstructed-measure force.
  const std::vector<Real>& G = pl.cumulative_force();
  A.handed_force_reconstructed = G.back() - G.front();

  // --- force against the TRUE surface element --------------------------------
  //
  // int p(tau) 2 pi r(s) ds with s = tau * L, by composite Gauss-Legendre with
  // two points per bin.  This is what the capillary solver actually integrates,
  // and it is reported next to the exact one instead of in place of it.
  {
    const Real L = fs.arclength;
    const int nb = ProjectedLoad::kBins, nq = 4;
    const Real h = 1.0 / static_cast<Real>(nb);
    // 4-point Gauss-Legendre on [0,1]
    const Real xg[4] = {0.0694318442029737, 0.3300094782075719, 0.6699905217924281,
                        0.9305681557970263};
    const Real wg[4] = {0.1739274225687269, 0.3260725774312731, 0.3260725774312731,
                        0.1739274225687269};
    Real acc = 0.0;
    for (int b = 0; b < nb; ++b) {
      for (int q = 0; q < nq; ++q) {
        const Real tau = (static_cast<Real>(b) + xg[q]) * h;
        const Real s = tau * L;
        // r(s) from the surface: invert by walking the stored arclength nodes.
        Real r = 0.0;
        {
          const std::size_t m = fs.nodes.size();
          const Real ds = (m > 1) ? L / static_cast<Real>(m - 1) : 0.0;
          const Real x = (ds > 0.0) ? s / ds : 0.0;
          std::size_t k = static_cast<std::size_t>(x);
          if (k + 1 >= m) k = m - 2;
          const Real t = x - static_cast<Real>(k);
          r = (1.0 - t) * fs.nodes[k].r + t * fs.nodes[k + 1].r;
        }
        acc += wg[q] * h * pl.at(tau) * 2.0 * pi * r * L;
      }
    }
    A.handed_force_true = acc;
  }

  auto rel = [](Real x, Real ref) {
    return (std::abs(ref) > 0.0) ? std::abs(x - ref) / std::abs(ref) : kNaN;
  };
  A.error_segment_vs_analytic = rel(A.segment_force, analytic_force);
  A.error_bin_vs_segment = rel(A.bin_force, A.segment_force);
  A.error_handed_reconstructed = rel(A.handed_force_reconstructed, A.segment_force);
  A.error_handed_true = rel(A.handed_force_true, A.segment_force);

  // --- continuity across the bin boundaries ---------------------------------
  {
    const Real h = 1.0 / static_cast<Real>(ProjectedLoad::kBins);
    Real lo = std::numeric_limits<Real>::max(), hi = -lo;
    for (int pass = 0; pass < 2; ++pass) {
      const Real d = LoadProjectionAudit::kProbeOffset * h * (pass == 0 ? 1.0 : 0.1);
      Real jh = 0.0, js = 0.0;
      for (int b = 1; b < ProjectedLoad::kBins; ++b) {
        const Real tb = static_cast<Real>(b) * h;
        const Real a_left = pl.at(tb - d), a_right = pl.at(tb + d);
        const Real s_left = pl.bin_pressure_at(tb - d), s_right = pl.bin_pressure_at(tb + d);
        jh = std::max(jh, std::abs(a_right - a_left));
        js = std::max(js, std::abs(s_right - s_left));
        if (pass == 0) {
          lo = std::min({lo, a_left, a_right});
          hi = std::max({hi, a_left, a_right});
        }
      }
      if (pass == 0) {
        A.handed_max_jump = jh;
        A.staircase_max_jump = js;
      } else {
        A.handed_max_jump_tenth = jh;
        A.staircase_max_jump_tenth = js;
      }
    }
    A.load_span = hi - lo;
    A.handed_jump_ratio = (A.load_span > 0.0) ? A.handed_max_jump / A.load_span : kNaN;
    A.staircase_jump_ratio = (A.load_span > 0.0) ? A.staircase_max_jump / A.load_span : kNaN;
    A.handed_jump_decay =
        (A.handed_max_jump > 0.0) ? A.handed_max_jump_tenth / A.handed_max_jump : kNaN;
    A.staircase_jump_decay =
        (A.staircase_max_jump > 0.0) ? A.staircase_max_jump_tenth / A.staircase_max_jump : kNaN;
  }

  // --- edge and coverage ----------------------------------------------------
  if (!raw.seg_tau0.empty()) {
    A.tau_first = raw.seg_tau0.front();
    A.tau_last = raw.seg_tau1.back();
  }
  A.area_gap = rel(A.bin_area, A.segment_area);
  {
    const std::size_t last = static_cast<std::size_t>(ProjectedLoad::kBins) - 1;
    const Real f_last = pl.p[last] * pl.area[last];
    A.last_bin_force_fraction =
        (std::abs(A.bin_force) > 0.0) ? f_last / A.bin_force : kNaN;
  }
  for (int b = 0; b < ProjectedLoad::kBins; ++b)
    A.max_bin_pressure = std::max(A.max_bin_pressure, pl.p[static_cast<std::size_t>(b)]);
  {
    const int m = 4096;
    for (int k = 0; k <= m; ++k)
      A.max_handed_pressure =
          std::max(A.max_handed_pressure, pl.at(static_cast<Real>(k) / static_cast<Real>(m)));
  }
  for (Real v : raw.node_pM) A.max_node_pressure = std::max(A.max_node_pressure, v);
  return A;
}

void LoadProjectionAudit::print(std::FILE* out) const {
  std::fprintf(out, "  Lastprojektion '%s': %lld Knoten, %lld Segmente\n", tag.c_str(),
               static_cast<long long>(n_nodes), static_cast<long long>(n_segments));
  std::fprintf(out, "    Kraft analytisch          : %.9e N\n", analytic_force);
  std::fprintf(out, "    Kraft Segmentprojektion   : %.9e N  (rel. %.3e)\n", segment_force,
               error_segment_vs_analytic);
  std::fprintf(out, "    Kraft Bins                : %.9e N  (rel. %.3e)\n", bin_force,
               error_bin_vs_segment);
  std::fprintf(out, "    Kraft uebergebene Last    : %.9e N gegen A' (rel. %.3e)\n",
               handed_force_reconstructed, error_handed_reconstructed);
  std::fprintf(out, "                                %.9e N gegen 2 pi r ds (rel. %.3e)\n",
               handed_force_true, error_handed_true);
  std::fprintf(out, "    Sprung uebergebene Last   : %.3e Pa  (%.3e der Spannweite),"
                    "  Abfall bei delta/10: %.4f\n",
               handed_max_jump, handed_jump_ratio, handed_jump_decay);
  std::fprintf(out, "    Sprung Treppenfunktion    : %.3e Pa  (%.3e der Spannweite),"
                    "  Abfall bei delta/10: %.4f\n",
               staircase_max_jump, staircase_jump_ratio, staircase_jump_decay);
  std::fprintf(out, "    tau-Abdeckung             : %.3e .. %.6f, leere Bins %lld\n", tau_first,
               tau_last, static_cast<long long>(empty_bins));
  std::fprintf(out, "    Kraft im letzten Bin      : %.3e der Gesamtkraft\n",
               last_bin_force_fraction);
  std::fprintf(out, "    Maximum uebergeben/Bin/Knoten: %.4e / %.4e / %.4e Pa\n",
               max_handed_pressure, max_bin_pressure, max_node_pressure);
}

// ===========================================================================
// Richardson
// ===========================================================================

RichardsonEstimate richardson(const std::vector<Real>& values, Real ratio) {
  RichardsonEstimate e;
  e.n_levels = static_cast<int>(values.size());
  e.ratio = ratio;
  e.observed_order = kNaN;
  e.extrapolated = kNaN;
  e.relative_error_finest = kNaN;
  e.last_relative_change = kNaN;
  if (values.size() < 2) {
    e.note = "Weniger als zwei Stufen.";
    return e;
  }
  const std::size_t n = values.size();
  const Real f_fine = values[n - 1], f_mid = values[n - 2];
  e.last_relative_change =
      (std::abs(f_fine) > 0.0) ? std::abs(f_fine - f_mid) / std::abs(f_fine) : kNaN;
  if (values.size() < 3) {
    e.note = "Weniger als drei Stufen: nur die Aenderung ist bekannt.";
    return e;
  }
  const Real f_coarse = values[n - 3];
  const Real d1 = f_mid - f_coarse, d2 = f_fine - f_mid;
  e.monotone = (d1 * d2 > 0.0);
  if (!e.monotone || !(std::abs(d2) > 0.0)) {
    e.note = "Die Folge ist nicht monoton oder steht still: keine Ordnung beobachtbar, "
             "nur die Aenderung zwischen den beiden feinsten Stufen ist bekannt.";
    return e;
  }
  const Real p = std::log(std::abs(d1 / d2)) / std::log(ratio);
  const Real rp = std::pow(ratio, p);
  if (!(std::isfinite(p)) || !(rp > 1.0 + 1e-12)) {
    e.note = "Die beobachtete Ordnung ist nicht auswertbar (r^p <= 1).";
    return e;
  }
  e.observed_order = p;
  e.extrapolated = f_fine + d2 / (rp - 1.0);
  e.relative_error_finest = (std::abs(e.extrapolated) > 0.0)
                                ? std::abs(e.extrapolated - f_fine) / std::abs(e.extrapolated)
                                : kNaN;
  e.usable = std::isfinite(e.relative_error_finest);
  e.note = e.usable ? "Richardson auf den drei feinsten Stufen." : "Extrapolation nicht endlich.";
  return e;
}

const char* to_string(DiscretizationVerdict v) {
  switch (v) {
    case DiscretizationVerdict::NotAttempted: return "NotAttempted";
    case DiscretizationVerdict::Converged: return "Converged";
    case DiscretizationVerdict::DiscretizationNotConverged: return "DiscretizationNotConverged";
    case DiscretizationVerdict::NotInAsymptoticRange: return "NotInAsymptoticRange";
    case DiscretizationVerdict::InsufficientLevels: return "InsufficientLevels";
  }
  return "?";
}

const char* explain(DiscretizationVerdict v) {
  switch (v) {
    case DiscretizationVerdict::NotAttempted:
      return "Es wurde keine Netzstudie gerechnet.";
    case DiscretizationVerdict::Converged:
      return "Der geschaetzte Diskretisierungsfehler der feinsten Stufe liegt unter dem vorab "
             "festgelegten Ziel von 1 %.";
    case DiscretizationVerdict::DiscretizationNotConverged:
      return "Der geschaetzte Diskretisierungsfehler der feinsten Stufe verfehlt das vorab "
             "festgelegte Ziel von 1 %.  Das Ergebnis ist nur qualitativ zu lesen.";
    case DiscretizationVerdict::NotInAsymptoticRange:
      return "Die Folge ueber die Netzstufen ist nicht monoton; eine Ordnung ist nicht "
             "beobachtbar und ein Fehler daher nicht schaetzbar.  Bekannt ist nur die Aenderung "
             "zwischen den beiden feinsten Stufen, und das ist kein Fehler.";
    case DiscretizationVerdict::InsufficientLevels:
      return "Weniger als drei Netzstufen: eine Extrapolation ist nicht moeglich.";
  }
  return "?";
}

DiscretizationVerdict verdict_of(const RichardsonEstimate& e) {
  if (e.n_levels < 3) return DiscretizationVerdict::InsufficientLevels;
  if (!e.usable) return DiscretizationVerdict::NotInAsymptoticRange;
  return (e.relative_error_finest < kDiscretizationTarget)
             ? DiscretizationVerdict::Converged
             : DiscretizationVerdict::DiscretizationNotConverged;
}

}  // namespace es
