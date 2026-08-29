#include "es/bem.hpp"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

#include "es/constants.hpp"
#include "es/elliptic.hpp"

namespace es {
namespace {

using constants::eps0;
using constants::pi;

// --- Gauss-Legendre nodes/weights on [-1,1], generated once on first use ----
struct GaussRule {
  std::vector<Real> x, w;
};

GaussRule make_gauss(int n) {
  GaussRule g;
  g.x.resize(static_cast<std::size_t>(n));
  g.w.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    // Chebyshev-like initial guess, then Newton on the Legendre polynomial.
    Real z = std::cos(pi * (static_cast<Real>(i) + 0.75) / (static_cast<Real>(n) + 0.5));
    Real pp = 0.0;
    for (int it = 0; it < 100; ++it) {
      Real p0 = 1.0, p1 = 0.0;
      for (int j = 0; j < n; ++j) {
        const Real p2 = p1;
        p1 = p0;
        p0 = ((2.0 * j + 1.0) * z * p1 - static_cast<Real>(j) * p2) / (j + 1.0);
      }
      pp = static_cast<Real>(n) * (z * p0 - p1) / (z * z - 1.0);
      const Real dz = p0 / pp;
      z -= dz;
      if (std::abs(dz) < 1e-16) break;
    }
    g.x[static_cast<std::size_t>(i)] = z;
    g.w[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - z * z) * pp * pp);
  }
  return g;
}

const GaussRule& gauss(int n) {
  static const GaussRule g4 = make_gauss(4);
  static const GaussRule g8 = make_gauss(8);
  static const GaussRule g12 = make_gauss(12);
  if (n <= 4) return g4;
  if (n <= 8) return g8;
  return g12;
}

// Innermost panel size as a fraction of the element length.  Set so that the
// complementary parameter mc ~ (eps/2r)^2 never drops below ~1e-13, keeping the
// AGM well away from its m -> 1 limit, while the log singularity contained in
// the untreated core contributes < 1e-6 relative.
constexpr Real kMinPanelFrac = 1e-6;
constexpr int kMaxPanels = 24;

/// Panel-refined quadrature over one element.  Calls f(x_quad, weight); the sum
/// of weight*integrand approximates \int_element (...) ds'.
template <class F>
void integrate_panels(const Element& el, Vec2 x, bool self, F&& f) {
  const Real L = el.len;
  if (!(L > 0.0)) return;
  const Vec2 d = el.b - el.a;

  Real tstar, dist;
  if (self) {
    tstar = 0.5;
    dist = 0.0;
  } else {
    tstar = std::clamp(dot(x - el.a, d) / (L * L), 0.0, 1.0);
    dist = norm(x - (el.a + tstar * d));
  }

  const Real ratio = dist / L;
  if (ratio > 5.0) {  // far field: a low-order rule on the whole element
    const GaussRule& g = gauss(4);
    for (std::size_t k = 0; k < g.x.size(); ++k) {
      const Real t = 0.5 * (g.x[k] + 1.0);
      f(el.a + t * d, 0.5 * g.w[k] * L);
    }
    return;
  }
  if (ratio > 1.0) {  // moderately near: higher order, still one panel
    const GaussRule& g = gauss(8);
    for (std::size_t k = 0; k < g.x.size(); ++k) {
      const Real t = 0.5 * (g.x[k] + 1.0);
      f(el.a + t * d, 0.5 * g.w[k] * L);
    }
    return;
  }

  // Near or singular: geometric panels marching away from the closest point in
  // both directions, so the log singularity is resolved to machine-ish accuracy.
  const GaussRule& g = gauss(8);
  const Real eps = std::max(dist, kMinPanelFrac * L);

  auto do_side = [&](Real t_from, Real t_to) {
    const Real span = std::abs(t_to - t_from) * L;
    if (span <= 0.0) return;
    const Real sign = (t_to > t_from) ? 1.0 : -1.0;
    Real s0 = 0.0;
    Real h = eps;
    for (int p = 0; p < kMaxPanels && s0 < span; ++p) {
      const Real s1 = std::min(span, s0 + h);
      const Real ta = t_from + sign * s0 / L;
      const Real tb = t_from + sign * s1 / L;
      const Real half = 0.5 * (tb - ta);
      const Real mid = 0.5 * (ta + tb);
      for (std::size_t k = 0; k < g.x.size(); ++k) {
        const Real t = mid + half * g.x[k];
        f(el.a + t * d, std::abs(half) * g.w[k] * L);
      }
      s0 = s1;
      h *= 2.0;
      if (p == kMaxPanels - 2) h = span;  // make sure the last panel closes
    }
  };

  do_side(tstar, 1.0);
  do_side(tstar, 0.0);
}

}  // namespace

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

Real kernel_G(Vec2 x, Vec2 xp) {
  const Real dz = x.z - xp.z;
  const Real sr = x.r + xp.r;
  const Real S2 = sr * sr + dz * dz;
  if (!(S2 > 0.0) || xp.r <= 0.0) return 0.0;
  const Real dr = x.r - xp.r;
  const Real mc = (dr * dr + dz * dz) / S2;  // = 1 - m, formed without cancellation
  Real K, E;
  ellipKE_mc(mc, K, E);
  return xp.r * K / (pi * eps0 * std::sqrt(S2));
}

Vec2 kernel_gradG(Vec2 x, Vec2 xp) {
  const Real dz = x.z - xp.z;
  const Real sr = x.r + xp.r;
  const Real S2 = sr * sr + dz * dz;
  if (!(S2 > 0.0) || xp.r <= 0.0) return {0.0, 0.0};
  const Real dr = x.r - xp.r;
  const Real mc = (dr * dr + dz * dz) / S2;
  const Real m = 1.0 - mc;
  const Real S = std::sqrt(S2);

  Real K, E;
  ellipKE_mc(mc, K, E);

  Real dKdm;
  if (m < 1e-6) {
    dKdm = 0.5 * pi * (0.25 + (9.0 / 32.0) * m + (75.0 / 256.0) * m * m);
  } else {
    dKdm = (E - mc * K) / (2.0 * m * mc);
  }

  const Real S4 = S2 * S2;
  const Real dm_dr = 4.0 * xp.r / S2 - 8.0 * x.r * xp.r * sr / S4;
  const Real dm_dz = -8.0 * x.r * xp.r * dz / S4;

  const Real C = xp.r / (pi * eps0);
  const Real gr = C * (dKdm * dm_dr / S - K * sr / (S2 * S));
  const Real gz = C * (dKdm * dm_dz / S - K * dz / (S2 * S));
  return {gr, gz};
}

Real integrate_G(const Element& el, Vec2 x, bool self) {
  Real acc = 0.0;
  integrate_panels(el, x, self, [&](Vec2 xp, Real w) { acc += w * kernel_G(x, xp); });
  return acc;
}

Vec2 integrate_gradG(const Element& el, Vec2 x) {
  Vec2 acc{0.0, 0.0};
  integrate_panels(el, x, false, [&](Vec2 xp, Real w) { acc += w * kernel_gradG(x, xp); });
  return acc;
}

// ---------------------------------------------------------------------------
// Solver
// ---------------------------------------------------------------------------

Electrode electrode_of(Tag t) {
  switch (t) {
    case Tag::Emitter:
    case Tag::FreeSurface:
      return Electrode::Emitter;
    case Tag::Extractor:
      return Electrode::Extractor;
    default:
      return Electrode::Collector;
  }
}

void BemSolver::assemble() {
  const Index n = mesh_.size();
  if (n == 0) throw std::runtime_error("BemSolver::assemble: empty mesh");
  lu_ = Matrix(n, n);

#ifdef ES_HAVE_OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
  for (Index i = 0; i < n; ++i) {
    const Vec2 xi = mesh_.elems[static_cast<std::size_t>(i)].mid;
    for (Index j = 0; j < n; ++j) {
      lu_(i, j) = integrate_G(mesh_.elems[static_cast<std::size_t>(j)], xi, i == j);
    }
  }

  if (!lu_factor(lu_, piv_)) throw std::runtime_error("BemSolver: singular influence matrix");
  factored_ = true;
  basis_.clear();
}

void BemSolver::solve_basis() {
  if (!factored_) assemble();
  const Index n = mesh_.size();
  const int ne = static_cast<int>(Electrode::Count);
  basis_.assign(static_cast<std::size_t>(ne) + 1, {});

  for (int k = 0; k < ne; ++k) {
    std::vector<Real> rhs(static_cast<std::size_t>(n), 0.0);
    for (Index i = 0; i < n; ++i) {
      if (electrode_of(mesh_.elems[static_cast<std::size_t>(i)].tag) == static_cast<Electrode>(k))
        rhs[static_cast<std::size_t>(i)] = 1.0;
    }
    lu_solve(lu_, piv_, rhs);
    basis_[static_cast<std::size_t>(k)] = std::move(rhs);
  }

  // Extra basis vector for the externally imposed potential (space charge).
  std::vector<Real> rhs(static_cast<std::size_t>(n), 0.0);
  if (!phi_ext_.empty()) {
    for (Index i = 0; i < n; ++i) rhs[static_cast<std::size_t>(i)] = -phi_ext_[static_cast<std::size_t>(i)];
    lu_solve(lu_, piv_, rhs);
  }
  basis_[static_cast<std::size_t>(ne)] = std::move(rhs);
}

std::vector<Real> BemSolver::sigma_for(const std::array<Real, 3>& V) const {
  if (basis_.empty()) throw std::runtime_error("BemSolver: call solve_basis() first");
  const Index n = mesh_.size();
  std::vector<Real> s(static_cast<std::size_t>(n), 0.0);
  for (int k = 0; k < 3; ++k) {
    if (V[static_cast<std::size_t>(k)] == 0.0) continue;
    const std::vector<Real>& b = basis_[static_cast<std::size_t>(k)];
    for (Index i = 0; i < n; ++i)
      s[static_cast<std::size_t>(i)] += V[static_cast<std::size_t>(k)] * b[static_cast<std::size_t>(i)];
  }
  const std::vector<Real>& be = basis_[3];
  for (Index i = 0; i < n; ++i) s[static_cast<std::size_t>(i)] += be[static_cast<std::size_t>(i)];
  return s;
}

void BemSolver::solve(const std::array<Real, 3>& V) {
  if (basis_.empty()) solve_basis();
  sigma_ = sigma_for(V);
  applied_ = V;
  // Keep the mesh in sync with what was actually applied, so that residual
  // checks, CSV dumps and the meniscus solver all see one consistent state
  // instead of whatever default the geometry builder happened to carry.
  for (Element& e : mesh_.elems)
    e.potential = V[static_cast<std::size_t>(electrode_of(e.tag))];
}

void BemSolver::set_external_potential(std::vector<Real> phi_ext) {
  phi_ext_ = std::move(phi_ext);
  if (factored_) solve_basis();
}

void BemSolver::clear_external_potential() {
  phi_ext_.clear();
  if (factored_) solve_basis();
}

Real BemSolver::En(Index i) const {
  return sigma_.empty() ? 0.0 : sigma_[static_cast<std::size_t>(i)] / eps0;
}

Real BemSolver::charge_on(Tag t) const {
  Real q = 0.0;
  for (Index i = 0; i < mesh_.size(); ++i)
    if (mesh_.elems[static_cast<std::size_t>(i)].tag == t)
      q += sigma_[static_cast<std::size_t>(i)] * mesh_.elems[static_cast<std::size_t>(i)].area;
  return q;
}

Real BemSolver::total_charge() const {
  Real q = 0.0;
  for (Index i = 0; i < mesh_.size(); ++i)
    q += sigma_[static_cast<std::size_t>(i)] * mesh_.elems[static_cast<std::size_t>(i)].area;
  return q;
}

Real BemSolver::capacitance(Electrode e) const {
  Real q = 0.0, v = 0.0;
  for (Index i = 0; i < mesh_.size(); ++i) {
    const Element& el = mesh_.elems[static_cast<std::size_t>(i)];
    if (electrode_of(el.tag) != e) continue;
    q += sigma_[static_cast<std::size_t>(i)] * el.area;
    v = el.potential;
  }
  return (v != 0.0) ? q / v : 0.0;
}

Real BemSolver::potential_at(Vec2 x) const {
  Real v = 0.0;
  for (Index i = 0; i < mesh_.size(); ++i)
    v += sigma_[static_cast<std::size_t>(i)] *
         integrate_G(mesh_.elems[static_cast<std::size_t>(i)], x, false);
  return v;
}

Vec2 BemSolver::field_at(Vec2 x) const {
  Vec2 g{0.0, 0.0};
  for (Index i = 0; i < mesh_.size(); ++i)
    g += sigma_[static_cast<std::size_t>(i)] *
         integrate_gradG(mesh_.elems[static_cast<std::size_t>(i)], x);
  return -1.0 * g;  // E = -grad V
}

Real BemSolver::peak_field(Tag t, Index* which) const {
  Real best = 0.0;
  Index arg = -1;
  for (Index i = 0; i < mesh_.size(); ++i) {
    if (mesh_.elems[static_cast<std::size_t>(i)].tag != t) continue;
    const Real e = std::abs(En(i));
    if (e > best) { best = e; arg = i; }
  }
  if (which) *which = arg;
  return best;
}

Real BemSolver::peak_emitter_field(Index* which) const {
  Real best = 0.0;
  Index arg = -1;
  for (Index i = 0; i < mesh_.size(); ++i) {
    if (electrode_of(mesh_.elems[static_cast<std::size_t>(i)].tag) != Electrode::Emitter) continue;
    const Real e = std::abs(En(i));
    if (e > best) { best = e; arg = i; }
  }
  if (which) *which = arg;
  return best;
}

void BemSolver::write_surface_csv(const std::string& path,
                                  const std::string& header) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  if (!header.empty()) std::fputs(header.c_str(), f);
  std::fprintf(f, "i,r,z,arclen_hint,sigma,En,p_electric,tag\n");
  Real s = 0.0;
  for (Index i = 0; i < mesh_.size(); ++i) {
    const Element& e = mesh_.elems[static_cast<std::size_t>(i)];
    s += e.len;
    const Real En_i = En(i);
    std::fprintf(f, "%td,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%s\n", i, e.mid.r, e.mid.z, s,
                 sigma_.empty() ? 0.0 : sigma_[static_cast<std::size_t>(i)], En_i,
                 0.5 * eps0 * En_i * En_i, tag_name(e.tag));
  }
  std::fclose(f);
}

}  // namespace es
