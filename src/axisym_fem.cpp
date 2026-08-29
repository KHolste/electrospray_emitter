#include "es/axisym_fem.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "es/constants.hpp"
#include "es/linalg.hpp"

namespace es {

using constants::eps0;
using constants::pi;

const char* to_string(FarField f) {
  switch (f) {
    case FarField::Asymptotic: return "asymptotic";
    case FarField::Grounded: return "grounded";
    case FarField::Insulated: return "insulated";
  }
  return "asymptotic";
}

// ---------------------------------------------------------------------------
// Mesh
// ---------------------------------------------------------------------------

std::array<Index, 4> QuadMesh::cell_nodes(Index i, Index j) const {
  return {node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)};
}

Real QuadMesh::cell_meridian_area(Index i, Index j) const {
  const auto c = cell_nodes(i, j);
  Real a = 0.0;
  for (int k = 0; k < 4; ++k) {
    const Vec2& p = nodes[static_cast<std::size_t>(c[k])];
    const Vec2& q = nodes[static_cast<std::size_t>(c[(k + 1) % 4])];
    a += p.r * q.z - q.r * p.z;
  }
  return 0.5 * a;
}

Real QuadMesh::cell_revolved_volume(Index i, Index j) const {
  // int_A 2 pi r dA = oint pi r^2 dz, exact over straight edges.
  const auto c = cell_nodes(i, j);
  Real v = 0.0;
  for (int k = 0; k < 4; ++k) {
    const Vec2& p = nodes[static_cast<std::size_t>(c[k])];
    const Vec2& q = nodes[static_cast<std::size_t>(c[(k + 1) % 4])];
    v += pi * (p.r * p.r + p.r * q.r + q.r * q.r) / 3.0 * (q.z - p.z);
  }
  return v;
}

Real QuadMesh::total_revolved_volume() const {
  Real v = 0.0;
  for (Index j = 0; j + 1 < nz; ++j)
    for (Index i = 0; i + 1 < nr; ++i) v += cell_revolved_volume(i, j);
  return v;
}

bool QuadMesh::validate_level_rows(Real tol_rel) const {
  Real span = 0.0;
  for (Index j = 0; j < nz; ++j) span = std::max(span, std::abs(z_of_row(j)));
  const Real tol = tol_rel * std::max(span, 1e-30);
  for (Index j = 0; j < nz; ++j) {
    const Real z0 = z_of_row(j);
    for (Index i = 0; i < nr; ++i)
      if (std::abs(at(i, j).z - z0) > tol) return false;
  }
  return true;
}

namespace {

// Q1 on the unit square, nodes ordered (0,0) (1,0) (1,1) (0,1).
inline void shape_derivatives(Real xi, Real eta, Real dxi[4], Real deta[4]) {
  dxi[0] = -(1.0 - eta);
  dxi[1] = (1.0 - eta);
  dxi[2] = eta;
  dxi[3] = -eta;
  deta[0] = -(1.0 - xi);
  deta[1] = -xi;
  deta[2] = xi;
  deta[3] = (1.0 - xi);
}

inline void shape(Real xi, Real eta, Real N[4]) {
  N[0] = (1.0 - xi) * (1.0 - eta);
  N[1] = xi * (1.0 - eta);
  N[2] = xi * eta;
  N[3] = (1.0 - xi) * eta;
}

/// 2x2 Gauss on [0,1]^2.
constexpr Real kG0 = 0.5 - 0.5 / 1.7320508075688772935;
constexpr Real kG1 = 0.5 + 0.5 / 1.7320508075688772935;
const Real kGaussXi[4] = {kG0, kG1, kG0, kG1};
const Real kGaussEta[4] = {kG0, kG0, kG1, kG1};
constexpr Real kGaussW = 0.25;

/// Geometry of one Gauss point of a cell: radius, |J|, and the four gradients.
struct GaussPoint {
  Real r{0.0}, detJ{0.0};
  Real dNdr[4]{}, dNdz[4]{};
};

GaussPoint gauss_point(const QuadMesh& m, Index i, Index j, Real xi, Real eta) {
  const auto c = m.cell_nodes(i, j);
  Vec2 p[4];
  for (int k = 0; k < 4; ++k) p[k] = m.nodes[static_cast<std::size_t>(c[k])];

  Real dxi[4], deta[4], N[4];
  shape_derivatives(xi, eta, dxi, deta);
  shape(xi, eta, N);

  Real drdxi = 0, drdeta = 0, dzdxi = 0, dzdeta = 0;
  GaussPoint g;
  for (int k = 0; k < 4; ++k) {
    drdxi += p[k].r * dxi[k];
    drdeta += p[k].r * deta[k];
    dzdxi += p[k].z * dxi[k];
    dzdeta += p[k].z * deta[k];
    g.r += p[k].r * N[k];
  }
  g.detJ = drdxi * dzdeta - drdeta * dzdxi;
  if (!(g.detJ > 0.0))
    throw std::runtime_error("axisym_fem: Zelle mit nicht positiver Jacobi-Determinante");
  for (int k = 0; k < 4; ++k) {
    g.dNdr[k] = (dzdeta * dxi[k] - dzdxi * deta[k]) / g.detJ;
    g.dNdz[k] = (-drdeta * dxi[k] + drdxi * deta[k]) / g.detJ;
  }
  return g;
}

/// Element stiffness  K_ab = \int eps grad Na . grad Nb  2 pi r dr dz  [F].
void element_matrix(const QuadMesh& m, Index i, Index j, Real eps_abs, Real Ke[4][4]) {
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) Ke[a][b] = 0.0;
  for (int q = 0; q < 4; ++q) {
    const GaussPoint g = gauss_point(m, i, j, kGaussXi[q], kGaussEta[q]);
    const Real w = eps_abs * 2.0 * pi * g.r * g.detJ * kGaussW;
    for (int a = 0; a < 4; ++a)
      for (int b = 0; b < 4; ++b)
        Ke[a][b] += w * (g.dNdr[a] * g.dNdr[b] + g.dNdz[a] * g.dNdz[b]);
  }
}

/// Element load vector of a source that is constant within the cell:
/// fe[a] = int N_a s 2 pi r dA, on the same four Gauss points.
void element_source(const QuadMesh& m, Index i, Index j, Real s, Real fe[4]) {
  for (int a = 0; a < 4; ++a) fe[a] = 0.0;
  if (s == 0.0) return;
  for (int q = 0; q < 4; ++q) {
    const GaussPoint g = gauss_point(m, i, j, kGaussXi[q], kGaussEta[q]);
    const Real xi = kGaussXi[q], eta = kGaussEta[q];
    const Real N[4] = {(1 - xi) * (1 - eta), xi * (1 - eta), xi * eta, (1 - xi) * eta};
    const Real w = s * 2.0 * pi * g.r * g.detJ * kGaussW;
    for (int a = 0; a < 4; ++a) fe[a] += w * N[a];
  }
}

// ---------------------------------------------------------------------------
// Symmetric band matrix, lower triangle.  store[i*(b+1) + (i-j)] = A(i,j).
// ---------------------------------------------------------------------------
class BandMatrix {
 public:
  BandMatrix(Index n, Index b) : n_(n), b_(b), a_(static_cast<std::size_t>(n) * (b + 1), 0.0) {}

  std::size_t bytes() const { return a_.size() * sizeof(Real); }

  Real& at(Index i, Index j) {  // requires i >= j and i - j <= b
    return a_[static_cast<std::size_t>(i) * (b_ + 1) + static_cast<std::size_t>(i - j)];
  }
  Real at(Index i, Index j) const {
    return a_[static_cast<std::size_t>(i) * (b_ + 1) + static_cast<std::size_t>(i - j)];
  }
  void add(Index i, Index j, Real v) {
    if (i < j) std::swap(i, j);
    at(i, j) += v;
  }
  Real get(Index i, Index j) const {
    if (i < j) std::swap(i, j);
    return (i - j <= b_) ? at(i, j) : 0.0;
  }
  void set(Index i, Index j, Real v) {
    if (i < j) std::swap(i, j);
    at(i, j) = v;
  }

  Index n() const { return n_; }
  Index bandwidth() const { return b_; }

  /// In-place Cholesky of the band.  False if not positive definite.
  bool factor() {
    for (Index j = 0; j < n_; ++j) {
      Real s = at(j, j);
      for (Index k = std::max<Index>(0, j - b_); k < j; ++k) s -= at(j, k) * at(j, k);
      if (!(s > 0.0)) return false;
      const Real d = std::sqrt(s);
      at(j, j) = d;
      const Index imax = std::min<Index>(n_ - 1, j + b_);
      for (Index i = j + 1; i <= imax; ++i) {
        Real t = at(i, j);
        for (Index k = std::max<Index>(0, i - b_); k < j; ++k) t -= at(i, k) * at(j, k);
        at(i, j) = t / d;
      }
    }
    return true;
  }

  /// Solve with the factor produced by factor(); x is overwritten.
  void solve(std::vector<Real>& x) const {
    for (Index i = 0; i < n_; ++i) {
      Real s = x[static_cast<std::size_t>(i)];
      for (Index k = std::max<Index>(0, i - b_); k < i; ++k)
        s -= at(i, k) * x[static_cast<std::size_t>(k)];
      x[static_cast<std::size_t>(i)] = s / at(i, i);
    }
    for (Index i = n_ - 1; i >= 0; --i) {
      Real s = x[static_cast<std::size_t>(i)];
      const Index kmax = std::min<Index>(n_ - 1, i + b_);
      for (Index k = i + 1; k <= kmax; ++k) s -= at(k, i) * x[static_cast<std::size_t>(k)];
      x[static_cast<std::size_t>(i)] = s / at(i, i);
    }
  }

 private:
  Index n_, b_;
  std::vector<Real> a_;
};

/// Outward unit normal of a boundary segment, oriented away from `inside`.
Vec2 outward_normal(Vec2 a, Vec2 b, Vec2 inside) {
  Vec2 n = normalized(perp(b - a));
  const Vec2 mid = 0.5 * (a + b);
  if (dot(n, mid - inside) < 0.0) n = -1.0 * n;
  return n;
}

}  // namespace

// ---------------------------------------------------------------------------

void AxisymProblem::check() const {
  if (!mesh) throw std::runtime_error("AxisymProblem: kein Netz");
  const QuadMesh& m = *mesh;
  if (m.nr < 2 || m.nz < 2) throw std::runtime_error("AxisymProblem: Netz zu klein");
  if (static_cast<Index>(m.nodes.size()) != m.n_nodes())
    throw std::runtime_error("AxisymProblem: Knotenzahl passt nicht zu (nr, nz)");
  if (static_cast<Index>(eps_r.size()) != m.n_cells() ||
      static_cast<Index>(active.size()) != m.n_cells())
    throw std::runtime_error("AxisymProblem: eps_r/active haben nicht die Zellenzahl");
  if (static_cast<Index>(fixed.size()) != m.n_nodes() ||
      static_cast<Index>(fixed_value.size()) != m.n_nodes())
    throw std::runtime_error("AxisymProblem: fixed/fixed_value haben nicht die Knotenzahl");
  if (require_level_rows && !m.validate_level_rows())
    throw std::runtime_error("AxisymProblem: die Netzzeilen liegen nicht auf konstantem z; "
                             "die exakte Punktlokalisierung setzt das voraus.  Ein Aufrufer, "
                             "der die Zeilen bewusst verformt hat, setzt require_level_rows = "
                             "false und lokalisiert Punkte selbst.");
  if (!cell_source.empty() && static_cast<Index>(cell_source.size()) != m.n_cells())
    throw std::runtime_error("AxisymProblem: cell_source hat nicht die Zellenzahl");
  if (!node_source_density.empty() &&
      static_cast<Index>(node_source_density.size()) != m.n_nodes())
    throw std::runtime_error("AxisymProblem: node_source_density hat nicht die Knotenzahl");
  if (!node_charge.empty() && static_cast<Index>(node_charge.size()) != m.n_nodes())
    throw std::runtime_error("AxisymProblem: node_charge hat nicht die Knotenzahl");
  // The coefficient bound depends on what the coefficient MEANS.  A relative
  // permittivity below one is unphysical; a conductivity of 1.5 S/m or a
  // viscosity of 0.036 Pa s are perfectly ordinary.  So the strict bound is
  // applied only to the electrostatic use, which is the one with the default
  // scale, and the general requirement is that the operator stay positive.
  const bool electrostatic = (coefficient_scale == eps0);
  for (Index c = 0; c < m.n_cells(); ++c) {
    if (!active[static_cast<std::size_t>(c)]) continue;
    const Real k = eps_r[static_cast<std::size_t>(c)];
    if (electrostatic && !(k >= 1.0))
      throw std::runtime_error("AxisymProblem: aktive Zelle mit eps_r < 1");
    if (!(k > 0.0))
      throw std::runtime_error("AxisymProblem: aktive Zelle mit nicht positivem Koeffizienten");
  }
  if (!(coefficient_scale > 0.0))
    throw std::runtime_error("AxisymProblem: coefficient_scale muss positiv sein");
  bool any_fixed = false;
  for (char f : fixed) any_fixed |= (f != 0);
  if (!any_fixed && far_field != FarField::Grounded)
    throw std::runtime_error("AxisymProblem: keine Dirichlet-Bedingung -- das Problem waere "
                             "nur bis auf eine Konstante bestimmt");
  if (far_field == FarField::Asymptotic && far_edges.empty())
    throw std::runtime_error("AxisymProblem: asymptotischer Fernrand ohne Randsegmente");
}

// ---------------------------------------------------------------------------

AxisymSolution solve_axisym(const AxisymProblem& p, LinearSolver which,
                            std::size_t memory_cap_bytes) {
  p.check();
  const QuadMesh& m = *p.mesh;
  const Index N = m.n_nodes();

  // Numbering: the narrow direction runs fastest, so the half bandwidth is
  // min(nr, nz) + 1 rather than whichever of the two happens to be first.
  const bool r_fastest = (m.nr <= m.nz);
  auto idx = [&](Index i, Index j) { return r_fastest ? j * m.nr + i : i * m.nz + j; };
  /// Mesh node index (always j*nr + i) -> internal index.
  auto to_internal = [&](Index n) { return idx(n % m.nr, n / m.nr); };
  const Index band = (r_fastest ? m.nr : m.nz) + 1;

  // Two band matrices are held at once: the one that gets factorised, and an
  // untouched copy of the assembled operator, which is what turns nodal
  // reactions into conductor charges without a hand-rolled surface integral.
  const std::size_t need = 2 * static_cast<std::size_t>(N) *
                           static_cast<std::size_t>(band + 1) * sizeof(Real);
  if (which == LinearSolver::Band && need > memory_cap_bytes) {
    char buf[360];
    std::snprintf(buf, sizeof buf,
                  "axisym_fem: die Bandfaktorisierung braeuchte %.2f GiB (N = %lld, "
                  "Halbbandbreite %lld, zwei Bandmatrizen) und ueberschreitet die Grenze "
                  "von %.2f GiB. Netzstufe verringern.",
                  need / 1073741824.0, static_cast<long long>(N),
                  static_cast<long long>(band), memory_cap_bytes / 1073741824.0);
    throw std::runtime_error(buf);
  }

  // --- assemble the unreduced system ---------------------------------------
  BandMatrix K(N, band);
  std::vector<Real> f(static_cast<std::size_t>(N), 0.0);

  for (Index j = 0; j + 1 < m.nz; ++j) {
    for (Index i = 0; i + 1 < m.nr; ++i) {
      const Index c = m.cell(i, j);
      if (!p.active[static_cast<std::size_t>(c)]) continue;
      Real Ke[4][4];
      element_matrix(m, i, j, p.coefficient_scale * p.eps_r[static_cast<std::size_t>(c)], Ke);
      const Index gi[4] = {idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)};
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
          if (gi[a] >= gi[b]) K.add(gi[a], gi[b], Ke[a][b]);
      if (!p.cell_source.empty()) {
        Real fe[4];
        element_source(m, i, j, p.cell_source[static_cast<std::size_t>(c)], fe);
        for (int a = 0; a < 4; ++a) f[static_cast<std::size_t>(gi[a])] += fe[a];
      }
      if (!p.node_source_density.empty()) {
        // int rho N_a 2 pi r dA with rho bilinear inside the element.
        const std::array<Index, 4> cn = m.cell_nodes(i, j);
        const Real rho[4] = {p.node_source_density[static_cast<std::size_t>(cn[0])],
                             p.node_source_density[static_cast<std::size_t>(cn[1])],
                             p.node_source_density[static_cast<std::size_t>(cn[2])],
                             p.node_source_density[static_cast<std::size_t>(cn[3])]};
        for (int q = 0; q < 4; ++q) {
          const GaussPoint g = gauss_point(m, i, j, kGaussXi[q], kGaussEta[q]);
          const Real xi = kGaussXi[q], eta = kGaussEta[q];
          const Real N[4] = {(1 - xi) * (1 - eta), xi * (1 - eta), xi * eta, (1 - xi) * eta};
          Real r_here = 0.0;
          for (int a = 0; a < 4; ++a) r_here += N[a] * rho[a];
          const Real w = r_here * 2.0 * pi * g.r * g.detJ * kGaussW;
          for (int a = 0; a < 4; ++a) f[static_cast<std::size_t>(gi[a])] += w * N[a];
        }
      }
    }
  }

  // --- discrete nodal charges ----------------------------------------------
  //
  // sum_p q_p N_a(x_p) has already been accumulated per node by the caller;
  // here it only has to reach the load vector.  No 2 pi r: see the header.
  if (!p.node_charge.empty())
    for (Index n = 0; n < N; ++n)
      f[static_cast<std::size_t>(to_internal(n))] += p.node_charge[static_cast<std::size_t>(n)];

  // --- far-field Robin term -------------------------------------------------
  //
  //   \int eps grad phi . grad v dV + \oint alpha phi v dS = 0,
  //   alpha = eps0 * (n . (x - x0)) / |x - x0|^2   [F/m^2]
  //
  // which is the monopole decay written as a natural condition.  alpha > 0 for
  // every boundary point of a box containing x0, so the system stays positive
  // definite.
  std::vector<std::array<Index, 2>> robin_edges;
  if (p.far_field == FarField::Asymptotic) {
    robin_edges = p.far_edges;
    const Real gq[2] = {0.5 - 0.5 / 1.7320508075688772935, 0.5 + 0.5 / 1.7320508075688772935};
    for (const auto& e : robin_edges) {
      const Vec2 a = m.nodes[static_cast<std::size_t>(e[0])];
      const Vec2 b = m.nodes[static_cast<std::size_t>(e[1])];
      const Real len = norm(b - a);
      if (!(len > 0.0)) continue;
      const Vec2 n = outward_normal(a, b, p.far_field_origin);
      Real Me[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
      for (int q = 0; q < 2; ++q) {
        const Real t = gq[q];
        const Vec2 x = a + t * (b - a);
        const Vec2 d = x - p.far_field_origin;
        const Real R2 = norm2(d);
        if (!(R2 > 0.0)) continue;
        const Real alpha = eps0 * dot(n, d) / R2;
        const Real Nf[2] = {1.0 - t, t};
        const Real w = alpha * 2.0 * pi * x.r * len * 0.5;  // 2-pt Gauss weight = 1/2
        for (int u = 0; u < 2; ++u)
          for (int v = 0; v < 2; ++v) Me[u][v] += w * Nf[u] * Nf[v];
      }
      // far_edges carries MESH node indices; the matrix uses the internal
      // numbering, which differs whenever the fast direction is not r.  Mixing
      // the two scatters the boundary term into unrelated matrix entries -- it
      // is invisible when nr <= nz, because the two numberings then coincide,
      // and it destroys the solution when they do not.
      const Index gi[2] = {to_internal(e[0]), to_internal(e[1])};
      for (int u = 0; u < 2; ++u)
        for (int v = 0; v < 2; ++v)
          if (gi[u] >= gi[v]) K.add(gi[u], gi[v], Me[u][v]);
    }
  }

  std::vector<char> fixed_g(static_cast<std::size_t>(N), 0);
  std::vector<Real> value_g(static_cast<std::size_t>(N), 0.0);
  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i) {
      const Index n_mesh = m.node(i, j);
      if (p.fixed[static_cast<std::size_t>(n_mesh)]) {
        fixed_g[static_cast<std::size_t>(idx(i, j))] = 1;
        value_g[static_cast<std::size_t>(idx(i, j))] =
            p.fixed_value[static_cast<std::size_t>(n_mesh)];
      }
    }

  // A free node that touches no active cell has no equation at all, and the
  // factorisation would report a meaningless singularity twenty lines later.
  // It means the caller deactivated a region -- the interior of an ideal
  // conductor -- without fixing its nodes.  Say so where it can be fixed.
  for (Index i = 0; i < N; ++i)
    if (!fixed_g[static_cast<std::size_t>(i)] && K.at(i, i) == 0.0)
      throw std::runtime_error(
          "axisym_fem: es gibt freie Knoten, die an keine aktive Zelle grenzen. Ein Gebiet "
          "wurde als inaktiv markiert (Inneres eines idealen Leiters), ohne seine Knoten "
          "festzuhalten. Das Gleichungssystem waere dort leer.");

  // Keep a copy of the unreduced operator for the reaction / residual.  It is
  // the same band matrix; copying it costs one more factor-sized allocation,
  // which is the price of getting the charge from the operator rather than from
  // a hand-rolled surface integral.
  BandMatrix K0 = K;

  // --- Dirichlet elimination, symmetric --------------------------------------

  for (Index d = 0; d < N; ++d) {
    if (!fixed_g[static_cast<std::size_t>(d)]) continue;
    const Real g = value_g[static_cast<std::size_t>(d)];
    if (g == 0.0) continue;
    const Index lo = std::max<Index>(0, d - band), hi = std::min<Index>(N - 1, d + band);
    for (Index i = lo; i <= hi; ++i) {
      if (fixed_g[static_cast<std::size_t>(i)]) continue;
      f[static_cast<std::size_t>(i)] -= K.get(i, d) * g;
    }
  }
  for (Index d = 0; d < N; ++d) {
    if (!fixed_g[static_cast<std::size_t>(d)]) continue;
    const Index lo = std::max<Index>(0, d - band), hi = std::min<Index>(N - 1, d + band);
    for (Index i = lo; i <= hi; ++i)
      if (i != d) K.set(i, d, 0.0);
    K.set(d, d, 1.0);
    f[static_cast<std::size_t>(d)] = value_g[static_cast<std::size_t>(d)];
  }

  // --- solve -----------------------------------------------------------------
  std::vector<Real> x = f;
  if (which == LinearSolver::Dense) {
    Matrix A(N, N);
    for (Index i = 0; i < N; ++i)
      for (Index j2 = 0; j2 < N; ++j2) A(i, j2) = K.get(i, j2);
    if (!solve_dense(A, x))
      throw std::runtime_error("axisym_fem: dichte Faktorisierung singulaer");
  } else {
    if (!K.factor())
      throw std::runtime_error("axisym_fem: Bandfaktorisierung nicht positiv definit -- "
                               "das deutet auf ein degeneriertes Netz oder eine fehlende "
                               "Dirichlet-Bedingung hin");
    K.solve(x);
  }

  // --- unpack, reactions, residual -------------------------------------------
  AxisymSolution s;
  s.n_nodes = N;
  s.half_bandwidth = band;
  s.factor_bytes = need;
  s.phi.assign(static_cast<std::size_t>(N), 0.0);
  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i)
      s.phi[static_cast<std::size_t>(m.node(i, j))] = x[static_cast<std::size_t>(idx(i, j))];

  s.n_free = 0;
  for (char c : p.fixed)
    if (!c) ++s.n_free;

  // R = K0 * phi in the internal numbering, mapped back to mesh order.
  std::vector<Real> R(static_cast<std::size_t>(N), 0.0);
  for (Index i = 0; i < N; ++i) {
    Real acc = 0.0;
    const Index lo = std::max<Index>(0, i - band), hi = std::min<Index>(N - 1, i + band);
    for (Index k = lo; k <= hi; ++k) acc += K0.get(i, k) * x[static_cast<std::size_t>(k)];
    R[static_cast<std::size_t>(i)] = acc;
  }
  s.reaction.assign(static_cast<std::size_t>(N), 0.0);
  s.residual_inf = 0.0;
  for (Index j = 0; j < m.nz; ++j)
    for (Index i = 0; i < m.nr; ++i) {
      const Index n_mesh = m.node(i, j);
      const Real v = R[static_cast<std::size_t>(idx(i, j))];
      s.reaction[static_cast<std::size_t>(n_mesh)] = v;
      if (!p.fixed[static_cast<std::size_t>(n_mesh)])
        s.residual_inf = std::max(s.residual_inf, std::abs(v));
    }
  return s;
}

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

bool locate(const QuadMesh& m, Vec2 x, Index* io, Index* jo, Real* xio, Real* etao) {
  if (m.nz < 2 || m.nr < 2) return false;
  const Real z0 = m.z_of_row(0), z1 = m.z_of_row(m.nz - 1);
  if (x.z < z0 || x.z > z1) return false;
  Index lo = 0, hi = m.nz - 1;
  while (hi - lo > 1) {
    const Index mid = (lo + hi) / 2;
    (m.z_of_row(mid) <= x.z ? lo : hi) = mid;
  }
  const Index j = lo;
  const Real za = m.z_of_row(j), zb = m.z_of_row(j + 1);
  Real eta = (zb > za) ? (x.z - za) / (zb - za) : 0.0;
  eta = std::min(1.0, std::max(0.0, eta));

  auto radius = [&](Index i) { return (1.0 - eta) * m.at(i, j).r + eta * m.at(i, j + 1).r; };
  if (x.r < radius(0) - 1e-15 || x.r > radius(m.nr - 1) + 1e-15) return false;
  Index rlo = 0, rhi = m.nr - 1;
  while (rhi - rlo > 1) {
    const Index mid = (rlo + rhi) / 2;
    (radius(mid) <= x.r ? rlo : rhi) = mid;
  }
  const Index i = rlo;
  const Real ra = radius(i), rb = radius(i + 1);
  Real xi = (rb > ra) ? (x.r - ra) / (rb - ra) : 0.0;
  xi = std::min(1.0, std::max(0.0, xi));

  *io = i;
  *jo = j;
  *xio = xi;
  *etao = eta;
  return true;
}

Real potential_at(const QuadMesh& m, const std::vector<Real>& phi, Vec2 x) {
  Index i, j;
  Real xi, eta;
  if (!locate(m, x, &i, &j, &xi, &eta))
    throw std::runtime_error("potential_at: Punkt liegt ausserhalb des Netzes");
  Real N[4];
  shape(xi, eta, N);
  const auto c = m.cell_nodes(i, j);
  Real v = 0.0;
  for (int k = 0; k < 4; ++k) v += N[k] * phi[static_cast<std::size_t>(c[k])];
  return v;
}

Real potential_in_cell(const QuadMesh& m, const std::vector<Real>& phi, Index i, Index j,
                       Real xi, Real eta) {
  Real N[4];
  shape(xi, eta, N);
  const auto c = m.cell_nodes(i, j);
  Real v = 0.0;
  for (int k = 0; k < 4; ++k) v += N[k] * phi[static_cast<std::size_t>(c[k])];
  return v;
}

Vec2 field_in_cell(const QuadMesh& m, const std::vector<Real>& phi, Index i, Index j, Real xi,
                   Real eta) {
  const GaussPoint g = gauss_point(m, i, j, xi, eta);
  const auto c = m.cell_nodes(i, j);
  Real dr = 0.0, dz = 0.0;
  for (int k = 0; k < 4; ++k) {
    const Real v = phi[static_cast<std::size_t>(c[k])];
    dr += g.dNdr[k] * v;
    dz += g.dNdz[k] * v;
  }
  return {-dr, -dz};
}

Vec2 field_at(const QuadMesh& m, const std::vector<Real>& phi, Vec2 x) {
  Index i, j;
  Real xi, eta;
  if (!locate(m, x, &i, &j, &xi, &eta))
    throw std::runtime_error("field_at: Punkt liegt ausserhalb des Netzes");
  return field_in_cell(m, phi, i, j, xi, eta);
}

Real QuadMesh::min_jacobian() const {
  Real lo = 0.0;
  bool first = true;
  for (Index j = 0; j + 1 < nz; ++j)
    for (Index i = 0; i + 1 < nr; ++i)
      for (int q = 0; q < 4; ++q) {
        const GaussPoint g = gauss_point(*this, i, j, kGaussXi[q], kGaussEta[q]);
        if (first || g.detJ < lo) {
          lo = g.detJ;
          first = false;
        }
      }
  return lo;
}

Vec2 field_recovered_at_node(const QuadMesh& m, const std::vector<Real>& phi,
                             const std::vector<Real>& eps_r, const std::vector<char>& active,
                             Index i, Index j, Real eps_select) {
  Vec2 acc{0.0, 0.0};
  Real w_sum = 0.0;
  for (Index dj = -1; dj <= 0; ++dj)
    for (Index di = -1; di <= 0; ++di) {
      const Index ci = i + di, cj = j + dj;
      if (ci < 0 || cj < 0 || ci + 1 >= m.nr || cj + 1 >= m.nz) continue;
      const Index c = m.cell(ci, cj);
      if (!active[static_cast<std::size_t>(c)]) continue;
      if (eps_r[static_cast<std::size_t>(c)] != eps_select) continue;
      const Real w = m.cell_revolved_volume(ci, cj);
      acc += w * field_in_cell(m, phi, ci, cj, 0.5, 0.5);
      w_sum += w;
    }
  if (!(w_sum > 0.0))
    throw std::runtime_error("field_recovered_at_node: kein passendes aktives Nachbarelement");
  return acc / w_sum;
}

Vec2 field_recovered_at(const QuadMesh& m, const std::vector<Real>& phi,
                        const std::vector<Real>& eps_r, const std::vector<char>& active,
                        Vec2 x) {
  Index i, j;
  Real xi, eta;
  if (!locate(m, x, &i, &j, &xi, &eta))
    throw std::runtime_error("field_recovered_at: Punkt liegt ausserhalb des Netzes");
  const Index c = m.cell(i, j);
  if (!active[static_cast<std::size_t>(c)])
    throw std::runtime_error("field_recovered_at: Punkt liegt in einer inaktiven Zelle "
                             "(Inneres eines idealen Leiters)");
  const Real e = eps_r[static_cast<std::size_t>(c)];
  Real N[4];
  shape(xi, eta, N);
  const Index ii[4] = {i, i + 1, i + 1, i};
  const Index jj[4] = {j, j, j + 1, j + 1};
  Vec2 E{0.0, 0.0};
  for (int k = 0; k < 4; ++k)
    E += N[k] * field_recovered_at_node(m, phi, eps_r, active, ii[k], jj[k], e);
  return E;
}

Real charge_of(const AxisymSolution& s, const std::vector<char>& mask) {
  if (mask.size() != s.reaction.size())
    throw std::runtime_error("charge_of: Maske hat nicht die Knotenzahl");
  Real q = 0.0;
  for (std::size_t k = 0; k < mask.size(); ++k)
    if (mask[k]) q += s.reaction[k];
  return q;
}

}  // namespace es
