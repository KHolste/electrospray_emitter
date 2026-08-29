#include "es/linalg.hpp"

#include <cmath>

namespace es {

bool lu_factor(Matrix& A, std::vector<Index>& piv) {
  const Index n = A.rows();
  piv.resize(static_cast<std::size_t>(n));
  for (Index i = 0; i < n; ++i) piv[static_cast<std::size_t>(i)] = i;

  for (Index k = 0; k < n; ++k) {
    Index p = k;
    Real best = std::abs(A(k, k));
    for (Index i = k + 1; i < n; ++i) {
      const Real v = std::abs(A(i, k));
      if (v > best) { best = v; p = i; }
    }
    if (best == 0.0) return false;
    if (p != k) {
      for (Index j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
      std::swap(piv[static_cast<std::size_t>(k)], piv[static_cast<std::size_t>(p)]);
    }
    const Real inv = 1.0 / A(k, k);
    for (Index i = k + 1; i < n; ++i) {
      const Real f = A(i, k) * inv;
      A(i, k) = f;
      if (f == 0.0) continue;
      Real* ai = A.row(i);
      const Real* ak = A.row(k);
      for (Index j = k + 1; j < n; ++j) ai[j] -= f * ak[j];
    }
  }
  return true;
}

void lu_solve(const Matrix& LU, const std::vector<Index>& piv, std::vector<Real>& b) {
  const Index n = LU.rows();
  std::vector<Real> x(static_cast<std::size_t>(n));
  for (Index i = 0; i < n; ++i) x[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(piv[static_cast<std::size_t>(i)])];
  // forward substitution (unit lower triangle)
  for (Index i = 1; i < n; ++i) {
    Real s = x[static_cast<std::size_t>(i)];
    const Real* li = LU.row(i);
    for (Index j = 0; j < i; ++j) s -= li[j] * x[static_cast<std::size_t>(j)];
    x[static_cast<std::size_t>(i)] = s;
  }
  // back substitution
  for (Index i = n - 1; i >= 0; --i) {
    Real s = x[static_cast<std::size_t>(i)];
    const Real* li = LU.row(i);
    for (Index j = i + 1; j < n; ++j) s -= li[j] * x[static_cast<std::size_t>(j)];
    x[static_cast<std::size_t>(i)] = s / li[i];
  }
  b.swap(x);
}

std::vector<Real> matvec(const Matrix& A, const std::vector<Real>& x) {
  std::vector<Real> y(static_cast<std::size_t>(A.rows()), 0.0);
  for (Index i = 0; i < A.rows(); ++i) {
    const Real* ai = A.row(i);
    Real s = 0.0;
    for (Index j = 0; j < A.cols(); ++j) s += ai[j] * x[static_cast<std::size_t>(j)];
    y[static_cast<std::size_t>(i)] = s;
  }
  return y;
}

bool solve_dense(const Matrix& A, std::vector<Real>& b) {
  Matrix LU = A;
  std::vector<Index> piv;
  if (!lu_factor(LU, piv)) return false;

  const std::vector<Real> rhs = b;
  lu_solve(LU, piv, b);

  // one step of iterative refinement
  std::vector<Real> res = matvec(A, b);
  for (std::size_t i = 0; i < res.size(); ++i) res[i] = rhs[i] - res[i];
  lu_solve(LU, piv, res);
  for (std::size_t i = 0; i < b.size(); ++i) b[i] += res[i];
  return true;
}

}  // namespace es
