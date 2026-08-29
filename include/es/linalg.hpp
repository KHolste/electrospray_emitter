#pragma once
#include <vector>

#include "es/types.hpp"

namespace es {

/// Dense row-major matrix.  Sized for axisymmetric BEM (N ~ 10^3), so a plain
/// LU factorisation with partial pivoting is entirely adequate.
class Matrix {
 public:
  Matrix() = default;
  Matrix(Index n, Index m) : n_(n), m_(m), a_(static_cast<std::size_t>(n * m), 0.0) {}

  Real& operator()(Index i, Index j) { return a_[static_cast<std::size_t>(i * m_ + j)]; }
  Real operator()(Index i, Index j) const { return a_[static_cast<std::size_t>(i * m_ + j)]; }

  Index rows() const { return n_; }
  Index cols() const { return m_; }
  Real* data() { return a_.data(); }
  const Real* data() const { return a_.data(); }
  Real* row(Index i) { return a_.data() + i * m_; }
  const Real* row(Index i) const { return a_.data() + i * m_; }

 private:
  Index n_{0}, m_{0};
  std::vector<Real> a_;
};

/// In-place LU with partial pivoting.  Returns false on (numerically) singular.
bool lu_factor(Matrix& A, std::vector<Index>& piv);

/// Solve A x = b using a factorisation produced by lu_factor (b is overwritten).
void lu_solve(const Matrix& LU, const std::vector<Index>& piv, std::vector<Real>& b);

/// Convenience: factor a copy and solve.  One step of iterative refinement is
/// applied, which costs one extra mat-vec and typically recovers the digits lost
/// to the ill-conditioning of graded BEM meshes.
bool solve_dense(const Matrix& A, std::vector<Real>& b);

std::vector<Real> matvec(const Matrix& A, const std::vector<Real>& x);

}  // namespace es
