#include "es/elliptic.hpp"

#include <limits>

#include "es/constants.hpp"

namespace es {
namespace {

constexpr int kMaxAgm = 40;
constexpr Real kAgmTol = 1e-17;

}  // namespace

void ellipKE_mc(Real mc, Real& K, Real& E) {
  if (mc > 1.0) mc = 1.0;
  if (mc <= 0.0) {
    K = std::numeric_limits<Real>::infinity();
    E = 1.0;
    return;
  }
  const Real m = 1.0 - mc;

  // Standard AGM recursion.  a_{n+1} = (a+b)/2, b_{n+1} = sqrt(a b),
  // c_{n+1} = (a-b)/2, with a0 = 1, b0 = sqrt(1-m), c0 = sqrt(m).
  //   K = pi / (2 a_N)
  //   E = K * (1 - sum_{n>=0} 2^{n-1} c_n^2)
  Real a = 1.0;
  Real b = std::sqrt(mc);
  Real c = std::sqrt(m);
  Real sum = 0.5 * m;  // n = 0 term: 2^{-1} c_0^2
  Real pow2 = 1.0;     // 2^{n-1} for n = 1 is 2^0

  for (int n = 1; n < kMaxAgm; ++n) {
    const Real an = 0.5 * (a + b);
    const Real bn = std::sqrt(a * b);
    c = 0.5 * (a - b);
    a = an;
    b = bn;
    sum += pow2 * c * c;
    pow2 *= 2.0;
    if (c < kAgmTol * a) break;
  }

  K = constants::pi / (2.0 * a);
  E = K * (1.0 - sum);
}

void ellipKE(Real m, Real& K, Real& E) {
  if (m < 0.0) m = 0.0;
  ellipKE_mc(1.0 - m, K, E);
}

Real ellipK(Real m) {
  Real K, E;
  ellipKE(m, K, E);
  return K;
}

Real ellipE(Real m) {
  Real K, E;
  ellipKE(m, K, E);
  return E;
}

void ellipKE_deriv(Real m, Real& dKdm, Real& dEdm) {
  // Exact:  dK/dm = [E - (1-m) K] / (2 m (1-m)),   dE/dm = (E - K) / (2 m).
  // Both are 0/0 at m = 0, so switch to the hypergeometric series there.
  constexpr Real kSmall = 1e-6;
  if (m < kSmall) {
    // K = (pi/2)(1 + m/4 + 9m^2/64 + 25m^3/256 + ...)
    // E = (pi/2)(1 - m/4 - 3m^2/64 -  5m^3/256 - ...)
    const Real hp = 0.5 * constants::pi;
    dKdm = hp * (0.25 + (9.0 / 32.0) * m + (75.0 / 256.0) * m * m);
    dEdm = hp * (-0.25 - (3.0 / 32.0) * m - (15.0 / 256.0) * m * m);
    return;
  }
  Real K, E;
  ellipKE(m, K, E);
  const Real om = 1.0 - m;
  dKdm = (E - om * K) / (2.0 * m * om);
  dEdm = (E - K) / (2.0 * m);
}

}  // namespace es
