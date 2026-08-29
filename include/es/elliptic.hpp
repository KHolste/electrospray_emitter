#pragma once
#include "es/types.hpp"

namespace es {

/// Complete elliptic integrals in the PARAMETER convention m = k^2:
///   K(m) = \int_0^{pi/2} dt / sqrt(1 - m sin^2 t)
///   E(m) = \int_0^{pi/2} sqrt(1 - m sin^2 t) dt
/// Both are evaluated with the arithmetic-geometric mean, which stays accurate
/// right up to m -> 1 (the regime that dominates near-singular BEM kernels).
Real ellipK(Real m);
Real ellipE(Real m);

/// K and E in one AGM pass (the BEM field kernel needs both).
void ellipKE(Real m, Real& K, Real& E);

/// Same, but parameterised by the COMPLEMENTARY parameter mc = 1 - m.  BEM
/// kernels naturally produce mc = d^2/S^2 (d = meridian distance between the
/// two rings), so passing mc through directly avoids the catastrophic
/// cancellation of forming 1 - m when m -> 1.  This is what keeps near-singular
/// panels accurate.
void ellipKE_mc(Real mc, Real& K, Real& E);

/// dK/dm and dE/dm, with series expansions used near m = 0 where the closed
/// forms are 0/0.
void ellipKE_deriv(Real m, Real& dKdm, Real& dEdm);

}  // namespace es
