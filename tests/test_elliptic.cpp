#include <cmath>
#include <cstdio>

#include "es/constants.hpp"
#include "es/elliptic.hpp"

using namespace es;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("%-42s got=%-22.15g want=%-22.15g relerr=%-10.3g %s\n", what, got, want, err,
              ok ? "OK" : "FAIL");
}

/// Brute-force reference: adaptive Simpson on the defining integrals.
static Real ref_K(Real m) {
  const int n = 2000001;
  const Real hstep = 0.5 * constants::pi / (n - 1);
  Real s = 0.0;
  for (int i = 0; i < n; ++i) {
    const Real t = i * hstep;
    const Real f = 1.0 / std::sqrt(1.0 - m * std::sin(t) * std::sin(t));
    const Real w = (i == 0 || i == n - 1) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
    s += w * f;
  }
  return s * hstep / 3.0;
}

static Real ref_E(Real m) {
  const int n = 2000001;
  const Real hstep = 0.5 * constants::pi / (n - 1);
  Real s = 0.0;
  for (int i = 0; i < n; ++i) {
    const Real t = i * hstep;
    const Real f = std::sqrt(1.0 - m * std::sin(t) * std::sin(t));
    const Real w = (i == 0 || i == n - 1) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
    s += w * f;
  }
  return s * hstep / 3.0;
}

int main() {
  std::printf("=== complete elliptic integrals (parameter convention m = k^2) ===\n");

  // Exact special values
  check("K(0)", ellipK(0.0), 0.5 * constants::pi, 1e-15);
  check("E(0)", ellipE(0.0), 0.5 * constants::pi, 1e-15);
  check("E(1)", ellipE(1.0), 1.0, 1e-12);
  // K(1/2) = Gamma(1/4)^2 / (4 sqrt(pi))
  check("K(1/2)", ellipK(0.5), 1.8540746773013719, 1e-14);
  check("E(1/2)", ellipE(0.5), 1.3506438810476755, 1e-14);

  // Against numerical quadrature
  for (Real m : {0.1, 0.35, 0.7, 0.9, 0.99}) {
    char buf[80];
    std::snprintf(buf, sizeof buf, "K(%.2f) vs Simpson", m);
    check(buf, ellipK(m), ref_K(m), 5e-10);
    std::snprintf(buf, sizeof buf, "E(%.2f) vs Simpson", m);
    check(buf, ellipE(m), ref_E(m), 5e-10);
  }

  // Near-singular limit: K(m) -> ln(4/sqrt(1-m)) as m -> 1.  This is the regime
  // BEM near-panels live in, so it must stay accurate.
  for (Real mc : {1e-6, 1e-10, 1e-14}) {
    Real K, E;
    ellipKE_mc(mc, K, E);
    const Real asym = std::log(4.0 / std::sqrt(mc));
    char buf[80];
    std::snprintf(buf, sizeof buf, "K(mc=%.0e) vs ln(4/k')", mc);
    // The asymptote itself is only good to O(mc log mc), so loosen accordingly.
    check(buf, K, asym, 10.0 * mc + 1e-13);
  }

  // Derivatives against central differences
  std::printf("\n=== derivatives ===\n");
  for (Real m : {1e-8, 0.01, 0.3, 0.8, 0.95}) {
    Real dK, dE;
    ellipKE_deriv(m, dK, dE);
    const Real hh = 1e-6 * std::max(m, 1e-3);
    const Real dKfd = (ellipK(m + hh) - ellipK(m - hh)) / (2 * hh);
    const Real dEfd = (ellipE(m + hh) - ellipE(m - hh)) / (2 * hh);
    char buf[80];
    std::snprintf(buf, sizeof buf, "dK/dm(%.0e) vs FD", m);
    check(buf, dK, dKfd, 1e-6);
    std::snprintf(buf, sizeof buf, "dE/dm(%.0e) vs FD", m);
    check(buf, dE, dEfd, 1e-6);
  }

  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
