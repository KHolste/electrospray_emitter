#include "es/capillary.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/constants.hpp"

namespace es {

using constants::pi;

// ---------------------------------------------------------------------------
// Status
// ---------------------------------------------------------------------------

const char* to_string(CapillaryStatus s) {
  switch (s) {
    case CapillaryStatus::Solved: return "Solved";
    case CapillaryStatus::PressureOutsideCapillaryRange: return "PressureOutsideCapillaryRange";
    case CapillaryStatus::HemisphericalLimit: return "HemisphericalLimit";
    case CapillaryStatus::InvalidGeometry: return "InvalidGeometry";
    case CapillaryStatus::InvalidLiquid: return "InvalidLiquid";
    case CapillaryStatus::ContactAngleAndPinningBothPrescribed:
      return "ContactAngleAndPinningBothPrescribed";
    case CapillaryStatus::AccuracyNotReached: return "AccuracyNotReached";
    case CapillaryStatus::ArclengthNotBracketed: return "ArclengthNotBracketed";
    case CapillaryStatus::NotAttempted: return "NotAttempted";
  }
  return "?";
}

const char* explain(CapillaryStatus s) {
  switch (s) {
    case CapillaryStatus::Solved:
      return "Statisches Kapillargleichgewicht gefunden und bis zur geforderten Genauigkeit "
             "aufgeloest.";
    case CapillaryStatus::PressureOutsideCapillaryRange:
      return "Zu diesem Druck gibt es bei diesem Pinningradius keine glatte Oberflaeche: die "
             "Meridiankurve wird senkrecht, bevor sie den Kontaktradius erreicht (|Pi| > 2). "
             "Es wird keine Ersatzform und kein letzter Iterationsstand zurueckgegeben.";
    case CapillaryStatus::HemisphericalLimit:
      return "Der Druck liegt auf der hemisphaerischen Grenze |Pi| = 2. Die Halbkugel ist eine "
             "gueltige Loesung und geschlossen bekannt, aber die Schiessbedingung r(L) = a wird "
             "dort tangential (doppelte Nullstelle). Die numerische Form wird deshalb nicht als "
             "geloest ausgegeben.";
    case CapillaryStatus::InvalidGeometry:
      return "Der Pinningradius aus der Geraetegeometrie ist nicht verwendbar.";
    case CapillaryStatus::InvalidLiquid:
      return "Der Stoffdatensatz ist nicht verwendbar.";
    case CapillaryStatus::ContactAngleAndPinningBothPrescribed:
      return "Kontaktwinkel und gepinnte Kontaktlinie wurden gleichzeitig vorgeschrieben. Das "
             "ueberbestimmt die Kante; an einer scharfen, gepinnten Austrittskante gilt kein "
             "zusaetzlicher Young-Winkel.";
    case CapillaryStatus::AccuracyNotReached:
      return "Die automatische Verfeinerung hat ihre Obergrenze erreicht, ohne die geforderte "
             "Genauigkeit zu unterschreiten.";
    case CapillaryStatus::ArclengthNotBracketed:
      return "Die Bogenlaengensuche hat innerhalb ihrer Schranken keine Einschliessung "
             "gefunden -- ein numerisches, kein physikalisches Scheitern.";
    case CapillaryStatus::NotAttempted:
      return "Nicht gerechnet.";
  }
  return "?";
}

namespace capillary {

Real pi_from_pressure(Real delta_p, Real a, Real gamma) {
  if (!(gamma > 0.0)) throw std::runtime_error("gamma muss positiv sein.");
  return delta_p * a / gamma;
}

Real pressure_from_pi(Real Pi, Real a, Real gamma) {
  if (!(a > 0.0)) throw std::runtime_error("Der Pinningradius muss positiv sein.");
  return Pi * gamma / a;
}

}  // namespace capillary

// ---------------------------------------------------------------------------
// Closed-form spherical cap
// ---------------------------------------------------------------------------

bool spherical_cap_exists(Real a, Real delta_p, Real gamma) {
  if (!(a > 0.0) || !(gamma > 0.0) || !std::isfinite(delta_p)) return false;
  return std::abs(capillary::pi_from_pressure(delta_p, a, gamma)) <= capillary::kMaxAbsPi;
}

SphericalCap spherical_cap(Real a, Real delta_p, Real gamma) {
  if (!(a > 0.0)) throw std::runtime_error("Der Pinningradius muss positiv sein.");
  if (!(gamma > 0.0)) throw std::runtime_error("Die Oberflaechenspannung muss positiv sein.");
  if (!std::isfinite(delta_p)) throw std::runtime_error("delta_p muss endlich sein.");

  SphericalCap c;
  c.a = a;
  c.gamma = gamma;
  c.delta_p = delta_p;
  c.Pi = capillary::pi_from_pressure(delta_p, a, gamma);
  c.curvature = delta_p / gamma;

  if (std::abs(c.Pi) > capillary::kMaxAbsPi)
    throw std::runtime_error(
        "Zu Pi = " + std::to_string(c.Pi) +
        " existiert keine Kugelkappe: der Kugelradius 2*gamma/delta_p waere kleiner als der "
        "Pinningradius. Der darstellbare Bereich ist |Pi| <= 2.");

  if (delta_p == 0.0) {
    c.sphere_radius = std::numeric_limits<Real>::infinity();
    c.apex_height = 0.0;
    c.arclength = a;
    c.revolved_area = pi * a * a;
    c.revolved_volume = 0.0;
    c.contact_tangent_angle = 0.0;
    return c;
  }

  const Real R = 2.0 * gamma / delta_p;   // signed
  const Real absR = std::abs(R);
  const Real sg = (R > 0.0) ? 1.0 : -1.0;
  // Clamped only against round-off at exactly |Pi| = 2, where a == absR.
  const Real ratio = std::min(1.0, a / absR);
  const Real q = absR * std::sqrt(std::max(0.0, 1.0 - ratio * ratio));  // sqrt(R^2 - a^2)
  const Real h_abs = absR - q;

  c.sphere_radius = R;
  c.apex_height = sg * h_abs;
  c.arclength = absR * std::asin(ratio);
  c.revolved_area = 2.0 * pi * absR * h_abs;
  c.revolved_volume = sg * pi * h_abs * h_abs * (3.0 * absR - h_abs) / 3.0;
  c.contact_tangent_angle = sg * std::asin(ratio);
  return c;
}

Real SphericalCap::z_at_radius(Real r) const {
  if (delta_p == 0.0) return 0.0;
  const Real R = sphere_radius;
  const Real absR = std::abs(R);
  const Real sg = (R > 0.0) ? 1.0 : -1.0;
  const Real z_centre = -sg * std::sqrt(std::max(0.0, absR * absR - a * a));
  return z_centre + sg * std::sqrt(std::max(0.0, absR * absR - r * r));
}

// ---------------------------------------------------------------------------
// The solver
// ---------------------------------------------------------------------------

namespace {

/// State integrated along the meridian.  area and volume ride along so that
/// both are fourth-order accurate instead of being read off the polyline.
struct State {
  Real r{0}, z{0}, psi{0}, area{0}, vol{0};
};

State operator+(const State& a, const State& b) {
  return {a.r + b.r, a.z + b.z, a.psi + b.psi, a.area + b.area, a.vol + b.vol};
}
State operator*(Real s, const State& a) {
  return {s * a.r, s * a.z, s * a.psi, s * a.area, s * a.vol};
}

/// Right-hand side per unit arclength.
///
/// The axis is the only special case and it is ANALYTIC: regularity makes the
/// two principal curvatures equal at r = 0, so sin(psi)/r -> dpsi/ds there and
/// dpsi/ds = kappa/2.  Nothing is divided by zero and nothing is guessed.
State rhs(const State& y, Real kappa) {
  State d;
  const Real c = std::cos(y.psi), s = std::sin(y.psi);
  d.r = c;
  d.z = -s;
  d.psi = (y.r > 0.0) ? kappa - s / y.r : 0.5 * kappa;
  d.area = 2.0 * pi * y.r;
  d.vol = pi * y.r * y.r * s;
  return d;
}

/// Classical RK4 in the normalised arclength tau = s/L, n uniform steps.
/// If `trace` is non-null it receives the n+1 states.
State integrate(Real L, Real kappa, int n, std::vector<State>* trace) {
  const Real h = 1.0 / static_cast<Real>(n);
  State y;  // apex: r = z = psi = 0, area = 0, volume = 0
  if (trace) {
    trace->clear();
    trace->reserve(static_cast<std::size_t>(n) + 1);
    trace->push_back(y);
  }
  for (int i = 0; i < n; ++i) {
    const State k1 = L * rhs(y, kappa);
    const State k2 = L * rhs(y + (0.5 * h) * k1, kappa);
    const State k3 = L * rhs(y + (0.5 * h) * k2, kappa);
    const State k4 = L * rhs(y + h * k3, kappa);
    y = y + (h / 6.0) * (k1 + (2.0 * k2 + (2.0 * k3 + k4)));
    if (trace) trace->push_back(y);
  }
  return y;
}

struct ShootResult {
  CapillaryStatus status{CapillaryStatus::NotAttempted};
  Real L{0};
  State end;
};

/// Find the arclength L at which the meridian first reaches r = a.
///
/// The bracket is found by marching, WITHOUT using the closed-form sphere
/// radius: if the tangent turns vertical while r is still short of a, then no
/// smooth surface pinned at a exists for this pressure, and that is reported as
/// such.  Inside the bracket a safeguarded Newton iteration (dr_end/dL =
/// cos psi_end) is used, falling back to bisection whenever the Newton step
/// leaves the bracket.
ShootResult shoot(Real a, Real kappa, int n) {
  ShootResult out;
  auto eval = [&](Real L) { return integrate(L, kappa, n, nullptr); };

  Real L_lo = a;                      // r_end(L) <= L always, so F(a) <= 0
  Real L_hi = L_lo;
  Real L_below_vertical = L_lo;       // largest L known to have |psi| < pi/2
  bool bracketed = false;
  for (int k = 0; k < 60 && !bracketed; ++k) {
    L_hi *= 1.25;
    if (L_hi > 8.0 * a) break;
    const State y = eval(L_hi);
    if (y.r >= a) {
      bracketed = true;
      break;
    }
    if (std::abs(y.psi) >= 0.5 * pi) {
      // The march has stepped past the point where the meridian stands
      // vertical.  That point carries the LARGEST radius the surface attains,
      // so it, and not the step size of the march, decides whether a shape
      // pinned at a exists.  Locate it by bisection on psi and compare.
      Real lo = L_below_vertical, hi = L_hi;
      for (int it = 0; it < 200; ++it) {
        if (hi - lo <= 4.0 * std::numeric_limits<Real>::epsilon() * hi) break;
        const Real mid = 0.5 * (lo + hi);
        if (std::abs(eval(mid).psi) < 0.5 * pi) lo = mid; else hi = mid;
      }
      if (eval(hi).r >= a) {
        L_hi = hi;                    // the crossing lies in [L_lo, L_hi]
        bracketed = true;
      } else {
        // The meridian turns vertical while still short of the contact radius:
        // no smooth surface pinned at a exists for this pressure.  Decided from
        // the integrated shape, not from the closed form.
        out.status = CapillaryStatus::PressureOutsideCapillaryRange;
        return out;
      }
      break;
    }
    L_below_vertical = L_hi;
    L_lo = L_hi;                      // r_end < a here, so F(L_lo) < 0 still
  }
  if (!bracketed) {
    out.status = CapillaryStatus::ArclengthNotBracketed;
    return out;
  }

  Real L = 0.5 * (L_lo + L_hi);
  State y = eval(L);
  const Real tol_r = 8.0 * std::numeric_limits<Real>::epsilon() * a;
  for (int it = 0; it < 200; ++it) {
    const Real F = y.r - a;
    if (F <= 0.0) L_lo = L; else L_hi = L;
    if (std::abs(F) <= tol_r) break;
    if (L_hi - L_lo <= 4.0 * std::numeric_limits<Real>::epsilon() * L_hi) break;

    const Real slope = std::cos(y.psi);
    Real L_next = (slope > 1.0e-12) ? L - F / slope : 0.5 * (L_lo + L_hi);
    if (!(L_next > L_lo && L_next < L_hi)) L_next = 0.5 * (L_lo + L_hi);
    L = L_next;
    y = eval(L);
  }
  out.status = CapillaryStatus::Solved;
  out.L = L;
  out.end = y;
  return out;
}

/// One solve at a fixed resolution.  Fills everything except the refinement
/// bookkeeping.
CapillaryMeniscus solve_at(Real a, Real z_contact, Real gamma, Real delta_p, int n) {
  CapillaryMeniscus m;
  m.contact_radius = a;
  m.contact_z = z_contact;
  m.delta_p_exit = delta_p;
  m.gamma = gamma;
  m.Pi = capillary::pi_from_pressure(delta_p, a, gamma);
  m.n_intervals = n;

  const Real kappa = delta_p / gamma;
  const ShootResult sh = shoot(a, kappa, n);
  if (sh.status != CapillaryStatus::Solved) {
    m.status = sh.status;
    m.message = explain(sh.status);
    return m;
  }

  std::vector<State> trace;
  integrate(sh.L, kappa, n, &trace);
  const Real z_end = trace.back().z;

  m.nodes.reserve(trace.size());
  m.psi.reserve(trace.size());
  for (const State& y : trace) {
    m.nodes.push_back({y.r, y.z - z_end + z_contact});
    m.psi.push_back(y.psi);
  }
  m.arclength = sh.L;
  m.apex_height = m.nodes.front().z - z_contact;
  m.revolved_area = trace.back().area;
  m.revolved_volume = trace.back().vol;
  m.contact_tangent_angle = trace.back().psi;

  // Independent second evaluation from the polyline, with the tested helpers.
  m.polyline_area = revolved_area(m.nodes);
  {
    // Closed by the flat disc at z = z_contact.  The loop apex -> contact ->
    // axis runs clockwise for a bulge, so Green's formula returns the negative
    // of the enclosed volume; the sign therefore comes from the ORIENTATION of
    // the computed shape, not from the sign of the prescribed pressure.
    std::vector<Vec2> loop = m.nodes;                 // apex -> contact
    loop.push_back({0.0, z_contact});                 // close along z = z_contact
    m.polyline_volume = -revolved_volume(loop);
  }

  m.status = CapillaryStatus::Solved;
  return m;
}

/// Largest normalised change between a run at n and one at 2n intervals.
Real refinement_change(const CapillaryMeniscus& coarse, const CapillaryMeniscus& fine, Real a) {
  Real e = 0.0;
  for (std::size_t i = 0; i < coarse.nodes.size(); ++i) {
    const Vec2 d = fine.nodes[2 * i] - coarse.nodes[i];
    e = std::max(e, norm(d) / a);
  }
  e = std::max(e, std::abs(fine.apex_height - coarse.apex_height) / a);
  e = std::max(e, std::abs(fine.arclength - coarse.arclength) / a);
  e = std::max(e, std::abs(fine.revolved_area - coarse.revolved_area) / (a * a));
  e = std::max(e, std::abs(fine.revolved_volume - coarse.revolved_volume) / (a * a * a));
  return e;
}

constexpr int kStartIntervals = 64;

}  // namespace

// ---------------------------------------------------------------------------

CapillaryMeniscus solve_capillary_meniscus(Real contact_radius, Real contact_z,
                                           const LiquidProperties& liquid,
                                           const CapillaryRequest& q) {
  CapillaryMeniscus m;
  m.contact_radius = contact_radius;
  m.contact_z = contact_z;
  m.delta_p_exit = q.delta_p_exit;
  m.gamma = liquid.gamma;

  if (q.contact_angle_prescribed) {
    m.status = CapillaryStatus::ContactAngleAndPinningBothPrescribed;
    m.message = explain(m.status);
    return m;
  }
  if (!(contact_radius > 0.0) || !std::isfinite(contact_radius) || !std::isfinite(contact_z)) {
    m.status = CapillaryStatus::InvalidGeometry;
    m.message = "Der Pinningradius muss endlich und positiv sein, ist aber " +
                std::to_string(contact_radius) + " m.";
    return m;
  }
  {
    const std::string why = liquid.why_unusable();
    if (!why.empty()) {
      m.status = CapillaryStatus::InvalidLiquid;
      m.message = why;
      return m;
    }
  }
  if (!std::isfinite(q.delta_p_exit)) {
    m.status = CapillaryStatus::InvalidGeometry;
    m.message = "delta_p_exit muss endlich sein.";
    return m;
  }

  const Real Pi = capillary::pi_from_pressure(q.delta_p_exit, contact_radius, liquid.gamma);
  m.Pi = Pi;
  if (std::abs(Pi) > capillary::kMaxAbsPi) {
    m.status = CapillaryStatus::PressureOutsideCapillaryRange;
    m.message = std::string(explain(m.status)) + "  Pi = " + std::to_string(Pi) + ".";
    return m;
  }
  if (std::abs(std::abs(Pi) - capillary::kMaxAbsPi) <= capillary::kHemisphereBand) {
    m.status = CapillaryStatus::HemisphericalLimit;
    m.message = std::string(explain(m.status)) + "  Pi = " + std::to_string(Pi) + ".";
    return m;
  }

  if (q.forced_intervals > 0) {
    CapillaryMeniscus r = solve_at(contact_radius, contact_z, liquid.gamma, q.delta_p_exit,
                                  q.forced_intervals);
    r.discretisation_was_forced = true;
    r.estimated_relative_error = std::numeric_limits<Real>::quiet_NaN();
    return r;
  }

  CapillaryMeniscus coarse =
      solve_at(contact_radius, contact_z, liquid.gamma, q.delta_p_exit, kStartIntervals);
  if (coarse.status != CapillaryStatus::Solved) return coarse;

  for (int n = 2 * kStartIntervals;; n *= 2) {
    CapillaryMeniscus fine =
        solve_at(contact_radius, contact_z, liquid.gamma, q.delta_p_exit, n);
    if (fine.status != CapillaryStatus::Solved) return fine;
    const Real e = refinement_change(coarse, fine, contact_radius);
    fine.estimated_relative_error = e;
    if (e <= q.target_relative_accuracy) return fine;
    if (2 * n > q.max_intervals) {
      fine.status = CapillaryStatus::AccuracyNotReached;
      fine.message = std::string(explain(CapillaryStatus::AccuracyNotReached)) +
                     "  Erreicht wurde " + std::to_string(e) + " bei " + std::to_string(n) +
                     " Intervallen; gefordert war " +
                     std::to_string(q.target_relative_accuracy) + ".";
      return fine;
    }
    coarse = fine;
  }
}

CapillaryMeniscus solve_capillary_meniscus(const DeviceGeometry& geometry,
                                           const LiquidProperties& liquid,
                                           const CapillaryRequest& q) {
  // The exit edge comes from the device geometry and from nowhere else.
  const Vec2 edge = geometry.feature(FeatureId::PinnedContactEdge);
  if (std::abs(edge.r - geometry.contact_radius()) > 0.0) {
    CapillaryMeniscus m;
    m.status = CapillaryStatus::InvalidGeometry;
    m.message = "Die gepinnte Austrittskante und phi_2/2 stimmen nicht ueberein. Die "
                "Kapillarrechnung nimmt keine zweite Geometriebeschreibung an.";
    return m;
  }
  return solve_capillary_meniscus(edge.r, edge.z, liquid, q);
}

// ---------------------------------------------------------------------------
// Independent checks on the produced polyline
// ---------------------------------------------------------------------------

ResidualProfile young_laplace_residual(const CapillaryMeniscus& m) {
  ResidualProfile out;
  const std::size_t N = m.nodes.size();
  if (N < 4 || !(m.gamma > 0.0) || !(m.contact_radius > 0.0)) return out;

  const std::size_t nseg = N - 1;
  std::vector<Real> seg_psi(nseg), seg_len(nseg), s_node(N, 0.0);
  for (std::size_t i = 0; i < nseg; ++i) {
    const Vec2 d = m.nodes[i + 1] - m.nodes[i];
    seg_psi[i] = std::atan2(-d.z, d.r);
    seg_len[i] = norm(d);
    s_node[i + 1] = s_node[i] + seg_len[i];
  }

  const Real scale = m.contact_radius / m.gamma;  // (gamma*kappa - dp) * a / gamma
  out.s.resize(N);
  out.residual.resize(N);
  Real sum2 = 0.0;
  for (std::size_t i = 0; i < N; ++i) {
    Real kappa_total = 0.0;
    if (i == 0) {
      // Apex: both principal curvatures are equal by regularity.  The chord
      // angle of the first segment is psi at half a step, so dpsi/ds is twice
      // that over the step, and the total curvature is twice again.
      kappa_total = 4.0 * seg_psi[0] / seg_len[0];
    } else if (i + 1 < N) {
      const Real psi_node = 0.5 * (seg_psi[i - 1] + seg_psi[i]);
      const Real kappa_m = (seg_psi[i] - seg_psi[i - 1]) / (0.5 * (seg_len[i - 1] + seg_len[i]));
      kappa_total = kappa_m + std::sin(psi_node) / m.nodes[i].r;
    } else {
      const Real kappa_m =
          (seg_psi[nseg - 1] - seg_psi[nseg - 2]) / (0.5 * (seg_len[nseg - 2] + seg_len[nseg - 1]));
      const Real psi_node = seg_psi[nseg - 1] + 0.5 * (seg_psi[nseg - 1] - seg_psi[nseg - 2]);
      kappa_total = kappa_m + std::sin(psi_node) / m.nodes[i].r;
    }
    const Real res = (m.gamma * kappa_total - m.delta_p_exit) * scale;
    out.s[i] = s_node[i];
    out.residual[i] = res;
    out.max_abs = std::max(out.max_abs, std::abs(res));
    sum2 += res * res;
  }
  out.rms = std::sqrt(sum2 / static_cast<Real>(N));
  return out;
}

Real profile_error_against_cap(const CapillaryMeniscus& m) {
  const SphericalCap c = spherical_cap(m.contact_radius, m.delta_p_exit, m.gamma);
  Real e = 0.0;
  if (m.delta_p_exit == 0.0) {
    for (const Vec2& p : m.nodes) e = std::max(e, std::abs(p.z - m.contact_z));
    return e / m.contact_radius;
  }
  const Real absR = std::abs(c.sphere_radius);
  const Real sg = (c.sphere_radius > 0.0) ? 1.0 : -1.0;
  const Real z_centre = m.contact_z - sg * std::sqrt(std::max(0.0, absR * absR -
                                                              m.contact_radius * m.contact_radius));
  for (const Vec2& p : m.nodes) {
    const Real d = std::hypot(p.r, p.z - z_centre) - absR;
    e = std::max(e, std::abs(d));
  }
  return e / m.contact_radius;
}

Real profile_z_error_against_cap(const CapillaryMeniscus& m) {
  const SphericalCap c = spherical_cap(m.contact_radius, m.delta_p_exit, m.gamma);
  Real e = 0.0;
  for (const Vec2& p : m.nodes) {
    const Real r = std::min(p.r, m.contact_radius);
    e = std::max(e, std::abs((p.z - m.contact_z) - c.z_at_radius(r)));
  }
  return e / m.contact_radius;
}

// ---------------------------------------------------------------------------

void print(const CapillaryMeniscus& m, std::FILE* out) {
  std::fprintf(out, "Statischer Kapillarmeniskus (P3a, kein elektrisches Feld)\n");
  std::fprintf(out, "  Status            : %s\n", to_string(m.status));
  if (!m.message.empty()) std::fprintf(out, "  Begruendung       : %s\n", m.message.c_str());
  std::fprintf(out, "  Pinningradius a   : %.9e m  (aus phi_2/2)\n", m.contact_radius);
  std::fprintf(out, "  delta_p_exit      : %.9e Pa (= p_fluessig - p_vakuum)\n", m.delta_p_exit);
  std::fprintf(out, "  gamma             : %.9e N/m\n", m.gamma);
  std::fprintf(out, "  Pi = dp*a/gamma   : %.9e\n", m.Pi);
  if (!is_usable(m.status)) return;
  std::fprintf(out, "  Intervalle        : %d%s\n", m.n_intervals,
               m.discretisation_was_forced ? " (vorgegeben, Studienbetrieb)" : " (automatisch)");
  std::fprintf(out, "  Schaetzfehler     : %.3e (Aenderung zur halben Aufloesung)\n",
               m.estimated_relative_error);
  std::fprintf(out, "  Apexhoehe h       : %.9e m   (h/a = %.9e)\n", m.apex_height,
               m.apex_height / m.contact_radius);
  std::fprintf(out, "  Bogenlaenge       : %.9e m\n", m.arclength);
  std::fprintf(out, "  Rotationsflaeche  : %.9e m^2\n", m.revolved_area);
  std::fprintf(out, "  Rotationsvolumen  : %.9e m^3\n", m.revolved_volume);
  std::fprintf(out, "  psi(Kontaktlinie) : %.9e rad = %.4f deg\n", m.contact_tangent_angle,
               m.contact_tangent_angle * 180.0 / pi);
  std::fprintf(out, "  Kontaktradius ge-  \n");
  std::fprintf(out, "  troffen mit       : |r_end - a|/a = %.3e\n",
               std::abs(m.contact().r - m.contact_radius) / m.contact_radius);
}

}  // namespace es
