#include "es/emission.hpp"

#include <algorithm>
#include <cmath>

#include "es/constants.hpp"

namespace es {

using constants::e;
using constants::eps0;
using constants::g0;
using constants::h_planck;
using constants::kB;
using constants::pi;
using constants::taylor_angle;

// ---------------------------------------------------------------------------
// Onset
// ---------------------------------------------------------------------------

Real onset_field_hemisphere(Real r, Real gamma) {
  if (!(r > 0.0)) return 0.0;
  return 2.0 * std::sqrt(gamma / (eps0 * r));
}

Real onset_voltage_taylor(Real r_c, Real d, Real gamma) {
  if (!(r_c > 0.0) || !(d > r_c)) return 0.0;
  return std::sqrt(2.0 * gamma * r_c * std::cos(taylor_angle) / eps0) * std::log(4.0 * d / r_c);
}

Real rayleigh_charge(Real d, Real gamma) {
  const Real R = 0.5 * d;
  if (!(R > 0.0)) return 0.0;
  return 8.0 * pi * std::sqrt(eps0 * gamma * R * R * R);
}

// ---------------------------------------------------------------------------
// Ion evaporation
// ---------------------------------------------------------------------------

Real schottky_lowering(Real E) {
  if (E <= 0.0) return 0.0;
  return std::sqrt(e * e * e * E / (4.0 * pi * eps0));
}

Real ion_current_density(Real E, const Fluid& f, Real T) {
  if (E <= 0.0 || T <= 0.0) return 0.0;
  const Real barrier = f.dG_solvation - schottky_lowering(E);
  const Real x = -barrier / (kB * T);
  // Guard the exponential: above the barrier the rate saturates at the attempt
  // frequency, which the Iribarne-Thomson form does not model.  Clamping keeps
  // parameter sweeps finite instead of producing inf.
  const Real ex = std::exp(std::min(x, 30.0));
  return f.evap_prefactor * (kB * T / h_planck) * eps0 * E * ex;
}

Real field_for_current_density(Real j_target, const Fluid& f, Real T) {
  if (!(j_target > 0.0)) return 0.0;
  Real lo = 1e6, hi = 1e11;
  if (ion_current_density(hi, f, T) < j_target) return 0.0;
  if (ion_current_density(lo, f, T) > j_target) return lo;
  for (int i = 0; i < 200; ++i) {
    const Real mid = std::sqrt(lo * hi);  // geometric bisection: spans decades
    if (ion_current_density(mid, f, T) < j_target) lo = mid; else hi = mid;
  }
  return std::sqrt(lo * hi);
}

Real characteristic_evaporation_field(const Fluid& f) {
  // G(E) = dG  =>  E = 4 pi eps0 dG^2 / e^3
  return 4.0 * pi * eps0 * f.dG_solvation * f.dG_solvation / (e * e * e);
}

IonEmission integrate_ion_emission(const BemSolver& bem, const Fluid& f, Real T,
                                   bool include_wetted_metal) {
  IonEmission out;
  const Mesh& m = bem.mesh();
  std::vector<std::pair<Real, Real>> contrib;  // (current, area), descending
  contrib.reserve(static_cast<std::size_t>(m.size()));

  for (Index i = 0; i < m.size(); ++i) {
    const Element& el = m.elems[static_cast<std::size_t>(i)];
    const bool emits = (el.tag == Tag::FreeSurface) ||
                       (include_wetted_metal && el.tag == Tag::Emitter);
    if (!emits) continue;
    // Ions leave along the outward normal, so only an outward-pulling field
    // (positive sigma for a positively biased emitter) drives emission.  Use
    // |E_n|: the polarity is set by the applied voltage sign, not by the model.
    const Real E = std::abs(bem.En(i));
    const Real j = ion_current_density(E, f, T);
    const Real dI = j * el.area;
    out.current += dI;
    if (j > out.peak_j) { out.peak_j = j; out.peak_E = E; }
    contrib.emplace_back(dI, el.area);
  }

  // "Emitting area" = the smallest area collecting 99% of the current.  Field
  // emission is so steep that this is typically a tiny fraction of the wetted
  // surface -- worth reporting explicitly, because it drives the local heat
  // load and the evaporative cooling of the apex.
  std::sort(contrib.begin(), contrib.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });
  Real acc = 0.0;
  for (const auto& c : contrib) {
    if (acc >= 0.99 * out.current) break;
    acc += c.first;
    out.emitting_area += c.second;
  }

  // Mass flow of the emitted clusters, from the mean q/m.
  const Real qm = f.qm_cluster();
  out.mdot = (qm > 0.0) ? out.current / qm : 0.0;
  return out;
}

// ---------------------------------------------------------------------------
// Cone-jet
// ---------------------------------------------------------------------------

ConeJetState cone_jet(const Fluid& f, Real Q, const ConeJetModel& model) {
  ConeJetState s;
  s.Q = Q;
  if (!(Q > 0.0)) return s;

  s.current = model.f_current * std::sqrt(f.gamma * f.K * Q / f.eps_r);
  s.d_jet = model.d_prefactor * std::cbrt(Q * eps0 * f.eps_r / f.K);
  s.d_droplet = model.jet_to_drop * s.d_jet;
  s.mdot = f.rho * Q;
  s.qm = (s.mdot > 0.0) ? s.current / s.mdot : 0.0;

  const Real vol = pi / 6.0 * s.d_droplet * s.d_droplet * s.d_droplet;
  s.q_droplet = s.qm * f.rho * vol;
  const Real qR = rayleigh_charge(s.d_droplet, f.gamma);
  s.fissility = (qR > 0.0) ? s.q_droplet / qR : 0.0;

  const Real qmin = f.q_min();
  s.q_over_qmin = (qmin > 0.0) ? Q / qmin : 0.0;
  s.extrapolated = (f.eps_r < 40.0);
  return s;
}

// ---------------------------------------------------------------------------
// Beam figures
// ---------------------------------------------------------------------------

BeamFigures beam_figures(const std::vector<Species>& species, Real V_accel) {
  BeamFigures b;
  if (V_accel <= 0.0) return b;
  Real sum_mv = 0.0;   // thrust
  Real sum_mq = 0.0;   // current
  Real sum_m = 0.0;    // mass flow
  for (const Species& s : species) {
    if (s.mdot <= 0.0 || s.qm <= 0.0) continue;
    const Real v = std::sqrt(2.0 * s.qm * V_accel);
    sum_mv += s.mdot * v;
    sum_mq += s.mdot * s.qm;
    sum_m += s.mdot;
  }
  b.thrust = sum_mv;
  b.current = sum_mq;
  b.mdot = sum_m;
  b.beam_power = b.current * V_accel;
  b.mean_qm = (sum_m > 0.0) ? sum_mq / sum_m : 0.0;
  b.Isp = (sum_m > 0.0) ? b.thrust / (sum_m * g0) : 0.0;
  // eta_pol = F^2 / (2 mdot P).  Exactly 1 for a single species; the deficit is
  // the energy wasted accelerating a mixture through one common potential.
  b.eta_polydispersity =
      (sum_m > 0.0 && b.beam_power > 0.0) ? b.thrust * b.thrust / (2.0 * sum_m * b.beam_power) : 0.0;
  return b;
}

void print_operating_point(std::FILE* out, const Fluid& f, const ConeJetState* cj,
                           const IonEmission* ion, const BeamFigures* fig) {
  f.print(out);
  if (cj) {
    std::fprintf(out, "\ncone-jet (droplet) mode\n");
    std::fprintf(out, "  flow rate           : %10.4g m^3/s (= %.4g nL/s)\n", cj->Q, cj->Q * 1e12);
    std::fprintf(out, "  Q / Q_min           : %10.2f %s\n", cj->q_over_qmin,
                 cj->q_over_qmin < 1.0 ? "  <-- below the stability floor" : "");
    std::fprintf(out, "  current             : %10.4g A  (= %.4g nA)\n", cj->current,
                 cj->current * 1e9);
    std::fprintf(out, "  jet diameter        : %10.4g m  (= %.4g nm)\n", cj->d_jet,
                 cj->d_jet * 1e9);
    std::fprintf(out, "  droplet diameter    : %10.4g m  (= %.4g nm)\n", cj->d_droplet,
                 cj->d_droplet * 1e9);
    std::fprintf(out, "  droplet q/m         : %10.4g C/kg\n", cj->qm);
    std::fprintf(out, "  fissility q/q_R     : %10.3f %s\n", cj->fissility,
                 cj->fissility > 0.8 ? "  <-- Coulomb fission expected" : "");
    if (cj->extrapolated)
      std::fprintf(out,
                   "  NOTE: eps_r = %.1f < 40, so f_current is extrapolated far outside the\n"
                   "        range where the correlation was established.  Calibrate against\n"
                   "        your own I(Q) data before believing the absolute current.\n",
                   f.eps_r);
  }
  if (ion) {
    std::fprintf(out, "\nion evaporation (Iribarne-Thomson)\n");
    std::fprintf(out, "  ion current         : %10.4g A  (= %.4g nA)\n", ion->current,
                 ion->current * 1e9);
    std::fprintf(out, "  peak current density: %10.4g A/m^2\n", ion->peak_j);
    std::fprintf(out, "  field at the peak   : %10.4g V/m (= %.3f V/nm)\n", ion->peak_E,
                 ion->peak_E * 1e-9);
    std::fprintf(out, "  99%%-current area    : %10.4g m^2\n", ion->emitting_area);
    std::fprintf(out, "  ion mass flow       : %10.4g kg/s\n", ion->mdot);
  }
  if (fig) {
    std::fprintf(out, "\nbeam figures of merit\n");
    std::fprintf(out, "  total current       : %10.4g A  (= %.4g nA)\n", fig->current,
                 fig->current * 1e9);
    std::fprintf(out, "  total mass flow     : %10.4g kg/s\n", fig->mdot);
    std::fprintf(out, "  thrust              : %10.4g N  (= %.4g uN)\n", fig->thrust,
                 fig->thrust * 1e6);
    std::fprintf(out, "  specific impulse    : %10.1f s\n", fig->Isp);
    std::fprintf(out, "  beam power          : %10.4g W\n", fig->beam_power);
    std::fprintf(out, "  mean q/m            : %10.4g C/kg\n", fig->mean_qm);
    std::fprintf(out, "  eta_polydispersity  : %10.4f\n", fig->eta_polydispersity);
  }
}

}  // namespace es
