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
// Literature closed forms
// ---------------------------------------------------------------------------

Real hemisphere_balance_field(Real r, Real gamma) {
  if (!(r > 0.0)) return 0.0;
  return 2.0 * std::sqrt(gamma / (eps0 * r));
}

Real literature_onset_voltage_smith(Real r_c, Real d, Real gamma) {
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
  // Above the barrier the rate saturates at the attempt frequency, which the
  // Iribarne-Thomson form does not model.  Clamp so sweeps stay finite.
  const Real ex = std::exp(std::min(x, 30.0));
  return f.evap_prefactor * (kB * T / h_planck) * eps0 * E * ex;
}

Real field_for_current_density(Real j_target, const Fluid& f, Real T) {
  if (!(j_target > 0.0)) return 0.0;
  Real lo = 1e6, hi = 1e11;
  if (ion_current_density(hi, f, T) < j_target) return 0.0;
  if (ion_current_density(lo, f, T) > j_target) return lo;
  for (int i = 0; i < 200; ++i) {
    const Real mid = std::sqrt(lo * hi);
    if (ion_current_density(mid, f, T) < j_target) lo = mid; else hi = mid;
  }
  return std::sqrt(lo * hi);
}

Real characteristic_evaporation_field(const Fluid& f) {
  return 4.0 * pi * eps0 * f.dG_solvation * f.dG_solvation / (e * e * e);
}

void require_modelled_polarity(const BemSolver& bem) {
  const std::array<Real, 3>& V = bem.applied();
  require_modelled_polarity(V[static_cast<std::size_t>(Electrode::Emitter)] -
                            V[static_cast<std::size_t>(Electrode::Extractor)]);
}

void require_modelled_polarity(Real U) {
  if (U < 0.0) {
    throw NotImplementedInThisPhase(
        "Negative Emissionspolaritaet (Anionenemission)",
        "Phase P4 (Speziesbehandlung) und P5 (emittierender Betrieb)",
        "Modelliert ist ausschliesslich die Kationenreihe. Anionen haben andere Masse, "
        "andere Solvatationsenergie und eine andere Speziesverteilung; das Vorzeichen von "
        "E_n entscheidet ausserdem, ob ueberhaupt emittiert wird. Der Prototyp benutzte "
        "|E_n| und Kationenmassen und lieferte deshalb fuer beide Polaritaeten identische "
        "Stroeme -- eine Zahl, die fuer eine der beiden Polaritaeten falsch ist.");
  }
}

IonEmission integrate_ion_emission(const BemSolver& bem, const Fluid& f, Real T,
                                   bool include_wetted_metal) {
  require_modelled_polarity(bem);

  IonEmission out;
  const Mesh& m = bem.mesh();

  // A_eff = (int j dA)^2 / int j^2 dA -- a smooth functional of j, unlike the
  // element-quantised "99% area" it replaces.
  Real sum_j2A = 0.0;

  for (Index i = 0; i < m.size(); ++i) {
    const Element& el = m.elems[static_cast<std::size_t>(i)];
    const bool emits = (el.tag == Tag::FreeSurface) ||
                       (include_wetted_metal && el.tag == Tag::Emitter);
    if (!emits) continue;
    const Real E = std::abs(bem.En(i));
    const Real j = ion_current_density(E, f, T);
    out.current += j * el.area;
    sum_j2A += j * j * el.area;
    if (j > out.peak_j) { out.peak_j = j; out.peak_E = E; }
  }

  out.effective_area = (sum_j2A > 0.0) ? (out.current * out.current) / sum_j2A : 0.0;

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
  Real sum_mv = 0.0, sum_mq = 0.0, sum_m = 0.0;
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
  b.eta_polydispersity =
      (sum_m > 0.0 && b.beam_power > 0.0) ? b.thrust * b.thrust / (2.0 * sum_m * b.beam_power) : 0.0;
  return b;
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

void print_diagnostic_estimate(std::FILE* out, const Fluid& f, const IonEmission* ion,
                               const BeamFigures* fig) {
  std::fprintf(out,
      "\n=====================================================================\n"
      "NICHT GEKOPPELTE DIAGNOSTISCHE ABSCHAETZUNG -- kein Betriebspunkt\n"
      "=====================================================================\n"
      "Die folgenden Zahlen entstehen, indem die Iribarne-Thomson-Rate auf das\n"
      "Feld eines STATISCHEN, perfekt leitenden, NICHT emittierenden Meniskus\n"
      "angewandt wird. Im rein ionischen Betrieb ist die Stroemung\n"
      "viskositaetsdominiert und der Strom wird von der endlichen Leitfaehigkeit\n"
      "kontrolliert (Higuera 2008, Phys. Rev. E 77, 026308) -- das hier benutzte\n"
      "Feld ist also nicht das Feld, das tatsaechlich anlaege.\n"
      "Es handelt sich NICHT um eine Stromvorhersage. Das selbstkonsistente\n"
      "Modell ist fuer Phase P5 vorgesehen.\n"
      "=====================================================================\n");
  f.print(out);
  if (ion) {
    std::fprintf(out, "\nIonenverdampfung (Abschaetzung)\n");
    std::fprintf(out, "  Strom               : %10.4g A  (= %.4g nA)\n", ion->current,
                 ion->current * 1e9);
    std::fprintf(out, "  Spitzenstromdichte  : %10.4g A/m^2\n", ion->peak_j);
    std::fprintf(out, "  Feld am Maximum     : %10.4g V/m (= %.4f V/nm)\n", ion->peak_E,
                 ion->peak_E * 1e-9);
    std::fprintf(out, "  effektive Flaeche   : %10.4g m^2   (A_eff = (int j dA)^2 / int j^2 dA)\n",
                 ion->effective_area);
    std::fprintf(out, "  Ionenmassenstrom    : %10.4g kg/s\n", ion->mdot);
    std::fprintf(out,
        "  Empfindlichkeit     : d ln I / d dG = %.3g pro eV -- eine Aenderung von\n"
        "                        0,1 eV in dG aendert den Strom um Faktor %.2g.\n"
        "                        Literaturwerte fuer EMI-BF4 streuen 1,0-1,4 eV.\n",
        -1.0 / (constants::kB * 298.15) * constants::eV,
        std::exp(0.1 * constants::eV / (constants::kB * 298.15)));
  }
  if (fig) {
    std::fprintf(out, "\nDaraus abgeleitete Kennzahlen (erben die obige Einschraenkung)\n");
    std::fprintf(out, "  Strom               : %10.4g A\n", fig->current);
    std::fprintf(out, "  Massenstrom         : %10.4g kg/s\n", fig->mdot);
    std::fprintf(out, "  Schub               : %10.4g N  (= %.4g uN)\n", fig->thrust,
                 fig->thrust * 1e6);
    std::fprintf(out, "  spezifischer Impuls : %10.1f s\n", fig->Isp);
    std::fprintf(out, "  mittleres q/m       : %10.4g C/kg\n", fig->mean_qm);
    std::fprintf(out, "  eta_polydispersity  : %10.4f\n", fig->eta_polydispersity);
  }
}

void print_cone_jet_correlation(std::FILE* out, const Fluid& f, const ConeJetState& cj) {
  std::fprintf(out,
      "\n---------------------------------------------------------------------\n"
      "EMPIRISCHE KORRELATION (empirical = true) -- nicht an Geometrie, Feld\n"
      "oder Meniskus gekoppelt. Reine Formelauswertung.\n"
      "Quelle: Fernandez de la Mora & Loscertales (1994), JFM 260, 155-184.\n"
      "---------------------------------------------------------------------\n");
  std::fprintf(out, "  Flussrate           : %10.4g m^3/s (= %.4g nL/s)\n", cj.Q, cj.Q * 1e12);
  std::fprintf(out, "  Q / Q_min           : %10.2f %s\n", cj.q_over_qmin,
               cj.q_over_qmin < 1.0 ? "  <-- unterhalb der Stabilitaetsgrenze" : "");
  std::fprintf(out, "  Strom               : %10.4g A  (= %.4g nA)\n", cj.current,
               cj.current * 1e9);
  std::fprintf(out, "  Jetdurchmesser      : %10.4g m  (= %.4g nm)\n", cj.d_jet, cj.d_jet * 1e9);
  std::fprintf(out, "  Tropfendurchmesser  : %10.4g m  (= %.4g nm)\n", cj.d_droplet,
               cj.d_droplet * 1e9);
  std::fprintf(out, "  Tropfen-q/m         : %10.4g C/kg\n", cj.qm);
  std::fprintf(out, "  Fissilitaet q/q_R   : %10.3f %s\n", cj.fissility,
               cj.fissility > 0.8 ? "  <-- Coulomb-Fission zu erwarten" : "");
  if (cj.extrapolated)
    std::fprintf(out,
                 "  ACHTUNG: eps_r = %.1f < 40. Die Stromkorrelation wurde an\n"
                 "           Flüssigkeiten mit eps_r >~ 40 etabliert; f_current ist hier\n"
                 "           weit ausserhalb des Etablierungsbereichs extrapoliert.\n",
                 f.eps_r);
  std::fprintf(out,
               "  Diese Werte gehen NICHT in den Strahltransport ein. Die Kopplung ist\n"
               "  fuer Phase P6 vorgesehen.\n");
}

}  // namespace es
