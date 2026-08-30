// es_transport -- P3: finite conductivity and feed flow, minimal honest step.
//
//   es_transport <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN COMPUTES
//
//   1. The steady axisymmetric pipe flow, reduced exactly to a scalar problem
//      for a straight pipe, solved on the existing FEM and compared with
//      Hagen-Poiseuille and with the hydraulic resistance P1 asserts.
//   2. The charge-transport contract: the relaxation time tau = eps/sigma, the
//      closed-form decay, the steady conduction current in a cylinder, and the
//      check that no current crosses a surface carrying the zero-flux
//      condition -- which is the condition a NON-EMITTING free surface has.
//   3. The perfect-conductor limit as a measurable ratio, judged from the
//      SOURCED material data, which for EMI-BF4 means it is not computable.
//
// WHAT THIS RUN IS NOT.  No general Stokes solver: no entrance flow, no
// pressure-velocity coupling, no curved geometry, no free surface.  No coupled
// finite-conductivity meniscus.  No emission, and therefore no emitted current.
// No time stepping of a shape.
//
// Exit code 2 means a declared check failed.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/feed.hpp"
#include "es/liquid.hpp"
#include "es/material_data.hpp"
#include "es/transport.hpp"

using namespace es;
using constants::eps0;
using constants::pi;

namespace {
void put(std::FILE* f, Real v) {
  if (std::isfinite(v))
    std::fprintf(f, ",%.9e", v);
  else
    std::fprintf(f, ",nan");
}

/// A quoted CSV field.  These notes contain commas and would otherwise split
/// the row; the field is quoted rather than the text mutilated.
std::string csv_quote(const std::string& in) {
  std::string out = "\"";
  for (char c : in) {
    if (c == '"')
      out += "\"\"";
    else if (c == '\n' || c == '\r')
      out += ' ';
    else
      out += c;
  }
  out += '"';
  return out;
}
}  // namespace

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.empty()) {
    std::printf("es_transport -- P3: endliche Leitfaehigkeit und Zulaufstroemung\n\n"
                "  es_transport <ausgabeverzeichnis> [key=value ...]\n");
    return 1;
  }
  Config cfg;
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  const MaterialDataset& mat = emibf4_sourced();
  const Real T = cfg.num("material.temperature", 298.15);
  const Real R = cfg.num("feed.channel_radius", 5.0e-6);
  const Real L = cfg.num("feed.channel_length", 300.0e-6);
  const Real a = 0.5 * cfg.num("device.phi_2", 10.0e-6);

  const MaterialValue mu_v = material_value(mat, PropertyKind::DynamicViscosity, T);
  const MaterialValue rho_v = material_value(mat, PropertyKind::Density, T);
  const MaterialValue gam_v = material_value(mat, PropertyKind::SurfaceTension, T);
  const MaterialValue sig_v = material_value(mat, PropertyKind::ElectricalConductivity, T);
  const MaterialValue eps_v = material_value(mat, PropertyKind::RelativePermittivity, T);

  int exit_code = 0;
  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& s) {
    std::printf("%s\n", s.c_str());
    std::fprintf(log, "%s\n", s.c_str());
    std::fflush(log);
  };

  say("P3 -- endliche Leitfaehigkeit und Zulaufstroemung, minimale ehrliche Stufe");
  say("");
  say("NICHT gerechnet: kein allgemeiner Stokes-Loeser (keine Einlaufstroemung, keine");
  say("Druck-Geschwindigkeits-Kopplung, keine gekruemmte Geometrie, keine freie");
  say("Oberflaeche), keine gekoppelte finite-conductivity-Meniskuskopplung, keine");
  say("Emission, keine Zeitintegration einer Form.");
  say("");
  for (const auto& pr : {std::pair<const char*, const MaterialValue*>{"mu", &mu_v},
                         {"rho", &rho_v},
                         {"gamma", &gam_v},
                         {"sigma", &sig_v},
                         {"eps_r", &eps_v}})
    say(std::string("  ") + pr.first + ": " + to_string(pr.second->status) +
        (pr.second->usable() ? " = " + std::to_string(pr.second->value) : ""));

  if (!mu_v.usable()) {
    say("  Ohne Viskositaet gibt es keine Stroemungsrechnung.  Geschlossener Abbruch.");
    std::fclose(log);
    return 2;
  }
  const Real mu = mu_v.value;

  // =========================================================================
  // 1.  Pipe flow: the profile, and the mesh study
  // =========================================================================
  {
    const Real dp = cfg.num("flow.pressure_drop", 1.0e4);
    const Real dpdz = -dp / L;
    std::FILE* f = std::fopen((outdir + "/pipe_profile.csv").c_str(), "w");
    std::fprintf(f, "# Geloestes Geschwindigkeitsprofil in der Mittelebene gegen die\n"
                    "# geschlossene Parabel.  Die Reduktion ist EXAKT fuer das gerade Rohr:\n"
                    "# Kontinuitaet erzwingt du_z/dz = 0, der konvektive Term verschwindet\n"
                    "# identisch, und die z-Impulsgleichung IST der achsensymmetrische\n"
                    "# Laplace-Operator.  Fuer alles andere gilt das nicht.\n");
    std::fprintf(f, "r_over_R,u_solved_m_per_s,u_closed_form_m_per_s,difference\n");
    const PipeFlowSolution s = solve_pipe_flow(R, L, mu, dpdz, 81, 9);
    const Real u0 = -dpdz * R * R / (4.0 * mu);
    for (int k = 0; k <= 80; ++k) {
      const Real x = static_cast<Real>(k) / 80.0;
      const Real want = u0 * (1.0 - x * x);
      const Real got = s.u[static_cast<std::size_t>(4 * 81 + k)];
      std::fprintf(f, "%.9e,%.9e,%.9e,%.9e\n", x, got, want, got - want);
    }
    std::fclose(f);

    std::FILE* g = std::fopen((outdir + "/pipe_convergence.csv").c_str(), "w");
    std::fprintf(g, "# Netzkonvergenz der geloesten Rohrstroemung gegen die geschlossene\n"
                    "# Form, und der daraus gewonnene hydraulische Widerstand gegen den,\n"
                    "# den P1 behauptet.  Die beiden teilen KEINEN Code.\n");
    std::fprintf(g, "nr,n_nodes,flow_rate_m3_per_s,flow_rate_closed_form,relative_error,"
                    "R_h_solved,R_h_p1,R_h_relative_error,profile_error\n");
    const Real Rh_p1 = 8.0 * mu * L / (pi * std::pow(R, 4.0));
    for (Index nr : {11, 21, 41, 81, 161}) {
      const PipeFlowSolution q = solve_pipe_flow(R, L, mu, dpdz, nr, 9);
      const Real Rh = dp / q.flow_rate;
      std::fprintf(g, "%lld,%lld", static_cast<long long>(nr),
                   static_cast<long long>(q.n_nodes));
      put(g, q.flow_rate);
      put(g, q.flow_rate_closed_form);
      put(g, q.flow_rate_error);
      put(g, Rh);
      put(g, Rh_p1);
      put(g, std::abs(Rh - Rh_p1) / Rh_p1);
      put(g, q.max_profile_error);
      std::fprintf(g, "\n");
    }
    std::fclose(g);
    say("");
    say("  Rohrstroemung: Volumenstrom trifft Hagen-Poiseuille auf " +
        std::to_string(s.flow_rate_error) + " relativ; der daraus gewonnene hydraulische");
    say("  Widerstand trifft die geschlossene Form von P1 auf " +
        std::to_string(std::abs(dp / s.flow_rate - Rh_p1) / Rh_p1) + " relativ.");
    if (!(s.flow_rate_error < 1.0e-2)) exit_code = 2;
  }

  // =========================================================================
  // 2.  Charge relaxation and the perfect-conductor limit
  //
  // WHICH permittivity belongs in tau_q is settled here rather than assumed.
  // An earlier version of this run demanded a single "DC permittivity" and
  // rejected the 1-18 GHz measurement as "not DC".  The free charge decays on
  // the time scale tau itself, so the permittivity that belongs in tau is the
  // one at f* = 1/(2 pi tau) -- which for this liquid falls in the low GHz
  // range, i.e. exactly where that measurement was taken.  The equation for
  // tau is therefore implicit, and it is solved on the MEASURED curve.
  //
  // A single eps_r is STILL not selected: no source states purity and water
  // content, and the single-value query still fails closed.  What is
  // established instead is a justified BAND plus the sensitivity over it, and
  // the verdict is taken at the corner of that band which is worst for the
  // approximation.
  // =========================================================================
  {
    const Real t_cap = (rho_v.usable() && gam_v.usable())
                           ? std::sqrt(rho_v.value * a * a * a / gam_v.value)
                           : std::numeric_limits<Real>::quiet_NaN();
    const Real t_vis = gam_v.usable() ? mu * a / gam_v.value
                                      : std::numeric_limits<Real>::quiet_NaN();

    const PermittivityBand band = permittivity_band(mat, T);
    const SelfConsistentRelaxation sc = self_consistent_relaxation(mat, T);
    band.print(stdout);  band.print(log);
    sc.print(stdout);    sc.print(log);

    // --- the five things that get called "the permittivity" ----------------
    {
      std::FILE* f = std::fopen((outdir + "/permittivity_concepts.csv").c_str(), "w");
      std::fprintf(f,
                   "# Fuenf verschiedene Groessen heissen bei einer leitfaehigen ionischen\n"
                   "# Fluessigkeit 'die Permittivitaet'.  Sie sind keine Varianten einer\n"
                   "# Zahl.  Welche in tau_q = eps0 eps_r / K gehoert, folgt aus der\n"
                   "# Herleitung und nicht aus einer Konvention: die freie Ladung zerfaellt\n"
                   "# auf der Zeitskala tau selbst, ihr Spektrum liegt also bei\n"
                   "# f* = 1/(2 pi tau), und dort ist eps_r abzulesen.\n");
      std::fprintf(f, "concept,enters_tau_q,explanation\n");
      const PermittivityConcept all[] = {PermittivityConcept::StaticApparentLowFrequency,
                                         PermittivityConcept::IntrinsicStatic,
                                         PermittivityConcept::FrequencyResolved,
                                         PermittivityConcept::ElectrodePolarisation,
                                         PermittivityConcept::DcConductivity};
      for (PermittivityConcept c : all)
        std::fprintf(f, "%s,%s,%s\n", to_string(c),
                     c == PermittivityConcept::FrequencyResolved ? "yes_at_f_star" : "no",
                     csv_quote(explain(c)).c_str());
      std::fclose(f);
      say("  permittivity_concepts.csv geschrieben");
    }

    // --- every permittivity datum, with its frequency and admissibility ----
    {
      std::FILE* f = std::fopen((outdir + "/permittivity_points.csv").c_str(), "w");
      std::fprintf(f,
                   "# Jeder Permittivitaetsmesspunkt bei dieser Temperatur.\n"
                   "# frequency_Hz = 0 heisst 'von der Quelle als statisch berichtet' -- fuer\n"
                   "# ionische Fluessigkeiten ist das kein Messwert bei null Hertz, sondern\n"
                   "# ein aus Mikrowellenspektren extrapolierter Grenzwert.\n"
                   "# admissible=no hiesse: unterhalb der Elektrodenpolarisationsschwelle und\n"
                   "# damit eine Eigenschaft der Messzelle statt der Fluessigkeit.  Der\n"
                   "# Datensatz enthaelt keinen solchen Punkt; die Spalte steht trotzdem da,\n"
                   "# damit die Pruefung sichtbar ist und nicht nur behauptet.\n");
      std::fprintf(f, "frequency_Hz,eps_r,uncertainty,concept,admissible,reference,method\n");
      const MaterialProperty* pp = mat.find(PropertyKind::RelativePermittivity);
      std::size_t n_written = 0;
      for (std::size_t k = 0; pp != nullptr && k < pp->n_sources; ++k) {
        const PropertySource& src = pp->sources[k];
        for (std::size_t j = 0; j < src.n_points; ++j) {
          const PropertyPoint& q = src.points[j];
          if (!q.ambient() || std::abs(q.T - T) > 2.0) continue;
          const bool below = q.frequency_resolved() &&
                             q.frequency_Hz < transport::kElectrodePolarisationFloor;
          const PermittivityConcept c =
              below ? PermittivityConcept::ElectrodePolarisation
                    : (q.frequency_resolved() ? PermittivityConcept::FrequencyResolved
                                              : PermittivityConcept::IntrinsicStatic);
          std::fprintf(f, "%.9e,%.9e,%.9e,%s,%s,%s,%s\n", q.frequency_Hz, q.value,
                       q.uncertainty, to_string(c), below ? "no" : "yes",
                       csv_quote(src.reference).c_str(), csv_quote(src.method).c_str());
          ++n_written;
        }
      }
      std::fclose(f);
      say("  permittivity_points.csv geschrieben (" + std::to_string(n_written) + " Punkte)");
    }

    // --- the implicit equation, as two curves whose intersection IS tau ----
    {
      std::FILE* f = std::fopen((outdir + "/self_consistency.csv").c_str(), "w");
      std::fprintf(f,
                   "# Die implizite Gleichung fuer tau, als zwei Kurven ueber der Frequenz:\n"
                   "#   tau_from_eps(f) = eps0 eps_r(f) / K   -- aus der GEMESSENEN Kurve\n"
                   "#   tau_from_f(f)   = 1 / (2 pi f)        -- die Definition von f*\n"
                   "# Ihr Schnittpunkt ist die Loesung.  Es wird nichts angepasst und keine\n"
                   "# Dispersionsfunktion erfunden; zwischen Messfrequenzen wird logarithmisch\n"
                   "# interpoliert und ausserhalb gar nicht erst gerechnet.\n");
      std::fprintf(f, "frequency_Hz,eps_r,tau_from_eps_s,tau_from_f_s,is_solution\n");
      if (sc.ok && sig_v.usable()) {
        const int n = 400;
        const Real lf0 = std::log10(sc.f_measured_min), lf1 = std::log10(sc.f_measured_max);
        for (int k = 0; k <= n; ++k) {
          const Real fq = std::pow(10.0, lf0 + (lf1 - lf0) * static_cast<Real>(k) / n);
          bool ex = false;
          const Real e = dielectric_permittivity_at(mat, T, fq, &ex);
          std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,no\n", fq, e,
                       constants::eps0 * e / sc.sigma, 1.0 / (2.0 * constants::pi * fq));
        }
        std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,yes\n", sc.f_star, sc.eps_r, sc.tau, sc.tau);
      }
      std::fclose(f);
      say("  self_consistency.csv geschrieben");
    }

    // --- the solution, and every corner of what it rests on ----------------
    {
      std::FILE* f = std::fopen((outdir + "/relaxation.csv").c_str(), "w");
      std::fprintf(f,
                   "# Ladungsrelaxation.  KEIN Einzelwert von eps_r ist ausgewaehlt -- keine\n"
                   "# Quelle nennt Reinheit und Wassergehalt -- und die Einzelwertabfrage\n"
                   "# meldet weiterhin MissingMaterialData.  Statt eines Ersatzwertes stehen\n"
                   "# hier die selbstkonsistente Loesung auf der GEMESSENEN Dispersionskurve\n"
                   "# und die vier Ecken des begruendeten Bandes aus eps_r und K.  Das Urteil\n"
                   "# wird an der fuer die Naeherung UNGUENSTIGSTEN Ecke gefaellt.\n");
      std::fprintf(f, "case,eps_r,eps_r_basis,sigma_S_per_m,tau_s,process_time_s,"
                      "process_name,ratio,verdict\n");
      auto write_row = [&](const char* tag, Real e, const char* basis, Real sg, Real tau,
                           Real tp, const char* tpname) {
        const Real ratio = tp / tau;
        std::fprintf(f, "%s,%.9e,%s,%.9e,%.9e,%.9e,%s,%.9e,%s\n", tag, e,
                     csv_quote(basis).c_str(), sg, tau, tp, tpname, ratio,
                     ratio > transport::kPerfectConductorMargin
                         ? "PerfectConductorJustified"
                         : "FiniteConductivityRequired");
      };
      const struct { const char* name; Real t; } procs[] = {{"t_capillary_inertial", t_cap},
                                                            {"t_visco_capillary", t_vis}};
      for (const auto& pr : procs) {
        if (!std::isfinite(pr.t)) continue;
        if (sc.ok)
          write_row("selbstkonsistent", sc.eps_r,
                    "eps_r(f*) auf der gemessenen Kurve, f* = 1/(2 pi tau)", sc.sigma, sc.tau,
                    pr.t, pr.name);
        const BandedRelaxationVerdict bv = judge_conductor_limit_over_band(mat, T, pr.t);
        bv.print(stdout);
        bv.print(log);
        if (bv.limit != ConductorLimit::PerfectConductorJustified) exit_code = 2;
        const Real eps[2] = {bv.eps_lo, bv.eps_hi};
        const Real sig[2] = {bv.sigma_lo, bv.sigma_hi};
        const char* en[2] = {"eps_lo", "eps_hi"};
        const char* sn[2] = {"K_lo", "K_hi"};
        for (int ie = 0; ie < 2; ++ie)
          for (int is = 0; is < 2; ++is) {
            char tag[64];
            std::snprintf(tag, sizeof tag, "Bandecke_%s_%s", en[ie], sn[is]);
            write_row(tag, eps[ie], "Ecke des begruendeten Bandes", sig[is],
                      constants::eps0 * eps[ie] / sig[is], pr.t, pr.name);
          }
      }
      std::fclose(f);
      say("  relaxation.csv geschrieben");
    }

    // --- the decay itself, on the self-consistent tau ----------------------
    if (sc.ok) {
      std::FILE* g = std::fopen((outdir + "/decay.csv").c_str(), "w");
      std::fprintf(g, "# rho(t) = rho0 exp(-t/tau) mit dem SELBSTKONSISTENTEN tau.\n"
                      "# Die Kurve zeigt die Form des Zerfalls und die Lage von tau.\n");
      std::fprintf(g, "t_over_tau,rho_over_rho0\n");
      for (int k = 0; k <= 200; ++k) {
        const Real x = 6.0 * static_cast<Real>(k) / 200.0;
        std::fprintf(g, "%.9e,%.9e\n", x, relaxed_charge_density(1.0, x * sc.tau, sc.tau));
      }
      std::fclose(g);
    }

    // --- the time scales, with tau as a BAND and not as a point ------------
    {
      std::FILE* h = std::fopen((outdir + "/time_scales.csv").c_str(), "w");
      std::fprintf(h, "# Zeitskalen nebeneinander.  Der Perfect-Conductor-Grenzfall ist ein\n"
                      "# VERHAELTNIS und keine Eigenschaft von K allein: welche Prozesszeit\n"
                      "# gilt, ist eine Modellentscheidung.  tau steht als BAND, weil kein\n"
                      "# einzelner eps_r-Wert belegt ist; lo/hi sind die Ecken des\n"
                      "# begruendeten Bandes, value_s ist der selbstkonsistente Wert.\n");
      std::fprintf(h, "scale,value_s,lo_s,hi_s,status,note\n");
      const BandedRelaxationVerdict bv = judge_conductor_limit_over_band(mat, T, t_cap);
      char note[512];
      std::snprintf(note, sizeof note,
                    "eps0 eps_r(f*)/K mit f* = %.4g Hz; Band ueber eps_r %.3g..%.3g und "
                    "K %.3g..%.3g S/m",
                    sc.ok ? sc.f_star : 0.0, bv.eps_lo, bv.eps_hi, bv.sigma_lo, bv.sigma_hi);
      std::fprintf(h, "tau_charge_relaxation,%.9e,%.9e,%.9e,%s,%s\n",
                   sc.ok ? sc.tau : std::numeric_limits<Real>::quiet_NaN(), bv.tau_min,
                   bv.tau_max, sc.ok ? "self_consistent_band" : "missing",
                   csv_quote(note).c_str());
      std::fprintf(h, "t_capillary_inertial,%.9e,nan,nan,%s,\"sqrt(rho a^3 / gamma)\"\n", t_cap,
                   (rho_v.usable() && gam_v.usable()) ? "sourced" : "missing");
      std::fprintf(h, "t_visco_capillary,%.9e,nan,nan,%s,\"mu a / gamma\"\n", t_vis,
                   gam_v.usable() ? "sourced" : "missing");
      std::fclose(h);
      say("  time_scales.csv geschrieben");
    }

    // The single-value query must STILL fail closed.  Written out and checked,
    // so that nobody reads the band as if a value had been selected.
    const RelaxationVerdict single = judge_conductor_limit(mat, T, t_cap);
    single.print(stdout);
    single.print(log);
    if (single.limit != ConductorLimit::MissingMaterialData) {
      say("  FEHLER: die Einzelwertabfrage haette geschlossen fehlschlagen muessen.");
      exit_code = 2;
    }
  }

  // =========================================================================
  // 3.  Steady conduction in a cylinder
  // =========================================================================
  if (sig_v.usable()) {
    const Real sigma = sig_v.value;
    const Real V = cfg.num("transport.test_voltage", 1.0);
    std::FILE* f = std::fopen((outdir + "/conduction.csv").c_str(), "w");
    std::fprintf(f, "# Stationaerer Leitungsstrom im Zylinder: div(sigma grad phi) = 0 ist\n"
                    "# dieselbe Gleichung wie die Elektrostatik, die Knotenreaktion ist dann\n"
                    "# ein STROM in Ampere.  Die Mantelflaeche traegt die natuerliche\n"
                    "# Bedingung j.n = 0 -- genau die, die eine NICHT emittierende freie\n"
                    "# Oberflaeche tragen muss, weil die Ladung sonst nirgendwohin koennte.\n");
    std::fprintf(f, "nr,nz,n_nodes,current_A,current_closed_form_A,relative_error,"
                    "resistance_Ohm,resistance_closed_form_Ohm,lateral_leakage,"
                    "potential_error,fem_residual\n");
    Real worst_leak = 0.0;
    for (Index n : {11, 21, 41}) {
      const ConductionSolution c = solve_cylinder_conduction(R, L, sigma, V, n, 2 * n - 1);
      std::fprintf(f, "%lld,%lld,%lld", static_cast<long long>(n),
                   static_cast<long long>(2 * n - 1), static_cast<long long>(c.n_nodes));
      put(f, c.current);
      put(f, c.current_closed_form);
      put(f, c.current_error);
      put(f, c.resistance);
      put(f, c.resistance_closed_form);
      put(f, c.lateral_leakage);
      put(f, c.max_potential_error);
      put(f, c.fem_residual);
      std::fprintf(f, "\n");
      worst_leak = std::max(worst_leak, c.lateral_leakage);
    }
    std::fclose(f);
    say("");
    say("  Leitungsstrom: groesste seitliche Leckstromdichte relativ zur axialen: " +
        std::to_string(worst_leak) + ".  Ohne Emission darf sie null sein und ist es.");
    if (!(worst_leak < 1.0e-10)) exit_code = 2;
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_transport (P3)\nphase=P3\ncommit=%s\n",
                 cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "channel_radius_m=%.9e\nchannel_length_m=%.9e\ncontact_radius_m=%.9e\n", R,
                 L, a);
    std::fprintf(f, "T_K=%.9e\n", T);
    std::fprintf(f, "mu_status=%s\nrho_status=%s\ngamma_status=%s\nsigma_status=%s\n"
                    "eps_r_status=%s\n",
                 to_string(mu_v.status), to_string(rho_v.status), to_string(gam_v.status),
                 to_string(sig_v.status), to_string(eps_v.status));
    if (mu_v.usable()) std::fprintf(f, "mu=%.9e\n", mu_v.value);
    if (sig_v.usable()) std::fprintf(f, "sigma=%.9e\n", sig_v.value);
    std::fprintf(f, "perfect_conductor_margin=%.9e\n", transport::kPerfectConductorMargin);
    std::fprintf(f, "eps_r_unsourced=%.9e\n", cfg.num("transport.eps_r_unsourced", 12.8));
    std::fprintf(f, "status=validated_subset\n");
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  say("");
  say(exit_code == 0 ? "Alle deklarierten Pruefungen dieses Laufs bestanden."
                     : "MINDESTENS EINE DEKLARIERTE PRUEFUNG IST FEHLGESCHLAGEN.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_transport: %s\n", e.what());
  return 2;
}
