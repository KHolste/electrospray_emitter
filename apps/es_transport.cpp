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
  // =========================================================================
  {
    std::FILE* f = std::fopen((outdir + "/relaxation.csv").c_str(), "w");
    std::fprintf(f, "# Ladungsrelaxation.  tau = eps0 eps_r / sigma ist NUR berechenbar, wenn\n"
                    "# BEIDE Zahlen belegt sind.  Fuer EMI-BF4 fehlt eps_r (siehe\n"
                    "# docs/13_material_data.md), also steht der belegte Fall als nan.\n"
                    "# Der Fall 'unbelegt' benutzt einen ausdruecklich NICHT belegten Wert\n"
                    "# und ist ein Rechenbeispiel, keine Aussage ueber den Stoff.\n");
    std::fprintf(f, "case,eps_r,eps_r_status,sigma_S_per_m,sigma_status,tau_s,"
                    "process_time_s,ratio,verdict\n");

    // The capillary and visco-capillary times, from the SOURCED data.
    const Real t_cap = (rho_v.usable() && gam_v.usable())
                           ? std::sqrt(rho_v.value * a * a * a / gam_v.value)
                           : std::numeric_limits<Real>::quiet_NaN();
    const Real t_vis = gam_v.usable() ? mu * a / gam_v.value
                                      : std::numeric_limits<Real>::quiet_NaN();

    auto row = [&](const char* tag, const RelaxationVerdict& v) {
      std::fprintf(f, "%s", tag);
      put(f, v.eps_r > 0.0 ? v.eps_r : std::numeric_limits<Real>::quiet_NaN());
      std::fprintf(f, ",%s", to_string(v.eps_status));
      put(f, v.sigma > 0.0 ? v.sigma : std::numeric_limits<Real>::quiet_NaN());
      std::fprintf(f, ",%s", to_string(v.sigma_status));
      put(f, v.tau);
      put(f, v.process_time);
      put(f, v.ratio);
      std::fprintf(f, ",%s\n", to_string(v.limit));
    };
    const RelaxationVerdict sourced = judge_conductor_limit(mat, T, t_cap);
    row("belegt_kapillar", sourced);
    sourced.print(stdout);
    sourced.print(log);
    // An explicitly UNSOURCED what-if, labelled as such, so that the order of
    // magnitude is visible without being claimed.
    const Real eps_guess = cfg.num("transport.eps_r_unsourced", 12.8);
    if (sig_v.usable()) {
      row("unbelegt_kapillar",
          judge_conductor_limit_explicit(eps_guess, sig_v.value, t_cap));
      row("unbelegt_viskokapillar",
          judge_conductor_limit_explicit(eps_guess, sig_v.value, t_vis));
    }
    std::fclose(f);

    // The decay itself, drawn.
    if (sig_v.usable()) {
      std::FILE* g = std::fopen((outdir + "/decay.csv").c_str(), "w");
      std::fprintf(g, "# Geschlossene Ladungsrelaxation rho(t) = rho0 exp(-t/tau) mit einem\n"
                      "# ausdruecklich UNBELEGTEN eps_r.  Die Kurve zeigt die Form des\n"
                      "# Zerfalls und die Lage von tau, nicht einen Stoffwert.\n");
      std::fprintf(g, "t_over_tau,rho_over_rho0\n");
      const Real tau = charge_relaxation_time(eps_guess, sig_v.value);
      for (int k = 0; k <= 200; ++k) {
        const Real x = 6.0 * static_cast<Real>(k) / 200.0;
        std::fprintf(g, "%.9e,%.9e\n", x, relaxed_charge_density(1.0, x * tau, tau));
      }
      std::fclose(g);
    }

    // The time scales, side by side.
    std::FILE* h = std::fopen((outdir + "/time_scales.csv").c_str(), "w");
    std::fprintf(h, "# Zeitskalen nebeneinander.  Der Perfect-Conductor-Grenzfall ist ein\n"
                    "# VERHAELTNIS und keine Eigenschaft von sigma allein: welche Prozesszeit\n"
                    "# gilt, ist eine Modellentscheidung.\n");
    std::fprintf(h, "scale,value_s,status,note\n");
    std::fprintf(h, "tau_charge_relaxation,%.9e,%s,eps0 eps_r / sigma -- eps_r fehlt\n",
                 sourced.tau, to_string(sourced.limit));
    std::fprintf(h, "tau_charge_relaxation_unsourced,%.9e,unsourced,"
                    "mit eps_r = %.3g (NICHT belegt)\n",
                 sig_v.usable() ? charge_relaxation_time(eps_guess, sig_v.value)
                                : std::numeric_limits<Real>::quiet_NaN(),
                 eps_guess);
    std::fprintf(h, "t_capillary_inertial,%.9e,%s,sqrt(rho a^3 / gamma)\n", t_cap,
                 (rho_v.usable() && gam_v.usable()) ? "sourced" : "missing");
    std::fprintf(h, "t_visco_capillary,%.9e,%s,mu a / gamma\n", t_vis,
                 gam_v.usable() ? "sourced" : "missing");
    std::fclose(h);
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
