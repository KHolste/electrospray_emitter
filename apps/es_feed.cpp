// es_feed -- P1: the pressure budget at the exit plane.
//
//   es_feed <geometrie.cfg> [<feed.cfg> ...] <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS
//
//     delta_p_exit = (p_reservoir - p_vacuum) - delta_p_hydrostatic
//                                             - delta_p_viscous
//
// with the viscous term the laminar Hagen-Poiseuille loss of a straight,
// COMPLETELY FILLED circular feed channel.  It turns the free input
// delta_p_exit of P3a/P3b into a budget whose every term is either an explicit
// input or a closed form of explicit inputs, and it writes the resulting
// pressure into the same form the meniscus solver takes.
//
// WHAT THIS RUN IS NOT.  No flow solver: one closed-form resistance, steady,
// laminar, fully developed, incompressible.  No capillary rise and no moving
// contact line -- the channel is full, so it has no free surface of its own and
// no Young angle to insert; the only free surface is the pinned meniscus, whose
// capillary pressure P3a/P3b computes and which is NOT added here again.  The
// reservoir is a boundary condition, not a volume.  The reservoir pressure and
// the volumetric flow rate stay INPUTS; nothing here derives them.
//
// Exit code 2 means a declared check failed.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/capillary.hpp"
#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/feed.hpp"
#include "es/liquid.hpp"

using namespace es;
using constants::pi;

namespace {

constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();

void put(std::FILE* f, Real v) {
  if (std::isfinite(v))
    std::fprintf(f, ",%.9e", v);
  else
    std::fprintf(f, ",nan");
}

}  // namespace

// ===========================================================================

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.size() < 2) {
    std::printf(
        "es_feed -- P1: Druckhaushalt am Austritt\n\n"
        "  es_feed <geometrie.cfg> [<feed.cfg> ...] <ausgabeverzeichnis> [key=value ...]\n\n"
        "Gerechnet wird delta_p_exit = (p_res - p_vak) - rho g H - 8 mu L Q / (pi R^4).\n"
        "Kein Stroemungsloeser, kein kapillarer Aufstieg, keine bewegliche Kontaktlinie.\n"
        "Reservoirdruck und Volumenstrom bleiben Eingaben.\n");
    return 1;
  }

  Config cfg;
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) cfg.load(pos[i]);
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);
  for (std::size_t i = 0; i + 1 < pos.size(); ++i) {
    const std::filesystem::path src(pos[i]);
    std::filesystem::copy_file(src, std::filesystem::path(outdir) / src.filename(),
                               std::filesystem::copy_options::overwrite_existing);
  }

  // A contact angle anywhere is a modelling error in a filled channel.
  for (const char* key : {"feed.contact_angle_deg", "capillary.contact_angle_deg",
                          "wetting.contact_angle_deg", "liquid.contact_angle_deg"})
    if (cfg.has(key))
      throw std::runtime_error(
          std::string(key) +
          " ist gesetzt.  Der Zulaufkanal ist voll gefuellt und hat dort keine freie "
          "Oberflaeche; die einzige liegt an der gepinnten Austrittskante und wird von "
          "P3a/P3b gerechnet.");

  LiquidProperties liquid = liquid_data_by_name(cfg.str("liquid.name", "emi-bf4"));
  liquid.validate_or_throw();
  const Real a = 0.5 * cfg.num("device.phi_2", 10.0e-6);
  const Real gamma_over_a = liquid.gamma / a;

  FeedRequest base;
  base.mode = (cfg.str("feed.mode", "budget") == "direct") ? PressureMode::Direct
                                                           : PressureMode::Budget;
  base.delta_p_exit_direct = cfg.num("feed.delta_p_exit", 0.0);
  base.p_reservoir = cfg.num("feed.p_reservoir", 0.0);
  base.p_vacuum = cfg.num("feed.p_vacuum", 0.0);
  base.z_exit = cfg.num("feed.z_exit", 0.0);
  base.z_reservoir = cfg.num("feed.z_reservoir", 0.0);
  base.gravity_axial = cfg.num("feed.gravity_axial", 0.0);
  base.Q = cfg.num("feed.flow_rate", 0.0);
  base.channel.radius = cfg.num("feed.channel_radius", a);
  base.channel.length = cfg.num("feed.channel_length",
                                cfg.num("reservoir.feed_channel_length", 300.0e-6));

  int exit_code = 0;
  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& s) {
    std::printf("%s\n", s.c_str());
    std::fprintf(log, "%s\n", s.c_str());
    std::fflush(log);
  };

  say("P1 -- Druckhaushalt am Austritt");
  say("  Stoffdatensatz: " + liquid.substance + " [" + to_string(liquid.status) + "]");
  say("  a = " + std::to_string(a) + " m, gamma/a = " + std::to_string(gamma_over_a) + " Pa");
  say("  Kanal: R = " + std::to_string(base.channel.radius) + " m, L = " +
      std::to_string(base.channel.length) + " m");
  say("");
  say("WAS NICHT GERECHNET WIRD: kein Stroemungsloeser, kein kapillarer Aufstieg, keine");
  say("bewegliche Kontaktlinie, keine Reservoirentleerung.  p_reservoir und Q sind");
  say("Eingaben.  Der Kapillardruck des Meniskus steckt in P3a/P3b, nicht hier.");

  // --- the nominal point ----------------------------------------------------
  {
    const PressureBudget b = solve_pressure_budget(base, liquid, a);
    say("");
    b.print(stdout);
    b.print(log);
    if (!is_usable(b.status)) {
      say("  Der Nennpunkt ist nicht auswertbar; die Studien unten laufen trotzdem und");
      say("  markieren jeden nicht auswertbaren Punkt als nan.");
    }
    std::FILE* f = std::fopen((outdir + "/budget_nominal.csv").c_str(), "w");
    std::fprintf(f, "# Der Nennpunkt des Druckhaushalts, Term fuer Term.\n");
    std::fprintf(f, "term,value_Pa,note\n");
    std::fprintf(f, "driving,%.9e,p_reservoir - p_vacuum\n", b.driving);
    std::fprintf(f, "hydrostatic,%.9e,-rho g_z (z_exit - z_reservoir); wird abgezogen\n",
                 b.hydrostatic);
    std::fprintf(f, "viscous,%.9e,8 mu L Q / (pi R^4); wird abgezogen\n", b.viscous);
    std::fprintf(f, "delta_p_exit,%.9e,Ergebnis\n", b.delta_p_exit);
    std::fprintf(f, "gamma_over_a,%.9e,Kapillardruckskala\n", b.gamma_over_a);
    std::fprintf(f, "Pi,%.9e,delta_p_exit / (gamma/a)\n", b.Pi);
    std::fclose(f);
  }

  // --- 1. the budget against the flow rate ----------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/budget_vs_flow.csv").c_str(), "w");
    std::fprintf(f, "# Druckhaushalt gegen den Volumenstrom.  Q ist eine EINGABE; dieser\n"
                    "# Lauf berechnet ihn nicht.  Ein nicht auswertbarer Punkt steht als nan.\n");
    std::fprintf(f, "Q_m3_per_s,status,driving_Pa,hydrostatic_Pa,viscous_Pa,delta_p_exit_Pa,"
                    "Pi,reynolds,entrance_fraction,wall_shear_Pa,mean_velocity_m_per_s\n");
    const int n = 121;
    const Real Qmax = cfg.num("feed.flow_rate_max", 2.0e-14);
    for (int k = 0; k < n; ++k) {
      FeedRequest q = base;
      q.mode = PressureMode::Budget;
      q.Q = Qmax * static_cast<Real>(k) / static_cast<Real>(n - 1);
      const PressureBudget b = solve_pressure_budget(q, liquid, a);
      const bool ok = is_usable(b.status);
      std::fprintf(f, "%.9e,%s", q.Q, to_string(b.status));
      put(f, b.driving);
      put(f, b.hydrostatic);
      put(f, ok ? b.viscous : kNaN);
      put(f, ok ? b.delta_p_exit : kNaN);
      put(f, ok ? b.Pi : kNaN);
      put(f, b.reynolds);
      put(f, b.entrance_fraction);
      put(f, ok ? b.wall_shear_stress : kNaN);
      put(f, b.mean_velocity);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
    say("  budget_vs_flow.csv geschrieben");
  }

  // --- 2. sensitivity to channel length and radius --------------------------
  {
    std::FILE* f = std::fopen((outdir + "/sensitivity.csv").c_str(), "w");
    std::fprintf(f, "# Empfindlichkeit des viskosen Terms gegen Kanallaenge und -radius, bei\n"
                    "# festem Volumenstrom.  R_h ~ L und R_h ~ R^-4; das steht hier als\n"
                    "# Rechnung, nicht als Behauptung.\n");
    std::fprintf(f, "sweep,value_m,R_h_Pa_s_per_m3,viscous_Pa,viscous_over_gamma_over_a,"
                    "delta_p_exit_Pa,status\n");
    const Real Q = cfg.num("feed.flow_rate_reference", 1.0e-15);
    for (int k = 0; k < 61; ++k) {
      const Real L = 50.0e-6 * std::pow(20.0, static_cast<Real>(k) / 60.0);
      FeedRequest q = base;
      q.mode = PressureMode::Budget;
      q.Q = Q;
      q.channel.length = L;
      const PressureBudget b = solve_pressure_budget(q, liquid, a);
      std::fprintf(f, "length,%.9e", L);
      put(f, b.hydraulic_resistance);
      put(f, b.viscous);
      put(f, b.viscous / gamma_over_a);
      put(f, is_usable(b.status) ? b.delta_p_exit : kNaN);
      std::fprintf(f, ",%s\n", to_string(b.status));
    }
    for (int k = 0; k < 61; ++k) {
      const Real R = 1.0e-6 * std::pow(20.0, static_cast<Real>(k) / 60.0);
      FeedRequest q = base;
      q.mode = PressureMode::Budget;
      q.Q = Q;
      q.channel.radius = R;
      const PressureBudget b = solve_pressure_budget(q, liquid, a);
      std::fprintf(f, "radius,%.9e", R);
      put(f, b.hydraulic_resistance);
      put(f, b.viscous);
      put(f, b.viscous / gamma_over_a);
      put(f, is_usable(b.status) ? b.delta_p_exit : kNaN);
      std::fprintf(f, ",%s\n", to_string(b.status));
    }
    std::fclose(f);
    say("  sensitivity.csv geschrieben");
  }

  // --- 3. the pressure scales side by side ----------------------------------
  {
    std::FILE* f = std::fopen((outdir + "/pressure_scales.csv").c_str(), "w");
    std::fprintf(f, "# Einordnung der Druckskalen.  gamma/a ist die Skala, gegen die der\n"
                    "# Haushalt zaehlt; alles darunter aendert die Meniskusform kaum, alles\n"
                    "# darueber verlaesst den Bereich, in dem eine gepinnte statische Form\n"
                    "# ueberhaupt existiert (|delta_p| <= 2 gamma/a).\n");
    std::fprintf(f, "scale,value_Pa,relative_to_gamma_over_a,note\n");
    const Real Q = cfg.num("feed.flow_rate_reference", 1.0e-15);
    FeedRequest q = base;
    q.mode = PressureMode::Budget;
    q.Q = Q;
    const PressureBudget b = solve_pressure_budget(q, liquid, a);
    auto row = [&](const char* name, Real v, const char* note) {
      std::fprintf(f, "%s,%.9e,%.9e,%s\n", name, v, v / gamma_over_a, note);
    };
    row("gamma_over_a", gamma_over_a, "Kapillardruckskala");
    row("capillary_limit", 2.0 * gamma_over_a, "Grenze der gepinnten statischen Form");
    row("viscous", b.viscous, "8 mu L Q / (pi R^4) beim Referenzstrom");
    row("hydrostatic_1mm", liquid.rho * constants::g0 * 1.0e-3,
        "rho g h fuer 1 mm auf der Erde; im Orbit null");
    row("hydrostatic_bore", liquid.rho * constants::g0 * a,
        "rho g a -- die Bondzahl mal gamma/a");
    std::fclose(f);
    say("  pressure_scales.csv geschrieben");
    say("");
    say("  Einordnung: der viskose Term beim Referenzstrom ist " +
        std::to_string(b.viscous / gamma_over_a) + " gamma/a, die Hydrostatik ueber einen");
    say("  Bohrungsradius " + std::to_string(liquid.rho * constants::g0 * a / gamma_over_a) +
        " gamma/a.");
  }

  // --- 4. the coupling: what the meniscus solver makes of the budget --------
  {
    std::FILE* f = std::fopen((outdir + "/coupled_to_p3a.csv").c_str(), "w");
    std::fprintf(f, "# Der Haushalt, an den P3a-Kapillarloeser uebergeben.  KEINE neue\n"
                    "# Physik: es ist derselbe Loeser, nur mit einem delta_p_exit, das aus\n"
                    "# dem Haushalt kommt statt aus der Hand.  Ein Punkt ausserhalb\n"
                    "# |delta_p| <= 2 gamma/a hat keine gepinnte statische Form und steht\n"
                    "# als nan mit seinem Status.\n");
    std::fprintf(f, "p_reservoir_Pa,Q_m3_per_s,feed_status,delta_p_exit_Pa,Pi,"
                    "capillary_status,h_over_a,apex_curvature_1_per_m\n");
    const Real Q = cfg.num("feed.flow_rate_reference", 1.0e-15);
    for (int k = 0; k <= 40; ++k) {
      const Real p_res = -2.5 * gamma_over_a + 5.0 * gamma_over_a *
                                                   static_cast<Real>(k) / 40.0;
      FeedRequest q = base;
      q.mode = PressureMode::Budget;
      q.p_reservoir = p_res;
      q.Q = Q;
      const PressureBudget b = solve_pressure_budget(q, liquid, a);
      std::fprintf(f, "%.9e,%.9e,%s", p_res, Q, to_string(b.status));
      put(f, is_usable(b.status) ? b.delta_p_exit : kNaN);
      put(f, is_usable(b.status) ? b.Pi : kNaN);
      if (is_usable(b.status)) {
        CapillaryRequest cr;
        cr.delta_p_exit = b.delta_p_exit;
        cr.target_relative_accuracy = 1.0e-10;
        const CapillaryMeniscus m = solve_capillary_meniscus(a, 0.0, liquid, cr);
        std::fprintf(f, ",%s", to_string(m.status));
        put(f, is_usable(m.status) ? m.apex_height / a : kNaN);
        put(f, is_usable(m.status) ? 2.0 * b.delta_p_exit / liquid.gamma : kNaN);
      } else {
        std::fprintf(f, ",NotAttempted,nan,nan");
      }
      std::fprintf(f, "\n");
    }
    std::fclose(f);
    say("  coupled_to_p3a.csv geschrieben");
  }

  {
    std::FILE* f = std::fopen((outdir + "/parameters.csv").c_str(), "w");
    std::fprintf(f, "name,value_SI,unit,role\n");
    std::fprintf(f, "contact_radius,%.9e,m,gepinnter Kontaktradius a\n", a);
    std::fprintf(f, "gamma,%.9e,N/m,Oberflaechenspannung (%s)\n", liquid.gamma,
                 to_string(liquid.status));
    std::fprintf(f, "rho,%.9e,kg/m^3,Dichte (%s)\n", liquid.rho, to_string(liquid.status));
    std::fprintf(f, "mu,%.9e,Pa s,dynamische Viskositaet (%s)\n", liquid.documented_only.mu,
                 to_string(liquid.status));
    std::fprintf(f, "gamma_over_a,%.9e,Pa,Kapillardruckskala\n", gamma_over_a);
    std::fprintf(f, "channel_radius,%.9e,m,Zulaufkanal\n", base.channel.radius);
    std::fprintf(f, "channel_length,%.9e,m,Zulaufkanal\n", base.channel.length);
    std::fprintf(f, "gravity_axial,%.9e,m/s^2,Komponente von g entlang +z\n",
                 base.gravity_axial);
    std::fprintf(f, "reynolds_laminar,%.9e,-,Gueltigkeitsgrenze\n", feed::kReynoldsLaminar);
    std::fclose(f);
  }
  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "run=es_feed\nphase=P1\ncommit=%s\nliquid_status=%s\nexit_code=%d\n",
                 cfg.str("meta.commit", "unbekannt").c_str(), to_string(liquid.status),
                 exit_code);
    std::fclose(f);
  }
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_feed: %s\n", e.what());
  return 2;
}
