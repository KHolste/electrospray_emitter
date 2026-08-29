// es_kinematics -- P4: the kinematic boundary condition, and what is missing.
//
//   es_kinematics <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS.  The kinematic condition dx/dt . n = u . n integrated on a
// PRESCRIBED velocity field, against two fields whose exact Lagrangian map is
// known in closed form.  That is KINEMATICS: the velocity comes from the
// caller, nothing is solved for it, and no force balance is evaluated anywhere.
//
// WHAT THIS RUN IS NOT.  There is no dynamic meniscus.  It is not implemented
// because the foundation is missing, not for lack of time -- the run writes the
// exact reason out and calls solve_dynamic_meniscus() to show that it fails
// closed.  There is NO mobility, NO artificial damping and NO stability
// statement anywhere.
//
// Exit code 2 means a declared check failed.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/surface_kinematics.hpp"

using namespace es;
using constants::pi;

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.empty()) {
    std::printf("es_kinematics -- P4: kinematische Randbedingung\n\n"
                "  es_kinematics <ausgabeverzeichnis> [key=value ...]\n");
    return 1;
  }
  Config cfg;
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  const Real R = cfg.num("kinematics.radius", 5.0e-6);
  const Real alpha = cfg.num("kinematics.alpha", 1.0e6);
  const Real T = cfg.num("kinematics.time", 5.0e-7);
  int exit_code = 0;

  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& s) {
    std::printf("%s\n", s.c_str());
    std::fprintf(log, "%s\n", s.c_str());
  };

  say("P4 -- kinematische Randbedingung auf vorgeschriebenem Feld");
  say("");
  say("DIE VORBEDINGUNG IST NICHT ERFUELLT.  Ein zeitabhaengiger Meniskus braucht ein");
  say("Geschwindigkeitsfeld mit freier Oberflaeche und einen Oberflaechenladungs-");
  say("transport; P3 liefert beides nicht.  Der dynamische Loeser schlaegt deshalb");
  say("geschlossen fehl:");
  try {
    solve_dynamic_meniscus();
    say("  FEHLER: er hat NICHT fehlgeschlagen.");
    exit_code = 2;
  } catch (const NotImplementedInThisPhase& e) {
    say(std::string("  ") + e.what());
  }
  say("");

  // --- the two validated cases ---------------------------------------------
  struct Case {
    const char* tag;
    VelocityField u;
    Vec2 (*exact)(Vec2, Real, Real);
    const char* note;
  };
  const Case cases[] = {
      {"dilatation", dilation_field(alpha), dilation_exact,
       "div u = 3 alpha; Volumen waechst wie e^{3 alpha t} -- eine BEKANNTE Aenderung"},
      {"squeeze", squeeze_field(alpha), squeeze_exact,
       "div u = 0 exakt; Volumen exakt erhalten, Form aendert sich stark"}};

  {
    std::FILE* f = std::fopen((outdir + "/shapes.csv").c_str(), "w");
    std::fprintf(f, "# Advektierte Oberflaeche gegen die exakte Lagrange-Abbildung.\n"
                    "# Das Feld ist VORGESCHRIEBEN; nichts wird dafuer geloest.\n");
    std::fprintf(f, "case,t_over_T,node,r_m,z_m,r_exact_m,z_exact_m\n");
    for (const Case& c : cases) {
      const SurfacePolyline s0 = hemisphere(R, 81);
      for (int frame = 0; frame <= 4; ++frame) {
        const Real t = T * static_cast<Real>(frame) / 4.0;
        const int steps = 200 * frame;
        AdvectionResult a;
        if (frame == 0) {
          a.surface = s0;
          a.status = StepStatus::Ok;
        } else {
          a = advect_surface(s0, c.u, t / steps, steps, KinematicMode::Lagrangian,
                             ContactLine::Free, 0.0);
        }
        if (a.status != StepStatus::Ok) {
          std::fprintf(f, "# %s bei t/T = %.2f: %s\n", c.tag, t / T, to_string(a.status));
          exit_code = 2;
          continue;
        }
        for (std::size_t k = 0; k < a.surface.nodes.size(); ++k) {
          const Vec2 x = a.surface.nodes[k];
          const Vec2 e = c.exact(s0.nodes[k], alpha, t);
          std::fprintf(f, "%s,%.4f,%zu,%.9e,%.9e,%.9e,%.9e\n", c.tag, t / T, k, x.r, x.z, e.r,
                       e.z);
        }
      }
    }
    std::fclose(f);
    say("  shapes.csv geschrieben");
  }

  {
    std::FILE* f = std::fopen((outdir + "/convergence.csv").c_str(), "w");
    std::fprintf(f, "# Zeitkonvergenz der Advektion gegen die exakte Abbildung, und die\n"
                    "# Volumenbilanz.  Fuer squeeze ist die exakte Volumenaenderung NULL,\n"
                    "# fuer dilatation ist sie e^{3 alpha T} -- die beiden Faelle trennen\n"
                    "# 'das Volumen stimmt' von 'die Form stimmt'.\n");
    std::fprintf(f, "case,steps,dt_s,status,shape_error_over_R,volume_change,"
                    "volume_change_exact,max_node_motion_fraction\n");
    for (const Case& c : cases) {
      const SurfacePolyline s0 = hemisphere(R, 81);
      const Real exact_change =
          (std::string(c.tag) == "dilatation") ? std::exp(3.0 * alpha * T) - 1.0 : 0.0;
      for (int n : {75, 150, 300, 600, 1200}) {
        const AdvectionResult a = advect_surface(s0, c.u, T / n, n, KinematicMode::Lagrangian,
                                                ContactLine::Free, 0.0);
        Real se = 0.0;
        for (std::size_t k = 0; k < a.surface.nodes.size(); ++k)
          se = std::max(se, norm(a.surface.nodes[k] - c.exact(s0.nodes[k], alpha,
                                                              a.surface.time)) / R);
        std::fprintf(f, "%s,%d,%.9e,%s,%.9e,%.9e,%.9e,%.9e\n", c.tag, n, T / n,
                     to_string(a.status), se, a.volume_change, exact_change,
                     a.max_node_motion_fraction);
      }
    }
    std::fclose(f);
    say("  convergence.csv geschrieben");
  }

  {
    // The redistribution: mesh motion, and how exactly it preserves the surface.
    std::FILE* f = std::fopen((outdir + "/redistribution.csv").c_str(), "w");
    std::fprintf(f, "# Die tangentiale Umverteilung ist NETZBEWEGUNG: sie schiebt Knoten\n"
                    "# entlang der Flaeche und kann die Flaeche deshalb nicht aendern.  Auf\n"
                    "# einem POLYGONZUG stimmt das nur bis auf die Diskretisierung, und wie\n"
                    "# genau, steht hier.\n");
    std::fprintf(f, "n_nodes,volume_lagrangian,volume_normal_only,relative_difference,"
                    "spacing_lagrangian,spacing_normal_only\n");
    for (std::size_t n : {41u, 81u, 161u, 321u, 641u}) {
      const SurfacePolyline sn = hemisphere(R, n);
      const AdvectionResult l = advect_surface(sn, squeeze_field(alpha), T / 900, 900,
                                              KinematicMode::Lagrangian, ContactLine::Free, 0.0);
      const AdvectionResult q = advect_surface(sn, squeeze_field(alpha), T / 900, 900,
                                              KinematicMode::NormalOnly, ContactLine::Free, 0.0);
      if (l.status != StepStatus::Ok || q.status != StepStatus::Ok) {
        std::fprintf(f, "# n = %zu: %s / %s\n", n, to_string(l.status), to_string(q.status));
        continue;
      }
      std::fprintf(f, "%zu,%.9e,%.9e,%.9e,%.9e,%.9e\n", n, l.volume_final, q.volume_final,
                   std::abs(q.volume_final - l.volume_final) / l.volume_final,
                   l.spacing_nonuniformity, q.spacing_nonuniformity);
    }
    std::fclose(f);
    say("  redistribution.csv geschrieben");
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_kinematics (P4)\nphase=P4\nstatus=infrastructure_only\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "radius_m=%.9e\nalpha_per_s=%.9e\ntime_s=%.9e\n", R, alpha, T);
    std::fprintf(f, "dynamic_solver=NotImplementedInThisPhase\n");
    std::fprintf(f, "max_node_motion=%.9e\n", kinematics::kMaxNodeMotion);
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  say("");
  say(exit_code == 0 ? "Alle deklarierten Pruefungen dieses Laufs bestanden."
                     : "MINDESTENS EINE DEKLARIERTE PRUEFUNG IST FEHLGESCHLAGEN.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_kinematics: %s\n", e.what());
  return 2;
}
