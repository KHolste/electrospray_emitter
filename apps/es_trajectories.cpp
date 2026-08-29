// es_trajectories -- P7: electrostatic particle transport and the balance.
//
//   es_trajectories <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS.  A TRANSPORT RESPONSE: a PRESCRIBED launch distribution is
// pushed through the field of the P6 problem, once without and once with a
// PRESCRIBED space charge, and the trajectories, impact sites, energies and the
// current balance are written out.
//
// WHAT IT IS NOT.  Not a current prediction.  P5 is blocked, there is no
// physical particle source, and the launched current is an input.  No droplet
// beam.  No self-consistent space-charge loop: the particles do not update the
// charge while they fly, which is what makes the with/without comparison a
// clean difference instead of an entangled one.
//
// Exit code 2 means a declared check failed.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/particle_transport.hpp"
#include "es/space_charge.hpp"

using namespace es;
using constants::e;

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
    std::printf("es_trajectories -- P7: Teilchentransport\n\n"
                "  es_trajectories <ausgabeverzeichnis> [key=value ...]\n");
    return 1;
  }
  Config cfg;
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  const Real R = cfg.num("box.radius", 2.0e-5);
  const Real Z = cfg.num("box.length", 1.0e-4);
  const Real V = cfg.num("box.voltage", -1000.0);
  const Real ap = cfg.num("box.aperture", 0.35) * R;
  const Index nr = cfg.integer("box.nr", 61), nz = cfg.integer("box.nz", 121);
  int exit_code = 0;

  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& s) {
    std::printf("%s\n", s.c_str());
    std::fprintf(log, "%s\n", s.c_str());
  };
  say("P7 -- elektrostatischer Teilchentransport");
  say("");
  say("TRANSPORTANTWORT AUF EINE VORGESCHRIEBENE STARTVERTEILUNG.  P5 ist");
  say("blockiert; es gibt keine physikalische Teilchenquelle.  Der gestartete");
  say("Strom ist eine EINGABE, und keine Zahl dieses Laufs ist eine");
  say("Stromvorhersage.  Kein Tropfenstrahl, keine selbstkonsistente");
  say("Raumladungsschleife.");

  // --- the mesh, the electrodes and the two fields --------------------------
  QuadMesh mesh;
  mesh.nr = nr;
  mesh.nz = nz;
  mesh.nodes.resize(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nz));
  for (Index j = 0; j < nz; ++j)
    for (Index i = 0; i < nr; ++i)
      mesh.nodes[static_cast<std::size_t>(j * nr + i)] =
          Vec2{R * static_cast<Real>(i) / static_cast<Real>(nr - 1),
               Z * static_cast<Real>(j) / static_cast<Real>(nz - 1)};
  std::vector<Real> eps_r(static_cast<std::size_t>(mesh.n_cells()), 1.0);
  std::vector<char> active(static_cast<std::size_t>(mesh.n_cells()), 1);
  std::vector<char> fixed(static_cast<std::size_t>(mesh.n_nodes()), 0);
  std::vector<Real> fixed_value(static_cast<std::size_t>(mesh.n_nodes()), 0.0);
  for (Index i = 0; i < nr; ++i) {
    fixed[static_cast<std::size_t>(mesh.node(i, 0))] = 1;          // emitter, 0 V
    const Index c = mesh.node(i, nz - 1);
    fixed[static_cast<std::size_t>(c)] = 1;                        // extractor, V
    fixed_value[static_cast<std::size_t>(c)] = V;
  }
  // The side is NOT held at a potential.  It carries the natural zero-flux
  // condition, which makes the field the clean parallel-plate one the transport
  // test needs.  A first version grounded the whole wall; with R/Z = 0.2 that
  // screens the applied field to a few thousand V/m near the emitter and the
  // particles crept instead of flying -- correct for THAT geometry, and a
  // useless one for testing transport.  The wall is still a SURFACE a particle
  // can hit; that is the domain classifier's job, not the field's.

  // A PRESCRIBED space charge in the beam channel.
  std::vector<Macroparticle> cloud;
  {
    const int np = cfg.integer("charge.n", 600);
    const Real q = cfg.num("charge.total", 3.0e-14) / static_cast<Real>(np);
    const Real rb = cfg.num("charge.radius", 0.3) * R;
    auto vdc = [](int n) {
      Real x = 0.0, b = 0.5;
      while (n > 0) { x += b * (n % 2); n /= 2; b *= 0.5; }
      return x;
    };
    for (int k = 0; k < np; ++k)
      cloud.push_back({{rb * std::sqrt(vdc(k + 1)),
                        Z * (static_cast<Real>(k) + 0.5) / static_cast<Real>(np)},
                       q});
  }
  const DepositionResult dep = deposit(mesh, cloud);
  const SpaceChargeSolution sol =
      solve_with_space_charge(mesh, eps_r, active, fixed, fixed_value, dep.node_charge, {});
  say("");
  {
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "  vorgeschriebene Raumladung: %.4e C, Erhaltungsfehler %.3e; "
                  "Potentialverschiebung %.4f V",
                  dep.total_particles, dep.conservation_error, sol.max_potential_shift);
    say(buf);
  }
  if (!(dep.conservation_error < 1.0e-13)) exit_code = 2;

  // --- the species, explicit and with signs ---------------------------------
  const TransportSpecies cation{"EMI+ (Testspezies)", +e, 111.17e-3 / 6.02214076e23};
  const TransportSpecies anion{"BF4- (Testspezies)", -e, 86.81e-3 / 6.02214076e23};
  const TransportDomain domain{R, Z, ap};

  // --- the prescribed launch distribution -----------------------------------
  std::vector<TracedParticle> launch;
  {
    const int np = cfg.integer("beam.n", 41);
    const Real I = cfg.num("beam.current", 1.0e-7);
    const Real rb = cfg.num("beam.radius", 0.5) * R;
    for (int k = 0; k < np; ++k) {
      TracedParticle p;
      const Real t = (static_cast<Real>(k) + 0.5) / static_cast<Real>(np);
      p.x = {rb * t, 1.0e-7};
      p.v = {0.0, 0.0};
      p.current = I / static_cast<Real>(np);
      launch.push_back(p);
    }
  }
  const Real dt = cfg.num("beam.dt", 2.0e-13);
  const int steps = cfg.integer("beam.max_steps", 200000);

  struct Case {
    const char* tag;
    const std::vector<Real>* phi;
    const TransportSpecies* sp;
    Real polarity_voltage;
  };

  // The anion case needs the opposite applied voltage, so it gets its own
  // solve.  The two polarities are never taken from one field.
  std::vector<Real> fixed_value_neg = fixed_value;
  for (Index i = 0; i < nr; ++i)
    fixed_value_neg[static_cast<std::size_t>(mesh.node(i, nz - 1))] = -V;
  const SpaceChargeSolution sol_neg = solve_with_space_charge(
      mesh, eps_r, active, fixed, fixed_value_neg, dep.node_charge, {});

  const Case cases[] = {
      {"kation_ohne_raumladung", &sol.phi_no_charge, &cation, V},
      {"kation_mit_raumladung", &sol.phi, &cation, V},
      {"anion_ohne_raumladung", &sol_neg.phi_no_charge, &anion, -V},
  };

  std::FILE* ft = std::fopen((outdir + "/trajectories.csv").c_str(), "w");
  std::fprintf(ft, "# Bahnen im r-z-Schnitt.  Jede Zeile ist ein Bahnpunkt; die Bahnen sind\n"
                   "# durch Neustarten der Integration mit wachsender Schrittzahl abgetastet,\n"
                   "# damit die Datei nicht die ganze Historie tragen muss.\n");
  std::fprintf(ft, "case,particle,step,r_m,z_m\n");
  std::FILE* fi = std::fopen((outdir + "/impacts.csv").c_str(), "w");
  std::fprintf(fi, "# Auftrefforte, Schicksale, Energien und der getragene Strom.\n"
                   "# Der gestartete Strom ist eine EINGABE.\n");
  std::fprintf(fi, "case,particle,r_start_m,z_start_m,r_end_m,z_end_m,fate,current_A,"
                   "time_s,energy_gain_J,energy_gain_eV,dV_V,q_dV_J,speed_m_per_s\n");
  std::FILE* fb = std::fopen((outdir + "/balance.csv").c_str(), "w");
  std::fprintf(fb, "# Strom- und Teilchenbilanz.  extrahiert + abgefangen + noch fliegend\n"
                   "# = gestartet, exakt: jedes Teilchen hat genau ein Schicksal.\n");
  std::fprintf(fb, "case,species,launched_A,extracted_A,intercepted_A,flying_A,emitter_A,"
                   "polymer_A,extractor_A,n_launched,n_extracted,n_intercepted,n_flying,"
                   "closure_error,max_energy_error,transmission\n");

  for (const Case& c : cases) {
    const TransportResult r = transport_particles(mesh, *c.phi, eps_r, active, domain, *c.sp,
                                                  launch, dt, steps);
    r.balance.print(stdout);
    r.balance.print(log);
    if (!(r.balance.closure_error < 1.0e-13)) exit_code = 2;

    for (std::size_t k = 0; k < r.particles.size(); ++k) {
      const TracedParticle& q = r.particles[k];
      const Real dV = q.potential_start - q.potential_end;
      std::fprintf(fi, "%s,%zu,%.9e,%.9e,%.9e,%.9e,%s", c.tag, k, q.x_start.r, q.x_start.z,
                   q.x_end.r, q.x_end.z, to_string(q.fate));
      put(fi, q.current);
      put(fi, q.time);
      put(fi, q.energy_gain);
      put(fi, q.energy_gain / constants::e);
      put(fi, dV);
      put(fi, c.sp->charge * dV);
      put(fi, norm(q.v));
      std::fprintf(fi, "\n");
    }
    // Sample the trajectories by re-integrating with a growing step budget.
    for (std::size_t k = 0; k < launch.size(); k += 4) {
      std::vector<TracedParticle> one = {launch[k]};
      const int total = std::max(1, r.particles[k].steps);
      for (int n = 0; n <= 60; ++n) {
        const int budget = 1 + (total - 1) * n / 60;
        const TransportResult t =
            transport_particles(mesh, *c.phi, eps_r, active, domain, *c.sp, one, dt, budget);
        std::fprintf(ft, "%s,%zu,%d,%.9e,%.9e\n", c.tag, k, budget, t.particles[0].x.r,
                     t.particles[0].x.z);
      }
    }
    const Real trans =
        (r.balance.launched > 0.0) ? r.balance.extracted / r.balance.launched : 0.0;
    std::fprintf(fb, "%s,%s", c.tag, c.sp->name);
    put(fb, r.balance.launched);
    put(fb, r.balance.extracted);
    put(fb, r.balance.intercepted);
    put(fb, r.balance.still_flying);
    put(fb, r.balance.hit_emitter);
    put(fb, r.balance.hit_polymer);
    put(fb, r.balance.hit_extractor);
    std::fprintf(fb, ",%lld,%lld,%lld,%lld", static_cast<long long>(r.balance.n_launched),
                 static_cast<long long>(r.balance.n_extracted),
                 static_cast<long long>(r.balance.n_intercepted),
                 static_cast<long long>(r.balance.n_flying));
    put(fb, r.balance.closure_error);
    put(fb, r.balance.max_energy_error);
    put(fb, trans);
    std::fprintf(fb, "\n");
    {
      char buf[256];
      std::snprintf(buf, sizeof buf,
                    "  %-24s Transmission %.4f, Schliessfehler %.2e, groesster "
                    "Energiefehler %.2e",
                    c.tag, trans, r.balance.closure_error, r.balance.max_energy_error);
      say(buf);
    }
  }
  std::fclose(ft);
  std::fclose(fi);
  std::fclose(fb);

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_trajectories (P7)\nphase=P7\nstatus=validated_subset\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "box_radius_m=%.9e\nbox_length_m=%.9e\nvoltage_V=%.9e\n"
                    "aperture_radius_m=%.9e\n", R, Z, V, ap);
    std::fprintf(f, "dt_s=%.9e\nmax_steps=%d\n", dt, steps);
    std::fprintf(f, "prescribed_space_charge_C=%.9e\n", dep.total_particles);
    std::fprintf(f, "source=prescribed launch distribution; P5 is blocked, no physical "
                    "particle source; the launched current is an INPUT\n");
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  say("");
  say(exit_code == 0 ? "Alle deklarierten Pruefungen dieses Laufs bestanden."
                     : "MINDESTENS EINE DEKLARIERTE PRUEFUNG IST FEHLGESCHLAGEN.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_trajectories: %s\n", e.what());
  return 2;
}
