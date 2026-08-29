// es_space_charge -- P6: Poisson with free space charge on the axisymmetric mesh.
//
//   es_space_charge <ausgabeverzeichnis> [key=value ...]
//
// WHAT THIS RUN IS.  The Poisson problem with a PRESCRIBED charge distribution
// -- called a TEST SOURCE everywhere, because P5 is blocked and there is no
// physical source of particles.  It writes: the volume charge, the potential
// change against rho = 0, the field change, the mesh convergence against a
// manufactured solution, and the charge conservation of the deposition.
//
// WHAT IT IS NOT.  No emission, no self-consistent PIC loop, no particle
// motion (P7), no ring singularity.  The charge is deposited with the element
// shape functions and enters the FEM load vector; the solution is piecewise
// bilinear and has no singularity to regularise.
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
#include "es/space_charge.hpp"

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

struct Box {
  QuadMesh mesh;
  std::vector<Real> eps_r;
  std::vector<char> active, fixed;
  std::vector<Real> fixed_value;
};

Box grounded_box(Real R, Real L, Index nr, Index nz, Real V) {
  Box b;
  b.mesh = cylinder_mesh_symmetric(R, L, nr, nz);
  b.eps_r.assign(static_cast<std::size_t>(b.mesh.n_cells()), 1.0);
  b.active.assign(static_cast<std::size_t>(b.mesh.n_cells()), 1);
  b.fixed.assign(static_cast<std::size_t>(b.mesh.n_nodes()), 0);
  b.fixed_value.assign(static_cast<std::size_t>(b.mesh.n_nodes()), 0.0);
  for (Index j = 0; j < b.mesh.nz; ++j)
    b.fixed[static_cast<std::size_t>(b.mesh.node(b.mesh.nr - 1, j))] = 1;
  for (Index i = 0; i < b.mesh.nr; ++i) {
    b.fixed[static_cast<std::size_t>(b.mesh.node(i, 0))] = 1;
    const Index c = b.mesh.node(i, b.mesh.nz - 1);
    b.fixed[static_cast<std::size_t>(c)] = 1;
    b.fixed_value[static_cast<std::size_t>(c)] = V;
  }
  return b;
}
}  // namespace

int main(int argc, char** argv) try {
  const std::vector<std::string> pos = Config::positional_args(argc, argv);
  if (pos.empty()) {
    std::printf("es_space_charge -- P6: Poisson mit freier Raumladung\n\n"
                "  es_space_charge <ausgabeverzeichnis> [key=value ...]\n");
    return 1;
  }
  Config cfg;
  cfg.apply_cli(argc, argv);
  const std::string outdir = pos.back();
  std::filesystem::create_directories(outdir);

  const Real R = cfg.num("box.radius", 2.0e-5);
  const Real L = cfg.num("box.half_length", 4.0e-5);
  const Real V = cfg.num("box.voltage", 250.0);
  const Index nr = cfg.integer("box.nr", 61);
  const Index nz = cfg.integer("box.nz", 121);
  int exit_code = 0;

  std::FILE* log = std::fopen((outdir + "/run.log").c_str(), "w");
  auto say = [&](const std::string& s) {
    std::printf("%s\n", s.c_str());
    std::fprintf(log, "%s\n", s.c_str());
  };
  say("P6 -- Poisson mit freier Raumladung auf dem achsensymmetrischen Volumennetz");
  say("");
  say("DIE LADUNGSVERTEILUNG IST EINE VORGESCHRIEBENE TESTQUELLE.  P5 ist");
  say("blockiert, es gibt keine physikalische Teilchenquelle, und aus keiner Zahl");
  say("dieses Laufs darf ein emittierter Strom gelesen werden.");
  say("");
  say("Kein Ringmodell: die Ladung wird mit den Formfunktionen auf die Knoten");
  say("deponiert und geht in den FEM-Lastvektor.  Die Loesung ist stueckweise");
  say("bilinear und hat gar keine Singularitaet, die regularisiert werden muesste.");

  // --- 1. the manufactured solution and its convergence ---------------------
  {
    std::FILE* f = std::fopen((outdir + "/manufactured.csv").c_str(), "w");
    std::fprintf(f, "# Hergestellte Poisson-Loesung:\n"
                    "#   phi = phi0 (R^2-r^2)(L^2-z^2)/(R^2 L^2),\n"
                    "#   rho = -eps0 lap phi = (2 eps0 phi0/(R^2 L^2)) [2(L^2-z^2)+(R^2-r^2)].\n"
                    "# Der Feldfehler ist NACH RADIUS getrennt: die volumengewichtete\n"
                    "# Knotenrekonstruktion traegt in Achsensymmetrie einen Faktor 2 pi r und\n"
                    "# ist deshalb achsennah nur erster Ordnung.  Ein einziger Wert fuer das\n"
                    "# ganze Gebiet wuerde das verstecken.\n");
    std::fprintf(f, "nr,nz,n_nodes,h_over_R,phi_error,E_error_far,E_error_near\n");
    const Real phi0 = 250.0;
    for (Index n : {21, 41, 81, 161}) {
      Box b = grounded_box(R, L, n, 2 * n - 1, 0.0);
      std::vector<Real> rho(static_cast<std::size_t>(b.mesh.n_nodes()), 0.0);
      for (Index j = 0; j < b.mesh.nz; ++j)
        for (Index i = 0; i < b.mesh.nr; ++i)
          rho[static_cast<std::size_t>(b.mesh.node(i, j))] =
              manufactured_charge_density(b.mesh.at(i, j), R, L, phi0);
      const SpaceChargeSolution s =
          solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed, b.fixed_value, {}, rho);
      Real wp = 0.0, wf = 0.0, wn = 0.0;
      for (Index j = 1; j + 1 < b.mesh.nz; ++j)
        for (Index i = 0; i + 1 < b.mesh.nr; ++i) {
          const Vec2 x = b.mesh.at(i, j);
          wp = std::max(wp, std::abs(s.phi[static_cast<std::size_t>(b.mesh.node(i, j))] -
                                     manufactured_potential(x, R, L, phi0)));
          const Vec2 Ee{2.0 * phi0 * x.r * (L * L - x.z * x.z) / (R * R * L * L),
                        2.0 * phi0 * (R * R - x.r * x.r) * x.z / (R * R * L * L)};
          const Real d =
              norm(field_recovered_at_node(b.mesh, s.phi, b.eps_r, b.active, i, j, 1.0) - Ee);
          if (x.r < 0.25 * R)
            wn = std::max(wn, d);
          else
            wf = std::max(wf, d);
        }
      std::fprintf(f, "%lld,%lld,%lld", static_cast<long long>(n),
                   static_cast<long long>(2 * n - 1),
                   static_cast<long long>(b.mesh.n_nodes()));
      put(f, 1.0 / static_cast<Real>(n - 1));
      put(f, wp / phi0);
      put(f, wf / (phi0 / R));
      put(f, wn / (phi0 / R));
      std::fprintf(f, "\n");
    }
    std::fclose(f);
    say("  manufactured.csv geschrieben");
  }

  // --- 2. the test source, and what it does to the field --------------------
  Box b = grounded_box(R, L, nr, nz, V);
  std::vector<Macroparticle> parts;
  {
    // A PRESCRIBED beam-like column of macroparticles on and near the axis.
    const int np = cfg.integer("source.n_particles", 400);
    const Real q = cfg.num("source.total_charge", 2.0e-14) / static_cast<Real>(np);
    const Real rb = cfg.num("source.radius", 0.25) * R;
    // A prescribed column: uniform in z over the middle half, and uniform over
    // the disc of radius rb in the cross section.  The radial coordinate uses a
    // van der Corput sequence so that the sampling is even and does not print
    // as diagonal streaks -- which is cosmetic, but a figure that looks like a
    // pattern invites being read as one.
    auto van_der_corput = [](int n) {
      Real x = 0.0, base = 0.5;
      while (n > 0) {
        x += base * (n % 2);
        n /= 2;
        base *= 0.5;
      }
      return x;
    };
    for (int k = 0; k < np; ++k) {
      const Real t = (static_cast<Real>(k) + 0.5) / static_cast<Real>(np);
      const Real rr = rb * std::sqrt(van_der_corput(k + 1));
      const Real zz = -0.5 * L + L * t;
      parts.push_back({{rr, zz}, q});
    }
  }
  const DepositionResult d = deposit(b.mesh, parts);
  say("");
  say("  Testquelle: " + std::to_string(parts.size()) + " Makropartikel, Summe " +
      std::to_string(d.total_particles) + " C, deponiert " +
      std::to_string(d.total_deposited) + " C, Erhaltungsfehler " +
      std::to_string(d.conservation_error));
  if (!(d.conservation_error < 1.0e-13) || d.n_outside != 0) exit_code = 2;

  const SpaceChargeSolution s =
      solve_with_space_charge(b.mesh, b.eps_r, b.active, b.fixed, b.fixed_value, d.node_charge,
                              {});
  say("  Potentialverschiebung gegen rho = 0: " + std::to_string(s.max_potential_shift) +
      " V; Feldverschiebung " + std::to_string(s.max_field_shift) + " V/m");

  {
    std::FILE* f = std::fopen((outdir + "/fields.csv").c_str(), "w");
    std::fprintf(f, "# Volumenladung, Potential mit und ohne Ladung, und die Differenz.\n"
                    "# Die Ladungsdichte ist die deponierte Knotenladung geteilt durch das\n"
                    "# zugehoerige Knotenvolumen -- eine ABGELEITETE Groesse zur Anzeige,\n"
                    "# nicht die Groesse, mit der gerechnet wurde (das ist die Knotenladung).\n");
    std::fprintf(f, "r_m,z_m,node_charge_C,phi_V,phi_no_charge_V,dphi_V,Er_V_per_m,"
                    "Ez_V_per_m,dEr_V_per_m,dEz_V_per_m\n");
    for (Index j = 0; j < b.mesh.nz; j += std::max<Index>(1, b.mesh.nz / 240))
      for (Index i = 0; i < b.mesh.nr; i += std::max<Index>(1, b.mesh.nr / 120)) {
        const Vec2 x = b.mesh.at(i, j);
        const std::size_t n = static_cast<std::size_t>(b.mesh.node(i, j));
        const Vec2 Ea = field_recovered_at_node(b.mesh, s.phi, b.eps_r, b.active, i, j, 1.0);
        const Vec2 Eb =
            field_recovered_at_node(b.mesh, s.phi_no_charge, b.eps_r, b.active, i, j, 1.0);
        std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", x.r, x.z,
                     d.node_charge[n], s.phi[n], s.phi_no_charge[n],
                     s.phi[n] - s.phi_no_charge[n], Ea.r, Ea.z, Ea.r - Eb.r, Ea.z - Eb.z);
      }
    std::fclose(f);
    say("  fields.csv geschrieben");
  }

  // --- 3. conservation and the self-field scaling ---------------------------
  {
    std::FILE* f = std::fopen((outdir + "/conservation.csv").c_str(), "w");
    std::fprintf(f, "# Ladungserhaltung der Deposition ueber die Netzstufen, und das\n"
                    "# Spitzenpotential eines EINZELNEN Makropartikels.  Letzteres waechst\n"
                    "# unter Verfeinerung: das Selbstfeld ist eine NETZGROESSE.  Eine\n"
                    "# PIC-Schleife darf ein Teilchen sein eigenes deponiertes Feld nicht\n"
                    "# ungefiltert fuehlen lassen.\n");
    std::fprintf(f, "nr,nz,n_nodes,deposited_C,particles_C,conservation_error,"
                    "partition_of_unity_error,single_particle_peak_phi_V\n");
    const std::vector<Macroparticle> one = {{{0.35 * R, 0.0}, 1.0e-15}};
    for (Index n : {21, 41, 81, 161}) {
      Box c = grounded_box(R, L, n, 2 * n - 1, 0.0);
      const DepositionResult dd = deposit(c.mesh, parts);
      const DepositionResult d1 = deposit(c.mesh, one);
      const SpaceChargeSolution s1 = solve_with_space_charge(
          c.mesh, c.eps_r, c.active, c.fixed, c.fixed_value, d1.node_charge, {});
      Real peak = 0.0;
      for (Real v : s1.phi) peak = std::max(peak, std::abs(v));
      std::fprintf(f, "%lld,%lld,%lld", static_cast<long long>(n),
                   static_cast<long long>(2 * n - 1),
                   static_cast<long long>(c.mesh.n_nodes()));
      put(f, dd.total_deposited);
      put(f, dd.total_particles);
      put(f, dd.conservation_error);
      put(f, dd.partition_of_unity_error);
      put(f, peak);
      std::fprintf(f, "\n");
    }
    std::fclose(f);

    // The approach study: potential and field along a line through a particle.
    std::FILE* g = std::fopen((outdir + "/approach.csv").c_str(), "w");
    std::fprintf(g, "# Annaeherung an ein einzelnes Makropartikel.  Ein Ringmodell wuerde\n"
                    "# hier logarithmisch divergieren; die FEM-Loesung eines Knotenlast-\n"
                    "# vektors ist stueckweise bilinear und bleibt beschraenkt.\n");
    std::fprintf(g, "d_over_R,phi_V,E_magnitude_V_per_m\n");
    Box c = grounded_box(R, L, 81, 161, 0.0);
    const DepositionResult d1 = deposit(c.mesh, one);
    const SpaceChargeSolution s1 = solve_with_space_charge(
        c.mesh, c.eps_r, c.active, c.fixed, c.fixed_value, d1.node_charge, {});
    for (int k = 0; k < 14; ++k) {
      const Real dr = 0.3 * R * std::pow(0.5, k);
      const Vec2 x{0.35 * R + dr, 0.0};
      std::fprintf(g, "%.9e,%.9e,%.9e\n", dr / R, potential_at(c.mesh, s1.phi, x),
                   norm(interpolated_field(c.mesh, s1.phi, c.eps_r, c.active, x)));
    }
    std::fclose(g);
    say("  conservation.csv und approach.csv geschrieben");
  }

  {
    std::FILE* f = std::fopen((outdir + "/meta.txt").c_str(), "w");
    std::fprintf(f, "app=es_space_charge (P6)\nphase=P6\nstatus=validated_subset\n");
    std::fprintf(f, "commit=%s\n", cfg.str("meta.commit", "unbekannt").c_str());
    std::fprintf(f, "box_radius_m=%.9e\nbox_half_length_m=%.9e\nbox_voltage_V=%.9e\n", R, L, V);
    std::fprintf(f, "n_particles=%zu\ntotal_charge_C=%.9e\n", parts.size(), d.total_particles);
    std::fprintf(f, "conservation_error=%.9e\n", d.conservation_error);
    std::fprintf(f, "max_potential_shift_V=%.9e\nmax_field_shift_V_per_m=%.9e\n",
                 s.max_potential_shift, s.max_field_shift);
    std::fprintf(f, "source=prescribed test source; P5 is blocked, there is no physical "
                    "particle source\n");
    std::fprintf(f, "exit_code=%d\n", exit_code);
    std::fclose(f);
  }
  say("");
  say(exit_code == 0 ? "Alle deklarierten Pruefungen dieses Laufs bestanden."
                     : "MINDESTENS EINE DEKLARIERTE PRUEFUNG IST FEHLGESCHLAGEN.");
  std::fclose(log);
  return exit_code;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_space_charge: %s\n", e.what());
  return 2;
}
