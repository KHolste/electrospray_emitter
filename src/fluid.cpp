#include "es/fluid.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {
namespace {

using constants::eV;
using constants::N_A;

constexpr Real g_per_mol = 1e-3 / N_A;  // g/mol -> kg per particle

std::string canon(const std::string& s) {
  std::string o;
  for (char c : s) {
    if (c == '-' || c == '_' || c == ' ' || c == '.' || c == '+') continue;
    o.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  }
  return o;
}

struct Entry {
  const char* key;
  Fluid f;
};

std::vector<Entry> build_table() {
  std::vector<Entry> t;

  {  // 1-ethyl-3-methylimidazolium tetrafluoroborate -- the workhorse of PIR
     // electrospray thrusters.  High conductivity, low viscosity.
    Fluid f;
    f.name = "EMI-BF4";
    f.rho = 1279.0;
    f.gamma = 0.0452;
    f.K = 1.36;
    f.mu = 0.0371;
    f.eps_r = 12.8;
    f.M_cation = 111.17 * g_per_mol;  // EMI+
    f.M_anion = 86.81 * g_per_mol;    // BF4-
    f.dG_solvation = 1.09 * eV;
    t.push_back({"emibf4", f});
  }
  {  // EMI bis(trifluoromethylsulfonyl)imide, also written EMI-TFSI / EMIM-Tf2N.
     // Lower conductivity and surface tension, much heavier anion -- markedly
     // lower Isp in the negative-emission polarity.
    Fluid f;
    f.name = "EMI-Im";
    f.rho = 1518.0;
    f.gamma = 0.0348;
    f.K = 0.88;
    f.mu = 0.0340;
    f.eps_r = 12.0;
    f.M_cation = 111.17 * g_per_mol;  // EMI+
    f.M_anion = 280.15 * g_per_mol;   // Tf2N-
    f.dG_solvation = 1.05 * eV;
    t.push_back({"emiim", f});
    t.push_back({"emitfsi", f});
    t.push_back({"emimtf2n", f});
  }
  {  // 1-butyl-3-methylimidazolium tetrafluoroborate.  Roughly 4x lower
     // conductivity and 3x higher viscosity than EMI-BF4 -- useful for probing
     // how the operating point moves with K.
    Fluid f;
    f.name = "BMI-BF4";
    f.rho = 1201.0;
    f.gamma = 0.0444;
    f.K = 0.35;
    f.mu = 0.104;
    f.eps_r = 11.7;
    f.M_cation = 139.22 * g_per_mol;  // BMI+
    f.M_anion = 86.81 * g_per_mol;
    f.dG_solvation = 1.10 * eV;
    t.push_back({"bmibf4", f});
  }
  {  // Classic droplet-mode reference liquid: formamide doped with NaI.  Very
     // high permittivity, so it sits in the regime where the de la Mora
     // correlations were originally established.  Conductivity is set by the
     // salt loading -- override it.
    Fluid f;
    f.name = "formamide+NaI";
    f.rho = 1130.0;
    f.gamma = 0.058;
    f.K = 0.1;
    f.mu = 0.0033;
    f.eps_r = 111.0;
    f.M_cation = 22.99 * g_per_mol;   // Na+
    f.M_anion = 126.90 * g_per_mol;   // I-
    f.dG_solvation = 2.0 * eV;        // effectively no field evaporation
    f.mean_solvation_n = 0.0;
    t.push_back({"formamidenai", f});
    t.push_back({"formamide", f});
  }
  {  // Glycerol + NaI: the other classic, three orders of magnitude more
     // viscous.  Handy for checking that nothing in the model secretly assumes
     // an inviscid jet.
    Fluid f;
    f.name = "glycerol+NaI";
    f.rho = 1260.0;
    f.gamma = 0.063;
    f.K = 0.01;
    f.mu = 1.0;
    f.eps_r = 42.0;
    f.M_cation = 22.99 * g_per_mol;
    f.M_anion = 126.90 * g_per_mol;
    f.dG_solvation = 2.0 * eV;
    f.mean_solvation_n = 0.0;
    t.push_back({"glycerolnai", f});
    t.push_back({"glycerol", f});
  }
  return t;
}

const std::vector<Entry>& table() {
  static const std::vector<Entry> t = build_table();
  return t;
}

}  // namespace

Fluid Fluid::at_temperature(Real T) const {
  Fluid f = *this;
  f.T_ref = T;
  const Real dT = T - T_ref;
  f.gamma = std::max(1e-6, gamma + dgamma_dT * dT);
  f.rho = std::max(1.0, rho + drho_dT * dT);
  if (T > vft_T0 + 1.0 && T_ref > vft_T0 + 1.0) {
    f.K = K * std::exp(-vft_B_K / (T - vft_T0) + vft_B_K / (T_ref - vft_T0));
    f.mu = mu * std::exp(vft_B_mu / (T - vft_T0) - vft_B_mu / (T_ref - vft_T0));
  }
  return f;
}

Real Fluid::charge_relaxation_time() const { return constants::eps0 * eps_r / K; }

Real Fluid::ehd_length() const {
  // r* = (Q_min eps0 eps_r / K)^(1/3) evaluated at Q_min = gamma eps0 eps_r/(rho K),
  // i.e. the cone-jet radius scale at the stability floor:
  //     r* = ( gamma eps0^2 eps_r^2 / (rho K^2) )^(1/3)
  // For EMI-BF4 this is a few nanometres -- if you get micrometres, a factor of
  // eps0 has gone missing somewhere.
  return std::cbrt(gamma * constants::eps0 * constants::eps0 * eps_r * eps_r / (rho * K * K));
}

Real Fluid::q_min() const { return gamma * constants::eps0 * eps_r / (rho * K); }

Real Fluid::qm_bare_cation() const { return constants::e / M_cation; }

Real Fluid::qm_cluster() const {
  // A cluster [EMI]_{n+1}[BF4]_n carries one net charge and n ion pairs on top
  // of the bare cation.
  const Real n = std::max(0.0, mean_solvation_n);
  return constants::e / (M_cation + n * (M_cation + M_anion));
}

void Fluid::print(std::FILE* f) const {
  std::fprintf(f, "fluid                 : %s  (properties at T = %.2f K)\n", name.c_str(), T_ref);
  std::fprintf(f, "  density             : %10.1f kg/m^3\n", rho);
  std::fprintf(f, "  surface tension     : %10.4f N/m\n", gamma);
  std::fprintf(f, "  conductivity        : %10.4f S/m\n", K);
  std::fprintf(f, "  viscosity           : %10.4g Pa s\n", mu);
  std::fprintf(f, "  rel. permittivity   : %10.2f\n", eps_r);
  std::fprintf(f, "  dG (ion evaporation): %10.3f eV\n", dG_solvation / constants::eV);
  std::fprintf(f, "  charge relaxation   : %10.3g s\n", charge_relaxation_time());
  std::fprintf(f, "  EHD length r*       : %10.3g m\n", ehd_length());
  std::fprintf(f, "  Q_min (cone-jet)    : %10.3g m^3/s  (= %.3g nL/s)\n", q_min(), q_min() * 1e12);
  std::fprintf(f, "  q/m bare cation     : %10.3g C/kg\n", qm_bare_cation());
  std::fprintf(f, "  q/m mean cluster    : %10.3g C/kg  (n = %.2f)\n", qm_cluster(),
               mean_solvation_n);
}

Fluid fluid_by_name(const std::string& name) {
  const std::string c = canon(name);
  for (const Entry& e : table())
    if (c == e.key) return e.f;
  std::string msg = "unknown fluid '" + name + "'; known:";
  for (const std::string& n : fluid_names()) msg += " " + n;
  throw std::runtime_error(msg);
}

std::vector<std::string> fluid_names() {
  std::vector<std::string> v;
  for (const Entry& e : table())
    if (std::find(v.begin(), v.end(), e.f.name) == v.end()) v.push_back(e.f.name);
  return v;
}

}  // namespace es
