#include "es/cone_jet_contract.hpp"

#include <cmath>
#include <limits>

#include "es/constants.hpp"
#include "es/status.hpp"

namespace es {

using constants::eps0;
using constants::pi;

namespace {
constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();
}

const char* to_string(ConeJetStatus s) {
  switch (s) {
    case ConeJetStatus::Blocked: return "Blocked";
    case ConeJetStatus::MissingMaterialData: return "MissingMaterialData";
    case ConeJetStatus::MissingOperatingPoint: return "MissingOperatingPoint";
    case ConeJetStatus::DiagnosisAvailable: return "DiagnosisAvailable";
  }
  return "?";
}

const char* explain(ConeJetStatus s) {
  switch (s) {
    case ConeJetStatus::Blocked:
      return "Der Cone-Jet-Modus ist blockiert.  Er braucht eine Zweiphasenstroemung mit "
             "freier Oberflaeche und endlicher Leitfaehigkeit; P3 und P4 haben festgestellt, "
             "dass es beides hier nicht gibt.  Die empirische Skalierung wird ebenfalls "
             "nicht uebernommen, weil Originalgleichung, Erratum, Vorfaktor und "
             "Gueltigkeitsbereich in diesem Lauf an keiner Quelle geprueft werden konnten.";
    case ConeJetStatus::MissingMaterialData:
      return "Eine Stoffgroesse, die die Diagnose braucht, ist MissingMaterialData.  Es wird "
             "kein Ersatzwert gesetzt; die betroffenen Kennzahlen bleiben nan.";
    case ConeJetStatus::MissingOperatingPoint:
      return "Volumenstrom oder Kanalradius fehlen.  Beide sind Eingaben und werden hier "
             "nicht bestimmt.";
    case ConeJetStatus::DiagnosisAvailable:
      return "Die dimensionslose Diagnose ist auswertbar.  Sie ist eine DIAGNOSE: sie sagt, "
             "welche Physik dominieren wuerde, und sagt KEINEN Betriebsmodus voraus.";
  }
  return "?";
}

// ---------------------------------------------------------------------------

std::vector<ConeJetRequirement> cone_jet_requirements() {
  return {
      {"Zweiphasenstroemung mit freier Oberflaeche",
       "Ein Cone-Jet HAT einen Jet: ein schlankes geladenes Filament, das in Tropfen "
       "zerfaellt.  Ohne einen Loeser fuer die Innenstroemung mit beweglicher Oberflaeche "
       "gibt es weder Jet noch Zerfall.",
       false, "fehlt; P4 spezifiziert den Vertrag (docs/15)"},
      {"endliche Leitfaehigkeit mit Oberflaechenladungstransport",
       "Der Strom eines Cone-Jets wird vom Leitungsstrom in den Jet hinein und von der "
       "Konvektion der Flaechenladung getragen.  Im Perfect-Conductor-Grenzfall existiert "
       "die dafuer noetige tangentiale Traktion q_s E_t nicht einmal.",
       false, "fehlt; P3 spezifiziert den Vertrag (docs/14, 14.5)"},
      {"belegte relative Permittivitaet",
       "Sie geht in die Ladungsrelaxationszeit und in jede elektrohydrodynamische Laenge "
       "ein.",
       false, "MissingMaterialData (P2, docs/13)"},
      {"belegte Leitfaehigkeit", "dieselbe Rolle wie die Permittivitaet.", true,
       "P2, gewaehlte Quelle mit Methode, Reinheit und Wassergehalt"},
      {"belegte Oberflaechenspannung, Dichte, Viskositaet",
       "Sie setzen alle kapillaren und viskosen Zeitskalen.", true, "P2"},
      {"Zerfallsmodell des Jets",
       "Der Tropfendurchmesser folgt aus einer Instabilitaet des Jets, nicht aus dem Jet "
       "selbst.  Ohne Zeitintegration einer freien Oberflaeche gibt es das nicht.",
       false, "fehlt; P4 ist infrastructure_only"},
      {"gepruefte empirische Skalierung",
       "Zulaessig nur nach Pruefung von Originalgleichung, Erratum, Vorfaktor und "
       "Gueltigkeitsbereich an der Quelle.  In diesem Lauf war weder Ganan-Calvo (1997) "
       "noch das Erratum (2000) im Volltext erreichbar.",
       false, "nicht geprueft"},
  };
}

// ---------------------------------------------------------------------------

Real charge_relaxation_time_cj(Real eps_r, Real K) {
  if (!(eps_r > 0.0) || !(K > 0.0)) return kNaN;
  return eps0 * eps_r / K;
}

Real capillary_inertial_time(Real rho, Real gamma, Real a) {
  if (!(rho > 0.0) || !(gamma > 0.0) || !(a > 0.0)) return kNaN;
  return std::sqrt(rho * a * a * a / gamma);
}

Real visco_capillary_time(Real mu, Real gamma, Real a) {
  if (!(mu > 0.0) || !(gamma > 0.0) || !(a > 0.0)) return kNaN;
  return mu * a / gamma;
}

Real electrohydrodynamic_length(Real gamma, Real eps_r, Real rho, Real K) {
  if (!(gamma > 0.0) || !(eps_r > 0.0) || !(rho > 0.0) || !(K > 0.0)) return kNaN;
  // From  eps0 eps_r / K = sqrt(rho r^3 / gamma):
  //   (eps0 eps_r / K)^2 = rho r^3 / gamma
  //   r^3 = gamma eps0^2 eps_r^2 / (rho K^2)
  return std::cbrt(gamma * eps0 * eps0 * eps_r * eps_r / (rho * K * K));
}

Real ohnesorge(Real mu, Real rho, Real gamma, Real a) {
  if (!(mu > 0.0) || !(rho > 0.0) || !(gamma > 0.0) || !(a > 0.0)) return kNaN;
  return mu / std::sqrt(rho * gamma * a);
}

Real electric_bond(Real E, Real a, Real gamma) {
  if (!(a > 0.0) || !(gamma > 0.0)) return kNaN;
  return eps0 * E * E * a / (2.0 * gamma);
}

Real feed_reynolds(Real rho, Real Q, Real R, Real mu) {
  if (!(rho > 0.0) || !(R > 0.0) || !(mu > 0.0)) return kNaN;
  return 2.0 * rho * std::abs(Q) / (pi * R * mu);
}

Real feed_capillary(Real mu, Real Q, Real R, Real gamma) {
  if (!(mu > 0.0) || !(R > 0.0) || !(gamma > 0.0)) return kNaN;
  return mu * std::abs(Q) / (pi * R * R) / gamma;
}

// ---------------------------------------------------------------------------

void ConeJetDiagnosis::print(std::FILE* out) const {
  std::fprintf(out, "  Cone-Jet-Diagnose: %s\n", to_string(status));
  std::fprintf(out, "    Stoffwerte: gamma %s, rho %s, mu %s, K %s, eps_r %s\n",
               to_string(gamma_status), to_string(rho_status), to_string(mu_status),
               to_string(K_status), to_string(eps_r_status));
  std::fprintf(out, "    Zeitskalen [s]: tau_e %.4e, kapillar %.4e, viskos %.4e, "
                    "Verweilzeit %.4e\n",
               tau_e, t_capillary, t_viscous, t_residence);
  std::fprintf(out, "    Laengen [m]: a %.4e, r* %.4e\n", a, r_star);
  std::fprintf(out, "    dimensionslos: Oh %.4e, Bo_E %.4e, Re %.4e, Ca %.4e\n", Oh, Bo_E, Re,
               Ca);
  std::fprintf(out, "    %s\n", message.c_str());
}

ConeJetDiagnosis diagnose_cone_jet(const MaterialDataset& d, Real T, Real a, Real Q,
                                   Real R_channel, Real E_surface) {
  ConeJetDiagnosis c;
  c.a = a;
  c.Q = Q;
  c.R_channel = R_channel;
  c.E_surface = E_surface;

  const MaterialValue g = material_value(d, PropertyKind::SurfaceTension, T);
  const MaterialValue r = material_value(d, PropertyKind::Density, T);
  const MaterialValue m = material_value(d, PropertyKind::DynamicViscosity, T);
  const MaterialValue k = material_value(d, PropertyKind::ElectricalConductivity, T);
  const MaterialValue e = material_value(d, PropertyKind::RelativePermittivity, T);
  c.gamma_status = g.status;
  c.rho_status = r.status;
  c.mu_status = m.status;
  c.K_status = k.status;
  c.eps_r_status = e.status;
  c.gamma = g.usable() ? g.value : 0.0;
  c.rho = r.usable() ? r.value : 0.0;
  c.mu = m.usable() ? m.value : 0.0;
  c.K = k.usable() ? k.value : 0.0;
  c.eps_r = e.usable() ? e.value : 0.0;

  c.tau_e = charge_relaxation_time_cj(c.eps_r, c.K);
  c.t_capillary = capillary_inertial_time(c.rho, c.gamma, a);
  c.t_viscous = visco_capillary_time(c.mu, c.gamma, a);
  c.t_residence = (Q != 0.0 && a > 0.0) ? (4.0 / 3.0 * pi * a * a * a) / std::abs(Q) : kNaN;
  c.r_star = electrohydrodynamic_length(c.gamma, c.eps_r, c.rho, c.K);
  c.Oh = ohnesorge(c.mu, c.rho, c.gamma, a);
  c.Bo_E = electric_bond(E_surface, a, c.gamma);
  c.Re = feed_reynolds(c.rho, Q, R_channel, c.mu);
  c.Ca = feed_capillary(c.mu, Q, R_channel, c.gamma);

  if (!e.usable() || !k.usable() || !g.usable() || !r.usable() || !m.usable()) {
    c.status = ConeJetStatus::MissingMaterialData;
    c.message = std::string(explain(c.status)) + "  Fehlend: ";
    if (!g.usable()) c.message += "gamma ";
    if (!r.usable()) c.message += "rho ";
    if (!m.usable()) c.message += "mu ";
    if (!k.usable()) c.message += "K ";
    if (!e.usable()) c.message += "eps_r ";
    c.message += " -- die davon abhaengigen Kennzahlen bleiben nan.";
    return c;
  }
  if (!(Q != 0.0) || !(R_channel > 0.0)) {
    c.status = ConeJetStatus::MissingOperatingPoint;
    c.message = explain(c.status);
    return c;
  }
  c.status = ConeJetStatus::DiagnosisAvailable;
  c.message = explain(c.status);
  return c;
}

void cone_jet_current() {
  throw NotImplementedInThisPhase(
      "Der Cone-Jet-Strom",
      "nach einer Zweiphasenstroemung mit freier Oberflaeche, einem "
      "Oberflaechenladungstransport und einem geprueften Zerfallsmodell; der Vertrag steht "
      "in docs/19_cone_jet.md",
      "Zwei getrennte Gruende, und keiner davon ist Zeitmangel.  ERSTENS ist der Modus "
      "physikalisch nicht gerechnet: ein Cone-Jet hat einen Jet, und dafuer braucht es die "
      "Zweiphasenstroemung mit freier Oberflaeche und den Oberflaechenladungstransport, die "
      "P3 und P4 als fehlend festgestellt haben.  ZWEITENS wird auch die empirische "
      "Skalierung nicht uebernommen: der Auftrag erlaubt das nur nach Pruefung von "
      "Originalgleichung, Erratum, Vorfaktor und Gueltigkeitsbereich an der Quelle, und in "
      "diesem Lauf war weder Ganan-Calvo, Phys. Rev. Lett. 79, 217 (1997), noch das Erratum, "
      "Phys. Rev. Lett. 85, 4193 (2000), im Volltext erreichbar.  Zahlen aus "
      "Suchtreffer-Schnipseln gelten in diesem Projekt nicht als Quelle.  Der alte "
      "ConeJetModel in emission.hpp bleibt unangetastet und wird von hier aus nicht benutzt.");
}

}  // namespace es
