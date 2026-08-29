#include "es/liquid.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {
namespace {

std::string canon(const std::string& s) {
  std::string o;
  for (char c : s) {
    if (c == '-' || c == '_' || c == ' ' || c == '.' || c == '+') continue;
    o.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  }
  return o;
}

}  // namespace

const char* to_string(LiquidDataStatus s) {
  switch (s) {
    case LiquidDataStatus::Verified: return "verified";
    case LiquidDataStatus::Provisional: return "provisional";
    case LiquidDataStatus::Illustrative: return "illustrative";
    case LiquidDataStatus::Unknown: return "unknown";
  }
  return "?";
}

const char* explain(LiquidDataStatus s) {
  switch (s) {
    case LiquidDataStatus::Verified:
      return "Stoff, Temperatur und Primaerquelle sind eindeutig belegt.";
    case LiquidDataStatus::Provisional:
      return "Benannte, aber nicht im Volltext gepruefte Quelle. Traegt eine "
             "Empfindlichkeitsstudie, keine Absolutaussage.";
    case LiquidDataStatus::Illustrative:
      return "Beispielwert ohne Primaerquelle. Traegt ausschliesslich eine "
             "dimensionslose Demonstration.";
    case LiquidDataStatus::Unknown:
      return "Registriert, aber ohne Wert. Die Verwendung ist ein Fehler.";
  }
  return "?";
}

LiquidDataStatus liquid_status_from_string(const std::string& s) {
  const std::string c = canon(s);
  if (c == "verified") return LiquidDataStatus::Verified;
  if (c == "provisional") return LiquidDataStatus::Provisional;
  if (c == "illustrative") return LiquidDataStatus::Illustrative;
  if (c == "unknown") return LiquidDataStatus::Unknown;
  throw std::runtime_error(
      "liquid.status muss verified, provisional, illustrative oder unknown sein; "
      "angegeben war: '" + s + "'");
}

// ---------------------------------------------------------------------------

Real LiquidProperties::bond_number(Real a) const {
  if (!(gamma > 0.0) || !(a > 0.0)) return 0.0;
  return rho * constants::g0 * a * a / gamma;
}

Real LiquidProperties::capillary_pressure_scale(Real a) const {
  if (!(a > 0.0)) return 0.0;
  return gamma / a;
}

std::string LiquidProperties::why_unusable() const {
  if (status == LiquidDataStatus::Unknown)
    return "Der Stoffdatensatz '" + substance +
           "' ist registriert, traegt aber keinen Status und damit keinen Wert.";
  if (!std::isfinite(gamma) || gamma <= 0.0)
    return "Die Oberflaechenspannung muss endlich und positiv sein, ist aber " +
           std::to_string(gamma) + " N/m.";
  if (!std::isfinite(rho) || rho <= 0.0)
    return "Die Dichte muss endlich und positiv sein, ist aber " + std::to_string(rho) +
           " kg/m^3.  Sie geht in P3a nicht in das Gleichgewicht ein, wohl aber in die "
           "Bond-Zahl, mit der die Vernachlaessigung der Schwerkraft begruendet wird.";
  if (!std::isfinite(T) || T <= 0.0)
    return "Die Bezugstemperatur muss endlich und positiv sein, ist aber " +
           std::to_string(T) + " K.";
  if (status != LiquidDataStatus::Illustrative && source.empty())
    return "Der Datensatz '" + substance + "' hat den Status " + to_string(status) +
           ", nennt aber keine Fundstelle.  Ein Stoffwert ohne Herkunft ist keine "
           "Eingabe, sondern eine Behauptung.";
  return {};
}

void LiquidProperties::validate_or_throw() const {
  const std::string why = why_unusable();
  if (!why.empty()) throw std::runtime_error(why);
}

void LiquidProperties::print(std::FILE* out) const {
  std::fprintf(out, "Fluessigkeit  : %s\n", substance.c_str());
  std::fprintf(out, "  Status      : %s -- %s\n", to_string(status), explain(status));
  std::fprintf(out, "  Temperatur  : %.4g K\n", T);
  std::fprintf(out, "  gamma       : %.6g N/m   (geht in das Gleichgewicht ein)\n", gamma);
  std::fprintf(out, "  rho         : %.6g kg/m^3 (nur Bond-Zahl, nicht im Gleichgewicht)\n", rho);
  std::fprintf(out, "  Fundstelle  : %s\n", source.empty() ? "(keine)" : source.c_str());
  if (!caveat.empty()) std::fprintf(out, "  Vorbehalt   : %s\n", caveat.c_str());
  std::fprintf(out, "  NUR VORGEMERKT, in P3a nicht verwendet: mu = %.6g Pa s, K = %.6g S/m, "
                    "eps_r = %.6g\n",
               documented_only.mu, documented_only.K, documented_only.eps_r);
}

// ---------------------------------------------------------------------------

LiquidProperties emibf4_illustrative() {
  LiquidProperties L;
  L.substance = "EMI-BF4 (1-Ethyl-3-methylimidazolium-tetrafluoroborat)";
  L.T = 298.15;
  // These two numbers are carried over UNCHANGED from the unsourced table in
  // src/fluid.cpp.  They are not measurements and no primary source for them is
  // in hand, which is exactly what the status says.
  L.gamma = 0.0452;
  L.rho = 1279.0;
  L.status = LiquidDataStatus::Illustrative;
  L.source =
      "Stoffidentitaet belegt: KunzeFynn-2024-12-10.pdf, Abschnitt 2.3.2, gedruckte Seite 28 "
      "(PDF-Seite 36), sowie Tabelle 'List of Publications', gedruckte Seite 30 "
      "(PDF-Seite 38), Publikationen I-IV.  Die ZAHLENWERTE stammen NICHT daraus: die "
      "Dissertation enthaelt keine Stoffwerttabelle und keinen numerischen Wert fuer "
      "Oberflaechenspannung oder Dichte.  gamma und rho sind unveraendert aus der "
      "quellenlosen Tabelle in src/fluid.cpp uebernommen.";
  L.caveat =
      "Nur fuer ein gekennzeichnetes Stoffbeispiel verwendbar.  Kein absoluter Druck, keine "
      "absolute Apexhoehe und keine Betriebsspannung darf auf diesen Datensatz gestuetzt "
      "werden.  Die Pruefung des Loesers laeuft unabhaengig davon dimensionslos.";
  L.documented_only.mu = 0.0371;
  L.documented_only.K = 1.36;
  L.documented_only.eps_r = 12.8;
  return L;
}

LiquidProperties unit_liquid() {
  LiquidProperties L;
  L.substance = "Einheitsfluessigkeit (kein Stoff)";
  L.T = 298.15;
  L.gamma = 1.0;
  L.rho = 1000.0;
  L.status = LiquidDataStatus::Illustrative;
  L.source = "Definition, kein Stoff: gamma = 1 N/m exakt, fuer dimensionslose Pruefungen.";
  L.caveat = "Beschreibt keine reale Fluessigkeit.";
  return L;
}

LiquidProperties liquid_data_by_name(const std::string& name) {
  const std::string c = canon(name);
  if (c == "emibf4") return emibf4_illustrative();
  if (c == "unit" || c == "einheit") return unit_liquid();
  std::string known;
  for (const std::string& n : liquid_data_names()) known += (known.empty() ? "" : ", ") + n;
  throw std::runtime_error("Unbekannter Stoffdatensatz '" + name + "'.  Bekannt sind: " +
                           known + ".  Ein anderer Stoff wird ueber liquid.* in der "
                           "Konfiguration vollstaendig angegeben, mit Status und Fundstelle.");
}

std::vector<std::string> liquid_data_names() { return {"emi-bf4", "unit"}; }

}  // namespace es
