#include "es/material_data.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace es {

namespace {
constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();
}

const char* to_string(PropertyKind k) {
  switch (k) {
    case PropertyKind::SurfaceTension: return "surface_tension";
    case PropertyKind::Density: return "density";
    case PropertyKind::DynamicViscosity: return "dynamic_viscosity";
    case PropertyKind::KinematicViscosity: return "kinematic_viscosity";
    case PropertyKind::ElectricalConductivity: return "electrical_conductivity";
    case PropertyKind::RelativePermittivity: return "relative_permittivity";
  }
  return "?";
}

const char* si_unit(PropertyKind k) {
  switch (k) {
    case PropertyKind::SurfaceTension: return "N/m";
    case PropertyKind::Density: return "kg/m^3";
    case PropertyKind::DynamicViscosity: return "Pa s";
    case PropertyKind::KinematicViscosity: return "m^2/s";
    case PropertyKind::ElectricalConductivity: return "S/m";
    case PropertyKind::RelativePermittivity: return "-";
  }
  return "?";
}

const char* to_string(MaterialDataStatus s) {
  switch (s) {
    case MaterialDataStatus::Measured: return "measured";
    case MaterialDataStatus::ManufacturerSpec: return "manufacturer_spec";
    case MaterialDataStatus::Literature: return "literature";
    case MaterialDataStatus::Derived: return "derived";
    case MaterialDataStatus::Illustrative: return "illustrative";
    case MaterialDataStatus::MissingMaterialData: return "MissingMaterialData";
  }
  return "?";
}

const char* explain(MaterialDataStatus s) {
  switch (s) {
    case MaterialDataStatus::Measured:
      return "Eine Primaerquelle ist an der Stelle gelesen, an der der Wert steht; Stoff, "
             "Temperatur und Messmethode sind eindeutig.";
    case MaterialDataStatus::ManufacturerSpec:
      return "Ein Herstellerdatenblatt liegt vor.  Der Wert steht mit seiner Temperatur, "
             "Messmethode und Wassergehalt aber meist nicht.";
    case MaterialDataStatus::Literature:
      return "Eine benannte Quelle, die an der Wertstelle nicht gelesen wurde.";
    case MaterialDataStatus::Derived:
      return "NICHT gemessen.  Aus anderen ausgewaehlten Groessen desselben Datensatzes ueber "
             "eine genannte Identitaet berechnet, nachdem die dokumentierten Bedingungen "
             "beider Quellen auf Vertraeglichkeit geprueft wurden; die Unsicherheit ist "
             "fortgepflanzt.  Der Wert darf nicht so dargestellt werden, als sei die Groesse "
             "selbst gemessen worden.";
    case MaterialDataStatus::Illustrative:
      return "Ein Beispielwert ohne Primaerquelle.  Er traegt eine dimensionslose "
             "Vorfuehrung und sonst nichts.";
    case MaterialDataStatus::MissingMaterialData:
      return "Der Wert fehlt.  Es gibt hier KEINEN Ersatzwert: eine Rechnung, die ihn "
             "braucht, muss geschlossen fehlschlagen statt eine Vorgabe zu benutzen.";
  }
  return "?";
}

// ---------------------------------------------------------------------------

std::size_t PropertySource::n_non_ambient() const {
  std::size_t n = 0;
  for (std::size_t k = 0; k < n_points; ++k)
    if (!points[k].ambient()) ++n;
  return n;
}

bool PropertySource::has_ambient_points() const {
  for (std::size_t k = 0; k < n_points; ++k)
    if (points[k].ambient()) return true;
  return false;
}

Real PropertySource::T_min() const {
  Real m = std::numeric_limits<Real>::max();
  bool any = false;
  for (std::size_t k = 0; k < n_points; ++k)
    if (points[k].ambient()) { m = std::min(m, points[k].T); any = true; }
  return any ? m : kNaN;
}

Real PropertySource::T_max() const {
  Real m = -std::numeric_limits<Real>::max();
  bool any = false;
  for (std::size_t k = 0; k < n_points; ++k)
    if (points[k].ambient()) { m = std::max(m, points[k].T); any = true; }
  return any ? m : kNaN;
}

bool PropertySource::covers(Real T) const {
  if (!has_ambient_points()) return false;
  const Real lo = T_min(), hi = T_max();
  if (lo == hi) return std::abs(lo - T) <= 1.0;  // a single ambient point
  return T >= lo - 1e-12 && T <= hi + 1e-12;
}

Real PropertySource::value_at(Real T) const {
  // AMBIENT points only.  The others stay in the record, but they measure the
  // substance in a state this project does not operate it in, and quoting one
  // of them for an ambient number is exactly the kind of silent substitution
  // this file exists to prevent.
  if (!has_ambient_points()) return kNaN;
  const Real tlo = T_min(), thi = T_max();
  if (tlo == thi) {
    if (std::abs(tlo - T) > 1.0) return kNaN;
    for (std::size_t k = 0; k < n_points; ++k)
      if (points[k].ambient() && points[k].T == tlo) return points[k].value;
    return kNaN;
  }
  if (T < tlo - 1e-12 || T > thi + 1e-12) return kNaN;
  // The points are stored in the order the source lists them, which is not
  // necessarily sorted; find the bracketing pair explicitly.
  std::size_t lo = 0, hi = 0;
  Real best_lo = -std::numeric_limits<Real>::max(), best_hi = std::numeric_limits<Real>::max();
  bool have_lo = false, have_hi = false;
  for (std::size_t k = 0; k < n_points; ++k) {
    if (!points[k].ambient()) continue;
    const Real t = points[k].T;
    if (t <= T && t > best_lo) { best_lo = t; lo = k; have_lo = true; }
    if (t >= T && t < best_hi) { best_hi = t; hi = k; have_hi = true; }
  }
  if (!have_lo || !have_hi) return kNaN;
  if (best_hi == best_lo) return points[lo].value;
  const Real w = (T - best_lo) / (best_hi - best_lo);
  return (1.0 - w) * points[lo].value + w * points[hi].value;
}

Real PropertySource::uncertainty_at(Real T) const {
  // Same bracketing as value_at(), but a missing figure poisons the result: a
  // point with uncertainty == 0 did not state one, and interpolating between a
  // stated and an unstated figure would manufacture an error bar.
  if (!has_ambient_points()) return kNaN;
  const Real tlo = T_min(), thi = T_max();
  if (T < tlo - 1e-12 || T > thi + 1e-12) return kNaN;
  if (tlo == thi) {
    if (std::abs(tlo - T) > 1.0) return kNaN;
    for (std::size_t k = 0; k < n_points; ++k)
      if (points[k].ambient() && points[k].T == tlo)
        return points[k].has_uncertainty() ? points[k].uncertainty : kNaN;
    return kNaN;
  }
  std::size_t lo = 0, hi = 0;
  Real best_lo = -std::numeric_limits<Real>::max(), best_hi = std::numeric_limits<Real>::max();
  bool have_lo = false, have_hi = false;
  for (std::size_t k = 0; k < n_points; ++k) {
    if (!points[k].ambient()) continue;
    const Real t = points[k].T;
    if (t <= T && t > best_lo) { best_lo = t; lo = k; have_lo = true; }
    if (t >= T && t < best_hi) { best_hi = t; hi = k; have_hi = true; }
  }
  if (!have_lo || !have_hi) return kNaN;
  if (!points[lo].has_uncertainty() || !points[hi].has_uncertainty()) return kNaN;
  if (best_hi == best_lo) return points[lo].uncertainty;
  const Real w = (T - best_lo) / (best_hi - best_lo);
  return (1.0 - w) * points[lo].uncertainty + w * points[hi].uncertainty;
}

bool PropertySource::is_frequency_resolved() const {
  for (std::size_t k = 0; k < n_points; ++k)
    if (points[k].frequency_resolved()) return true;
  return false;
}

int PropertySource::provenance_completeness() const {
  return (states_method() ? 1 : 0) + (states_purity() ? 1 : 0) +
         (states_water_content() ? 1 : 0);
}

// ---------------------------------------------------------------------------

const PropertySource& MaterialProperty::selection() const {
  if (!has_selection())
    throw std::runtime_error(std::string("MissingMaterialData: fuer ") + to_string(kind) +
                             " ist keine Quelle ausgewaehlt.  " +
                             explain(MaterialDataStatus::MissingMaterialData));
  return sources[static_cast<std::size_t>(selected)];
}

Real MaterialProperty::min_at(Real T, Real tol, bool include_fr) const {
  Real m = std::numeric_limits<Real>::max();
  bool any = false;
  for (std::size_t k = 0; k < n_sources; ++k) {
    const PropertySource& s = sources[k];
    if (!include_fr && s.is_frequency_resolved()) continue;
    for (std::size_t j = 0; j < s.n_points; ++j)
      if (s.points[j].ambient() && std::abs(s.points[j].T - T) <= tol) {
        m = std::min(m, s.points[j].value);
        any = true;
      }
  }
  return any ? m : kNaN;
}

Real MaterialProperty::max_at(Real T, Real tol, bool include_fr) const {
  Real m = -std::numeric_limits<Real>::max();
  bool any = false;
  for (std::size_t k = 0; k < n_sources; ++k) {
    const PropertySource& s = sources[k];
    if (!include_fr && s.is_frequency_resolved()) continue;
    for (std::size_t j = 0; j < s.n_points; ++j)
      if (s.points[j].ambient() && std::abs(s.points[j].T - T) <= tol) {
        m = std::max(m, s.points[j].value);
        any = true;
      }
  }
  return any ? m : kNaN;
}

std::size_t MaterialProperty::n_sources_at(Real T, Real tol) const {
  std::size_t n = 0;
  for (std::size_t k = 0; k < n_sources; ++k) {
    const PropertySource& s = sources[k];
    if (s.is_frequency_resolved()) continue;
    for (std::size_t j = 0; j < s.n_points; ++j)
      if (s.points[j].ambient() && std::abs(s.points[j].T - T) <= tol) { ++n; break; }
  }
  return n;
}

Real MaterialProperty::relative_spread_at(Real T, Real tol) const {
  if (!has_selection()) return kNaN;
  const Real v = selection().value_at(T);
  const Real lo = min_at(T, tol), hi = max_at(T, tol);
  if (!(std::isfinite(v) && std::isfinite(lo) && std::isfinite(hi)) || !(std::abs(v) > 0.0))
    return kNaN;
  return (hi - lo) / std::abs(v);
}

// ---------------------------------------------------------------------------

const MaterialProperty* MaterialDataset::find(PropertyKind k) const {
  for (std::size_t i = 0; i < n_properties; ++i)
    if (properties[i].kind == k) return &properties[i];
  return nullptr;
}

void MaterialDataset::print(std::FILE* out) const {
  std::fprintf(out, "Stoffdatensatz %s (CAS %s), M = %.6g kg/mol\n", name, cas, molar_mass);
  std::fprintf(out, "  Stoffidentitaet: %s\n", identity_source);
  for (std::size_t i = 0; i < n_properties; ++i) {
    const MaterialProperty& p = properties[i];
    std::fprintf(out, "  %-24s %2zu Quellen", to_string(p.kind), p.n_sources);
    if (p.has_selection()) {
      const PropertySource& s = p.selection();
      const Real v = s.value_at(298.15);
      std::fprintf(out, ", gewaehlt: %.6g %s bei 298,15 K\n", v, si_unit(p.kind));
      std::fprintf(out, "      %s\n", s.reference);
      std::fprintf(out, "      Methode %s | Reinheit %s | Wasser %s\n", s.method,
                   s.states_purity() ? s.purity : "NICHT ANGEGEBEN",
                   s.states_water_content() ? s.water_content : "NICHT ANGEGEBEN");
      const Real spread = p.relative_spread_at(298.15);
      if (std::isfinite(spread))
        std::fprintf(out, "      Literaturstreuung bei 298,15 K: %.1f %% ueber %zu Quellen\n",
                     100.0 * spread, p.n_sources_at(298.15));
    } else {
      std::fprintf(out, ", KEINE ausgewaehlt -> MissingMaterialData\n");
    }
  }
}

// ---------------------------------------------------------------------------

Real MaterialValue::value_or_throw() const {
  if (!usable()) throw std::runtime_error("MissingMaterialData: " + message);
  return value;
}

void MaterialValue::print(std::FILE* out) const {
  if (!usable()) {
    std::fprintf(out, "  MissingMaterialData: %s\n", message.c_str());
    return;
  }
  std::fprintf(out, "  %.9g bei %.2f K [%s]\n", value, T, to_string(status));
  if (source) std::fprintf(out, "    %s\n", source->reference);
  if (relative_spread > 0.0)
    std::fprintf(out, "    Literaturstreuung %.1f %%\n", 100.0 * relative_spread);
}

// ---------------------------------------------------------------------------
// What a change of gamma does
// ---------------------------------------------------------------------------

const char* to_string(GammaScalingCategory c) {
  switch (c) {
    case GammaScalingCategory::DimensionlessEquilibrium:
      return "scaling_of_dimensionless_equilibrium";
    case GammaScalingCategory::DimensionlessGroup:
      return "dimensionless_group_at_fixed_field";
    case GammaScalingCategory::InvariantAtFixedState:
      return "invariant_at_fixed_geometry_voltage_permittivity";
  }
  return "?";
}

namespace {
const GammaScalingRow kGammaRows[] = {
    {"capillary_pressure_scale_gamma_over_a", GammaScalingCategory::DimensionlessEquilibrium,
     "gamma^1", 1.0, false, "die dimensionslose Form und der Kruemmungsradius a",
     "Die Druckskala des Gleichgewichts ist gamma/a.  Bei fester dimensionsloser Form "
     "wandert jeder Kapillardruck proportional zu gamma."},
    {"mechanical_load_for_same_dimensionless_shape",
     GammaScalingCategory::DimensionlessEquilibrium, "gamma^1", 1.0, false,
     "die dimensionslose Form",
     "Der mechanische Druck beziehungsweise die Gleichgewichtslast, die DIESELBE "
     "dimensionslose Form traegt, skaliert mit derselben Druckskala gamma/a."},
    {"voltage_for_same_dimensionless_shape", GammaScalingCategory::DimensionlessEquilibrium,
     "sqrt(gamma)", 0.5, false, "die dimensionslose Form und die Geometrie",
     "Elektrokapillares Aehnlichkeitsargument: bei fester dimensionsloser Form ist die "
     "elektrische Bondzahl Gamma = eps0 E^2 a / (2 gamma) fest, also E ~ sqrt(gamma) und "
     "damit U ~ sqrt(gamma)."},
    {"electric_bond_number_at_fixed_field", GammaScalingCategory::DimensionlessGroup,
     "1/gamma", -1.0, false, "das Feld E und die Geometrie",
     "Bei festgehaltenem FELD steht gamma im Nenner: Gamma = eps0 E^2 a / (2 gamma).  Ein "
     "groesseres gamma macht dieselbe Feldstaerke elektrisch weniger wirksam."},
    {"maxwell_traction_at_fixed_geometry_and_voltage",
     GammaScalingCategory::InvariantAtFixedState, "gamma^0", 0.0, false,
     "Geometrie, angelegte Spannung und Permittivitaetsverteilung",
     "KEINE Skalierung mit gamma.  Bei festgehaltener Geometrie, Spannung und "
     "Permittivitaetsverteilung folgt das Feld aus der Laplace- beziehungsweise "
     "Poisson-Gleichung, in der gamma nirgends vorkommt; die Maxwell-Traktion eps0 E^2 / 2 "
     "ist ein Funktional allein dieses Feldes und bleibt unveraendert.  gamma entscheidet, "
     "WELCHE Form im Gleichgewicht steht -- nicht, welche Kraft ein gegebenes Feld auf eine "
     "gegebene Form ausuebt.  Diese Zeile stand frueher faelschlich als 'linear in gamma' "
     "in impact.csv."},
};
}  // namespace

const GammaScalingRow* gamma_scaling_rows(std::size_t& n) {
  n = sizeof(kGammaRows) / sizeof(kGammaRows[0]);
  return kGammaRows;
}

// ---------------------------------------------------------------------------
// nu = mu / rho
// ---------------------------------------------------------------------------

void KinematicViscosityDerivation::print(std::FILE* out) const {
  if (!ok) {
    std::fprintf(out, "  nu ist NICHT ableitbar: %s\n", blocker.c_str());
    return;
  }
  std::fprintf(out, "  nu = %.9g m^2/s bei %.2f K [abgeleitet, %s]\n", value, T, identity);
  std::fprintf(out, "    mu  = %.9g", mu);
  if (std::isfinite(mu_uncertainty)) std::fprintf(out, " +- %.3g", mu_uncertainty);
  std::fprintf(out, " Pa s  -- %s\n", mu_source ? mu_source->reference : "?");
  std::fprintf(out, "    rho = %.9g", rho);
  if (std::isfinite(rho_uncertainty)) std::fprintf(out, " +- %.3g", rho_uncertainty);
  std::fprintf(out, " kg/m^3 -- %s\n", rho_source ? rho_source->reference : "?");
  std::fprintf(out, "    Bedingungen: %s\n", conditions.c_str());
  if (uncertainty_propagated)
    std::fprintf(out,
                 "    fortgepflanzte Unsicherheit: %.3g m^2/s (%.2f %%), linear %.3g m^2/s\n",
                 uncertainty, 100.0 * relative_uncertainty, uncertainty_linear);
  else
    std::fprintf(out, "    Unsicherheit NICHT fortpflanzbar: eine Elternquelle gibt keine an\n");
  std::fprintf(out, "    dieselbe Publikation fuer mu und rho: %s\n",
               same_publication ? "ja" : "nein");
  std::fprintf(out, "    NICHT gemessen -- die Groesse selbst wurde nicht bestimmt.\n");
}

KinematicViscosityDerivation derive_kinematic_viscosity(const MaterialDataset& d, Real T) {
  KinematicViscosityDerivation r;
  r.T = T;

  const MaterialProperty* pm = d.find(PropertyKind::DynamicViscosity);
  const MaterialProperty* pr = d.find(PropertyKind::Density);

  // (C1) both parents must have a selection.
  if (pm == nullptr || !pm->has_selection()) {
    r.blocker = "C1 verletzt: fuer die dynamische Viskositaet ist keine Quelle ausgewaehlt, "
                "die die Auswahlregel erfuellt.";
    return r;
  }
  if (pr == nullptr || !pr->has_selection()) {
    r.blocker = "C1 verletzt: fuer die Dichte ist keine Quelle ausgewaehlt, die die "
                "Auswahlregel erfuellt.";
    return r;
  }
  const PropertySource& sm = pm->selection();
  const PropertySource& sr = pr->selection();
  r.mu_source = &sm;
  r.rho_source = &sr;

  // (C2) both must cover T without extrapolation.
  const Real mu = sm.value_at(T);
  const Real rho = sr.value_at(T);
  char buf[512];
  if (!std::isfinite(mu)) {
    std::snprintf(buf, sizeof buf,
                  "C2 verletzt: die gewaehlte mu-Quelle deckt %.2f K nicht ab (gemessen "
                  "%.2f .. %.2f K).  Es wird nicht extrapoliert.",
                  T, sm.T_min(), sm.T_max());
    r.blocker = buf;
    return r;
  }
  if (!std::isfinite(rho)) {
    std::snprintf(buf, sizeof buf,
                  "C2 verletzt: die gewaehlte rho-Quelle deckt %.2f K nicht ab (gemessen "
                  "%.2f .. %.2f K).  Es wird nicht extrapoliert.",
                  T, sr.T_min(), sr.T_max());
    r.blocker = buf;
    return r;
  }
  if (!(rho > 0.0)) {
    r.blocker = "C2 verletzt: die Dichte bei dieser Temperatur ist nicht positiv.";
    return r;
  }

  // (C3) ambient on both sides is guaranteed by value_at(), which never reads a
  // non-ambient point.  It is asserted here so that the condition is visible in
  // the code rather than only in a comment.
  if (!sm.has_ambient_points() || !sr.has_ambient_points()) {
    r.blocker = "C3 verletzt: eine der beiden Quellen hat keine Punkte bei Umgebungsdruck.";
    return r;
  }
  if (sm.is_frequency_resolved() || sr.is_frequency_resolved()) {
    r.blocker = "C3 verletzt: eine der beiden Quellen ist frequenzaufgeloest.";
    return r;
  }

  // (C4) the documented sample conditions must agree verbatim.
  auto same = [](const char* a, const char* b) { return std::string(a) == std::string(b); };
  if (!same(sm.purity, sr.purity)) {
    std::snprintf(buf, sizeof buf,
                  "C4 verletzt: die Reinheit ist verschieden angegeben -- mu-Quelle '%s', "
                  "rho-Quelle '%s'.  Zwei verschiedene Proben ergeben keinen gemeinsamen "
                  "Quotienten.",
                  sm.purity, sr.purity);
    r.blocker = buf;
    return r;
  }
  if (!same(sm.water_content, sr.water_content)) {
    std::snprintf(buf, sizeof buf,
                  "C4 verletzt: der Wassergehalt ist verschieden angegeben -- mu-Quelle '%s', "
                  "rho-Quelle '%s'.",
                  sm.water_content, sr.water_content);
    r.blocker = buf;
    return r;
  }
  if (!same(sm.sample_source, sr.sample_source)) {
    std::snprintf(buf, sizeof buf,
                  "C4 verletzt: die Probenherkunft ist verschieden angegeben -- mu-Quelle "
                  "'%s', rho-Quelle '%s'.",
                  sm.sample_source, sr.sample_source);
    r.blocker = buf;
    return r;
  }

  // All conditions hold.
  r.ok = true;
  r.mu = mu;
  r.rho = rho;
  r.value = mu / rho;
  r.mu_uncertainty = sm.uncertainty_at(T);
  r.rho_uncertainty = sr.uncertainty_at(T);
  r.same_publication = same(sm.reference, sr.reference);

  if (std::isfinite(r.mu_uncertainty) && std::isfinite(r.rho_uncertainty)) {
    const Real em = r.mu_uncertainty / mu;
    const Real er = r.rho_uncertainty / rho;
    r.relative_uncertainty = std::sqrt(em * em + er * er);
    r.uncertainty = r.value * r.relative_uncertainty;
    r.uncertainty_linear = r.value * (em + er);
    r.uncertainty_propagated = true;
  } else {
    r.relative_uncertainty = kNaN;
    r.uncertainty = kNaN;
    r.uncertainty_linear = kNaN;
    r.uncertainty_propagated = false;
  }

  std::snprintf(buf, sizeof buf,
                "T = %.2f K, Umgebungsdruck, Reinheit '%s', Wassergehalt '%s', Probe '%s'; "
                "Methoden: mu %s, rho %s",
                T, sm.purity, sm.water_content, sm.sample_source, sm.method, sr.method);
  r.conditions = buf;
  return r;
}

MaterialValue derived_kinematic_viscosity(const MaterialDataset& d, Real T) {
  const KinematicViscosityDerivation r = derive_kinematic_viscosity(d, T);
  MaterialValue v;
  v.T = T;
  if (!r.ok) {
    v.message = "kinematic_viscosity: " + r.blocker;
    return v;
  }
  v.value = r.value;
  v.status = MaterialDataStatus::Derived;
  // The parent that dominates the error bar is the one a reader should look up
  // first; naming it here keeps a single-source field usable without pretending
  // the other parent does not exist -- the full record is the derivation.
  v.source = r.mu_source;
  v.relative_spread = std::isfinite(r.relative_uncertainty) ? r.relative_uncertainty : 0.0;
  v.message = std::string("abgeleitet: ") + r.identity + "; " + r.conditions;
  return v;
}

MaterialValue material_value(const MaterialDataset& d, PropertyKind kind, Real T) {
  MaterialValue v;
  v.T = T;
  const MaterialProperty* p = d.find(kind);
  if (p == nullptr) {
    v.message = std::string("Der Datensatz '") + d.name + "' fuehrt " + to_string(kind) +
                " gar nicht.";
    return v;
  }
  if (!p->has_selection()) {
    v.message = std::string(to_string(kind)) + ": keine Quelle erfuellt die Auswahlregel " +
                "(Methode, Reinheit UND Wassergehalt angegeben).  Es wird kein Ersatzwert "
                "gesetzt.";
    if (kind == PropertyKind::KinematicViscosity)
      v.message += "  Es gibt hier also keine DIREKTE Messung, die der Vertrag zulaesst; "
                   "ein Wert kann aber aus mu und rho ABGELEITET werden -- siehe "
                   "derive_kinematic_viscosity(), das seine Bedingungen einzeln prueft.";
    return v;
  }
  const PropertySource& s = p->selection();
  const Real value = s.value_at(T);
  if (!std::isfinite(value)) {
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "%s: die gewaehlte Quelle deckt %.2f K nicht ab (gemessen %.2f .. %.2f K).  "
                  "Es wird NICHT extrapoliert.",
                  to_string(kind), T, s.T_min(), s.T_max());
    v.message = buf;
    return v;
  }
  v.value = value;
  v.status = s.status;
  v.source = &s;
  const Real spread = p->relative_spread_at(T);
  v.relative_spread = std::isfinite(spread) ? spread : 0.0;
  v.message = "ok";
  return v;
}

}  // namespace es
