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
