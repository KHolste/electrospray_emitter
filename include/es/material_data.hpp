#pragma once
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "es/types.hpp"

namespace es {

// ===========================================================================
// P2 -- a material-data contract with provenance and a validity range
// ===========================================================================
//
// WHY THE EXISTING DATA SET IS NOT ENOUGH
//
// es::Fluid carries numbers without a source.  es::LiquidProperties carries a
// status with them and says `Illustrative`, which is honest but still leaves
// the solver with a number that nothing supports.  Surface tension enters the
// P3a/P3b equilibrium LINEARLY: every pressure scale, every apex height and --
// through V^2 ~ gamma -- every voltage of that phase is proportional to it.  A
// value that is wrong by fifteen per cent makes every number of those phases
// wrong by fifteen per cent, and no amount of mesh refinement helps.
//
// WHAT THIS FILE IS
//
// A typed record of MEASUREMENTS, one per source, each carrying:
//
//   * the full primary citation, so the value can be found again;
//   * the measurement method, verbatim;
//   * the sample provenance -- source, purification, purity, WATER CONTENT --
//     verbatim, including the fact that a source did not state it;
//   * the temperature at which each point was measured, and the range the
//     source covers;
//   * the reported uncertainty where the source gave one;
//   * the measurement FREQUENCY where the quantity is frequency dependent, so
//     that a value measured at 18 GHz cannot be quoted for a DC problem.
//
// WHAT THIS FILE IS NOT.  It does not average.  Sources that disagree stand
// next to each other, and the spread is a reportable quantity, not something to
// be hidden in a mean.  The one value a solver uses is SELECTED by a stated
// rule and named in every output.
//
// FAIL CLOSED.  A property with no selected source yields MissingMaterialData.
// There is no silent default anywhere in this file.

// ---------------------------------------------------------------------------

enum class PropertyKind {
  SurfaceTension = 0,     ///< [N/m]  -- enters the P3a/P3b equilibrium
  Density,                ///< [kg/m^3]
  DynamicViscosity,       ///< [Pa s]
  KinematicViscosity,     ///< [m^2/s] -- reported separately, never converted
                          ///< silently: converting needs a density, and which
                          ///< density the original author used is not recorded.
  ElectricalConductivity, ///< [S/m]
  RelativePermittivity,   ///< [-]
};
const char* to_string(PropertyKind k);
const char* si_unit(PropertyKind k);

enum class MaterialDataStatus {
  /// A primary source has been read at the place the value is taken from, and
  /// substance, temperature and method are unambiguous.
  Measured = 0,
  /// A manufacturer data sheet: the document is in hand and the value is
  /// stated with its temperature, but the method and the water content are
  /// usually not.
  ManufacturerSpec,
  /// A named source that has not been read at the value.
  Literature,
  /// An example value with no primary source.  It may carry a dimensionless
  /// demonstration and nothing else.
  Illustrative,
  /// Registered without a number.  Using it is an error, not a default.
  MissingMaterialData,
};
const char* to_string(MaterialDataStatus s);
const char* explain(MaterialDataStatus s);
inline bool carries_quantitative_claim(MaterialDataStatus s) {
  return s == MaterialDataStatus::Measured || s == MaterialDataStatus::ManufacturerSpec;
}

// ---------------------------------------------------------------------------

/// Ambient pressure, and the band around it a point may lie in and still count
/// as an ambient measurement.  A density measured at 60 MPa is not the density
/// at 1 atm, and ILThermo carries both.
inline constexpr Real kAmbientPressure = 101325.0;   ///< [Pa]
inline constexpr Real kAmbientTolerance = 0.05;      ///< relative

struct PropertyPoint {
  Real T{0};             ///< [K]
  Real value{0};         ///< in si_unit(kind)
  Real uncertainty{0};   ///< as reported; 0 means "not reported", never "zero"
  Real frequency_Hz{0};  ///< 0 means "not frequency resolved" (a static value)
  /// Measurement pressure [Pa].  0 means the source stated none, which for
  /// these compilations means ambient -- so 0 counts as ambient, and a stated
  /// 60 MPa does not.
  Real pressure_Pa{0};
  bool has_uncertainty() const { return uncertainty > 0.0; }
  bool frequency_resolved() const { return frequency_Hz > 0.0; }
  bool ambient() const {
    return pressure_Pa <= 0.0 ||
           (pressure_Pa > (1.0 - kAmbientTolerance) * kAmbientPressure &&
            pressure_Pa < (1.0 + kAmbientTolerance) * kAmbientPressure);
  }
};

/// One source's measurements of one property.  Everything a reader would have
/// to look up is here, verbatim; nothing is paraphrased and nothing is filled
/// in.  An empty string means the source did not state it.
struct PropertySource {
  PropertyKind kind{PropertyKind::SurfaceTension};
  const char* database_id{""};   ///< ILThermo set id, or another retrieval key
  const char* reference{""};     ///< full primary citation
  const char* paper_title{""};
  const char* method{""};        ///< measurement method, verbatim
  const char* sample_source{""};
  const char* purity{""};        ///< verbatim; empty = not stated
  const char* water_content{""}; ///< verbatim; empty = not stated
  const char* constraints{""};   ///< e.g. "Pressure of 1 atm"
  MaterialDataStatus status{MaterialDataStatus::Measured};
  const PropertyPoint* points{nullptr};
  std::size_t n_points{0};

  /// Over the AMBIENT points only; see PropertyPoint::ambient().
  Real T_min() const;
  Real T_max() const;
  bool covers(Real T) const;
  /// Number of points measured away from ambient pressure.  They are kept in
  /// the record -- they are real measurements -- but they never back an
  /// ambient-pressure number.
  std::size_t n_non_ambient() const;
  bool has_ambient_points() const;
  /// Linear interpolation between the two neighbouring points.  Returns NaN
  /// outside [T_min, T_max]: extrapolating a fit through measurements that were
  /// not taken is exactly what a provenance record is supposed to prevent.
  Real value_at(Real T) const;
  /// True where at least one point carries a non-zero frequency, i.e. the
  /// source measured a frequency-resolved quantity.  Such a source must not
  /// back a DC number without saying so.
  bool is_frequency_resolved() const;
  bool states_purity() const { return purity[0] != '\0'; }
  bool states_water_content() const { return water_content[0] != '\0'; }
  bool states_method() const { return method[0] != '\0'; }
  /// How completely the provenance is documented: method, purity and water
  /// content, each worth one point.  Used by the SELECTION RULE below, which is
  /// stated rather than applied silently.
  int provenance_completeness() const;
};

/// All the sources for one property, plus the one that is selected.
struct MaterialProperty {
  PropertyKind kind{PropertyKind::SurfaceTension};
  const PropertySource* sources{nullptr};
  std::size_t n_sources{0};
  /// Index into `sources` of the value a computation uses, or -1 for none.
  ///
  /// THE SELECTION RULE, stated once and applied mechanically:
  ///   1. the source must state method, purity AND water content;
  ///   2. it must have AMBIENT-pressure points and must not be frequency
  ///      resolved -- a density at 60 MPa and a permittivity at 18 GHz are
  ///      real measurements of something else;
  ///   3. among those, the one whose ambient temperature range contains
  ///      298.15 K and which has the most ambient points;
  ///   4. ties broken by the smaller reported uncertainty.
  /// If no source satisfies (1), nothing is selected and the property is
  /// MissingMaterialData.  The rule is not tuned to produce a particular
  /// number; it is applied by tools/fetch_material_data.py and re-checked in
  /// tests/test_material_data.cpp.
  int selected{-1};

  bool has_selection() const { return selected >= 0; }
  const PropertySource& selection() const;   ///< throws if none

  /// Smallest and largest value any source reports within `tol` of `T`.
  /// Sources that are frequency resolved are EXCLUDED from this band unless
  /// `include_frequency_resolved` -- otherwise a GHz measurement would widen
  /// the band of a DC quantity.
  Real min_at(Real T, Real tol = 2.0, bool include_frequency_resolved = false) const;
  Real max_at(Real T, Real tol = 2.0, bool include_frequency_resolved = false) const;
  /// (max - min) / selected value at T.  This is the LITERATURE SCATTER and it
  /// is the honest uncertainty of the number, far larger than any single
  /// source's own error bar.
  Real relative_spread_at(Real T, Real tol = 2.0) const;
  std::size_t n_sources_at(Real T, Real tol = 2.0) const;
};

/// A substance: the properties, and what is known about the substance itself.
struct MaterialDataset {
  const char* name{""};
  const char* cas{""};
  Real molar_mass{0};          ///< [kg/mol], 0 = not recorded
  const char* identity_source{""};  ///< where the substance identity comes from
  const char* note{""};
  const MaterialProperty* properties{nullptr};
  std::size_t n_properties{0};

  const MaterialProperty* find(PropertyKind k) const;
  void print(std::FILE* out) const;
};

// ---------------------------------------------------------------------------
// Querying, fail-closed
// ---------------------------------------------------------------------------

struct MaterialValue {
  Real value{0};
  MaterialDataStatus status{MaterialDataStatus::MissingMaterialData};
  const PropertySource* source{nullptr};
  Real T{0};
  Real relative_spread{0};   ///< literature scatter at T, or 0 if unknown
  std::string message;

  bool usable() const { return status != MaterialDataStatus::MissingMaterialData; }
  /// Throws std::runtime_error naming the property and the reason.  For callers
  /// that must stop rather than continue with a number they do not have.
  Real value_or_throw() const;
  void print(std::FILE* out) const;
};

/// The selected value of `kind` at temperature `T`.  Returns
/// MissingMaterialData -- never a default -- if the property has no selection,
/// if the selected source does not cover T, or if the dataset does not carry
/// the property at all.
MaterialValue material_value(const MaterialDataset& d, PropertyKind kind, Real T);

// ---------------------------------------------------------------------------
// The datasets
// ---------------------------------------------------------------------------

/// EMI-BF4 from the NIST ILThermo compilation (Standard Reference Database
/// #147) plus one manufacturer data sheet.  Generated by
/// tools/fetch_material_data.py into src/material_data_emibf4.cpp; the audit
/// trail is docs/13_material_data.md.
const MaterialDataset& emibf4_sourced();

}  // namespace es
