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
  KinematicViscosity,     ///< [m^2/s] -- kept as its own record of DIRECT
                          ///< measurements.  It is never silently converted
                          ///< from a foreign mu and rho, because which density
                          ///< the original author divided by is not recorded.
                          ///< A value may however be DERIVED from this data
                          ///< set's own selected mu and rho -- see
                          ///< derive_kinematic_viscosity() below, which states
                          ///< its compatibility conditions and fails closed.
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
  /// NOT measured.  Computed from other SELECTED properties of THIS data set by
  /// a stated identity, from sources whose documented conditions were checked
  /// for compatibility, with the uncertainty propagated.  It carries both
  /// parent citations and both parent values and must never be presented as if
  /// the quantity itself had been measured.
  Derived,
  /// An example value with no primary source.  It may carry a dimensionless
  /// demonstration and nothing else.
  Illustrative,
  /// Registered without a number.  Using it is an error, not a default.
  MissingMaterialData,
};
const char* to_string(MaterialDataStatus s);
const char* explain(MaterialDataStatus s);
inline bool carries_quantitative_claim(MaterialDataStatus s) {
  return s == MaterialDataStatus::Measured || s == MaterialDataStatus::ManufacturerSpec ||
         s == MaterialDataStatus::Derived;
}
/// True only where the quantity ITSELF was measured.  A derived value carries a
/// quantitative claim but is not a measurement of its own quantity, and any
/// output that distinguishes the two must ask this and not the above.
inline bool is_direct_measurement(MaterialDataStatus s) {
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
  /// The reported uncertainty at T, interpolated between the two neighbouring
  /// AMBIENT points exactly as value_at() interpolates the value.  Returns NaN
  /// outside the measured range, and NaN if either neighbouring point reports
  /// no uncertainty -- because uncertainty == 0 means "the source did not state
  /// one", never "the source stated zero", and averaging a stated figure with a
  /// missing one would invent an error bar.
  Real uncertainty_at(Real T) const;
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
// Deriving nu = mu / rho -- allowed, but only under stated conditions
// ---------------------------------------------------------------------------
//
// The kinematic viscosity is not a second, independent material property: it is
// mu / rho by definition.  Refusing to form it at all would be over-strict.
// Forming it from any mu and any rho would be worse: the two would then be a
// viscosity of one sample at one temperature divided by a density of a
// DIFFERENT sample, possibly at a different temperature, water content or
// pressure, and the quotient would describe no liquid that exists.
//
// So the derivation is permitted exactly when the two parent values can be
// shown to describe the same liquid in the same state.  The conditions are
// checked mechanically and named individually when one of them fails:
//
//   (C1) mu and rho each have a SELECTED source (i.e. each satisfies the
//        selection rule: method, purity and water content stated, ambient,
//        not frequency resolved);
//   (C2) both selected sources cover T without extrapolation;
//   (C3) the values used are both ambient-pressure values -- guaranteed by
//        value_at(), which never touches a non-ambient point;
//   (C4) the documented sample conditions agree VERBATIM between the two
//        sources: purity, water content and sample origin.  Two different
//        strings are not proof of incompatibility, but they are the absence of
//        proof of compatibility, and this contract fails closed on that.
//
// It is deliberately NOT required that the two sources be the same publication.
// Whether they are is recorded (`same_publication`) and reported, because it
// makes the claim stronger, not because the derivation would otherwise be void.
//
// The uncertainty is propagated through nu = mu/rho as independent relative
// errors.  mu and rho from one and the same instrument run are not strictly
// independent, and the quadratic sum is then the smaller of the two plausible
// combinations; the linear sum is also reported so that a reader can take the
// conservative one.  The propagated figure is dominated by mu either way.

struct KinematicViscosityDerivation {
  bool ok{false};              ///< false: nothing below except `blocker` is valid
  const char* identity{"nu = mu / rho"};
  Real T{0};                   ///< [K], the one temperature both values are at

  Real mu{0};                  ///< [Pa s]   the value actually used
  Real mu_uncertainty{0};      ///< [Pa s]   as reported at T; NaN = not reported
  Real rho{0};                 ///< [kg/m^3]
  Real rho_uncertainty{0};     ///< [kg/m^3]
  const PropertySource* mu_source{nullptr};
  const PropertySource* rho_source{nullptr};

  Real value{0};               ///< [m^2/s]  = mu / rho
  Real uncertainty{0};         ///< [m^2/s]  quadratic propagation; NaN if a
                               ///< parent reported none
  Real uncertainty_linear{0};  ///< [m^2/s]  linear propagation (conservative)
  Real relative_uncertainty{0};

  bool same_publication{false};
  bool uncertainty_propagated{false};
  std::string conditions;      ///< purity / water / pressure, verbatim
  std::string blocker;         ///< empty iff ok; names the failing condition

  void print(std::FILE* out) const;
};

/// nu at T, derived from the SELECTED dynamic viscosity and the SELECTED
/// density of `d`.  Fails closed -- ok == false and `blocker` naming the one
/// condition that prevents the derivation -- rather than returning a number.
KinematicViscosityDerivation derive_kinematic_viscosity(const MaterialDataset& d, Real T);

/// The same result in the common MaterialValue shape, with status Derived (or
/// MissingMaterialData with the blocker as its message).  Callers that only
/// want a number and a status use this; callers that must report HOW the number
/// arose use derive_kinematic_viscosity() and get both parents.
MaterialValue derived_kinematic_viscosity(const MaterialDataset& d, Real T);

// ---------------------------------------------------------------------------
// What a change of gamma does -- and what it demonstrably does NOT do
// ---------------------------------------------------------------------------
//
// A corrected surface tension changes P3a/P3b numbers.  How it changes them is
// a physical statement, and an earlier version of this table got one row of it
// wrong, so the rows now live here, in the library, where a test can check them
// instead of in an application's printf.
//
// THE DISTINCTION THIS TABLE EXISTS TO MAKE
//
//   (1) SCALING A DIMENSIONLESS EQUILIBRIUM.  The P3a/P3b solver computes a
//       nondimensional shape in which gamma appears only through the capillary
//       pressure scale gamma/a and the electric Bond number
//       Gamma = eps0 E^2 a / (2 gamma).  Holding the DIMENSIONLESS solution
//       fixed and changing gamma therefore moves the dimensional pressures and
//       voltages by an exact power of gamma.  Nothing is re-solved.  This is a
//       unit conversion of a result already computed, and it is exact.
//
//   (2) A RECOMPUTED COUPLED SOLUTION at fixed voltage and fixed geometry.  That
//       is a different question and this table does not answer it.  No row here
//       is a new simulation, and the `recomputed` field says so on every row.
//
// THE ROW THAT WAS WRONG, AND WHY
//
// A previous version listed "maxwell_force_for_fixed_shape" as "linear in
// gamma".  That is false.  With the geometry, the applied voltage and the
// permittivity distribution all held fixed, the field follows from Laplace's
// (or Poisson's) equation, in which gamma does not appear at ALL -- not in the
// operator, not in the boundary conditions, not in the coefficients.  The
// Maxwell traction eps0 E^2 / 2 is then a functional of that field alone and is
// INVARIANT under a change of gamma: exponent 0, factor exactly 1.
//
// gamma decides WHICH shape stands in equilibrium.  It does not decide the
// force a given field exerts on a given shape.  Conflating the two turned a
// statement about the equilibrium into a false statement about the electric
// stress, and it is kept in the table -- as an invariant, with its exponent of
// zero -- rather than deleted, so that the correction is visible.

enum class GammaScalingCategory {
  /// An exact scaling of the nondimensional P3a/P3b solution.  Not a new solve.
  DimensionlessEquilibrium = 0,
  /// A dimensionless group in which gamma appears explicitly, evaluated at a
  /// fixed FIELD rather than at a fixed shape.
  DimensionlessGroup,
  /// Does not scale with gamma at all.  Listed so the invariance is stated
  /// rather than left to be assumed.
  InvariantAtFixedState,
};
const char* to_string(GammaScalingCategory c);

struct GammaScalingRow {
  const char* quantity{""};
  GammaScalingCategory category{GammaScalingCategory::DimensionlessEquilibrium};
  const char* law{""};        ///< "gamma^1", "sqrt(gamma)", "gamma^0", "1/gamma"
  Real exponent{0};           ///< factor = (gamma_new / gamma_old)^exponent
  /// false on every row: nothing in this table is a newly computed coupled
  /// solution, and an output that presented one as such would be wrong.
  bool recomputed{false};
  const char* what_is_held_fixed{""};
  const char* note{""};
};

/// The complete table.  `n` receives the row count.
const GammaScalingRow* gamma_scaling_rows(std::size_t& n);

// ---------------------------------------------------------------------------
// The datasets
// ---------------------------------------------------------------------------

/// EMI-BF4 from the NIST ILThermo compilation (Standard Reference Database
/// #147) plus one manufacturer data sheet.  Generated by
/// tools/fetch_material_data.py into src/material_data_emibf4.cpp; the audit
/// trail is docs/13_material_data.md.
const MaterialDataset& emibf4_sourced();

}  // namespace es
