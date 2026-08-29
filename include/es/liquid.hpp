#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/status.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P3a -- working liquid as a value PLUS its provenance
// ===========================================================================
//
// WHY THIS IS NOT es::Fluid.
//
// The table in fluid.hpp carries numbers without a source.  Its own header says
// so ("literature-typical room-temperature values ... NOT a substitute for
// characterising the batch you actually spray"), but nothing in the type stops
// a caller from quoting one of those numbers as if it had been measured.  The
// static meniscus is the first place where a surface tension enters an answer
// that is supposed to be quantitative, so the status travels WITH the value and
// every output has to print it.
//
// WHAT THE KUNZE DOCUMENTS ACTUALLY SUPPLY -- checked, not remembered:
//
//   * The liquid is named.  KunzeFynn-2024-12-10.pdf, section 2.3.2 "Ionic
//     liquids for electrospray thrusters", printed page 28 (PDF page 36):
//     "EMI-BF4 and EMI-Im were also selected for this project."  The table
//     "List of Publications", printed page 30 (PDF page 38), assigns EMI-BF4 as
//     the propellant to publications I-IV, i.e. to exactly the straight
//     10 um capillaries in SU-8 and IP-Q that the reference geometry models.
//     The IDENTITY of the liquid is therefore documented.
//
//   * No property value is.  The full text was searched for numeric surface
//     tension, density, viscosity and conductivity of either liquid; it
//     contains none, and no property table for the propellants exists in the
//     document.  Section 2.3.2 compares the two liquids only qualitatively
//     ("The surface tension of EMI-Im is slightly lower than that of EMI-BF4").
//
// Consequence, and it is deliberate: the built-in EMI-BF4 entry below carries
// status Illustrative.  Nothing was added from memory, no digit was invented,
// and the numbers are the ones already in es::Fluid -- carried over unchanged,
// with the fact that they have no primary source written into the record
// instead of left out of it.  Every quantitative statement made in P3a is
// therefore either dimensionless or explicitly marked as an example.

/// How much weight a property value carries.  Ordered strongest to weakest.
enum class LiquidDataStatus {
  /// Substance, temperature and primary source are unambiguous, and the source
  /// has been read at the place the value is taken from.
  Verified = 0,
  /// A named but unchecked source, or a value measured on a nominally equal
  /// substance in an unstated state.  Usable for a sensitivity study.
  Provisional,
  /// An example value with no primary source.  It may carry a dimensionless
  /// demonstration and nothing else.
  Illustrative,
  /// Registered without a number.  Using it is an error, not a default.
  Unknown,
};
const char* to_string(LiquidDataStatus s);
const char* explain(LiquidDataStatus s);
LiquidDataStatus liquid_status_from_string(const std::string& s);  ///< throws if unknown

/// True where a value may back a quantitative claim without a caveat.
inline bool carries_quantitative_claim(LiquidDataStatus s) {
  return s == LiquidDataStatus::Verified;
}

// ---------------------------------------------------------------------------

struct LiquidProperties {
  std::string substance{"unbenannt"};
  Real T{298.15};        ///< temperature the values below apply at [K]
  Real gamma{0.0};       ///< surface tension [N/m]     -- enters the equilibrium
  Real rho{0.0};         ///< density [kg/m^3]          -- Bond number only
  LiquidDataStatus status{LiquidDataStatus::Unknown};
  std::string source;    ///< precise enough to find the value again
  std::string caveat;    ///< what must NOT be concluded from it

  /// Properties that are documented here so they need not be looked up twice,
  /// and that are NOT inputs of the P3a equilibrium problem.  Nothing in this
  /// struct is read by the solver.  A value of 0 means "not supplied"; a value
  /// present here must never be reported as physics that has been accounted
  /// for.
  struct NotUsedInP3a {
    Real mu{0.0};     ///< dynamic viscosity [Pa s] -- needs flow (later phase)
    Real K{0.0};      ///< electrical conductivity [S/m] -- needs charge transport
    Real eps_r{0.0};  ///< relative permittivity [-] -- needs the field coupling
  } documented_only;

  /// Bond number rho g a^2 / gamma for a pinning radius a [m].  The ratio of
  /// the hydrostatic pressure over one bore radius to the capillary pressure
  /// scale gamma/a.  Gravity is NOT part of the solved equilibrium; this number
  /// is what justifies leaving it out, and it is reported, never assumed.
  Real bond_number(Real a) const;

  /// Capillary pressure scale gamma / a [Pa].
  Real capillary_pressure_scale(Real a) const;

  /// Empty unless the data set cannot be used: the returned string names the
  /// first defect in plain language.  Fail-closed, and no exception, so a
  /// solver can turn it into a status instead of a crash.
  std::string why_unusable() const;
  bool usable() const { return why_unusable().empty(); }
  /// Same test, as an exception, for applications that must stop at once.
  void validate_or_throw() const;

  void print(std::FILE* out) const;
};

// ---------------------------------------------------------------------------

/// EMI-BF4 as an EXPLICITLY ILLUSTRATIVE data set.  See the header comment: the
/// substance is documented in the Kunze dissertation, the numbers are not.
LiquidProperties emibf4_illustrative();

/// A liquid whose surface tension is exactly 1 N/m, for dimensionless tests.
/// Status Illustrative -- it is not a substance.
LiquidProperties unit_liquid();

/// Built-in data sets by name ("emi-bf4", "unit").  Throws if unknown.
LiquidProperties liquid_data_by_name(const std::string& name);
std::vector<std::string> liquid_data_names();

}  // namespace es
