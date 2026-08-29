#pragma once
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "es/types.hpp"

namespace es {

// ===========================================================================
// P5 -- the ion-emission contract.  DISABLED BY DEFAULT, and blocked.
// ===========================================================================
//
// WHAT THIS FILE IS
//
// An API for ion emission with every physical input made explicit -- species,
// charge SIGN, mass, activation barrier, temperature, field direction -- and a
// fail-closed path for every one of them.  It ships with NO enabled model.
//
// WHY IT IS BLOCKED, and this is a finding, not an omission
//
// The task allows an emission model only from a FULLY CHECKED equation with
// SOURCED parameters.  Neither was obtainable in this run:
//
//   * THE EQUATION.  The Iribarne-Thomson ion-evaporation rate is cited
//     everywhere as
//         j = (kT/h) sigma_s exp[ -(dG - sqrt(e^3 E /(4 pi eps0))) / (kT) ] ,
//     and that form is what src/emission.cpp already computes.  The primary
//     sources are Iribarne & Thomson (1976) J. Chem. Phys. 64, 2287-2294 and
//     Thomson & Iribarne (1979) J. Chem. Phys. 71, 4451-4463.  Neither full
//     text was reachable in this run, so the PREFACTOR, the exact definition of
//     sigma_s and the stated validity range were NOT read at the source.  A
//     secondary source that was read (Wiley, "A Brief Overview of the
//     Mechanisms Involved in Electrospray Mass Spectrometry", section 1.2.8)
//     cites both papers, describes the model qualitatively, and states that it
//     "is experimentally well supported for small ions of the kind that one
//     encounters in inorganic and organic chemistry" -- with an explicit
//     reservation for larger ions, which is precisely the ionic-liquid case.
//     It prints no equation.  So the equation is UNVERIFIED here.
//
//   * THE PARAMETERS.  No sourced activation barrier dG for EMI-BF4 was found.
//     The literature that was reached says the thermochemical data needed to
//     evaluate the theoretical rate "are of insufficient accuracy".  The value
//     1.09 eV that sits in src/fluid.cpp has no source at all and its own
//     header calls it "the least certain quantity here by far", spanning
//     roughly 1.0-1.4 eV depending on extraction method.  Since j depends on dG
//     EXPONENTIALLY, a 0.4 eV span at 298 K is a factor of e^(0.4/0.0257) =
//     1e7 in the current.  A model with that parameter unsourced does not
//     predict a current; it reports the parameter.
//
// CONSEQUENCE, applied mechanically: emission_status() returns Blocked, the
// model is disabled, and emitted_current_density() returns
// MissingEmissionParameters.  There is no default barrier, no default species
// and no fallback to the diagnostic numbers of the old es::Fluid path.
//
// THREE FURTHER PROHIBITIONS, honoured by there being no code for them:
//   * no cone-jet correlation is used as an ion-emission current.  Those are
//     different physics and the old block stays disabled (see P8);
//   * no current computed from a PERFECT-CONDUCTOR field is called a
//     prediction.  Without the finite-conductivity back-reaction of P3/P4 the
//     surface field at the emission site is not the field the emission would
//     see;
//   * the two polarities are never computed from the same species data.

// ---------------------------------------------------------------------------

enum class Polarity { Positive = 0, Negative };
const char* to_string(Polarity p);

enum class EmissionStatus {
  /// The model is switched off.  This is the DEFAULT and the shipped state.
  Disabled = 0,
  /// A required parameter is absent.  No number is produced and no default is
  /// substituted.
  MissingEmissionParameters,
  /// The rate equation itself has not been verified against a primary source in
  /// this project.  Even with complete parameters, nothing is evaluated.
  EquationNotValidated,
  /// The species sign does not match the requested polarity, or the field
  /// points the wrong way to emit it.
  PolarityMismatch,
  /// A parameter lies outside the range the equation was stated for.
  ParameterOutOfRange,
  /// Everything checked out.  Unreachable with the shipped data, by design.
  Ok,
};
const char* to_string(EmissionStatus s);
const char* explain(EmissionStatus s);
inline bool is_usable(EmissionStatus s) { return s == EmissionStatus::Ok; }

// ---------------------------------------------------------------------------

/// One emitted species.  Every field is explicit; an empty source string means
/// the value is NOT sourced, and a value of zero means it is absent.
struct EmittedSpecies {
  const char* name{""};
  /// Signed charge number.  +1 for a cation, -1 for an anion.  The sign is
  /// carried, never taken from a magnitude: the prototype used |E_n| and cation
  /// masses for both polarities and therefore reported identical currents for
  /// signs that are physically different.
  int charge_number{0};
  Real mass{0};                    ///< [kg]
  const char* mass_source{""};
  Real activation_barrier{0};      ///< dG [J]
  const char* barrier_source{""};
  Real barrier_uncertainty{0};     ///< [J]; 0 means "not stated"
  const char* note{""};

  bool has_charge() const { return charge_number != 0; }
  bool has_mass() const { return mass > 0.0 && mass_source[0] != '\0'; }
  bool has_barrier() const { return activation_barrier > 0.0 && barrier_source[0] != '\0'; }
  bool complete() const { return has_charge() && has_mass() && has_barrier(); }
  Polarity polarity() const { return charge_number > 0 ? Polarity::Positive : Polarity::Negative; }
};

struct EmissionModel {
  /// OFF.  Turning it on is not enough to get a number; the checks below still
  /// have to pass, and with the shipped data they cannot.
  bool enabled{false};
  /// Set only by a source audit that read the primary paper.  It is false for
  /// this project and the reason is in the header comment.
  bool equation_validated{false};
  Real temperature{0};             ///< [K]
  Polarity polarity{Polarity::Positive};
  const EmittedSpecies* species{nullptr};
  std::size_t n_species{0};

  /// The first reason this model cannot produce a number, or Ok.
  EmissionStatus status() const;
  void print(std::FILE* out) const;
};

struct EmissionResult {
  EmissionStatus status{EmissionStatus::Disabled};
  Real current_density{0};   ///< [A/m^2]; 0 ONLY when the status says so
  Real barrier_lowering{0};  ///< [J]
  Real effective_barrier{0}; ///< [J]
  std::string message;
  bool usable() const { return is_usable(status); }
};

/// Emitted current density at a surface with outward normal field `E_n`.
/// Fails closed.  With the shipped data it always returns Disabled or
/// MissingEmissionParameters -- never a number.
EmissionResult emitted_current_density(const EmissionModel& m, Real E_n);

// ---------------------------------------------------------------------------
// The mathematical kernel, separated from the physical model on purpose
// ---------------------------------------------------------------------------
//
// The functions below are PURE FUNCTIONS of their arguments.  They implement
// the functional form that is commonly quoted for the Iribarne-Thomson rate.
// They are separated from EmissionModel so that the FORM can be tested --
// dimensions, monotonicity, limits, the sign of the field dependence -- without
// any material data and without anything being called a prediction.
//
// Testing the form is not validating the model.  What the tests below establish
// is that the code computes the function it says it computes; whether that
// function describes ionic-liquid emission is exactly the open question.

/// Schottky-type barrier lowering by a field of magnitude E [J]:
///     G(E) = sqrt( e^3 E / (4 pi eps0) ) .
/// Zero at E = 0, monotone increasing, and proportional to sqrt(E).
Real schottky_barrier_lowering(Real E);

/// The field at which the lowering equals the barrier, G(E*) = dG [V/m]:
///     E* = 4 pi eps0 dG^2 / e^3 .
/// Above it the barrier is gone and the rate expression has left the regime it
/// was derived for -- which is why it is reported, not clipped.
Real barrier_free_field(Real dG);

/// The quoted Iribarne-Thomson rate, as a function of its arguments alone:
///     j = (kT/h) * sigma_s * exp[ -(dG - G(E)) / (kT) ],   sigma_s = eps0 E .
/// [A/m^2].  Returns NaN for non-positive T or a negative field magnitude.
///
/// NOT VALIDATED: see the header.  Nothing in this project may present its
/// output as a predicted current.
Real iribarne_thomson_rate(Real E, Real dG, Real T);

/// The dimensionless form, so that a sensitivity study need not quote any
/// unsourced number: with x = E / E* and b = dG / (kT),
///     j / j0 = exp[ -b (1 - sqrt(x)) ],     j0 = (kT/h) eps0 E* .
/// This is the ONLY form any figure of this phase may use.
Real iribarne_thomson_dimensionless(Real x, Real b);

// ---------------------------------------------------------------------------
// The shipped data: a species table with the gaps left open
// ---------------------------------------------------------------------------
//
// The masses ARE sourced -- they follow from the molar masses of the ions,
// which the IoLiTec data sheet and the ILThermo component record give for the
// compound.  The BARRIERS are not, and are left at zero with an empty source.
// That is the blocker, in the data rather than in a comment.

/// EMI+ and BF4-, with masses from the compound's molar mass and the barriers
/// deliberately absent.
const EmissionModel& emibf4_emission_blocked(Polarity p);

/// A model whose parameters are all present but whose equation is still not
/// validated.  It exists ONLY so that the tests can reach the
/// EquationNotValidated branch and show that a complete parameter set is still
/// not enough.  It is not a substance.
///
/// IT OWNS ITS SPECIES.  A first version returned an EmissionModel pointing at
/// a shared static buffer, so building a second model silently rewrote the
/// first one's species -- and a test that compared the two polarities compared
/// one model with itself.  The test caught it because it printed both currents
/// and one of them was zero.  Ownership is the fix, not a comment.
struct SyntheticEmissionModel {
  EmittedSpecies species;
  EmissionModel model;
  const EmissionModel& operator*() const { return model; }
};
SyntheticEmissionModel synthetic_complete_model(Polarity p, Real dG, Real T);

}  // namespace es
