#pragma once
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "es/axisym_fem.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P9 -- axisymmetric versus 3D, and the contract for comparing with data
// ===========================================================================
//
// WHAT THIS FILE IS FOR.  Three things that have to exist before a 3D solver or
// a comparison with measurements can be honest:
//
//   1. A geometry kind that TRAVELS WITH a result, so that an axisymmetric
//      number cannot be reported as a three-dimensional one.  That is not
//      pedantry: an axisymmetric emitter and a real one differ by exactly the
//      thing an axisymmetric code cannot see, and the failure mode is a number
//      that looks right and answers a different question.
//   2. A REVOLUTION REFERENCE: the axisymmetric result, evaluated as a genuine
//      three-dimensional field and integrated by explicit 3D quadrature.  It is
//      not a 3D solver -- it solves nothing new -- but it is a real check that
//      the 2 pi r weighting of every axisymmetric integral in this project is
//      the three-dimensional integral it claims to be.
//   3. An IMPORT CONTRACT for measured data: value, unit, uncertainty with its
//      type, provenance, and the geometry the measurement belongs to.  Any of
//      them missing fails closed.
//
// WHAT IT IS NOT.  There is no 3D mesh, no 3D solver and no 3D result anywhere
// in this project.  ThreeDimensional exists as a LABEL that nothing can
// currently produce, which is the point: the test shows that the promotion path
// is closed.

// ---------------------------------------------------------------------------

enum class GeometryKind {
  /// The meridian half-plane with rotational symmetry.  Everything this project
  /// computes is of this kind.
  Axisymmetric = 0,
  /// A genuinely three-dimensional geometry.  NOTHING here produces one.
  ThreeDimensional,
  /// An axisymmetric result evaluated as a 3D field by revolution.  It is
  /// STILL axisymmetric physics; the label exists so that a revolved reference
  /// cannot be passed off as a 3D computation.
  RevolvedAxisymmetric,
};
const char* to_string(GeometryKind k);
const char* explain(GeometryKind k);
/// True only for a genuinely three-dimensional computation.
inline bool is_three_dimensional(GeometryKind k) { return k == GeometryKind::ThreeDimensional; }

/// A number with the geometry it came from attached.  The constructor is the
/// only way to set the kind, and there is no setter: a result cannot be
/// relabelled after the fact.
class LabelledResult {
 public:
  LabelledResult(Real value, GeometryKind kind, std::string quantity, std::string unit)
      : value_(value), kind_(kind), quantity_(std::move(quantity)), unit_(std::move(unit)) {}

  Real value() const { return value_; }
  GeometryKind kind() const { return kind_; }
  const std::string& quantity() const { return quantity_; }
  const std::string& unit() const { return unit_; }

  /// Returns the value ONLY if the result really is three-dimensional.
  /// Throws otherwise, naming the quantity and its actual kind.  This is the
  /// single gate through which a 3D claim has to pass.
  Real value_as_three_dimensional() const;

  void print(std::FILE* out) const;

 private:
  Real value_;
  GeometryKind kind_;
  std::string quantity_, unit_;
};

// ---------------------------------------------------------------------------
// The revolution reference
// ---------------------------------------------------------------------------
//
// An axisymmetric field phi(r, z) IS a three-dimensional field
// phi3(x, y, z) = phi(sqrt(x^2+y^2), z), and an axisymmetric vector field
// (E_r, E_z) IS the 3D field (E_r cos t, E_r sin t, E_z) with t the azimuth.
// Both statements are trivially true and both are worth checking in code,
// because the place they are usually got wrong is the INTEGRAL.

struct RevolvedField {
  const QuadMesh* mesh{nullptr};
  const std::vector<Real>* phi{nullptr};
  const std::vector<Real>* eps_r{nullptr};
  const std::vector<char>* active{nullptr};

  /// Potential at a genuine 3D point.
  Real potential(Real x, Real y, Real z) const;
  /// Field at a genuine 3D point, as three Cartesian components.
  void field(Real x, Real y, Real z, Real* Ex, Real* Ey, Real* Ez) const;
  /// Largest deviation from rotational invariance over `n_azimuth` azimuths at
  /// a given (r, z).  It must be at round-off: the field IS invariant by
  /// construction, and measuring it checks the evaluation, not the physics.
  Real azimuthal_variation(Real r, Real z, int n_azimuth) const;
};

/// Integral of `f(r, z)` over a surface of revolution, computed in TWO ways:
///   * the axisymmetric form  int f 2 pi r ds  along the meridian;
///   * an explicit 3D quadrature  int int f r dt ds  over `n_azimuth` azimuths.
/// The two must agree to the quadrature order.  This is the check that the
/// 2 pi r weighting used everywhere in this project really is the 3D integral.
struct RevolutionCheck {
  Real axisymmetric{0};
  Real three_dimensional{0};
  Real relative_difference{0};
  int n_azimuth{0}, n_meridian{0};
};

/// `meridian` gives the (r, z) points of the curve, in order; `f` is evaluated
/// at each of them.
RevolutionCheck revolution_surface_integral(const std::vector<Vec2>& meridian,
                                            const std::function<Real(Vec2)>& f,
                                            int n_azimuth);

/// The same for a VOLUME of revolution bounded by the meridian polyline and the
/// plane z = z_base.
RevolutionCheck revolution_volume(const std::vector<Vec2>& meridian, Real z_base,
                                  int n_azimuth);

// ---------------------------------------------------------------------------
// The import contract for measured data
// ---------------------------------------------------------------------------
//
// A measurement without a unit is a number; without a provenance it is a
// rumour; and without the geometry it belongs to it cannot be compared with
// anything.  Those three are HARD requirements and a set that breaks one of
// them is rejected whole.
//
// THE UNCERTAINTY IS DIFFERENT, and an earlier version of this contract got it
// wrong by treating it the same way.  A publication that reports a current of
// 210 nA without an error bar has not produced a broken record; it has produced
// an incomplete one.  Refusing to import it throws away a real measurement.
// Importing it as if it had an uncertainty would be worse.  So it is imported
// with the explicit state NotReported, it is archived, it may be DRAWN with its
// state visible -- and usable_quantitatively() is false for it, so no deviation,
// no chi-square and no pass/fail can be computed from it.

enum class UncertaintyType {
  /// Evaluated by statistical analysis of repeated observations (GUM type A).
  TypeA = 0,
  /// Evaluated by other means -- instrument specification, calibration
  /// certificate, judgement (GUM type B).
  TypeB,
  /// THE PUBLICATION DOES NOT STATE ONE.  This is a fact about the source, not
  /// a defect of the import, and the two must not be confused: a missing unit
  /// is a broken record, a missing uncertainty is an incomplete publication.
  /// Such a point may be archived and shown QUALITATIVELY; it may never carry a
  /// quantitative validation.  An earlier version of this file called this
  /// NotStated and made it a hard import failure, which threw away real data.
  NotReported,
};
const char* to_string(UncertaintyType u);
const char* explain(UncertaintyType u);

enum class ImportStatus {
  /// Complete: unit, uncertainty, provenance and geometry all present.  Only
  /// this status may back a quantitative comparison.
  Ok = 0,
  /// Everything present EXCEPT the uncertainty, which the publication does not
  /// state.  The point is imported, archived and may be drawn -- with its
  /// status visible -- but no quantitative claim may rest on it.
  OkUncertaintyNotReported,
  /// The hard errors.  Each of them means the record cannot be interpreted at
  /// all, and each rejects the whole set.
  MissingUnit,
  MissingProvenance,
  MissingGeometryKind,
  UnitMismatch,
};
const char* to_string(ImportStatus s);
const char* explain(ImportStatus s);
/// True if the point may be kept and displayed at all.
inline bool is_usable(ImportStatus s) {
  return s == ImportStatus::Ok || s == ImportStatus::OkUncertaintyNotReported;
}
/// True ONLY where a quantitative comparison is permitted.  Every place that
/// computes a deviation, a chi-square or a pass/fail must ask this one.
inline bool usable_quantitatively(ImportStatus s) { return s == ImportStatus::Ok; }
/// True for the failures that reject the whole set.
inline bool is_hard_error(ImportStatus s) { return !is_usable(s); }

struct MeasuredPoint {
  std::string quantity;
  Real value{0};
  std::string unit;              ///< SI, spelled out; empty = missing
  Real uncertainty{0};           ///< in the same unit; <= 0 = not reported
  UncertaintyType uncertainty_type{UncertaintyType::NotReported};
  Real coverage_factor{0};       ///< k; 0 = not stated
  std::string provenance;        ///< where the number can be found again
  GeometryKind geometry{GeometryKind::Axisymmetric};
  bool geometry_stated{false};
  std::string conditions;        ///< temperature, polarity, flow rate, ...

  ImportStatus check() const;
  void print(std::FILE* out) const;
};

/// Import a set of measured points.  Returns the first failing status and the
/// index of the point that failed, or Ok.  NOTHING is imported partially: a set
/// with one bad point is rejected whole, because a comparison that silently
/// dropped a point would be a comparison with a different data set.
struct ImportResult {
  /// The worst status in the set: a hard error if any point has one, otherwise
  /// OkUncertaintyNotReported if any point lacks an uncertainty, otherwise Ok.
  ImportStatus status{ImportStatus::Ok};
  Index first_bad{-1};                 ///< index of the first HARD failure
  std::string message;
  /// Empty on a hard error.  Otherwise every point, each carrying its own
  /// status -- the set is not homogeneous and pretending it is would be the
  /// same mistake in a new place.
  std::vector<MeasuredPoint> points;
  Index n_quantitative{0};             ///< points that may back a number
  Index n_qualitative_only{0};         ///< points archived without uncertainty
};
ImportResult import_measurements(const std::vector<MeasuredPoint>& raw);

// ---------------------------------------------------------------------------
// The validation matrix
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// SIX QUESTIONS, NOT ONE
// ---------------------------------------------------------------------------
//
// An earlier version of this matrix carried a single Comparability per row and
// the figure coloured the row by it.  A quantity that is comparable IN
// PRINCIPLE but blocked, or comparable and not converged, therefore appeared
// GREEN -- and green reads as success.  The total current is the clearest case:
// it is directly comparable between an axisymmetric model and a real device,
// and this project cannot compute it at all.
//
// Comparability and validation are different questions and are now answered
// separately.  A row is a set of independent verdicts:
//
//   1. comparable_geometry  -- can the quantity be compared between an
//                              axisymmetric computation and a 3D device at all?
//   2. implemented          -- does something in this project compute it?
//   3. converged            -- is the numerical result converged, by a stated
//                              criterion that was fixed BEFORE the measurement?
//   4. comparable_with_data -- could it be compared with a measurement, given
//                              an import that satisfies the contract?
//   5. validated            -- has it ACTUALLY been compared with measured data
//                              and agreed within the stated uncertainties?
//   6. blocked + reason     -- is it blocked, and by what?
//
// AND THE INVARIANT, checked in code rather than promised in prose:
// validated == Yes requires implemented == Yes AND converged == Yes AND
// comparable_with_data == Yes AND not blocked.  A quantity can never be shown
// as validated because it is merely comparable.

/// A per-axis answer.  There is no single overall verdict, on purpose.
enum class Assessment {
  No = 0,          ///< the axis is not satisfied
  Partial,         ///< satisfied only under a restriction that is stated
  Yes,             ///< satisfied
  NotApplicable,   ///< the axis does not apply to this quantity
};
const char* to_string(Assessment a);
const char* symbol(Assessment a);

/// Whether a quantity can be compared between an axisymmetric computation and a
/// three-dimensional device at all -- and if so, under what condition.
enum class Comparability {
  /// The same number in both, by symmetry or by construction.
  Direct = 0,
  /// Comparable only after an explicit reduction is stated (an average over the
  /// azimuth, a total instead of a local value, ...).
  AfterStatedReduction,
  /// Not comparable: the quantity does not exist in one of the two, or the
  /// axisymmetric model cannot represent what sets it.
  NotComparable,
};
const char* to_string(Comparability c);
/// The comparability axis as an Assessment, so that the figure never has to
/// map a Comparability onto a colour that could read as an overall verdict.
Assessment as_assessment(Comparability c);

struct ValidationEntry {
  const char* quantity{""};
  const char* unit{""};

  // --- axis 1: geometric / physical comparability -------------------------
  Comparability comparability{Comparability::NotComparable};
  const char* condition{""};       ///< what has to hold, or why it does not

  // --- axes 2 to 5 --------------------------------------------------------
  Assessment implemented{Assessment::No};
  Assessment converged{Assessment::No};
  Assessment comparable_with_data{Assessment::No};
  Assessment validated{Assessment::No};

  // --- axis 6 -------------------------------------------------------------
  bool blocked{false};
  const char* blocked_reason{""};   ///< non-empty exactly when blocked

  // --- the evidence behind axes 2 to 5 ------------------------------------
  const char* computed_by{""};      ///< which phase computes it here, or "-"
  const char* phase_status{""};     ///< that phase's own status word
  const char* convergence_note{""}; ///< what is known about convergence
  const char* validation_note{""};  ///< why it is not validated
  bool measurable{false};           ///< can it be measured on a bench at all
};

/// The invariant above, evaluated for one row.  Returns an empty string if the
/// row is consistent, and otherwise names what is wrong with it.
std::string inconsistency(const ValidationEntry& e);

/// The matrix.  It is data, not code: a reader should be able to check it
/// against the phase documents line by line.
std::vector<ValidationEntry> validation_matrix();

/// How many rows reach each verdict on the validation axis.  Used by the run
/// and by the test so that "nothing is validated" is a measured statement.
struct ValidationTally {
  Index n_rows{0};
  Index n_comparable{0};       ///< Yes or Partial on axis 1
  Index n_implemented{0};
  Index n_converged{0};
  Index n_comparable_with_data{0};
  Index n_validated{0};
  Index n_blocked{0};
};
ValidationTally tally(const std::vector<ValidationEntry>& m);

}  // namespace es
