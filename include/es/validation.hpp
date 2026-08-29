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
// A measurement without a unit is a number; without an uncertainty it is an
// anecdote; without a provenance it is a rumour; and without the geometry it
// belongs to it cannot be compared with anything.  All four are required.

enum class UncertaintyType {
  /// Evaluated by statistical analysis of repeated observations (GUM type A).
  TypeA = 0,
  /// Evaluated by other means -- instrument specification, calibration
  /// certificate, judgement (GUM type B).
  TypeB,
  /// Not stated by the source.  A dataset with this fails the import.
  NotStated,
};
const char* to_string(UncertaintyType u);

enum class ImportStatus {
  Ok = 0,
  MissingUnit,
  MissingUncertainty,
  MissingProvenance,
  MissingGeometryKind,
  UnitMismatch,
};
const char* to_string(ImportStatus s);
const char* explain(ImportStatus s);
inline bool is_usable(ImportStatus s) { return s == ImportStatus::Ok; }

struct MeasuredPoint {
  std::string quantity;
  Real value{0};
  std::string unit;              ///< SI, spelled out; empty = missing
  Real uncertainty{0};           ///< in the same unit; <= 0 = missing
  UncertaintyType uncertainty_type{UncertaintyType::NotStated};
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
  ImportStatus status{ImportStatus::Ok};
  Index first_bad{-1};
  std::string message;
  std::vector<MeasuredPoint> points;   ///< empty unless status is Ok
};
ImportResult import_measurements(const std::vector<MeasuredPoint>& raw);

// ---------------------------------------------------------------------------
// The validation matrix
// ---------------------------------------------------------------------------

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

struct ValidationEntry {
  const char* quantity{""};
  const char* unit{""};
  Comparability comparability{Comparability::NotComparable};
  const char* condition{""};       ///< what has to hold, or why it does not
  const char* computed_by{""};     ///< which phase computes it here, or "-"
  const char* status{""};          ///< the phase status
  bool measurable{false};          ///< can it be measured on a bench at all
};

/// The matrix.  It is data, not code: a reader should be able to check it
/// against the phase documents line by line.
std::vector<ValidationEntry> validation_matrix();

}  // namespace es
