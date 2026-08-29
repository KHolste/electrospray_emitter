#pragma once
#include <stdexcept>
#include <string>

namespace es {

// ---------------------------------------------------------------------------
// Result status
// ---------------------------------------------------------------------------
//
// A single `converged` flag cannot distinguish "the iteration reached its
// tolerance" from "the iteration reached its tolerance but at a different
// operating point than the one that was requested".  The prototype conflated
// the two and reported success for a shape belonging to 872 V when 500 V had
// been asked for.  Every solver entry point therefore returns an explicit
// status, and callers must not treat any value other than Converged as usable.

enum class SolveStatus {
  Converged = 0,          ///< converged AND, where a target was given, it was met
  NotConverged,           ///< iteration did not reach its tolerance
  VoltageNotBracketed,    ///< requested voltage lies outside the reachable branch
  VoltageMismatch,        ///< search terminated away from the requested voltage
  ShapeIntegrationFailed, ///< the Young-Laplace march did not reach the contact line
  NoStaticFoldFound,      ///< branch has no interior turning point to search below
  AmbiguousBranch,        ///< several solutions exist; none may be chosen silently
  NotAttempted,           ///< default-constructed; nothing was solved
};

const char* to_string(SolveStatus s);

/// Human-readable cause, suitable for printing to a user.
const char* explain(SolveStatus s);

inline bool is_usable(SolveStatus s) { return s == SolveStatus::Converged; }

// ---------------------------------------------------------------------------

/// Why a branch does not have a usable static turning point.
enum class FoldStatus {
  Found = 0,
  TooFewPoints,        ///< fewer than three converged points
  Monotone,            ///< voltage rises (or falls) throughout: no turning point
  MaximumAtBoundary,   ///< largest value is the first or last point, so not interior
  NotAttempted,
};

const char* to_string(FoldStatus s);
const char* explain(FoldStatus s);

// ---------------------------------------------------------------------------
// Fail-closed marker for models that are not implemented yet
// ---------------------------------------------------------------------------
//
// Options whose physics is not implemented must not silently produce numbers
// that look like results.  They throw this, and the applications turn it into a
// clear message naming the phase in which the model is due.

class NotImplementedInThisPhase : public std::runtime_error {
 public:
  NotImplementedInThisPhase(const std::string& what_feature, const std::string& phase,
                            const std::string& reason)
      : std::runtime_error(what_feature + " ist in dieser Phase nicht implementiert und wurde "
                           "deshalb abgelehnt.\n  Grund : " + reason +
                           "\n  Geplant: " + phase +
                           "\n  Siehe docs/05_implementation_plan.md."),
        feature_(what_feature), phase_(phase), reason_(reason) {}

  const std::string& feature() const { return feature_; }
  const std::string& phase() const { return phase_; }
  const std::string& reason() const { return reason_; }

 private:
  std::string feature_, phase_, reason_;
};

}  // namespace es
