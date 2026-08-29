#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/axisym_fem.hpp"
#include "es/space_charge.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P7 -- electrostatic particle transport on the field interface of P6
// ===========================================================================
//
// WHAT THIS IS.  Newton's equation for an explicit species (q, m) in the field
// interpolated from the axisymmetric volume mesh, with the surface a particle
// hits classified, and with a current and particle balance that closes exactly.
//
// WHAT IT IS NOT, and this decides what the numbers mean:
//
//   * THERE IS NO PHYSICAL SOURCE.  P5 is blocked, so particles are launched
//     from an explicitly PRESCRIBED distribution.  The result is therefore a
//     TRANSPORT RESPONSE -- what the device does to a given beam -- and NOT a
//     current prediction.  Every output says so.
//   * no droplet-beam weighting.  The prototype's droplet branch stays
//     disabled and nothing here reactivates it.
//   * no magnetic field, no collisions, no secondary emission, no sputtering.
//   * no self-consistent space-charge loop.  The field may INCLUDE a prescribed
//     space charge (P6), but the particles do not update it while they fly.
//     That is a deliberate separation: it makes the two effects comparable
//     (with and without prescribed space charge) instead of entangled.
//
// THE SPECIES IS EXPLICIT.  Charge and mass are given by the caller with their
// signs; there is no default species, no "the ion", and no magnitude-only path.
// Reversing the polarity means giving a species with the opposite charge sign
// AND reversing the applied potentials -- the two together leave the trajectory
// invariant, and that is a test.

// ---------------------------------------------------------------------------

struct TransportSpecies {
  const char* name{""};
  Real charge{0};   ///< [C], SIGNED
  Real mass{0};     ///< [kg]
  bool complete() const { return charge != 0.0 && mass > 0.0; }
  Real charge_to_mass() const { return charge / mass; }
};

/// What a particle ran into.  `Flying` means it was still inside when the step
/// budget ran out, which is a result and not a failure -- but it is counted
/// separately, because a balance that folded it into "lost" would hide it.
enum class Fate {
  Flying = 0,
  HitEmitter,
  HitPolymer,
  HitExtractor,
  LeftThroughAperture,
  LeftDomain,
};
const char* to_string(Fate f);
/// True for the fates that mean the particle left usefully.
inline bool is_extracted(Fate f) { return f == Fate::LeftThroughAperture; }
/// True for the fates that mean the particle was intercepted by a surface.
inline bool is_intercepted(Fate f) {
  return f == Fate::HitEmitter || f == Fate::HitPolymer || f == Fate::HitExtractor;
}

// ---------------------------------------------------------------------------
// The geometry a particle can hit
// ---------------------------------------------------------------------------
//
// A rectangular meridian domain r in [0, R], z in [0, Z] with named surfaces.
// It is deliberately simple: the point of P7 is the transport, the hit
// classification and the balance, not a new device geometry.  The classifier
// takes a position and says which surface it crossed, and the four kinds are
// all reachable in this one geometry, so the balance can actually be exercised.
struct TransportDomain {
  Real R{0}, Z{0};
  Real aperture_radius{0};   ///< r <= this at z = Z is open; beyond it is the extractor
  /// Classify a position that has left the domain.  Returns Flying if it has
  /// not.  The axis is a symmetry line, not a wall: r < 0 is reflected by the
  /// caller, never counted as a hit.
  Fate classify(Vec2 x) const;
  bool inside(Vec2 x) const;
};

// ---------------------------------------------------------------------------

struct TracedParticle {
  Vec2 x, v;
  /// Where it was launched.  Kept separately because `x` is advanced in place
  /// and a report that printed `x` as the start would print the end.
  Vec2 x_start;
  /// The current this macroparticle stands for [A].  For a steady beam this is
  /// the natural bookkeeping quantity: the balance below sums it, and the sum
  /// is exactly the launched total because nothing is created or destroyed.
  Real current{0};
  Fate fate{Fate::Flying};
  Real time{0};
  Real path_length{0};
  Vec2 x_end;                ///< where it ended
  Real energy_gain{0};       ///< [J], measured as the kinetic energy change
  Real potential_start{0}, potential_end{0};   ///< [V]
  int steps{0};
};

struct TransportBalance {
  Real launched{0}, extracted{0}, intercepted{0}, still_flying{0};   ///< [A]
  Real hit_emitter{0}, hit_polymer{0}, hit_extractor{0};             ///< [A]
  Index n_launched{0}, n_extracted{0}, n_intercepted{0}, n_flying{0};
  /// |launched - (extracted + intercepted + still_flying)| / launched.
  /// It must be at machine level: every particle has exactly one fate.
  Real closure_error{0};
  /// Largest relative violation of  1/2 m v^2 - q phi = const  along a
  /// trajectory.  It is the integrator's energy error, not a physical loss.
  Real max_energy_error{0};
  void print(std::FILE* out) const;
};

struct TransportResult {
  std::vector<TracedParticle> particles;
  TransportBalance balance;
  Real dt{0};
  int max_steps{0};
  std::string message;
};

/// Push the given particles through the field of `phi` on `mesh` until they
/// leave the domain or the step budget runs out.
///
/// Integrator: velocity Verlet, which for a force that depends on position only
/// is time-reversible and conserves energy to second order without drift.  The
/// field is the interpolated recovered field of P6, so deposition and
/// interpolation use the same basis.
TransportResult transport_particles(const QuadMesh& mesh, const std::vector<Real>& phi,
                                    const std::vector<Real>& eps_r,
                                    const std::vector<char>& active,
                                    const TransportDomain& domain,
                                    const TransportSpecies& species,
                                    const std::vector<TracedParticle>& launch, Real dt,
                                    int max_steps);

}  // namespace es
