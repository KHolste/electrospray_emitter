#pragma once
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/emission.hpp"
#include "es/fluid.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// Beam transport
// ===========================================================================
//
// Particles are tracked in the meridian half-plane.  In axisymmetry a
// macroparticle is not a point but a RING, which means its own field is exactly
// the axisymmetric Green's function the BEM already uses -- space charge needs
// no second discretisation and no artificial field boundary.
//
// Self-consistency loop, when space charge is enabled:
//   1. trace rays in the Laplace field of the electrodes
//   2. convert each trajectory into ring macroparticles, Q_i = I_ray dl_i / v_i
//      (the steady-state charge held in that piece of the beam)
//   3. sample the ring potential at the BEM collocation points and re-solve, so
//      the electrodes carry the charge induced by the beam
//   4. re-trace in (induced BEM field + direct ring field), under-relaxed
//
// Space charge is off by default: at the 10-100 nA of a single emitter it is a
// small correction, and each iteration costs roughly as much as the whole
// Laplace trace.  It matters for emitter arrays and for the near-apex region.

/// One population launched from the surface, as a fraction of the local current.
struct BeamSpecies {
  std::string name{"ion"};
  Real qm{5.0e5};      ///< charge-to-mass ratio [C/kg]
  Real fraction{1.0};  ///< share of the emitted CURRENT carried by this species
};

enum class RayStatus : int {
  Flying = 0,       ///< still in flight when the step budget ran out
  Transmitted = 1,  ///< reached the downstream plane
  Intercepted = 2,  ///< hit an electrode
  Escaped = 3,      ///< left the radial domain
};

const char* ray_status_name(RayStatus s);

struct Ray {
  Vec2 x0{};               ///< launch point
  Vec2 x{};                ///< current / final position
  Vec2 v{};                ///< current / final velocity [m/s]
  Real qm{0};              ///< [C/kg]
  Real current{0};         ///< current carried by this ray [A]
  std::string species{};
  RayStatus status{RayStatus::Flying};
  Tag hit_tag{Tag::Other};  ///< which surface it struck, if intercepted
  Real angle{0};            ///< final half-angle from the axis [rad]
  Real energy_eV{0};        ///< final kinetic energy per charge [eV]
  int steps{0};
  std::vector<Vec2> path;   ///< sampled trajectory (for plotting / space charge)
  std::vector<Real> speed;  ///< |v| at each sampled point
};

struct BeamParams {
  Real z_end{2.0e-3};       ///< downstream plane where tracing stops [m]
  Real r_max{1.0e-2};       ///< radial cut-off [m]
  int max_steps{20000};
  Real cfl{0.05};           ///< step limiter: fraction of the local scale per step
  Real launch_offset{0.05}; ///< launch height above the surface, in element lengths
  int path_samples{80};     ///< trajectory points stored per ray
  int space_charge_iters{0};
  Real space_charge_relax{0.5};
  bool include_wetted_metal{false};
  bool verbose{false};
};

struct BeamResult {
  std::vector<Ray> rays;
  Real current_launched{0};
  Real current_transmitted{0};
  Real current_intercepted{0};
  Real interception_fraction{0};
  Real half_angle_50{0};   ///< half-angle containing 50% of the transmitted current
  Real half_angle_95{0};   ///< ... 95%
  Real mean_energy_eV{0};  ///< current-weighted
  int space_charge_iterations{0};

  void write_rays_csv(const std::string& path) const;
  void write_paths_csv(const std::string& path) const;
  void print(std::FILE* out) const;
};

/// Launch weights taken from the Iribarne-Thomson rate on the solved surface.
/// `bem` must already carry a solution.
BeamResult trace_beam(BemSolver& bem, const Fluid& f, Real T,
                      const std::vector<BeamSpecies>& species, const BeamParams& p);

/// Launch with an explicitly prescribed current per surface element (same order
/// as bem.mesh().elems; zero entries are skipped).  Use this to drive the trace
/// from a cone-jet model, from measured data, or from any custom distribution.
BeamResult trace_beam_with_weights(BemSolver& bem, const std::vector<Real>& element_current,
                                   const std::vector<BeamSpecies>& species, const BeamParams& p);

/// Potential and field of a single charged ring of total charge Q at xp.
/// V = Q K(m) / (2 pi^2 eps0 S) -- the same kernel the BEM uses, rescaled from
/// "per unit surface density" to "per unit total ring charge".
Real ring_potential(Real Q, Vec2 x, Vec2 xp);
Vec2 ring_field(Real Q, Vec2 x, Vec2 xp);

}  // namespace es
