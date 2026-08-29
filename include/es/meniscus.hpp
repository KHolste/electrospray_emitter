#pragma once
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/fluid.hpp"
#include "es/geometry.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// Static electrified meniscus
// ===========================================================================
//
// The free surface of a perfectly conducting liquid pinned at a contact line of
// radius r_c satisfies the Young-Laplace equation with a Maxwell traction:
//
//     gamma ( dphi/ds + sin(phi)/r ) = dp - rho g z + (eps0/2) E_n^2
//
// with s the arclength from the apex, phi the angle of the surface tangent
// below the horizontal, and the two terms on the left the meridional and
// azimuthal curvatures.  The electric term always pulls outward, so raising the
// voltage inflates the meniscus.
//
// WHAT IS AN INDEPENDENT VARIABLE, AND WHY IT MATTERS NUMERICALLY
//
// dp (the feed pressure) and U (the applied voltage) are the physical inputs;
// the apex height h comes out.  But h(U) has a vertical tangent at the onset of
// emission -- that IS the onset: a saddle-node bifurcation beyond which no
// static meniscus exists and the liquid erupts into a jet.  Marching in U can
// therefore never reach it; the iteration simply stops converging somewhere
// short of the fold and the result looks like a numerical failure rather than
// the physics it is.
//
// So the solver inverts the roles: h is prescribed and U is solved for.  That
// parameterisation is regular through the fold, the whole branch including its
// turning point is traced, and the onset voltage is read off as the maximum of
// U(h).  Everything past that maximum is the unstable branch -- real, but not
// realisable at fixed voltage.
//
// ASSUMPTIONS
//  * static: no flow.  Valid near onset and at the low flow rates of the pure
//    ionic regime.  At cone-jet flow rates the viscous pressure drop along the
//    cone is not negligible and this model overestimates the apex height.
//  * perfectly conducting liquid.  Justified by the charge relaxation time,
//    ~1e-10 s for ionic liquids -- see Fluid::charge_relaxation_time().
//  * pinned contact line at r_c.  Correct for a capillary with a sharp edge.
//    For externally wetted or porous emitters the contact line is not pinned;
//    treating the tip as pinned at its radius is a first-order stand-in only.
//  * the free surface is monotone in r.  Overhanging menisci -- a pendant drop
//    bulged past a hemisphere -- are legitimate Young-Laplace solutions but are
//    not on the path to a Taylor cone, and the solver reports them as
//    unreachable rather than tracking them.  In practice this bounds the
//    reachable apex height at roughly r_c in the field-free limit; with a field
//    the cone-like branch goes far higher without ever turning vertical.

struct MeniscusParams {
  Real r_contact{1.0e-5};  ///< pinning radius [m]
  Real z_contact{0.0};     ///< axial position of the contact line [m]
  Real delta_p{0.0};       ///< feed pressure at the apex, above ambient [Pa]
  Real gamma{0.0452};      ///< surface tension [N/m]
  Real rho{0.0};           ///< density, for the hydrostatic term; 0 disables it

  int n_nodes{81};             ///< nodes on the free surface
  Real apex_clustering{1.8};   ///< >1 clusters nodes toward the apex
  Real h_far{0};               ///< element size on the electrodes (0 = builder default)

  int max_outer{40};       ///< shape <-> field iterations
  Real relax{0.6};         ///< under-relaxation of the shape update
  Real tol{2.0e-4};        ///< convergence: max node motion / r_contact
  bool verbose{false};
};

struct MeniscusShape {
  std::vector<Vec2> nodes;  ///< contact line -> apex (BEM traversal order)
  Real height{0};           ///< apex above the contact line [m]
  Real arclength{0};        ///< meridian length of the free surface [m]
  Real apex_radius{0};      ///< radius of curvature at the apex [m]
  Real half_angle{0};       ///< local cone half-angle at mid-arc [rad]

  /// Free-surface mesh in BEM orientation (outward normal into the vacuum).
  Mesh to_mesh(Real potential) const;
};

/// Initial guess: an ellipse arc of height h pinned at r_c.  Reduces to a
/// spherical cap when h == r_c.
MeniscusShape initial_shape(Real r_c, Real z_c, Real h, int n_nodes, Real clustering);

struct MeniscusSolution {
  MeniscusShape shape;
  Real voltage{0};        ///< emitter-to-extractor voltage sustaining this shape [V]
  Real apex_field{0};     ///< |E_n| at the apex [V/m]
  Real peak_field{0};     ///< max |E_n| on the free surface [V/m]
  Real delta_p{0};
  bool converged{false};
  /// Set when the feed pressure alone already lifts the meniscus past the
  /// requested height, so no voltage was needed.  The achieved height is then
  /// shape.height, which will not equal the h that was asked for.
  bool voltage_clamped{false};
  int iterations{0};
  Real residual{0};       ///< final max node motion / r_contact
};

/// Couples the Young-Laplace shape solver to the BEM field solver.
class MeniscusSolver {
 public:
  /// `electrodes` must contain everything except the free surface: the emitter
  /// body (see make_capillary_open) and any extractor.  The emitter is held at
  /// U and the extractor at 0, so U is the emitter-to-extractor voltage.
  MeniscusSolver(Mesh electrodes, MeniscusParams params);

  /// Solve for the voltage that sustains an apex height h.  Warm-starts from
  /// the shape left by the previous call, which is what makes continuation
  /// cheap; pass a shape explicitly to override.
  MeniscusSolution solve_at_height(Real h, const MeniscusShape* start = nullptr);

  /// Solve for the STABLE meniscus that a given voltage sustains.  Runs a
  /// coarse continuation to locate the fold, then bisects on apex height along
  /// the stable side of the branch.  Returns converged == false if U exceeds
  /// the onset voltage, because then no static meniscus exists at all.
  MeniscusSolution solve_at_voltage(Real U, Real h_max, int scout_steps = 14);

  /// Trace the branch h = h_min ... h_max.  The onset voltage is the maximum of
  /// `voltage` over the returned branch; everything beyond it is unstable.
  std::vector<MeniscusSolution> continuation(Real h_min, Real h_max, int n_steps);

  /// Convenience: run a continuation and return the peak of U(h), refined by a
  /// parabolic fit through the three points around the maximum.
  struct Onset {
    Real voltage{0};
    Real height{0};
    Real apex_field{0};
    Real apex_radius{0};
    Real half_angle{0};
    bool found{false};
  };
  static Onset find_onset(const std::vector<MeniscusSolution>& branch);

  /// The solver keeps the last assembled system, so callers can post-process
  /// (ion emission, beam launch) without rebuilding it.
  const BemSolver& bem() const { return bem_; }
  BemSolver& bem() { return bem_; }

  const MeniscusParams& params() const { return params_; }

  static void write_branch_csv(const std::string& path,
                               const std::vector<MeniscusSolution>& branch);

 private:
  Mesh electrodes_;
  MeniscusParams params_;
  BemSolver bem_;
  MeniscusShape last_;
  bool have_last_{false};
};

}  // namespace es
