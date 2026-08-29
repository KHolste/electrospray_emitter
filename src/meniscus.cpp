#include "es/meniscus.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {
namespace {

using constants::eps0;
using constants::g0;
using constants::pi;

constexpr Real kInf = std::numeric_limits<Real>::infinity();

/// Normal field per volt on the free surface, as a function of radius.
/// Sampled from the BEM basis solution; the shape ODE reads it back while it
/// marches outward, so it must be queryable at any r in [0, r_c].
class FieldProfile {
 public:
  FieldProfile() = default;
  void set(std::vector<Real> r, std::vector<Real> e) {
    r_ = std::move(r);
    e_ = std::move(e);
    // Sort by radius; the mesh runs contact line -> apex, i.e. backwards.
    std::vector<std::size_t> idx(r_.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return r_[a] < r_[b]; });
    std::vector<Real> rr, ee;
    rr.reserve(r_.size());
    ee.reserve(e_.size());
    for (std::size_t i : idx) { rr.push_back(r_[i]); ee.push_back(e_[i]); }
    r_ = std::move(rr);
    e_ = std::move(ee);
  }
  bool empty() const { return r_.size() < 2; }

  Real at(Real r) const {
    if (r_.empty()) return 0.0;
    if (r <= r_.front()) return e_.front();
    if (r >= r_.back()) return e_.back();
    const auto it = std::lower_bound(r_.begin(), r_.end(), r);
    const std::size_t j = static_cast<std::size_t>(it - r_.begin());
    const Real t = (r - r_[j - 1]) / std::max(r_[j] - r_[j - 1], 1e-300);
    return e_[j - 1] + t * (e_[j] - e_[j - 1]);
  }

 private:
  std::vector<Real> r_, e_;
};

struct RawShape {
  std::vector<Vec2> pts;  ///< apex (0,0) outward, z measured DOWN from the apex
  std::vector<Real> s;    ///< arclength
  std::vector<Real> phi;  ///< tangent angle below horizontal
  Real height{0};
  Real arclength{0};
  Real apex_curvature{0};
  bool ok{false};
};

struct OdeCtx {
  Real gamma;
  Real delta_p;
  Real rho;
  Real U;
  const FieldProfile* prof;
};

inline Real dphi_ds(const OdeCtx& c, Real r, Real z, Real phi) {
  const Real E = c.U * c.prof->at(r);
  const Real p = c.delta_p - c.rho * g0 * z + 0.5 * eps0 * E * E;
  return p / c.gamma - std::sin(phi) / std::max(r, 1e-300);
}

/// March the Young-Laplace ODE from the apex until r reaches r_c.
/// State: r (outward), z (measured downward from the apex, so z <= 0 going
/// down is expressed as z increasing downward -- here z is the DROP, positive
/// downward, to keep signs readable), phi (tangent angle below horizontal).
RawShape integrate_shape(const OdeCtx& c, Real r_c) {
  RawShape out;
  // Regular series start at the apex, where sin(phi)/r is 0/0.  Both principal
  // curvatures are equal there, so each carries half the driving pressure.
  const Real E0 = c.U * c.prof->at(0.0);
  const Real kappa0 = 0.5 * (c.delta_p + 0.5 * eps0 * E0 * E0) / c.gamma;
  out.apex_curvature = 2.0 * kappa0;  // total curvature = 2 * kappa0

  Real s = 1e-7 * r_c;
  Real r = s;
  Real zdrop = 0.5 * kappa0 * s * s;
  Real phi = kappa0 * s;

  out.pts.push_back({0.0, 0.0});
  out.s.push_back(0.0);
  out.phi.push_back(0.0);

  const Real s_max = 200.0 * r_c;
  const Real ds_max = r_c / 40.0;

  auto deriv = [&](Real rr, Real zz, Real pp, Real& dr, Real& dz, Real& dp) {
    dr = std::cos(pp);
    dz = std::sin(pp);  // z is the drop below the apex, so it grows with sin(phi)
    dp = dphi_ds(c, rr, -zz, pp);  // hydrostatic term wants the signed height
  };

  int guard = 0;
  while (r < r_c && s < s_max && ++guard < 2000000) {
    // Step control: never turn more than ~1 degree, never move more than a few
    // percent of the local radius (the sin(phi)/r term stiffens near the axis).
    Real d1, d2, d3;
    deriv(r, zdrop, phi, d1, d2, d3);
    // ~0.3 deg of turning per step and 2% of the local radius.  The resampling
    // that follows interpolates linearly between raw points, so the chord error
    // ds^2/(8R) is what sets the accuracy of the returned shape; this keeps it
    // near 1e-6 of the meniscus size.
    const Real ds_curv = 0.005 / std::max(std::abs(d3), 1e-300);
    Real ds = std::min({ds_max, ds_curv, 0.02 * std::max(r, 1e-12) + 1e-9 * r_c});
    ds = std::max(ds, 1e-12 * r_c);

    // classic RK4
    Real k1r, k1z, k1p, k2r, k2z, k2p, k3r, k3z, k3p, k4r, k4z, k4p;
    deriv(r, zdrop, phi, k1r, k1z, k1p);
    deriv(r + 0.5 * ds * k1r, zdrop + 0.5 * ds * k1z, phi + 0.5 * ds * k1p, k2r, k2z, k2p);
    deriv(r + 0.5 * ds * k2r, zdrop + 0.5 * ds * k2z, phi + 0.5 * ds * k2p, k3r, k3z, k3p);
    deriv(r + ds * k3r, zdrop + ds * k3z, phi + ds * k3p, k4r, k4z, k4p);

    const Real rn = r + ds / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r);
    const Real zn = zdrop + ds / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z);
    const Real pn = phi + ds / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p);

    // The surface has turned vertical: dr/ds = cos(phi) has reached zero, so r
    // stops growing and the contact line is never reached on this branch.
    // Overhanging (past-hemisphere) menisci do exist as Young-Laplace solutions
    // but they are not the electrospray regime -- an electrified meniscus on
    // its way to a Taylor cone is monotone in r all the way out.  Report the
    // failure rather than stalling: the bisection reads it as "too much
    // voltage" and backs off, which is exactly the right response.
    if (pn >= 0.5 * pi) { out.ok = false; return out; }
    if (!std::isfinite(rn) || !std::isfinite(zn) || !std::isfinite(pn)) { out.ok = false; return out; }

    if (rn >= r_c) {
      // Land exactly on the contact line by linear interpolation in this step.
      const Real t = (r_c - r) / std::max(rn - r, 1e-300);
      out.pts.push_back({r_c, -(zdrop + t * (zn - zdrop))});
      out.s.push_back(s + t * ds);
      out.phi.push_back(phi + t * (pn - phi));
      out.height = zdrop + t * (zn - zdrop);
      out.arclength = s + t * ds;
      out.ok = true;
      return out;
    }

    r = rn; zdrop = zn; phi = pn; s += ds;
    out.pts.push_back({r, -zdrop});
    out.s.push_back(s);
    out.phi.push_back(phi);
  }
  out.ok = false;
  return out;
}

/// Resample a raw shape onto n nodes, clustered toward the apex, and return it
/// in BEM order (contact line first, apex last).
MeniscusShape resample(const RawShape& raw, Real z_contact, int n, Real clustering) {
  MeniscusShape m;
  m.height = raw.height;
  m.arclength = raw.arclength;
  m.apex_radius = (raw.apex_curvature > 0.0) ? 2.0 / raw.apex_curvature : kInf;

  // Cluster: s(t) = S * t^p, t in [0,1] measured from the apex.
  const Real p = std::max(1.0, clustering);
  std::vector<Vec2> apex_to_rim(static_cast<std::size_t>(n));
  std::size_t k = 0;
  for (int i = 0; i < n; ++i) {
    const Real t = static_cast<Real>(i) / (n - 1);
    const Real starget = raw.arclength * std::pow(t, p);
    while (k + 2 < raw.s.size() && raw.s[k + 1] < starget) ++k;
    const Real s0 = raw.s[k], s1 = raw.s[k + 1];
    const Real u = (s1 > s0) ? std::clamp((starget - s0) / (s1 - s0), 0.0, 1.0) : 0.0;
    apex_to_rim[static_cast<std::size_t>(i)] = raw.pts[k] + u * (raw.pts[k + 1] - raw.pts[k]);
  }
  apex_to_rim.front() = raw.pts.front();
  apex_to_rim.back() = raw.pts.back();

  // Local cone half-angle at mid-arc: 90 deg minus the tangent angle.  For a
  // fully developed Taylor cone this approaches 49.3 deg.
  {
    const Real starget = 0.5 * raw.arclength;
    std::size_t j = 0;
    while (j + 2 < raw.s.size() && raw.s[j + 1] < starget) ++j;
    m.half_angle = 0.5 * pi - raw.phi[j];
  }

  m.nodes.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const Vec2 q = apex_to_rim[static_cast<std::size_t>(n - 1 - i)];
    m.nodes[static_cast<std::size_t>(i)] = {q.r, q.z + z_contact + raw.height};
  }
  return m;
}

}  // namespace

// ---------------------------------------------------------------------------

Mesh MeniscusShape::to_mesh(Real potential) const {
  Mesh m;
  m.begin_body(Tag::FreeSurface, potential);
  for (const Vec2& p : nodes) m.add_node(p);
  m.end_body(false);
  return m;
}

MeniscusShape initial_shape(Real r_c, Real z_c, Real h, int n_nodes, Real clustering) {
  MeniscusShape m;
  m.height = h;
  m.nodes.resize(static_cast<std::size_t>(n_nodes));
  const Real p = std::max(1.0, clustering);
  // Ellipse arc: r = r_c sin(t), drop = h (1 - cos t).  Circle when h == r_c.
  for (int i = 0; i < n_nodes; ++i) {
    const Real tt = static_cast<Real>(n_nodes - 1 - i) / (n_nodes - 1);  // 1 at rim
    const Real t = 0.5 * pi * std::pow(tt, p);
    m.nodes[static_cast<std::size_t>(i)] = {r_c * std::sin(t), z_c + h * std::cos(t)};
  }
  m.nodes.front() = {r_c, z_c};
  m.nodes.back() = {0.0, z_c + h};
  m.apex_radius = (h > 0.0) ? r_c * r_c / h : kInf;
  m.arclength = 0.0;
  for (std::size_t i = 1; i < m.nodes.size(); ++i)
    m.arclength += norm(m.nodes[i] - m.nodes[i - 1]);
  return m;
}

// ---------------------------------------------------------------------------

MeniscusSolver::MeniscusSolver(Mesh electrodes, MeniscusParams params)
    : electrodes_(std::move(electrodes)), params_(params) {}

MeniscusSolution MeniscusSolver::solve_at_height(Real h, const MeniscusShape* start) {
  MeniscusSolution sol;
  sol.delta_p = params_.delta_p;

  MeniscusShape shape;
  if (start) shape = *start;
  else if (have_last_) shape = last_;
  else shape = initial_shape(params_.r_contact, params_.z_contact, h, params_.n_nodes,
                             params_.apex_clustering);

  // If the warm start has the wrong height, stretch it so the first BEM solve
  // already sees roughly the right geometry.
  if (shape.height > 0.0 && std::abs(shape.height - h) > 1e-12) {
    const Real f = h / shape.height;
    for (Vec2& p : shape.nodes) p.z = params_.z_contact + f * (p.z - params_.z_contact);
    shape.height = h;
  }

  FieldProfile prof;
  Real U = 0.0;
  Real relax = params_.relax;
  Real prev_disp = kInf;

  for (int it = 0; it < params_.max_outer; ++it) {
    // ---- field on the current shape ---------------------------------------
    Mesh full = merge({electrodes_, shape.to_mesh(0.0)});
    bem_.set_mesh(std::move(full));
    bem_.solve_basis();
    // Emitter at 1 V, extractor at 0 V.  E_n per volt = sigma_basis / eps0.
    {
      const std::vector<Real>& b = bem_.sigma_for({1.0, 0.0, 0.0});
      std::vector<Real> rr, ee;
      for (Index i = 0; i < bem_.size(); ++i) {
        const Element& el = bem_.mesh().elems[static_cast<std::size_t>(i)];
        if (el.tag != Tag::FreeSurface) continue;
        rr.push_back(el.mid.r);
        ee.push_back(std::abs(b[static_cast<std::size_t>(i)]) / eps0);
      }
      if (rr.size() < 2) throw std::runtime_error("MeniscusSolver: free surface has < 2 elements");
      prof.set(std::move(rr), std::move(ee));
    }

    // ---- voltage that reproduces the prescribed apex height ---------------
    OdeCtx ctx{params_.gamma, params_.delta_p, params_.rho, 0.0, &prof};
    auto height_at = [&](Real u) -> Real {
      ctx.U = u;
      const RawShape rs = integrate_shape(ctx, params_.r_contact);
      return rs.ok ? rs.height : kInf;
    };

    Real lo = 0.0, hi = 0.0;
    const Real h_lo = height_at(0.0);
    if (!(h_lo < h)) {
      // Feed pressure alone already exceeds the target height: no voltage
      // needed, and h is not reachable from above by raising U.
      U = 0.0;
    } else {
      // Bracket by doubling.  A field scale that lifts the meniscus by ~r_c is
      // the natural starting guess.
      hi = std::sqrt(4.0 * params_.gamma / (eps0 * params_.r_contact)) /
           std::max(prof.at(0.0), 1e-300);
      if (!(hi > 0.0) || !std::isfinite(hi)) hi = 1.0;
      int guard = 0;
      while (height_at(hi) < h && ++guard < 60) { lo = hi; hi *= 1.6; }
      for (int b = 0; b < 80; ++b) {
        const Real mid = 0.5 * (lo + hi);
        if (height_at(mid) < h) lo = mid; else hi = mid;
      }
      U = 0.5 * (lo + hi);
    }

    // ---- new shape, under-relaxed -----------------------------------------
    ctx.U = U;
    const RawShape rs = integrate_shape(ctx, params_.r_contact);
    if (!rs.ok) {
      sol.status = SolveStatus::ShapeIntegrationFailed;
      sol.iterations = it + 1;
      break;
    }
    const MeniscusShape fresh =
        resample(rs, params_.z_contact, params_.n_nodes, params_.apex_clustering);

    Real disp = 0.0;
    for (std::size_t i = 0; i < shape.nodes.size(); ++i)
      disp = std::max(disp, norm(fresh.nodes[i] - shape.nodes[i]));
    disp /= params_.r_contact;

    // The shape-field coupling is positively fed back -- a sharper meniscus
    // raises the apex field, which sharpens it further -- so back off the
    // relaxation whenever an iteration fails to reduce the residual.
    if (disp > prev_disp) relax = std::max(0.15, 0.5 * relax);
    prev_disp = disp;

    MeniscusShape blended = fresh;
    for (std::size_t i = 0; i < shape.nodes.size(); ++i)
      blended.nodes[i] = shape.nodes[i] + relax * (fresh.nodes[i] - shape.nodes[i]);
    // Pin the contact line exactly -- that IS a boundary condition.  The apex is
    // NOT: its height is an output of the ODE.  Forcing it to the requested h
    // would leave a permanent mismatch whenever the voltage cannot deliver
    // exactly that height (e.g. when U is clamped at zero because the feed
    // pressure alone already overshoots), and the iteration would stall on that
    // artificial offset instead of converging.
    blended.nodes.front() = {params_.r_contact, params_.z_contact};
    blended.nodes.back().r = 0.0;
    shape = blended;

    sol.iterations = it + 1;
    sol.residual = disp;
    sol.status = SolveStatus::NotConverged;
    if (params_.verbose)
      std::printf("    outer %2d: U = %10.2f V, residual = %.3e\n", it + 1, U, disp);
    if (disp < params_.tol) { sol.status = SolveStatus::Converged; break; }
  }

  // Final field evaluation on the converged shape.
  Mesh full = merge({electrodes_, shape.to_mesh(0.0)});
  bem_.set_mesh(std::move(full));
  bem_.solve({U, 0.0, 0.0});

  sol.shape = shape;
  sol.voltage = U;
  sol.peak_field = bem_.peak_field(Tag::FreeSurface);
  for (Index i = 0; i < bem_.size(); ++i) {
    const Element& el = bem_.mesh().elems[static_cast<std::size_t>(i)];
    if (el.tag == Tag::FreeSurface && el.mid.r < 0.3 * params_.r_contact)
      sol.apex_field = std::max(sol.apex_field, std::abs(bem_.En(i)));
  }
  last_ = shape;
  have_last_ = true;
  return sol;
}

const char* to_string(BranchSide s) {
  switch (s) {
    case BranchSide::LowerHeight: return "lower_height";
    case BranchSide::UpperHeight: return "upper_height";
    default: return "unspecified";
  }
}

const char* to_string(BranchTermination s) {
  switch (s) {
    case BranchTermination::ReachedRequestedRange: return "reached_requested_range";
    case BranchTermination::SolverStopped: return "solver_stopped";
    default: return "not_traced";
  }
}

MeniscusSolution MeniscusSolver::solve_at_voltage(Real U, Real h_max, BranchSide side,
                                                  int scout_steps) {
  MeniscusSolution bad;
  bad.target_voltage = U;
  bad.delta_p = params_.delta_p;
  bad.side = side;

  const int requested = std::max(4, scout_steps);
  const std::vector<MeniscusSolution> scout =
      continuation(0.1 * params_.r_contact, h_max, requested);

  // The continuation stops early when the shape solver gives up.  That is a
  // different situation from having covered the whole requested range, and the
  // two must not be confused.
  bad.termination_reason = (static_cast<int>(scout.size()) >= requested && !scout.empty() &&
                            scout.back().ok())
                               ? BranchTermination::ReachedRequestedRange
                               : BranchTermination::SolverStopped;

  std::vector<const MeniscusSolution*> pts;
  for (const MeniscusSolution& m : scout)
    if (m.ok()) pts.push_back(&m);
  if (pts.size() < 2) {
    bad.status = SolveStatus::BranchCoverageIncomplete;
    return bad;
  }

  // --- crossings inside the traced range -----------------------------------
  // Each bracket carries the voltages that were already computed for its
  // endpoints.  Re-solving them here would throw away the continuation's warm
  // start and can fail near the far end of the branch, which would turn a
  // perfectly good bracket into a spurious NotConverged.
  // The seed is the converged shape at h_lo.  Without it the bisection would
  // warm-start from whatever the scout left behind -- after a continuation that
  // ended on a failed step that is a distorted shape, and the bisection then
  // fails on a bracket that is perfectly good.
  struct Bracket { Real h_lo, h_hi, v_lo, v_hi; MeniscusShape seed; };
  std::vector<Bracket> brackets;
  for (std::size_t i = 1; i < pts.size(); ++i) {
    const Real a0 = pts[i - 1]->voltage - U;
    const Real a1 = pts[i]->voltage - U;
    if (a0 == 0.0 || (a0 * a1) < 0.0)
      brackets.push_back({pts[i - 1]->shape.height, pts[i]->shape.height,
                          pts[i - 1]->voltage, pts[i]->voltage, pts[i - 1]->shape});
  }

  // The branch may already be above U at the smallest sampled height; the
  // crossing then lies below it.  March down rather than returning the smallest
  // sampled shape.
  bool low_end_covered = pts.front()->voltage <= U;
  if (!low_end_covered) {
    Real h = pts.front()->shape.height;
    for (int i = 0; i < 30 && h > 1e-4 * params_.r_contact; ++i) {
      const Real h_next = 0.5 * h;
      MeniscusSolution m = solve_at_height(h_next);
      if (m.ok() && m.voltage <= U) {
        brackets.insert(brackets.begin(),
                        {h_next, h, m.voltage, pts.front()->voltage, m.shape});
        low_end_covered = true;
        break;
      }
      h = h_next;
    }
  }

  // --- is the traced range sufficient to decide? ---------------------------
  // It is, for this target, when the branch has turned over and moved away
  // BELOW the target at the far end: the examined monotone segments then run
  // from below the target, up over the turning point, and back below it.
  // Anything else means we stopped looking, not that the branch stops.
  const bool descending_at_end = pts.back()->voltage < pts[pts.size() - 2]->voltage;
  const bool far_end_below = pts.back()->voltage < U;
  const bool coverage = low_end_covered && descending_at_end && far_end_below;

  bad.crossings_in_range = static_cast<int>(brackets.size());
  bad.coverage_complete = coverage;
  bad.additional_crossing_possible = !coverage;

  const int n = static_cast<int>(brackets.size());

  // --- decide what may be returned -----------------------------------------
  if (n == 0) {
    // No crossing here.  Only a fully investigated branch may claim there is
    // none at all.
    bad.status = coverage ? SolveStatus::VoltageNotBracketed
                          : SolveStatus::BranchCoverageIncomplete;
    return bad;
  }

  Bracket use{};
  switch (side) {
    case BranchSide::Unspecified:
      if (n >= 2) { bad.status = SolveStatus::AmbiguousBranch; return bad; }
      if (!coverage) { bad.status = SolveStatus::BranchCoverageIncomplete; return bad; }
      use = brackets.front();
      break;

    case BranchSide::LowerHeight:
      // The lowest crossing is bracketed and may be delivered.  Whether further
      // solutions exist is reported, not hidden.
      use = brackets.front();
      break;

    case BranchSide::UpperHeight:
      if (n >= 2) {
        use = brackets.back();
      } else if (coverage) {
        // Exactly one crossing on a fully investigated branch: it is the only
        // solution, so delivering it is not a substitution.
        use = brackets.front();
      } else {
        // The upper crossing was never bracketed.  Do not hand back the lower
        // one in its place.
        bad.status = SolveStatus::BranchCoverageIncomplete;
        return bad;
      }
      break;
  }

  // --- bisect inside the chosen bracket ------------------------------------
  Real lo = use.h_lo, hi = use.h_hi;
  const bool rising = use.v_hi > use.v_lo;

  MeniscusSolution best;
  Real best_err = 1e300;
  const Real vtol = params_.voltage_tol * std::max(std::abs(U), 1.0);
  for (int i = 0; i < 40; ++i) {
    const Real mid = 0.5 * (lo + hi);
    MeniscusSolution m = (i == 0) ? solve_at_height(mid, &use.seed) : solve_at_height(mid);
    if (!m.ok()) {
      // Cannot tell which side the crossing is on; give up the half that is
      // further from the best candidate so far.
      if (best_err < 1e299 && mid > best.shape.height) hi = mid; else lo = mid;
      if (hi - lo < 1e-9 * params_.r_contact) break;
      continue;
    }
    // Keep the closest candidate, not merely the last converged one: the last
    // one can be worse than an earlier midpoint if the interval later collapses.
    const Real err = std::abs(m.voltage - U);
    if (err < best_err) { best_err = err; best = m; }
    if (err <= vtol) break;
    const bool below = m.voltage < U;
    if (below == rising) lo = mid; else hi = mid;
    if (hi - lo < 1e-9 * params_.r_contact) break;
  }

  best.target_voltage = U;
  best.crossings_in_range = n;
  best.coverage_complete = coverage;
  best.additional_crossing_possible = !coverage;
  best.termination_reason = bad.termination_reason;
  best.side = side;
  if (best.status == SolveStatus::NotAttempted) best.status = SolveStatus::NotConverged;
  if (!best.ok()) return best;
  if (std::abs(best.voltage - U) > vtol) best.status = SolveStatus::VoltageMismatch;
  return best;
}

void MeniscusSolver::realize(const MeniscusSolution& sol) {
  if (sol.shape.nodes.size() < 2)
    throw std::runtime_error("MeniscusSolver::realize: solution carries no free surface");
  Mesh full = merge({electrodes_, sol.shape.to_mesh(0.0)});
  bem_.set_mesh(std::move(full));
  bem_.solve({sol.voltage, 0.0, 0.0});
  last_ = sol.shape;
  have_last_ = true;
}

std::vector<MeniscusSolution> MeniscusSolver::continuation(Real h_min, Real h_max, int n_steps) {
  std::vector<MeniscusSolution> branch;
  branch.reserve(static_cast<std::size_t>(n_steps));
  have_last_ = false;
  for (int i = 0; i < n_steps; ++i) {
    const Real t = (n_steps > 1) ? static_cast<Real>(i) / (n_steps - 1) : 0.0;
    const Real h = h_min + t * (h_max - h_min);
    if (params_.verbose) std::printf("  h = %.4g m (%.3f r_c)\n", h, h / params_.r_contact);
    MeniscusSolution s = solve_at_height(h);
    branch.push_back(s);
    if (!s.ok() && i > 2) {
      if (params_.verbose) std::printf("  branch stopped: no converged shape at h/r_c = %.3f\n",
                                       h / params_.r_contact);
      break;
    }
  }
  return branch;
}

MeniscusSolver::StaticFold MeniscusSolver::find_static_fold(
    const std::vector<MeniscusSolution>& branch) {
  StaticFold f;

  std::vector<std::size_t> idx;
  for (std::size_t i = 0; i < branch.size(); ++i)
    if (branch[i].ok()) idx.push_back(i);

  if (idx.size() < 3) {
    // One or two points cannot exhibit a maximum.  The prototype reported one
    // anyway, from a single converged point.
    f.status = FoldStatus::TooFewPoints;
    return f;
  }

  std::size_t kbest = 0;
  for (std::size_t k = 1; k < idx.size(); ++k)
    if (branch[idx[k]].voltage > branch[idx[kbest]].voltage) kbest = k;

  if (kbest == 0 || kbest + 1 == idx.size()) {
    // The largest value sits at the edge of the traced range: either the branch
    // is monotone, or it was not followed far enough.  Neither demonstrates an
    // interior turning point.
    const bool rising = branch[idx.back()].voltage > branch[idx.front()].voltage;
    const bool falling = branch[idx.back()].voltage < branch[idx.front()].voltage;
    bool monotone = true;
    for (std::size_t k = 1; k < idx.size(); ++k) {
      const Real d = branch[idx[k]].voltage - branch[idx[k - 1]].voltage;
      if ((rising && d < 0.0) || (falling && d > 0.0)) { monotone = false; break; }
    }
    f.status = monotone ? FoldStatus::Monotone : FoldStatus::MaximumAtBoundary;
    return f;
  }

  const MeniscusSolution& mm = branch[idx[kbest]];
  const MeniscusSolution& ml = branch[idx[kbest - 1]];
  const MeniscusSolution& mr = branch[idx[kbest + 1]];

  // Require a strict interior maximum: rising into it, falling out of it.
  if (!(ml.voltage < mm.voltage && mr.voltage < mm.voltage)) {
    f.status = FoldStatus::MaximumAtBoundary;
    return f;
  }

  f.status = FoldStatus::Found;
  f.index = idx[kbest];
  f.voltage = mm.voltage;
  f.height = mm.shape.height;
  f.apex_field = mm.apex_field;
  f.apex_radius = mm.shape.apex_radius;
  f.half_angle = mm.shape.half_angle;

  // Parabolic refinement through the three points around the maximum; without
  // it the fold voltage is biased low by up to half a continuation step.
  const Real h0 = ml.shape.height, v0 = ml.voltage;
  const Real h1 = mm.shape.height, v1 = mm.voltage;
  const Real h2 = mr.shape.height, v2 = mr.voltage;
  if (h1 > h0 && h2 > h1) {
    const Real d1 = (v1 - v0) / (h1 - h0);
    const Real d2 = (v2 - v1) / (h2 - h1);
    const Real a2 = (d2 - d1) / (h2 - h0);
    if (a2 < 0.0) {
      const Real hstar = 0.5 * (h0 + h1) - d1 / (2.0 * a2);
      if (hstar > h0 && hstar < h2) {
        f.height = hstar;
        f.voltage = v0 + d1 * (hstar - h0) + a2 * (hstar - h0) * (hstar - h1);
      }
    }
  }
  return f;
}

void MeniscusSolver::write_branch_csv(const std::string& path,
                                      const std::vector<MeniscusSolution>& branch,
                                      const std::string& header) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  if (!header.empty()) std::fputs(header.c_str(), f);
  std::fprintf(f, "# the maximum of `voltage` over converged rows is the STATIC FOLD,\n");
  std::fprintf(f, "# not an emission onset -- see docs/02_model_specification.md 2.4.\n");
  std::fprintf(f, "height,voltage,apex_field,peak_field,apex_radius,half_angle_deg,arclength,"
                  "status,iterations,residual\n");
  for (const MeniscusSolution& s : branch)
    std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,%.9e,%.6f,%.9e,%s,%d,%.3e\n", s.shape.height, s.voltage,
                 s.apex_field, s.peak_field, s.shape.apex_radius,
                 s.shape.half_angle * 180.0 / pi, s.shape.arclength, to_string(s.status),
                 s.iterations, s.residual);
  std::fclose(f);
}

}  // namespace es
