#include "es/surface_kinematics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

using constants::pi;

// ---------------------------------------------------------------------------

const char* to_string(KinematicMode m) {
  switch (m) {
    case KinematicMode::Lagrangian: return "Lagrangian";
    case KinematicMode::NormalOnly: return "NormalOnly";
  }
  return "?";
}

const char* to_string(ContactLine c) {
  switch (c) {
    case ContactLine::Pinned: return "Pinned";
    case ContactLine::Free: return "Free";
  }
  return "?";
}

const char* to_string(StepStatus s) {
  switch (s) {
    case StepStatus::Ok: return "Ok";
    case StepStatus::NodeCrossedAxis: return "NodeCrossedAxis";
    case StepStatus::SegmentCollapsed: return "SegmentCollapsed";
    case StepStatus::StepTooLarge: return "StepTooLarge";
    case StepStatus::NotAttempted: return "NotAttempted";
  }
  return "?";
}

const char* explain(StepStatus s) {
  switch (s) {
    case StepStatus::Ok:
      return "Der Schritt ist zulaessig: kein Knoten hat die Achse ueberquert, kein Segment "
             "ist kollabiert, und die Knotenbewegung bleibt unter der Diskretisierungs"
             "schranke.";
    case StepStatus::NodeCrossedAxis:
      return "Ein Knoten hat r < 0 erreicht.  Die Flaeche waere keine achsensymmetrische "
             "Flaeche mehr; der Schritt wird verworfen und nicht mit einer Warnung "
             "angenommen.";
    case StepStatus::SegmentCollapsed:
      return "Ein Segment hat seine Laenge verloren oder sich umgekehrt.  Der Polygonzug "
             "waere keine einfache Kurve mehr.";
    case StepStatus::StepTooLarge:
      return "Ein Knoten hat sich in einem Schritt weiter bewegt als die zugelassene "
             "Bruchteil der kuerzesten Segmentlaenge.  Das ist eine "
             "DISKRETISIERUNGSschranke, keine physikalische: hier wird keine Kraft "
             "integriert.";
    case StepStatus::NotAttempted:
      return "Es wurde kein Schritt versucht.";
  }
  return "?";
}

// ---------------------------------------------------------------------------

Vec2 SurfacePolyline::tangent_at(std::size_t k) const {
  const std::size_t n = nodes.size();
  if (n < 2) return {0, 0};
  Vec2 t;
  if (k == 0)
    t = nodes[1] - nodes[0];
  else if (k + 1 == n)
    t = nodes[n - 1] - nodes[n - 2];
  else
    t = nodes[k + 1] - nodes[k - 1];
  return normalized(t);
}

Vec2 SurfacePolyline::normal_at(std::size_t k) const {
  // On the axis, axisymmetry forces the normal to be axial: any radial
  // component would single out a direction the problem does not have.
  if (k == 0 && std::abs(nodes[0].r) < 1e-15) return {0.0, 1.0};
  const Vec2 t = tangent_at(k);
  // perp(t) = (t.z, -t.r); for a surface running apex -> contact line with the
  // liquid below, that points outward.
  return normalized(perp(t));
}

Real SurfacePolyline::arclength() const {
  Real s = 0.0;
  for (std::size_t k = 0; k + 1 < nodes.size(); ++k) s += norm(nodes[k + 1] - nodes[k]);
  return s;
}

Real SurfacePolyline::revolved_volume(Real z_base) const {
  // Exact revolved volume of the polyline over the plane z = z_base, by the
  // truncated-cone formula per segment: no quadrature error.
  Real v = 0.0;
  for (std::size_t k = 0; k + 1 < nodes.size(); ++k) {
    const Real r0 = nodes[k].r, r1 = nodes[k + 1].r;
    const Real z0 = nodes[k].z - z_base, z1 = nodes[k + 1].z - z_base;
    v += pi / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1) * (z0 - z1);
  }
  return v;
}

Real SurfacePolyline::revolved_area() const {
  Real a = 0.0;
  for (std::size_t k = 0; k + 1 < nodes.size(); ++k) {
    const Real r0 = nodes[k].r, r1 = nodes[k + 1].r;
    a += pi * (r0 + r1) * norm(nodes[k + 1] - nodes[k]);
  }
  return a;
}

Real SurfacePolyline::node_spacing_nonuniformity() const {
  const std::size_t n = nodes.size();
  if (n < 3) return 0.0;
  Real lo = std::numeric_limits<Real>::max(), hi = 0.0, sum = 0.0;
  for (std::size_t k = 0; k + 1 < n; ++k) {
    const Real d = norm(nodes[k + 1] - nodes[k]);
    lo = std::min(lo, d);
    hi = std::max(hi, d);
    sum += d;
  }
  const Real mean = sum / static_cast<Real>(n - 1);
  return (mean > 0.0) ? (hi - lo) / mean : 0.0;
}

// ---------------------------------------------------------------------------

void AdvectionResult::print(std::FILE* out) const {
  std::fprintf(out, "  Advektion: %s, %d Schritte a %.4e s\n", to_string(status), steps, dt);
  std::fprintf(out, "    Volumen %.9e -> %.9e, relative Aenderung %.3e\n", volume_initial,
               volume_final, volume_change);
  std::fprintf(out, "    groesste Knotenbewegung / kuerzestes Segment: %.3f (Grenze %.2f)\n",
               max_node_motion_fraction, kinematics::kMaxNodeMotion);
  std::fprintf(out, "    Ungleichmaessigkeit der Knotenabstaende: %.3e\n",
               spacing_nonuniformity);
  if (!message.empty()) std::fprintf(out, "    %s\n", message.c_str());
}

namespace {

/// Redistribute the nodes to equal arclength along the SAME polyline.  This is
/// mesh motion: it moves nodes along the curve and therefore cannot change the
/// curve.  The test measures that it does not.
void redistribute(SurfacePolyline& s, bool pin_ends) {
  const std::size_t n = s.nodes.size();
  if (n < 3) return;
  std::vector<Real> cum(n, 0.0);
  for (std::size_t k = 1; k < n; ++k) cum[k] = cum[k - 1] + norm(s.nodes[k] - s.nodes[k - 1]);
  const Real L = cum[n - 1];
  if (!(L > 0.0)) return;
  std::vector<Vec2> out(n);
  out[0] = s.nodes[0];
  out[n - 1] = s.nodes[n - 1];
  for (std::size_t k = 1; k + 1 < n; ++k) {
    const Real target = L * static_cast<Real>(k) / static_cast<Real>(n - 1);
    std::size_t j = 0;
    while (j + 2 < n && cum[j + 1] < target) ++j;
    const Real seg = cum[j + 1] - cum[j];
    const Real w = (seg > 0.0) ? (target - cum[j]) / seg : 0.0;
    out[k] = (1.0 - w) * s.nodes[j] + w * s.nodes[j + 1];
  }
  (void)pin_ends;  // the ends are kept in both cases; pinning is applied outside
  s.nodes = out;
}

}  // namespace

AdvectionResult advect_surface(const SurfacePolyline& initial, const VelocityField& u, Real dt,
                               int n_steps, KinematicMode mode, ContactLine contact,
                               Real z_base) {
  AdvectionResult out;
  out.dt = dt;
  out.surface = initial;
  out.volume_initial = initial.revolved_volume(z_base);
  if (n_steps < 0 || !(dt > 0.0)) {
    out.status = StepStatus::NotAttempted;
    out.message = "dt muss positiv und n_steps nicht negativ sein.";
    out.volume_final = out.volume_initial;
    return out;
  }

  const std::size_t n = out.surface.nodes.size();
  const Vec2 pinned_end = out.surface.nodes.back();

  for (int step = 0; step < n_steps; ++step) {
    const SurfacePolyline before = out.surface;
    const Real t = out.surface.time;

    // Shortest segment before the step: the yardstick of the CFL-like bound.
    Real shortest = std::numeric_limits<Real>::max();
    for (std::size_t k = 0; k + 1 < n; ++k)
      shortest = std::min(shortest, norm(before.nodes[k + 1] - before.nodes[k]));

    // Classical RK4 on the node trajectories.  In NormalOnly the velocity is
    // projected onto the normal of the CURRENT surface at every stage; the
    // normal is a property of the surface, so it is evaluated on `before`.
    auto vel = [&](std::size_t k, Vec2 x, Real tt) {
      const Vec2 v = u(x, tt);
      if (mode == KinematicMode::Lagrangian) return v;
      const Vec2 nrm = before.normal_at(k);
      return dot(v, nrm) * nrm;
    };

    std::vector<Vec2> moved(n);
    for (std::size_t k = 0; k < n; ++k) {
      const Vec2 x0 = before.nodes[k];
      const Vec2 k1 = vel(k, x0, t);
      const Vec2 k2 = vel(k, x0 + (0.5 * dt) * k1, t + 0.5 * dt);
      const Vec2 k3 = vel(k, x0 + (0.5 * dt) * k2, t + 0.5 * dt);
      const Vec2 k4 = vel(k, x0 + dt * k3, t + dt);
      moved[k] = x0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    if (contact == ContactLine::Pinned) moved.back() = pinned_end;
    // The apex sits on the axis and must stay there: axisymmetry, not a choice.
    if (std::abs(before.nodes[0].r) < 1e-15) moved[0].r = 0.0;

    // --- admissibility, checked BEFORE the step is accepted -----------------
    Real worst_motion = 0.0;
    for (std::size_t k = 0; k < n; ++k)
      worst_motion = std::max(worst_motion, norm(moved[k] - before.nodes[k]));
    out.max_node_motion_fraction =
        std::max(out.max_node_motion_fraction, (shortest > 0.0) ? worst_motion / shortest : 0.0);

    for (std::size_t k = 0; k < n; ++k)
      if (moved[k].r < 0.0) {
        out.status = StepStatus::NodeCrossedAxis;
        out.message = explain(out.status);
        out.volume_final = out.surface.revolved_volume(z_base);
        out.volume_change = (std::abs(out.volume_initial) > 0.0)
                                ? (out.volume_final - out.volume_initial) / out.volume_initial
                                : 0.0;
        out.spacing_nonuniformity = out.surface.node_spacing_nonuniformity();
        return out;
      }
    for (std::size_t k = 0; k + 1 < n; ++k)
      if (!(norm(moved[k + 1] - moved[k]) > kinematics::kMinSegmentFraction * shortest)) {
        out.status = StepStatus::SegmentCollapsed;
        out.message = explain(out.status);
        out.volume_final = out.surface.revolved_volume(z_base);
        out.volume_change = (std::abs(out.volume_initial) > 0.0)
                                ? (out.volume_final - out.volume_initial) / out.volume_initial
                                : 0.0;
        out.spacing_nonuniformity = out.surface.node_spacing_nonuniformity();
        return out;
      }
    if (shortest > 0.0 && worst_motion > kinematics::kMaxNodeMotion * shortest) {
      out.status = StepStatus::StepTooLarge;
      out.message = explain(out.status);
      out.volume_final = out.surface.revolved_volume(z_base);
      out.volume_change = (std::abs(out.volume_initial) > 0.0)
                              ? (out.volume_final - out.volume_initial) / out.volume_initial
                              : 0.0;
      out.spacing_nonuniformity = out.surface.node_spacing_nonuniformity();
      return out;
    }

    out.surface.nodes = moved;
    out.surface.time = t + dt;
    ++out.steps;
    if (mode == KinematicMode::NormalOnly) redistribute(out.surface, true);
  }

  out.status = StepStatus::Ok;
  out.message = explain(out.status);
  out.volume_final = out.surface.revolved_volume(z_base);
  out.volume_change = (std::abs(out.volume_initial) > 0.0)
                          ? (out.volume_final - out.volume_initial) / out.volume_initial
                          : 0.0;
  out.spacing_nonuniformity = out.surface.node_spacing_nonuniformity();
  return out;
}

// ---------------------------------------------------------------------------

VelocityField dilation_field(Real alpha) {
  return [alpha](Vec2 x, Real) { return Vec2{alpha * x.r, alpha * x.z}; };
}

VelocityField squeeze_field(Real alpha) {
  // div u = (1/r) d(r u_r)/dr + du_z/dz = -alpha + alpha = 0, exactly.
  return [alpha](Vec2 x, Real) { return Vec2{-0.5 * alpha * x.r, alpha * x.z}; };
}

Vec2 dilation_exact(Vec2 x0, Real alpha, Real t) {
  const Real e = std::exp(alpha * t);
  return {x0.r * e, x0.z * e};
}

Vec2 squeeze_exact(Vec2 x0, Real alpha, Real t) {
  return {x0.r * std::exp(-0.5 * alpha * t), x0.z * std::exp(alpha * t)};
}

SurfacePolyline hemisphere(Real R, std::size_t n) {
  if (n < 3) throw std::runtime_error("hemisphere: zu wenige Knoten");
  SurfacePolyline s;
  s.nodes.resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    const Real th = 0.5 * pi * static_cast<Real>(k) / static_cast<Real>(n - 1);
    s.nodes[k] = Vec2{R * std::sin(th), R * std::cos(th)};
  }
  s.nodes[0].r = 0.0;
  s.nodes[n - 1].z = 0.0;
  return s;
}

void solve_dynamic_meniscus() {
  throw NotImplementedInThisPhase(
      "Der zeitabhaengige Meniskus",
      "nach einem Stroemungsloeser mit freier Oberflaeche und einem "
      "Oberflaechenladungstransport; der vollstaendige Vertrag steht in "
      "docs/15_free_surface_dynamics.md",
      "Die Grundlage fehlt, und zwar strukturell, nicht aus Zeitmangel.  (a) Es gibt kein "
      "Geschwindigkeitsfeld im Inneren mit freier Oberflaeche: die Rohrstroemung aus P3 "
      "traegt das nicht, ihre Exaktheit beruht auf du_z/dz = 0, was eine sich verformende "
      "Oberflaeche zerstoert.  (b) Es gibt keinen Oberflaechenladungstransport mit "
      "Konvektion und keine TANGENTIALE Traktion q_s E_t -- diese existiert im "
      "Perfect-Conductor-Grenzfall gar nicht und ist genau das, was einen endlich "
      "leitenden Meniskus antreibt.  (c) eps_r ist nach dem Stoffdatenvertrag von P2 "
      "MissingMaterialData, also ist nicht einmal die Ladungsrelaxationszeit berechenbar.  "
      "Es gibt hier keine Mobilitaet, keine kuenstliche Daempfung und keine "
      "Stabilitaetsaussage; ein Ersatzmodell waere ein freier Parameter im Kostuem der "
      "Physik.");
}

}  // namespace es
