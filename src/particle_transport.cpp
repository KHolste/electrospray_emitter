#include "es/particle_transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace es {

const char* to_string(Fate f) {
  switch (f) {
    case Fate::Flying: return "Flying";
    case Fate::HitEmitter: return "HitEmitter";
    case Fate::HitPolymer: return "HitPolymer";
    case Fate::HitExtractor: return "HitExtractor";
    case Fate::LeftThroughAperture: return "LeftThroughAperture";
    case Fate::LeftDomain: return "LeftDomain";
  }
  return "?";
}

bool TransportDomain::inside(Vec2 x) const {
  return x.r <= R && x.z >= 0.0 && x.z <= Z;
}

Fate TransportDomain::classify(Vec2 x) const {
  // The order matters and is stated: the ends are checked before the side, so
  // a particle that leaves through the aperture is not counted as a wall hit
  // just because it also passed r = R in the same step.
  if (x.z <= 0.0) return Fate::HitEmitter;
  if (x.z >= Z) return (x.r <= aperture_radius) ? Fate::LeftThroughAperture : Fate::HitExtractor;
  if (x.r >= R) return Fate::HitPolymer;
  return Fate::Flying;
}

void TransportBalance::print(std::FILE* out) const {
  std::fprintf(out, "  Strombilanz [A]:\n");
  std::fprintf(out, "    gestartet     %.9e  (%lld Teilchen)\n", launched,
               static_cast<long long>(n_launched));
  std::fprintf(out, "    extrahiert    %.9e  (%lld)\n", extracted,
               static_cast<long long>(n_extracted));
  std::fprintf(out, "    abgefangen    %.9e  (%lld)   davon Emitter %.3e, Polymer %.3e, "
                    "Extraktor %.3e\n",
               intercepted, static_cast<long long>(n_intercepted), hit_emitter, hit_polymer,
               hit_extractor);
  std::fprintf(out, "    noch fliegend %.9e  (%lld)\n", still_flying,
               static_cast<long long>(n_flying));
  std::fprintf(out, "    Schliessfehler %.3e, groesster Energiefehler %.3e\n", closure_error,
               max_energy_error);
}

// ---------------------------------------------------------------------------

TransportResult transport_particles(const QuadMesh& mesh, const std::vector<Real>& phi,
                                    const std::vector<Real>& eps_r,
                                    const std::vector<char>& active,
                                    const TransportDomain& domain,
                                    const TransportSpecies& species,
                                    const std::vector<TracedParticle>& launch, Real dt,
                                    int max_steps) {
  if (!species.complete())
    throw std::runtime_error(
        "transport_particles: die Spezies ist unvollstaendig.  Ladung und Masse sind "
        "PFLICHTANGABEN mit Vorzeichen; es gibt keine Vorgabespezies.");
  if (!(dt > 0.0) || max_steps < 1)
    throw std::runtime_error("transport_particles: dt muss positiv und max_steps >= 1 sein");

  TransportResult out;
  out.dt = dt;
  out.max_steps = max_steps;
  out.particles = launch;
  const Real qm = species.charge_to_mass();

  auto accel = [&](Vec2 x) {
    const Vec2 E = interpolated_field(mesh, phi, eps_r, active, x);
    return qm * E;
  };
  auto pot = [&](Vec2 x) {
    Index i, j;
    Real xi, eta;
    if (!locate(mesh, x, &i, &j, &xi, &eta)) return 0.0;
    return potential_in_cell(mesh, phi, i, j, xi, eta);
  };

  for (std::size_t pk = 0; pk < out.particles.size(); ++pk) {
    TracedParticle& p = out.particles[pk];
    const Vec2 v_launch = launch[pk].v;
    p.fate = Fate::Flying;
    p.time = 0.0;
    p.steps = 0;
    p.path_length = 0.0;
    p.x_start = p.x;
    p.potential_start = pot(p.x);
    // THE ENERGY INVARIANT.  The force is F = q E = -q grad phi, so the
    // potential energy is U = +q phi and the constant of the motion is
    //     E = 1/2 m v^2 + q phi .
    // The first version of this file wrote a minus sign there; the energy check
    // then reported a relative error of 2, which is exactly what an inverted
    // potential energy produces, and the test caught it.
    const Real E0 = 0.5 * species.mass * norm2(p.v) + species.charge * p.potential_start;

    Vec2 a = accel(p.x);
    for (int s = 0; s < max_steps; ++s) {
      // Velocity Verlet.  For a position-dependent force this is time
      // reversible and does not drift in energy.
      const Vec2 x_new = p.x + dt * p.v + (0.5 * dt * dt) * a;
      Vec2 x_use = x_new;
      // The axis is a symmetry line, not a wall.  A particle that crosses it
      // in the meridian plane has passed through the axis in three dimensions;
      // reflecting r and v_r is the axisymmetric statement of that, and it is
      // NOT a hit.
      Vec2 v_half = p.v + (0.5 * dt) * a;
      if (x_use.r < 0.0) {
        x_use.r = -x_use.r;
        v_half.r = -v_half.r;
      }
      const Fate f = domain.classify(x_use);
      if (f != Fate::Flying) {
        // A PARTIAL LAST STEP, not a full one that is then clipped.  Bisect on
        // the step FRACTION s so that the Verlet position
        //     x(s) = x + s dt v + (s dt)^2 a / 2
        // is the last point still inside, then take a consistent Verlet step of
        // length s dt.  A first version clipped the position but kept the
        // velocity of the full step; the exit speed was then wrong by up to
        // 2.5e-3 and the measured time order came out negative, which is how
        // the test found it.
        Real lo_s = 0.0, hi_s = 1.0;
        auto step_to = [&](Real ss) { return p.x + (ss * dt) * p.v + (0.5 * ss * ss * dt * dt) * a; };
        for (int it = 0; it < 60; ++it) {
          const Real mid = 0.5 * (lo_s + hi_s);
          Vec2 xm = step_to(mid);
          if (xm.r < 0.0) xm.r = -xm.r;
          if (domain.classify(xm) == Fate::Flying)
            lo_s = mid;
          else
            hi_s = mid;
        }
        const Real dts = lo_s * dt;
        Vec2 x_exit = step_to(lo_s);
        Vec2 v_h = p.v + (0.5 * dts) * a;
        if (x_exit.r < 0.0) {
          x_exit.r = -x_exit.r;
          v_h.r = -v_h.r;
        }
        p.v = v_h + (0.5 * dts) * accel(x_exit);
        p.fate = f;
        p.x_end = x_exit;
        p.path_length += norm(x_exit - p.x);
        p.time += dts;
        p.x = x_exit;
        ++p.steps;
        break;
      }
      const Vec2 a_new = accel(x_use);
      p.v = v_half + (0.5 * dt) * a_new;
      p.path_length += norm(x_use - p.x);
      p.x = x_use;
      a = a_new;
      p.time += dt;
      ++p.steps;
      // Energy check along the way, on the CURRENT position.
      const Real E = 0.5 * species.mass * norm2(p.v) + species.charge * pot(p.x);
      const Real scale = std::max(std::abs(E0), 0.5 * species.mass * norm2(p.v));
      if (scale > 0.0)
        out.balance.max_energy_error =
            std::max(out.balance.max_energy_error, std::abs(E - E0) / scale);
    }
    if (p.fate == Fate::Flying) p.x_end = p.x;
    p.potential_end = pot(p.x_end);
    p.energy_gain = 0.5 * species.mass * (norm2(p.v) - norm2(v_launch));
  }

  // --- the balance, and it must close ---------------------------------------
  TransportBalance& b = out.balance;
  for (const TracedParticle& p : out.particles) {
    b.launched += p.current;
    ++b.n_launched;
    switch (p.fate) {
      case Fate::LeftThroughAperture:
        b.extracted += p.current;
        ++b.n_extracted;
        break;
      case Fate::HitEmitter:
        b.hit_emitter += p.current;
        b.intercepted += p.current;
        ++b.n_intercepted;
        break;
      case Fate::HitPolymer:
        b.hit_polymer += p.current;
        b.intercepted += p.current;
        ++b.n_intercepted;
        break;
      case Fate::HitExtractor:
        b.hit_extractor += p.current;
        b.intercepted += p.current;
        ++b.n_intercepted;
        break;
      case Fate::LeftDomain:
        b.extracted += p.current;
        ++b.n_extracted;
        break;
      case Fate::Flying:
        b.still_flying += p.current;
        ++b.n_flying;
        break;
    }
  }
  const Real sum = b.extracted + b.intercepted + b.still_flying;
  b.closure_error =
      (std::abs(b.launched) > 0.0) ? std::abs(sum - b.launched) / std::abs(b.launched) : 0.0;
  out.message = "Transportantwort auf eine VORGESCHRIEBENE Startverteilung.  Keine "
                "Stromvorhersage: es gibt keine physikalische Teilchenquelle (P5 blockiert).";
  return out;
}

}  // namespace es
