#include "es/beam.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {
namespace {

using constants::eps0;
using constants::pi;

/// A ring macroparticle standing in for a piece of the steady-state beam.
struct RingCharge {
  Vec2 x{};
  Real Q{0};
};

/// Does the segment a->b cross element e?  Returns the crossing parameter along
/// a->b, or -1.  Plain 2D segment intersection in the meridian plane.
Real segment_hit(Vec2 a, Vec2 b, const Element& e) {
  const Vec2 r = b - a;
  const Vec2 s = e.b - e.a;
  const Real denom = r.r * s.z - r.z * s.r;
  if (std::abs(denom) < 1e-30) return -1.0;
  const Vec2 q = e.a - a;
  const Real t = (q.r * s.z - q.z * s.r) / denom;
  const Real u = (q.r * r.z - q.z * r.r) / denom;
  if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) return t;
  return -1.0;
}

}  // namespace

const char* ray_status_name(RayStatus s) {
  switch (s) {
    case RayStatus::Flying: return "flying";
    case RayStatus::Transmitted: return "transmitted";
    case RayStatus::Intercepted: return "intercepted";
    default: return "escaped";
  }
}

Real ring_potential(Real Q, Vec2 x, Vec2 xp) {
  if (xp.r <= 0.0) {
    // Degenerate ring on the axis: a point charge.
    const Real d = norm(x - xp);
    return (d > 0.0) ? Q / (4.0 * pi * eps0 * d) : 0.0;
  }
  return Q * kernel_G(x, xp) / (2.0 * pi * xp.r);
}

Vec2 ring_field(Real Q, Vec2 x, Vec2 xp) {
  if (xp.r <= 0.0) {
    const Vec2 d = x - xp;
    const Real n = norm(d);
    if (!(n > 0.0)) return {0.0, 0.0};
    return (Q / (4.0 * pi * eps0 * n * n * n)) * d;
  }
  const Vec2 g = kernel_gradG(x, xp);
  return (-Q / (2.0 * pi * xp.r)) * g;
}

// ---------------------------------------------------------------------------

namespace {

/// Push one ray through the field with adaptive leapfrog-style steps.
void push_ray(Ray& ray, const BemSolver& bem, const std::vector<RingCharge>* rings,
              const BeamParams& p) {
  const Mesh& mesh = bem.mesh();
  const Real qm = ray.qm;

  auto field = [&](Vec2 x) {
    Vec2 E = bem.field_at(x);
    if (rings)
      for (const RingCharge& rc : *rings) E += ring_field(rc.Q, x, rc.x);
    return E;
  };

  // Length scale for the step limiter: the smallest element near the launch
  // point sets how fast the field varies there.
  Real scale = 1e30;
  for (const Element& e : mesh.elems) scale = std::min(scale, e.len);
  scale = std::max(scale, 1e-12);

  const int sample_every =
      std::max(1, p.max_steps / std::max(1, p.path_samples));
  ray.path.clear();
  ray.speed.clear();
  ray.path.push_back(ray.x);
  ray.speed.push_back(norm(ray.v));

  for (int step = 0; step < p.max_steps; ++step) {
    const Vec2 E = field(ray.x);
    const Vec2 acc = qm * E;
    const Real a = norm(acc);
    const Real v = norm(ray.v);

    // Step so that neither the velocity change nor the displacement exceeds the
    // CFL fraction of the local scales.
    Real dt = 1e30;
    if (a > 0.0) dt = std::min(dt, std::sqrt(2.0 * p.cfl * scale / a));
    if (v > 0.0) dt = std::min(dt, p.cfl * std::max(scale, 0.02 * norm(ray.x - ray.x0)) / v);
    if (!std::isfinite(dt) || dt <= 0.0) break;

    // velocity Verlet
    const Vec2 x_new = ray.x + dt * ray.v + (0.5 * dt * dt) * acc;
    const Vec2 acc_new = qm * field(x_new);
    const Vec2 v_new = ray.v + (0.5 * dt) * (acc + acc_new);

    // Collision with any electrode along this step.
    Real t_hit = 2.0;
    const Element* hit = nullptr;
    for (const Element& e : mesh.elems) {
      const Real t = segment_hit(ray.x, x_new, e);
      if (t >= 0.0 && t < t_hit) { t_hit = t; hit = &e; }
    }
    if (hit) {
      ray.x = ray.x + t_hit * (x_new - ray.x);
      ray.v = v_new;
      ray.status = RayStatus::Intercepted;
      ray.hit_tag = hit->tag;
      ray.steps = step + 1;
      ray.path.push_back(ray.x);
      ray.speed.push_back(norm(ray.v));
      return;
    }

    ray.x = x_new;
    ray.v = v_new;
    ray.steps = step + 1;
    if (step % sample_every == 0) {
      ray.path.push_back(ray.x);
      ray.speed.push_back(norm(ray.v));
    }

    if (ray.x.z >= p.z_end) { ray.status = RayStatus::Transmitted; break; }
    if (ray.x.r >= p.r_max || ray.x.z < -p.r_max) { ray.status = RayStatus::Escaped; break; }
    if (ray.x.r < 0.0) {  // crossed the axis: reflect, it is a 2D projection
      ray.x.r = -ray.x.r;
      ray.v.r = -ray.v.r;
    }
  }
  ray.path.push_back(ray.x);
  ray.speed.push_back(norm(ray.v));
}

std::vector<RingCharge> rings_from_rays(const std::vector<Ray>& rays) {
  std::vector<RingCharge> out;
  for (const Ray& r : rays) {
    for (std::size_t i = 1; i < r.path.size(); ++i) {
      const Real dl = norm(r.path[i] - r.path[i - 1]);
      const Real vbar = 0.5 * (r.speed[i] + r.speed[i - 1]);
      if (!(vbar > 0.0) || !(dl > 0.0)) continue;
      // Steady state: the charge sitting in this piece of the beam is the
      // current times the time a particle spends crossing it.
      out.push_back({0.5 * (r.path[i] + r.path[i - 1]), r.current * dl / vbar});
    }
  }
  return out;
}

void finalise(BeamResult& res) {
  res.current_launched = 0;
  res.current_transmitted = 0;
  res.current_intercepted = 0;
  Real energy_acc = 0.0;
  for (Ray& r : res.rays) {
    const Real v2 = norm2(r.v);
    r.energy_eV = 0.5 * v2 / std::max(r.qm, 1e-300);  // (1/2 m v^2)/q in volts
    r.angle = std::atan2(r.v.r, std::max(r.v.z, 1e-300));
    res.current_launched += r.current;
    if (r.status == RayStatus::Intercepted) res.current_intercepted += r.current;
    if (r.status == RayStatus::Transmitted) {
      res.current_transmitted += r.current;
      energy_acc += r.current * r.energy_eV;
    }
  }
  res.interception_fraction =
      (res.current_launched > 0) ? res.current_intercepted / res.current_launched : 0.0;
  res.mean_energy_eV = (res.current_transmitted > 0) ? energy_acc / res.current_transmitted : 0.0;

  // Current-weighted divergence quantiles over the transmitted beam.
  std::vector<std::pair<Real, Real>> ang;
  for (const Ray& r : res.rays)
    if (r.status == RayStatus::Transmitted) ang.emplace_back(std::abs(r.angle), r.current);
  std::sort(ang.begin(), ang.end());
  Real acc = 0.0;
  res.half_angle_50 = res.half_angle_95 = 0.0;
  bool got50 = false;
  for (const auto& a : ang) {
    acc += a.second;
    if (!got50 && acc >= 0.50 * res.current_transmitted) { res.half_angle_50 = a.first; got50 = true; }
    if (acc >= 0.95 * res.current_transmitted) { res.half_angle_95 = a.first; break; }
  }
}

}  // namespace

// ---------------------------------------------------------------------------

BeamResult trace_beam_with_weights(BemSolver& bem, const std::vector<Real>& element_current,
                                   const std::vector<BeamSpecies>& species,
                                   const BeamParams& p) {
  if (static_cast<Index>(element_current.size()) != bem.size())
    throw std::runtime_error("trace_beam_with_weights: weight vector length != element count");
  if (species.empty()) throw std::runtime_error("trace_beam_with_weights: no species");

  BeamResult res;
  const Mesh& mesh = bem.mesh();

  // --- launch ---------------------------------------------------------------
  for (Index i = 0; i < bem.size(); ++i) {
    const Real I = element_current[static_cast<std::size_t>(i)];
    if (!(I > 0.0)) continue;
    const Element& e = mesh.elems[static_cast<std::size_t>(i)];
    for (const BeamSpecies& sp : species) {
      if (!(sp.fraction > 0.0) || !(sp.qm > 0.0)) continue;
      Ray r;
      r.x0 = e.mid + (p.launch_offset * e.len) * e.normal;
      r.x = r.x0;
      r.v = {0.0, 0.0};
      r.qm = sp.qm;
      r.current = I * sp.fraction;
      r.species = sp.name;
      res.rays.push_back(std::move(r));
    }
  }
  if (res.rays.empty()) throw std::runtime_error("trace_beam: no emitting elements");

  // --- Laplace pass ---------------------------------------------------------
  const Index nr = static_cast<Index>(res.rays.size());
#ifdef ES_HAVE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for (Index k = 0; k < nr; ++k) push_ray(res.rays[static_cast<std::size_t>(k)], bem, nullptr, p);

  // --- self-consistent space charge ----------------------------------------
  std::vector<RingCharge> rings;
  for (int it = 0; it < p.space_charge_iters; ++it) {
    std::vector<RingCharge> fresh = rings_from_rays(res.rays);
    if (rings.empty()) {
      rings = std::move(fresh);
    } else {
      // Under-relax the ring charges.  The geometry has changed too, so blend
      // by replacing the set and scaling the charges -- crude but stable.
      for (RingCharge& rc : fresh) rc.Q *= p.space_charge_relax;
      for (RingCharge& rc : rings) rc.Q *= (1.0 - p.space_charge_relax);
      fresh.insert(fresh.end(), rings.begin(), rings.end());
      rings = std::move(fresh);
    }

    // Electrodes respond to the beam: impose the ring potential as an external
    // field on the BEM right-hand side and re-solve for the induced charge.
    std::vector<Real> phi(static_cast<std::size_t>(bem.size()), 0.0);
#ifdef ES_HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (Index i = 0; i < bem.size(); ++i) {
      const Vec2 xi = mesh.elems[static_cast<std::size_t>(i)].mid;
      Real v = 0.0;
      for (const RingCharge& rc : rings) v += ring_potential(rc.Q, xi, rc.x);
      phi[static_cast<std::size_t>(i)] = v;
    }
    const std::array<Real, 3> applied = bem.applied();
    bem.set_external_potential(std::move(phi));
    bem.solve(applied);

    for (Ray& r : res.rays) {
      r.x = r.x0;
      r.v = {0.0, 0.0};
      r.status = RayStatus::Flying;
      r.hit_tag = Tag::Other;
    }
#ifdef ES_HAVE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (Index k = 0; k < nr; ++k)
      push_ray(res.rays[static_cast<std::size_t>(k)], bem, &rings, p);

    res.space_charge_iterations = it + 1;
    if (p.verbose) {
      Real Qtot = 0.0;
      for (const RingCharge& rc : rings) Qtot += rc.Q;
      std::printf("  space charge iteration %d: %zu rings, Q_beam = %.4g C\n", it + 1,
                  rings.size(), Qtot);
    }
  }

  finalise(res);
  return res;
}

BeamResult trace_beam(BemSolver& bem, const Fluid& f, Real T,
                      const std::vector<BeamSpecies>& species, const BeamParams& p) {
  std::vector<Real> w(static_cast<std::size_t>(bem.size()), 0.0);
  for (Index i = 0; i < bem.size(); ++i) {
    const Element& e = bem.mesh().elems[static_cast<std::size_t>(i)];
    const bool emits =
        (e.tag == Tag::FreeSurface) || (p.include_wetted_metal && e.tag == Tag::Emitter);
    if (!emits) continue;
    w[static_cast<std::size_t>(i)] = ion_current_density(std::abs(bem.En(i)), f, T) * e.area;
  }
  return trace_beam_with_weights(bem, w, species, p);
}

// ---------------------------------------------------------------------------

void BeamResult::print(std::FILE* out) const {
  std::fprintf(out, "\nbeam transport\n");
  std::fprintf(out, "  rays                : %10zu\n", rays.size());
  std::fprintf(out, "  launched current    : %10.4g A  (= %.4g nA)\n", current_launched,
               current_launched * 1e9);
  std::fprintf(out, "  transmitted         : %10.4g A  (= %.2f %%)\n", current_transmitted,
               current_launched > 0 ? 100.0 * current_transmitted / current_launched : 0.0);
  std::fprintf(out, "  intercepted         : %10.4g A  (= %.2f %%)\n", current_intercepted,
               100.0 * interception_fraction);
  std::fprintf(out, "  half-angle (50%% I)  : %10.2f deg\n", half_angle_50 * 180.0 / constants::pi);
  std::fprintf(out, "  half-angle (95%% I)  : %10.2f deg\n", half_angle_95 * 180.0 / constants::pi);
  std::fprintf(out, "  mean beam energy    : %10.1f eV\n", mean_energy_eV);
  if (space_charge_iterations)
    std::fprintf(out, "  space-charge passes : %10d\n", space_charge_iterations);
}

void BeamResult::write_rays_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "i,species,current,qm,r0,z0,r_end,z_end,vr,vz,angle_deg,energy_eV,status,hit,steps\n");
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const Ray& r = rays[i];
    std::fprintf(f, "%zu,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.6f,%.6f,%s,%s,%d\n", i,
                 r.species.c_str(), r.current, r.qm, r.x0.r, r.x0.z, r.x.r, r.x.z, r.v.r, r.v.z,
                 r.angle * 180.0 / constants::pi, r.energy_eV, ray_status_name(r.status),
                 tag_name(r.hit_tag), r.steps);
  }
  std::fclose(f);
}

void BeamResult::write_paths_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "ray,k,r,z,speed\n");
  for (std::size_t i = 0; i < rays.size(); ++i)
    for (std::size_t k = 0; k < rays[i].path.size(); ++k)
      std::fprintf(f, "%zu,%zu,%.9e,%.9e,%.9e\n", i, k, rays[i].path[k].r, rays[i].path[k].z,
                   k < rays[i].speed.size() ? rays[i].speed[k] : 0.0);
  std::fclose(f);
}

}  // namespace es
