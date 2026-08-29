#include "es/io.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {
namespace {

using constants::pi;

/// Is the point enclosed by one of the meshed bodies?  Crossing number of a ray
/// shot in +r within the meridian half-plane, counted per body: bodies are
/// closed contours (possibly through the axis), so an odd count means inside.
bool inside_conductor(const Mesh& m, Vec2 x) {
  // Bodies that touch the axis are closed through it; extend each such contour
  // with the implicit axial segment by simply not shooting the ray along +r
  // from a point with r < 0 (which cannot happen).  Counting crossings of the
  // meridian contour is then exact.
  std::map<int, int> crossings;
  for (const Element& e : m.elems) {
    const Real z0 = e.a.z, z1 = e.b.z;
    if ((z0 > x.z) == (z1 > x.z)) continue;  // does not straddle the ray
    const Real t = (x.z - z0) / (z1 - z0);
    const Real rc = e.a.r + t * (e.b.r - e.a.r);
    if (rc > x.r) ++crossings[e.body];
  }
  for (const auto& c : crossings)
    if (c.second % 2 == 1) return true;
  return false;
}

}  // namespace

FieldGrid sample_field(const BemSolver& bem, int nr, int nz, Real r0, Real r1, Real z0, Real z1) {
  FieldGrid g;
  g.nr = std::max(2, nr);
  g.nz = std::max(2, nz);
  g.r0 = r0; g.r1 = r1; g.z0 = z0; g.z1 = z1;
  const std::size_t n = static_cast<std::size_t>(g.nr) * static_cast<std::size_t>(g.nz);
  g.V.assign(n, 0.0);
  g.Er.assign(n, 0.0);
  g.Ez.assign(n, 0.0);
  const Real nan = std::numeric_limits<Real>::quiet_NaN();

#ifdef ES_HAVE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for (Index iz = 0; iz < g.nz; ++iz) {
    for (Index ir = 0; ir < g.nr; ++ir) {
      const Vec2 x{g.r_at(static_cast<int>(ir)), g.z_at(static_cast<int>(iz))};
      const std::size_t k = static_cast<std::size_t>(iz) * static_cast<std::size_t>(g.nr) +
                            static_cast<std::size_t>(ir);
      if (inside_conductor(bem.mesh(), x)) {
        g.V[k] = g.Er[k] = g.Ez[k] = nan;
        continue;
      }
      g.V[k] = bem.potential_at(x);
      const Vec2 E = bem.field_at(x);
      g.Er[k] = E.r;
      g.Ez[k] = E.z;
    }
  }
  return g;
}

void write_grid_csv(const FieldGrid& g, const std::string& path) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "r,z,V,Er,Ez,Emag\n");
  for (int iz = 0; iz < g.nz; ++iz)
    for (int ir = 0; ir < g.nr; ++ir) {
      const std::size_t k = static_cast<std::size_t>(iz) * g.nr + ir;
      std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", g.r_at(ir), g.z_at(iz), g.V[k], g.Er[k],
                   g.Ez[k], std::hypot(g.Er[k], g.Ez[k]));
    }
  std::fclose(f);
}

void write_grid_vtk(const FieldGrid& g, const std::string& path) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  const Real dr = (g.nr > 1) ? (g.r1 - g.r0) / (g.nr - 1) : 1.0;
  const Real dz = (g.nz > 1) ? (g.z1 - g.z0) / (g.nz - 1) : 1.0;
  std::fprintf(f, "# vtk DataFile Version 3.0\n");
  std::fprintf(f, "electrospray axisymmetric field (meridian half-plane)\n");
  std::fprintf(f, "ASCII\nDATASET STRUCTURED_POINTS\n");
  std::fprintf(f, "DIMENSIONS %d %d 1\n", g.nr, g.nz);
  std::fprintf(f, "ORIGIN %.9e %.9e 0\n", g.r0, g.z0);
  std::fprintf(f, "SPACING %.9e %.9e 1\n", dr, dz);
  std::fprintf(f, "POINT_DATA %zu\n", g.V.size());
  std::fprintf(f, "SCALARS potential double 1\nLOOKUP_TABLE default\n");
  for (Real v : g.V) std::fprintf(f, "%.9e\n", v);
  std::fprintf(f, "VECTORS E double\n");
  for (std::size_t k = 0; k < g.V.size(); ++k)
    std::fprintf(f, "%.9e %.9e 0\n", g.Er[k], g.Ez[k]);
  std::fclose(f);
}

void write_shape_csv(const MeniscusShape& s, const std::string& path,
                     const std::string& header) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  if (!header.empty()) std::fputs(header.c_str(), f);
  std::fprintf(f, "i,r,z\n");
  for (std::size_t i = 0; i < s.nodes.size(); ++i)
    std::fprintf(f, "%zu,%.9e,%.9e\n", i, s.nodes[i].r, s.nodes[i].z);
  std::fclose(f);
}

// ---------------------------------------------------------------------------

std::string output_path(const std::string& prefix, const std::string& app,
                        const std::string& state, Real voltage, const std::string& what,
                        const std::string& ext) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.1f", voltage);
  std::string v(buf);
  for (char& c : v)
    if (c == '.') c = 'p';  // keep the name shell- and glob-friendly
  return prefix + "_" + app + "_" + state + "_U" + v + "V_" + what + "." + ext;
}

std::string meta_header(const std::string& app, const std::string& state, Real voltage,
                        const std::string& note) {
  char buf[1024];
  std::snprintf(buf, sizeof buf,
                "# application : %s\n"
                "# state       : %s\n"
                "# voltage     : %.6g V (emitter minus extractor)\n"
                "# note        : %s\n"
                "# WARNUNG     : Prototyp, keine validierte Simulation. Siehe docs/.\n",
                app.c_str(), state.c_str(), voltage,
                note.empty() ? "-" : note.c_str());
  return std::string(buf);
}

void Setup::print(std::FILE* out) const {
  std::fprintf(out, "geometry\n");
  std::fprintf(out, "  emitter type        : %s\n", emitter_type.c_str());
  std::fprintf(out, "  elements            : %zu\n", electrodes.elems.size());
  if (gap > 0) std::fprintf(out, "  tip-extractor gap   : %10.4g m\n", gap);
  if (r_contact > 0) std::fprintf(out, "  contact radius      : %10.4g m\n", r_contact);
  std::fprintf(out, "  applied voltage     : %10.2f V\n", voltage);
  std::fprintf(out, "\n");
  fluid.print(out);
}

Setup build_setup(const Config& cfg) {
  Setup s;
  s.emitter_type = cfg.str("emitter.type", "capillary");
  s.temperature = cfg.num("fluid.temperature", 298.15);
  s.voltage = cfg.num("solve.voltage", 1500.0);

  // --- fluid ---------------------------------------------------------------
  s.fluid = fluid_by_name(cfg.str("fluid.name", "EMI-BF4"));
  s.fluid = s.fluid.at_temperature(s.temperature);
  // Explicit overrides win over the table and over the temperature model.
  s.fluid.rho = cfg.num("fluid.rho", s.fluid.rho);
  s.fluid.gamma = cfg.num("fluid.gamma", s.fluid.gamma);
  s.fluid.K = cfg.num("fluid.conductivity", s.fluid.K);
  s.fluid.mu = cfg.num("fluid.viscosity", s.fluid.mu);
  s.fluid.eps_r = cfg.num("fluid.eps_r", s.fluid.eps_r);
  s.fluid.dG_solvation = cfg.num("fluid.dg_ev", s.fluid.dG_solvation / constants::eV) * constants::eV;
  s.fluid.mean_solvation_n = cfg.num("fluid.solvation_n", s.fluid.mean_solvation_n);
  s.fluid.evap_prefactor = cfg.num("fluid.evap_prefactor", s.fluid.evap_prefactor);

  // --- emitter -------------------------------------------------------------
  std::vector<Mesh> parts;
  Real z_tip = 0.0;
  if (s.emitter_type == "capillary" || s.emitter_type == "capillary_open") {
    OpenCapillaryParams p;
    p.r_bore = cfg.num("emitter.r_bore", 1.0e-5);
    p.r_outer = cfg.num("emitter.r_outer", 2.0 * p.r_bore);
    p.shank_length = cfg.num("emitter.shank_length", 1.0e-3);
    p.z_rim = cfg.num("emitter.z_tip", 0.0);
    p.h_rim = cfg.num("emitter.h_tip", 0.0);
    p.h_far = cfg.num("emitter.h_far", 0.0);
    parts.push_back(make_capillary_open(p));
    s.r_contact = p.r_bore;
    s.wetted = true;
    z_tip = p.z_rim;
  } else if (s.emitter_type == "capillary_dry") {
    CapillaryParams p;
    p.r_inner = cfg.num("emitter.r_bore", 1.0e-5);
    p.r_outer = cfg.num("emitter.r_outer", 2.0 * p.r_inner);
    p.rim_radius = cfg.num("emitter.rim_radius", 0.0);
    p.shank_length = cfg.num("emitter.shank_length", 1.0e-3);
    p.z_tip = cfg.num("emitter.z_tip", 0.0);
    p.h_tip = cfg.num("emitter.h_tip", 0.0);
    p.h_far = cfg.num("emitter.h_far", 0.0);
    parts.push_back(make_capillary(p));
    s.r_contact = p.r_inner;
    z_tip = p.z_tip;
  } else if (s.emitter_type == "needle") {
    NeedleParams p;
    p.tip_radius = cfg.num("emitter.tip_radius", 2.0e-6);
    p.half_angle = cfg.num("emitter.half_angle", 15.0 * pi / 180.0);
    p.shank_radius = cfg.num("emitter.shank_radius", 1.5e-4);
    p.length = cfg.num("emitter.length", 1.0e-3);
    p.z_tip = cfg.num("emitter.z_tip", 0.0);
    p.h_tip = cfg.num("emitter.h_tip", 0.0);
    p.h_far = cfg.num("emitter.h_far", 0.0);
    parts.push_back(make_needle(p));
    // An externally wetted tip has no pinned edge; the tip radius is the best
    // available stand-in for the contact radius (see meniscus.hpp).
    s.r_contact = p.tip_radius;
    z_tip = p.z_tip;
  } else if (s.emitter_type == "sphere") {
    parts.push_back(make_sphere(cfg.num("emitter.radius", 1.0e-3), 0.0,
                                cfg.integer("emitter.n_elem", 200)));
  } else if (s.emitter_type == "spheroid") {
    parts.push_back(make_prolate_spheroid(cfg.num("emitter.a", 2.0e-4),
                                          cfg.num("emitter.b", 1.0e-5), 0.0,
                                          cfg.integer("emitter.n_elem", 400)));
  } else {
    throw std::runtime_error("unknown emitter.type '" + s.emitter_type +
                             "'; expected capillary, capillary_dry, needle, sphere or spheroid");
  }

  // --- extractor -----------------------------------------------------------
  if (cfg.flag("extractor.enabled", true)) {
    ExtractorParams p;
    p.aperture_radius = cfg.num("extractor.aperture_radius", 2.0e-4);
    p.outer_radius = cfg.num("extractor.outer_radius", 3.0e-3);
    p.thickness = cfg.num("extractor.thickness", 1.0e-4);
    s.gap = cfg.num("extractor.gap", 5.0e-4);
    p.z_plate = z_tip + s.gap;
    p.edge_radius = cfg.num("extractor.edge_radius", 0.0);
    p.h_edge = cfg.num("extractor.h_edge", 0.0);
    p.h_far = cfg.num("extractor.h_far", 0.0);
    parts.push_back(make_extractor(p));
  }

  s.electrodes = merge(parts);
  return s;
}

MeniscusParams meniscus_params_from(const Config& cfg, const Setup& s) {
  MeniscusParams p;
  p.r_contact = cfg.num("meniscus.r_contact", s.r_contact);
  p.z_contact = cfg.num("emitter.z_tip", 0.0);
  p.gamma = s.fluid.gamma;
  p.delta_p = cfg.num("meniscus.delta_p", 0.0);
  p.rho = cfg.flag("meniscus.gravity", false) ? s.fluid.rho : 0.0;
  p.n_nodes = cfg.integer("meniscus.n_nodes", 81);
  p.apex_clustering = cfg.num("meniscus.apex_clustering", 1.8);
  p.max_outer = cfg.integer("meniscus.max_outer", 40);
  p.relax = cfg.num("meniscus.relax", 0.5);
  p.tol = cfg.num("meniscus.tol", 3.0e-4);
  p.verbose = cfg.flag("meniscus.verbose", false);
  if (!(p.r_contact > 0.0))
    throw std::runtime_error("meniscus.r_contact is zero: this emitter type has no bore, "
                             "set it explicitly");
  return p;
}

void print_key_reference(std::FILE* out) {
  std::fprintf(out, R"(configuration keys (SI unless a unit suffix is given; see examples/*.cfg)

[fluid]
  name              EMI-BF4 | EMI-Im | BMI-BF4 | formamide+NaI | glycerol+NaI
  temperature       operating temperature, e.g. 25C or 298.15K
  rho gamma conductivity viscosity eps_r    override the table
  dg_ev             ion-evaporation activation energy in eV (primary fit knob)
  solvation_n       mean cluster solvation number (sets q/m of the ion beam)
  evap_prefactor    multiplier on the Iribarne-Thomson prefactor

[emitter]
  type              capillary | capillary_dry | needle | sphere | spheroid
  r_bore r_outer shank_length z_tip rim_radius      (capillary)
  tip_radius half_angle shank_radius length         (needle)
  radius | a b n_elem                               (sphere / spheroid)
  h_tip h_far       element sizes; 0 selects the builder default

[extractor]
  enabled           true | false
  gap               emitter tip to extractor face, e.g. 500um
  aperture_radius outer_radius thickness edge_radius h_edge h_far

[meniscus]
  r_contact delta_p n_nodes apex_clustering max_outer relax tol gravity verbose
  branch            lower | upper | none -- which meniscus is meant when a
                    voltage admits more than one.  Refers to APEX HEIGHT only,
                    not to stability.  "none" makes an ambiguous request fail.

[solve]
  voltage           emitter-to-extractor voltage, e.g. 1.5kV

[beam]
  z_end r_max cfl max_steps path_samples space_charge_iters
  species           ion | droplet | both
  droplet_qm droplet_fraction

[output]
  prefix            file-name prefix for the CSV/VTK dumps
  grid              true | false
  grid_nr grid_nz grid_rmax grid_zmin grid_zmax
)");
}

}  // namespace es
