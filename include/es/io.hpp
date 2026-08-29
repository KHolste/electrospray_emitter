#pragma once
#include <string>
#include <vector>

#include "es/bem.hpp"
#include "es/config.hpp"
#include "es/fluid.hpp"
#include "es/geometry.hpp"
#include "es/meniscus.hpp"
#include "es/types.hpp"

namespace es {

// ---------------------------------------------------------------------------
// Field sampling and output
// ---------------------------------------------------------------------------

struct FieldGrid {
  int nr{0}, nz{0};
  Real r0{0}, r1{0}, z0{0}, z1{0};
  std::vector<Real> V, Er, Ez;  ///< row-major, index = iz * nr + ir

  Real r_at(int ir) const { return r0 + (r1 - r0) * ir / std::max(1, nr - 1); }
  Real z_at(int iz) const { return z0 + (z1 - z0) * iz / std::max(1, nz - 1); }
};

/// Sample potential and field on a uniform (r,z) grid.  Points that fall inside
/// a conductor are meaningless for a boundary-integral solution (the formula
/// evaluates the exterior field everywhere), so they are marked with NaN.
FieldGrid sample_field(const BemSolver& bem, int nr, int nz, Real r0, Real r1, Real z0, Real z1);

void write_grid_csv(const FieldGrid& g, const std::string& path);
/// Legacy ASCII VTK structured points -- opens directly in ParaView.
void write_grid_vtk(const FieldGrid& g, const std::string& path);
void write_shape_csv(const MeniscusShape& s, const std::string& path);

// ---------------------------------------------------------------------------
// Building a run from a config file
// ---------------------------------------------------------------------------

struct Setup {
  Mesh electrodes;        ///< everything except the free surface
  Fluid fluid;
  Real temperature{298.15};
  Real voltage{0};        ///< emitter-to-extractor voltage [V]
  Real gap{0};            ///< emitter tip to extractor face [m]
  Real r_contact{0};      ///< pinning radius, if the emitter has a bore
  bool wetted{false};     ///< emitter geometry leaves room for a meniscus
  std::string emitter_type{"capillary"};

  void print(std::FILE* out) const;
};

/// Assemble geometry and fluid from a Config.  Recognised keys are documented
/// in examples/*.cfg and listed by `--help` in each application.
Setup build_setup(const Config& cfg);

/// Fill MeniscusParams from the config and the setup.
MeniscusParams meniscus_params_from(const Config& cfg, const Setup& s);

void print_key_reference(std::FILE* out);

}  // namespace es
