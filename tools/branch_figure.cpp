// branch_figure -- generate the data behind the branch-ambiguity figures.
//
//   branch_figure <output-directory>
//
// Recomputes everything from scratch and writes CSV files that
// python/plot_branch.py turns into two PNGs.  No previously computed file is
// read, so the figures can be reproduced from source at any time.
//
// This is a diagnostic tool, not part of the physics.  It exists because the
// question "how many menisci belong to one voltage" is easier to settle by
// looking at the branch than by reading a status code.

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "es/constants.hpp"
#include "es/fluid.hpp"
#include "es/io.hpp"
#include "es/meniscus.hpp"

using namespace es;
using constants::pi;

namespace {

/// Target voltage for the figures [V].  Fixed, so the figures are stable.
constexpr Real kTargetVoltage = 1154.2;

Mesh electrodes() {
  OpenCapillaryParams cp;
  cp.r_bore = 1.0e-5;
  cp.r_outer = 2.0e-5;
  cp.shank_length = 1.0e-3;
  cp.z_rim = 0.0;
  cp.h_rim = 1.0e-5 / 14.0;
  ExtractorParams ep;
  ep.aperture_radius = 2.0e-4;
  ep.outer_radius = 3.0e-3;
  ep.thickness = 1.0e-4;
  ep.z_plate = 5.0e-4;
  ep.h_edge = 1.0e-5;
  return merge({make_capillary_open(cp), make_extractor(ep)});
}

MeniscusParams params() {
  MeniscusParams mp;
  mp.r_contact = 1.0e-5;
  mp.z_contact = 0.0;
  mp.gamma = fluid_by_name("EMI-BF4").gamma;
  mp.delta_p = 0.0;
  mp.n_nodes = 61;
  mp.max_outer = 60;
  mp.relax = 0.5;
  mp.tol = 3e-4;
  mp.voltage_tol = 1e-3;
  return mp;
}

void write_branch(const std::string& path, const std::vector<MeniscusSolution>& br, Real r_c,
                  const std::string& label) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# branch segment: %s\n", label.c_str());
  std::fprintf(f, "# r_contact (m) : %.9e\n", r_c);
  std::fprintf(f, "h_over_rc,height_m,voltage_V,apex_field_Vpm,apex_radius_m,status\n");
  for (const MeniscusSolution& m : br)
    std::fprintf(f, "%.9e,%.9e,%.9e,%.9e,%.9e,%s\n", m.shape.height / r_c, m.shape.height,
                 m.voltage, m.apex_field, m.shape.apex_radius, to_string(m.status));
  std::fclose(f);
}

}  // namespace

int main(int argc, char** argv) try {
  if (argc < 2) {
    std::fprintf(stderr, "usage: branch_figure <output-directory>\n");
    return 1;
  }
  const std::string dir = std::string(argv[1]) + "/";
  const Real r_c = params().r_contact;
  const Real U = kTargetVoltage;

  MeniscusSolver s(electrodes(), params());

  // --- full branch ---------------------------------------------------------
  std::printf("tracing the full branch ...\n");
  const std::vector<MeniscusSolution> full = s.continuation(0.10 * r_c, 2.2 * r_c, 22);
  const MeniscusSolver::StaticFold fold = MeniscusSolver::find_static_fold(full);
  if (!fold.found())
    throw std::runtime_error(std::string("no fold candidate: ") + to_string(fold.status));
  write_branch(dir + "branch_full.csv", full, r_c, "full traced range, h_max = 2.2 r_c");

  // --- two deliberately short ranges, for the coverage illustration --------
  const Real h_before = 0.90 * fold.height;
  const Real h_past = 1.25 * fold.height;
  write_branch(dir + "branch_before_fold.csv",
               s.continuation(0.10 * r_c, h_before, 12), r_c,
               "range ends before the turning point");
  write_branch(dir + "branch_past_fold.csv",
               s.continuation(0.10 * r_c, h_past, 12), r_c,
               "range ends past the turning point but above the target voltage");

  // --- the two solutions at the target voltage -----------------------------
  std::printf("solving both crossings at U = %.1f V ...\n", U);
  MeniscusSolution lo = s.solve_at_voltage(U, 2.2 * r_c, BranchSide::LowerHeight, 22);
  MeniscusSolution hi = s.solve_at_voltage(U, 2.2 * r_c, BranchSide::UpperHeight, 22);
  if (!lo.ok() || !hi.ok())
    throw std::runtime_error("could not bracket both crossings: lower=" +
                             std::string(to_string(lo.status)) + ", upper=" +
                             std::string(to_string(hi.status)));

  write_shape_csv(lo.shape, dir + "shape_lower_height.csv",
                  meta_header("branch_figure", "LowerHeight solution at the target voltage",
                              lo.voltage, "apex height only -- no stability claim"));
  write_shape_csv(hi.shape, dir + "shape_upper_height.csv",
                  meta_header("branch_figure", "UpperHeight solution at the target voltage",
                              hi.voltage, "apex height only -- no stability claim"));

  // --- geometry for context ------------------------------------------------
  s.realize(lo);
  s.bem().mesh().write_csv(dir + "mesh_electrodes.csv",
                           meta_header("branch_figure", "emitter and extractor geometry",
                                       lo.voltage, "meniscus in this mesh is the lower one"));

  // --- markers -------------------------------------------------------------
  {
    std::FILE* f = std::fopen((dir + "markers.csv").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open markers.csv");
    std::fprintf(f, "# points to be marked in the U(h) figure\n");
    std::fprintf(f, "name,h_over_rc,height_m,voltage_V,apex_radius_m,apex_field_Vpm\n");
    std::fprintf(f, "target_voltage,,,%.9e,,\n", U);
    std::fprintf(f, "fold_candidate,%.9e,%.9e,%.9e,%.9e,%.9e\n", fold.height / r_c, fold.height,
                 fold.voltage, fold.apex_radius, fold.apex_field);
    std::fprintf(f, "crossing_lower_height,%.9e,%.9e,%.9e,%.9e,%.9e\n", lo.shape.height / r_c,
                 lo.shape.height, lo.voltage, lo.shape.apex_radius, lo.apex_field);
    std::fprintf(f, "crossing_upper_height,%.9e,%.9e,%.9e,%.9e,%.9e\n", hi.shape.height / r_c,
                 hi.shape.height, hi.voltage, hi.shape.apex_radius, hi.apex_field);
    std::fprintf(f, "range_end_before_fold,%.9e,%.9e,,,\n", h_before / r_c, h_before);
    std::fprintf(f, "range_end_past_fold,%.9e,%.9e,,,\n", h_past / r_c, h_past);
    std::fprintf(f, "range_end_full,%.9e,%.9e,,,\n", 2.2, 2.2 * r_c);
    std::fclose(f);
  }

  // --- provenance ----------------------------------------------------------
  {
    std::FILE* f = std::fopen((dir + "meta.txt").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open meta.txt");
    std::fprintf(f,
        "branch_figure -- data for the branch-ambiguity figures\n"
        "=====================================================\n\n"
        "geometry   : capillary, bore radius %.4g m, outer radius %.4g m,\n"
        "             extractor gap %.4g m, aperture radius %.4g m\n"
        "fluid      : EMI-BF4, gamma = %.4f N/m, feed pressure %.4g Pa\n"
        "meniscus   : %d free-surface nodes, tol = %.1e, voltage_tol = %.1e\n\n"
        "target voltage        : %.4f V\n"
        "fold candidate        : %.4f V at h/r_c = %.6f\n"
        "  NOTE: a discrete maximum of the sampled branch, refined parabolically.\n"
        "  It is a CANDIDATE for a static turning point.  No dynamic stability,\n"
        "  no emission onset and no cone-jet transition follows from it.\n\n"
        "solution LowerHeight  : h/r_c = %.6f, U = %.4f V, R_apex = %.6g m, E_apex = %.6g V/m\n"
        "solution UpperHeight  : h/r_c = %.6f, U = %.4f V, R_apex = %.6g m, E_apex = %.6g V/m\n"
        "  The names refer to APEX HEIGHT only.  Neither is labelled stable or\n"
        "  unstable: no stability analysis is implemented.\n\n"
        "investigated ranges (h/r_c):\n"
        "  ends before the turning point : %.6f   coverage incomplete\n"
        "  ends past it, still above U   : %.6f   coverage incomplete\n"
        "  full                          : %.6f   coverage complete for this target\n",
        1.0e-5, 2.0e-5, 5.0e-4, 2.0e-4, params().gamma, params().delta_p, params().n_nodes,
        params().tol, params().voltage_tol, U, fold.voltage, fold.height / r_c,
        lo.shape.height / r_c, lo.voltage, lo.shape.apex_radius, lo.apex_field,
        hi.shape.height / r_c, hi.voltage, hi.shape.apex_radius, hi.apex_field,
        h_before / r_c, h_past / r_c, 2.2);
    std::fclose(f);
  }

  std::printf("\nwrote to %s\n", dir.c_str());
  std::printf("  fold candidate : %.2f V at h/r_c = %.4f\n", fold.voltage, fold.height / r_c);
  std::printf("  lower          : h/r_c = %.4f, R_apex/r_c = %.4f\n", lo.shape.height / r_c,
              lo.shape.apex_radius / r_c);
  std::printf("  upper          : h/r_c = %.4f, R_apex/r_c = %.4f\n", hi.shape.height / r_c,
              hi.shape.apex_radius / r_c);
  return 0;
} catch (const std::exception& e) {
  std::fprintf(stderr, "branch_figure: %s\n", e.what());
  return 1;
}
