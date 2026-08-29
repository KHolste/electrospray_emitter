// es_field -- electrostatics of an emitter/extractor pair.
//
// Solves the axisymmetric Laplace problem for the given geometry, reports the
// tip field enhancement and the onset estimates, and optionally sweeps the
// voltage or dumps the field on a grid.

#include <cstdio>
#include <stdexcept>

#include "es/bem.hpp"
#include "es/constants.hpp"
#include "es/emission.hpp"
#include "es/io.hpp"

using namespace es;

int main(int argc, char** argv) try {
  Config cfg;
  // Files first, command line second: an override on the command line must win.
  const std::vector<std::string> rest = Config::positional_args(argc, argv);
  for (const std::string& a : rest) {
    if (a == "--help" || a == "-h") {
      std::printf("usage: es_field [file.cfg] [key=value ...]\n\n");
      print_key_reference(stdout);
      return 0;
    }
    cfg.load(a);
  }
  cfg.apply_cli(argc, argv);

  Setup s = build_setup(cfg);
  BemSolver bem(s.electrodes);
  std::printf("assembling %d boundary elements ...\n", static_cast<int>(bem.size()));
  bem.solve_basis();
  bem.solve({s.voltage, 0.0, 0.0});

  s.print(stdout);

  Index which = -1;
  const Real Epeak = bem.peak_emitter_field(&which);
  std::printf("\nelectrostatics\n");
  std::printf("  peak |E_n| on emitter : %10.4g V/m  (= %.4f V/nm)\n", Epeak, Epeak * 1e-9);
  if (which >= 0) {
    const Element& e = bem.mesh().elems[static_cast<std::size_t>(which)];
    std::printf("  at (r,z)              : (%.4g, %.4g) m\n", e.mid.r, e.mid.z);
  }
  if (s.gap > 0.0)
    std::printf("  enhancement E/(U/gap) : %10.2f\n", Epeak / (s.voltage / s.gap));
  std::printf("  charge on emitter     : %10.4g C\n", bem.charge_on(Tag::Emitter));
  std::printf("  charge on extractor   : %10.4g C\n", bem.charge_on(Tag::Extractor));
  std::printf("  net charge            : %10.4g C\n", bem.total_charge());

  // Where does this field sit relative to the two thresholds that matter?
  const Real E_evap = characteristic_evaporation_field(s.fluid);
  std::printf("\nthresholds for %s\n", s.fluid.name.c_str());
  std::printf("  field for G(E) = dG   : %10.4g V/m  -> ratio %.3f\n", E_evap, Epeak / E_evap);
  if (s.r_contact > 0.0) {
    const Real E_hemi = hemisphere_balance_field(s.r_contact, s.fluid.gamma);
    std::printf("  hemispherical onset   : %10.4g V/m  -> ratio %.3f\n", E_hemi, Epeak / E_hemi);
    if (s.gap > 0.0)
      std::printf("  Taylor/Smith V_onset  : %10.1f V   (applied %.1f V)\n",
                  literature_onset_voltage_smith(s.r_contact, s.gap, s.fluid.gamma), s.voltage);
  }
  std::printf("\nNOTE: this is the DRY field of the metal geometry.  A wetted emitter\n"
              "      forms a meniscus that sharpens under the field and raises the apex\n"
              "      field far above this value -- run es_meniscus for that.\n");

  const std::string prefix = cfg.str("output.prefix", "field");
  bem.write_surface_csv(prefix + "_surface.csv");
  bem.mesh().write_csv(prefix + "_mesh.csv");
  std::printf("\nwrote %s_surface.csv, %s_mesh.csv\n", prefix.c_str(), prefix.c_str());

  // --- optional voltage sweep ----------------------------------------------
  if (cfg.has("sweep.voltage_max")) {
    const Real v0 = cfg.num("sweep.voltage_min", 0.0);
    const Real v1 = cfg.num("sweep.voltage_max", 2000.0);
    const int n = cfg.integer("sweep.steps", 21);
    const std::string path = prefix + "_sweep.csv";
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot open " + path);
    std::fprintf(f, "voltage,peak_field,emitter_charge,ion_current\n");
    std::printf("\nvoltage sweep (the field is linear in U; only the emission is not)\n");
    std::printf("  %10s %14s %14s\n", "U [V]", "E_peak [V/m]", "I_ion [A]");
    for (int i = 0; i < n; ++i) {
      const Real U = v0 + (v1 - v0) * i / std::max(1, n - 1);
      bem.solve({U, 0.0, 0.0});
      IonEmission ie;
      bool ion_available = true;
      try {
        ie = integrate_ion_emission(bem, s.fluid, s.temperature, true);
      } catch (const NotImplementedInThisPhase& ex) {
        // Do not report a silent zero: say that the quantity is unavailable.
        ion_available = false;
        if (i == 0) std::printf("\n  I_ion nicht verfuegbar: %s\n", ex.what());
      }
      if (ion_available) {
        std::fprintf(f, "%.9e,%.9e,%.9e,%.9e\n", U, bem.peak_emitter_field(),
                     bem.charge_on(Tag::Emitter), ie.current);
        std::printf("  %10.1f %14.4g %14.4g\n", U, bem.peak_emitter_field(), ie.current);
      } else {
        std::fprintf(f, "%.9e,%.9e,%.9e,nan\n", U, bem.peak_emitter_field(),
                     bem.charge_on(Tag::Emitter));
        std::printf("  %10.1f %14.4g %14s\n", U, bem.peak_emitter_field(), "n/a");
      }
    }
    std::fclose(f);
    std::printf("wrote %s\n", path.c_str());
    bem.solve({s.voltage, 0.0, 0.0});
  }

  // --- optional field grid -------------------------------------------------
  if (cfg.flag("output.grid", false)) {
    const Real rmax = cfg.num("output.grid_rmax", 5.0 * std::max(s.gap, 1e-4));
    const Real zmin = cfg.num("output.grid_zmin", -2.0 * std::max(s.gap, 1e-4));
    const Real zmax = cfg.num("output.grid_zmax", 3.0 * std::max(s.gap, 1e-4));
    const int nr = cfg.integer("output.grid_nr", 120);
    const int nz = cfg.integer("output.grid_nz", 200);
    std::printf("\nsampling %d x %d field grid ...\n", nr, nz);
    const FieldGrid g = sample_field(bem, nr, nz, 0.0, rmax, zmin, zmax);
    write_grid_csv(g, prefix + "_grid.csv");
    write_grid_vtk(g, prefix + "_grid.vtk");
    std::printf("wrote %s_grid.csv and %s_grid.vtk\n", prefix.c_str(), prefix.c_str());
  }

  cfg.warn_about_unused(stdout, {"beam.", "meniscus.", "operate."});
  return 0;
} catch (const std::exception& e) {
  std::fprintf(stderr, "es_field: %s\n", e.what());
  return 1;
}
