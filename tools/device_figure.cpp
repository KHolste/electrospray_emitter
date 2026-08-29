// device_figure -- build the P1 device geometry from a config file and write
// the data behind the geometry figures.
//
//   device_figure <config.cfg> <output-directory>
//
// Writes regions.csv, boundaries.csv, features.csv and parameters.csv, plus a
// meta.txt with the resolved parameter set and the derived measures.  Nothing
// is read back in; python/plot_device.py turns these into PNGs.

#include <cstdio>
#include <stdexcept>
#include <string>

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/device_geometry.hpp"

using namespace es;
using constants::pi;

namespace {

DeviceParameters from_config(const Config& c) {
  DeviceParameters p;
  p.phi_3 = c.num("device.phi_3", p.phi_3);
  p.phi_1 = c.num("device.phi_1", p.phi_1);
  p.phi_2 = c.num("device.phi_2", p.phi_2);
  p.emitter_height = c.num("device.emitter_height", p.emitter_height);
  p.extraction_distance = c.num("device.extraction_distance", p.extraction_distance);
  p.extractor_aperture_diameter =
      c.num("device.extractor_aperture_diameter", p.extractor_aperture_diameter);
  p.extractor_thickness = c.num("device.extractor_thickness", p.extractor_thickness);
  p.domain_radius = c.num("domain.radius", p.domain_radius);
  p.domain_z_min = c.num("domain.z_min", p.domain_z_min);
  p.domain_z_max = c.num("domain.z_max", p.domain_z_max);
  p.extractor_outer_radius = c.num("device.extractor_outer_radius", p.extractor_outer_radius);
  // Reserved keys are read so that a typo is not silently ignored; any value
  // other than the default makes build() refuse.
  p.reserved.edge_radius_inner = c.num("reserved.edge_radius_inner", 0.0);
  p.reserved.edge_radius_outer = c.num("reserved.edge_radius_outer", 0.0);
  p.reserved.contact_angle_deg = c.num("reserved.contact_angle_deg", 0.0);
  p.reserved.bore_diameter_at_inlet = c.num("reserved.bore_diameter_at_inlet", 0.0);
  p.reserved.porous_emitter = c.flag("reserved.porous_emitter", false);
  p.reserved.collector_enabled = c.flag("reserved.collector_enabled", false);
  return p;
}

}  // namespace

int main(int argc, char** argv) try {
  if (argc < 3) {
    std::fprintf(stderr, "usage: device_figure <config.cfg> <output-directory>\n");
    return 1;
  }
  Config cfg;
  cfg.load(argv[1]);
  const std::string dir = argv[2];

  const DeviceGeometry g = DeviceGeometry::build(from_config(cfg));
  g.print(stdout);
  g.write_csv(dir);

  {
    std::FILE* f = std::fopen((dir + "/meta.txt").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open meta.txt");
    std::fprintf(f, "parametrische P1-Geometrie\n"
                    "==========================\n\n"
                    "Erzeugt aus %s durch tools/device_figure.cpp.\n"
                    "Achsensymmetrisch in der Meridianhalbebene (r, z), SI-Einheiten.\n"
                    "z = 0 ist die Stirnebene des Emitters, z waechst zum Extraktor.\n\n"
                    "In dieser Phase wird KEINE Physik geloest: die Randkennungen sind\n"
                    "Bezeichner, keine Randbedingungen. Die Flaeche 'free_surface_reference'\n"
                    "ist die ebene Referenzflaeche am Bohrungsaustritt, nicht ein gerechneter\n"
                    "Meniskus.\n\n", argv[1]);
    g.print(f);
    std::fclose(f);
  }

  std::printf("\ngeschrieben nach %s\n", dir.c_str());
  cfg.warn_about_unused(stdout, {"reserved."});
  return 0;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\ndevice_figure: %s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "device_figure: %s\n", e.what());
  return 1;
}
