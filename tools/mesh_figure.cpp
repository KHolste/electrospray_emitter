// mesh_figure -- build the P1 device geometry and its automatic boundary mesh
// from a config file, validate the mesh, and write the data behind the mesh
// figures.
//
//   mesh_figure <config.cfg> <output-directory>
//
// Writes the geometry CSVs (regions, boundaries, features, parameters) so the
// figures have their background, plus mesh_nodes.csv, mesh_elements.csv,
// mesh_boundaries.csv, mesh_size_field.csv and mesh_report.txt.  Nothing is
// read back in; python/plot_mesh.py turns these into PNGs.
//
// Exit code 2 means the mesh failed one of its own checks.

#include <cstdio>
#include <stdexcept>
#include <string>

#include "es/boundary_mesh.hpp"
#include "es/config.hpp"
#include "es/device_geometry.hpp"

using namespace es;

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
  p.reserved.edge_radius_inner = c.num("reserved.edge_radius_inner", 0.0);
  p.reserved.edge_radius_outer = c.num("reserved.edge_radius_outer", 0.0);
  p.reserved.contact_angle_deg = c.num("reserved.contact_angle_deg", 0.0);
  p.reserved.bore_diameter_at_inlet = c.num("reserved.bore_diameter_at_inlet", 0.0);
  p.reserved.porous_emitter = c.flag("reserved.porous_emitter", false);
  p.reserved.collector_enabled = c.flag("reserved.collector_enabled", false);
  return p;
}

const char* kPreamble =
    "Automatischer achsensymmetrischer Randvernetzer (P1)\n"
    "====================================================\n\n"
    "Erzeugt aus %s durch tools/mesh_figure.cpp.\n"
    "Meridianhalbebene (r, z), SI-Einheiten. z = 0 ist die Stirnebene des\n"
    "Emitters, z waechst zum Extraktor.\n\n"
    "NUR der Rand ist diskretisiert -- keine Volumenvernetzung, keine\n"
    "Feldloesung, kein Meniskus, keine Stroemung, keine Emission, keine\n"
    "Raumladung. Die Randkennungen sind Bezeichner, keine Randbedingungen.\n"
    "Die P0-BEM benutzt dieses Netz nicht; sie behaelt ihre eigenen Netze.\n\n"
    "Die Flaeche bei z = 0 ist die anfaengliche ebene Fluessigkeitsoberflaeche\n"
    "-- noch kein berechneter Meniskus. Sie ist der geometrische Ausgangszustand\n"
    "und traegt keine physikalische Aussage.\n\n"
    "Scharfe Kanten (Austrittskante, Aperturkanten) werden bewusst verfeinert.\n"
    "Ein dort spaeter berechnetes maximales elektrisches Feld ist deshalb KEINE\n"
    "netzkonvergente Groesse: das Eckfeld einer unverrundeten Kante divergiert\n"
    "und folgt der oertlichen Elementgroesse.\n\n";

}  // namespace

int main(int argc, char** argv) try {
  if (argc < 3) {
    std::fprintf(stderr, "usage: mesh_figure <config.cfg> <output-directory>\n");
    return 1;
  }
  Config cfg;
  cfg.load(argv[1]);
  const std::string dir = argv[2];

  const DeviceGeometry g = DeviceGeometry::build(from_config(cfg));
  g.write_csv(dir);

  const BoundaryMesh mesh = BoundaryMesh::generate(g);
  mesh.write_csv(dir, g);
  mesh.print(stdout, g);

  std::printf("\nPruefungen des Randnetzes\n");
  const MeshReport rep = mesh.validate(g);
  rep.print(stdout);

  {
    std::FILE* f = std::fopen((dir + "/mesh_report.txt").c_str(), "w");
    if (!f) throw std::runtime_error("cannot open mesh_report.txt");
    std::fprintf(f, kPreamble, argv[1]);
    mesh.print(f, g);
    std::fprintf(f, "\nPruefungen des Randnetzes\n");
    rep.print(f);
    std::fclose(f);
  }

  std::printf("\ngeschrieben nach %s\n", dir.c_str());
  cfg.warn_about_unused(stdout, {"reserved."});
  return rep.all_passed() ? 0 : 2;
} catch (const NotImplementedInThisPhase& e) {
  std::fprintf(stderr, "\nmesh_figure: %s\n", e.what());
  return 3;
} catch (const std::exception& e) {
  std::fprintf(stderr, "mesh_figure: %s\n", e.what());
  return 1;
}
