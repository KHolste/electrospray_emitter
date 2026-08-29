#include <cstdio>
#include <fstream>

#include "es/config.hpp"
#include "es/constants.hpp"
#include "es/io.hpp"

using namespace es;
using constants::pi;

static int failures = 0;

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-44s got=%-15.7g want=%-15.7g %s\n", what, got, want, ok ? "OK" : "FAIL");
}

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-44s %s\n", what, ok ? "OK" : "FAIL");
}

static void test_units() {
  std::printf("\n=== unit suffixes ===\n");
  check("bare number is SI", Config::parse_value("1.5"), 1.5, 1e-15);
  check("10um", Config::parse_value("10um"), 1e-5, 1e-15);
  check("500nm", Config::parse_value("500nm"), 5e-7, 1e-15);
  check("2mm", Config::parse_value("2mm"), 2e-3, 1e-15);
  check("1.5kV", Config::parse_value("1.5kV"), 1500.0, 1e-15);
  check("80nA", Config::parse_value("80nA"), 8e-8, 1e-15);
  check("0.5nL/s", Config::parse_value("0.5nL/s"), 5e-13, 1e-15);
  check("49.3deg", Config::parse_value("49.3deg"), 49.3 * pi / 180.0, 1e-15);
  check("25C", Config::parse_value("25C"), 298.15, 1e-15);
  check("1bar", Config::parse_value("1bar"), 1e5, 1e-15);
  check("scientific 3.2e-5", Config::parse_value("3.2e-5"), 3.2e-5, 1e-15);
  check("negative -50Pa", Config::parse_value("-50Pa"), -50.0, 1e-15);
  check("space before unit: 10 um", Config::parse_value("10 um"), 1e-5, 1e-15);

  // "mm" must not be swallowed by "m", and an unknown unit must be an error,
  // not a silent zero.
  expect("mm is not parsed as m", Config::parse_value("1mm") != Config::parse_value("1m"));
  bool threw = false;
  try { Config::parse_value("10furlong"); } catch (const std::exception&) { threw = true; }
  expect("unknown unit is rejected", threw);
  threw = false;
  try { Config::parse_value("abc"); } catch (const std::exception&) { threw = true; }
  expect("non-numeric value is rejected", threw);
}

static void test_file_and_overrides() {
  std::printf("\n=== file parsing, sections and override order ===\n");
  const char* path = "test_config_tmp.cfg";
  {
    std::ofstream f(path);
    f << "# comment line\n";
    f << "[emitter]\n";
    f << "  type = capillary   ; trailing comment\n";
    f << "  r_bore = 10um\n";
    f << "\n";
    f << "[solve]\n";
    f << "voltage = 1.2kV\n";
    f << "emitter.r_outer = 25um   # already dotted, section must not double-prefix\n";
  }

  Config cfg;
  cfg.load(path);
  check("section prefixing: emitter.r_bore", cfg.num("emitter.r_bore", 0.0), 1e-5, 1e-15);
  check("dotted key inside a section", cfg.num("emitter.r_outer", 0.0), 2.5e-5, 1e-15);
  check("solve.voltage", cfg.num("solve.voltage", 0.0), 1200.0, 1e-15);
  expect("comments are stripped", cfg.str("emitter.type", "") == "capillary");
  check("missing key falls back", cfg.num("emitter.nonexistent", 42.0), 42.0, 1e-15);

  // Command line must win over the file.
  const char* argv[] = {"prog", path, "solve.voltage=900", "emitter.r_bore=5um"};
  Config cfg2;
  for (const std::string& a : Config::positional_args(4, const_cast<char**>(argv))) cfg2.load(a);
  cfg2.apply_cli(4, const_cast<char**>(argv));
  check("CLI overrides the file (voltage)", cfg2.num("solve.voltage", 0.0), 900.0, 1e-15);
  check("CLI overrides the file (r_bore)", cfg2.num("emitter.r_bore", 0.0), 5e-6, 1e-15);

  // Typo detection.
  Config cfg3;
  cfg3.set("emitter.r_bore", "10um");
  cfg3.set("emitter.tipradius", "2um");  // typo: should be tip_radius
  cfg3.set("beam.z_end", "2mm");
  (void)cfg3.num("emitter.r_bore", 0.0);
  const std::vector<std::string> unused = cfg3.unused_keys({"beam."});
  expect("unused key detection finds exactly the typo",
         unused.size() == 1 && unused[0] == "emitter.tipradius");

  std::remove(path);
}

static void test_setup_from_config() {
  std::printf("\n=== Setup construction ===\n");
  Config cfg;
  cfg.set("emitter.type", "capillary");
  cfg.set("emitter.r_bore", "8um");
  cfg.set("extractor.gap", "400um");
  cfg.set("fluid.name", "EMI-Im");
  cfg.set("fluid.temperature", "40C");
  cfg.set("solve.voltage", "1.3kV");
  const Setup s = build_setup(cfg);
  check("contact radius from r_bore", s.r_contact, 8e-6, 1e-15);
  check("gap", s.gap, 4e-4, 1e-15);
  check("voltage", s.voltage, 1300.0, 1e-15);
  check("temperature", s.temperature, 313.15, 1e-15);
  expect("fluid resolved", s.fluid.name == "EMI-Im");
  expect("temperature model applied (K rises above the 25 C table value)",
         s.fluid.K > fluid_by_name("EMI-Im").K);
  expect("capillary geometry is wetted", s.wetted);
  expect("mesh was built", s.electrodes.size() > 20);

  // An explicit fluid override must beat both the table and the T model.
  cfg.set("fluid.conductivity", "0.5");
  const Setup s2 = build_setup(cfg);
  check("explicit conductivity override wins", s2.fluid.K, 0.5, 1e-15);

  bool threw = false;
  try {
    Config bad;
    bad.set("emitter.type", "banana");
    build_setup(bad);
  } catch (const std::exception&) { threw = true; }
  expect("unknown emitter type is rejected", threw);
}

int main() {
  test_units();
  test_file_and_overrides();
  test_setup_from_config();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
