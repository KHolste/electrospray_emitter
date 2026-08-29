// Regression tests for the eleven confirmed defects of docs/01_gap_analysis.md.
//
// Each test is written against the behaviour that was measured on the prototype,
// so that the specific failure cannot come back unnoticed.  Test names carry the
// finding number.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "es/beam.hpp"
#include "es/bem.hpp"
#include "es/constants.hpp"
#include "es/emission.hpp"
#include "es/fluid.hpp"
#include "es/io.hpp"
#include "es/meniscus.hpp"
#include "es/status.hpp"

using namespace es;
using constants::pi;

static int failures = 0;

static void expect(const char* what, bool ok) {
  if (!ok) ++failures;
  std::printf("  %-62s %s\n", what, ok ? "OK" : "FAIL");
}

static void check(const char* what, Real got, Real want, Real rtol) {
  const Real err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  const bool ok = err <= rtol;
  if (!ok) ++failures;
  std::printf("  %-42s got=%-13.6g want=%-13.6g %s\n", what, got, want, ok ? "OK" : "FAIL");
}

// --------------------------------------------------------------------------
static Mesh standard_electrodes() {
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

static MeniscusParams standard_params() {
  MeniscusParams mp;
  mp.r_contact = 1.0e-5;
  mp.z_contact = 0.0;
  mp.gamma = fluid_by_name("EMI-BF4").gamma;
  mp.delta_p = 0.0;
  mp.n_nodes = 61;
  mp.max_outer = 40;
  mp.relax = 0.5;
  mp.tol = 3e-4;
  mp.voltage_tol = 1e-3;
  return mp;
}

static BemSolver solved_bem(Real U_emitter) {
  BemSolver bem(standard_electrodes());
  bem.solve({U_emitter, 0.0, 0.0});
  return bem;
}

/// Capture what a printing function emits, so that wording requirements can be
/// tested rather than merely intended.
static std::string capture(void (*fn)(std::FILE*)) {
  const char* tmp = "test_regressions_capture.tmp";
  std::FILE* f = std::fopen(tmp, "w+");
  fn(f);
  std::fflush(f);
  std::fseek(f, 0, SEEK_SET);
  std::string out;
  char buf[512];
  while (std::fgets(buf, sizeof buf, f)) out += buf;
  std::fclose(f);
  std::remove(tmp);
  return out;
}

// ==========================================================================
// B1  solve_at_voltage must never report success at a different voltage
// ==========================================================================
static void b1_solve_at_voltage() {
  std::printf("\n=== B1  solve_at_voltage: Status statt stiller Fehlanpassung ===\n");
  MeniscusSolver s(standard_electrodes(), standard_params());
  const Real h_max = 2.0e-5;

  // The prototype answered a 500 V request with a shape belonging to 871.8 V
  // and flagged it converged.
  for (Real U : {500.0, 800.0, 1000.0, 1150.0}) {
    MeniscusSolution m = s.solve_at_voltage(U, h_max);
    const Real tol = standard_params().voltage_tol * std::max(U, 1.0);
    const bool honest = !m.ok() || std::abs(m.voltage - U) <= tol;
    char buf[128];
    std::snprintf(buf, sizeof buf, "U = %7.1f V -> %s, geliefert %.1f V", U,
                  to_string(m.status), m.voltage);
    expect(buf, honest);
    if (!honest)
      std::printf("      Abweichung %.1f V bei Toleranz %.1f V\n", m.voltage - U, tol);
  }

  // Above the fold there is no static solution at all; that must be said, not
  // approximated.
  const std::vector<MeniscusSolution> branch = s.continuation(0.15e-5, h_max, 20);
  const MeniscusSolver::StaticFold fold = MeniscusSolver::find_static_fold(branch);
  expect("Referenzast besitzt einen inneren Umkehrpunkt", fold.found());
  if (fold.found()) {
    MeniscusSolution m = s.solve_at_voltage(fold.voltage * 1.2, h_max);
    expect("Spannung oberhalb der Falte -> VoltageNotBracketed",
           m.status == SolveStatus::VoltageNotBracketed);
    expect("und ok() ist falsch", !m.ok());
  }
}

// ==========================================================================
// B2  a single or monotone branch has no turning point
// ==========================================================================
static MeniscusSolution mk(Real h, Real v, SolveStatus st = SolveStatus::Converged) {
  MeniscusSolution m;
  m.shape.height = h;
  m.voltage = v;
  m.status = st;
  return m;
}

/// The prototype's algorithm, reproduced here so that the tests below can be
/// shown to discriminate: it reported a turning point whenever ANY converged
/// point existed, and took the largest voltage regardless of its position.
static bool prototype_find_onset(const std::vector<MeniscusSolution>& branch, Real* v) {
  std::size_t best = 0;
  bool any = false;
  for (std::size_t i = 0; i < branch.size(); ++i) {
    if (!branch[i].ok()) continue;
    if (!any || branch[i].voltage > branch[best].voltage) { best = i; any = true; }
  }
  if (any && v) *v = branch[best].voltage;
  return any;
}

static void b2_static_fold() {
  std::printf("\n=== B2  find_static_fold: kein Umkehrpunkt aus dem Nichts ===\n");

  // Prototype: a single converged point produced found = true.
  {
    std::vector<MeniscusSolution> b{mk(1e-6, 994.9)};
    const auto f = MeniscusSolver::find_static_fold(b);
    expect("ein einziger Punkt -> TooFewPoints", f.status == FoldStatus::TooFewPoints);
    expect("und found() ist falsch", !f.found());
  }
  {
    std::vector<MeniscusSolution> b{mk(1e-6, 900.0), mk(2e-6, 1000.0)};
    const auto f = MeniscusSolver::find_static_fold(b);
    expect("zwei Punkte -> TooFewPoints", f.status == FoldStatus::TooFewPoints);
  }
  // Prototype: a strictly rising branch produced found = true at the last point.
  {
    std::vector<MeniscusSolution> b;
    for (int i = 0; i < 6; ++i) b.push_back(mk(1e-6 * (i + 1), 500.0 + 100.0 * i));
    const auto f = MeniscusSolver::find_static_fold(b);
    expect("streng monoton steigend -> Monotone", f.status == FoldStatus::Monotone);
    expect("und found() ist falsch", !f.found());
  }
  {
    std::vector<MeniscusSolution> b;
    for (int i = 0; i < 6; ++i) b.push_back(mk(1e-6 * (i + 1), 1000.0 - 100.0 * i));
    const auto f = MeniscusSolver::find_static_fold(b);
    expect("streng monoton fallend -> Monotone", f.status == FoldStatus::Monotone);
  }
  // Show that these tests discriminate: the prototype's algorithm accepts
  // exactly the two cases the new one rejects.
  {
    std::vector<MeniscusSolution> single{mk(1e-6, 994.9)};
    std::vector<MeniscusSolution> rising;
    for (int i = 0; i < 6; ++i) rising.push_back(mk(1e-6 * (i + 1), 500.0 + 100.0 * i));
    Real v1 = 0, v2 = 0;
    const bool p1 = prototype_find_onset(single, &v1);
    const bool p2 = prototype_find_onset(rising, &v2);
    std::printf("  [Gegenprobe] alter Algorithmus: 1 Punkt -> found=%d (%.1f V), "
                "monoton -> found=%d (%.1f V)\n", p1 ? 1 : 0, v1, p2 ? 1 : 0, v2);
    expect("der alte Algorithmus haette beide akzeptiert (Test diskriminiert)", p1 && p2);
    expect("der neue lehnt beide ab",
           !MeniscusSolver::find_static_fold(single).found() &&
               !MeniscusSolver::find_static_fold(rising).found());
  }

  // A genuine interior maximum must be found, and refined.
  {
    std::vector<MeniscusSolution> b{mk(1e-6, 900.0), mk(2e-6, 1000.0), mk(3e-6, 1050.0),
                                    mk(4e-6, 1040.0), mk(5e-6, 980.0)};
    const auto f = MeniscusSolver::find_static_fold(b);
    expect("echter innerer Umkehrpunkt -> Found", f.status == FoldStatus::Found);
    expect("parabolisch verfeinerte Faltenspannung >= Stuetzwert", f.voltage >= 1050.0);
    expect("Faltenhoehe liegt zwischen den Nachbarn",
           f.height > 2e-6 && f.height < 4e-6);
  }
  // Non-converged points must not be counted towards the three.
  {
    std::vector<MeniscusSolution> b{mk(1e-6, 900.0), mk(2e-6, 1000.0, SolveStatus::NotConverged),
                                    mk(3e-6, 800.0, SolveStatus::NotConverged)};
    const auto f = MeniscusSolver::find_static_fold(b);
    expect("nicht konvergierte Punkte zaehlen nicht -> TooFewPoints",
           f.status == FoldStatus::TooFewPoints);
  }
}

// ==========================================================================
// B3  output files identify application, state and voltage; realize() binds
//     surface data to the reported state
// ==========================================================================
static void b3_output_identity() {
  std::printf("\n=== B3  Ausgabedateien: Anwendung, Zustand, Spannung ===\n");
  const std::string a = output_path("out", "meniscus", "fold", 1179.8, "surface");
  const std::string b = output_path("out", "beam", "traced", 1150.0, "surface");
  std::printf("  %s\n  %s\n", a.c_str(), b.c_str());
  expect("es_meniscus und es_beam erzeugen verschiedene Namen", a != b);
  expect("Name nennt die Anwendung", a.find("meniscus") != std::string::npos);
  expect("Name nennt den Zustand", a.find("fold") != std::string::npos);
  expect("Name nennt die Spannung", a.find("1179p8V") != std::string::npos);

  const std::string h = meta_header("es_meniscus", "static_fold", 1179.8, "Testlauf");
  expect("Kopfzeilen sind Kommentare", h.rfind("#", 0) == 0);
  expect("Kopf nennt die Anwendung", h.find("es_meniscus") != std::string::npos);
  expect("Kopf nennt den Zustand", h.find("static_fold") != std::string::npos);
  expect("Kopf nennt die Spannung", h.find("1179.8") != std::string::npos);

  // realize() must put the solver into exactly the state it is handed, so a
  // surface dump provably belongs to it.  The prototype dumped whatever the
  // last internal iteration left behind -- the last NON-converged branch point.
  MeniscusSolver s(standard_electrodes(), standard_params());
  const std::vector<MeniscusSolution> branch = s.continuation(0.15e-5, 2.0e-5, 20);
  const auto fold = MeniscusSolver::find_static_fold(branch);
  if (fold.found()) {
    const MeniscusSolution& at_fold = branch[fold.index];
    // Deliberately disturb the solver state first.
    (void)s.solve_at_height(1.9e-5);
    s.realize(at_fold);
    const Real peak = s.bem().peak_field(Tag::FreeSurface);
    check("Spitzenfeld nach realize() == Zustand am Umkehrpunkt", peak, at_fold.peak_field, 1e-9);
    const Real last_peak = branch.back().peak_field;
    expect("und es ist NICHT der letzte Astpunkt",
           std::abs(peak - last_peak) > 0.1 * std::max(peak, last_peak));
  }
}

// ==========================================================================
// B4 / B5  space charge fails closed
// ==========================================================================
static void b4_space_charge_closed() {
  std::printf("\n=== B4/B5  Raumladung: geschlossenes Fehlschlagen ===\n");
  BemSolver bem = solved_bem(1500.0);
  std::vector<Real> w(static_cast<std::size_t>(bem.size()), 0.0);
  for (Index i = 0; i < bem.size(); ++i)
    if (bem.mesh().elems[static_cast<std::size_t>(i)].tag == Tag::Emitter)
      w[static_cast<std::size_t>(i)] = 1e-9;

  BeamParams bp;
  bp.z_end = 1.0e-3;
  bp.max_steps = 2000;
  bp.space_charge_iters = 1;

  bool threw = false;
  std::string msg;
  try {
    (void)trace_beam_with_weights(bem, w, {{"ion", 5e5, 1.0, SpeciesKind::IonEvaporated}}, bp);
  } catch (const NotImplementedInThisPhase& e) {
    threw = true;
    msg = e.what();
  }
  expect("space_charge_iters > 0 wirft NotImplementedInThisPhase", threw);
  expect("Meldung nennt die vorgesehene Phase", msg.find("P4") != std::string::npos);
  expect("Meldung nennt den Grund", msg.find("Eigenfeld") != std::string::npos ||
                                    msg.find("nicht wohlgestellt") != std::string::npos);

  // The ring self-field is still singular -- that is why the model is closed.
  const Vec2 xp{5e-5, 1e-4};
  const Vec2 E0 = ring_field(1e-15, xp, xp);
  expect("Ring-Eigenfeld bleibt singulaer (Begruendung der Sperre)",
         !std::isfinite(E0.r) || !std::isfinite(E0.z));
}

// ==========================================================================
// B6  droplets must not be launched from the ion evaporation distribution
// ==========================================================================
static void b6_droplet_closed() {
  std::printf("\n=== B6  Tropfenstrahl: gesperrt bis zur Cone-Jet-Kopplung ===\n");
  BemSolver bem = solved_bem(1500.0);
  std::vector<Real> w(static_cast<std::size_t>(bem.size()), 0.0);
  for (Index i = 0; i < bem.size(); ++i)
    if (bem.mesh().elems[static_cast<std::size_t>(i)].tag == Tag::Emitter)
      w[static_cast<std::size_t>(i)] = 1e-9;
  BeamParams bp;
  bp.z_end = 1.0e-3;
  bp.max_steps = 2000;

  bool threw = false;
  std::string msg;
  try {
    (void)trace_beam_with_weights(bem, w, {{"droplet", 1e4, 1.0, SpeciesKind::Droplet}}, bp);
  } catch (const NotImplementedInThisPhase& e) {
    threw = true;
    msg = e.what();
  }
  expect("Droplet-Spezies wirft NotImplementedInThisPhase", threw);
  expect("Meldung nennt Phase P6", msg.find("P6") != std::string::npos);

  // A mixed request must fail too, not silently drop the droplets.
  threw = false;
  try {
    (void)trace_beam_with_weights(
        bem, w,
        {{"ion", 5e5, 0.5, SpeciesKind::IonEvaporated}, {"droplet", 1e4, 0.5, SpeciesKind::Droplet}},
        bp);
  } catch (const NotImplementedInThisPhase&) { threw = true; }
  expect("gemischte Anforderung faellt ebenfalls durch", threw);
}

// ==========================================================================
// B8  negative polarity fails closed
// ==========================================================================
static void b8_polarity_closed() {
  std::printf("\n=== B8  Negative Polaritaet: geschlossenes Fehlschlagen ===\n");
  const Fluid f = fluid_by_name("EMI-BF4");

  BemSolver pos = solved_bem(1500.0);
  bool ok_positive = true;
  try {
    (void)integrate_ion_emission(pos, f, 298.15, true);
  } catch (const NotImplementedInThisPhase&) { ok_positive = false; }
  expect("positive Polaritaet rechnet weiterhin", ok_positive);

  BemSolver neg = solved_bem(-1500.0);
  bool threw = false;
  std::string msg;
  try {
    (void)integrate_ion_emission(neg, f, 298.15, true);
  } catch (const NotImplementedInThisPhase& e) { threw = true; msg = e.what(); }
  expect("negative Polaritaet wirft NotImplementedInThisPhase", threw);
  expect("Meldung nennt Anionen", msg.find("Anion") != std::string::npos);

  // ... and the beam entry point must refuse as well, not just the integral.
  std::vector<Real> w(static_cast<std::size_t>(neg.size()), 1e-12);
  BeamParams bp;
  bp.z_end = 1.0e-3;
  bp.max_steps = 500;
  threw = false;
  try {
    (void)trace_beam(neg, f, 298.15, {{"ion", 5e5, 1.0, SpeciesKind::IonEvaporated}}, bp);
  } catch (const NotImplementedInThisPhase&) { threw = true; }
  expect("auch trace_beam verweigert die negative Polaritaet", threw);
}

// ==========================================================================
// B7 / B11  wording: diagnostic estimate, empirical correlation
// ==========================================================================
static void b7_b11_wording() {
  std::printf("\n=== B7/B11  Kennzeichnung von Abschaetzung und Korrelation ===\n");

  const std::string diag = capture([](std::FILE* f) {
    const Fluid fl = fluid_by_name("EMI-BF4");
    IonEmission ion;
    ion.current = 4.2e-16;
    ion.peak_E = 3.2e7;
    ion.effective_area = 3e-10;
    ion.mdot = 8e-22;
    print_diagnostic_estimate(f, fl, &ion, nullptr);
  });
  expect("Ausgabe sagt 'kein Betriebspunkt'", diag.find("kein Betriebspunkt") != std::string::npos);
  expect("Ausgabe sagt 'NICHT GEKOPPELT'",
         diag.find("NICHT GEKOPPELTE") != std::string::npos);
  expect("Ausgabe sagt, es sei keine Stromvorhersage",
         diag.find("NICHT um eine Stromvorhersage") != std::string::npos);
  expect("Ausgabe nennt Higuera als Grund", diag.find("Higuera") != std::string::npos);
  expect("Ausgabe nennt die vorgesehene Phase P5", diag.find("P5") != std::string::npos);
  expect("Ausgabe gibt die dG-Empfindlichkeit an",
         diag.find("Empfindlichkeit") != std::string::npos);

  const std::string cj = capture([](std::FILE* f) {
    const Fluid fl = fluid_by_name("EMI-BF4");
    print_cone_jet_correlation(f, fl, cone_jet(fl, 1e-13));
  });
  expect("Cone-Jet-Block ist als empirisch markiert",
         cj.find("empirical = true") != std::string::npos);
  expect("Cone-Jet-Block sagt 'nicht gekoppelt'",
         cj.find("nicht an Geometrie") != std::string::npos);
  expect("Cone-Jet-Block warnt bei eps_r < 40", cj.find("ACHTUNG") != std::string::npos);
  expect("ConeJetState traegt das Merkmal empirical", ConeJetState::empirical);
}

// ==========================================================================
// B9  emitting area is a continuous functional
// ==========================================================================
static void b9_effective_area() {
  std::printf("\n=== B9  Emittierende Flaeche: stetiges Funktional ===\n");
  // A sphere at fixed potential has a uniform normal field, hence a uniform
  // current density.  For uniform j the functional (int j dA)^2 / int j^2 dA is
  // exactly the total area -- an exact test the element-quantised "99% area"
  // could never pass.
  Fluid f = fluid_by_name("EMI-BF4");
  f.dG_solvation = 0.05 * constants::eV;  // keep j well away from underflow

  for (int n : {80, 160, 320}) {
    Mesh m = make_sphere(1e-6, 0.0, n);
    for (Element& e : m.elems) e.tag = Tag::FreeSurface;
    m.finalize();
    BemSolver bem(m);
    bem.solve({1000.0, 0.0, 0.0});
    const IonEmission ie = integrate_ion_emission(bem, f, 298.15, false);
    char buf[96];
    std::snprintf(buf, sizeof buf, "A_eff == Kugelflaeche (n = %d)", n);
    check(buf, ie.effective_area, m.total_area(), 2e-3);
  }
}

// ==========================================================================
// B10  the toolchain workaround is a build setting, recorded here for
//      completeness: ES_NATIVE_ARCH defaults to OFF (see CMakeLists.txt).
// ==========================================================================

int main() {
  b1_solve_at_voltage();
  b2_static_fold();
  b3_output_identity();
  b4_space_charge_closed();
  b6_droplet_closed();
  b8_polarity_closed();
  b7_b11_wording();
  b9_effective_area();
  std::printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
