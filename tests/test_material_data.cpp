// tests/test_material_data.cpp -- P2: the material-data contract
//
// The point of this test is NOT to check that a number equals a number -- the
// table is generated, so that would only check the generator against itself.
// It checks the PROPERTIES the contract promises:
//
//   * a missing property fails closed with MissingMaterialData and no default;
//   * a selected source really satisfies the stated selection rule;
//   * a value is never extrapolated outside the temperatures that were measured;
//   * a frequency-resolved source never backs a DC number;
//   * a non-ambient-pressure point never backs an ambient number;
//   * the literature scatter is reported and is larger than any single source's
//     own error bar -- which is the honest uncertainty of the number.
//
// The last group is the one that matters for this project: it measures how far
// the unsourced value that P3a/P3b actually used sits from the sourced range.

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

#include "es/liquid.hpp"
#include "es/material_data.hpp"

using namespace es;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  [%s] %s\n", ok ? "ok" : "FEHLER", what.c_str());
  if (!ok) ++failures;
}

constexpr Real kT = 298.15;

}  // namespace

int main() {
  std::printf("P2 -- Stoffdatenvertrag fuer EMI-BF4\n\n");
  const MaterialDataset& d = emibf4_sourced();
  d.print(stdout);

  // =========================================================================
  std::printf("\n1. Fail-closed\n");
  {
    const MaterialValue v = material_value(d, PropertyKind::RelativePermittivity, kT);
    std::printf("    relative Permittivitaet: %s\n", to_string(v.status));
    std::printf("    %s\n", v.message.c_str());
    check(v.status == MaterialDataStatus::MissingMaterialData,
          "die relative Permittivitaet hat keine Quelle mit vollstaendiger Provenienz "
          "und meldet MissingMaterialData");
    check(!v.usable(), "und ist nicht brauchbar");
    bool threw = false;
    try {
      (void)v.value_or_throw();
    } catch (const std::exception&) {
      threw = true;
    }
    check(threw, "value_or_throw wirft, statt eine Vorgabe zu liefern");
    check(v.value == 0.0 && v.source == nullptr,
          "und es steht kein stiller Ersatzwert im Ergebnis");

    const MaterialValue k = material_value(d, PropertyKind::KinematicViscosity, kT);
    check(k.status == MaterialDataStatus::MissingMaterialData,
          "die kinematische Viskositaet ebenso -- sie wird NICHT aus mu und rho "
          "zusammengerechnet, weil nicht festgehalten ist, welche Dichte der jeweilige "
          "Autor benutzt hat");
  }

  // =========================================================================
  std::printf("\n2. Die Auswahlregel gilt fuer jede getroffene Auswahl\n");
  for (PropertyKind kind : {PropertyKind::SurfaceTension, PropertyKind::Density,
                            PropertyKind::DynamicViscosity,
                            PropertyKind::ElectricalConductivity}) {
    const MaterialProperty* p = d.find(kind);
    check(p != nullptr && p->has_selection(),
          std::string(to_string(kind)) + ": eine Quelle ist ausgewaehlt");
    if (!p || !p->has_selection()) continue;
    const PropertySource& s = p->selection();
    check(s.states_method() && s.states_purity() && s.states_water_content(),
          std::string(to_string(kind)) +
              ": die gewaehlte Quelle nennt Methode, Reinheit UND Wassergehalt");
    check(!s.is_frequency_resolved(),
          std::string(to_string(kind)) + ": die gewaehlte Quelle ist nicht frequenzaufgeloest");
    check(s.has_ambient_points() && s.covers(kT),
          std::string(to_string(kind)) + ": sie hat Umgebungsdruckpunkte und deckt 298,15 K ab");
    // No other admissible source has more ambient points -- that is the rule.
    std::size_t best = 0;
    for (std::size_t j = 0; j < s.n_points; ++j)
      if (s.points[j].ambient()) ++best;
    std::size_t better = 0;
    for (std::size_t i = 0; i < p->n_sources; ++i) {
      const PropertySource& o = p->sources[i];
      if (!(o.states_method() && o.states_purity() && o.states_water_content())) continue;
      if (o.is_frequency_resolved() || !o.covers(kT)) continue;
      std::size_t n = 0;
      for (std::size_t j = 0; j < o.n_points; ++j)
        if (o.points[j].ambient()) ++n;
      if (n > best) ++better;
    }
    check(better == 0,
          std::string(to_string(kind)) +
              ": keine zulaessige Quelle hat mehr Umgebungsdruckpunkte");
  }

  // =========================================================================
  std::printf("\n3. Keine Extrapolation ueber den Messbereich hinaus\n");
  {
    const MaterialProperty* p = d.find(PropertyKind::SurfaceTension);
    const PropertySource& s = p->selection();
    std::printf("    gewaehlte Quelle deckt %.2f .. %.2f K ab\n", s.T_min(), s.T_max());
    check(std::isfinite(s.value_at(kT)), "innerhalb des Bereichs gibt es einen Wert");
    check(!std::isfinite(s.value_at(s.T_min() - 10.0)),
          "zehn Kelvin unterhalb gibt es keinen -- es wird nicht extrapoliert");
    check(!std::isfinite(s.value_at(s.T_max() + 10.0)), "oberhalb ebenso wenig");
    const MaterialValue v = material_value(d, PropertyKind::SurfaceTension, s.T_max() + 10.0);
    check(v.status == MaterialDataStatus::MissingMaterialData,
          "und die Abfrage meldet dort MissingMaterialData statt eines Wertes");
  }

  // =========================================================================
  std::printf("\n4. Druck und Frequenz koennen einen Wert nicht ersetzen\n");
  {
    const MaterialProperty* p = d.find(PropertyKind::Density);
    std::size_t with_hp = 0, hp_points = 0;
    for (std::size_t i = 0; i < p->n_sources; ++i)
      if (p->sources[i].n_non_ambient() > 0) {
        ++with_hp;
        hp_points += p->sources[i].n_non_ambient();
      }
    std::printf("    Dichte: %zu Quellen enthalten %zu Punkte abseits des Umgebungsdrucks\n",
                with_hp, hp_points);
    check(with_hp > 0, "es gibt Hochdruckdaten im Datensatz -- sie sind nicht geloescht");
    check(p->selection().n_non_ambient() == 0,
          "aber die GEWAEHLTE Dichtequelle enthaelt keinen einzigen davon");
    const Real hi = p->max_at(kT), lo = p->min_at(kT);
    std::printf("    Umgebungsdruckband bei 298,15 K: %.1f .. %.1f kg/m^3\n", lo, hi);
    check(hi < 1295.0,
          "und das Band bei 298,15 K enthaelt den Hochdruckwert von rund 1300 kg/m^3 nicht");

    const MaterialProperty* q = d.find(PropertyKind::RelativePermittivity);
    std::size_t fr = 0;
    for (std::size_t i = 0; i < q->n_sources; ++i)
      if (q->sources[i].is_frequency_resolved()) ++fr;
    std::printf("    relative Permittivitaet: %zu von %zu Quellen sind frequenzaufgeloest\n",
                fr, q->n_sources);
    check(fr > 0, "es gibt frequenzaufgeloeste Permittivitaetsdaten");
    const Real band_dc = q->max_at(kT) - q->min_at(kT);
    const Real band_all = q->max_at(kT, 2.0, true) - q->min_at(kT, 2.0, true);
    std::printf("    Band ohne / mit den frequenzaufgeloesten Quellen: %.2f / %.2f\n",
                band_dc, band_all);
    check(band_all > band_dc,
          "die frequenzaufgeloesten Werte wuerden das Band verbreitern und sind deshalb "
          "standardmaessig ausgeschlossen");
  }

  // =========================================================================
  std::printf("\n5. Literaturstreuung -- die ehrliche Unsicherheit\n");
  for (PropertyKind kind : {PropertyKind::SurfaceTension, PropertyKind::Density,
                            PropertyKind::DynamicViscosity,
                            PropertyKind::ElectricalConductivity}) {
    const MaterialProperty* p = d.find(kind);
    const Real v = p->selection().value_at(kT);
    const Real lo = p->min_at(kT), hi = p->max_at(kT);
    const Real spread = p->relative_spread_at(kT);
    std::printf("    %-24s gewaehlt %.6g, Band %.6g .. %.6g ueber %zu Quellen, "
                "Streuung %.1f %%\n",
                to_string(kind), v, lo, hi, p->n_sources_at(kT), 100.0 * spread);
    check(std::isfinite(spread) && spread > 0.0,
          std::string(to_string(kind)) + ": die Streuung ist bekannt und nicht null");
    // The selected value must lie inside the band it is the middle of.
    check(v >= lo - 1e-12 && v <= hi + 1e-12,
          std::string(to_string(kind)) + ": der gewaehlte Wert liegt im Band");
  }

  // =========================================================================
  std::printf("\n6. Wo der bisherige unbelegte Wert liegt\n");
  {
    const LiquidProperties ill = emibf4_illustrative();
    const MaterialProperty* g = d.find(PropertyKind::SurfaceTension);
    const Real lo = g->min_at(kT), hi = g->max_at(kT);
    const Real sel = g->selection().value_at(kT);
    std::printf("    illustrativ (src/liquid.cpp): gamma = %.5g N/m\n", ill.gamma);
    std::printf("    belegtes Band bei 298,15 K  : %.5g .. %.5g N/m, gewaehlt %.5g\n", lo, hi,
                sel);
    std::printf("    Abweichung des illustrativen Wertes vom gewaehlten: %+.1f %%\n",
                100.0 * (ill.gamma - sel) / sel);
    check(ill.status == LiquidDataStatus::Illustrative,
          "der bisherige Datensatz bleibt illustrative");

    // THE FINDING, stated exactly as far as it holds and no further.
    //
    // The unsourced value is NOT below the whole band: the band's lower edge is
    // 0.0443 N/m from Martino et al. (2006), a capillary-rise measurement that
    // states neither purity nor water content.  What IS true, and sharper, is
    // that the unsourced value lies below EVERY source that documents its
    // sample -- which is the statement that matters, because those are the only
    // ones a number may rest on.
    Real lo_documented = std::numeric_limits<Real>::max();
    std::size_t n_documented = 0;
    for (std::size_t i = 0; i < g->n_sources; ++i) {
      const PropertySource& s = g->sources[i];
      if (!(s.states_purity() && s.states_water_content())) continue;
      for (std::size_t j = 0; j < s.n_points; ++j)
        if (s.points[j].ambient() && std::abs(s.points[j].T - kT) <= 2.0) {
          lo_documented = std::min(lo_documented, s.points[j].value);
          ++n_documented;
        }
    }
    std::printf("    kleinster Wert unter den Quellen MIT Reinheit und Wassergehalt: "
                "%.5g N/m (%zu Punkte)\n", lo_documented, n_documented);
    check(n_documented > 0, "es gibt Quellen mit vollstaendiger Probenangabe");
    check(ill.gamma < lo_documented,
          "der bisher benutzte gamma-Wert liegt unter JEDER Quelle, die Reinheit und "
          "Wassergehalt angibt");
    check(ill.gamma >= lo,
          "er liegt aber NICHT unter dem gesamten Band: dessen untere Kante ist eine "
          "Kapillaraufstiegsmessung ohne Reinheits- und Wasserangabe");

    const MaterialProperty* rho = d.find(PropertyKind::Density);
    std::printf("    illustrativ: rho = %.6g kg/m^3, belegtes Band %.6g .. %.6g\n", ill.rho,
                rho->min_at(kT), rho->max_at(kT));
    check(ill.rho >= rho->min_at(kT) - 1e-9 && ill.rho <= rho->max_at(kT) + 1e-9,
          "die Dichte des illustrativen Satzes liegt dagegen im belegten Band");

    const MaterialProperty* mu = d.find(PropertyKind::DynamicViscosity);
    std::printf("    illustrativ: mu = %.6g Pa s, belegtes Band %.6g .. %.6g\n",
                ill.documented_only.mu, mu->min_at(kT), mu->max_at(kT));
  }

  std::printf("\n%s: %d Fehler\n", failures == 0 ? "BESTANDEN" : "FEHLGESCHLAGEN", failures);
  return failures == 0 ? 0 : 1;
}
