# electrospray — Simulation von Kolloid-/Electrospray-Emittern

> ## ⚠ Status: Prototyp, keine validierte Simulation
>
> Dieser Code ist ein **Prototyp in Neuausrichtung**. Er ist nicht geeignet, um
> Emitter auszulegen oder Messdaten zu interpretieren.
>
> Belastbar geprüft ist ausschließlich die **raumladungsfreie Elektrostatik**
> (gegen analytische Lösungen für Kugel und Rotationsellipsoid) sowie die
> **Bahnintegration** (Energieerhaltung). Die eigentliche Emitterphysik —
> Meniskus im Betrieb, Emission, Raumladung — ist **nicht validiert** und in
> Teilen nachweislich fehlerhaft.
>
> Die Neuausrichtung, die bestätigten Fehlerbefunde und der Umbauplan stehen in
> **[docs/](docs/)**. Vor jeder Benutzung bitte
> [docs/01_gap_analysis.md](docs/01_gap_analysis.md) lesen.

---

Axialsymmetrischer Simulationscode für Emitter ionischer Flüssigkeiten.
C++20, keine externen Abhängigkeiten außer CMake und einem C++-Compiler.

## 1. Was der Code enthält

| Baustein | Was er tut | Prüfstand |
|---|---|---|
| Axialsymmetrische BEM | Laplace-Gleichung für vorgegebene Elektrodenpotentiale | gegen analytische Lösungen geprüft |
| Statischer Meniskus | Young-Laplace mit Maxwell-Zug, perfekt leitend, ohne Strömung | im feldfreien Grenzfall geprüft; im Betriebsfall nicht |
| Emissionsformeln | Iribarne-Thomson-Rate; Cone-Jet-Korrelationen | Formeln korrekt programmiert; keine Validierung des Ergebnisses |
| Strahlverfolgung | Bahnintegration im gelösten Feld | Energieerhaltung geprüft; Raumladungsmodell unbrauchbar |

Die Sprachregelung folgt [docs/README.md](docs/README.md): *geprüft* heißt
Vergleich mit einer unabhängigen Referenz mit angegebener Abweichung;
*validiert* wird für keinen Teil der Emitterphysik beansprucht.

## 2. Was geprüft ist

| Prüfung | Referenz | Abweichung |
|---|---|---|
| Elliptische Integrale K, E | exakte Werte, unabhängige Quadratur | 10⁻¹⁴ |
| Kapazität isolierte Kugel | 4πε₀R | 2,5·10⁻⁶ |
| Spitzenfeld Rotationsellipsoid (Aspekt 20) | prolat-sphäroidale Analytik | 1,8·10⁻³ |
| Randbedingungsresiduum | Dirichletdaten an Kollokationspunkten | 1,5·10⁻¹⁵ |
| Fernfeld-Monopol | Q/4πε₀R | 10⁻³ |
| Meniskus ohne Feld | Kugelkappe R = 2γ/Δp | 4·10⁻⁶ |
| Energieerhaltung Trajektorie | E_kin/q = Potentialdifferenz | 7·10⁻⁶ |
| Ring-Kern (Fernfeld, Gradient) | Punktladung, E = −∇V | 10⁻⁵ |

Das ist die Elektrostatik und die Bahnintegration. Vollständige Matrix mit
offenen Punkten: [docs/06_validation_matrix.md](docs/06_validation_matrix.md).

Ein Hinweis zum Kegelhalbwinkel: der Meniskusast läuft gegen 49,16° und damit
in die Nähe von Taylors 49,29°. Das ist ein **Diskretisierungstest** des
Young-Laplace-Lösers, liegt auf dem *instabilen* Ast und ist kein Nachweis über
Emissionsphysik, Betriebspunkt oder Stabilität. In einer früheren Fassung dieses
README stand er in dieser Rolle — das war überzogen.

## 3. Was nicht geprüft ist — und was in Phase P0 korrigiert wurde

Ausführlich in [docs/01_gap_analysis.md](docs/01_gap_analysis.md). Die elf
bestätigten Befunde sind testgetrieben bearbeitet
(`tests/test_regressions.cpp`); die inhaltlichen Grenzen bleiben bestehen.

**Weiterhin gültig — daran hat P0 nichts geändert:**

* **Der emittierende Betrieb ist nicht modelliert.** Der Meniskus wird statisch
  und perfekt leitend gerechnet, die Ionenrate danach auf dessen Feld
  angewandt. Higuera (2008) zeigt, dass im rein ionischen Regime der Strom von
  der endlichen Leitfähigkeit kontrolliert wird; Collins et al. (2008) zeigen,
  dass Tip Streaming im perfekt leitenden Grenzfall gar nicht auftritt. Die
  ausgegebenen Ionenströme sind **keine Vorhersage** — `es_operate` sagt das
  jetzt in jedem Lauf ausdrücklich.
* **Der berechnete Umkehrpunkt ist kein Emissions-Onset.** Er heißt im Code
  jetzt „static fold". Eine Stabilitätsanalyse findet nicht statt; ein Onset
  darf erst nach physikalischem Nachweis ausgewiesen werden (Phase P3).
* **Stoffdaten ohne Einzelnachweis.** Die Werte in `src/fluid.cpp` sind
  literaturtypisch, aber nicht quellenbelegt und nicht zitierfähig (Phase P1).
* **Die statische Rechnung löst die nm-Emissionsstruktur nicht auf**, aus der
  im PIR tatsächlich emittiert wird.

**Bewusst deaktiviert — schlägt geschlossen fehl statt scheinbar zu rechnen:**

| Funktion | Verhalten | vorgesehen |
|---|---|---|
| Raumladung (`beam.space_charge_iters > 0`) | `NotImplementedInThisPhase`, Exit 3 | P4, Poisson-FEM/FVM mit PIC |
| Tropfenstrahl (`beam.species = droplet\|both`) | `NotImplementedInThisPhase`, Exit 3 | P6, Kopplung an das Cone-Jet-Modell |
| Negative Polarität (U < 0) | `NotImplementedInThisPhase`, Exit 3 | P4/P5, Anionenspezies |

Das Ring-Makroteilchenmodell wurde **nicht** notdürftig regularisiert. Sein
Eigenfeld divergiert weiterhin; genau das ist der Grund für die Sperre.

**In P0 korrigiert:**

* `solve_at_voltage` meldet Erfolg nur noch, wenn die gelieferte Spannung die
  angeforderte trifft; sonst ein eindeutiger Status mit Ursache.
* `find_static_fold` verlangt ein nachgewiesenes inneres Maximum. Ein einzelner
  oder monotoner Ast liefert keinen Umkehrpunkt mehr.
* Ausgabedateien tragen Anwendung, Zustand und Spannung im Namen plus einen
  Provenienzkopf; `es_meniscus` und `es_beam` überschreiben einander nicht mehr.
  Oberflächen- und Meshdaten gehören nachweislich zum ausgewiesenen Zustand.
* Die emittierende Fläche ist das stetige Funktional
  `A_eff = (∫j dA)² / ∫j² dA` statt einer Summe ganzer Elemente.
* Netzkonvergenz wird für Faltenspannung, Apexfeld, Apexradius, Ionenstrom und
  `A_eff` geprüft, nicht mehr nur für die Spannung.
* Neu gefunden und behoben: `BemSolver::set_mesh` invalidierte die
  Basislösung nicht, wodurch ein späteres `solve()` die Basis des alten Netzes
  auf die neue Geometrie superponierte.

## 4. Netzerzeugung

Die Zielfassung erzeugt das Netz **automatisch**: der Benutzer gibt Geometrie-
und Genauigkeitsparameter vor, kein Netz. Größenverteilung aus Krümmung und
Merkmalsgröße, danach a-posteriori-Verfeinerung an Meniskusspitze,
Kontaktlinie, Elektrodenkanten, Extraktoröffnung und Strahlkern, mit Abbruch
über Netzunabhängigkeit von Zielgrößen (Apexfeld, Strom, emittierende Fläche,
Divergenz). Spezifikation:
[docs/04_geometry_model.md](docs/04_geometry_model.md), Abschnitt 4.4.

**Der aktuelle Prototyp leistet das nicht.** Dort werden Elementgrößen
(`h_tip`, `h_far`) von Hand im Konfigurationsfile gesetzt; es gibt weder einen
Fehlerschätzer noch Verfeinerung.

## 5. Bauen

```sh
cmake -S . -B build -G Ninja
cmake --build build
cd build && ctest --output-on-failure
```

Getestet mit MinGW-w64 GCC 16.1 unter Windows 11, CMake 4.3, Ninja, sowie GCC
13.3 unter Ubuntu 24.04 (WSL2). OpenMP wird genutzt, wenn vorhanden.

### `-march=native` ist abgeschaltet — Grund

Mit AVX-Codegen (`-mavx`, `-mavx2`, und damit `-march=native` auf gängigen
CPUs) stürzt der Code unter MinGW-w64 GCC 16.1 **nichtdeterministisch in etwa
der Hälfte der Läufe** mit SIGSEGV ab. Ursache ist im Disassemblat belegt:

```
=> vmovdqa %ymm0,0x20(%rsp)      rsp = 0x5ffaf0  ->  Ziel 0x5ffb10
```

`vmovdqa` verlangt 32-Byte-Ausrichtung, die Zieladresse ist aber nur
16-Byte-ausgerichtet. Der Prolog von `main()` richtet den Stack nicht nach
(kein `and $-32,%rsp`), emittiert aber fünf alignment-pflichtige
256-Bit-Zugriffe. Die Win64-ABI garantiert nur 16 Byte; ob der Stack zufällig
32-Byte-ausgerichtet startet, entscheidet die Adressraum-Randomisierung — daher
die Nichtdeterminismus.

Nicht reproduzierbar unter Linux/GCC 13.3. ASan und UBSan melden nichts.
`-mstackrealign` hilft nicht. `-mprefer-vector-width=128` vermeidet es.

Wer die Host-CPU-Codegen trotzdem will:
`cmake -S . -B build -DES_NATIVE_ARCH=ON` — mit dem obigen Vorbehalt.

## 6. Anwendungen

Konfigurationsdatei plus `key=value`-Overrides auf der Kommandozeile
(Kommandozeile schlägt Datei). Werte tragen Einheiten: `10um`, `1.5kV`,
`500Pa`, `0.5nL/s`, `25C`, `15deg`. Nicht gelesene Schlüssel werden gemeldet.

| Programm | Zweck | Vertrauenswürdigkeit |
|---|---|---|
| `es_field` | Laplace-Lösung, Spitzenfeld, Feldüberhöhung, Feldgitter | Elektrostatik geprüft |
| `es_meniscus` | Fortsetzung des statischen Astes | Ast ja, „Onset"-Etikett nein |
| `es_operate` | Betriebspunkt, Ionenstrom, Schub/Isp | **nicht belastbar**, siehe Abschnitt 3 |
| `es_beam` | Strahltransport, Divergenz, Interzeption | ohne Raumladung brauchbar, mit nicht |
| `es_vacuum` | P2a: statische Vakuum-Elektrostatik auf der P1-Geometrie | **überholt** — behandelt den Emitterkörper als Metall. Bleibt als unabhängige BEM-Vergleichsrechnung für ε_r = 1 erhalten, siehe `docs/08_dielectric_model.md` |
| `es_dielectric` | P2b: dielektrische Elektrostatik des kapillaren Kunze-Emitters | kantenfernes Feld **netz**konvergent und gegen die BEM geprüft; Kantenfelder nicht; **trunkierungs**konvergent ist nichts — jede Größe gilt für das angegebene `feed.liquid_feed_z` |

```sh
./build/es_field --help          # Schlüsselreferenz
./build/es_meniscus examples/capillary_emibf4.cfg
python python/plot.py out_capillary

./build/es_vacuum examples/device_p1.cfg examples/vacuum_p2a.cfg \
    results/2026-08-29_p2a_vacuum_electrostatics meta.commit=$(git rev-parse HEAD)
python python/plot_vacuum.py results/2026-08-29_p2a_vacuum_electrostatics

./build/es_dielectric examples/device_p1.cfg examples/dielectric_p2b.cfg \
    results/2026-08-29_p2b_dielectric_electrostatics meta.commit=$(git rev-parse HEAD)
python python/plot_dielectric.py results/2026-08-29_p2b_dielectric_electrostatics
```

Bekanntes Problem: verschiedene Anwendungen überschreiben bei gleichem
`output.prefix` einander kommentarlos, und die geschriebenen Dateien gehören
nicht zwingend zu dem Zustand, den der Bericht ausweist.

## 7. Struktur

```
include/es/  types constants elliptic linalg geometry bem fluid
             emission meniscus beam config io status
             device_geometry boundary_mesh vacuum_bem
src/         Implementierungen
apps/        es_field es_meniscus es_operate es_beam es_vacuum
tests/       10 Testprogramme
examples/    Konfigurationsbeispiele
python/      Auswertung und Plots (keine Rechnung)
docs/        Neuausrichtung, Modellspezifikation, Umbauplan
```

## 8. Weiteres Vorgehen

Kein neuer Funktionsumfang, bis die Punkte aus
[docs/05_implementation_plan.md](docs/05_implementation_plan.md) Phase P0
abgearbeitet sind. Danach Geometrie und automatische Vernetzung (P1),
Elektrostatik konsolidieren (P2), Meniskus neu (P3), Volumenlöser und PIC (P4),
emittierender Betrieb im PIR (P5).
