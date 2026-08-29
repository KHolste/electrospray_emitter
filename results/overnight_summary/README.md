# Nachtlauf 2026-08-29/30 — Übersicht

21 lokale Commits, **nichts gepusht**. Vollständiger Build aus sauberem
Zustand: fehlerfrei. `ctest`: **26/26 grün**. Sanitizerlauf unter WSL
(GCC 13.3, ASan + UBSan, `detect_leaks=1`): **26/26 grün**, 385 s.

`origin/main` steht auf `a932ff1`; HEAD ist **23 Commits voraus** (die zwei
P3b-Commits `0677e3b` und `6dc833b` plus die 21 dieses Laufs).

## Die Punkte

| Punkt | Status | Code-Commit | Artefakt-Commit | Tests | wichtigste Abbildung | Blocker |
|---|---|---|---|---|---|---|
| 0 — P3b numerisch bereinigen | `validated_subset` + `DiscretizationNotConverged` | `713ab86` | `85e7096` | `test_load_projection`, 40+ Prüfungen | `fig3_exclusion.png` — `kTolExclusion` widerlegt | 1-%-Ziel der gekoppelten Netzkonvergenz **verfehlt** (4,6–6,1 %); bei 1000 V nicht einmal asymptotisch |
| 1 — Druck- und Zulaufmodell | `implemented` | `1cb6108` | `3ffd040` | `test_feed` | `fig1_budget_vs_flow.png` | keiner; `p_reservoir` und `Q` bleiben Eingaben |
| 2 — Reale Stoffdaten | `implemented` (γ, ρ, µ, K) + `MissingMaterialData` (ε_r, ν) | `4c4bbfe` | `42e24d4` | `test_material_data` | `fig2_scatter.png` — Literaturstreuung je Quelle | ε_r und kinematische Viskosität ohne Quelle mit Methode+Reinheit+Wassergehalt |
| 3 — Endliche Leitfähigkeit und Strömung | `validated_subset` | `2a7205d` | `26e8a27` | `test_transport` | `fig1_pipe_flow.png` | τ = ε₀ε_r/σ **nicht berechenbar**; keine gekoppelte finite-conductivity-Meniskusrechnung |
| 4 — Zeitabhängige freie Oberfläche | `infrastructure_only` | `c4ecbb7` | `5067cea` | `test_surface_kinematics` | `fig1_advection.png` | Vorbedingung aus Punkt 3 **nicht erfüllt**: kein Feld mit freier Oberfläche, keine tangentiale Traktion q_s E_t, kein ε_r |
| 5 — Ionenemission | `blocked` | `baaf6fd` | `114768e` | `test_emission_contract` | `fig1_blocker.png` — Hebelwirkung der Barriere | Ratengleichung an keiner Primärquelle geprüft; kein belegtes ΔG für EMI-BF4 |
| 6 — Raumladung, Poisson/FEM und PIC-Grundlage | `validated_subset` | `8de4558` | `93b3be6` | `test_space_charge` | `fig2_checks.png` | keine Emissions-PIC-Schleife (P5 blockiert); Selbstfeld ist eine Netzgröße |
| 7 — Ionenbahnen und messbare Größen | `validated_subset` | `cc8f52f` | `3e97de0` | `test_particle_transport` | `fig1_trajectories.png` | keine physikalische Teilchenquelle → **Transportantwort**, keine Stromvorhersage |
| 8 — Cone-Jet und Tropfenbetrieb | `blocked` | `4274453` | `be884ab` | `test_cone_jet_contract` | `fig1_status.png` — Gültigkeitskarte | fünf von sieben Teilmodellen fehlen; Gañán-Calvo 1997 und Erratum 2000 nicht im Volltext erreichbar |
| 9 — 3D-Erweiterung und Validierung | `infrastructure_only` | `e8713ee` | `f6fc6f2` | `test_validation` | `fig1_validation.png` | kein 3D-Netz, kein 3D-Löser; keine Messdaten importiert |

## Die sieben Befunde, die etwas ändern

1. **γ war um 19 % zu klein.** Der bisher benutzte Wert 0,0452 N/m liegt unter
   *jeder* Quelle, die Reinheit und Wassergehalt angibt (niedrigste 0,0501). Die
   gewählte Quelle (Souckova et al. 2011, zwei Methoden, 99,8 %, 0,0164 % Wasser)
   gibt 0,05401 N/m. **Faktor 1,195 auf jede Druckskala, 1,093 auf jede
   Spannung** — größer als jeder Diskretisierungsfehler dieses Projekts.
2. **`kTolExclusion` ist widerlegt, nicht getauscht.** Die geprüfte Größe hat für
   p = C d^β einen geschlossenen Grenzwert ohne Netzweite; Messung trifft Formel
   auf 2·10⁻⁴. Die 11,3 % aus P3b sind der von β festgelegte Wert.
3. **Das 1-%-Ziel der gekoppelten Netzkonvergenz ist nirgends erreicht.** Die
   P3b-Ergebnisse sind qualitativ; keine Apexhöhe trägt drei Stellen.
   Die Randsingularität setzt die Rate: Ordnung ≈ 1+β ≈ 0,55.
4. **τ = ε₀ε_r/σ ist mit den belegten Daten nicht berechenbar.** Der
   Äquipotentialansatz von P3b ist damit eine **Annahme**, kein nachgewiesener
   Grenzfall.
5. **Die vorhandene Feldrekonstruktion ist achsennah nur erster Ordnung.**
   Ursache: die 2πr-Volumengewichtung, Unsymmetrie h/r. Abseits der Achse ist sie
   zweiter Ordnung. Betrifft auch P2b, P2c und P3b.
6. **Die Emissionsbarriere ist der Hebel.** Über die Literaturspanne 1,0–1,4 eV
   ändert sich die Rate bei festem Feld um den Faktor **5,8·10⁶**.
7. **Die 2πr-Wichtung ist als 3D-Integral nachgewiesen** (≤10⁻¹² gegen explizite
   3D-Quadratur) — die eine Sache, die zwischen 2D und 3D wirklich stimmt.

## Fehler, die die Tests im Code dieses Laufs gefunden haben

| Fehler | wo | wie gefunden |
|---|---|---|
| Vorzeichen im Energieinvariant (E = ½mv² **+** qφ) | P7 | Energietest meldete relativen Fehler von genau 2 |
| φ am Austritt außerhalb des Netzes ausgewertet | P7 | Energiebilanz um die volle Beschleunigungsspannung falsch |
| inkonsistenter letzter Zeitschritt | P7 | gemessene Zeitordnung kam **negativ** heraus |
| geteilter statischer Puffer in `synthetic_complete_model` | P5 | einer der beiden Polaritätsströme war null |
| zu starke Behauptung über γ | P2 | Test schlug fehl; Aussage auf das Haltbare korrigiert |
| `CMakeLists.txt` in `713ab86` nennt Dateien, die erst mit P1 kommen | P0 | beim Artefaktlauf bemerkt, nach vorn behoben (siehe `1cb6108`) |

## Was NICHT gerechnet wurde

Keine Emission. Keine selbstkonsistente Emissions-PIC-Schleife. Keine
Zweiphasenströmung, kein Jet, kein Tropfenzerfall. Kein dynamischer Meniskus.
Keine gekoppelte finite-conductivity-Meniskusrechnung. Keine Stabilitätsaussage,
kein Taylor-Kegel-Onset. Kein 3D-Netz, kein 3D-Löser, kein 3D-Ergebnis. Keine
importierten Messdaten. Keine Stromvorhersage.

## Reproduktion

```sh
cmake -S . -B build -G Ninja && cmake --build build && (cd build && ctest)

./build/es_p3b_audit     examples/device_p1.cfg examples/electrocapillary_p3b.cfg \
                         results/2026-08-29_p0_p3b_audit    meta.commit=$(git rev-parse HEAD)
./build/es_feed          examples/device_p1.cfg examples/feed_p1.cfg \
                         results/2026-08-29_p1_pressure_budget meta.commit=$(git rev-parse HEAD)
./build/es_material      results/2026-08-29_p2_material_data  meta.commit=$(git rev-parse HEAD)
./build/es_transport     results/2026-08-29_p3_transport      meta.commit=$(git rev-parse HEAD)
./build/es_kinematics    results/2026-08-29_p4_kinematics     meta.commit=$(git rev-parse HEAD)
./build/es_emission_audit results/2026-08-29_p5_emission_audit meta.commit=$(git rev-parse HEAD)
./build/es_space_charge  results/2026-08-29_p6_space_charge   meta.commit=$(git rev-parse HEAD)
./build/es_trajectories  results/2026-08-29_p7_trajectories   meta.commit=$(git rev-parse HEAD)
./build/es_cone_jet      results/2026-08-29_p8_cone_jet       meta.commit=$(git rev-parse HEAD)
./build/es_validation    results/2026-08-29_p9_validation     meta.commit=$(git rev-parse HEAD)
```

Danach je Ordner das zugehörige `python/plot_*.py`. Der Punkt-0-Lauf braucht
rund 20 Minuten (gekoppelte Netzstufe 4 bei 1400 V und Kanten-Gate bis Stufe 5,
3 GiB Bandfaktorisierung); alle anderen laufen in Sekunden bis wenigen Minuten.

Die Stoffdatentabelle wird von `python tools/fetch_material_data.py` neu erzeugt;
ein leerer `git diff` danach ist die Reproduzierbarkeitsprüfung.
