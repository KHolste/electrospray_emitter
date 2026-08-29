# P2a — statische Vakuum-Elektrostatik auf der P1-Geometrie — 2026-08-29

Alle Daten und Abbildungen in diesem Ordner sind in diesem Lauf frisch erzeugt.

## Was gerechnet wurde — und was nicht

Genau ein Problem:

* ladungsfreies Vakuum, `rho = 0`;
* Emittermetall **und die anfängliche ebene Flüssigkeitsoberfläche bei z = 0**
  auf `V_emitter`, die ebene Fläche ausschließlich als
  Perfect-Conductor-Referenz;
* Extraktor auf `V_extractor`;
* `V -> 0` im Unendlichen.

**Nicht enthalten:** Meniskusverformung, endliche Flüssigkeitsleitfähigkeit,
Ionen- oder Tröpfchenemission, Strömung, Raumladung. Die Fläche bei z = 0 ist
die anfängliche ebene Flüssigkeitsoberfläche — kein berechneter Meniskus, kein
emittierender Meniskus und keine Aussage über das Pure-Ion-Regime.

Die P0-BEM bleibt unverändert und behält ihre eigenen Geometrien und Netze.
P2a benutzt denselben BEM-Kern über einen ausdrücklichen Adapter, nicht statt
ihm.

## Reproduktion

```sh
cmake -S . -B build -G Ninja
cmake --build build
./build/es_vacuum.exe examples/device_p1.cfg examples/vacuum_p2a.cfg \
    results/2026-08-29_p2a_vacuum_electrostatics \
    meta.commit=$(git rev-parse HEAD) \
    meta.config="examples/device_p1.cfg + examples/vacuum_p2a.cfg"
python python/plot_vacuum.py results/2026-08-29_p2a_vacuum_electrostatics
```

Die verwendeten Konfigurationsdateien liegen als Kopie daneben
(`device_p1.cfg`, `vacuum_p2a.cfg`); die aufgelösten Werte in
`parameters.csv`, die Herkunft in `meta.txt`, der volle Bericht in
`report.txt`. `es_vacuum` gibt 2 zurück, wenn eine Netzprüfung oder eine
numerische Prüfung fehlschlägt — dieser Lauf gab 0 zurück.

## Der BEM-Kern und was er mathematisch kann

Aus `src/bem.cpp` abgelesen, nicht angenommen:

| Frage | Antwort |
|---|---|
| Formulierung | indirekte Einfachschicht; Unbekannte ist `sigma`, Kollokation in den Elementmitten |
| Randbedingung im Unendlichen | `V -> 0`, in der azimutal integrierten Freiraum-Greensfunktion enthalten; keine Abschneidefläche |
| geschlossene Leiter nötig? | nein — die Einfachschicht ist auf offenen wie geschlossenen Flächen wohlgestellt. Auf einer offenen Fläche ist `sigma` die **Summe** beider Seiten, `E_n = sigma/eps0` dort also kein einseitiges Feld |
| mehrere Leiter | Gruppierung nach `es::Tag` in bis zu drei Elektroden, eine Einheitspotential-Basis je Elektrode, Superposition |
| Falle | `Tag::Other` bildet auf `Electrode::Collector` ab — ein unmarkiertes Panel würde stillschweigend zu einer dritten, geerdeten Elektrode. Der Adapter erzeugt `Tag::Other` nie |

Folge aus der Randbedingung im Unendlichen: die angelegten Potentiale sind
absolut gegen unendlich, und die Gesamtladung ist **nicht** null.

## Was in die Vakuum-BEM eingeht

Ausgewählt wird ausschließlich über Rand- und Gebietskennungen, nie über
Koordinaten oder Reihenfolge: ein Element geht ein, wenn genau eine seiner
beiden Seiten Vakuum ist und die andere ein Leiter (Emitterfestkörper,
Flüssigkeitssäule, Extraktorfestkörper). Die vollständige Entscheidungstabelle
steht in `bem_selection.csv`.

| Randkurve | Kennung | übernommen |
|---|---|---|
| `free_surface_reference` | `free_surface_reference` | ja → ebene Flüssigkeitsoberfläche |
| `emitter_tip_land` | `emitter_tip_land` | ja → Emitter |
| `emitter_outer_surface` | `emitter_outer_surface` | ja → Emitter |
| `extractor_surface.aperture` / `.front` / `.back` / `.rim` | `extractor_surface` | ja → Extraktor |
| `symmetry_axis.liquid` / `.vacuum` | `symmetry_axis` | **nein** — kein Interface, keine Rotationsfläche |
| `bore_wall` | `bore_wall` | **nein** — innere Materialgrenze, kein Vakuumkontakt |
| `liquid_inlet` | `liquid_inlet` | **nein** — offener Domänenrand |
| `open_boundary.*` | `open_boundary` | **nein** — kein Leiter, keine Randbedingung |

Die offene rechteckige Rechendomäne ist kein Leiter und wird an keiner Stelle
als BEM-Rand verwendet.

## Pflichtangabe `extractor_outer_radius`

Bisher bedeutete `0` „die Elektrode reicht bis zum Domänenrand". Das setzte
einen Leiter mit der offenen Kante der Rechenbox gleich, die keiner ist, und
versteckte damit eine fehlende physikalische Abmessung. Der Außenradius ist
jetzt Pflicht, `0` wird abgelehnt, und `domain_radius` muss **echt** größer
sein. Der Wert im Beispielparametersatz (2 mm) ist ausdrücklich ein
Beispielwert. Wie stark das Ergebnis davon abhängt, steht in
`extractor_radius_study.csv` und in Abbildung 3(c): zwischen 0.5 mm und 2.5 mm
ändert sich `C_m` um den Faktor 2.

## Der offene Emitterbogen

Der Emitterbogen ist am Schnitt durch die Domänensohle offen; der offene
Domänenrand wurde **nicht** als Deckel benutzt. Der dahinterliegende Hohlraum
ist ein Rohr vom Radius 20 µm und der Länge 400 µm und deshalb stark
abgeschirmt — gemessen: `max |V − V_emitter| = 3.9e-2 V` von 1500 V, also
2.6e-5 relativ. Erst das macht `sigma/eps0` auf Stirnfläche und ebener
Flüssigkeitsoberfläche zum einseitigen Vakuumfeld. Am offenen Bogenende selbst
hat die Einfachschichtdichte eine `1/sqrt(d)`-Singularität; diese Stelle ist
wie eine scharfe Kante markiert und aus jeder Feldauswertung ausgeschlossen.

## Kapazitätsgrößen — benannt, nicht „C"

Mit `V = 0` im Unendlichen hat ein Zweileitersystem eine 2×2-Maxwell-Matrix:

```
Q_E = c_EE V_E + c_EX V_X
Q_X = c_XE V_E + c_XX V_X
```

| Größe | Wert (Referenzstufe 2, 341 Panels) |
|---|---|
| `c_EE` Maxwell-Selbstkoeffizient Emitter | 8.665508961e-15 F |
| `c_EX = c_XE` Influenzkoeffizient | −6.846168706e-15 / −6.846168792e-15 F |
| Reziprozitätsfehler | 1.2e-08 |
| `c_XX` Maxwell-Selbstkoeffizient Extraktor | 1.537010769e-13 F |
| `C_m = −c_EX` Gegenkapazität | 6.846168749e-15 F |
| `C_E,inf = c_EE + c_EX` gegen unendlich | 1.819340255e-15 F |
| `Q_E/(V_E − V_X)` am Arbeitspunkt | 8.665508961e-15 F |

`Q_E/(V_E − V_X)` ist bei geerdetem Extraktor identisch `c_EE` und **nicht**
`C_m`, weil das System bei `V(unendlich) = 0` Nettoladung trägt
(`Q_E + Q_X = 2.73e-12 C ≠ 0`).

## Netzkonvergenz

Vier automatisch erzeugte Stufen; `size_scale` ist ein gleichmäßiger Faktor auf
das automatische Größenfeld und ausschließlich ein Konvergenzwerkzeug.

| Stufe | size_scale | Randelemente | BEM-Panels | `c_EE` (F) | `C_m` (F) | `E_z`(ref) (V/m) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1.0 | 404 | 161 | 8.661339e-15 | 6.842876e-15 | 4.933072e+07 |
| 1 | 0.5 | 593 | 228 | 8.664084e-15 | 6.845055e-15 | 4.932283e+07 |
| 2 | 0.25 | 942 | 341 | 8.665509e-15 | 6.846169e-15 | 4.931964e+07 |
| 3 | 0.125 | 1611 | 544 | 8.666204e-15 | 6.846704e-15 | 4.931834e+07 |

Relative Änderung zwischen den beiden feinsten Stufen: `c_EE` 8.0e-05,
`C_m` 7.8e-05, `E_z`(ref) 2.6e-05.

Referenzpunkt für das konvergierte Axialfeld: `r = 0`, `z = 0.5 µm`, also ein
Zehntel Bohrungsradius über der Mitte der ebenen Flüssigkeitsoberfläche. Der
Punkt ist durch die Geometrie festgelegt, auf jeder Netzstufe derselbe und
1.005 Bohrungsradien von der nächsten unverrundeten Kante entfernt — das
Vierfache des Markierungsradius jener Kantenzone.

## Was hier ausdrücklich **keine** konvergierte Größe ist

Austrittskante, äußere Stirnkante und beide Aperturkanten sind unverrundet, das
offene Bogenende ist ein Modellschnitt. Das Feld folgt dort der Elementgröße:
zwischen `size_scale` 0.5 und 0.125 wächst das Maximum am Kantenknoten um 85 %,
während der berichtete kantenferne Wert um 6 % wandert. Die markierten Zonen
(`edge_zones.csv`, Radius = lokale Merkmalsgröße / 4, rein geometrisch) sind in
Abbildung 2 weiß ausgespart. Verrundungsradien sind ein P3-Parameter.

## Inhalt

| Datei | Inhalt |
|---|---|
| `report.txt` | vollständiger Textbericht des Laufs |
| `meta.txt` | Herkunft: Commit, Konfiguration, Spannungen, Referenzstufe |
| `bem_selection.csv` | Entscheidung je Randkurve mit Begründung |
| `convergence.csv` | alle vier Netzstufen, alle ausgewerteten Größen |
| `extractor_radius_study.csv` | Einfluss des Elektrodenaußenradius |
| `linearity.csv` | Skalierung von `sigma`, `V`, `E`, `Q` mit der Spannung |
| `polarity.csv` | vollständige Polaritätsumkehr |
| `surface.csv` | `sigma`, `E_n`, Kennung und Kantenzonenflag je Panel |
| `edge_zones.csv` | markierte Stellen, an denen `\|E\|` nicht konvergiert |
| `axis_profile.csv` | `V` und `E_z` auf `r = 0` |
| `field_full.csv` | Feldgitter über die ganze Domäne |
| `field_tip.csv` | Feldgitter an der Emitterspitze |
| `field_aperture.csv` | Feldgitter an der Extraktoröffnung |
| `field_truncation.csv` | Feldgitter am offenen Bogenende |
| `parameters.csv`, `regions.csv`, `boundaries.csv`, `features.csv` | Geometrie als Hintergrund der Abbildungen |
| `device_p1.cfg`, `vacuum_p2a.cfg` | Kopien der verwendeten Eingaben |
| `fig1_potential.png` | Potentialkarte mit Equipotentiallinien |
| `fig2_field_magnitude.png` | `\|E\|`-Karten mit ausgesparten Kanten |
| `fig3_convergence.png` | Netzkonvergenz und Elektrodenradiusstudie |
| `fig4_surface_charge.png` | Flächenladungsdichte entlang der Leiter |
