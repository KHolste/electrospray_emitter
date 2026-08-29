# P2a — statische Vakuum-Elektrostatik auf der P1-Geometrie — 2026-08-29

Alle Daten und Abbildungen in diesem Ordner sind in diesem Lauf frisch erzeugt,
aus einem sauberen Arbeitsbaum, und jede Abbildung nennt den tatsächlichen
HEAD-Commit. Die Herkunft steht in `figures_provenance.txt`.

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
ctest --test-dir build
git status --porcelain          # muss leer sein
./build/es_vacuum.exe examples/device_p1.cfg examples/vacuum_p2a.cfg \
    results/2026-08-29_p2a_vacuum_electrostatics \
    meta.commit=$(git rev-parse HEAD) \
    meta.config="examples/device_p1.cfg + examples/vacuum_p2a.cfg"
python python/plot_vacuum.py results/2026-08-29_p2a_vacuum_electrostatics
```

`plot_vacuum.py` bestimmt den Commit selbst aus git, im Augenblick des
Zeichnens, und beendet sich mit Code 1, wenn der Arbeitsbaum schmutzig ist oder
die Daten nicht aus HEAD stammen; die Abbildungen tragen dann quer über den
Kopf **NICHT FREIGEGEBEN … DIRTY**. Die frühere Fassung dieser Abbildungen nannte
einen Commit, der aus einem Argument stammte und die gezeigte Implementierung
gar nicht enthielt — dafür ist die Messung da.

Die verwendeten Konfigurationsdateien liegen als Kopie daneben
(`device_p1.cfg`, `vacuum_p2a.cfg`); die aufgelösten Werte in
`parameters.csv`, die Herkunft in `meta.txt` und `figures_provenance.txt`, der
volle Bericht in `report.txt`. `es_vacuum` gibt 2 zurück, wenn eine Netzprüfung
oder eine numerische Prüfung fehlschlägt — dieser Lauf gab 0 zurück.

## Der geschlossene Emitterleiter

Der Leiter aus Emittermetall und Flüssigkeitssäule endet nicht mehr offen am
unteren Modellschnitt. Er hat eine rückwärtige Länge `emitter_back_length` und
ist dort durch eine **volle leitende Scheibe** von `r = 0` bis zum Fußradius auf
`V_emitter` geschlossen. Randkennung: `numerical_emitter_back_closure`.

| | |
|---|---|
| Art des Abschlusses | ebene, voll leitende Kreisscheibe auf `V_emitter`, eigene Randkennung |
| kein Bauteil | die Scheibe steht für die reale rückwärtige Terminierung, die das Modell nicht kennt |
| kein Domänenrand | die offene Rechendomäne ist weiterhin an keiner Stelle BEM-Rand |
| Kante | als `kind=numerical_closure` markiert und aus jeder Feldauswertung ausgeschlossen |
| Abstand zur ausgewerteten Region | siehe `meta.txt`, `back_closure_clearance_m` |

`vacuum_bem_mesh()` **lehnt** ein Netz mit offenem Leiter jetzt ab, statt es zu
dokumentieren: `open_arc_ends()` findet jedes freie Blechende aus der
Panelkonnektivität, und ein einziges genügt für einen Abbruch mit Hinweis auf
den Parameter, der es schließt. Die frühere Begründung „der Hohlraum dahinter
ist abgeschirmt" begrenzt den Fehler, beseitigt ihn aber nicht — und sie hat
genau das künstliche Feldmaximum erzeugt, um dessen Beseitigung es hier geht.

Im Inneren des geschlossenen Leiters ist das Potential das angelegte und das
Feld verschwindet; was übrig bleibt, ist Diskretisierungsfehler und fällt mit
der Verfeinerung (Spalten `interior_dV_V` und `interior_Emag_V_per_m` in
`convergence.csv`).

## Die Trunkierungsstudie — und ihr Ergebnis

Vier zunehmende rückwärtige Längen (`truncation.csv`), Domäne und `size_scale`
fest gehalten, damit sich nur das hintere Ende ändert. Dass die Vernetzung an
Spitze und Extraktor dabei **bitweise** dieselbe bleibt, wird geprüft und in der
Spalte `front_mesh_identical` mitgeschrieben; andernfalls bricht der Lauf ab.

Vorab festgelegte Toleranzen (`es::truncation` in `include/es/vacuum_bem.hpp`,
begründet dort): `E_z(ref)` 1e-3, `c_EX` 1e-3, Sondenpotentiale 1e-3 von
`|V_E − V_X|`.

**Die Toleranzen werden nicht erreicht.** Gemessen bei der letzten Verdopplung:

| Größe | Änderung pro Verdopplung | Toleranz | Ergebnis |
|---|---:|---:|---|
| `E_z(ref)` | 2.5e-02 | 1e-03 | nicht erreicht |
| `c_EX` | 3.3e-01 | 1e-03 | nicht erreicht |
| Sondenpotentiale | 1.4e-02 | 1e-03 | nicht erreicht |
| `c_EE` | 4.1e-01 | — | hängt zwangsläufig von der Länge ab |

Über den ganzen vermeßbaren Bereich (200 µm bis 3.2 mm) trägt `E_z(ref)` einen
`1/L`-Schwanz — 5.6, 4.0, 2.5, 1.5 Prozent je Verdopplung —, und `c_EX` wächst
um rund ein Drittel je Verdopplung, ohne sich zu setzen.

Der Grund liegt in der Randbedingung, nicht in der Diskretisierung: mit
`V -> 0` im Unendlichen und nur zwei Leitern trägt das System Nettoladung. Es
gibt im Modell keine Rückelektrode und keine Umhausung, also trägt ein längerer
Emitter schlicht mehr Ladung, und die wird an der Spitze ebenso gespürt wie am
Extraktor. Ob das lokale Feld unempfindlich wird, sobald eine geerdete
Umhausung im Modell ist, gehört in die Phase, die eine hinzufügt.

**Folge für die Ergebnissemantik.** Es gibt in P2a keine truncation-konvergierte
Größe. `emitter_back_length` ist deshalb eine **Pflichtangabe der Geometrie**,
kein Konvergenzparameter, und der benutzte Wert (800 µm) ist ein
**Beispielwert**, keine gemessene Abmessung — genau wie
`extractor_outer_radius`. Jeder Zahlenwert in diesem Ordner gilt für diese
Länge; `quantities.csv` sagt für jede Größe einzeln, wogegen sie konvergiert ist
und wovon sie noch abhängt.

## Was wovon abhängt

| Größe | Netzkonvergenz | Trunkierung | hängt ab von |
|---|---|---|---|
| `E_z(ref)`, Sondenpotentiale | konvergiert (≈ 2e-05) | **nicht** konvergiert | `emitter_back_length`, `extractor_outer_radius`, alle Geräteabmessungen |
| `c_EX`, `C_m` | konvergiert (≈ 1e-05) | **nicht** konvergiert | dieselben |
| `c_EE` | konvergiert | nicht anwendbar | Leiterlänge in erster Ordnung — **keine Geräteaussage** |
| `c_XX` | konvergiert | nicht anwendbar | Elektrodenabmessungen |
| `\|E\|` an unverrundeten Kanten | **nicht** konvergent | — | Elementgröße; Verrundung ist P3 |
| `\|E\|` auf der Rückschließung | — | — | keine Gerätefläche, wird nicht berichtet |

`extractor_radius_study.csv` und Abbildung 3(c) zeigen die zweite reale
Geometrieabhängigkeit: zwischen 0.5 mm und 2.5 mm Elektrodenaußenradius ändert
sich `C_m` um mehr als den Faktor 2. Das ist kein Netzfehler.

## Der BEM-Kern und was er mathematisch kann

Aus `src/bem.cpp` abgelesen, nicht angenommen:

| Frage | Antwort |
|---|---|
| Formulierung | indirekte Einfachschicht; Unbekannte ist `sigma`, Kollokation in den Elementmitten |
| Randbedingung im Unendlichen | `V -> 0`, in der azimutal integrierten Freiraum-Greensfunktion enthalten; keine Abschneidefläche |
| geschlossene Leiter nötig? | für die Lösbarkeit nein, für die Auswertung ja — auf einer offenen Fläche ist `sigma` die **Summe** beider Seiten, `E_n = sigma/eps0` dort also kein einseitiges Feld. P2a verlangt deshalb geschlossene Leiter |
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
| `numerical_emitter_back_closure.liquid` / `.solid` | `numerical_emitter_back_closure` | ja → Emitter, **numerisch**, nicht auswertbar |
| `extractor_surface.aperture` / `.front` / `.back` / `.rim` | `extractor_surface` | ja → Extraktor |
| `symmetry_axis.*` | `symmetry_axis` | **nein** — kein Interface, keine Rotationsfläche |
| `bore_wall` | `bore_wall` | **nein** — innere Materialgrenze, kein Vakuumkontakt |
| `open_boundary.*` | `open_boundary` | **nein** — kein Leiter, keine Randbedingung |

Die liquide Hälfte der Abschlussscheibe trägt `Tag::Emitter`, nicht
`Tag::FreeSurface`: sonst stünde eine numerische Kappe in jeder Zahl, die für
die ebene Flüssigkeitsoberfläche berichtet wird. Im geschlossenen Aufbau gibt
es keinen `liquid_inlet`; ein hydraulischer Zulauf gehört zum Strömungsmodell
und kommt mit ihm zurück.

Die offene rechteckige Rechendomäne ist kein Leiter und wird an keiner Stelle
als BEM-Rand verwendet.

## Wo `sigma/eps0` ein physikalisches Feld ist

Nach der Schließung ist `sigma/eps0` das einseitige Vakuum-Normalfeld — aber
nur dort, wo die Fläche die geschlossene Perfect-Conductor-Region gegen Vakuum
begrenzt. `surface.csv` führt das je Panel in der Spalte `evaluable`. Nicht
auswertbar sind:

* die numerische Rückschließung (`numerical = 1`);
* der Schaft hinter dem Kegelfuß, dessen Länge numerisch gewählt ist;
* jedes Panel innerhalb einer markierten Kantenzone.

## Pflichtangabe `extractor_outer_radius`

Bisher bedeutete `0` „die Elektrode reicht bis zum Domänenrand". Das setzte
einen Leiter mit der offenen Kante der Rechenbox gleich, die keiner ist, und
versteckte damit eine fehlende physikalische Abmessung. Der Außenradius ist
Pflicht, `0` wird abgelehnt, und `domain_radius` muss **echt** größer sein. Der
Wert im Beispielparametersatz (2 mm) ist ausdrücklich ein Beispielwert.

## Kapazitätsgrößen — benannt, nicht „C"

Mit `V = 0` im Unendlichen hat ein Zweileitersystem eine 2×2-Maxwell-Matrix:

```
Q_E = c_EE V_E + c_EX V_X
Q_X = c_XE V_E + c_XX V_X
```

Die Zahlenwerte der Referenzstufe stehen in `report.txt` und `quantities.csv`.
`Q_E/(V_E − V_X)` ist bei geerdetem Extraktor identisch `c_EE` und **nicht**
`C_m`, weil das System bei `V(unendlich) = 0` Nettoladung trägt.

## Was hier ausdrücklich **keine** konvergierte Größe ist

Austrittskante, äußere Stirnkante und beide Aperturkanten sind unverrundet. Das
Feld folgt dort der Elementgröße. Die markierten Zonen (`edge_zones.csv`,
Radius = lokale Merkmalsgröße / 4, rein geometrisch) sind in Abbildung 2 weiß
ausgespart. Verrundungsradien sind ein P3-Parameter. Ein offenes Bogenende gibt
es nicht mehr; `edge_zones.csv` würde es als `kind=truncation_end` führen, und
der Lauf käme nicht bis zur Auswertung.

## Inhalt

| Datei | Inhalt |
|---|---|
| `report.txt` | vollständiger Textbericht des Laufs |
| `meta.txt` | Herkunft, Toleranzen, gemessene Trunkierungsänderungen, Urteil |
| `figures_provenance.txt` | Stempel der Abbildungen, Sauberkeit des Arbeitsbaums |
| `quantities.csv` | je Größe: Netzkonvergenz, Trunkierung, Abhängigkeiten |
| `truncation.csv` | vier rückwärtige Längen, alle ausgewerteten Größen |
| `probe_points.csv` | die festen Punkte zwischen Emitter und Extraktor |
| `bem_selection.csv` | Entscheidung je Randkurve mit Begründung |
| `convergence.csv` | alle vier Netzstufen, alle ausgewerteten Größen |
| `extractor_radius_study.csv` | Einfluss des Elektrodenaußenradius |
| `linearity.csv` | Skalierung von `sigma`, `V`, `E`, `Q` mit der Spannung |
| `polarity.csv` | vollständige Polaritätsumkehr |
| `surface.csv` | `sigma`, `E_n`, Kennung, Kantenzone, `numerical`, `evaluable` je Panel |
| `edge_zones.csv` | markierte Stellen, an denen kein Feldwert abgelesen werden darf |
| `axis_profile.csv` | `V` und `E_z` auf `r = 0` |
| `field_full.csv` | Feldgitter über die ganze Domäne |
| `field_tip.csv` | Feldgitter an der Emitterspitze |
| `field_aperture.csv` | Feldgitter an der Extraktoröffnung |
| `field_back_closure.csv` | Feldgitter an der geschlossenen numerischen Rückfläche |
| `parameters.csv`, `regions.csv`, `boundaries.csv`, `features.csv` | Geometrie als Hintergrund der Abbildungen |
| `device_p1.cfg`, `vacuum_p2a.cfg` | Kopien der verwendeten Eingaben |
| `fig1_potential.png` | Potentialkarte mit Equipotentiallinien und Rückschließung |
| `fig2_field_magnitude.png` | `\|E\|`-Karten, darunter die geschlossene Rückfläche |
| `fig3_convergence.png` | Netzkonvergenz, Elektrodenradius- und Trunkierungsstudie |
| `fig4_surface_charge.png` | Flächenladungsdichte entlang der Leiter |
