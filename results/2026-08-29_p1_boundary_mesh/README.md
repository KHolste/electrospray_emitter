# Automatischer achsensymmetrischer Randvernetzer — 2026-08-29

Ergebnisordner der zweiten P1-Teilstufe: **Randvernetzung**.
Alle Daten in diesem Lauf frisch erzeugt.

**Keine Volumenvernetzung, keine Feldlösung, keine Meniskusrechnung, keine
Strömung, keine Emission, keine Raumladung.** Die Randkennungen sind
Bezeichner, keine Randbedingungen. Die P0-BEM benutzt dieses Netz nicht; sie
behält ihre eigenen Netze.

## Reproduktion

```sh
cmake -S . -B build -G Ninja
cmake --build build
./build/mesh_figure.exe examples/device_p1.cfg results/2026-08-29_p1_boundary_mesh
python python/plot_mesh.py results/2026-08-29_p1_boundary_mesh
```

Der Parametersatz steht vollständig in `../../examples/device_p1.cfg` (der
gleiche wie für die Geometriestufe); die aufgelösten Werte in
`parameters.csv`, das Netz und sein Prüfbericht in `mesh_report.txt`.

`mesh_figure` gibt 2 zurück, wenn eine Netzprüfung fehlschlägt — dieser Lauf
gab 0 zurück, alle 18 Prüfungen bestanden.

## Inhalt

| Datei | Inhalt |
|---|---|
| `mesh_nodes.csv` | Netzknoten (r, z), Achsen- und Eckenmarkierung, Merkmalszuordnung |
| `mesh_elements.csv` | Randelemente: Kennung, Kurve, Endpunkte, Normale, Nachbargebiete, Meridianlänge, Rotationsfläche |
| `mesh_boundaries.csv` | Elementstatistik je Randkurve, je Randkennung und gesamt |
| `mesh_size_field.csv` | Verfeinerungsquellen der Größenfunktion mit lokaler Merkmalsgröße |
| `mesh_report.txt` | Klartextbericht: Größenfunktion, Statistik, alle Prüfungen |
| `regions.csv`, `boundaries.csv`, `features.csv`, `parameters.csv` | Geometrie (Hintergrund der Abbildungen), aus P1a |
| `fig1_mesh_overview.png` | vollständiges Randnetz, maßstäblich, mit Statistiktafel |
| `fig2_mesh_tip_detail.png` | Bohrung, ebene Flüssigkeitsoberfläche, Austrittskante |
| `fig3_mesh_aperture.png` | Extraktoröffnung und Aperturkanten |

## Elementstatistik (Längen in µm)

| Rand | n | min | median | max |
|---|---:|---:|---:|---:|
| `symmetry_axis` | 78 | 0.692 | 9.007 | 86.789 |
| `emitter_outer_surface` | 49 | 0.170 | 4.381 | 39.745 |
| `emitter_tip_land` | 13 | 0.176 | 0.369 | 0.728 |
| `bore_wall` | 48 | 0.173 | 5.450 | 40.228 |
| `free_surface_reference` | 10 | 0.166 | 0.483 | 0.864 |
| `liquid_inlet` | 6 | 0.650 | 0.819 | 1.031 |
| `extractor_surface` | 105 | 3.504 | 68.002 | 87.759 |
| `open_boundary` | 125 | 0.664 | 87.005 | 88.235 |
| **gesamt** | **434** | **0.166** | **18.393** | **88.235** |

431 Knoten, größtes Verhältnis benachbarter Elementgrößen 1.283
(Schranke 1.5). Die Aufschlüsselung je einzelner Randkurve steht in
`mesh_boundaries.csv` und `mesh_report.txt`.

## Zur Lesart

**Die Fläche bei z = 0 ist die anfängliche ebene Flüssigkeitsoberfläche —
noch kein berechneter Meniskus.** Sie ist der geometrische Ausgangszustand.
Daraus darf keine physikalische Meniskuslösung abgeleitet werden; der Meniskus
gehört zu P3.

**Achsensymmetrisch, nicht planar.** Jedes Element ist eine Strecke in der
Meridianhalbebene; das Objekt, das es darstellt, ist die durch Rotation um
r = 0 erzeugte Fläche. `meridian_length` und `revolved_area` heißen deshalb
verschieden.

**Elemente auf r = 0 sind Symmetrieelemente, keine Ringelemente.** Ihre
Rotationsfläche ist exakt null, nicht klein. Sie werden in
`mesh_elements.csv` mit `kind = axis_symmetry` geführt.

**Keine Netzparameter.** Die Elementgröße folgt allein aus der Geometrie:
`h(x) = min(min_s [h_s + G·|x − x_s|], h_max)` mit Quellen auf den benannten
Merkmalen (`h = lfs/32`) und den übrigen Geometrieecken (`h = lfs/8`),
`G = 0.25`, `h_max = Domänendiagonale/40`. `lfs` ist die lokale
Merkmalsgröße. `h_tip` und `h_far` existieren nicht.

**Scharfe Kanten sind verfeinert — das macht das dortige Feld nicht
konvergent.** Das Eckfeld einer unverrundeten Kante divergiert und folgt der
örtlichen Elementgröße. Ein maximales elektrisches Feld an der Austrittskante
oder an einer Aperturkante darf später nicht als netzkonvergente Größe
behandelt werden. Endliche Kantenradien sind ein P3-Parameter.

**Die Maßstäbe der drei Abbildungen sind verschieden** (mm in Abb. 1, µm in
Abb. 2 und 3; die beiden Ausschnitte je Detailabbildung noch einmal
unterschiedlich). Das ist in den Titeln vermerkt.
