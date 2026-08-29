# Parametrische P1-Geometrie — 2026-08-29

Ergebnisordner der ersten P1-Teilstufe: **Geometrie und Randtopologie**.
Keine Physik, keine Randbedingungen, kein Netz. Alle Daten in diesem Lauf
frisch erzeugt.

## Reproduktion

```sh
cmake --build build
./build/device_figure.exe examples/device_p1.cfg results/2026-08-29_p1_device_geometry
python python/plot_device.py results/2026-08-29_p1_device_geometry
```

Der Parametersatz steht vollständig in `../../examples/device_p1.cfg`; die
aufgelösten Werte und alle abgeleiteten Maße in `meta.txt`.

## Inhalt

| Datei | Inhalt |
|---|---|
| `parameters.csv` | aufgelöster Parametersatz in SI, mit Einheit und Rolle |
| `regions.csv` | geschlossene Meridiankonturen der vier Gebiete (Loop 0 außen, >0 Löcher) |
| `boundaries.csv` | benannte Randpolylinien mit Kennung und angrenzenden Gebieten |
| `features.csv` | nulldimensionale Merkmale (Kreise in 3D), die ein Vernetzer auflösen muss |
| `meta.txt` | Klartextbericht: Parameter, Gebietsvolumina, Randflächen, Merkmale |
| `fig1_device_overview.png` | maßstäbliche Übersicht der r-z-Domäne |
| `fig2_device_detail.png` | Detail: Bohrung/Stirnfläche/Austrittskante und Extraktoröffnung |

## Zur Lesart

**Achsensymmetrisch, nicht planar.** Jede Fläche und jedes Volumen entsteht
durch Rotation um r = 0 und trägt den Faktor 2πr. Die Größen heißen deshalb
verschieden: `meridian_length`, `meridian_area` (in der r-z-Ebene) gegenüber
`revolved_area`, `revolved_volume` (nach Rotation).

Das äußere Rechteck ist die **offene Rechendomäne**, kein Leiter und kein
Gehäuse auf festem Potential.

`free_surface_reference` ist die ebene Referenzfläche am Bohrungsaustritt,
**nicht** ein gerechneter Meniskus. Der Meniskus gehört zu P3.

Die Randkennungen sind **Bezeichner**. Es ist in dieser Stufe keine
Randbedingung zugeordnet und keine Feldgleichung gelöst.
