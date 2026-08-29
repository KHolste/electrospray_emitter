# Zweig-Mehrdeutigkeit des Spannungs-Solvers — 2026-08-29

Ergebnisordner der P0-Nachprüfung. **Alle Daten sind in diesem Lauf frisch
berechnet**; keine älteren CSV-Dateien wurden wiederverwendet.

## Reproduktion

```sh
cmake --build build
./build/branch_figure.exe results/2026-08-29_branch_ambiguity
python python/plot_branch.py results/2026-08-29_branch_ambiguity
```

`tools/branch_figure.cpp` rechnet den Ast und beide Lösungen aus den
Quellparametern neu; `python/plot_branch.py` liest ausschließlich die hier
abgelegten CSV-Dateien.

## Inhalt

| Datei | Inhalt |
|---|---|
| `branch_full.csv` | Ast über den vollen angeforderten Bereich (h_max = 2,2 r_c) |
| `branch_before_fold.csv` | Ast, Bereich endet vor dem Umkehrpunkt |
| `branch_past_fold.csv` | Ast, Bereich endet dahinter, aber über der Zielspannung |
| `shape_lower_height.csv` | Meniskusprofil der LowerHeight-Lösung |
| `shape_upper_height.csv` | Meniskusprofil der UpperHeight-Lösung |
| `mesh_electrodes.csv` | Emitter- und Extraktorgeometrie |
| `markers.csv` | Zielspannung, Umkehrpunkt-Kandidat, beide Schnittpunkte, Bereichsenden |
| `meta.txt` | Parametersatz und Zahlenwerte im Klartext |
| `fig1_branch_U_vs_h.png` | U(h) mit Zielspannung, beiden Schnittpunkten, Kandidat und Bereichsenden |
| `fig2_meniscus_profiles.png` | beide Meniskusprofile bei 1154,2 V mit Emitter und Extraktor |

## Zur Lesart

`LowerHeight` und `UpperHeight` bezeichnen **ausschließlich die Apexhöhe**.
Keine der beiden Lösungen ist als stabil oder instabil gekennzeichnet — eine
Stabilitätsanalyse ist nicht implementiert.

Der markierte Umkehrpunkt ist ein **Kandidat**, ermittelt als diskretes Maximum
der abgetasteten Astpunkte. Daraus folgt weder dynamische Stabilität noch ein
Emissionsbeginn noch ein Cone-Jet-Übergang.

Abbildung 1a zeigt zusätzlich den Unterschied zwischen dem *angeforderten*
Bereich (h_max = 2,2 r_c) und dem *tatsächlich verfolgten*: die Fortsetzung
bricht bei h/r_c = 0,80 ab. Für die Zielspannung 1154,2 V liegen beide
Schnittpunkte dennoch innerhalb des verfolgten Bereichs, weshalb die Abdeckung
für genau diese Zielspannung als vollständig gilt.
