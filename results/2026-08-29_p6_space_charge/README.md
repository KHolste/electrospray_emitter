# P6 — Poisson mit freier Raumladung — 2026-08-29

Status: **`validated_subset`**.

```sh
./build/es_space_charge results/2026-08-29_p6_space_charge meta.commit=$(git rev-parse HEAD)
python python/plot_space_charge.py results/2026-08-29_p6_space_charge
```

Modellvertrag und Befunde: [`docs/17_space_charge.md`](../../docs/17_space_charge.md).

## Die Ladung ist eine VORGESCHRIEBENE Testquelle

P5 ist blockiert; es gibt keine physikalische Teilchenquelle. Aus keiner Zahl
dieses Laufs darf ein emittierter Strom gelesen werden.

## Kein neues Ringmodell

Das gesperrte Ringmodell divergierte, weil das analytische Potential eines
geladenen Rings auf dem Ring logarithmisch singulaer ist. Hier wird die Ladung
gar nicht als Singularitaet in die Feldloesung gesteckt: sie wird mit den
Formfunktionen auf die Knoten deponiert und geht als f_a = sum_p q_p N_a(x_p)
in den FEM-Lastvektor -- genau die schwache Form. Die Loesung ist stueckweise
bilinear, also endlich ueberall.

**Kein 2 pi r am Lastvektor**, und das ist keine Auslassung: das
Rotationsvolumenelement steckt bereits in dem Integral, aus dem die Ladung kam.

## Die sechs Pflichtpruefungen

| # | Pruefung | Ergebnis |
|---|---|---|
| 1 | rho = 0 reproduziert die Laplace-Loesung | **bitgleich** (0,000e+00), Feld ebenfalls |
| 2 | hergestellte Loesung phi = phi0 (R²−r²)(L²−z²)/(R²L²) | 4,5e-05 relativ, Ordnung 1,83 |
| 3 | deponierte Ladung = Summe der Makropartikelladungen | 2,4e-16 relativ |
| 4 | Annaeherung **bei festem Netz** (81x161) | beschraenkt ueber acht Halbierungen des Abstands -- und ueber h -> 0 sagt das nichts, siehe unten |
| 5 | Netzkonvergenz von Potential und Feld | siehe unten |
| 6 | Polaritaet | Vorzeichen richtig; Loesung **exakt ungerade** in der Ladung (1e-12) |

Zusaetzlich: Selbstfeld und Elektrodenfeld getrennt geloest, Ueberlagerung exakt
(2,6e-13 V von 250 V). Ein Teilchen ausserhalb des Netzes wird GEZAEHLT und
nicht deponiert -- die Erhaltungspruefung schlaegt dann aus, statt zu schweigen.

## Das Selbstfeld: die Deposition regularisiert NICHTS

Eine fruehere Fassung berichtete hier *„keine Divergenz beim Annaehern“*. Das
ist wahr und es ist eine Aussage ueber **ein festes Netz**. Ueber h -> 0 sagt es
nichts, und dort ist die Antwort fuer jede der fuenf beteiligten Groessen eine
andere. Sie werden jetzt getrennt gemessen.

| | Frage | bei festem Netz | bei h -> 0 |
|---|---|---|---|
| (a) | Ladungserhaltung der Deposition | exakt (2,4e-16) | bleibt exakt auf **jeder** Stufe |
| (b) | Fremdfeld bei **festem** Abstand | endlich | **konvergiert**: Ordnung 2,05 (phi) / 1,91 (E), letzte relative Aenderung 4e-05 |
| (c) | Selbstpotential am Teilchen | endlich | **divergiert** -- siehe unten |
| (d) | Breite der deponierten Wolke | eine Zelle | faellt **wie das Netz** gegen null |
| (e) | Anteil des Selbstfeldes | — | faellt asymptotisch wie N^-0,82 mit der Makropartikelzahl |

### (c) Das Selbstpotential divergiert, und das Gesetz haengt vom Ort ab

| Lage | ueber fuenf Stufen (h: 1e-6 → 6,25e-8 m) | angepasstes Gesetz |
|---|---|---|
| abseits der Achse, r = 0,35 R | 1,938 → 3,073 V, Faktor **1,59** | 0,409·ln(1/h), Rest 1,2e-04 |
| auf der Achse, r = 0 | 32,4 → 524,6 V, Faktor **16,18** | h^-1,004, Rest 3,0e-04 |

Beide Gesetze werden angepasst und beide Reste berichtet. Das Ergebnis ist genau
das physikalisch Erwartete: abseits der Achse ist ein Makropartikel ein **Ring**,
dessen exaktes Potential auf dem Ring logarithmisch singulaer ist -- bei h
abgeschnitten waechst es wie ln(1/h). Auf der Achse entartet der Ring zur
**Punktladung** mit 1/d-Singularitaet, also wie 1/h.

Der fruehere Header-Kommentar behauptete, es muesse „wie 1/h“ wachsen. Das ist
**nur** das Achsengesetz; abseits der Achse liefert der Potenzfit mit 0,166 einen
deutlich schlechteren Rest als der Logarithmus.

**Damit ist die Deposition keine netzunabhaengige Regularisierung.**

## Die drei Kandidaten, mit Urteil und mit der Messung dahinter

| Variante | Urteil | Begruendung |
|---|---|---|
| **Selbstfeldabzug** | `implemented` | Eine Ladung uebt auf sich selbst keine Kraft aus; das deponierte Selbstfeld ist reines Artefakt und traegt eine Richtung, die allein von der Lage im Zellinneren abhaengt. Die diskrete Aufgabe ist **linear**, also exakt abziehbar statt daempfbar. |
| feste physikalische Formbreite | `rejected_free_parameter` | Es gibt keine physikalische Laenge, die sie festlegte -- Messung (d) zeigt, dass in der Wolke nur h steckt. Sie waere ein **frei gewaehlter** Glaettungsparameter. |
| skalierte Makropartikelzahl | `measured_not_implemented` | Tragfaehig und der uebliche PIC-Weg, aber eine Eigenschaft einer **Schleife**, die es hier nicht gibt. Gemessen, nicht behauptet. |

Die Tabelle liegt als `es::pic_options()` in der Bibliothek und wird in
`tests/test_space_charge.cpp`, Abschnitt 9, geprueft -- unter anderem darauf,
dass **genau eine** Variante implementiert ist.

**Was der Abzug leistet.** Ein *einzelnes* Teilchen im sonst leeren, geerdeten
Kasten spuert ohne Abzug 1,13e5 bis 1,42e5 V/m, obwohl nichts anderes darin ist;
der Betrag haengt allein von seiner Lage in der Zelle ab. Mit Abzug spuert es
**exakt null**. Die Linearitaet wird gegen die volle Loesung geprueft:
Ueberlagerungsfehler 1,6e-15.

**Der Preis:** eine zusaetzliche Loesung je Teilchen. Fuer eine Diagnose richtig,
fuer eine Produktionsschleife untragbar.

## Die PIC-Schleife bleibt blockiert -- aus zwei unabhaengigen Gruenden

1. **Keine Quelle.** P5 ist blockiert; es gibt keine belegte Emissionsrate.
2. **Keine tragfaehige Kostenstruktur.** Der Abzug kostet eine Loesung je
   Teilchen.

Beide sind unabhaengig: faellt einer weg, bleibt der andere. P6 bleibt damit ein
**validierter Poisson- und Depositions-Teiltest** mit einer exakt begruendeten
Selbstfeldbehandlung fuer die Diagnose -- und keine PIC-Schleife.

## Der zweite Befund, der bleibt

**Die vorhandene Feldrekonstruktion ist achsennah nur ERSTER Ordnung.**

| Groesse | Ordnung |
|---|---|
| Potential | 1,83 |
| Feld, r > R/4 | **2,00** |
| Feld, r < R/4 | **1,00** |

Ursache: `field_recovered_at_node()` mittelt die superkonvergenten
Zellmittelpunkt-Gradienten mit dem ZELLVOLUMEN als Gewicht, und dieses traegt in
Achsensymmetrie einen Faktor 2 pi r. Die beiden Nachbarzellen eines Knotens bei
Radius r werden damit im Verhaeltnis (r+h/2)/(r-h/2) gewichtet -- eine
Unsymmetrie der Ordnung h/r.

Das betrifft mehr als P6: dieselbe Rekonstruktion liefert jedes Oberflaechenfeld
von P2b, P2c und P3b. Es entwertet diese Ergebnisse nicht, deren Groessen sind
Integrale oder liegen fern der Achse -- aber ein Feld erster Ordnung ist das, was
eine PIC-Schleife abtasten wuerde.

Gefunden wurde es, indem der Fehler NACH RADIUS getrennt berichtet wurde statt
als eine Zahl fuer das ganze Gebiet.

## NICHT gerechnet

Keine Emission, keine selbstkonsistente PIC-Schleife, keine Teilchenbewegung
(P7), kein Ringmodell, keine Magnetfelder, keine Stoesse, keine
raumladungsbegrenzte Emission.

Und ausdruecklich: **keine netzunabhaengige Regularisierung des Selbstfeldes.**
Der exakte Abzug entfernt die scheinbare Selbstkraft vollstaendig, macht das
Selbstpotential aber nicht netzunabhaengig -- es gibt es danach am Teilchen
schlicht nicht mehr. Fuer alles andere bleibt das deponierte Feld eine
Netzgroesse. Keine skalierte Makropartikelzahl und damit keine
PIC-Konvergenzstrategie: gemessen, nicht implementiert.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_space_charge.png` | deponierte Knotenladung, Potential mit Ladung, die Potentialaenderung gegen rho = 0 und die Feldaenderung |
| `fig2_checks.png` | Netzkonvergenz gegen die hergestellte Loesung nach Radius getrennt, Ladungserhaltung und Selbstfeldskalierung, und die Annaeherung an ein einzelnes Makropartikel **bei festem Netz** |
| `fig3_self_field.png` | was bei festem Netz beschraenkt ist, was bei h -> 0 divergiert (beide Lagen, mit angepasstem Gesetz), was bei h -> 0 konvergiert, und der Abzug, der die scheinbare Selbstkraft auf exakt null bringt |

## Dateien

`self_field_scaling.csv` -- das Selbstpotential ueber fuenf Netzstufen, auf und
abseits der Achse, mit beiden angepassten Gesetzen und beiden Resten.
`foreign_field.csv` -- das Fremdfeld bei festem Abstand und seine Konvergenzordnung.
`deposition_width.csv` -- die Breite der deponierten Wolke gegen die Zellgroesse.
`self_to_total.csv` -- der Selbstanteil gegen die Makropartikelzahl.
`self_field_exclusion.csv` -- die scheinbare Selbstkraft mit und ohne Abzug.
`pic_options.csv` -- die drei Kandidaten mit Urteil, Begruendung und Messung,
und der Status der PIC-Schleife.
