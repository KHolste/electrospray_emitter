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
| 4 | keine Divergenz beim Annaehern | beschraenkt ueber acht Halbierungen des Abstands |
| 5 | Netzkonvergenz von Potential und Feld | siehe unten |
| 6 | Polaritaet | Vorzeichen richtig; Loesung **exakt ungerade** in der Ladung (1e-12) |

Zusaetzlich: Selbstfeld und Elektrodenfeld getrennt geloest, Ueberlagerung exakt
(2,6e-13 V von 250 V). Ein Teilchen ausserhalb des Netzes wird GEZAEHLT und
nicht deponiert -- die Erhaltungspruefung schlaegt dann aus, statt zu schweigen.

## Zwei Befunde, die man behalten sollte

**1. Das Selbstfeld eines Makropartikels ist eine NETZGROESSE.** Sein
Spitzenpotential waechst unter Verfeinerung (1,94 -> 2,22 -> 2,51 V ueber drei
Stufen), weil ein feineres Netz die Ladung weniger verschmiert. Eine
PIC-Schleife darf ein Teilchen sein eigenes deponiertes Feld deshalb nicht
ungefiltert fuehlen lassen. Dieses Modul implementiert keine solche Korrektur;
es macht den Effekt messbar.

**2. Die vorhandene Feldrekonstruktion ist achsennah nur ERSTER Ordnung.**

| Groesse | Ordnung |
|---|---|
| Potential | 1,83 |
| Feld, r > R/4 | **2,00** |
| Feld, r < R/4 | **1,00** |

Ursache: `field_recovered_at_node()` mittelt die superkonvergenten
Zellmittelpunkt-Gradienten mit dem ZELLVOLUMEN als Gewicht, und dieses traegt in
Achsensymmetrie einen Faktor 2 pi r. Die beiden Nachbarzellen eines Knotens bei
Radius r werden damit im Verhaeltnis (r+h/2)/(r−h/2) gewichtet -- eine
Unsymmetrie der Ordnung h/r.

Das betrifft mehr als P6: dieselbe Rekonstruktion liefert jedes
Oberflaechenfeld von P2b, P2c und P3b. Es entwertet diese Ergebnisse nicht,
deren Groessen sind Integrale oder liegen fern der Achse -- aber ein Feld
erster Ordnung ist das, was eine PIC-Schleife abtasten wuerde.

Gefunden wurde es, indem der Fehler NACH RADIUS getrennt berichtet wurde statt
als eine Zahl fuer das ganze Gebiet.

## NICHT gerechnet

Keine Emission, keine selbstkonsistente PIC-Schleife, keine Teilchenbewegung
(P7), kein Ringmodell, keine Magnetfelder, keine Stoesse, keine
raumladungsbegrenzte Emission, keine Selbstfeldkorrektur.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_space_charge.png` | deponierte Knotenladung, Potential mit Ladung, die Potentialaenderung gegen rho = 0 und die Feldaenderung |
| `fig2_checks.png` | Netzkonvergenz gegen die hergestellte Loesung nach Radius getrennt, Ladungserhaltung und Selbstfeldskalierung, und die Annaeherung an ein einzelnes Makropartikel |
