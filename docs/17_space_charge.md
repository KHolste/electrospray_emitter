# 17 — P6: Poisson mit freier Raumladung und die PIC-Grundlage

**Stand: 2026-08-29.** Status: **`validated_subset`**. Die Poisson-Gleichung mit
freier Raumladung ist implementiert und gegen eine hergestellte Lösung geprüft;
die Ladungsverteilung ist eine **vorgeschriebene Testquelle**, weil P5 blockiert
ist.

Code: `include/es/space_charge.hpp`, `src/space_charge.cpp`,
`apps/es_space_charge.cpp`, `tests/test_space_charge.cpp`,
`python/plot_space_charge.py`, plus die Erweiterung von `axisym_fem`.
Ergebnisse: `results/2026-08-29_p6_space_charge/`.

---

## 17.1 Kein neues Ringmodell

Das gesperrte Ringmakropartikelmodell des Prototyps divergierte, weil das
analytische Potential eines geladenen Rings **auf** dem Ring logarithmisch
singulär ist: ein Teilchen auf seinem eigenen Ring sah ein unendliches Feld, und
Verfeinerung machte es schlimmer.

**Hier wird die Ladung gar nicht als Singularität in die Feldlösung gesteckt.**
Sie wird mit den Elementformfunktionen auf die Knoten **deponiert** und geht in
den FEM-Lastvektor,

$$
f_a = \sum_p q_p\,N_a(\mathbf x_p),
$$

was genau die schwache Form von $\nabla\cdot(\varepsilon\nabla\varphi)=-\rho$
für $\rho=\sum_p q_p\delta(\mathbf x-\mathbf x_p)$ ist. Die Lösung dieses
diskreten Systems ist eine **stückweise bilineare Funktion**: endlich überall,
mit endlichem Gradienten. Es gibt keine Singularität, die regularisiert werden
müsste, weil die Diskretisierung nie eine erzeugt.

**Kein $2\pi r$ am Lastvektor**, und das ist keine Auslassung: das
Rotationsvolumenelement steckt bereits in dem Integral, aus dem die Ladung
hervorging. Es ein zweites Mal anzuwenden hieße, jedes Makropartikel mit seinem
eigenen Radius zu multiplizieren.

**Der Preis, ausgesprochen statt versteckt.** Das Selbstfeld eines
Makropartikels ist das Feld einer über die Nachbarzellen verschmierten Ladung
und hängt damit vom Netz ab; sein Spitzenwert **wächst** unter Verfeinerung
(gemessen: 1,94 → 2,22 → 2,51 V über drei Stufen). Eine PIC-Schleife darf ein
Teilchen sein eigenes deponiertes Feld deshalb nicht ungefiltert fühlen lassen.
Dieses Modul implementiert keine solche Korrektur; es macht den Effekt messbar
und sagt das.

## 17.2 Was `axisym_fem` dazu bekommen hat

Drei Quellterme, alle mit klar getrennten Einheiten:

| Feld | integriert als | Einheit von φ |
|---|---|---|
| `cell_source` | $\int s\,N_a\,2\pi r\,\mathrm dA$, $s$ zellweise konstant | je nach Problem |
| `node_source_density` | $\int\rho\,N_a\,2\pi r\,\mathrm dA$, $\rho$ knotenweise bilinear | [V] mit $\rho$ in C/m³ |
| `node_charge` | **direkt** in $f_a$, ohne $2\pi r$ | [V] mit $q$ in C |

Die Vorgabe aller drei ist leer, also ändert sich für keinen bestehenden
Aufrufer etwas.

## 17.3 Die sechs Pflichtprüfungen

| # | Prüfung | Ergebnis |
|---|---|---|
| 1 | $\rho=0$ reproduziert die Laplace-Lösung | **bitgleich** (0,000e+00), Feld ebenfalls |
| 2 | hergestellte Poisson-Lösung | Potential 4,5·10⁻⁵ relativ bei 161×321, Ordnung **1,83** |
| 3 | deponierte Gesamtladung = Summe der Makropartikelladungen | 2,4·10⁻¹⁶ relativ |
| 4 | keine Divergenz beim Annähern | Potential und Feld bleiben beschränkt über acht Halbierungen des Abstands |
| 5 | Netzkonvergenz von Potential und Feld | siehe 17.4 |
| 6 | Polaritätsvorzeichen | positive Ladung hebt das Potential, negative senkt es; die Lösung ist **exakt ungerade** in der Ladung (10⁻¹²) |

### Die hergestellte Lösung

$$
\varphi = \varphi_0\frac{(R^2-r^2)(L^2-z^2)}{R^2L^2},\qquad
\rho = -\varepsilon_0\nabla^2\varphi
= \frac{2\varepsilon_0\varphi_0}{R^2L^2}\bigl[2(L^2-z^2)+(R^2-r^2)\bigr].
$$

Sie verschwindet auf $r=R$ und $z=\pm L$ — genau dort, wo die Testbox geerdet
ist —, ist auf der Achse regulär, und beide Ausdrücke sind unabhängig
hingeschrieben. Der Test vergleicht den Löser also mit einer geschlossenen Form
und nicht mit sich selbst.

### Selbstfeld und Elektrodenfeld getrennt

Beide werden einzeln gelöst und ihre Überlagerung geprüft: das Problem ist
linear, und die Diskretisierung bricht das nicht.

$$
\lvert\varphi_\text{beide}-(\varphi_\text{selbst}+\varphi_\text{Elektroden})\rvert
= 2{,}6\cdot10^{-13}\ \text{V von } 250\ \text{V}.
$$

Außerdem: die Verschiebung gegen $\rho=0$ ist mit und ohne Elektrodenspannung
**dieselbe Zahl** (auf $2\cdot10^{-16}$) — sie hängt nur an der Ladung.

## 17.4 Ein gemessener Befund über die vorhandene Feldrekonstruktion

Die Netzkonvergenz des **Feldes** wurde zunächst als ein Wert für das ganze
Gebiet gemessen und ergab Ordnung 1,00 — obwohl `field_recovered_at_node()` als
superkonvergent dokumentiert ist. Aufgeteilt nach Radius:

| Größe | Ordnung |
|---|---|
| Potential | 1,83 |
| Feld, $r>R/4$ | **2,00** |
| Feld, $r<R/4$ | **1,00** |

**Die Ursache.** Die Rekonstruktion mittelt die superkonvergenten
Zellmittelpunkt-Gradienten mit dem **Zellvolumen** als Gewicht. In
Achsensymmetrie trägt dieses Gewicht einen Faktor $2\pi r$, so dass die beiden
Nachbarzellen eines Knotens bei Radius $r$ im Verhältnis
$(r+h/2)/(r-h/2)$ gewichtet werden — eine Unsymmetrie der Ordnung $h/r$. Weit
von der Achse verschwindet sie und die Rekonstruktion ist zweiter Ordnung;
achsennah dominiert sie und die Ordnung fällt auf eins.

**Das betrifft mehr als P6.** Dieselbe Rekonstruktion liefert jedes
Oberflächenfeld von P2b, P2c und P3b. Es entwertet diese Ergebnisse nicht —
deren Größen sind Integrale oder liegen fern der Achse —, aber ein Feld erster
Ordnung ist das, was eine PIC-Schleife abtasten würde, und das muss man wissen,
bevor man es benutzt.

Gefunden wurde es, indem der Fehler **nach Radius getrennt** statt als eine Zahl
für das ganze Gebiet berichtet wurde.

## 17.5 Was P6 ausdrücklich nicht enthält

* **Keine Emission** und keine selbstkonsistente Emissions-PIC-Schleife. P5 ist
  blockiert; jede hier verwendete Ladungsverteilung ist **vorgeschrieben** und
  heißt in jeder Ausgabe Testquelle.
* **Keine Teilchenbewegung** — das ist P7.
* Keine Magnetfelder, keine Stöße, keine raumladungsbegrenzte Emission.
* Keine Korrektur des Selbstfeldes eines Makropartikels.
