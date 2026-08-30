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

**Der Preis, ausgesprochen statt versteckt — und diesmal richtig.** Die
Finite-Elemente-Lösung ist **auf einem festen Netz** beschränkt. Sie ist damit
**nicht regularisiert**: für $h\to0$ wächst das Selbstpotential eines
Makropartikels unbeschränkt. Siehe 17.4a, wo das getrennt gemessen wird.

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
| 4 | Annäherung **bei festem Netz** | Potential und Feld bleiben über acht Halbierungen des Abstands beschränkt — auf dem Netz 81×161 und über $h\to0$ sagt das nichts, siehe 17.4a |
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

## 17.4a Das Selbstfeld, in fünf getrennte Fragen zerlegt

Eine frühere Fassung dieses Abschnitts berichtete *„keine Divergenz beim
Annähern — beschränkt über acht Halbierungen des Abstands“*. Das ist wahr und
es ist eine Aussage über **ein festes Netz**. Über $h\to0$ sagt es nichts, und
dort ist die Antwort für jede der folgenden fünf Größen eine andere. Deshalb
werden sie getrennt gemessen.

| | Frage | Antwort bei festem Netz | Antwort bei $h\to0$ |
|---|---|---|---|
| (a) | Ladungserhaltung der Deposition | exakt bis auf Rundung ($2\cdot10^{-16}$) | bleibt exakt auf **jeder** Stufe |
| (b) | Fremdfeld bei **festem** Abstand | endlich | **konvergiert**, Ordnung 2,05 (φ) bzw. 1,91 (E), letzte relative Änderung $4\cdot10^{-5}$ |
| (c) | Selbstpotential am Teilchen | endlich | **divergiert** — siehe unten |
| (d) | Breite der deponierten Wolke | eine Zelle breit | $\to0$ **wie das Netz**; keine andere Länge steckt darin |
| (e) | Anteil des Selbstfeldes | — | fällt mit der Makropartikelzahl wie $N^{-0{,}82}$ (asymptotisch) |

### (c) Das Selbstpotential divergiert — und das Gesetz hängt vom Ort ab

| Lage | Wachstum über fünf Stufen (h: 1e-6 → 6,25e-8 m) | angepasstes Gesetz |
|---|---|---|
| abseits der Achse, $r=0{,}35R$ | 1,938 → 3,073 V, Faktor **1,59** | $0{,}409\cdot\ln(1/h)$, Rest $1{,}2\cdot10^{-4}$ |
| auf der Achse, $r=0$ | 32,4 → 524,6 V, Faktor **16,18** | $h^{-1{,}004}$, Rest $3{,}0\cdot10^{-4}$ |

Beide Gesetze werden angepasst und beide Reste berichtet, damit die Antwort aus
den Daten kommt und nicht aus der Erwartung. Sie ist physikalisch genau die,
die man erwarten muss:

* **Abseits der Achse ist ein Makropartikel ein RING.** Das exakte Potential
  eines geladenen Rings ist auf dem Ring logarithmisch singulär. Bei $h$
  abgeschnitten wächst es wie $\ln(1/h)$ — langsam, aber ohne Grenze.
* **Auf der Achse entartet der Ring zur Punktladung** mit $1/d$-Singularität.
  Bei $h$ abgeschnitten wächst es wie $1/h$; der gemessene Exponent 1,004 ist
  genau das.

Eine frühere Fassung des Header-Kommentars behauptete, das Selbstpotential
müsse „wie $1/h$“ wachsen. Das ist nur das **Achsen**gesetz; abseits der Achse
zeigt die Messung den Logarithmus, und der Potenzfit liefert dort mit 0,166
einen deutlich schlechteren Rest.

**Damit ist die Deposition keine netzunabhängige Regularisierung**, und keine
Formulierung in diesem Projekt darf mehr so klingen.

### (d) Warum eine „feste physikalische Formbreite“ ausscheidet

Die deponierte Wolke eines Teilchens erreicht höchstens die vier Knoten seiner
Zelle; ihr größter Abstand bleibt unter der Zelldiagonale, und ihre
ladungsgewichtete RMS-Breite fällt mit dem Netz gegen null. **In ihr steckt
keine andere Länge als $h$.** Es gibt also nichts, woraus eine physikalische
Breite folgen würde: sie wäre ein frei gewählter Glättungsparameter, der genau
die Größe festlegte, die er verbergen soll. Beschränktheit, die man sich
aussucht, ist keine Konvergenz.

## 17.4b Die drei Kandidaten, mit Urteil und mit der Messung dahinter

| Variante | Urteil | Begründung |
|---|---|---|
| **Selbstfeldabzug** | `implemented` | Physikalisch übt eine Ladung auf sich selbst keine Kraft aus; das deponierte Selbstfeld ist reines Diskretisierungsartefakt und trägt zudem eine Richtung, die allein von der Lage im Zellinneren abhängt. Die diskrete Aufgabe ist **linear**, also lässt sich das Artefakt **exakt abziehen** statt dämpfen — keine Glättungsbreite, kein Filter, kein freier Parameter. |
| feste physikalische Formbreite | `rejected_free_parameter` | Es gibt keine physikalische Länge, die sie festlegte (siehe 17.4a (d)). Sie wäre ein frei gewählter Parameter. |
| skalierte Makropartikelzahl | `measured_not_implemented` | Mathematisch tragfähig und der übliche PIC-Konvergenzweg, aber eine Eigenschaft einer **Schleife**, die es hier nicht gibt. Deshalb gemessen und nicht implementiert. |

Die Tabelle liegt als `es::pic_options()` in der Bibliothek und wird in
`tests/test_space_charge.cpp`, Abschnitt 9, geprüft — unter anderem darauf,
dass **genau eine** Variante implementiert ist.

### Was der Abzug leistet, gemessen

Ein **einzelnes** Teilchen im sonst leeren, geerdeten Kasten spürt ohne Abzug
ein Feld zwischen $1{,}13\cdot10^5$ und $1{,}42\cdot10^5$ V/m, obwohl nichts
anderes darin ist — der Betrag hängt allein davon ab, wo in der Zelle es sitzt.
Mit Abzug spürt es **exakt null** (nicht „klein“: der zurückgegebene Vektor ist
die Differenz zweier Lösungen, die identisch sind). Die Linearität wird dabei
gegen die volle Lösung geprüft; der Überlagerungsfehler beträgt
$1{,}6\cdot10^{-15}$.

**Der Preis wird genannt:** eine zusätzliche Lösung je Teilchen. Für eine
Diagnose ist das richtig, für eine Produktionsschleife untragbar.

## 17.4c Die PIC-Schleife bleibt blockiert — aus zwei unabhängigen Gründen

`es::pic_loop_status()` führt beide, damit das Wegfallen eines einzelnen sie
nicht stillschweigend freigibt:

1. **Keine Quelle.** P5 ist blockiert: es gibt keine belegte Emissionsrate und
   damit keine physikalische Teilchenquelle.
2. **Keine tragfähige Kostenstruktur.** Die einzige Selbstfeldbehandlung, die
   dieser Punkt rechtfertigen kann, kostet eine Lösung je Teilchen.

P6 bleibt damit das, was es ist: ein **validierter Poisson- und
Depositions-Teiltest** mit einer exakt begründeten Selbstfeldbehandlung für die
Diagnose — und keine PIC-Schleife.

## 17.5 Was P6 ausdrücklich nicht enthält

* **Keine Emission** und keine selbstkonsistente Emissions-PIC-Schleife. P5 ist
  blockiert; jede hier verwendete Ladungsverteilung ist **vorgeschrieben** und
  heißt in jeder Ausgabe Testquelle.
* **Keine Teilchenbewegung** — das ist P7.
* Keine Magnetfelder, keine Stöße, keine raumladungsbegrenzte Emission.
* **Keine netzunabhängige Regularisierung des Selbstfeldes.** Der exakte Abzug
  (17.4b) entfernt die scheinbare Selbstkraft vollständig, macht aber das
  Selbstpotential nicht netzunabhängig — es gibt es danach am Teilchen schlicht
  nicht mehr. Für alles andere bleibt das deponierte Feld eine Netzgröße.
* Keine skalierte Makropartikelzahl und damit keine PIC-Konvergenzstrategie —
  gemessen, nicht implementiert.
