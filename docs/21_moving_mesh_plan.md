# 21. Bewegliches Netz — Ist-Stand und Entwicklungsplan

**Stand: 2026-08-30.
Status: Planungsdokument. In diesem Lauf wurde nichts davon implementiert.**

Dieses Dokument trennt drei Dinge, die im Projekt bisher nebeneinander stehen
und leicht verwechselt werden: ein Netz, dessen Knoten sich bewegen; eine
Bewegung, die ein statisches Gleichgewicht sucht; und eine physikalische
Zeitentwicklung. Nur das erste und das zweite gibt es. Es legt außerdem fest,
was erfüllt sein muss, bevor dieses Projekt das Wort *Taylor-Kegel* in einer
Ergebnisaussage benutzen darf.

Ergänzt [10_electrocapillary_model.md](10_electrocapillary_model.md)
Abschnitt 10.2 (das P3b-Netz), [15_free_surface_dynamics.md](15_free_surface_dynamics.md)
(der Vertrag für die Dynamik) und [07_mesher_decision.md](07_mesher_decision.md)
(Vernetzungsentscheidungen).

---

## 21.1 Was das P3b-Netz ist

**Feste Topologie, bewegliche Knotenkoordinaten.**

Das Netz von P3b ist das P2c-Tensorgitter mit einem zusätzlichen axialen Warp
auf den Spalten innerhalb der Bohrung (`build_meniscus_mesh` in
`include/es/electrocapillary.hpp`). Bewegt werden ausschließlich
**Koordinaten**:

* Zellzahl, Zellzuordnung und Konnektivität sind die von P2c und bleiben über
  jeden Lauf hinweg unverändert;
* die Zeile $j_\mathrm{tip}$ wird von $z=0$ auf $z_\mathrm{Meniskus}(r_i)$
  gezogen, dazwischen ist die Abbildung stückweise linear in $z$;
* die Bandränder sind **vorhandene** Gitterzeilen bei $\pm 1{,}5\,a$, keine neu
  eingefügten;
* außerhalb der Bohrung bewegt sich kein Knoten.

Es gibt daher **keine** Knoteneinfügung, keine Knotenentfernung, keinen
Kantentausch und keine Änderung der Elementanzahl. `MeniscusMeshQuality` misst,
was diese Konstruktion an Qualität übriglässt — kleinste Jacobi-Determinante,
Zellstreckung, Scherung, Zahl invertierter Zellen —, aber sie kann nichts
reparieren: `build_meniscus_mesh` wirft, wenn die Form nicht ins Band passt oder
eine Zelle invertieren würde. Es gibt keinen Pfad, auf dem sich das Netz an eine
Form anpasst, die es nicht tragen kann.

## 21.2 Wozu sich dieses Netz bewegt

**Zur Suche eines statischen Gleichgewichts, nicht zu einer Zeitentwicklung.**

`solve_coupled` iteriert einen Fixpunkt aus Form und Last: Feld lösen, Maxwell-
Last projizieren, Young-Laplace-Gleichgewichtsform dazu bestimmen, Netz auf
diese Form ziehen, wiederholen — unterrelaxiert mit `kRelaxation = 0.5`, bis
`kTolShape`, `kTolLoad` und `kTolMechanical` erfüllt sind.

Die Iterationsschritte sind **keine Zeitschritte**. Es gibt keine Zeitvariable,
keine Trägheit, keine Viskosität in der Oberflächenbewegung und keine
kinematische Randbedingung; die Zwischenformen haben keine physikalische
Bedeutung, sondern sind Iterierte eines Nullstellenproblems. Ein
Unterrelaxationsfaktor ist ein Konvergenzmittel und keine Dämpfung eines
physikalischen Vorgangs.

`continue_over_voltage` setzt darauf eine **Fortsetzung über die Spannung**: der
Zweig, der die feldfreie P3a-Lösung enthält, wird durch Erhöhen von $|V|$
verfolgt, ein gescheiterter Schritt wird halbiert, und die Fortsetzung endet,
wenn der Schritt unter `kMinStep` fiele. Auch das ist eine Folge **statischer**
Gleichgewichte, keine Bahn in der Zeit.

## 21.3 Was P4 beisteuert — und was nicht

**Nur die Kinematik für ein vorgeschriebenes Geschwindigkeitsfeld.**

`include/es/surface_kinematics.hpp` implementiert und validiert genau ein Stück:
die kinematische Randbedingung $\dot{\mathbf{x}}\cdot\mathbf{n} =
\mathbf{u}\cdot\mathbf{n}$, in der Zeit integriert auf einem **vom Aufrufer
vorgegebenen** Geschwindigkeitsfeld. Das ist Kinematik, nicht Dynamik: für
$\mathbf{u}$ wird nichts gelöst, und nirgends in dieser Datei wird eine
Kraftbilanz ausgewertet. `solve_dynamic_meniscus()` wirft
`NotImplementedInThisPhase` und benennt, was fehlt. Der Status von P4 ist
`infrastructure_only`.

P4 und P3b sind außerdem **nicht gekoppelt**: P4 bewegt einen Polygonzug, nicht
das P3b-Volumennetz. Es gibt keinen Code, der beide verbindet.

## 21.4 Warum zugespitzte Formen das jetzige Netz sprengen

Ein Warp mit fester Topologie verteilt die vorhandenen Knoten um. Je stärker die
Form von der ebenen abweicht, desto ungleichmäßiger. Für einen schwach gewölbten
Meniskus reicht das; für eine stark zugespitzte Form nicht, aus vier
unabhängigen Gründen:

1. **Auflösung folgt der Topologie, nicht der Krümmung.** Die Knotenzahl entlang
   der Oberfläche ist durch das Tensorgitter festgelegt. Ein Apex mit kleinem
   Krümmungsradius bekommt dieselbe Knotenzahl wie ein flacher — die
   Auflösung *am Apex* wird mit wachsender Zuspitzung schlechter, nicht besser.
2. **Zellqualität verfällt.** Streckung und Scherung im Warpband wachsen mit der
   Apexhöhe. `MeniscusMeshQuality` misst das; das Band ist mit `kBandFactor
   = 1.5` bemessen, und eine Form, die nicht hineinpasst, wird abgewiesen statt
   vernetzt.
3. **Die Kantensingularität setzt die Ordnung.** An der ungerundeten Kontaktkante
   ist die integrierte Maxwell-Kraft der ebenen Form mit der Ordnung 0,541
   konvergent ($1+\beta = 0{,}554$) — zweite Ordnung ist dort nicht zu haben.
   Eine gleichmäßige Verfeinerung zahlt an dieser Stelle also schlecht; nötig
   ist gezielte Verfeinerung.
4. **Das 1-%-Ziel ist gemessen verfehlt.** P0 hat es an genau diesem Netz
   gemessen: bei 1400 V zwischen 4,59 % und 5,67 % geschätzter
   Diskretisierungsfehler, alles `DiscretizationNotConverged`; bei 1000 V sind
   $h/a$ und $E_n$ nicht einmal im asymptotischen Bereich
   (`NotInAsymptoticRange`). Das ist der Ausgangspunkt, nicht eine Sorge.

**Folgerung.** Taylor-Cone-nahe Formen brauchen **adaptive Verfeinerung**,
gesteuert von einem Fehlerindikator, und voraussichtlich **Neuvernetzung** —
also eine Änderung der Topologie, nicht nur der Koordinaten. Ein reines
Warpverfahren mit fester Konnektivität wird dafür nicht ausreichen.

## 21.5 Was vor einer Taylor-Cone-Aussage konvergiert sein muss

Vorab festgelegt, damit die Schranke nicht nachträglich zur Messung passend
gewählt wird. Eine Aussage über einen Taylor-Kegel setzt voraus, dass **alle
vier** folgenden Größen unter Verfeinerung *und* unter Neuvernetzung
konvergieren — jede mit gemessener Ordnung und geschätztem Fehler, wie in
[11_p0_p3b_audit.md](11_p0_p3b_audit.md):

| | Größe | warum sie einzeln nötig ist |
|---|---|---|
| (a) | **Apexkrümmungsradius** | die eigentliche Zielgröße der Zuspitzung; sie ist genau die Größe, deren Auflösung mit der Zuspitzung schlechter wird |
| (b) | **Apexhöhe** | die Formgröße, die die Fortsetzung verfolgt; ohne sie ist (a) nicht einzuordnen |
| (c) | **integrierte Maxwell-Kraft** | die Last, die die Form erzeugt; sie ist die einzige Feldgröße, die die Kantensingularität überlebt (punktweise Werte an der Kante sind nicht netzkonvergent) |
| (d) | **Netzqualität** | eine konvergierte Zahl auf einem entartenden Netz ist kein Ergebnis; kleinste Jacobi-Determinante, Zellstreckung, Scherung und invertierte Zellen müssen über die Stufen beschränkt bleiben |

Konvergenz **unter Neuvernetzung** ist dabei nicht dasselbe wie Konvergenz unter
Verfeinerung derselben Topologie: erst wenn zwei unabhängig erzeugte Netzfolgen
denselben Grenzwert liefern, ist ausgeschlossen, dass der Grenzwert eine
Eigenschaft des Warps ist.

## 21.6 Was kein Taylor-Cone-Onset ist

**Ein Abbruch des Fortsetzungs- oder Netzverfahrens ist kein physikalisches
Ereignis.**

`continue_over_voltage` endet, wenn der Spannungsschritt unter `kMinStep`
fiele; `build_meniscus_mesh` wirft, wenn eine Form nicht ins Warpband passt oder
eine Zelle invertieren würde; `solve_coupled` meldet Nichtkonvergenz nach
`kMaxIterations`. Jeder dieser drei Abbrüche ist eine Eigenschaft **dieses
Lösers und dieses Netzes**. Der Header sagt es bereits im Code:

> The end of the branch is where THIS SOLVER stopped. It is not an emission
> onset, not a Taylor-cone onset and not a stability limit; none of those is
> computed anywhere in this phase.

Ein Onset zu behaupten verlangt eine **Stabilitätsaussage** — einen
Eigenwertwechsel oder einen gerechneten Umkehrpunkt des Zweigs mit
netzkonvergenten Zielgrößen auf beiden Seiten. Nichts davon wird derzeit
gerechnet. Zusätzlich ist der Zweig selbst nur lokal eindeutig: die Nachprüfung
zur Zweig-Mehrdeutigkeit hat festgehalten, dass ein Schnittpunkt im verfolgten
Bereich keine Eindeutigkeit beweist.

## 21.7 Reihenfolge, wenn es angegangen wird

Ohne Terminzusage; jede Stufe endet mit einer Messung, nicht mit einer
Behauptung.

| Stufe | Inhalt | Abnahme |
|---|---|---|
| M1 | Fehlerindikator auf dem vorhandenen Netz, ohne jede Netzänderung | Indikator korreliert über mindestens zwei Formen mit dem gemessenen Diskretisierungsfehler |
| M2 | Anisotrope Verfeinerung entlang der Oberfläche bei fester Topologie im Volumen | (a)–(d) aus 21.5 unter dieser Verfeinerung konvergent |
| M3 | Neuvernetzung mit Übertragung des Zustands | zwei unabhängig erzeugte Netzfolgen liefern denselben Grenzwert; Übertragungsfehler beziffert |
| M4 | Fortsetzung mit Neuvernetzung je Schritt | Zweig unabhängig von der Netzfolge; Abbruchgrund benannt und als Löser- oder Netzabbruch klassifiziert |
| M5 | Stabilitätsaussage | erst hier darf das Wort Onset fallen |

**In diesem Lauf wurde keine adaptive Neuvernetzung implementiert.** Dieses
Dokument beschreibt ausschließlich den Ist-Stand und die Bedingungen; der Code
ist unverändert.
