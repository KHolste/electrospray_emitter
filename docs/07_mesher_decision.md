# 7. Vernetzung — getroffene und offene Entscheidungen

**Stand: 2026-08-29.
Status: Randvernetzung entschieden und gebaut. Volumenvernetzung offen, Frist
gesetzt.**

## 7.1 Was bereits festgelegt war

[04_geometry_model.md](04_geometry_model.md), Abschnitt 4.4, legt die
**Strategie** fest:

* Größenfunktion aus Krümmung, Merkmalsgröße und Wandabstand, mit beschränkter
  Wachstumsrate,
* a-posteriori-Verfeinerung an Meniskusspitze, Kontaktlinie, Elektrodenkanten,
  Extraktoröffnung und Strahlkern,
* Abbruch über Netzunabhängigkeit von Zielgrößen, nicht über eine Elementzahl,
* keinerlei Netzparameter in der Benutzerschnittstelle.

Offen waren Verfahren und Implementierung: Elementtyp, Algorithmus, Eigenbau
oder Bibliothek, und die Reihenfolge.

---

## 7.2 Entscheidung — Randvernetzer jetzt, im Eigenbau

**Getroffen. Umgesetzt in `include/es/boundary_mesh.hpp` und
`src/boundary_mesh.cpp`.**

Ein *Rand*vernetzer für diese Geometrieklasse ist keine offene
Forschungsfrage: die Ränder sind stückweise geradlinige Polylinien in der
Meridianhalbebene. Das ist einige hundert Zeilen, vollständig testbar gegen
geschlossene Formeln, und es begründet keine Abhängigkeit. Die Argumente, die
gegen einen selbstgebauten *Volumen*vernetzer sprechen — Wochen Arbeit,
eigener Wartungsgegenstand, schwer zu verifizieren — gelten hier nicht.

Er wird jetzt gebraucht, nicht später:

| Braucht ihn | wofür |
|---|---|
| P2 (Elektrostatik konsolidieren) | die BEM braucht **nur** eine Randdiskretisierung |
| P3 (Meniskus) | die bewegliche freie Oberfläche wird als Randkurve neu vernetzt |
| Gate P1b | „die Referenzgeometrie läuft ohne jede Netzangabe durch" |

Die Komponente ist bewusst **solverunabhängig**: sie hängt nur von
`device_geometry.hpp` ab, kennt keine Randbedingung, kein Potential und keine
Ladung, und bietet **keine** Umwandlung in das `es::Mesh` der P0-BEM an. Die
P0-BEM behält ihre eigenen Netze; sie umzustellen ist ein eigener, expliziter
Schritt in P2 und passiert nicht stillschweigend.

### Verfahren

Das Netz ist eine reine Funktion der `DeviceParameters`. Größenfunktion:

$$h(x) \;=\; \min\!\Big(\;\min_s \big[\,h_s + G\,|x - x_s|\,\big],\; h_{\max}\Big)$$

ein Minimum über Kegel auf endlich vielen Verfeinerungsquellen. Ein Minimum
G-Lipschitz-stetiger Funktionen ist G-Lipschitz-stetig, also beschränkt das
Feld selbst die Wachstumsrate; es braucht keine Nachglättung, und benachbarte
Elementgrößen können sich höchstens um $(1+G/2)/(1-G/2)$ unterscheiden.

Quellen und Konstanten (dokumentiert in `namespace es::mesher`):

| Konstante | Wert | Bedeutung |
|---|---|---|
| `kFeatureDivisions` | 32 | an einem benannten Merkmal ist `h = lfs / 32` |
| `kCornerDivisions` | 8 | an jeder übrigen Geometrieecke ist `h = lfs / 8` |
| `kDomainDivisions` | 40 | `h_max = Domänendiagonale / 40` |
| `kGradation` | 0.25 | Lipschitz-Konstante *G* |
| `kMinElementsPerSegment` | 4 | jedes gerade Geometriesegment trägt ≥ 4 Elemente |
| `kNodeSnapRelative` | 1e-12 | Knotenidentifikation, relativ zur Diagonale |
| `kMaxNeighbourRatio` | 1.5 | Prüfgrenze für benachbarte Elementgrößen |

`lfs` ist die **lokale Merkmalsgröße** an einem Polylinienpunkt: das Kleinere
aus (a) der kürzesten dort anliegenden Randkante und (b) dem Abstand zur
nächsten Randkante, die dort *nicht* anliegt. Beides folgt allein aus der
Geometrie. Damit verfeinert sich die Austrittskante von selbst, wenn die
Bohrung kleiner wird, und die Aperturkanten folgen der Elektrodendicke statt
dem Aperturdurchmesser — ohne eine einzige Benutzereingabe. `h_tip` und
`h_far` existieren nicht.

Die Knoten eines Segments liegen bei gleichen Zuwächsen des Dichteintegrals
$\int \mathrm{d}s / h(s)$; die Segmentenden sind immer Knoten, sodass jede
Geometrieecke bitgenau erhalten bleibt.

### Was der Randvernetzer ausdrücklich nicht behauptet

Scharfe Kanten werden verfeinert, und das ist beabsichtigt. Es macht ein dort
später berechnetes maximales elektrisches Feld **nicht** zu einer
netzkonvergenten Größe: das Eckfeld einer unverrundeten Kante divergiert und
folgt der örtlichen Elementgröße. Netzkonvergenz darf nur für Größen
beansprucht werden, die nicht an der Singularität ausgewertet werden
(Kapazität, Gesamtladung, Feld in definiertem Abstand). Endliche Kantenradien
sind ein P3-Parameter (`reserved.edge_radius_*`).

Die Fläche bei z = 0 heißt im Modell wie in den Abbildungen **anfängliche
ebene Flüssigkeitsoberfläche — noch kein berechneter Meniskus**
(`boundary_long_name(BoundaryId::FreeSurfaceReference)`). Sie ist der
geometrische Ausgangszustand, aus dem keine physikalische Meniskuslösung
abgeleitet werden darf.

---

## 7.3 Entscheidung — Volumenvernetzung: Frist vor P3, nicht P4

> **Nachtrag 2026-08-29 (P2b): das Verfahren für die STATISCHE Geometrie ist
> getroffen, die Frist bleibt.** Gebaut ist ein *blockstrukturierter, radial
> gewarpter* Vernetzer, `include/es/volume_mesh.hpp`, der eine reine Funktion
> der Geräteparameter ist. Er verletzt keine der beiden Bindungen unten: er ist
> **kein** allgemeiner unstrukturierter Vernetzer (kein Delaunay, kein Octree,
> keine Qualitätsgarantien für beliebige Eingaben — er kennt genau diese
> Geometrieklasse), und er führt **keine** externe Abhängigkeit ein. Weil die
> einzige nicht achsenparallele Grenze dieses Geräts die gerade Kegelflanke ist,
> genügt ein Tensorproduktgitter, dessen Radien so verzerrt sind, dass eine
> Gitterlinie exakt auf der Flanke liegt; jede Materialgrenze fällt damit exakt
> auf Gitterlinien, und jedes Gebietsvolumen wird exakt reproduziert.
>
> **Die Frist aus diesem Abschnitt ist damit nicht erledigt.** Dieser Vernetzer
> adaptiert nicht, verfeinert nicht a posteriori und führt keine bewegte freie
> Oberfläche. P3 braucht alle drei. Die Entscheidung unten ist weiterhin vor
> Beginn von P3 zu treffen; der hier gebaute Vernetzer ist billig wegzuwerfen.
>
> Begründung, Prüfliste und die beiden beim Bau gemessenen Fehler:
> [08_dielectric_model.md](08_dielectric_model.md), Abschnitt 8.5.

**Getroffen ist die Frist, nicht das Verfahren.**

Die frühere Fassung dieses Dokuments schlug vor, die Volumenfrage bis P4
offenzulassen. Das ist falsch und wird hiermit korrigiert: **bereits das
gekoppelte Flüssigkeitsmodell in P3 braucht eine Volumendiskretisierung** —
die Strömung in Bohrung und Meniskus ist ein Feldproblem im Flüssigkeitsgebiet
und nicht mit einer Randmethode allein zu behandeln. Die Entscheidung ist
deshalb **vor dem Beginn von P3 zu treffen**, nicht vor P4.

Bindend für diese Entscheidung:

* **Kein allgemeiner Volumenvernetzer im Eigenbau.** Der Randvernetzer war
  vertretbar, ein Delaunay- oder Octree-Vernetzer mit Qualitätsgarantien ist
  es nicht.
* **In diesem Lauf keine schwere externe Abhängigkeit.** Das Projekt hat außer
  CMake und einem C++-Compiler keine; eine Netzbibliothek wäre die erste. Sie
  wird bewusst und geprüft eingeführt, nicht beiläufig.

Was die Entscheidung treibt, ist nicht die einfache Ausgangsgeometrie, sondern:

| Anforderung | Konsequenz |
|---|---|
| Die freie Oberfläche bewegt sich in jeder Newton-Iteration (P3) | Neuvernetzung je Iteration oder mitbewegtes Gitter |
| Längenskalen von 10 nm bis 3 mm | fünf bis sechs Dekaden lokale Verfeinerung |
| Die Kontaktlinie ist ein Punkt mit Feldsingularität | starke, gerichtete Verfeinerung genau dort |
| PIC braucht brauchbare Zellstatistik im Strahl (P4) | Zellen im Strahlbereich dürfen nicht beliebig verzerrt sein |
| Querprüfung gegen die BEM bei ρ_f = 0 | das Volumenverfahren muss die BEM-Genauigkeit an der Spitze erreichen |

Der hier gebaute Randvernetzer liefert für jede der Kandidatenlösungen die
Eingabe, die sie braucht — eine konforme, orientierte, gebietstreue
Randdiskretisierung mit erhaltenen Merkmalen — und bleibt gültig, wie die
Volumenentscheidung auch ausfällt.

---

## 7.4 Was in diesem Lauf gebaut wurde

Randvernetzer, Prüfungen, Werkzeug, Abbildungen, Tests:

* `include/es/boundary_mesh.hpp`, `src/boundary_mesh.cpp` — Größenfunktion,
  Netzerzeugung, Topologieprüfungen, Statistik, CSV-Ausgabe.
* `tools/mesh_figure.cpp` — erzeugt Netz und Prüfbericht aus einer
  Konfigurationsdatei; Rückgabewert 2, wenn eine Prüfung fehlschlägt.
* `python/plot_mesh.py` — drei Abbildungen aus den geschriebenen CSVs.
* `tests/test_boundary_mesh.cpp` — Größenfeld gegen geschlossene Formeln,
  achsensymmetrischer Vertrag, Topologie, Gebietsumläufe, Determinismus,
  vier Parametervarianten.

Die Prüfliste, die `BoundaryMesh::validate()` abarbeitet und die im
Ergebnisordner mitgeschrieben wird:

1. keine Null­längen-Elemente, keine Nullradien-Ringelemente,
2. Achsenelemente exakt bei r = 0, Rotationsfläche exakt null,
3. keine doppelten Elemente,
4. kein Knoten mit weniger als zwei Elementen (keine Lücken),
5. jede Randkurve ist eine ununterbrochene Kette mit den richtigen Enden,
6. alle Geometrieecken und alle benannten Merkmale sind exakt Netzknoten,
7. kein Element überspannt eine Geometrieecke oder Materialgrenze,
8. keine ungewollten Überlappungen oder Kreuzungen,
9. jedes Element trennt genau die angegebenen Gebiete (geometrisch nachgeprüft,
   nicht aus der Punktreihenfolge übernommen),
10. Gebietsränder geschlossen und mathematisch positiv orientiert; Meridian­fläche
    und Rotationsvolumen jedes Gebiets werden aus dem Netz reproduziert,
11. Rotationsflächen und Meridianlängen bleiben exakt erhalten,
12. Elementgrößenverhältnis benachbarter Elemente beschränkt,
13. bitweise reproduzierbare Knotenanordnung.

Ergebnisse des P1-Beispielparametersatzes:
[`results/2026-08-29_p1_boundary_mesh/`](../results/2026-08-29_p1_boundary_mesh/).
