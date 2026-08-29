# 7. Vernetzung — offene Entscheidung

**Stand: 2026-08-29. Status: Entscheidung erforderlich, bevor P1 fortgesetzt wird.**

## Was bereits festgelegt ist

[04_geometry_model.md](04_geometry_model.md), Abschnitt 4.4, legt die
**Strategie** fest:

* Größenfunktion aus Krümmung, Merkmalsgröße und Wandabstand, mit beschränkter
  Wachstumsrate,
* a-posteriori-Verfeinerung an Meniskusspitze, Kontaktlinie, Elektrodenkanten,
  Extraktoröffnung und Strahlkern,
* Abbruch über Netzunabhängigkeit von Zielgrößen, nicht über eine Elementzahl,
* keinerlei Netzparameter in der Benutzerschnittstelle.

## Was nicht festgelegt ist

Weder ein **Verfahren** noch eine **Implementierung**. Konkret offen:

1. **Elementtyp und Zellform.** Dreiecke, Vierecke oder ein Blockgitter mit
   lokaler Verfeinerung.
2. **Algorithmus.** Delaunay mit Größenfunktion, Advancing Front,
   Quadtree/Octree mit Schnittzellen, oder ein strukturiertes Mehrblockgitter,
   das dieser Geometrie angepasst wird.
3. **Eigenbau oder Bibliothek.** Der Stufenplan notiert nur „zunächst
   strukturiert; unstrukturiert nur, falls die Geometrie es erzwingt". Das ist
   keine Entscheidung.

Deshalb wurde in diesem Lauf **kein Vernetzer eingebaut**. Ein provisorischer
Mesher wäre schwer wieder zu entfernen und würde die Entscheidung faktisch
vorwegnehmen.

## Die eigentliche Frage

Die Geometrie selbst ist einfach: achsensymmetrisch, stückweise geradlinig, vier
Gebiete, keine krummen Ränder außer der später beweglichen freien Oberfläche.
Ein strukturiertes Gitter käme damit zurecht. Was die Entscheidung treibt, ist
**nicht** die Ausgangsgeometrie, sondern:

| Anforderung | Konsequenz |
|---|---|
| Die freie Oberfläche bewegt sich in jeder Newton-Iteration (P3) | Entweder Neuvernetzung je Iteration oder ein mitbewegtes Gitter |
| Längenskalen von 10 nm (Emissionsstruktur) bis 3 mm (Domäne) | Fünf bis sechs Dekaden lokale Verfeinerung |
| Die Kontaktlinie ist ein Punkt mit Feldsingularität | Starke, gerichtete Verfeinerung genau dort |
| PIC braucht eine brauchbare Zellstatistik im Strahl (P4) | Zellen im Strahlbereich dürfen nicht beliebig verzerrt sein |
| Querprüfung gegen die vorhandene BEM bei ρ_f = 0 | Das Volumenverfahren muss die BEM-Genauigkeit an der Spitze erreichen |

## Zu entscheiden

**Entscheidung 1 — Eigenbau oder externe Bibliothek.**
Das Projekt hat bisher **keine** externen Abhängigkeiten außer CMake und einem
C++-Compiler; das ist ein Wert, der nicht beiläufig aufgegeben werden sollte.
Eine Netzbibliothek wäre die erste schwere Abhängigkeit.

* *Eigenbau*, z. B. Quadtree mit Schnittzellen: keine Abhängigkeit, auf diese
  Geometrieklasse zugeschnitten, dafür mehrere Wochen Arbeit und ein eigener
  Wartungsgegenstand.
* *Bibliothek*: schneller am Ziel, aber Bindung an deren Datenmodell,
  Lizenz, Bauweise und Plattformverfügbarkeit — letzteres ist auf diesem
  MinGW-Windows-Aufbau nicht selbstverständlich.

**Entscheidung 2 — strukturiert oder unstrukturiert.**
Ein achsensymmetrisches Blockgitter mit lokaler Verfeinerung ist einfacher zu
verifizieren und passt zur Rechteckdomäne. Die bewegliche freie Oberfläche und
die Kontaktliniensingularität sprechen für unstrukturiert.

**Entscheidung 3 — Reihenfolge.**
Der Volumenlöser wird erst in P4 gebraucht. Für P2 und P3 genügt weiterhin die
BEM, die nur eine Randvernetzung braucht — und die ist aus der hier
implementierten Randtopologie unmittelbar erzeugbar. Eine mögliche Antwort auf
Entscheidung 1 und 2 lautet daher: **jetzt nur den Randvernetzer bauen** und die
Volumenfrage bis P4 offen lassen, wenn die Anforderungen aus PIC und
Poisson-Löser konkret vorliegen.

## Was in diesem Lauf gebaut wurde

Nur Geometrie und Randtopologie: benannte Gebiete, benannte Randkurven,
benannte nulldimensionale Merkmale, Rotationsmaße mit analytischen Tests. Das
ist genau die Eingabe, die jeder der genannten Vernetzer braucht, und sie
bleibt unabhängig davon gültig, wie die Entscheidung ausfällt.
