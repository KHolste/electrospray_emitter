# Projektdokumentation — Neuausrichtung

**Stand: 2026-08-29. Status des Codes: Prototyp. Keine validierte Gesamtsimulation.**

Diese Dokumente ersetzen die Darstellung im Haupt-README, soweit sie den
Reifegrad betrifft. Das Haupt-README beschreibt, was der Prototyp *tut*; die
Dokumente hier beschreiben, was fehlt, was falsch ist und was gebaut werden soll.

| Dokument | Inhalt |
|---|---|
| [01_gap_analysis.md](01_gap_analysis.md) | Ist-Stand gegen Ziel; Verifikation der gemeldeten Fehler |
| [02_model_specification.md](02_model_specification.md) | Gleichungen, Unbekannte, Randbedingungen, Kopplungen |
| [03_module_decisions.md](03_module_decisions.md) | Erhalten / korrigieren / ersetzen, je Modul begründet |
| [04_geometry_model.md](04_geometry_model.md) | Parametrisches Geometriemodell zur Skizze |
| [05_implementation_plan.md](05_implementation_plan.md) | Stufenplan mit Freigabe-Gates |
| [06_validation_matrix.md](06_validation_matrix.md) | Analytische, numerische, experimentelle Referenzen |
| [07_mesher_decision.md](07_mesher_decision.md) | Vernetzung: Randvernetzer entschieden und gebaut, Volumenfrage vor P3 |
| [references.md](references.md) | Literaturverzeichnis mit Prüfstatus |

---

## Sprachregelung in diesen Dokumenten

Um genau die Vermischung zu vermeiden, die zur Neuausrichtung geführt hat:

| Begriff | Bedeutung hier |
|---|---|
| **implementiert** | Code existiert und läuft. Keine Aussage über Richtigkeit. |
| **gegen X geprüft** | Numerisches Ergebnis wurde mit X verglichen; Abweichung angegeben. |
| **validiert** | Gegen eine *unabhängige* Referenz geprüft, die nicht aus diesem Code stammt. |
| **offen** | Nicht geprüft, oder Prüfung nicht abgeschlossen. |

Ein interner Konsistenztest und die Reproduktion einer einprogrammierten
Skalierung sind ausdrücklich **keine** Validierung.

---

## Festgelegte Rahmenbedingungen

Diese Punkte waren offen und sind entschieden. Sie sind für den gesamten
Entwurf bindend.

### Geometrie (Referenzfall)

| Punkt | Festlegung |
|---|---|
| Bohrung | **zylindrisch**, ∅₂ über die ganze Länge konstant. Die grüne Verjüngung in der Skizze ist der Meniskus, nicht die Bohrung. |
| Äußeres Rechteck | **Simulationsdomäne**, kein leitendes Gehäuse auf festem Potential. |
| Stromab des Extraktors | zunächst **verlängerte offene Domäne**. Ein Kollektor ist als optionale Elektrode vorzusehen, aber nicht im Referenzfall. |
| Austrittskante | zunächst **scharf**, Kontaktlinie **gepinnt** (Gibbs). Ein endlicher Kantenradius mit Kontaktwinkelbedingung ist als spätere Erweiterung strukturell vorzusehen. |
| Emitterkörper | **massiv mit zentraler Kapillare**. Speisung über Hagen-Poiseuille bzw. die daraus abgeleitete hydraulische Impedanz. Poröse Emitter mit Darcy-Strömung sind ein **separates späteres Modell**. |

### Physik und Ausbaustufen

| Punkt | Festlegung |
|---|---|
| Betriebsregime | Langfristig PIR **und** Cone-Jet/Mischbetrieb, aber als **getrennte Modelle**. Priorität: ein wissenschaftlich belastbares **PIR**-Modell. |
| Zeitabhängigkeit | Zunächst **quasistationär**. Zeitabhängigkeit erst für eine spätere Ausbaustufe (Instabilitäten, Einschwingvorgänge, Regimewechsel). |
| Validierungsdaten | Zunächst **keine eigenen Messdaten**. Validierung gegen unabhängige Literaturdaten; ein Datenimport für eigene Messungen ist strukturell vorzusehen. |

### Netzerzeugung — bindende Anforderung

**Die Vernetzung erfolgt automatisch. Der Benutzer gibt ausschließlich
Geometrie- und Genauigkeitsparameter vor, niemals ein Netz, Elementgrößen oder
Verfeinerungszonen.**

Konkret bedeutet das:

* Aus den Geometrieparametern (∅₁, ∅₂, ∅₃, Höhe, Profil, Extraktorabstand,
  Apertur, Dicke) wird das Rand- **und** das Volumennetz ohne weitere Eingaben
  erzeugt.
* Die Elementgröße folgt einer Größenfunktion aus lokaler Krümmung und
  Merkmalsgröße mit beschränkter Wachstumsrate — nicht aus Benutzerangaben.
* Verfeinert wird a posteriori an Meniskusspitze, Kontaktlinie,
  Elektrodenkanten, Extraktoröffnung und Strahlkern, gesteuert von einem
  Fehlerindikator.
* Abbruchkriterium ist die **Netzunabhängigkeit von Zielgrößen** (Apexfeld,
  Apexkrümmungsradius, Emissionsstrom, emittierende Fläche, Divergenzwinkel),
  nicht eine Elementzahl.
* Der Benutzer steuert nur über eine Genauigkeitsvorgabe (z. B. zulässige
  relative Abweichung der Zielgrößen) und optional eine Obergrenze für den
  Rechenaufwand.

Der derzeitige Prototyp erfüllt das **nicht**: dort werden `h_tip` und `h_far`
von Hand gesetzt, es gibt keinen Fehlerschätzer und keine Verfeinerung.
Ausführung in [04_geometry_model.md](04_geometry_model.md), Abschnitt 4.4;
Umsetzung in Phase P1 des [Stufenplans](05_implementation_plan.md).

---

## Verbleibende offene Punkte

Keine, die den Entwurf blockieren. Zwei Punkte sind zum Zeitpunkt der
Implementierung zu klären:

1. Genauigkeitsvorgabe für die adaptive Verfeinerung: welche relative
   Abweichung der Zielgrößen soll die Voreinstellung sein?
2. Für die spätere Kollektorelektrode: Position, Potential und ob sie als
   Absorber oder als abbildendes Element modelliert werden soll.
