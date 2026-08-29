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
Skalierung sind ausdrücklich **keine** Validierung. Wo im bisherigen Haupt-README
„verifiziert" steht, ist im Sinne dieser Tabelle „gegen eine analytische Lösung
geprüft" gemeint — das trifft für die Elektrostatik zu und für nichts sonst.

---

## Offene Fragen, die vor Implementierungsbeginn zu klären sind

Diese Punkte blockieren Teile des Entwurfs. Sie sind in
[04_geometry_model.md](04_geometry_model.md) und
[02_model_specification.md](02_model_specification.md) an der jeweiligen Stelle
ausführlicher begründet.

**Geometrie (aus der Skizze nicht eindeutig ablesbar)**

1. Ist die Bohrung zylindrisch (∅₂ über die ganze Länge) oder verjüngt sie sich
   von einem größeren Wert an der Basis auf ∅₂ am Austritt? Die Skizze zeigt die
   beiden inneren Begrenzungslinien annähernd parallel; die grüne Füllung ist
   dagegen konisch gezeichnet. Der Entwurf nimmt vorläufig **zylindrisch** an
   und bietet die Verjüngung als optionalen Parameter an.
2. Ist das äußere Rechteck („mesh") ein physikalisches Gehäuse auf definiertem
   Potential, oder nur die Berechnungsdomäne? Für die BEM ist es entbehrlich,
   für den Poisson-Löser der Raumladung nicht.
3. Existiert stromab der Extraktionselektrode ein Kollektor bzw. eine
   Beschleunigungselektrode auf eigenem Potential?
4. Ist die Stirnfläche am Emitteraustritt scharfkantig oder verrundet? Das
   entscheidet, ob die Kontaktlinie gepinnt ist (Gibbs) oder über einen
   Kontaktwinkel läuft (Young) — zwei verschiedene Randbedingungen mit
   verschiedenen Lösungsmengen.
5. Ist der Emitterkörper ein massives, innen benetztes Röhrchen (wie gezeichnet)
   oder poröses Material? Der Speiseweg unterscheidet sich (Hagen-Poiseuille
   gegen Darcy).

**Physik / Betriebsbereich**

6. Welcher Betriebsbereich ist das Ziel: reines PIR, Cone-Jet, oder beides? Der
   Aufwand unterscheidet sich erheblich; PIR ist mit stationären Rechnungen
   erreichbar, der Cone-Jet-Übergang nicht.
7. Liegen eigene Messdaten vor (I–U-Kennlinien, TOF-Spektren, Strahlprofile)?
   Ohne solche Daten bleibt Validierungsstufe C
   ([06_validation_matrix.md](06_validation_matrix.md)) auf publizierte Daten
   Dritter beschränkt.
8. Soll die zeitabhängige Dynamik (Startvorgang, Instabilität, Tröpfchenabriss)
   erfasst werden, oder genügen stationäre Betriebspunkte?
