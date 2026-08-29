# 3. Modulentscheidungen

Bewertung jedes vorhandenen Moduls gegen die Spezifikation in
[02_model_specification.md](02_model_specification.md). „Validiert" ist im Sinne
der Tabelle in [README.md](README.md) gemeint.

| Modul | Zeilen | Entscheidung | Begründung |
|---|---:|---|---|
| `elliptic.{hpp,cpp}` | ~120 | **erhalten** | AGM-Auswertung von K und E; gegen exakte Werte und unabhängige Quadratur auf 10⁻¹⁴ geprüft. Die Formulierung über das komplementäre Argument $m_c=d^2/S^2$ vermeidet die Auslöschung bei $m\to1$ und ist genau das, was die Nahfeldpanels brauchen. Kein Änderungsbedarf. |
| `linalg.{hpp,cpp}` | ~110 | **erhalten, später ergänzen** | Dichte LU mit Teilpivotisierung und einem Schritt Nachiteration. Für $N\sim10^3$ ausreichend. Ergänzungsbedarf entsteht erst mit dem Newton-Verfahren (Stufe 2/3): dort werden Blockstrukturen und ein iterativer Löser sinnvoll. |
| `bem.{hpp,cpp}` — Kern und Quadratur | ~200 | **erhalten** | Der axialsymmetrische Kern und die geometrisch verfeinerte Panel-Quadratur sind gegen Kugel- und Rotationsellipsoid-Analytik geprüft (2,5·10⁻⁶ bzw. 0,18 %) und das Randbedingungsresiduum liegt bei 1,5·10⁻¹⁵. Das ist der belastbarste Teil des Prototyps. |
| `bem.{hpp,cpp}` — `BemSolver` | ~250 | **erhalten, erweitern** | Basislösung je Elektrode ist gut (Spannungssweeps kosten nichts). Fehlt: mehr als drei Elektroden, korrekte Polaritätsbehandlung, saubere Schnittstelle zum Volumenlöser. Die vorhandene `set_external_potential`-Kopplung ist konzeptionell richtig (induzierte Ladung auf den Elektroden), aber die Quelle der Außenpotentiale wird ersetzt. |
| `geometry.{hpp,cpp}` | ~330 | **ersetzen** | Die vorhandenen Bausteine bilden die Skizzengeometrie nicht ab: es gibt keinen konisch verjüngten Emitter mit ∅₁/∅₂/∅₃. Die Elementgrößen werden von Hand gesetzt, es gibt keinen Fehlerschätzer und keine Adaptivität. Die Verifikationsgeometrien (`make_sphere`, `make_prolate_spheroid`) und die Grading-Hilfsfunktion `graded_parameters` werden übernommen. |
| `fluid.{hpp,cpp}` | ~200 | **umbauen** | Als Datencontainer brauchbar, inhaltlich unzureichend: keine Quellenangaben, keine Unsicherheiten, nur eine Ionenspezies, keine Anionenseite, generische statt stofflicher VFT-Parameter. Umzubauen auf: Wert + Quelle + Unsicherheit + Gültigkeitstemperaturbereich je Größe; Speziesliste statt einer mittleren Solvatationszahl. |
| `emission.{hpp,cpp}` | ~230 | **aufteilen** | Drei Dinge in einer Datei, die getrennt gehören: (a) `ion_current_density` / `schottky_lowering` — Formel korrekt, erhalten, aber um Speziessumme und Polarität erweitern; (b) `cone_jet` — empirische Korrelation, in ein eigenes, als solches gekennzeichnetes Modul; (c) `onset_voltage_taylor`, `onset_field_hemisphere` — nur noch als Vergleichswerte, nicht als Modellaussagen. |
| `meniscus.{hpp,cpp}` | ~440 | **neu schreiben** | Vier unabhängige Gründe: `solve_at_voltage` meldet Konvergenz bei bis zu 74 % Spannungsabweichung (Befund 1); `find_onset` meldet Umkehrpunkte, wo keine sind (Befund 2); es gibt keine Stabilitätsanalyse, nur ein Maximum; das Modell ist auf perfekte Leitfähigkeit und $\mathbf{u}=0$ beschränkt und deckt damit den emittierenden Betrieb nicht ab. Erhalten wird die ODE-Formulierung in $(r,z,\varphi)$ mit Reihenstart an der Achse — sie ist gegen die Kugelkappe auf 4·10⁻⁶ geprüft. |
| `beam.{hpp,cpp}` | ~380 | **ersetzen** | Das Ring-Makroteilchenmodell ist kein wohlgestelltes Modell (Abschnitt 4.3 der Spezifikation): singuläres Eigenfeld, kein Teilchenmaß, kein definierter Abschneideradius. Zusätzlich falsche Tropfengewichtung und fehlende Polarität. Erhalten werden der Verlet-Integrator mit adaptiver Schrittweite (Energieerhaltung auf 7·10⁻⁶ geprüft) und der Segment-Kollisionstest. |
| `config.{hpp,cpp}` | ~190 | **erhalten** | Einheitensuffixe, Sektionen, CLI-Vorrang und die Warnung über ungenutzte Schlüssel sind für Parameterstudien genau richtig. Erweitern um Parametersätze für die neue Geometrie und um die Speziesliste. |
| `io.{hpp,cpp}` | ~260 | **erhalten, erweitern** | Feldgitter, CSV/VTK, Setup-Aufbau. Zwei Korrekturen nötig: eindeutige Ausgabepräfixe je Anwendung (Zusatzbefund in 01), und die Ausgabedateien müssen zu dem Zustand gehören, der im Bericht ausgewiesen wird (Befund 3). |
| `apps/*` | ~450 | **überarbeiten** | Die Aufteilung in vier Programme bleibt sinnvoll. `es_meniscus` und `es_operate` benutzen die neu zu schreibende Meniskusschicht und müssen entsprechend nachgezogen werden. |
| `tests/*` | ~990 | **erhalten, stark erweitern** | Die analytischen Tests bleiben. Zu ergänzen: alle in [06_validation_matrix.md](06_validation_matrix.md) genannten Fälle, insbesondere Netzkonvergenz für Apexfeld, Emissionsstrom und emittierende Fläche (nicht nur für die Onset-Spannung), sowie ein Test, der die in Befund 1–3 gezeigten Fehlverhalten festnagelt. |

## Zusammenfassung

| Entscheidung | Module |
|---|---|
| erhalten | `elliptic`, `linalg`, `bem` (Kern), `config` |
| erhalten und erweitern | `bem` (Solver), `io`, `apps`, `tests` |
| umbauen / aufteilen | `fluid`, `emission` |
| neu schreiben | `meniscus` |
| ersetzen | `geometry`, `beam` |

Rund 40 % des Codes ist tragfähig, im Wesentlichen die Elektrostatik und die
Infrastruktur. Die Emitterphysik im engeren Sinne — Meniskus im Betrieb,
Emission, Strahl — wird neu aufgebaut.

## Was neu hinzukommt und im Prototyp keine Entsprechung hat

| Neues Modul | Zweck |
|---|---|
| `mesh/` | automatische, adaptiv verfeinerte Vernetzung (Rand und Volumen) |
| `poisson/` | FEM- oder FVM-Löser für die Poisson-Gleichung mit Raumladung |
| `pic/` | Teilchen-in-Zelle mit Formfunktionen, ersetzt die Ringsummation |
| `ehd/` | Leaky-Dielectric-Modell: Strömung, Ladungstransport, Sprungbedingungen |
| `hydraulics/` | Speiseimpedanz, Hagen-Poiseuille bzw. Darcy |
| `species/` | Ionen- und Clusterspezies, Polarität, q/m-Verteilungen |
| `stability/` | Eigenwertanalyse des Newton-Systems entlang des Lösungsastes |
| `provenance/` | Quellen- und Unsicherheitsangaben, die in jede Ausgabedatei mitwandern |
