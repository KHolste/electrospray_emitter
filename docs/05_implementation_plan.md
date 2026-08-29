# 5. Stufenplan

Jede Phase endet an einem **Gate**. Ein Gate ist bestanden, wenn die genannten
Nachweise vorliegen — nicht, wenn der Code läuft. Ohne bestandenes Gate wird die
nächste Phase nicht begonnen. Die Nachweise verweisen auf
[06_validation_matrix.md](06_validation_matrix.md).

Aufwandsangaben sind grobe Größenordnungen in Personentagen und dienen der
Priorisierung, nicht der Terminplanung.

---

## P0 — Absichern, bevor irgendetwas gebaut wird

**Aufwand: 3–5 d**

1. Alle elf Befunde aus [01_gap_analysis.md](01_gap_analysis.md) als
   **fehlschlagende Tests** festschreiben, bevor sie repariert werden. Diese
   Tests bleiben dauerhaft im Bestand.
2. **Erledigt.** Der Absturz mit `-march=native` ist aufgeklärt: AVX-Codegen
   emittiert alignment-pflichtige 256-Bit-Stackzugriffe in einem Frame, den GCC
   nicht auf 32 Byte nachrichtet; die Win64-ABI garantiert nur 16 Byte. Nachweis
   über Minimalfall, Faktorexperiment und Disassemblat in
   [01_gap_analysis.md](01_gap_analysis.md) Abschnitt 1.2. ASan und UBSan unter
   WSL2/Ubuntu melden für alle sechs Tests **und** alle vier Anwendungen keinen
   Befund. Verbleibend: ein Reproduktionsfall ohne Projektabhängigkeit für einen
   GCC-Bugreport (niedrige Priorität, blockiert nichts).
2b. Den Sanitizer-Lauf unter WSL2 als wiederholbaren Schritt einrichten
   (Skript oder CI-Job), damit er bei jeder Phase erneut läuft und nicht einmalig
   bleibt.
3. Ausgabepräfixe je Anwendung eindeutig machen; jede Ausgabedatei bekommt einen
   Kopf mit Programmversion, Commit-Hash, vollständigem Parametersatz und dem
   Zustand, zu dem sie gehört.
4. Haupt-README auf den in [README.md](README.md) festgelegten Sprachgebrauch
   umstellen. Insbesondere die Verifikationstabelle: was dort steht, gilt für
   die Elektrostatik und die Bahnintegration, nicht für die Emitterphysik.

**Gate P0.** Sanitizer-Lauf ohne Befund **(erreicht)**; alle elf Befunde durch
fehlschlagende Tests abgedeckt **(offen)**; Ausgabedateien selbstbeschreibend
und je Anwendung eindeutig benannt **(offen)**; Haupt-README auf die
Sprachregelung umgestellt **(erreicht)**.

---

## P1 — Geometrie und automatische Vernetzung

**Aufwand: 8–12 d. Die Geometriefestlegungen liegen vor (siehe [README.md](README.md)).**

1. Parametrisches Modell nach [04_geometry_model.md](04_geometry_model.md),
   einschließlich der Zwangsprüfungen aus 4.3.
2. Randvernetzung mit krümmungs- und merkmalsgesteuerter Größenfunktion,
   beschränkte Wachstumsrate.
3. Volumenvernetzung des Gebiets (für P4), zunächst strukturiert; unstrukturiert
   nur, falls die Geometrie es erzwingt.
4. A-posteriori-Verfeinerung mit den Indikatoren aus 4.4, Abbruch über
   Zielgrößen.
5. **Bindend:** kein Netz-, Elementgrößen- oder Verfeinerungsparameter in der
   Benutzerschnittstelle. Der Benutzer gibt Geometrie und eine Genauigkeits-
   vorgabe vor, sonst nichts. Die Schlüssel `h_tip` und `h_far` des Prototyps
   entfallen ersatzlos.

**Gate P1.** Netzunabhängigkeit von Kapazität und Spitzenfeld über drei
Verfeinerungsstufen; die Verifikationsgeometrien Kugel und Rotationsellipsoid
werden vom neuen Vernetzer automatisch mit derselben Genauigkeit aufgelöst wie
bisher mit manueller Größenvorgabe; **die Referenzgeometrie der Skizze läuft
ohne jede Netzangabe durch**.

---

## P2 — Elektrostatik konsolidieren

**Aufwand: 4–6 d**

1. `BemSolver` auf beliebig viele Elektroden erweitern.
2. Kapazitätsmatrix und Feldüberhöhung als Standardausgabe.
3. Das Raumladungskriterium aus Spezifikation 1.1 bei jedem Lauf mit ausgeben.
4. Verifikationssatz erweitern: koaxiale Zylinder, Plattenkondensator mit Rand,
   Kugel im homogenen Feld.

**Gate P2.** Alle Fälle der Kategorie A der Validierungsmatrix bestanden, mit
angegebenen Toleranzen.

---

## P3 — Statischer Meniskus, korrekt

**Aufwand: 12–18 d**

1. Newton-Verfahren auf dem gekoppelten Residuum (Form + Randintegral), mit
   Formableitung des Operators. Ersetzt die Fixpunktiteration.
2. Kontaktlinie: beide Fälle (gepinnt mit **geprüfter** Gibbs-Bedingung;
   Kontaktwinkel).
3. `solve_at_voltage` neu: Rückgabe muss die angeforderte Spannung innerhalb
   einer angegebenen Toleranz treffen, sonst `converged = false` mit Angabe der
   erreichten Spannung. Befund 1 darf nicht wieder auftreten können.
4. `find_onset` ersetzen durch eine Fortsetzung mit **Bogenlängenparametrisierung**
   und Stabilitätsauswertung über die Eigenwerte der Jacobi-Matrix. Ein
   Umkehrpunkt wird nur gemeldet, wenn er als innerer Extremwert mit
   Vorzeichenwechsel eines Eigenwerts nachgewiesen ist.
5. Getrennte, klar benannte Ausgabegrößen für: Stabilitätsgrenze,
   Strom-Nachweisgrenze (Eingabeparameter), Cone-Jet-Übergang (in P6).
6. Kriterium von Gallud & Lozano (2022) als unabhängige Zweitbewertung
   implementieren und mit der Eigenwertanalyse vergleichen.

**Gate P3.** Kugelkappe im feldfreien Grenzfall; Wiedergabe der
Stabilitätsgrenzen von Wohlhuter & Basaran (1992) für den dielektrischen
Tropfen; Netzkonvergenz für Apexfeld **und** Apexkrümmungsradius; Befunde 1–3
durch Tests abgedeckt und behoben.

---

## P4 — Volumenlöser und PIC

**Aufwand: 15–20 d**

1. FEM- oder FVM-Poisson-Löser auf dem axialsymmetrischen Gebiet.
2. **Querprüfung gegen die BEM**: bei $\rho_f=0$ müssen beide Verfahren
   dasselbe Feld liefern. Das ist der wertvollste verfügbare Test und der Grund,
   die BEM zu behalten.
3. PIC mit Formfunktionen (cloud-in-cell), Ladungsdeposition, Feldinterpolation.
   Ersetzt die Ringsummation vollständig.
4. Getrennte Spezies mit eigenem $q/m$ und eigener Startverteilung.
5. Polarität durchgängig: Vorzeichen von $E_n$, Speziesauswahl, Massen.

**Gate P4.** Child-Langmuir-Gesetz für raumladungsbegrenzte Emission zwischen
Planparallelplatten auf besser als 2 %; Übereinstimmung mit der BEM bei
$\rho_f=0$; Energieerhaltung wie bisher; Konvergenz gegen die
Teilchenzahl je Zelle gezeigt.

---

## P5 — Emittierender Betrieb, PIR

**Aufwand: 20–30 d. Die inhaltlich schwierigste Phase.**

1. Leaky-Dielectric-Modell nach Spezifikation 3.1/3.2: Strömung,
   Ladungstransport, alle Sprungbedingungen einschließlich der
   **tangentialen** Maxwell-Spannung.
2. Ionenverdampfung mit Speziessumme und polaritätsabhängigen Parametern.
3. Hydraulische Speisung, Kopplung $\Delta p \leftrightarrow Q$.
4. Gekoppeltes Newton-System über alle Unbekannten auf $\Sigma$.
5. $\Delta G$ und der Verdampfungsvorfaktor als deklarierte Anpassungsparameter
   mit ausgegebener Empfindlichkeit $\partial\ln I/\partial\Delta G$.

**Gate P5.** Reproduktion der Meniskusfamilien von Coffman et al. (2016) und der
Emissionskurven von Gallud & Lozano (2022) innerhalb der dort angegebenen
Streuung; Wiedergabe der von Higuera (2008) berichteten Skalierungen
(viskositätsdominierte Strömung, leitfähigkeitskontrollierter Strom);
Netzkonvergenz für Strom und emittierende Fläche.

---

## P6 — Cone-Jet als getrenntes, gekennzeichnetes Modul

**Aufwand: 5–8 d**

1. Korrelationen in ein eigenes Modul, jede mit Gültigkeitsbereich und
   `empirical = true` in der Ausgabe.
2. Kopplung an den Strahltransport über den **Cone-Jet-Strom und die
   Tropfengrößenverteilung**, nicht über die Ionenverdampfungsrate (Befund 6).
3. Warnung bei Anwendung außerhalb des Etablierungsbereichs
   ($\varepsilon_r \lesssim 40$).

**Gate P6.** Die Skalierungen $I\propto\sqrt{Q}$ und $d\propto Q^{1/3}$ gegen
publizierte Messdaten geprüft, nicht gegen die eigene Formel. Getrennte
Ausgabeblöcke nachgewiesen.

---

## P7 — Gekoppelter Betriebspunkt

**Aufwand: 10–15 d**

Iteration Meniskus $\leftrightarrow$ Emission $\leftrightarrow$ Raumladung bis
zur Selbstkonsistenz; Ausgabe von Strom, Schub, $I_\mathrm{sp}$, Divergenz,
Interzeption, Wirkungsgraden, jeweils mit Angabe, welche Modellstufe den Wert
liefert.

**Gate P7.** Konvergenz der Kopplungsschleife nachgewiesen; die von Higuera
(2008) angegebene Vernachlässigbarkeit der Raumladung im PIR wird durch die
eigene Rechnung bestätigt oder begründet widerlegt.

---

## P8 — Validierung gegen Experiment

**Aufwand: laufend**

Vergleich mit publizierten I–U-Kennlinien, TOF-Spektren und Strahlprofilen; wenn
vorhanden, mit eigenen Messdaten (zunächst nicht verfügbar; Datenimport ist
strukturell vorzusehen). Erst nach dieser Phase ist
die Bezeichnung „validiert" für die Emitterphysik zulässig.

---

## Reihenfolge und Abhängigkeiten

```
P0 ──> P1 ──> P2 ──> P3 ──────────────> P5 ──> P7 ──> P8
              │                          ^      ^
              └──> P4 ───────────────────┘      │
                                    P6 ─────────┘
```

P4 kann parallel zu P3 laufen. P6 ist unabhängig und jederzeit einschiebbar,
weil es keine Kopplung an die Feldrechnung hat — was zugleich seine
Beschränkung ist.

## Was in diesem Plan bewusst nicht vorkommt

* **Zeitabhängige Dynamik** (Startvorgang, Tröpfchenabriss, Instabilitäts-
  entwicklung). Offene Frage 8. Wäre eine eigene Ausbaustufe mit
  Grenzflächenverfolgung.
* **Emitterarrays.** Die Axialsymmetrie schließt sie aus; nötig wären 3D-BEM
  oder ein Einheitszellenmodell.
* **Elektrochemie an der Elektrode** (Zersetzung bei Gleichspannungsbetrieb).
* **Wärmehaushalt.** Verdampfungskühlung an der Spitze und die
  Temperaturabhängigkeit von $K$ und $\mu$ koppeln zurück; im Plan nicht
  enthalten, als möglicher Erweiterungspunkt vermerkt.
