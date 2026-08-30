# P9 — Achsensymmetrie, 3D-Erweiterung und Validierung — 2026-08-29

Status: **`infrastructure_only`**. Es gibt **kein 3D-Netz, keinen 3D-Loeser und
kein 3D-Ergebnis**.

```sh
./build/es_validation results/2026-08-29_p9_validation meta.commit=$(git rev-parse HEAD)
python python/plot_validation.py results/2026-08-29_p9_validation
```

Vollstaendig: [`docs/20_validation_3d.md`](../../docs/20_validation_3d.md).

## Die Geometrieart reist mit dem Ergebnis

`LabelledResult` traegt seine `GeometryKind`, sie wird im Konstruktor gesetzt,
und es gibt KEINEN Setter. `value_as_three_dimensional()` ist das einzige
Gatter, durch das eine 3D-Behauptung muss, und es wirft fuer alles ausser einer
echten 3D-Rechnung -- auch fuer eine Rotationsreferenz.

Das ist keine Pedanterie: ein achsensymmetrisches Modell sieht genau das nicht,
worin sich eine reale Anordnung von ihm unterscheidet. Die Fehlerart ist eine
Zahl, die richtig aussieht und eine andere Frage beantwortet.

## Die Rotationsreferenz

Dieselbe Groesse zweimal: achsensymmetrisch als int f 2 pi r ds, und als
EXPLIZITE 3D-Quadratur ueber den Azimut.

| Groesse | Unterschied | geschlossene Form |
|---|---|---|
| Flaeche der Halbkugel | <= 1e-12 | 2 pi R^2, getroffen auf 8e-8 |
| int (z/R) dA | <= 1e-12 | pi R^2, getroffen auf 1,3e-7 |
| Rotationsvolumen | <= 1e-12 | (2/3) pi R^3, getroffen auf 1,5e-7 |

Der Rest ist **Summationsrundung**, nicht ein Modellunterschied: die
3D-Quadratur addiert 720 x 2000 Terme, wo die andere 2000 addiert. Die
Azimutzahl aendert nichts -- 3 und 997 geben dasselbe, weil die Mittelpunktregel
fuer einen azimutunabhaengigen Integranden exakt ist.

**Damit ist die 2 pi r-Wichtung, die jedes Integral dieses Projekts traegt, als
das dreidimensionale Integral NACHGEWIESEN.** Das ist kein 3D-Loeser: es wird
nichts Neues geloest, es wird eine Wichtung geprueft.

## Der Importvertrag

**Harte Bedingungen** -- eine fehlende schlaegt geschlossen fehl und lehnt den
Satz als GANZES ab, denn einen Punkt wegzulassen waere ein Vergleich mit einem
anderen Datensatz:

| Feld | fehlt → |
|---|---|
| `unit` | `MissingUnit` |
| `provenance` | `MissingProvenance` |
| `geometry_stated` | `MissingGeometryKind` |
| zwei Einheiten fuer dieselbe Groesse | `UnitMismatch` |

**Eine nicht angegebene Unsicherheit ist KEINE davon.** Eine fruehere Fassung
behandelte sie wie eine fehlende Einheit und warf damit echte Messungen weg. Die
beiden Faelle sind nicht dasselbe: eine fehlende Einheit macht die Zahl
unlesbar; eine fehlende Unsicherheit macht sie unbrauchbar fuer einen
*quantitativen* Vergleich -- die Messung selbst existiert aber, ist zitierbar und
traegt eine Aussage.

Dafuer gibt es den ausdruecklichen Zustand `UncertaintyType::NotReported` und
den Importstatus `ImportStatus::OkUncertaintyNotReported`:

* der Punkt wird **importiert und archiviert**;
* er darf **qualitativ dargestellt** werden, mit sichtbarem Status;
* `usable_quantitatively()` ist fuer ihn **falsch** -- jede Abweichung, jedes
  Chi-Quadrat und jedes Bestanden/Durchgefallen muss genau das abfragen.

Ein **gemischter Satz** bleibt ganz: `n_quantitative` und `n_qualitative_only`
werden getrennt gezaehlt, und jeder Punkt traegt seinen eigenen Status.

**Keiner der Punkte in `import_contract.csv` ist eine Messung.** Sie fuehren den
Vertrag vor.

## Die Validierungsmatrix: sechs Fragen, nicht eine

**Was daran falsch war.** Die fruehere Matrix trug *eine* Vergleichbarkeit je
Zeile, und die Abbildung faerbte die Zeile danach. Eine Groesse, die im Prinzip
vergleichbar, aber blockiert oder nicht konvergiert ist, erschien damit **gruen**
-- und gruen liest sich als Erfolg. Der klarste Fall ist der **Gesamtstrom**:
direkt vergleichbar, und von diesem Projekt ueberhaupt nicht rechenbar, weil P5
blockiert ist.

Jede Zeile traegt jetzt **sechs unabhaengige Urteile** (`yes` / `partial` /
`no` / `n/a`):

| Achse | Frage |
|---|---|
| `comparable_geometry` | laesst sich die Groesse zwischen achsensymmetrisch und 3D ueberhaupt vergleichen? |
| `implemented` | rechnet dieses Projekt sie? |
| `converged` | ist das Ergebnis nach einem **vorab** festgelegten Kriterium konvergiert? |
| `comparable_with_data` | liesse sie sich mit einer Messung vergleichen? |
| `validated` | ist sie **tatsaechlich** verglichen worden und hat innerhalb der angegebenen Unsicherheiten uebereingestimmt? |
| `blocked` + Grund | ist sie blockiert, und wodurch? |

In der Abbildung heisst gruen nur, dass **genau diese eine Achse** erfuellt ist.
Die Blockiert-Spalte steht **neben** den anderen, nicht statt ihrer.

**Die Invariante steht im Code**: `validated=yes` verlangt `implemented=yes` UND
`converged` in {yes, n/a} UND `comparable_with_data=yes` UND nicht blockiert;
und blockiert genau dann, wenn ein Grund genannt ist. Der Test baut zusaetzlich
von Hand Zeilen, die sie verletzen, und prueft, dass sie abgefangen werden.

### Was dabei herauskommt

| Achse | von 13 erreicht |
|---|---|
| vergleichbar (ganz oder nach Reduktion) | **9** |
| implementiert | **8** |
| numerisch konvergiert | **0** |
| mit Messdaten vergleichbar | **6** |
| **tatsaechlich validiert** | **0** |
| blockiert | **5** |

**Der Abstand zwischen der ersten und der vorletzten Zahl ist der ganze Punkt
dieser Tabelle.** Dass nichts validiert ist, ist eine **gepruefte** Aussage: es
sind ueberhaupt keine Messdaten importiert, und der Lauf setzt `exit_code=2`,
sobald eine Zeile etwas anderes behauptet.

Die **vier grundsaetzlich nicht vergleichbaren** Groessen -- azimutale
Asymmetrie, Versatz von Emitter und Blende, Neigung des Emitters,
Array-Uebersprechen -- sind genau das, was ein achsensymmetrisches Modell nicht
hat. Sie stehen benannt in der Matrix statt zu fehlen, und sie sind die Antwort
auf die Frage, wofuer ein 3D-Loeser gebraucht wuerde.

Jede Zeile traegt ausserdem je einen Satz zu Konvergenz und Validierung statt
einer Farbe. **Eine Zeile hat sich durch die P3-Korrektur geaendert:** die
Ladungsrelaxationszeit stand als `MissingMaterialData` und steht jetzt als
„Band belegt, Einzelwert fehlt" -- implementiert, geschlossen loesbar, mit
Messdaten vergleichbar, und weiterhin **nicht validiert**.

## Kunze-Geometrien und Messdaten

Die Geraetegeometrie modelliert die geraden 10-um-Kapillaren aus SU-8 und IP-Q,
fuer die die Kunze-Dissertation EMI-BF4 ausweist. **Messdaten daraus sind nicht
importiert**: der Vertrag verlangt zu jedem Punkt Einheit, Unsicherheit mit Typ,
Fundstelle und Geometrieart, und diese Zuordnung ist eine Arbeit am Dokument,
die dieser Lauf nicht geleistet hat. Das Schema steht bereit; die Daten fehlen.

## Abbildung

| Datei | Inhalt |
|---|---|
| `fig1_validation.png` | die Validierungsmatrix mit einer Spalte je Achse und einer eindeutigen Legende, die Rotationsreferenz ueber der Azimutzahl, und der Importvertrag mit seinen drei Ausgaengen |

## Dateien

`validation_matrix.csv` -- alle sechs Achsen je Groesse, der Blockergrund und je
ein Satz zu Konvergenz und Validierung. `validation_tally.csv` -- wie viele
Zeilen jede Achse erreichen; die vorletzte Zahl ist null.
`import_contract.csv` -- die drei Ausgaenge des Vertrags und ein gemischter Satz.
`revolution.csv` -- die Rotationsreferenz.
