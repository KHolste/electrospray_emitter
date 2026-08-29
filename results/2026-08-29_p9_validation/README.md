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

Einheit, Unsicherheit MIT TYP (GUM A oder B), Fundstelle und Geometrieart sind
Pflicht; dazu Ueberdeckungsfaktor und Bedingungen. Jede fehlende Angabe schlaegt
geschlossen fehl, und ein Satz mit einem unvollstaendigen Punkt wird als GANZES
abgelehnt -- einen Punkt wegzulassen waere ein Vergleich mit einem anderen
Datensatz.

**Keiner der Punkte in `import_contract.csv` ist eine Messung.** Sie fuehren den
Vertrag vor.

## Die Validierungsmatrix

13 Groessen: **6 direkt vergleichbar**, **3 erst nach ausgesprochener
Reduktion**, **4 grundsaetzlich nicht** -- azimutale Asymmetrie, Versatz von
Emitter und Blende, Neigung des Emitters, Array-Uebersprechen. Diese vier sind
genau das, was ein achsensymmetrisches Modell nicht hat; sie stehen benannt in
der Matrix statt zu fehlen, und sie sind die Antwort auf die Frage, wofuer ein
3D-Loeser gebraucht wuerde.

Von den 13 rechnet dieses Projekt 8 -- und **keine davon traegt den Status
"validiert"**.

## Kunze-Geometrien und Messdaten

Die Geraetegeometrie modelliert die geraden 10-um-Kapillaren aus SU-8 und IP-Q,
fuer die die Kunze-Dissertation EMI-BF4 ausweist. **Messdaten daraus sind nicht
importiert**: der Vertrag verlangt zu jedem Punkt Einheit, Unsicherheit mit Typ,
Fundstelle und Geometrieart, und diese Zuordnung ist eine Arbeit am Dokument,
die dieser Lauf nicht geleistet hat. Das Schema steht bereit; die Daten fehlen.

## Abbildung

| Datei | Inhalt |
|---|---|
| `fig1_validation.png` | die Validierungsmatrix nach Vergleichbarkeit, die Rotationsreferenz ueber der Azimutzahl, und der Importvertrag an einem absichtlich unvollstaendigen Beispielsatz |
