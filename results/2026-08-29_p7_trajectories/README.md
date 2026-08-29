# P7 — Ionenbahnen und messbare Größen — 2026-08-29

Status: **`validated_subset`**.

```sh
./build/es_trajectories results/2026-08-29_p7_trajectories meta.commit=$(git rev-parse HEAD)
python python/plot_trajectories.py results/2026-08-29_p7_trajectories
```

Modellvertrag und Befunde:
[`docs/18_particle_transport.md`](../../docs/18_particle_transport.md).

## Das Ergebnis ist eine TRANSPORTANTWORT

P5 ist blockiert; es gibt keine physikalische Teilchenquelle. Die
Startverteilung ist vorgeschrieben und der gestartete Strom ist eine EINGABE.
Keine Zahl dieses Laufs ist eine Stromvorhersage.

## Geprueft

| Pruefung | Ergebnis |
|---|---|
| gleichfoermiges Feld gegen die geschlossene Form | Rundungsboden (1e-15 … 1e-14) |
| radiale Koordinate im gleichfoermigen Feld | bewegt sich nicht (< 1e-12 R) |
| Energiegewinn = q dV | 3,1e-13 relativ |
| Energieerhaltung entlang der Bahn | 3,1e-13 |
| Polaritaetsumkehr BEIDER Vorzeichen | Bahn identisch, Flugzeit auf 1e-12 gleich |
| nur die Spezies umgekehrt | laeuft auf den Emitter zurueck |
| Treffererkennung | alle vier Flaechen einzeln getroffen |
| Strombilanz | Schliessfehler **exakt null** |
| Zeitordnung | **2,13**, gemessen im NICHT gleichfoermigen Feld |
| Netzunabhaengigkeit im exakt dargestellten Feld | < 1e-6 |

## Drei Fehler, die der Test gefunden hat

**1. Falsches Vorzeichen im Energieinvariant.** Die Erhaltungsgroesse ist
E = ½mv² + q phi; die erste Fassung schrieb ein Minus. Der Energietest meldete
daraufhin einen relativen Fehler von genau 2 -- was eine invertierte potentielle
Energie erzeugt.

**2. Potential ausserhalb des Netzes ausgewertet.** Der Austrittspunkt lag hinter
dem Rand, `locate()` lieferte null, und die Energiebilanz war um die volle
Beschleunigungsspannung falsch.

**3. Inkonsistenter letzter Schritt.** Eine erste Fassung schnitt die Position
auf den Rand zurueck, behielt aber die Geschwindigkeit des vollen Schritts; die
Austrittsgeschwindigkeit war dann um bis zu 2,5e-3 falsch und die gemessene
Zeitordnung kam negativ heraus. Jetzt wird auf den Schrittanteil bisektiert und
ein konsistenter Verlet-Teilschritt ausgefuehrt.

**Ausserdem gemessen:** im gleichfoermigen Feld ist velocity Verlet EXAKT, weil
die Beschleunigung konstant ist. Eine dort abgelesene Zeitordnung misst nur die
Rundung; die Ordnung wird deshalb in einem nicht gleichfoermigen Feld durch
Selbstkonvergenz gemessen.

## Die Bilanz

extrahiert + abgefangen + noch fliegend = gestartet, **exakt**. Jedes Teilchen
hat genau ein Schicksal, und "noch fliegend" wird getrennt gezaehlt statt als
verloren verbucht.

| Fall | Transmission |
|---|---|
| Kation ohne Raumladung | 70,7 % |
| Kation mit vorgeschriebener Raumladung (3e-14 C) | 31,7 % |
| Anion ohne Raumladung (umgekehrte Spannung) | 70,7 % |

Die vorgeschriebene positive Raumladung defokussiert den positiven Strahl. Das
ist eine Transportantwort auf eine gesetzte Ladung, kein Betriebspunkt.

Dass Kation und Anion dieselbe Transmission haben, ist KEINE Aussage ueber
Symmetrie der Physik: hier wurden nur die Vorzeichen von Ladung und Spannung
zugleich umgekehrt, waehrend Masse, Barriere und Speziesmischung in Wirklichkeit
verschieden sind.

## NICHT gerechnet

Keine Emission, kein Tropfenstrahl, keine selbstkonsistente Raumladungsschleife,
keine Magnetfelder, keine Stoesse, keine Sekundaeremission. Keine
Geraetegeometrie: die Testdomaene ist ein Rechteck mit benannten Flaechen.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_trajectories.png` | Bahnen im r-z-Schnitt mit den Auftrefforten, fuer beide Polaritaeten und ohne/mit vorgeschriebener Raumladung |
| `fig2_balance.png` | Energiegewinn gegen q dV, die Strombilanz als Balken, und die Transmission ohne/mit Raumladung |
