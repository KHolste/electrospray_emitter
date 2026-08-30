# P2 — Stoffdaten für EMI-BF4 — 2026-08-29

```sh
./build/es_material results/2026-08-29_p2_material_data meta.commit=$(git rev-parse HEAD)
python python/plot_material.py results/2026-08-29_p2_material_data
```

Audit, Quellenlage und Auswahlregel: [`docs/13_material_data.md`](../../docs/13_material_data.md).

## Woher die Zahlen kommen

**NIST ILThermo v2.0**, Standard Reference Database #147. Eine Zusammenstellung:
jeder Datensatz traegt die vollstaendige Zitation der PRIMAERpublikation samt
Messmethode und Probenbeschreibung. 154 Datensaetze, 1561 Messpunkte fuer
EMI-BF4. Abgefragt von `tools/fetch_material_data.py`, das die Zahlen ohne
Zwischenschritt in die versionierte Datei `src/material_data_emibf4.cpp`
schreibt.

**IoLiTec Technical Data Sheet IL-0006** (Revision 2/13/2012, ausgegeben
11/2/2012). Das PDF wurde geoeffnet und gelesen; uebernommen sind Dichte
(24 °C), Viskositaet (25 °C), Leitfaehigkeit (25 °C) und die Reinheit >98 %.
Das Blatt nennt KEINE Oberflaechenspannung, KEINE Permittivitaet, KEINEN
Wassergehalt und KEINE Messmethode -- diese Felder bleiben leer.

Keine Zahl stammt aus einem Suchtreffer-Snippet und keine aus dem Gedaechtnis.

## Kernzahlen bei 298,15 K

| Groesse | gewaehlt | Band (±2 K) | Streuung | Quellen |
|---|---|---|---|---|
| Oberflaechenspannung | 0,05401 N/m | 0,0443 … 0,0546 | 19,1 % | 11 |
| Dichte | 1280,9 kg/m³ | 1277,6 … 1290,0 | 1,0 % | 60 |
| dyn. Viskositaet | 0,03637 Pa·s | 0,0252 … 0,0450 | 54,4 % | 22 |
| el. Leitfaehigkeit | 1,5584 S/m | 0,91 … 1,63 | 46,2 % | 21 |
| kin. Viskositaet | **28,3941 mm²/s** (`derived`) | 27,1 … 27,5 (direkt gemessen) | ± 2,09 % (fortgepflanzt) | 2 Eltern |
| rel. Permittivitaet | **MissingMaterialData** | — | — | 4 |

Gewaehlt wird nach einer vorab festgelegten, mechanisch angewandten Regel:
Methode, Reinheit UND Wassergehalt muessen angegeben sein; Umgebungsdruck; nicht
frequenzaufgeloest; danach die meisten Punkte. Erfuellt keine Quelle die Regel,
wird NICHTS gewaehlt. **Es wird nicht gemittelt.**

## Der Befund, der P3a und P3b betrifft

gamma geht LINEAR in das statische Gleichgewicht ein; die Spannung skaliert mit
sqrt(gamma).

* bisher benutzt (ohne Quelle, `illustrative`): **0,0452 N/m**
* gewaehlte Quelle (Souckova, Klomfar, Patek 2011, *Fluid Phase Equilib.*
  303(2), 184-190; 99,8 mass %, 0,0164 mass % Wasser, Wilhelmy-Platte UND
  du-Noüy-Ring, 14 Punkte 288-356 K): **0,05401 N/m**
* Der bisherige Wert liegt **unter jeder Quelle, die Reinheit und Wassergehalt
  angibt** (die niedrigste davon ist 0,0501 N/m).
* Er liegt NICHT unter dem gesamten Band: dessen untere Kante 0,0443 N/m ist
  eine einzelne Kapillaraufstiegsmessung (Martino, de la Mora et al. 2006) ohne
  Reinheits- und Wasserangabe. Sie wird weder geloescht noch gemittelt.

**Faktor auf jede Druckskala 1,195, auf jede Spannung 1,093.** Das ist groesser
als jeder Diskretisierungsfehler, den P0 gemessen hat (4,6 bis 6,1 %).

Der bisherige Datensatz `emibf4_illustrative()` bleibt unveraendert und bleibt
`illustrative`. Der belegte ist ein zweiter, `emibf4_sourced()`.

## Was fehlt

* **relative Permittivitaet**: keine der vier Quellen nennt Reinheit und
  Wassergehalt; die einzige mehrpunktige (Bennett et al. 2019) misst von 1 bis
  18 GHz, wobei eps_r von 11,7 auf 7,7 faellt -- eine Dispersionsmessung und
  kein Gleichstromwert. Fuer P3b ist das folgenlos (die Fluessigkeit ist dort
  ein idealer Leiter), fuer jedes Modell mit endlicher Leitfaehigkeit nicht.
* **kinematische Viskositaet als DIREKTE Messung**: keine der drei Quellen
  nennt Methode, Reinheit und Wassergehalt zugleich, also waehlt die Regel keine
  aus. Ein Wert wird aber ABGELEITET -- siehe unten.

## nu = mu/rho -- abgeleitet, nicht gemessen

nu ist keine unabhaengige Stoffgroesse. Sie gar nicht zu bilden waere zu streng;
sie aus irgendeinem mu und irgendeinem rho zu bilden waere schlimmer -- der
Quotient waere die Viskositaet einer Probe geteilt durch die Dichte einer
ANDEREN. Die Ableitung ist deshalb genau dann erlaubt, wenn beide Elternwerte
nachweislich dieselbe Fluessigkeit im selben Zustand beschreiben. Die vier
Bedingungen stehen in `include/es/material_data.hpp` und werden einzeln
geprueft; ist eine verletzt, wird kein Wert geliefert, sondern die Bedingung
genannt.

Bei 298,15 K sind beide Eltern sogar **dieselbe Publikation und dieselbe Probe**
(Miranda & Santos 2025, *Fluid Phase Equilib.* 597, 114458; 99 mass %,
0,0100 water mass %, dried by vacuum heating):

| | Wert | Methode |
|---|---|---|
| mu | (36,370 ± 0,760) mPa·s | Concentric cylinders viscometry |
| rho | (1280,9 ± 1,1) kg/m³ | Vibrating tube method |
| **nu** | **(28,3941 ± 0,5938) mm²/s** | **abgeleitet, nicht gemessen** |

Die Unsicherheit ist **fortgepflanzt** (quadratische Summe der relativen
Anteile, 2,09 %, dominiert von mu); die konservativere lineare Addition
(± 0,6177 mm²/s) wird mitgeliefert, weil mu und rho aus demselben Messgang
stammen und damit nicht streng unabhaengig sind.

Der Status ist `derived`, nicht `measured`. In `summary.csv` steht dafuer die
Spalte `value_kind`, im Referenzfeld stehen **beide** Eltern, und
`relative_spread` ist dort die fortgepflanzte Unsicherheit und keine
Literaturstreuung. Vollstaendiger Nachweis: `kinematic_viscosity.csv`.

**Quervergleich.** Die drei direkt gemessenen Quellen (Kapillarviskosimetrie,
ohne Reinheits- und Wasserangabe) geben 27,1 … 27,5 mm²/s. Der abgeleitete Wert
liegt **3,3 % darueber**. Das wird berichtet und nicht angepasst: die Abweichung
ist klein gegen die Literaturstreuung von mu (54 %) und damit durch
Probenunterschiede erklaerbar. Die direkt gemessenen Werte gehen in keine
Rechnung ein.

## Was eine Aenderung von gamma tut -- und was nicht

`impact.csv` fuehrte frueher eine Zeile
`maxwell_force_for_fixed_shape (linear in gamma)`. **Das war falsch.** Bei
festgehaltener Geometrie, festgehaltener angelegter Spannung und festgehaltener
Permittivitaetsverteilung folgt das Feld aus der Laplace- bzw.
Poisson-Gleichung, in der gamma nirgends vorkommt; die Maxwell-Traktion
eps0·E²/2 ist dann ein Funktional allein dieses Feldes und **invariant** unter
einer Aenderung von gamma.

gamma entscheidet, WELCHE Form im Gleichgewicht steht -- nicht, welche Kraft ein
gegebenes Feld auf eine gegebene Form ausuebt.

| Groesse | festgehalten | Gesetz | Exponent |
|---|---|---|---|
| Kapillardruckskala gamma/a | dimensionslose Form, a | gamma¹ | +1 |
| mechanische Last fuer dieselbe dimensionslose Form | dimensionslose Form | gamma¹ | +1 |
| Spannung fuer dieselbe dimensionslose Form | dimensionslose Form, Geometrie | √gamma | +½ |
| elektrische Bondzahl bei festem **Feld** | E, Geometrie | 1/gamma | −1 |
| **Maxwell-Traktion bei festem Geometrie/U/eps** | Geometrie, U, eps-Verteilung | gamma⁰ | **0** |

Die Zeile wird nicht geloescht, sondern als Invariante mit Exponent 0 gefuehrt,
damit die Korrektur sichtbar bleibt. Jede Zeile traegt `recomputed=no`: **keine**
ist eine neu gerechnete gekoppelte Simulation, sondern die Umskalierung einer
bereits gerechneten dimensionslosen Loesung. Die Tabelle steht als
`es::gamma_scaling_rows()` in `src/material_data.cpp` und wird in
`tests/test_material_data.cpp`, Abschnitt 5b, geprueft.

## Widersprueche, nebeneinander stehen gelassen

* Das IoLiTec-Blatt nennt 25,2 mPa·s; die Primaerliteratur liegt zwischen 31,7
  und 45 mPa·s. Das Blatt liegt unter dem gesamten Literaturband und nennt
  keinen Wassergehalt.
* Ein Leitfaehigkeitsdatensatz (Rilo et al. 2009) berichtet 5e-4 S/m, mehr als
  drei Groessenordnungen unter allen anderen. Er wird nicht geloescht.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_properties_vs_T.png` | jeder gemessene Punkt jeder Quelle gegen die Temperatur; gewaehlte Quelle hervorgehoben, Nicht-Umgebungsdruck- und frequenzaufgeloeste Punkte getrennt markiert |
| `fig2_scatter.png` | Literaturstreuung im Fenster 298,15 K ± 2 K, Quelle fuer Quelle, mit der Provenienz, die jede angibt |
| `fig3_impact.png` | wie eine Aenderung von gamma die P3a/P3b-Zahlen verschiebt -- und, als eigene Zeile, was sie nachweislich NICHT verschiebt |
| `fig4_kinematic_viscosity.png` | nu = mu/rho: der abgeleitete Wert mit beiden Eltern, der fortgepflanzten Unsicherheit und dem direkt gemessenen Band |

## Dateien

`kinematic_viscosity.csv` -- die Ableitung nu = mu/rho mit beiden Eltern, beiden
Einzelwerten, den geprueften Bedingungen und der fortgepflanzten Unsicherheit.
`sources.csv` -- jede Quelle mit Provenienz. `points.csv` -- jeder einzelne
Messpunkt mit Temperatur, Druck und Frequenz. `summary.csv` -- was eine Rechnung
benutzen darf und wie unsicher es ist. `impact.csv` -- die Skalierungen.
