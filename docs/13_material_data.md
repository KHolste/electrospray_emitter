# 13 — P2: Stoffdatenvertrag für EMI-BF4

**Stand: 2026-08-29.** Dieses Dokument ist das Audit der Stoffdaten: was gesucht
wurde, was gefunden wurde, woher jede Zahl kommt, und was fehlt.

Code: `include/es/material_data.hpp`, `src/material_data.cpp`,
`src/material_data_emibf4.cpp` (**erzeugt**), `tools/fetch_material_data.py`,
`apps/es_material.cpp`, `tests/test_material_data.cpp`,
`python/plot_material.py`.
Ergebnisse: `results/2026-08-29_p2_material_data/`.

---

## 13.1 Was gesucht wurde

Dichte, Oberflächenspannung, dynamische Viskosität, elektrische Leitfähigkeit,
relative Permittivität, jeweils mit Temperaturabhängigkeit und Messbedingungen.

## 13.2 Wo gesucht wurde, und was das für Quellen sind

**NIST ILThermo v2.0**, NIST Standard Reference Database #147,
<https://ilthermo.boulder.nist.gov/>. Das ist eine *Zusammenstellung*: jeder
Datensatz darin trägt die vollständige Zitation der **Primärpublikation**, aus
der die Daten stammen, zusammen mit der Messmethode und der Probenbeschreibung,
die diese Publikation berichtet hat. Abgefragt wird die Zusammenstellung;
aufgezeichnet und zu zitieren ist die Primärquelle.

Die Datenbank ist frei abfragbar. `tools/fetch_material_data.py` holt sie über
die dokumentierte Abfrageschnittstelle und schreibt die Zahlen ohne
Zwischenschritt in die erzeugte C++-Datei. **Keine Zahl ist von Hand getippt,
keine stammt aus einem Suchtreffer-Snippet, keine aus dem Gedächtnis.**

**IoLiTec, Technical Data Sheet IL-0006.** Das Herstellerdatenblatt wurde als
PDF geöffnet und gelesen (Revision Date 2/13/2012, Date Issued 11/2/2012).
Übernommen sind daraus, mit ihren Temperaturen: Dichte 1,282 g/cm³ (24 °C),
Viskosität 25,2 mPa·s (25 °C), Leitfähigkeit 14,1 mS/cm (25 °C), Reinheit
>98 %. Das Blatt nennt **keine** Oberflächenspannung, **keine** Permittivität,
**keinen** Wassergehalt und **keine** Messmethode; diese Felder bleiben deshalb
leer statt geraten.

### Was NICHT als Quelle gilt

Suchtreffer-Zusammenfassungen, Datenbanken hinter Bezahlschranken, deren
Wertstelle nicht gelesen werden konnte, und alles aus dem Gedächtnis. Die
Recherche hat mehrere solcher Treffer geliefert (ResearchGate-Vorschauen,
pubs.acs.org-Abstracts); aus keinem davon ist eine Zahl übernommen worden.

## 13.3 Was gefunden wurde

| Größe | Quellen | Punkte | Band bei 298,15 K | Streuung |
|---|---|---|---|---|
| Oberflächenspannung | 14 | 87 | 0,0443 … 0,0546 N/m | 19 % |
| Dichte | 80 | 849 | 1277,6 … 1290,0 kg/m³ | 1,0 % |
| dynamische Viskosität | 32 | 256 | 0,0252 … 0,0450 Pa·s | 54 % |
| elektrische Leitfähigkeit | 24 | 326 | 0,91 … 1,63 S/m (bei 298,15 K; ein Ausreißer 5·10⁻⁴ liegt knapp ausserhalb des 2-K-Fensters) | 46 % |
| kinematische Viskosität | 3 | 22 | — | — |
| relative Permittivität | 4 | 21 | 12,8 … 13,6 (statisch) | — |

Die genauen Zahlen stehen in `results/2026-08-29_p2_material_data/summary.csv`;
sie werden hier nicht doppelt gepflegt.

## 13.4 Die Auswahlregel

Ein Löser braucht **eine** Zahl. Sie wird nach einer vorab festgelegten,
mechanisch angewandten Regel gewählt — nicht gemittelt:

1. Die Quelle muss **Methode, Reinheit UND Wassergehalt** angeben.
2. Sie muss Punkte bei **Umgebungsdruck** haben und darf **nicht
   frequenzaufgelöst** sein.
3. Unter den verbleibenden: diejenige, deren Umgebungsdruck-Temperaturbereich
   298,15 K enthält und die die meisten Umgebungsdruckpunkte hat.
4. Gleichstand: die kleinere angegebene Unsicherheit.

Erfüllt keine Quelle Punkt 1, wird **nichts** gewählt und die Größe ist
`MissingMaterialData`. Es gibt keinen stillen Ersatzwert.

Die Regel ist nicht darauf abgestimmt, eine bestimmte Zahl zu liefern. Sie steht
in `include/es/material_data.hpp`, wird von `tools/fetch_material_data.py`
angewandt und in `tests/test_material_data.cpp` für jede getroffene Auswahl
unabhängig nachgeprüft.

### Zwei Fallen, die die Regel abfängt

**Druck.** ILThermo enthält für die Dichte auch isochore PVT-Messungen bei
1–61 MPa (Klomfar et al. 2012). Ohne Druckprüfung gewinnt dieser Datensatz die
Auswahl, weil er die meisten Punkte hat — und liefert 1299,6 kg/m³ statt der
rund 1281 kg/m³ bei Umgebungsdruck, also 1,5 % daneben. Punkte abseits des
Umgebungsdrucks bleiben im Datensatz (sie sind echte Messungen), tragen aber
keinen Umgebungsdruckwert und gehen auch nicht in das Streuband ein.

**Frequenz.** Die einzige mehrpunktige Permittivitätsquelle (Bennett et al.
2019) misst von 1 GHz bis 18 GHz, wobei $\varepsilon_r$ von 11,7 auf 7,7 fällt.
Das ist eine dielektrische Dispersionsmessung und **kein Gleichstromwert**.
Frequenzaufgelöste Quellen sind deshalb von der Auswahl und vom Streuband
ausgeschlossen, und `PropertyPoint` trägt die Messfrequenz mit, damit das
sichtbar bleibt.

## 13.5 Der wichtigste Befund: γ

$\gamma$ ist die einzige Stoffgröße, die in das P3a/P3b-Gleichgewicht eingeht,
und sie geht **linear** ein. Jede Druckskala ist proportional zu $\gamma$, jede
Spannung zu $\sqrt\gamma$.

**Der bisher benutzte Wert 0,0452 N/m liegt unter jeder Quelle, die Reinheit und
Wassergehalt angibt** — die niedrigste davon ist 0,0501 N/m. Er liegt *nicht*
unter dem gesamten Band: dessen untere Kante ist 0,0443 N/m aus einer einzelnen
Kapillaraufstiegsmessung ohne Reinheits- und Wasserangabe (siehe unten). Die
schärfere Aussage ist die zutreffende, und sie ist die, auf die es ankommt: nur
Quellen mit dokumentierter Probe können eine Zahl tragen. Der Wert stammt
unverändert aus der quellenlosen Tabelle in `src/fluid.cpp`. Die bestbelegte Quelle — Souckova, Klomfar, Patek (2011),
*Fluid Phase Equilib.* 303(2), 184–190, 99,8 mass % Reinheit, 0,0164 mass %
Wasser (Karl-Fischer), zwei unabhängige Methoden (Wilhelmy-Platte und
du-Noüy-Ring), 14 Punkte von 288 bis 356 K, angegebene Unsicherheit 5–6·10⁻⁵
N/m — liefert 0,05397 bzw. 0,05401 N/m bei 298 K. Die beiden Methoden derselben
Arbeit stimmen auf 4·10⁻⁵ N/m überein.

Der einzige Literaturwert in der Nähe von 0,0452 ist 0,0443 N/m aus Martino, de
la Mora, Yoshida, Saito, Wilkes (2006) — einer Kapillaraufstiegsmessung ohne
angegebene Reinheit und ohne Wassergehalt, aus der Electrospray-Literatur
selbst. Sie wird **nicht** weggelassen und **nicht** gemittelt; sie steht als
untere Kante des Bandes im Datensatz, und sie ist auch der Grund, warum die
Streuung mit 19 % so groß ausfällt. Ohne sie liegen die elf übrigen Quellen
zwischen 0,0501 und 0,0546 N/m.

**Konsequenz, quantitativ:** `impact.csv` rechnet die Skalierung aus. Der
Faktor auf jede Druckskala ist rund 1,19, auf jede Spannung rund 1,09. Das ist
größer als jeder Diskretisierungsfehler, den P0 gemessen hat.

**Der bisherige Datensatz bleibt `illustrative`.** `emibf4_illustrative()` in
`src/liquid.cpp` ist unverändert; er wird nicht heimlich auf den belegten Wert
gesetzt. Der belegte Datensatz ist ein zweiter, `emibf4_sourced()`, und jede
Rechnung, die ihn benutzt, sagt das in ihrem Ausgabekopf.

## 13.6 Was fehlt, und was deshalb geschlossen fehlschlägt

| Größe | Status | Grund |
|---|---|---|
| relative Permittivität | `MissingMaterialData` | keine der vier Quellen nennt Reinheit und Wassergehalt; die mehrpunktige ist zudem frequenzaufgelöst (1–18 GHz) |
| kinematische Viskosität | `MissingMaterialData` | keine Quelle nennt Methode, Reinheit und Wassergehalt zugleich. Sie wird **nicht** aus $\mu$ und $\rho$ zusammengerechnet: welche Dichte der jeweilige Autor benutzt hat, ist nicht festgehalten |

Für P3b heißt das: die relative Permittivität der **Flüssigkeit** ist nicht
belegt. Sie geht dort auch nicht ein — die Flüssigkeit ist ein idealer Leiter —
aber jedes Modell mit endlicher Leitfähigkeit (Punkt 3) braucht sie, und dort
wird sie geschlossen fehlschlagen, bis eine Quelle sie trägt.

Die relativen Permittivitäten der **Polymere** (SU-8, IP-Q) sind eine andere
Baustelle und liegen in `es::MaterialLibrary` mit ihrem eigenen Status; P2 fasst
sie nicht an.

## 13.7 Widersprüche, nebeneinander stehen gelassen

* **Viskosität.** Das IoLiTec-Datenblatt nennt 25,2 mPa·s bei 25 °C. Die
  Primärliteratur bei 298,15 K liegt zwischen 31,7 und 45 mPa·s. Das
  Datenblatt liegt also **unter** dem gesamten Literaturband. Es nennt keinen
  Wassergehalt, und Wasser senkt die Viskosität ionischer Flüssigkeiten stark;
  mehr lässt sich dazu ohne die Angabe nicht sagen. Beide Werte stehen im
  Datensatz, mit ihrem jeweiligen Status.
* **Leitfähigkeit.** Ein Datensatz (Rilo et al. 2009, „crison
  conductivimeter") berichtet 5·10⁻⁴ S/m, mehr als drei Größenordnungen unter
  allen anderen. Der Wert wird **nicht gelöscht**; er steht im Datensatz und
  vergrößert das gemeldete Streuband, was genau die richtige Wirkung ist:
  ein Band, das ihn verschweigt, wäre schmaler, als der Kenntnisstand hergibt.
* **Oberflächenspannung.** Siehe 13.5.

## 13.8 Reproduzierbarkeit

`src/material_data_emibf4.cpp` ist **erzeugt und versioniert**. Der Lauf

```sh
python tools/fetch_material_data.py
```

lädt die Daten neu und schreibt die Datei neu; findet `git diff` danach keinen
Unterschied, ist die Zusammenstellung unverändert. Der Zwischenspeicher
`tools/_ilthermo_cache/` ist nicht versioniert — jede Zahl steht mit ihrer
Quelle in der erzeugten Datei, und eine fremde Datenbank gehört nicht als
Ganzes in dieses Repositorium.

Die Bibliothek bleibt abhängigkeitsfrei: die erzeugte Datei ist reines C++ ohne
Laufzeit-Dateizugriff.
