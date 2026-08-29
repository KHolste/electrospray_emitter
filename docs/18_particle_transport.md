# 18 — P7: Ionenbahnen und messbare Größen

**Stand: 2026-08-29.** Status: **`validated_subset`**. Die Bahnintegration, die
Treffererkennung und die Bilanz sind implementiert und geprüft. Die
Startverteilung ist **vorgeschrieben**, weil P5 blockiert ist; das Ergebnis ist
eine **Transportantwort** und keine Stromvorhersage.

Code: `include/es/particle_transport.hpp`, `src/particle_transport.cpp`,
`apps/es_trajectories.cpp`, `tests/test_particle_transport.cpp`,
`python/plot_trajectories.py`.
Ergebnisse: `results/2026-08-29_p7_trajectories/`.

---

## 18.1 Was gerechnet wird

Newtons Gleichung für eine **explizite** Spezies $(q,m)$ im Feld, das aus dem
achsensymmetrischen Volumennetz interpoliert wird — mit derselben bilinearen
Basis, mit der P6 die Ladung deponiert. Integrator: **velocity Verlet**,
zeitumkehrbar und ohne Energiedrift.

**Die Spezies ist Pflicht.** Ladung und Masse werden mit ihren Vorzeichen
übergeben; es gibt keine Vorgabespezies und keinen Pfad über einen Betrag. Eine
Spezies ohne Ladung wird abgelehnt.

**Die Achse ist eine Symmetrielinie, kein Wand.** Ein Teilchen, das $r<0$
erreicht, ist in drei Dimensionen durch die Achse gegangen; $r$ und $v_r$ werden
gespiegelt, und das zählt **nicht** als Treffer.

## 18.2 Die Pflichtprüfungen

| Prüfung | Ergebnis |
|---|---|
| gleichförmiges Feld gegen die geschlossene Form | Fehler auf dem **Rundungsboden** ($10^{-15}$…$10^{-14}$) |
| radiale Koordinate im gleichförmigen Feld | bewegt sich **nicht** ($<10^{-12}R$) |
| Energiegewinn $=q\,\Delta V$ | $3{,}1\cdot10^{-13}$ relativ |
| Energieerhaltung entlang der Bahn | $3{,}1\cdot10^{-13}$ |
| Polaritätsumkehr **beider** Vorzeichen | Bahn **identisch**, Flugzeit auf $10^{-12}$ gleich |
| nur die Spezies umgekehrt | läuft auf den Emitter zurück |
| Treffererkennung | alle vier Flächen einzeln getroffen |
| Bilanz | Schließfehler **exakt null** |
| Zeitordnung | **2,13** (im nicht gleichförmigen Feld) |
| Netzunabhängigkeit im exakt dargestellten Feld | $<10^{-6}$ |

### Zwei Befunde aus dem Test

**Im gleichförmigen Feld ist velocity Verlet exakt.** Die Beschleunigung ist
konstant, also *ist* der Schritt $x+\Delta t\,v+\tfrac12\Delta t^2a$ die
geschlossene Parabel. Der gemessene Fehler liegt beim Rundungsboden, und eine
dort abgelesene „Zeitordnung" misst die Akkumulation der Rundung — die erste
Fassung des Tests las so eine **negative** Ordnung ab. Die Ordnung wird deshalb
in einem **nicht** gleichförmigen Feld gemessen, durch Selbstkonvergenz gegen
den feinsten Schritt: 2,13.

**Der letzte Schritt muss ein Teilschritt sein.** Eine erste Fassung schnitt die
Position auf den Rand zurück, behielt aber die Geschwindigkeit des vollen
Schritts; die Austrittsgeschwindigkeit war dann um bis zu $2{,}5\cdot10^{-3}$
falsch. Jetzt wird auf den **Schrittanteil** $s$ bisektiert und ein
konsistenter Verlet-Schritt der Länge $s\,\Delta t$ ausgeführt. Ebenso wurde das
Potential am Austritt zunächst **außerhalb** des Netzes ausgewertet, wo
`locate()` null liefert — die Energiebilanz war dadurch um die volle
Beschleunigungsspannung falsch. Beide Fehler hat der Test gefunden.

**Und ein Vorzeichenfehler.** Die Erhaltungsgröße ist
$E=\tfrac12mv^2+q\varphi$; die erste Fassung schrieb ein Minuszeichen. Der
Energietest meldete daraufhin einen relativen Fehler von genau 2 — was eine
invertierte potentielle Energie erzeugt.

## 18.3 Treffererkennung und Bilanz

Vier Flächen: Emitter ($z=0$), Polymerwand ($r=R$), Extraktor ($z=Z$, außerhalb
der Blende) und die offene Blende ($z=Z$, $r\le r_\text{ap}$). Die Reihenfolge
der Prüfung ist festgelegt und begründet: die Stirnflächen zuerst, damit ein
Teilchen, das durch die Blende geht, nicht als Wandtreffer gezählt wird, nur
weil es im selben Schritt auch $r=R$ passiert hat.

$$
I_\text{extrahiert} + I_\text{abgefangen} + I_\text{noch fliegend}
= I_\text{gestartet}
$$

**exakt** — jedes Teilchen hat genau ein Schicksal. „Noch fliegend" wird
**getrennt** gezählt und nicht als verloren verbucht; der Test läuft absichtlich
einmal mit knappem Schrittbudget, damit diese Kategorie besetzt ist.

## 18.4 Was das Ergebnis ist und was nicht

`results/.../balance.csv` meldet für den Beispielaufbau (parallele Platten,
$U=-1000$ V, Blende $0{,}35R$):

| Fall | Transmission |
|---|---|
| Kation ohne Raumladung | 70,7 % |
| Kation mit vorgeschriebener Raumladung ($3\cdot10^{-14}$ C) | 31,7 % |
| Anion ohne Raumladung (bei umgekehrter Spannung) | 70,7 % |

**Das ist eine Transportantwort auf eine gesetzte Ladung, kein Betriebspunkt.**
Der gestartete Strom ist eine Eingabe, die Raumladung ist vorgeschrieben und
wird während des Fluges nicht aktualisiert — genau das macht den Vergleich der
ersten beiden Zeilen zu einer sauberen Differenz statt zu einer verschränkten.

Dass Kation und Anion dieselbe Transmission haben, ist **keine** Aussage über
Symmetrie der Physik: es ist die Folge davon, dass hier nur die Vorzeichen von
Ladung und Spannung zugleich umgekehrt wurden, während Masse, Barriere und
Speziesmischung in Wirklichkeit verschieden sind (P5, 16.2).

## 18.5 Was P7 ausdrücklich nicht enthält

* **Keine Emission.** P5 ist blockiert; es gibt keine physikalische
  Teilchenquelle und keine Zahl hier ist eine Stromvorhersage.
* **Kein Tropfenstrahl.** Der Tropfenzweig des Prototyps bleibt deaktiviert und
  wird hier nicht reaktiviert.
* **Keine selbstkonsistente Raumladungsschleife.** Die Teilchen ändern die
  Ladung während des Fluges nicht.
* Keine Magnetfelder, keine Stöße, keine Sekundäremission, kein Sputtern.
* Keine Gerätegeometrie: die Testdomäne ist ein Rechteck mit benannten Flächen.
  Der Anschluss an das P1/P2c-Gerät ist offen.
