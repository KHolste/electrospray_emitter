# 12 — P1: Druck- und Zulaufmodell

**Stand: 2026-08-29.** Dieses Dokument beschreibt genau, was `es::FeedRequest`,
`es::PressureBudget` und `apps/es_feed.cpp` rechnen — und, ebenso genau, was
sie **nicht** rechnen.

Code: `include/es/feed.hpp`, `src/feed.cpp`, `apps/es_feed.cpp`,
`tests/test_feed.cpp`, `python/plot_feed.py`.
Ergebnisse: `results/2026-08-29_p1_pressure_budget/`.

---

## 12.1 Modellvertrag

$$
\Delta p_\mathrm{exit}
= \underbrace{(p_\mathrm{reservoir}-p_\mathrm{vacuum})}_{\text{Antrieb}}
- \Delta p_\mathrm{hydrostatisch}
- \Delta p_\mathrm{viskos}
$$

mit

$$
\Delta p_\mathrm{hydrostatisch} = -\rho\,g_z\,(z_\mathrm{exit}-z_\mathrm{reservoir}),
\qquad
\Delta p_\mathrm{viskos} = \frac{8\mu L}{\pi R^4}\,Q .
$$

`PressureMode::Direct` behält die freie Eingabe von $\Delta p_\mathrm{exit}$ als
kontrollierten Modus. Das ist keine Bequemlichkeit: P3a und P3b wurden damit
gerechnet, und jedes vorhandene Ergebnis muss reproduzierbar bleiben. Im
direkten Modus sind die Einzelterme **nicht bekannt** und stehen als `nan`,
nicht als null.

## 12.2 Vorzeichen, Referenzhöhe, Strömungsrichtung — einmal festgelegt

**Geometrie.** $z$ läuft entlang der Emitterachse und wächst zum Extraktor hin.
Die Austrittsebene — die Ebene, in der die Kontaktlinie gepinnt ist — ist die
Referenzhöhe $z_\mathrm{exit}$. $z_\mathrm{reservoir}$ ist die Höhe, auf die
sich $p_\mathrm{reservoir}$ bezieht.

**Schwerkraft.** `gravity_axial` $=g_z$ ist die Komponente der
Erdbeschleunigung entlang $+z$ und ist **null als Vorgabe**, weil das
modellierte Gerät ein Raumfahrttriebwerk ist. Auf dem Tisch mit nach oben
zeigendem Emitter ist sie $-9{,}80665\ \mathrm{m/s^2}$. Aus
$p(z)=p_\mathrm{res}+\rho g_z(z-z_\mathrm{res})$ folgt der obige Ausdruck; für
einen nach oben zeigenden Emitter ist er $+\rho g H$, also ein **Verlust**, wie
er es sein muss.

**Strömungsrichtung.** $Q>0$ heißt: die Flüssigkeit strömt in $+z$, also vom
Vorrat zum Austritt — die Richtung, in die sie strömt, wenn der Emitter gespeist
wird. Der viskose Term ist dann ein Verlust. $Q<0$ ist zugelassen und kehrt das
Vorzeichen um; nichts an der Poiseuille-Beziehung stört sich daran.

**Jeder Term ist ein Abzug vom Antrieb.** Ein Haushalt, der einen davon
addierte, sähe genauso aus wie ein richtiger. `tests/test_feed.cpp` prüft
deshalb jeden Term **einzeln** auf sein Vorzeichen.

## 12.3 Was ausdrücklich nicht gerechnet wird

* **Kein Strömungslöser.** Eine geschlossene Widerstandsformel, stationär,
  laminar, ausgebildet, inkompressibel, konstanter Radius.
* **Kein kapillarer Aufstieg und keine bewegliche Kontaktlinie.** Beides wird
  nicht behauptet und nicht gerechnet.
* **Kein Young-Winkel im Kanal.** Der Kanal ist **voll gefüllt**; dort gibt es
  gar keine zweite freie Oberfläche, deren Kontaktwinkel man einsetzen könnte.
  Die einzige freie Oberfläche des Modells ist der an der Austrittskante
  gepinnte Meniskus, und dessen Kapillardruck rechnet P3a/P3b — er darf hier
  nicht ein zweites Mal addiert werden. Ein gesetzter Kontaktwinkel führt
  deshalb zu `MissingFeedInput` mit dieser Begründung und wird nicht
  stillschweigend übergangen.
* **Der Vorrat ist eine Randbedingung, kein Volumen.** Seine Geometrie geht
  nicht ein; es gibt keine Entleerung, keine freie Oberfläche darin und kein
  Gaspolster.
* **$p_\mathrm{reservoir}$ und $Q$ bleiben Eingaben.** Nichts hier berechnet
  sie. Ein Volumenstrom folgte aus der Emission, die diese Phase nicht hat.

## 12.4 Gültigkeitsgrenzen, geschlossen fehlschlagend

| Prüfung | Grenze | Status bei Verletzung |
|---|---|---|
| Reynoldszahl $\rho\bar u\,2R/\mu$ | $<2300$ | `NotLaminar` |
| Einlauflänge $0{,}06\,\mathrm{Re}\,2R$ | $<0{,}05\,L$ | `EntranceLengthNotShort` |
| Kanalradius und -länge | $>0$ | `ChannelGeometryInvalid` |
| $\mu$, $\rho$ vorhanden, wenn gebraucht | — | `MissingLiquidProperty` |

Fehlende Pflichtangaben liefern `nan`, **nie** null. Ein Nullwert wäre eine
Behauptung.

## 12.5 Prüfungen

`tests/test_feed.cpp` prüft gegen geschlossene Formeln und einen unabhängigen
Weg:

| Prüfung | Referenz | Abweichung |
|---|---|---|
| $R_h = 8\mu L/(\pi R^4)$ | die Definition, zweitgeschrieben | $0$ |
| **Simpson-Integral des Profils** $u(r)=\frac{\Delta p}{4\mu L}(R^2-r^2)$ über den Querschnitt | ergibt $Q$ | $<10^{-9}$ |
| Wandschubspannung $4\mu\bar u/R$ | gegen $\Delta p R/(2L)$ | $0$ |
| $R\to2R$ | $R_h/16$ | $0$ |
| $L\to2L$ | $2R_h$ | $0$ |
| $\rho g H$ beim Anheben um 10 mm | geschlossen | $0$ |
| Vorrat oberhalb des Austritts | Vorzeichenumkehr | $0$ |
| ohne Schwerkraft / ohne Strömung | Term **exakt** null | exakt |
| Haushalt affin in $Q$ | Verdopplung | $0$ |
| $\Pi=\Delta p_\mathrm{exit}/(\gamma/a)$ | geschlossen | $0$ |

Die zweite Zeile ist die eigentliche Prüfung: sie schließt den Kreis zwischen
dem Geschwindigkeitsfeld, aus dem der Widerstand hergeleitet ist, und dem
Widerstand selbst, ohne eine der beiden Formeln zweimal zu benutzen.

## 12.6 Kopplung an P3a/P3b

`coupled_to_p3a.csv` übergibt das gerechnete $\Delta p_\mathrm{exit}$ an den
**unveränderten** P3a-Kapillarlöser. Das ist keine neue Physik; es ist dieselbe
Rechnung mit einer anderen Herkunft der rechten Seite. Wo
$|\Delta p_\mathrm{exit}|>2\gamma/a$ liegt, existiert keine gepinnte statische
Form; der Punkt steht dann als Lücke mit seinem Status und nicht als Nullwert.
`PressureBudget::within_capillary_range` meldet das, **bevor** der Meniskuslöser
überhaupt gefragt wird.

## 12.7 Einordnung der Skalen (Rechenbeispiel, Stoffstatus `illustrative`)

Mit $a=5\ \mu$m, $\gamma/a = 9040$ Pa, $\mu = 37{,}1$ mPa·s, $\rho = 1279$
kg/m³, Kanal $R=5\ \mu$m, $L=300\ \mu$m:

| Skala | Wert | in $\gamma/a$ |
|---|---|---|
| $\gamma/a$ | 9040 Pa | 1 |
| Grenze der gepinnten Form $2\gamma/a$ | 18 080 Pa | 2 |
| viskoser Verlust bei $Q=10^{-13}$ m³/s | 4535 Pa | 0,502 |
| Hydrostatik über 1 mm (auf der Erde) | 12,5 Pa | 0,00139 |
| Hydrostatik über einen Bohrungsradius | 0,063 Pa | $6{,}9\cdot10^{-6}$ |

Ablesbar, und mehr behauptet dieser Abschnitt nicht: der **hydraulische
Widerstand** ist die Größe, die den Austrittsdruck auf der Kapillarskala bewegt;
die **Schwerkraft** ist es selbst im Laborfall nicht. Der hydraulische
Widerstand dieses Kanals ist $4{,}53\cdot10^{16}$ Pa·s/m³, so dass schon
$Q\approx4\cdot10^{-13}$ m³/s die ganze Spanne $|\Pi|\le2$ durchfährt.

Alle absoluten Drücke stützen sich auf $\mu$ und $\rho$ mit Status
`illustrative` und sind damit Rechenbeispiele. Die Verhältnisse und die
Exponenten sind es nicht.
