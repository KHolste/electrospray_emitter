# 09 — P3a: statischer Kapillarmeniskus ohne elektrisches Feld

**Stand: 2026-08-29.** Dieses Dokument beschreibt genau den Schritt P3a: die
statische Form der freien Flüssigkeitsoberfläche an der Austrittskante, gehalten
**allein** von der Oberflächenspannung gegen einen vorgegebenen
Flüssigkeitsdruck. Alles, was darüber hinausgeht, ist hier ausdrücklich **nicht**
enthalten und darf aus keiner Zahl dieses Schrittes gelesen werden.

Code: `include/es/capillary.hpp`, `src/capillary.cpp`,
`include/es/liquid.hpp`, `src/liquid.cpp`, `apps/es_capillary.cpp`,
`tests/test_capillary.cpp`, `python/plot_capillary.py`.
Ergebnisse: `results/2026-08-29_p3a_capillary_meniscus/`.

---

## 9.1 Modellvertrag

Achsensymmetrische Young-Laplace-Gleichung, Meridiankurve in
Bogenlängenparametrisierung:

$$
\gamma\,\kappa = \Delta p_\mathrm{exit},
\qquad
\kappa = \frac{\mathrm{d}\psi}{\mathrm{d}s} + \frac{\sin\psi}{r},
$$

$$
\frac{\mathrm{d}r}{\mathrm{d}s} = \cos\psi,
\qquad
\frac{\mathrm{d}z}{\mathrm{d}s} = -\sin\psi .
$$

$\kappa$ ist die **Summe** der beiden Hauptkrümmungen, $s$ die Bogenlänge vom
Apex auf der Achse, $\psi$ der Tangentenwinkel. Die äußere Normale ist
$\mathbf n = (\sin\psi, \cos\psi)$ und zeigt aus der Flüssigkeit ins Vakuum.

**Warum keine Höhenfunktion $z(r)$.** Sie hat eine unendliche Ableitung, sobald
die Oberfläche senkrecht steht. Das passiert bereits an der hemisphärischen
Grenze dieses Problems und in jedem Taylor-Kegel. Eine $z(r)$-Formulierung
müsste dafür später vollständig umgebaut werden; sie wird deshalb gar nicht
erst gebaut.

### Vorzeichen

`delta_p_exit = p_liquid − p_vacuum` an der Austrittsebene $z = 0$.

| Vorzeichen | Wölbung | Apexhöhe |
|---|---|---|
| $\Delta p > 0$ | nach $+z$, zum Extraktor hinaus | $h > 0$ |
| $\Delta p = 0$ | **exakt eben** | $h = 0$ |
| $\Delta p < 0$ | in die Bohrung hinein | $h < 0$ |

Das ist getestet, nicht nur kommentiert (`test_capillary`, Abschnitte 1 und 4).

### Randbedingungen

* **Achse $r = 0$:** Regularität, $\psi(0) = 0$.
* **Kontaktlinie:** an der scharfen Austrittskante **gepinnt**,
  $r = \phi_2/2$, $z = 0$, direkt aus `DeviceParameters` über das Feature
  `PinnedContactEdge`. Kein Radius ist im Kapillarmodul fest codiert, und es
  gibt keine zweite Geometriebeschreibung.
* **Kein Kontaktwinkel.** Pinning und ein vorgeschriebener Young-Winkel sind
  zwei einander ausschließende Beschreibungen derselben Kante; beides zugleich
  überbestimmt das Problem. Der Löser lehnt die Kombination mit dem eigenen
  Status `ContactAngleAndPinningBothPrescribed` ab, und die Anwendung bricht bei
  `capillary.contact_angle_deg`, `wetting.contact_angle_deg` oder
  `liquid.contact_angle_deg` mit einer Erklärung ab. Es wird keine der beiden
  Beschreibungen stillschweigend bevorzugt.

### Achsengrenzwert — analytisch, nicht geraten

Bei $r = 0$ ist $\sin\psi/r$ ein $0/0$-Ausdruck. Regularität einer glatten
Fläche macht dort beide Hauptkrümmungen gleich, also

$$
\lim_{r\to0}\frac{\sin\psi}{r} = \frac{\mathrm{d}\psi}{\mathrm{d}s}
\quad\Longrightarrow\quad
\left.\frac{\mathrm{d}\psi}{\mathrm{d}s}\right|_\mathrm{Apex} = \frac{\kappa}{2}.
$$

Genau dieser Grenzwert steht im Integrator. Es wird nirgends durch $r = 0$
geteilt und es gibt keinen Kleinradius-Korrekturfaktor. Geprüft wird der
Grenzwert **an der numerischen Lösung**, nicht an der Formel des Lösers
(`test_capillary`, Abschnitt 3).

---

## 9.2 Gültiger Druckbereich

Dimensionsloser Druck

$$
\Pi = \frac{\Delta p\,a}{\gamma},\qquad a = \frac{\phi_2}{2}.
$$

Ohne Feld und ohne Schwerkraft ist die Lösung eine **Kugelkappe** mit
Kugelradius $R = 2\gamma/\Delta p$ und konstanter mittlerer Krümmung. Eine an
$a$ gepinnte glatte Kappe existiert genau für

$$
|\Pi| \le 2 \qquad\Longleftrightarrow\qquad |\Delta p| \le \frac{2\gamma}{a},
$$

denn sonst wäre $|R| < a$. $|\Pi| = 2$ ist die Halbkugel; dort steht die
Meridiankurve an der Kante senkrecht.

**Wie der Löser das feststellt.** Nicht aus der geschlossenen Formel, sondern an
der integrierten Form: Er marschiert in der Gesamtbogenlänge, und wenn die
Tangente senkrecht steht, bevor der Kontaktradius erreicht ist, wird der Punkt
der senkrechten Tangente durch Bisektion auf $\psi$ genau lokalisiert. Dort ist
der Radius maximal; ist dieses Maximum kleiner als $a$, existiert keine Lösung.
Damit hängt die Entscheidung nicht von der Schrittweite des Marsches ab.

### Statusvertrag

Kein einzelnes `converged`. `CapillaryStatus` unterscheidet:

| Status | Bedeutung |
|---|---|
| `Solved` | gelöst und bis zur geforderten Genauigkeit aufgelöst |
| `PressureOutsideCapillaryRange` | $\vert\Pi\vert > 2$: keine glatte gepinnte Fläche. **Keine Ersatzform, kein letzter Iterationsstand** |
| `HemisphericalLimit` | $\vert\Pi\vert = 2$ bis auf $10^{-9}$: die Halbkugel ist gültig und geschlossen bekannt, aber die Schießbedingung $r(L)=a$ wird dort tangential (doppelte Nullstelle); die numerische Form wird deshalb nicht als gelöst ausgegeben |
| `InvalidGeometry` | Pinningradius unbrauchbar |
| `InvalidLiquid` | Stoffdatensatz unbrauchbar; `why_unusable()` benennt den Grund |
| `ContactAngleAndPinningBothPrescribed` | siehe oben |
| `AccuracyNotReached` | Verfeinerung an der Obergrenze, geforderte Genauigkeit nicht erreicht |
| `ArclengthNotBracketed` | numerisches Scheitern der Bogenlängensuche, ausdrücklich getrennt vom physikalischen |
| `NotAttempted` | nichts gerechnet |

Im Sweep über $\Pi \in [-2{,}10; +2{,}10]$ wurden 22 von 421 Punkten so
abgelehnt; der größte gelöste $|\Pi|$ war 1,99.

---

## 9.3 Diskretisierung

Es gibt **keine Netzeingabe**. Einzige Eingabe ist die geforderte relative
Profilgenauigkeit (`capillary.accuracy`, bezogen auf den Bohrungsradius). Die
Intervallzahl der Bogenlängenintegration wird durch Verdopplung bestimmt, bis
die Änderung zur halben Auflösung darunter liegt; die tatsächlich verwendete
Zahl und der Schätzfehler stehen in jedem Ergebnis. Kein `h_tip`, kein `h_far`,
keine Knotenzahl.

Bei $10^{-11}$ wählt der Löser für die Referenzgeometrie 128 Intervalle
($\Pi = 0$) bis 2048 Intervalle ($|\Pi| = 1{,}98$).

Für die Netzstudie kann die Intervallzahl über `CapillaryRequest::forced_intervals`
erzwungen werden. Das ist ein **Studienschalter der Programmierschnittstelle**
und in keiner Konfigurationsdatei erreichbar.

---

## 9.4 Stoffdaten

Ein Wert ohne Herkunft ist keine Eingabe, sondern eine Behauptung. `LiquidProperties`
trägt deshalb Status und Fundstelle mit sich, und jede Ausgabe druckt sie.

| Status | Bedeutung |
|---|---|
| `verified` | Stoff, Temperatur und Primärquelle eindeutig belegt |
| `provisional` | benannte, aber nicht im Volltext geprüfte Quelle; trägt eine Empfindlichkeitsstudie |
| `illustrative` | Beispielwert ohne Primärquelle; trägt ausschließlich eine dimensionslose Demonstration |
| `unknown` | registriert, ohne Wert; die Verwendung ist ein Fehler |

### Was die Kunze-Unterlagen wirklich hergeben

**Belegt ist die Stoffidentität.** `KunzeFynn-2024-12-10.pdf`, Abschnitt 2.3.2
„Ionic liquids for electrospray thrusters“, gedruckte Seite 28 (PDF-Seite 36):
„EMI-BF4 and EMI-Im were also selected for this project.“ Die Tabelle
„List of Publications“, gedruckte Seite 30 (PDF-Seite 38), ordnet **EMI-BF4** den
Publikationen I–IV zu — also genau den geraden 10-µm-Kapillaren aus SU-8 bzw.
IP-Q, die die Referenzgeometrie modelliert. EMI-Im wird in den Publikationen V
und VI verwendet (Spiral- und extern benetzte Emitter).

**Nicht belegt sind die Zahlenwerte.** Der Volltext wurde auf numerische Werte
für Oberflächenspannung, Dichte, Viskosität und Leitfähigkeit durchsucht: die
Arbeit enthält keine Stoffwerttabelle für die Treibstoffe und keinen einzigen
solchen Zahlenwert. Abschnitt 2.3.2 vergleicht die beiden Flüssigkeiten
ausschließlich qualitativ („The surface tension of EMI-Im is slightly lower than
that of EMI-BF4“).

**Konsequenz.** Der eingebaute Datensatz `emibf4_illustrative()` trägt den Status
`illustrative`. $\gamma = 0{,}0452$ N/m und $\rho = 1279$ kg/m³ sind
**unverändert** aus der quellenlosen Tabelle in `src/fluid.cpp` übernommen;
nichts wurde aus dem Gedächtnis ergänzt und keine Stelle erfunden. Der Datensatz
trägt genau eine Art von Aussage: ein gekennzeichnetes Rechenbeispiel. Die
Prüfung des Lösers läuft unabhängig davon **dimensionslos** (Einheitsflüssigkeit
$\gamma = 1$ N/m, exakt).

**Viskosität, Leitfähigkeit und relative Permittivität** stehen in
`LiquidProperties::documented_only` — vorgemerkt, damit sie nicht zweimal
nachgeschlagen werden müssen, und vom Löser nicht gelesen. Sie dürfen nicht als
bereits berücksichtigte Physik ausgegeben werden: Viskosität braucht eine
Strömung, Leitfähigkeit einen Ladungstransport, $\varepsilon_r$ eine Feldkopplung.

---

## 9.5 Schwerkraft und Bond-Zahl

Die Schwerkraft ist **nicht** gekoppelt. Ob das zulässig ist, wird ausgerechnet
und berichtet, nicht unterstellt:

$$
\mathrm{Bo} = \frac{\rho\,g\,a^2}{\gamma}
= \frac{1279 \cdot 9{,}80665 \cdot (5{,}0\cdot10^{-6})^2}{0{,}0452}
= 6{,}94\cdot10^{-6}.
$$

Der hydrostatische Druck über eine Bohrungsradiushöhe ist damit $6{,}9\cdot10^{-6}$
der Kapillardruckskala $\gamma/a = 9040$ Pa, also 0,063 Pa. Über die größte hier
gerechnete Apexhöhe $|h| = 4{,}34$ µm sind es 0,054 Pa bzw. $6{,}0\cdot10^{-6}$
der Druckskala. Die Vernachlässigung ist für **diese** Geometrie damit
quantitativ gerechtfertigt; für eine millimetergroße Geometrie wäre sie es
nicht, denn Bo skaliert mit $a^2$.

Nicht behauptet wird damit der kapillare Aufstieg vom Reservoir bis zur
Austrittskante. Der **gefüllte Zulauf ist eine Voraussetzung** von P3a; der
Zustand stromauf geht ausschließlich über `delta_p_exit` ein. Die
P2c-Vorratsgeometrie bleibt unverändert und wird mechanisch nicht gelöst.

---

## 9.6 Prüfungen

`tests/test_capillary` (alle Toleranzen im `namespace tol` **vor** der ersten
Auswertung festgelegt und danach nicht angepasst):

| # | Prüfung | Ergebnis |
|---|---|---|
| 1 | $\Delta p = 0$ ergibt exakt ebene Fläche | $z \equiv 0$ und $\psi \equiv 0$ **bitgenau**, Toleranz 0 |
| 2 | gepinnte Kante wird getroffen | $\vert r_\mathrm{Kontakt}-a\vert/a \le 2\cdot10^{-15}$; $z$ exakt |
| 3 | Symmetrie und Krümmungsgrenzwert bei $r=0$ | $r(\mathrm{Apex})=0$, $\psi(\mathrm{Apex})=0$ exakt; $2\,\mathrm{d}\psi/\mathrm{d}s = \Delta p/\gamma$ auf $10^{-8}$ |
| 4 | Profil gegen Kugelkappe, $\Pi = \pm1{,}98 \dots \pm0{,}05$ | Normalabstand $\le 3{,}4\cdot10^{-13}\,a$ |
| 5 | Apexhöhe, Bogenlänge, Fläche, Volumen, $\psi$ an der Kante | relative Abweichung $\le 4\cdot10^{-13}$ |
| 6 | Young-Laplace-Residuum über die ganze Fläche | $\le 2{,}2\cdot10^{-8}$ (Grenze $10^{-4}$) |
| 7 | Netzkonvergenz, 5 Auflösungen (16 … 256) | beobachtete Ordnung Profil 3,8 (Grenze 3,5), Residuum 2,00 (Grenze 1,5) |
| 8 | Skalierung, $a$ um Faktor 74 und anderes $\gamma$ | $h/a$, $s/a$, $A/a^2$, $V/a^3$ gleich auf $10^{-14}$ |
| 9 | ungültige Stoffdaten, nicht darstellbarer Druck, Kontaktwinkel+Pinning | jeweils eigener Status, **keine Form zurückgegeben** |
| 10 | vollständige Baseline | 15/15 Tests grün; P0–P2c-Quellen unverändert |

Das Residuum wird **allein aus den Knotenkoordinaten** gebildet: Tangentenwinkel
aus den Sehnen, Meridiankrümmung aus deren Drehung. Es teilt keinen Code mit dem
Integrator, konvergiert deshalb zweiter statt vierter Ordnung — und ein Löser,
der nur seine eigenen diskreten Gleichungen erfüllt, kann es nicht bestehen.

Zusätzlich werden Fläche und Volumen **zweimal unabhängig** ausgewertet: mit der
ODE mitintegriert (vierte Ordnung) und aus der Polylinie mit den geprüften
Helfern aus `device_geometry.hpp` (zweite Ordnung).

---

## 9.7 Was P3a ausdrücklich nicht enthält

* kein elektrisches Feld, kein Maxwell-Druck, keine Kopplung an den
  FEM-Elektrostatiksolver, keine Betriebsspannung;
* keine Emission, keine Raumladung, kein Strahl;
* keine Strömung: kein viskoser Druckabfall, keine Speiseimpedanz;
* keine Zeitabhängigkeit und **keine Stabilitätsaussage** — dass eine Lösung hier
  existiert, sagt nichts darüber, ob sie stabil ist;
* kein Taylor-Kegel, kein Cone-Jet;
* keine Schwerkraft (siehe 9.5);
* keine Verrundung der Austrittskante (`reserved.edge_radius_inner` bleibt
  abgelehnt);
* keine mechanische Rechnung im Vorratsraum.

Der nächste Schritt (P3b) ist nicht begonnen.
