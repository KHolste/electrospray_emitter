# 11 — P0: numerische Bereinigung von P3b

**Stand: 2026-08-29.** Dieses Dokument beschreibt genau, was in P0 an P3b
korrigiert und gemessen wurde. Es enthält **keine neue Physik**. Aus keiner Zahl
darf mehr gelesen werden als „so genau ist die vorhandene Rechnung, und das ist
der Grund".

Code: `include/es/load_projection.hpp`, `src/load_projection.cpp`,
`apps/es_p3b_audit.cpp`, `tests/test_load_projection.cpp`,
`python/plot_p3b_audit.py`.
Ergebnisse: `results/2026-08-29_p0_p3b_audit/`.

---

## 11.1 Drei Dinge hießen „die Last"

In P3b hießen drei verschiedene Objekte „die Last", und nur eines davon ist das,
was der Kapillarlöser integriert. Sie sind jetzt getrennte, benannte, von außen
messbare Typen:

| Name | Was es ist | Wofür es verwendet wird |
|---|---|---|
| **Rohlast** | punktweiser Knotenwert $p_\mathrm{M}=\varepsilon_0E_n^2/2$ | für **nichts**; an der Kante nicht konvergent |
| **Segmentlast** | eine Konstante je Oberflächensegment, Kraft durch Rotationsfläche | die Kraftbuchhaltung; eine **Treppenfunktion** |
| **übergebene Last** | $p(\tau)=G'(\tau)/A'(\tau)$ aus den monotonen kubischen Interpolanten | **das**, was `solve_capillary_meniscus` bekommt |

Der Typ, der die dritte trägt, hieß `TauLoad` und lag in einem anonymen
Namensraum in `src/electrocapillary.cpp` — weshalb seine zentrale Zusage
(„stetig, und trägt die integrierte Maxwell-Kraft") aus keinem Test heraus
prüfbar war. Er heißt jetzt `es::ProjectedLoad`, liegt in
`include/es/load_projection.hpp` und ist inhaltlich **unverändert**: dieselben
Bins, dieselben Fritsch-Carlson-Steigungen, dieselbe Arithmetik in `at()`.

## 11.2 Geprüft an hergestellten Lasten mit bekanntem Integral

Eine Projektion ist Quadratur und Buchhaltung. Um sie zu prüfen, muss die
Eingabe eine sein, deren Integral man kennt. Also wird die Last **vorgegeben**
und kein Feld gelöst; die Segmentquadratur ist dabei buchstäblich die, die
`maxwell_load()` benutzt, weil beide denselben Helfer
`assemble_load_segments()` aufrufen — sonst prüfte das Audit eine Kopie des
Codes und nicht den Code.

Zwei Lasten, wie gefordert:

* glatt: $p=p_0\,(1+(r/a)^2)$ auf der ebenen Scheibe, $F=\tfrac32\pi p_0a^2$;
* integrable Singularität: $p=C\,d^{\beta}$ mit $-1<\beta<0$,
  $F = 2\pi C\,a^{2+\beta}/((1+\beta)(2+\beta))$.

**Die Kante der singulären Last.** Bei $d=0$ ist der punktweise Wert unendlich,
und keine Diskretisierung kann ihn tragen. Der **hergestellten Eingabe** wird
deshalb an genau diesem einen Knoten das lokale Mittel über das letzte halbe
Segment gegeben,
$p_\text{Kante} = (2/h)\int_0^{h/2} C d^\beta\,\mathrm dd = C(h/2)^\beta/(1+\beta)$.
Das ist eine **ausgesprochene Eigenschaft der Eingabe** und keine
Regularisierung der Projektion; es spiegelt, was das rekonstruierte FEM-Feld
dort ohnehin tut (ein Zellmittel, kein Punktwert). Sonst wird nirgends etwas
regularisiert.

### Was gemessen wurde

| Aussage | Messung |
|---|---|
| Segmentprojektion, glatte Last | Ordnung **1,98** (erwartet 2) |
| Segmentprojektion, $p=Cd^{\beta}$ | Ordnung **0,805 / 0,616 / 0,278** bei $\beta=-0{,}25/-0{,}44/-0{,}75$; erwartet $1+\beta = 0{,}75/0{,}56/0{,}25$ |
| Binkraft gegen Segmentkraft | $\le 2\cdot10^{-16}$ relativ |
| übergebene Last gegen $A'$ | $\le 2\cdot10^{-16}$ relativ — **exakt erhalten** |
| übergebene Last gegen $2\pi r\,\mathrm ds$ | $4\cdot10^{-6}$ relativ (Interpolationsordnung) |

Die zweite Zeile ist eine ehrliche Aussage und keine Schwäche: eine
Trapezregel über ein Segment, das eine integrable Singularität enthält, ist um
$O(h^{1+\beta})$ falsch. Genau diese Rate zeigt die Messung.

Beide Kraftzahlen stehen im Ergebnis, nicht nur die günstigere. Die Erhaltung
ist **exakt gegenüber dem rekonstruierten Flächenmaß** $A'$ und gilt gegenüber
dem echten Oberflächenelement nur zur Interpolationsordnung.

### Stetigkeit — ohne erfundene Schranke

Die Stetigkeit wird nicht an einer Schranke gemessen, sondern am **Abfall**: die
Last wird um einen Probenabstand $\delta$ links und rechts jedes Binrandes
ausgewertet, dann noch einmal mit $\delta/10$. Eine stetige Funktion mit Knick
gibt dann etwa $0{,}1$, eine Treppenfunktion etwa $1$. Gemessen: übergebene Last
$\approx 0{,}1$, Segmenttreppe $\approx 1{,}0$; der Sprung der übergebenen Last
ist zugleich mehr als vier Größenordnungen kleiner als der der Treppe.

**Die Rekonstruktion ist also tatsächlich stetig, nicht stufenförmig.** Der
in P3b als D2 dokumentierte Fehler ist damit nachgewiesen behoben und nicht nur
behauptet.

### Kontaktkante, Abdeckung, keine versteckte Kappung

* Die Segmente reichen bis $\tau=1$, also bis an die Kontaktlinie.
* Kein Bin ist leer; die Binflächen summieren sich auf $10^{-13}$ relativ zur
  Segmentfläche.
* Das kantennächste Bin trägt Kraft; es gibt **keine Ausschlusszone**.
* Das Maximum der übergebenen Last erreicht das Binmaximum; es gibt **keine
  Kappung**.
* Der punktweise Kantenwert wächst mit jeder Verfeinerung
  ($1{,}11\to1{,}51\to2{,}04\to2{,}77\cdot10^{5}$ Pa über vier Stufen) und der
  Kraftanteil des letzten Segments **fällt** dabei mit der Ordnung 0,556
  gegenüber $1+\beta=0{,}560$. Beides steht nebeneinander, weil es zusammen die
  Aussage ist: die Last ist integrabel, obwohl ihr Punktwert divergiert.

---

## 11.3 `kTolExclusion` — aufgelöst, nicht getauscht

P3b berichtete: das Kriterium `edge_gate::kTolExclusion` (5 %) werde um 11,3 %
verfehlt, ein anderes Kriterium entscheide. Zwei Kriterien nebeneinander
stehenzulassen ist keine Auflösung. P0 leitet die Sache her.

**Die Rechnung.** Das Kriterium misst

$$
\frac{\lvert F(d_0/2)-F(d_0)\rvert}{\lvert F(d_0/2)\rvert},
\qquad d_0=\texttt{kExclusionMid}=0{,}05\,a .
$$

Für $p=C\,d^{\beta}$ auf der ebenen Scheibe ist $F(d_0)$ in geschlossener Form
bekannt,

$$
F(d_0)=2\pi C\left[\frac{a\,(a^{1+\beta}-d_0^{1+\beta})}{1+\beta}
-\frac{a^{2+\beta}-d_0^{2+\beta}}{2+\beta}\right],
$$

und der Quotient oben enthält die **Netzweite überhaupt nicht**. Er hat einen
endlichen Grenzwert ungleich null, der allein von $\beta$ und $d_0$ abhängt.
Ein Kriterium, das verlangt, dass er klein wird, ist damit **kein
Konvergenztest**: es fordert von einer geometrischen Größe, dass sie
verschwindet.

**Der Nachweis, nicht die Behauptung.** `tests/test_load_projection.cpp`,
Abschnitt 4, prüft für $\beta=-0{,}25/-0{,}44/-0{,}75$:

| $\beta$ | geschlossene Form | gemessen (feinste Stufe) | Abweichung |
|---|---|---|---|
| $-0{,}25$ | 0,08105 | 0,08106 | $1{,}1\cdot10^{-4}$ |
| $-0{,}44$ | 0,11230 | 0,11232 | $2{,}3\cdot10^{-4}$ |
| $-0{,}75$ | 0,17929 | 0,17937 | $4{,}7\cdot10^{-4}$ |

und zusätzlich, dass die Messgröße unter Verfeinerung **stehen bleibt** (Drift
$\le 2{,}4\cdot10^{-3}$ von Stufe 256 auf 512), während der
Diskretisierungsfehler **derselben Last** weiter mit der Rate $h^{1+\beta}$
fällt. Das ist der entscheidende Gegensatz.

**Am realen Feld.** Für die ebene Meniskusform misst der Lauf $\beta=-0{,}4456$
und eine Halbierungsänderung von 0,1110; die geschlossene Form für dieses
$\beta$ ergibt 0,1133. Die 11,3 % aus P3b sind also **kein numerischer Zufall,
sondern der von $\beta$ festgelegte Wert.** Für die nach außen gewölbte Form
($\beta=+0{,}041$) ist die Größe 0,0364 und liegt damit unter der Schranke —
weil $\beta$ größer ist, nicht weil das Netz feiner wäre. Ein Kriterium, dessen
Erfüllung so von der Form abhängt, kann keine Netzkonvergenz prüfen.

**Was nicht geändert wurde.** Die Schranke `kTolExclusion` steht unverändert im
Code, wird gemessen und wird weiter berichtet. Geändert hat sich nur, dass
jetzt eine Rechnung mit Regressionstest dahintersteht, warum sie nicht
entscheidet.

**Das ersetzende Kriterium wurde gegen eine bekannte Wahrheit geprüft.**
`kTolLimitAgreement` vergleicht zwei Extrapolationen desselben Grenzwertes. Auf
einer hergestellten Last ist dieser Grenzwert **bekannt**: die Extrapolation
über die Ausschlussdistanz trifft ihn auf 1,05 % (Abschnitt 5 desselben Tests).
Damit ist das ersetzende Kriterium nicht nur gegen sich selbst geprüft.

---

## 11.4 Kanten-Gate mit zusätzlichen Netzstufen

Geprüft wurden nur die **ebene** und eine **nach außen gewölbte** Form
($\Pi=+1{,}5$), auf den Netzstufen 1 bis 5. Stufe 5 ist die feinste, die die
Bandfaktorisierung noch trägt (3,01 GiB); dafür wurde die Speichergrenze
`DielectricSetup::memory_cap_bytes` konfigurierbar gemacht — Vorgabewert
unverändert. Eine Stufe, die nicht rechenbar ist, wird **gemeldet** und nicht
still weggelassen.

Die drei Sichten stehen in `fig4_gate.png` nebeneinander: die direkte
Kraftfolge über die Netzstufen, die Extrapolation daraus (Aitken) und die
Extrapolation über die Ausschlussdistanz.

| Form | $\beta$ | Urteil | Ordnung der Gesamtkraft | geschätzter Fehler der feinsten Stufe |
|---|---|---|---|---|
| eben ($\Pi=0$) | $-0{,}446$ | Passed | **0,541** | 6,1 % |
| gewölbt ($\Pi=+1{,}5$) | $+0{,}041$ | Passed | 0,800 | 1,8 % |

**Die wichtigste Zahl dieses Abschnitts:** die integrierte Maxwell-Kraft
konvergiert für die ebene Form mit der Ordnung 0,541, und $1+\beta = 0{,}554$.
Das ist keine Enttäuschung, sondern die Bestätigung, dass die Randsingularität
die Konvergenzrate der **gesamten** gekoppelten Rechnung setzt. Zweite Ordnung
ist mit einer unverrundeten Kante nicht zu haben.

---

## 11.5 Gekoppelte Netzkonvergenz — Ziel verfehlt, und das steht so da

**Kein neuer Spannungssweep.** Zwei Betriebspunkte, 1000 V und 1400 V, die drei
vorhandenen Netzstufen plus eine feinere (1, 2, 3, 4). Für $h/a$, das
kantenferne $E_n$ und die integrierte Maxwell-Kraft wird eine
Richardson-Extrapolation auf dem Verfeinerungsverhältnis $\sqrt2$ gerechnet, und
zwar auf **zwei Fenstern** (Stufen 1–3 und 2–4), damit eine zufällig
auswertbare Ordnung eine nicht-asymptotische Folge nicht verdecken kann.

**Ziel, vor der Messung festgelegt:** geschätzter Diskretisierungsfehler der
feinsten Stufe unter 1 %.

| Punkt | Größe | Ordnung | geschätzter Fehler | Urteil |
|---|---|---|---|---|
| 1000 V | $h/a$ | — | — | `NotInAsymptoticRange` |
| 1000 V | $E_n$ kantenfern | — | — | `NotInAsymptoticRange` |
| 1000 V | Gesamtkraft | 0,581 | 6,07 % | `DiscretizationNotConverged` |
| 1400 V | $h/a$ | 0,953 | 5,67 % | `DiscretizationNotConverged` |
| 1400 V | $E_n$ kantenfern | 0,231 | 4,83 % | `DiscretizationNotConverged` |
| 1400 V | Gesamtkraft | 0,910 | 4,59 % | `DiscretizationNotConverged` |

**Das Ziel ist nirgends erreicht.** Bei 1000 V sind $h/a$ und $E_n$ nicht einmal
im asymptotischen Bereich: die Differenzen zwischen den Stufen schrumpfen nicht
monoton, also ist keine Ordnung beobachtbar und **kein Fehler schätzbar**. Was
dort bekannt ist, ist die bloße Änderung zwischen den beiden feinsten Stufen
(0,32 % für $h/a$), und das ist eine Änderung und keine Fehlerschranke. Der
Unterschied steht in den Feldnamen (`last_relative_change` gegen
`relative_error_finest`) und in der Abbildung.

**Konsequenz für alles, was P3b berichtet hat.** Die gekoppelten Ergebnisse sind
**qualitativ**. Keine Apexhöhe und keine Kraft aus diesem Ast trägt drei
Stellen. Die Enden der Äste (1439,5 V; 1236,3 V) sind damit erst recht keine
physikalische Aussage — sie waren es schon vorher nicht, aber jetzt ist auch die
Zahl selbst mit mehreren Prozent behaftet.

---

## 11.6 Darstellung

Drei Regeln, weil P3b alle drei verletzt hat:

* Eine **gescheiterte Rechnung** steht als `nan` und erscheint in der Abbildung
  als Lücke mit ihrem Status. In `coupled_convergence.csv` von P3b stand der
  fehlgeschlagene Unterrelaxationslauf als lauter Nullen; null ist ein
  physikalischer Wert und wurde als einer gelesen.
* Eine Abbildung heißt nur dann **„Konvergenz"**, wenn das 1-%-Ziel erreicht
  ist. `fig5_coupled.png` heißt deshalb „gekoppelte Netzstudie: das 1-%-Ziel ist
  NICHT erreicht".
* **Rohlast und übergebene Last** werden mit verschiedenen Namen, Farben und
  Linienarten gezeichnet und nie als dieselbe Kurve.

---

## 11.7 Was P0 ausdrücklich nicht enthält

Keine neue Physik. Keine Emission, keine endliche Leitfähigkeit, keine Strömung,
keine Raumladung, keine Zeitabhängigkeit, keine Stabilitätsaussage, kein
Taylor-Kegel. $\Delta p_\mathrm{exit}$ bleibt in diesem Lauf eine Eingabe (P1
ändert das getrennt). Die Stoffwerte bleiben `illustrative`.
