# 10 — P3b: selbstkonsistentes statisches Elektro-Kapillargleichgewicht

**Stand: 2026-08-29.** Dieses Dokument beschreibt genau den Schritt P3b: die
Kopplung des geprüften P3a-Kapillarsolvers an die dielektrische P2c-Elektrostatik
über den Maxwell-Druck auf der freien Oberfläche. Alles darüber hinaus ist
ausdrücklich **nicht** enthalten und darf aus keiner Zahl gelesen werden.

Code: `include/es/electrocapillary.hpp`, `src/electrocapillary.cpp`,
`apps/es_electrocapillary.cpp`, `tests/test_electrocapillary.cpp`,
`python/plot_electrocapillary.py`.
Ergebnisse: `results/2026-08-29_p3b_electrocapillary/`.

---

## 10.1 Modellvertrag

$$
\gamma\,\kappa(s) \;=\; \Delta p_\mathrm{exit} \;+\; p_\mathrm{M}(s),
\qquad
p_\mathrm{M} \;=\; \tfrac12\,\varepsilon_0 E_n^2 ,
$$

mit $\kappa$ der Summe der Hauptkrümmungen wie in P3a und $E_n$ dem
**vakuumseitigen** Normalfeld auf der freien Oberfläche der Flüssigkeit, die als
idealer Leiter auf $V_\mathrm{emitter}$ behandelt wird. Die Elektrostatik ist das
unveränderte dielektrische P2c-Problem: die Polymerkörper polarisieren und sind
nie Elektroden, die Metallisierung des Extraktorträgers liegt auf
$V_\mathrm{extraktor}$.

### Vorzeichen — hergeleitet, nicht behauptet

Auf der Oberfläche eines Leiters steht das Feld senkrecht,
$\mathbf E = E_n\mathbf n$ mit $\mathbf n$ der ins Vakuum zeigenden Normalen, und
die Flächenladung ist $\sigma = \varepsilon_0 E_n$. Die Traktion ist das Produkt
aus Flächenladung und **Mittelwert** des Feldes auf beiden Seiten — innen null,
außen $E_n$ —, also

$$
\mathbf t \;=\; \sigma\,\frac{E_n}{2}\,\mathbf n
\;=\; \frac{\varepsilon_0 E_n^2}{2}\,\mathbf n
\;=\; p_\mathrm{M}\,\mathbf n,
\qquad p_\mathrm{M}\ge 0 .
$$

Die Last zieht die Oberfläche also **immer nach außen**, unabhängig vom
Vorzeichen von $E_n$, und geht deshalb wie eine Erhöhung des Innendrucks in die
Krümmungsgleichung ein. Zwei Folgerungen werden geprüft statt behauptet:

* bei $\Delta p_\mathrm{exit}=0$ und $V\ne0$ wölbt sich der Meniskus nach außen
  ($h>0$) — gemessen;
* die Umkehr **aller** angelegten Potentiale kehrt $E_n$ im Vorzeichen um und
  lässt $p_\mathrm{M}$ und damit die Form unverändert — gemessen: die
  Formunterschiede sind exakt null.

Weitere unabhängige Prüfungen der Feld- und Flussauswertung, alle am ebenen
Referenzfall:

| Prüfung | Messung |
|---|---|
| $E_n(+V)+E_n(-V)$, relativ | 0 (exakt) |
| $F(+V)-F(-V)$, relativ | 0 (exakt) |
| $F(2V)/F(V)$ gegen 4 | 0 (exakt) |
| Ladung aus $\sigma=\varepsilon_0E_n$ gegen die FEM-Knotenreaktionen | 1,4·10⁻³ |
| $\lvert E_t\rvert/\lvert E\rvert$ auf der Äquipotentialfläche, kantenfern | 3,0·10⁻² |
| $\varphi=\mathrm{const}$ bei gleichen Potentialen und natürlichem Rand | 1,7·10⁻¹¹ relativ |

Der Ladungsvergleich lässt den **Kontaktknoten** in beiden Summen weg, und das
muss er: sein Element-Patch reicht auch in das Emitterdielektrikum, seine
Knotenreaktion ist also die Ladung zweier verschiedener Flächen zusammen,
während das Flächenintegral nur den Meniskus kennt.

Der Test mit gleichen Potentialen braucht den **natürlichen** Fernrand. Am
offenen Rand gilt $\varphi\to0$ im Unendlichen; zwei Elektroden auf demselben
Potential ungleich null tragen dort Ladung und erzeugen sehr wohl ein Feld. Mit
der asymptotischen Bedingung hätte dieser Test nichts geprüft.

---

## 10.2 Das bewegliche Netz

Das P2c-Netz ist ein Tensorgitter in $(r_\mathrm{ref}, z)$ mit einem rein
**radialen** Warp, der innerhalb der Bohrung die Identität ist. Die Spalten
$i\le i_\mathrm{bore}$ sind daher senkrechte Geraden und die Zeile $j_\mathrm{tip}$
ist die ebene Flüssigkeitsoberfläche. P3b fügt einen **axialen** Warp hinzu, und
nur auf diesen Spalten:

* die Zeile $j_\mathrm{tip}$ wird von $z=0$ auf $z_\mathrm{Meniskus}(r_i)$ gezogen;
* zwei vorhandene Gitterzeilen, die erste unterhalb $-1{,}5a$ und die erste
  oberhalb $+1{,}5a$, werden festgehalten und bilden die Bandränder;
* dazwischen ist die Abbildung stückweise linear in $z$ und stetig.

Außerhalb der Bohrung bewegt sich **kein** Knoten, und bei $r=r_\mathrm{bore}$ ist
die Meniskushöhe wegen des Pinnings null, so dass beide Hälften dort exakt
zusammenpassen. Die Kontaktspalte wird gar nicht erst gerechnet, sondern
übersprungen: der Warp wäre dort die Identität, als Arithmetik aber nur bis auf
Rundung, und die Zusage „außerhalb der Bohrung bewegt sich nichts“ wäre dann auf
10⁻²¹ statt exakt richtig.

**Zellen können nicht invertieren.** In den verformten Spalten hängt $r$ nur von
$i$ ab, die Jacobi-Determinante ist also $(\partial r/\partial\xi)\,
(\partial z/\partial\eta)$ mit zwei positiven Faktoren. Gemessen wird es
trotzdem. Über alle geprüften Formen und Netzstufen: kleinste
Jacobi-Determinante > 0, **0 invertierte Zellen**, Kontaktlinie und Apex auf
Rundungsgenauigkeit getroffen, Fehler des Flüssigkeitsvolumens gegen die
geschlossene Form (Bohrungssäule + Rotationsvolumen des Meniskus) ≤ 2,3·10⁻¹¹.

Die **Zellzuordnung** bleibt unverändert: P2c weist die Materialien aus Indizes
allein zu, die Zellen unterhalb der Zeile $j_\mathrm{tip}$ innerhalb der Bohrung
sind also Flüssigkeit, wohin diese Zeile auch gelegt wird. Randbedingungen,
Audit und Assemblierung sind damit die P2c-Ursprünglichen, unangetastet.

Die Bandränder sind **vorhandene** Zeilen, keine neuen. Deshalb ist das Netz für
die ebene Oberfläche bitgleich das P2c-Netz — was die Nullfeld-Rückprüfung erst
zu einer Prüfung macht statt zu einem Vergleich zweier verschiedener Netze.

**Punktlokalisierung.** Die generische `locate()` in `axisym_fem.hpp` setzt
Zeilen auf konstantem $z$ voraus, die das verformte Netz innerhalb der Bohrung
nicht hat. `locate_meniscus()` macht dasselbe für dieses Netz und ist exakt:
außerhalb der Bohrung **sind** die Zeilen eben, innerhalb sind die Spalten
senkrecht, beide Fälle reduzieren sich auf eine Bisektion plus eine lineare
Abbildung. Geprüft durch Rücktransformation: der bilinear zurückgerechnete Punkt
trifft den abgefragten auf 10⁻¹² des Bohrungsradius.

`AxisymProblem::require_level_rows` ist neu und standardmäßig `true`; die
Assemblierung braucht ebene Zeilen nie, nur die Punktlokalisierung. Jeder
bestehende Aufrufer bleibt damit unverändert.

---

## 10.3 Das Kanten-Gate

Die Kontaktlinie sitzt auf einer unverrundeten Leiter/Dielektrikum/Vakuum-Kante.
Dort ist das Feld singulär und **kein punktweiser Wert netzkonvergent**. Ob die
Flächenlast überhaupt eine brauchbare Größe ist, wird deshalb **vor** jeder
Kopplung entschieden, an Messungen.

### Grenzen, vor der Messung festgelegt

| Größe | Grenze |
|---|---|
| kantenferne Last ($d\ge0{,}1a$), zwei feinste Netzstufen | 2·10⁻² |
| integrierte Gesamtkraft, zwei feinste Netzstufen | 5·10⁻² |
| Singularitätsexponent $\beta$ in $p_\mathrm{M}\sim d^\beta$ | $>-1$ (integrierbar) |
| Abgleich zweier unabhängiger Extrapolationen der Grenzkraft | 5·10⁻² |

### Ergebnis: das Gate hängt von der **Form** ab

| $\Pi$ | $\psi$ an der Kante | $\beta$ | Urteil |
|---|---|---|---|
| $+1{,}90$ | $+71{,}8°$ | $+0{,}308$ | Passed |
| $+1{,}50$ | $+48{,}6°$ | $+0{,}054$ | Passed |
| $+0{,}50$ | $+14{,}5°$ | $-0{,}292$ | Passed |
| $0$ | $0°$ | $-0{,}444$ | Passed |
| $-0{,}50$ | $-14{,}5°$ | $-0{,}611$ | **FailedExclusion** |
| $-1{,}50$ | $-48{,}6°$ | $-1{,}095$ | **FailedNotIntegrable** |

Das ist keine numerische Laune, sondern die Geometrie der Kante. Bei einer
ebenen oder nach außen gewölbten Oberfläche ist die Leiterkante von der
Vakuumseite aus konvex oder flach; die Singularität ist schwach oder
verschwindet ganz ($\beta>0$ bei $\Pi\ge1{,}5$: die Last **fällt** zur Kante hin
ab). Wird die Oberfläche in die Bohrung gezogen, wird die Leiterkante
**einspringend**, die Singularität verstärkt sich, und bei $\Pi=-1{,}5$ ist der
gemessene Exponent nicht mehr integrierbar.

**Konsequenz, mechanisch angewandt:** gekoppelt wird ausschließlich für
Kontaktwinkel im Bereich, in dem das Gate bestanden wurde, also
$\psi\in[0°,\;71{,}8°]$. Jeder Ast, dessen feldfreie Startform außerhalb liegt,
wird mit `EdgeLoadNotWellPosed` abgelehnt, und jeder konvergierte Punkt eines
Astes wird nachträglich gegen dieselbe Schranke geprüft.

### Was die Kante punktweise tut, und dass man es sieht

Der punktweise Knotenwert $p_\mathrm{M}$ an der Kontaktlinie **wächst mit jeder
Verfeinerung** (ebene Form: 1,66 → 1,99 → 2,36 → 2,77·10⁴ Pa über die Stufen
1 bis 4). Er wird berichtet, damit man das sieht, und er wird **nirgends
verwendet**. Nichts wird abgeschnitten, kein Maximalwert ist codiert, keine
Ausschlusszone ist frei gewählt, kein Kantenradius erfunden.

Die Kraft **außerhalb einer festen Ausschlussdistanz** konvergiert dagegen
ausgezeichnet: bei $d\ge0{,}1a$ ändert sie sich von Stufe 3 auf Stufe 4 um
1·10⁻⁴ relativ.

---

## 10.4 Die Lastprojektion

**Konservativ, segmentweise.** Für jedes Oberflächensegment wird die
Normalkraft mit dem einseitigen Feld integriert,

$$
\Delta F_k=\int_{\text{Segment}} \frac{\varepsilon_0E_n^2}{2}\,2\pi r\,\mathrm ds ,
$$

und der Segmentdruck ist $\Delta F_k$ geteilt durch die Rotationsfläche des
Segments. Die Summe über die Segmente **ist** die integrierte Maxwell-Kraft, per
Konstruktion. Nahe der Kante bleibt das Segmentmittel einer integrierbaren
Singularität endlich, während der punktweise Wert es nicht tut.

**Stetig, für die Übergabe an den Kapillarlöser.** Aus den Bindaten werden die
kumulierte Kraft $G$ und die kumulierte Fläche $A$ monoton kubisch interpoliert
und $p(\tau)=G'(\tau)/A'(\tau)$ gesetzt. Damit gilt für jedes Bin
$\int p\,\mathrm dA=\int G'\,\mathrm d\tau = G(\tau_{b+1})-G(\tau_b)$, die
integrierte Kraft bleibt also **exakt** erhalten, und die rechte Seite der ODE
ist stetig.

**Parametrisierung über die normierte Bogenlänge** $\tau=s/L$. Die Last lebt auf
der Oberfläche, und die Oberfläche ist das Gesuchte; nur eine Parametrisierung,
die das Ändern der Form übersteht, kann die Last mitführen. $s$ tut das nicht
(sein Bereich verschiebt sich), $r$ auch nicht (seine Abbildung auf die Fläche
wird singulär, sobald die Tangente senkrecht steht).

---

## 10.5 Der gekoppelte Fixpunkt und die Fortsetzung

Je Betriebspunkt: Form vorgeben → konformes Volumennetz → dielektrische
Elektrostatik → vakuumseitiges $E_n$ → konservative Projektion → Young-Laplace
mit ortsabhängiger Last → wiederholen. Die Last wird unterrelaxiert
($\omega=0{,}5$).

Grenzen, vor der Messung festgelegt und getrennt ausgegeben:

| Größe | Grenze |
|---|---|
| maximale Formänderung / $a$ | 10⁻⁶ |
| Änderung des Maxwell-Drucks / $(\gamma/a)$ | 10⁻⁴ |
| Treffer der gepinnten Kontaktlinie | 10⁻¹² |
| Young-Laplace-Residuum, **kantenfern** | 10⁻³ |
| FEM-Residuum, Netzqualität | berichtet |

Kein einzelnes boolesches `converged`; `CouplingStatus` unterscheidet neun
Fälle.

**Rückwärtskompatibilität als Pflicht-Gate.** Bei
$V_\mathrm{emitter}=V_\mathrm{extraktor}=0$ ist die Lasttabelle identisch null;
dann wird sie gar nicht erst übergeben und das Problem als das P3a-Problem
gelöst — gleiche rechte Seite, gleiche geforderte Genauigkeit, gleiche Antwort
Bit für Bit. Gemessen über vier Drücke: Knotenabstand **exakt null**.

**Fortsetzung.** Verfolgt wird ausschließlich der Ast, der die feldfreie
P3a-Lösung enthält. $|V|$ wird in Schritten erhöht, jeder Schritt startet von der
vorherigen konvergierten Form; ein gescheiterter Schritt wird halbiert, und die
Fortsetzung endet, wenn die Schrittweite unter ihre Untergrenze fiele. **Ein
gescheiterter Punkt wird nie als Lösung gespeichert.**

> **Was das Ende eines Astes ist und was nicht.** Es ist die Stelle, an der
> *dieser Löser* mit dieser Schrittweite und diesen vorab festgelegten Grenzen
> stehen bleibt. Es ist **kein** Emissionsbeginn, **kein** Taylor-Kegel-Onset
> und **keine** Aussage über dynamische Stabilität. Nichts davon wird in P3b
> gerechnet.

---

## 10.6 Gefundene Fehler

Die Aufgabe verlangt, gefundene Fehler zu dokumentieren statt still zu
reparieren. Vier Stück, alle mit Regressionstest in
`tests/test_electrocapillary.cpp`:

**D1 — `force_beyond` war durch ganze Segmente quantisiert.** Die Kraft jenseits
einer Ausschlussdistanz zählte Segmente nach ihrem Mittelpunkt, ganz oder gar
nicht. Nahe der Kante trägt ein Segment mehrere Prozent der Gesamtkraft, die
Größe sprang also mit dem Netz aus einem Grund, der nichts mit dem Feld zu tun
hatte. Das Segment, das die Ausschlussdistanz überspannt, wird jetzt nach dem
Anteil seiner Bogenlänge geteilt. Regressionstest: `force_beyond` ist **stetig**
in der Ausschlussdistanz (gemessen: größter Sprung 1,6·10⁻³ der Gesamtkraft) und
fällt monoton.

**D2 — die projizierte Last war unstetig.** Die stückweise konstante
Segmentprojektion ist die konservative Größe, als rechte Seite aber
unstetig — und gegen eine unstetige rechte Seite kann kein Integrator
konvergieren. Jeder gekoppelte Punkt kam mit `AccuracyNotReached` zurück, und
zwar bei jeder Auflösung. Der Fehler lag in der Projektion, nicht im Integrator.
Ersetzt durch die stetige Rekonstruktion aus 10.4, die die Binintegrale exakt
erhält. Regressionstest: die übergebene Last trägt die integrierte Maxwell-Kraft
(1,5·10⁻⁴ relativ) und der Kapillarlöser erreicht seine geforderte Genauigkeit.

**D3 — das mechanische Residuum wurde punktweise an der Kante bewertet.** Das
widerspricht dem, was das Gate zuvor festgestellt hatte: dort konvergiert nichts
punktweise, also auch kein Residuum. Die Fortsetzung blieb bei 733 V stehen, aus
einer Eigenschaft der Lastdarstellung statt der Lösung. Die Schranke ist
unverändert; geändert hat sich, dass sie auf die Größe angewandt wird, für die
sie gedacht war — das Residuum außerhalb $d\ge0{,}1a$. Der Wert über die ganze
Oberfläche wird weiterhin gemessen und berichtet. Nach der Korrektur läuft der
Ast bis 1439,5 V.

**D4 — ein Gate-Kriterium maß nicht, was es messen sollte.** `kTolExclusion`
verlangte, dass das Halbieren der Ausschlussdistanz die integrierte Kraft um
weniger als 5 % ändert. Was diese Zahl vergleicht, ist der **Kraftinhalt zweier
verschiedener Gebiete** der Oberfläche; er ist eine physikalische Größe, wird
unter Verfeinerung nicht klein und darf es nicht. Mit dem gemessenen Exponenten
verhält sich die fehlende Kraft wie $d^{1+\beta}=d^{0{,}56}$, die Änderung pro
Halbierung kann also gar nicht unter etwa ein Drittel der vorigen fallen.

*Die Grenze wurde nicht gelockert.* Sie steht unverändert im Code, wird gemessen
(11,3 % für die ebene Form) und als **nicht eingehalten** berichtet. Die
Entscheidung trägt stattdessen ein Kriterium, das Summierbarkeit so prüft, wie
sie prüfbar ist, und das dieselbe schon vorhandene Toleranz benutzt: zwei
Extrapolationen desselben Grenzwertes aus Daten, die nichts teilen —

* über die **Netzstufen**: Aitken auf die Gesamtkraft der drei feinsten Stufen;
* über die **Ausschlussdistanz**: Ausgleichsgerade von $F(d_0)$ gegen
  $d_0^{1+\beta}$ mit dem unabhängig gefitteten $\beta$.

Wäre der Kantenbeitrag nicht summierbar, gäbe es die zweite gar nicht. Für die
ebene Form: 7,79·10⁻⁷ N gegen 7,60·10⁻⁷ N, Abweichung 2,5 %.

Zusätzlich zwei kleinere, ebenfalls dokumentierte Korrekturen: der
Ladungsvergleich ließ den Kontaktknoten zunächst mitlaufen und lag deshalb um
18 % daneben (jetzt 1,4·10⁻³); und eine Testschranke war absolut in Volt
formuliert für eine Größe, die proportional zur angelegten Spannung ist — sie
ist jetzt relativ.

---

## 10.7 Was P3b ausdrücklich nicht enthält

* **keine Emission** und kein emittierter Strom;
* **keine endliche Leitfähigkeit** der Flüssigkeit: sie ist ein idealer
  Äquipotentialleiter, der P2b/P2c-Vertrag, und es wird keine Ladungsrelaxation
  modelliert;
* keine Strömung, keine Viskosität, keine Speiseimpedanz;
* keine Raumladung — das Vakuum trägt keine freie Ladung;
* keine Zeitabhängigkeit und **keine dynamische Stabilität**;
* **kein Taylor-Kegel-Onset und kein Cone-Jet**;
* keine Schwerkraft (Bond-Zahl in 9.5);
* keine Verrundung der Austrittskante — ein physischer Kantenradius oder ein
  anderes Kontaktlinienmodell wäre eine neue Modellentscheidung;
* $\Delta p_\mathrm{exit}$ bleibt eine **Eingabe**; P3b bestimmt weder den
  Vorrats- noch den Speisedruck;
* die Stoffwerte für EMI-BF4 bleiben `illustrative` (siehe 9.4), solange
  $\gamma$ und $\rho$ nicht mit Temperatur und Quelle belegt sind. Jede
  physikalische Abbildung unterscheidet deshalb dimensionslose
  Solvervalidierung, illustratives Rechenbeispiel und Vorhersage — die dritte
  Kategorie ist leer.
