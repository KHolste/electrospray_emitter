# 14 — P3: endliche Leitfähigkeit und Flüssigkeitsströmung

**Stand: 2026-08-29.** Status: **`validated_subset`**. Zwei Teilprobleme sind
implementiert und gegen geschlossene Formen geprüft; die gekoppelte
finite-conductivity-Meniskusrechnung ist **nicht** implementiert und ihr
Vertrag steht in 14.5 als Spezifikation.

Code: `include/es/transport.hpp`, `src/transport.cpp`, `apps/es_transport.cpp`,
`tests/test_transport.cpp`, `python/plot_transport.py`.
Ergebnisse: `results/2026-08-29_p3_transport/`.

---

## 14.1 Ein Operator, drei Probleme

`AxisymProblem` assembliert

$$
-\nabla\cdot\bigl(\text{scale}\cdot k\,\nabla\varphi\bigr) = s .
$$

Drei Probleme dieses Projekts sind **genau** diese Gleichung:

| Problem | scale | $k$ | $s$ | $\varphi$ | Knotenreaktion |
|---|---|---|---|---|---|
| Elektrostatik | $\varepsilon_0$ | $\varepsilon_r$ | 0 | Potential [V] | Ladung [C] |
| Leitung | 1 | $\sigma$ [S/m] | 0 | Potential [V] | **Strom [A]** |
| ausgebildete Rohrströmung | 1 | $\mu$ [Pa s] | $-\mathrm dp/\mathrm dz$ | $u_z$ [m/s] | — |

An der Assemblierung ändert sich **nichts**; nur die Bedeutung der Zahlen. Neu
sind `coefficient_scale` (Vorgabe $\varepsilon_0$, also unverändert für jeden
bestehenden Aufrufer) und `cell_source` (leer = null). Die Prüfung
$\varepsilon_r\ge1$ gilt weiterhin, aber nur für den elektrostatischen Fall —
eine Leitfähigkeit von 1,5 S/m oder eine Viskosität von 0,036 Pa·s sind völlig
gewöhnlich.

## 14.2 Warum die Strömungsreduktion exakt ist, und wo sie aufhört

Für stationäre inkompressible Strömung im **geraden** Kreisrohr mit
$\mathbf u = u_z(r,z)\,\mathbf e_z$ erzwingt die Kontinuität sofort
$\partial u_z/\partial z=0$, also $u_z=u_z(r)$. Damit verschwindet der
konvektive Term $(\mathbf u\cdot\nabla)\mathbf u$ **identisch** — nicht wegen
kleiner Reynoldszahl, sondern weil $\mathbf u$ senkrecht auf seinem eigenen
Gradienten steht. Die $z$-Impulsgleichung ist dann

$$
0 = -\frac{\mathrm dp}{\mathrm dz} + \mu\,\frac1r\frac{\partial}{\partial r}
\Bigl(r\frac{\partial u_z}{\partial r}\Bigr),
$$

und die rechte Seite **ist** der achsensymmetrische Laplace-Operator von $u_z$.

**Das ist kein Stokes-Löser.** Es kann keine Einlaufströmung, keine Verengung,
keine Krümmung und keine freie Oberfläche — in all diesen Fällen ist
$\partial u_z/\partial z\ne0$ und der Druck keine reine Funktion von $z$.

**Was es leistet:** eine unabhängige Prüfung des hydraulischen Widerstands aus
P1. P1 schreibt $8\mu L Q/(\pi R^4)$ hin; hier wird das Feld gelöst und das
Profil integriert. Die beiden teilen **keinen Code**.

Gemessen: Volumenstrom gegen Hagen-Poiseuille **2,6·10⁻⁵ relativ** bei 81
radialen Knoten, Netzkonvergenz **zweiter Ordnung** (beobachtet 2,00), und der
so gewonnene Widerstand trifft die Formel von P1 auf dieselben 2,6·10⁻⁵.
Zusätzlich exakt geprüft: Linearität in $\mathrm dp/\mathrm dz$, Umkehrung in
$\mu$, Haftbedingung, und dass ohne Druckgradient die Geschwindigkeit **exakt**
null ist.

Randbedingungen: Haftbedingung an der Wand; Achse und beide Stirnflächen tragen
die **natürliche** Bedingung $\partial u/\partial n=0$ — die die exakte Lösung
erfüllt, so dass nichts erzwungen wird, was sie nicht ohnehin tut.

## 14.3 Der Ladungstransportvertrag

$$
\mathbf j = \sigma\mathbf E = -\sigma\nabla\varphi, \qquad
\frac{\partial\rho}{\partial t} + \nabla\cdot\mathbf j = 0, \qquad
\nabla\cdot(\varepsilon\nabla\varphi) = -\rho .
$$

Für homogene $\sigma,\varepsilon$ folgt

$$
\frac{\partial\rho}{\partial t} + \frac{\sigma}{\varepsilon}\rho = 0
\;\Longrightarrow\;
\rho(t)=\rho(0)\,e^{-t/\tau}, \qquad
\tau = \frac{\varepsilon}{\sigma} = \frac{\varepsilon_0\varepsilon_r}{\sigma}.
$$

Freie Ladung kann im **Inneren** eines homogenen Leiters nicht leben; sie zerfällt
mit $\tau$ und landet auf dem Rand. Das ist der ganze Inhalt des
Perfect-Conductor-Grenzfalls.

**Stationär** wird die Ladungserhaltung zu $\nabla\cdot(\sigma\nabla\varphi)=0$
— dieselbe Gleichung wie die Elektrostatik. Geprüft im Zylinder gegen
$I=\sigma A V/L$ und $R=L/(\sigma A)$: **3,7·10⁻¹¹ relativ**, Potentialfehler
3,9·10⁻¹² (die Q1-Lösung einer linearen exakten Lösung *ist* diese Lösung; der
Rest ist die Rundung des direkten Lösers). Ohmsches Gesetz Term für Term
geprüft: linear in $V$, in $\sigma$, in $A$, umgekehrt in $L$, unabhängig von
der Polarität, und bei $V=0$ **exakt** kein Strom.

### Ohne Emission kein Normalstrom durch die freie Oberfläche

Das ist eine physikalische Aussage, keine numerische Feinheit: Ladung, die die
freie Oberfläche überquert, hätte nirgendwohin zu gehen, also verbietet die
Erhaltung sie. Die richtige Randbedingung ist $\mathbf j\cdot\mathbf n=0$ — die
**natürliche** Bedingung desselben Operators. Ein Löser, der dort stattdessen
ein Potential vorschriebe, triebe einen Strom aus der Flüssigkeit heraus und
nennte das Physik.

Gemessen an der Mantelfläche des Zylinders, die genau diese Bedingung trägt:
größte radiale Stromdichte relativ zur mittleren axialen **< 10⁻¹²**.

## 14.4 Der Grenzfall ist ein Verhältnis — und welche Permittivität darin steht

$\tau\ll t_\text{Prozess}$. **Welche** Prozesszeit gilt, ist eine
Modellentscheidung: für eine statische Form die Einstellzeit der Form, für
emittierenden Betrieb die Transitzeit durch die Emissionszone. Der Code nimmt
$t_\text{Prozess}$ deshalb als **explizite Eingabe** und berichtet das
Verhältnis; er wählt keine.

### 14.4.1 Eine frühere Fassung dieses Abschnitts hat zweimal falsch gelegen

Sie verlangte einen einzelnen „DC-Wert“ von $\varepsilon_r$ und verwarf die
vorhandene 1–18-GHz-Messung mit der Begründung, sie sei „nicht DC“. Beides war
falsch, und zwar in entgegengesetzte Richtungen: die Forderung nach einem
Gleichstromwert war unbegründet, und die Ablehnung der GHz-Daten verwarf genau
die Messung, die die Formel verlangt.

### 14.4.2 Fünf Größen heißen „die Permittivität“

Bei einer **leitfähigen** ionischen Flüssigkeit sind das keine Varianten einer
Zahl, sondern verschiedene Dinge:

| | Größe | was sie ist | im Datensatz |
|---|---|---|---|
| (1) | **statische bzw. niederfrequente Scheinpermittivität** | was eine Kapazitätsmessbrücke bei kHz anzeigt. Bei $K\sim1$ S/m vom Leitungsbeitrag und von der Elektrodenpolarisation beherrscht, Werte bis $10^4$ und darüber. Eine Eigenschaft der **Messzelle** | **kein einziger Punkt** — geprüft, nicht behauptet |
| (2) | **intrinsische statische Permittivität** $\varepsilon_s$ | Grenzwert der *dielektrischen* Dispersion für $f\to0$, nachdem der Leitungsbeitrag abgetrennt wurde. Wird bei ionischen Flüssigkeiten nicht bei null Hertz gemessen, sondern aus Mikrowellenspektren extrapoliert | 3 Punkte: 12,8 / 12,9 / 13,6 |
| (3) | **frequenzabhängige komplexe Permittivität** $\varepsilon^*(f)$ | an einzelnen Frequenzen gemessen | 18 Punkte, 1–18 GHz, mit Unsicherheiten |
| (4) | **Elektrodenpolarisation** | Messartefakt: Ionen schirmen an den Elektroden das Feld ab und blähen die scheinbare Kapazität auf. Der Grund, warum für diese Flüssigkeiten überhaupt Mikrowellenspektroskopie benutzt wird | — |
| (5) | **DC-Leitfähigkeit** $K$ | eigene Größe mit eigener Quelle. Sie als imaginäre Permittivität $K/(\varepsilon_0\omega)$ zu schreiben ist Buchhaltung, keine zweite Dielektrizitätszahl | 1,5584 S/m, `measured` |

Die Schwelle, unterhalb derer ein Punkt als (1) gilt, ist in
`es::transport::kElectrodePolarisationFloor` **vorab** auf 1 MHz festgelegt.
`permittivity_points.csv` führt jeden Punkt mit seiner Einordnung, und der Test
prüft, dass der Datensatz keinen Punkt unterhalb der Schwelle enthält — die
Aussage „hier steckt keine kHz-Scheinpermittivität drin“ ist damit nachgeprüft
und nicht angenommen.

### 14.4.3 Welche davon in $\tau_q$ gehört — hergeleitet, nicht gewählt

Aus Ladungserhaltung und Gauß folgt

$$
\frac{\partial\rho_f}{\partial t}+\frac{K}{\varepsilon_0\varepsilon_r}\rho_f=0 ,
$$

wobei $\varepsilon_r$ die **gebundene** Ladungsantwort ist, also die
Polarisation, die dem Feld folgt, *während* die freie Ladung zerfällt. Die
freie Ladung zerfällt auf der Zeitskala $\tau$ selbst; ihr Spektrum liegt damit
bei $\omega\sim1/\tau$, also bei

$$
f^{*}=\frac{1}{2\pi\tau}.
$$

Da $\varepsilon_r$ dispersiv ist, ist die Gleichung für $\tau$ **implizit**:

$$
\boxed{\;\tau=\frac{\varepsilon_0\,\varepsilon_r(f^{*})}{K},\qquad
f^{*}=\frac{1}{2\pi\tau}\;}
$$

Für diese Flüssigkeit liegt $f^{*}$ bei einigen GHz — **genau dort, wo Bennett
et al. gemessen haben**. Die 1–18-GHz-Daten sind also nicht die falsche
Frequenz, sondern die richtige. Sie umgekehrt als Gleichstromwert zu übernehmen
wäre ebenso falsch; die Kurve wird an der Stelle ausgewertet, die die Physik
auswählt.

### 14.4.4 Die Lösung auf der gemessenen Kurve

`self_consistent_relaxation()` löst die implizite Gleichung durch
Fixpunktiteration auf den **Messpunkten**, logarithmisch zwischen benachbarten
Messfrequenzen interpoliert. Es wird **keine Dispersionsfunktion angepasst** —
eine anzupassen hieße, Daten zu erfinden — und außerhalb des Messbereichs wird
nicht extrapoliert; ob $f^{*}$ innerhalb liegt, wird berichtet.

Bei 298,15 K:

| | |
|---|---|
| $f^{*}$ | $2{,}627\;\mathrm{GHz}$ — **innerhalb** des gemessenen Bereichs 1–18 GHz |
| $\varepsilon_r(f^{*})$ | $10{,}664$ |
| $K$ | $1{,}5584\;\mathrm{S/m}$ (ausgewählte Quelle) |
| $\tau$ | $6{,}059\cdot10^{-11}\;\mathrm{s}$ |
| Iterationen | 18, Selbstkonsistenz-Residuum $0$ |
| zum Vergleich mit $\varepsilon_s=13{,}6$ | $\tau=7{,}727\cdot10^{-11}$ s, also **+27,5 %** |

Die letzte Zeile ist der ganze Unterschied zwischen „statisch“ und „bei
$f^{*}$“: knapp 28 %. Er ist berichtenswert und für das Urteil unten belanglos.

Die beiden definierenden Identitäten $\tau=\varepsilon_0\varepsilon_r/K$ und
$f^{*}=1/(2\pi\tau)$ gelten in der Rückgabe **exakt**; das Residuum misst die
eine Identität, die ein Fixpunkt nur bis auf seine Konvergenz erfüllen kann,
nämlich $\varepsilon_r=\varepsilon_r(f^{*})$.

### 14.4.5 Ein Einzelwert bleibt fehlend — ein Band ist belegt

Keine der vier Permittivitätsquellen nennt Reinheit **und** Wassergehalt. Die
Auswahlregel von P2 wählt deshalb weiterhin **keine** aus, und
`material_value(RelativePermittivity)` meldet unverändert
`MissingMaterialData`. `judge_conductor_limit()` — die Einzelwertabfrage —
schlägt weiterhin geschlossen fehl. Der Lauf prüft das ausdrücklich und setzt
sonst `exit_code = 2`.

Was **statt** eines Ersatzwertes belegt werden kann, ist ein Band und die
Empfindlichkeit darüber:

$$
\varepsilon_r\in[7{,}7,\;13{,}6],\qquad K\in[0{,}91,\;1{,}63]\;\mathrm{S/m}
$$

Das $\varepsilon_r$-Band umfasst alle als dielektrisch zulässigen Punkte, also
die drei extrapolierten statischen Werte **und** die ganze gemessene
Dispersion. Es wird nicht gemittelt.

### 14.4.6 Das Urteil, an der ungünstigsten Ecke gefällt

$\tau$ wächst mit $\varepsilon_r$ und fällt mit $K$; die für die Näherung
schlechteste Ecke ist also $(\varepsilon_\text{hi},K_\text{lo})$. Das wird
**gerechnet und geprüft**, nicht behauptet — ein Vorzeichenfehler in dieser
Überlegung wäre sonst unsichtbar.

| Ecke | $\tau$ [s] | $t_\text{kap}/\tau$ | $t_\text{vis}/\tau$ |
|---|---|---|---|
| $\varepsilon_\text{lo},K_\text{hi}$ | $4{,}183\cdot10^{-11}$ | $4{,}12\cdot10^{4}$ | $8{,}05\cdot10^{4}$ |
| $\varepsilon_\text{lo},K_\text{lo}$ | $7{,}492\cdot10^{-11}$ | $2{,}30\cdot10^{4}$ | $4{,}49\cdot10^{4}$ |
| $\varepsilon_\text{hi},K_\text{hi}$ | $7{,}388\cdot10^{-11}$ | $2{,}33\cdot10^{4}$ | $4{,}56\cdot10^{4}$ |
| **$\varepsilon_\text{hi},K_\text{lo}$** (schlechteste) | $1{,}323\cdot10^{-10}$ | $\mathbf{1{,}30\cdot10^{4}}$ | $\mathbf{2{,}54\cdot10^{4}}$ |
| selbstkonsistent | $6{,}059\cdot10^{-11}$ | $2{,}84\cdot10^{4}$ | $5{,}56\cdot10^{4}$ |

mit $t_\text{kap}=\sqrt{\rho a^3/\gamma}=1{,}722\cdot10^{-6}$ s und
$t_\text{vis}=\mu a/\gamma=3{,}367\cdot10^{-6}$ s aus den **belegten**
Stoffdaten. Die geforderte Schranke ist 100.

**Damit ändert sich der Befund von P3.** Der Äquipotentialansatz von P3b ist
für die **statische** Form nicht mehr eine Annahme, sondern über das gesamte
begründete Band belegt — und zwar an der ungünstigsten Ecke noch um mehr als
zwei Größenordnungen über der Schranke. Es wird dafür **kein einziger
unbelegter $\varepsilon_r$-Wert** benutzt: das Ergebnis hängt nicht davon ab,
welchen Wert im Band man nähme, weil *jeder* Wert im Band dasselbe Urteil
liefert. Genau das ist der Sinn einer Empfindlichkeitsrechnung.

**Was das NICHT sagt.** Es sagt nichts über emittierenden Betrieb: dort ist die
Prozesszeit die Transitzeit durch die Emissionszone, und die ist hier nicht
gerechnet. Es sagt nichts über die tangentiale Traktion $q_sE_t$, die im
Perfect-Conductor-Grenzfall gar nicht existiert (14.5). Und es macht
$\varepsilon_r$ nicht zu einer belegten Zahl — für jede Rechnung, die einen
*Wert* statt eines *Verhältnisses* braucht, fehlt er weiterhin.

## 14.5 Was NICHT implementiert ist — der Vertrag, damit es später geht

Die gekoppelte finite-conductivity-Meniskusrechnung fehlt. Ihre
Randbedingungen sind hier festgehalten, damit klar ist, was fehlt:

Auf der freien Oberfläche mit Flächenladungsdichte $q_s$ und
Oberflächengeschwindigkeit $\mathbf u_s$:

* **Ladungserhaltung auf der Fläche**
  $$\frac{\partial q_s}{\partial t} + \nabla_s\cdot(q_s\mathbf u_s)
  = \sigma E_n^\text{innen} - j_\text{Emission},$$
  d. h. der Leitungsstrom aus dem Inneren speist die Flächenladung, die
  Oberflächenkonvektion transportiert sie, und die Emission entfernt sie. Ohne
  Emission muss der stationäre Zustand $\sigma E_n^\text{innen}=0$ erfüllen —
  genau der in 14.3 geprüfte Fall.
* **Sprung der Normalkomponente**
  $\varepsilon_0E_n^\text{aussen}-\varepsilon_0\varepsilon_r E_n^\text{innen}=q_s$.
* **Tangentiale Spannung.** Bei endlicher Leitfähigkeit ist $E_t\ne0$ auf der
  Oberfläche, und $q_sE_t$ ist eine **tangentiale** Traktion, die die
  Flüssigkeit in Bewegung setzt. Diese Traktion existiert im
  Perfect-Conductor-Grenzfall **gar nicht** — sie ist genau das, was P3b
  strukturell fehlt, und der Grund, weshalb ein endlich leitender Meniskus keine
  bloße Korrektur des P3b-Ergebnisses ist.
* **Normale Spannungsbilanz** wie P3b, aber mit dem vollständigen
  Maxwell-Spannungstensor beider Seiten statt nur $\varepsilon_0E_n^2/2$.

**Was dafür fehlt:** ein Löser für das Geschwindigkeitsfeld im Inneren mit
freier Oberfläche (die Reduktion aus 14.2 trägt das nicht), ein
Oberflächenladungstransport mit einer bewegten Fläche, und $\varepsilon_r$ als
belegte Zahl. Solange das fehlt, ist hier `infrastructure_only`, und der Code
enthält **keine** Näherung, die so täte, als wäre es gelöst.

## 14.6 Was P3 ausdrücklich nicht enthält

Kein allgemeiner Stokes-Löser (keine Einlaufströmung, keine
Druck-Geschwindigkeits-Kopplung, keine gekrümmte Geometrie, keine freie
Oberfläche). Keine gekoppelte finite-conductivity-Meniskusrechnung. Keine
Emission und kein emittierter Strom. Keine Zeitintegration einer Form —
integriert wird nur die geschlossene Relaxations-ODE, die eine Lösung hat.
Keine Raumladung. Keine Stabilitätsaussage.
