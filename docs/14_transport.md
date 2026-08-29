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

## 14.4 Der Grenzfall ist ein Verhältnis — und mit unseren Daten unbekannt

$\tau\ll t_\text{Prozess}$. **Welche** Prozesszeit gilt, ist eine
Modellentscheidung: für eine statische Form die Einstellzeit der Form, für
emittierenden Betrieb die Transitzeit durch die Emissionszone. Der Code nimmt
$t_\text{Prozess}$ deshalb als **explizite Eingabe** und berichtet das
Verhältnis; er wählt keine.

**Und hier steht der Befund.** $\tau=\varepsilon_0\varepsilon_r/\sigma$ braucht
$\varepsilon_r$, und $\varepsilon_r$ ist für EMI-BF4 nach dem Stoffdatenvertrag
von P2 **`MissingMaterialData`**: keine der vier Quellen nennt Reinheit und
Wassergehalt, und die einzige mehrpunktige misst von 1 bis 18 GHz. Also ist
$\tau$ **mit den belegten Daten dieses Projekts nicht berechenbar**, und
`judge_conductor_limit()` liefert `MissingMaterialData` statt einer Zahl.

Damit ist der Äquipotentialansatz von P3b **eine Annahme und kein
nachgewiesener Grenzfall**. Mit einem ausdrücklich unbelegten
$\varepsilon_r=12{,}8$ ergäbe sich $\tau\approx7{,}3\cdot10^{-11}$ s gegen
$t_\text{kap}=\sqrt{\rho a^3/\gamma}\approx1{,}7\cdot10^{-6}$ s und
$t_\text{vis}=\mu a/\gamma\approx3{,}4\cdot10^{-6}$ s, also Verhältnisse von
$2\cdot10^4$ bis $5\cdot10^4$ — die Annahme ist plausibel, aber sie ist mit
diesem Datenstand nicht belegt, und das steht in jeder Ausgabe.

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
