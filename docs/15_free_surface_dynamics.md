# 15 — P4: zeitabhängige freie Oberfläche

**Stand: 2026-08-29.** Status: **`infrastructure_only`**. Ein dynamischer
Meniskus ist **nicht** implementiert. Implementiert und gegen exakte Lösungen
geprüft ist genau ein Baustein: die **kinematische Randbedingung** auf einem
vorgeschriebenen Geschwindigkeitsfeld.

Code: `include/es/surface_kinematics.hpp`, `src/surface_kinematics.cpp`,
`apps/es_kinematics.cpp`, `tests/test_surface_kinematics.cpp`,
`python/plot_kinematics.py`.
Ergebnisse: `results/2026-08-29_p4_kinematics/`.

---

## 15.1 Die Vorbedingung — geprüft, und nicht erfüllt

Der Auftrag verlangt, zuerst zu prüfen, ob Punkt 3 eine physikalisch
geschlossene Strömungs- und Ladungstransportgrundlage liefert. **Er tut es
nicht**, und zwar aus drei strukturellen Gründen:

1. **Die Strömung von P3 hat keine freie Oberfläche.** Ihre Exaktheit beruht
   darauf, dass die Kontinuität im geraden Rohr $\partial u_z/\partial z=0$
   erzwingt und der konvektive Term deshalb identisch verschwindet. Eine sich
   verformende Oberfläche zerstört genau diese Voraussetzung. Es gibt also kein
   Geschwindigkeitsfeld, das man in eine kinematische Bedingung einsetzen
   könnte.
2. **Der Ladungstransport von P3 ist stationär und kennt keine Flächenladung.**
   Es fehlen die Oberflächenladungsbilanz mit Konvektion und — strukturell
   entscheidend — die **tangentiale Traktion** $q_sE_t$. Diese existiert im
   Perfect-Conductor-Grenzfall gar nicht und ist genau das, was einen endlich
   leitenden Meniskus antreibt.
3. **$\varepsilon_r$ ist `MissingMaterialData`** (P2). Damit ist nicht einmal
   die Ladungsrelaxationszeit berechenbar.

`solve_dynamic_meniscus()` wirft deshalb `NotImplementedInThisPhase` und nennt
alle drei Gründe. Der Test prüft, dass die Meldung sie wirklich nennt.

## 15.2 Drei Verbote, im Code als Abwesenheit durchgesetzt

* **Keine Mobilität.** Es gibt nirgends einen Koeffizienten, der eine Form mit
  einer gewählten Rate ins Gleichgewicht relaxiert. Ein solcher Koeffizient
  wäre ein freier Parameter im Kostüm der Physik.
* **Keine künstliche Dämpfung.** Es wird nichts geglättet, gefiltert oder
  relaxiert. Die tangentiale Umverteilung schiebt Knoten **entlang** der Fläche
  und kann die Fläche deshalb nicht ändern; wie genau das auf einem Polygonzug
  stimmt, wird gemessen (15.5).
* **Keine Stabilitätsaussage.** Das Ende eines statischen Astes ist keine
  dynamische Instabilität, und nichts hier macht eine daraus.

## 15.3 Was implementiert und validiert ist

$$
\frac{\mathrm d\mathbf x}{\mathrm dt}\cdot\mathbf n = \mathbf u\cdot\mathbf n
\qquad\text{auf der freien Oberfläche.}
$$

**Nur die Normalkomponente ist Physik.** Die Tangentialbewegung eines Knotens
ist eine Eigenschaft der Parametrisierung: ein Knoten, der entlang der Fläche
gleitet, lässt die Fläche unverändert. Zwei Modi, ausdrücklich **nicht**
austauschbar:

| Modus | $\mathrm d\mathbf x/\mathrm dt$ | wofür |
|---|---|---|
| `Lagrangian` | $\mathbf u$ | Knoten sind materielle Punkte; nur hier gibt es eine exakte Lösung, weil die Knotenbahnen die Charakteristiken sind |
| `NormalOnly` | $(\mathbf u\cdot\mathbf n)\mathbf n$ + Umverteilung | dieselbe Fläche, andere Parametrisierung; die Umverteilung ist **Netzbewegung** |

Zeitintegration: klassisches RK4 auf den Knotenbahnen.

### Die zwei analytisch kontrollierbaren Fälle

Beide sind vorgeschriebene Felder, deren **Lagrange-Abbildung in geschlossener
Form** bekannt ist — nicht nur das Volumen, sondern jede Knotenposition.

| Feld | $\nabla\cdot\mathbf u$ | exakte Abbildung | was es prüft |
|---|---|---|---|
| $\mathbf u=\alpha\mathbf x$ | $3\alpha$ | $\mathbf x\mapsto\mathbf x e^{\alpha t}$ | das Volumen muss einer **bekannten, von null verschiedenen** Änderung folgen |
| $\mathbf u=(-\alpha r/2,\ \alpha z)$ | $0$ **exakt** | $(r,z)\mapsto(re^{-\alpha t/2},\ ze^{\alpha t})$ | das Volumen muss **exakt erhalten** bleiben, während sich die Form stark ändert |

Zusammen trennen sie „das Volumen stimmt" von „die Form stimmt"; ein einzelner
Fall kann das nicht.

**Gemessen** (81 Knoten, $\alpha T=0{,}3$ bzw. $0{,}5$):

| Prüfung | Ergebnis |
|---|---|
| Nullfeld bewegt nichts | **exakt null**, bitgenau |
| Dilatation, Volumen gegen $V_0e^{3\alpha T}$ | $1{,}1\cdot10^{-15}$ |
| Dilatation, jeder Knoten gegen die exakte Abbildung | $3{,}9\cdot10^{-15}$ von $R$ |
| beobachtete Zeitordnung (80 → 160 Schritte) | **4,0** (RK4) |
| Squeeze, Volumenerhaltung | $2{,}1\cdot10^{-15}$ |
| Squeeze, Form gegen die exakte Abbildung | $6{,}8\cdot10^{-15}$ von $R$ |
| gepinnte Kontaktlinie | **bitgenau** unverändert |
| Apex auf der Achse | **exakt** $r=0$ |

Bei 640 Schritten liegt der Formfehler bei einigen $10^{-15}$ von $R$, also bei
etwa $10^{-20}$ m — das ist der Rundungsboden, nicht die Integration. Die
Ordnung wird deshalb auf den **groben** Stufen abgelesen und der Boden
angegeben, statt eine bedeutungslose Ordnung vom feinsten Paar zu nehmen.

## 15.4 Der Zeitschrittvertrag

Ein Schritt ist zulässig, wenn die entstehende Kurve noch eine Fläche ist.
Geprüft **nach** jedem Schritt, und ein verletzender Schritt wird **verworfen**,
nie mit einer Warnung angenommen:

1. kein Knoten überquert die Achse ($r<0$) → `NodeCrossedAxis`;
2. kein Segment kollabiert oder kehrt sich um → `SegmentCollapsed`;
3. kein Knoten bewegt sich weiter als `kMaxNodeMotion` = 0,25 der **kürzesten**
   Segmentlänge → `StepTooLarge`.

Die dritte ist eine **Diskretisierungs**schranke, keine physikalische: hier wird
keine Kraft integriert, also gilt keine kapillare oder akustische
Zeitschrittgrenze. Diese gehörten zu dem dynamischen Löser, den es nicht gibt.

Bemerkenswert und gemessen: beim Squeeze **verschärft** sich die Schranke im
Lauf, weil die Segmente am Äquator schrumpfen. 200 Schritte bleiben bei Schritt
180 mit `StepTooLarge` stehen — das ist der Vertrag bei der Arbeit, nicht ein
Fehler.

## 15.5 Eine gemessene Grenze der Netzbewegung

Im Kontinuum kann eine tangentiale Umverteilung die Fläche nicht ändern. Auf
einem **Polygonzug** schneidet das Neuabtasten Ecken ab. Gemessen wird der
Volumenunterschied zwischen `Lagrangian` und `NormalOnly` nach 900 Schritten:

| Knoten | 41 | 81 | 161 | 321 |
|---|---|---|---|---|
| $\lvert\Delta V\rvert/V$ | $4{,}1\cdot10^{-3}$ | $2{,}0\cdot10^{-3}$ | $9{,}1\cdot10^{-4}$ | $3{,}5\cdot10^{-4}$ |

**Beobachtete Ordnung 1,18 — nicht 2.** Der Verlust pro Umverteilung ist
$O(h^2)$, aber er akkumuliert über die Schritte, und die Akkumulation ist nicht
$h$-unabhängig. Das ist eine echte Grenze dieser Netzbewegung und steht hier,
statt weggeschärft zu werden: **ein dynamischer Löser, der in jedem Schritt
umverteilt, braucht eine gekrümmte Rekonstruktion oder muss den Verlust
bilanzieren.**

## 15.6 Der Vertrag des fehlenden Lösers

Damit die Lücke spezifiziert ist statt nur zugegeben:

### Zustandsgrößen

* die freie Oberfläche $\Gamma(t)$ als Polygonzug (`SurfacePolyline`), Apex auf
  der Achse, Kontaktlinie gepinnt oder mit einem Kontaktlinienmodell;
* das Geschwindigkeitsfeld $\mathbf u(\mathbf x,t)$ im **Volumen** der
  Flüssigkeit;
* der Druck $p(\mathbf x,t)$;
* die Flächenladungsdichte $q_s(s,t)$ auf $\Gamma$;
* das Potential $\varphi(\mathbf x,t)$ in Flüssigkeit, Dielektrikum und Vakuum.

### Randbedingungen auf $\Gamma$

* **kinematisch** — implementiert und geprüft:
  $\mathrm d\mathbf x/\mathrm dt\cdot\mathbf n=\mathbf u\cdot\mathbf n$;
* **Flächenladungsbilanz:**
  $$\frac{\partial q_s}{\partial t}+\nabla_s\cdot(q_s\mathbf u_s)
  =\sigma E_n^\text{innen}-j_\text{Emission};$$
  ohne Emission verlangt der stationäre Zustand $\sigma E_n^\text{innen}=0$ —
  genau der in P3 geprüfte Fall;
* **Normalsprung:**
  $\varepsilon_0E_n^\text{aussen}-\varepsilon_0\varepsilon_rE_n^\text{innen}=q_s$;
* **Normalspannung:**
  $\gamma\kappa = \Delta p + [\![T_{nn}]\!]$ mit dem **vollständigen**
  Maxwell-Spannungstensor beider Seiten, nicht nur $\varepsilon_0E_n^2/2$;
* **Tangentialspannung:**
  $\mu\bigl(\partial u_t/\partial n+\partial u_n/\partial t\bigr)=q_sE_t$.
  **Dieser Term fehlt im Perfect-Conductor-Grenzfall vollständig.** Ein endlich
  leitender Meniskus ist deshalb keine Korrektur des P3b-Ergebnisses, sondern
  ein anderes Problem.

### Zeitschrittvertrag des dynamischen Lösers

Zusätzlich zu 15.4 kämen die Schranken, die eine Kraftbilanz mitbringt:

* kapillar: $\Delta t \lesssim \sqrt{\rho\,h^3/(2\pi\gamma)}$ (Brackbill-Kothe-Zemach);
* viskos: $\Delta t\lesssim \rho h^2/(2\mu)$;
* elektrisch: $\Delta t\ll\tau=\varepsilon/\sigma$, sonst ist die
  Ladungsrelaxation nicht aufgelöst;
* Advektion: die CFL-Bedingung aus 15.4.

Der bindende ist der kleinste, und mit $h\sim0{,}1\,a$, $a=5\ \mu$m wäre
$\tau\approx10^{-10}$ s der bindende — die Rechnung wäre also
ladungsrelaxationsbegrenzt, was allein schon eine Aussage über die nötige
Rechenzeit ist.

### Diskretisierung

* **Volumen:** ein Löser für Stokes oder Navier-Stokes mit freier Oberfläche auf
  einem beweglichen achsensymmetrischen Netz. Q1/Q1 für Geschwindigkeit und
  Druck ist **instabil** (LBB); es braucht Taylor-Hood (Q2/Q1) oder eine
  Stabilisierung. Das vorhandene `axisym_fem` ist ein **Skalar**-Q1-Löser und
  trägt das nicht.
* **Oberfläche:** die Krümmung braucht mehr als einen Polygonzug — eine
  Spline- oder Bogenlängenrekonstruktion wie in P3a, sonst ist $\kappa$ nicht
  einmal stetig.
* **Netzbewegung:** siehe 15.5.

### Stabilitäts- und Energieprüfungen, die ein solcher Löser bestehen müsste

* **Volumenerhaltung** ohne Emission, bis auf die Zeitordnung;
* **Energiebilanz:** $\mathrm d/\mathrm dt(E_\text{kin}+\gamma A+E_\text{el})
  = -\Phi_\text{visc} - \Phi_\text{Joule} + P_\text{Zulauf}$, wobei beide
  Dissipationsterme **nicht negativ** sein müssen;
* **Ladungserhaltung:** Flächenladung plus Volumenladung plus emittierte Ladung
  konstant;
* **linearisierte Schwingung** einer gepinnten Kappe gegen eine Referenz —
  aber es gibt für die gepinnte achsensymmetrische Kappe keine geschlossene
  Eigenfrequenz, so dass die Referenz selbst eine Rechnung wäre. Das ist ein
  bekannter Validierungsengpass und kein Detail.

## 15.7 Was P4 ausdrücklich nicht enthält

Kein dynamischer Meniskus. Keine Kraftbilanz. Kein Strömungslöser mit freier
Oberfläche. Kein Oberflächenladungstransport. Keine Emission. Keine
Stabilitätsaussage, keine Taylor-Kegel-Aussage. Keine Mobilität, keine
künstliche Dämpfung. Das Geschwindigkeitsfeld ist **vorgeschrieben**; nichts
wird dafür gelöst.
