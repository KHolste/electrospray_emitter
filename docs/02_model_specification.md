# 2. Modellspezifikation

Jede Gleichung trägt Annahmen, Gültigkeitsbereich und Primärquelle. Der
Prüfstatus der Quellen steht in [references.md](references.md); in diesem
Session-Lauf bibliographisch verifizierte Stellen sind dort mit ✓ markiert,
alle übrigen sind aus dem Gedächtnis wiedergegeben und **vor Zitation in einer
Publikation zu prüfen**.

## 2.0 Gebiete, Symmetrie, Notation

Axialsymmetrische Meridianhalbebene, Koordinaten $(r,z)$, $r \ge 0$. Die Achse
$r=0$ ist Symmetrierand, kein physikalischer Rand.

| Gebiet | Bedeutung |
|---|---|
| $\Omega_\mathrm{v}$ | Vakuum zwischen Emitter, Extraktor und Gehäuse |
| $\Omega_\ell$ | Flüssigkeit (Speisekanal und Meniskus) |
| $\Sigma$ | freie Oberfläche des Meniskus |
| $\Gamma_k$ | Elektrodenoberfläche $k$ (Emitter, Extraktor, Gehäuse, Kollektor) |
| $L$ | Kontaktlinie, $\Sigma \cap \Gamma_\mathrm{emitter}$ |

Feldgrößen: Potential $\phi$, elektrisches Feld $\mathbf{E}=-\nabla\phi$,
Geschwindigkeit $\mathbf{u}$, Druck $p$, freie Raumladungsdichte $\rho_f$,
Flächenladungsdichte $\sigma_s$ auf $\Sigma$. Normale $\hat{n}$ zeigt aus der
Flüssigkeit in das Vakuum. Hochindizes $^\mathrm{v}$ und $^\ell$ bezeichnen die
Vakuum- bzw. Flüssigkeitsseite von $\Sigma$.

---

## Stufe 1 — Elektrostatik vor der Emission

### 1.1 Feldgleichung

$$\nabla^2\phi = 0 \quad \text{in } \Omega_\mathrm{v}, \qquad
\phi = U_k \ \text{auf } \Gamma_k, \qquad
\phi \to 0 \ \text{für } |x|\to\infty$$

**Unbekannte:** $\phi$ in $\Omega_\mathrm{v}$. In der Randintegralformulierung
äquivalent die Flächenladungsdichte $\sigma$ auf $\bigcup_k\Gamma_k$.

**Randbedingungen:** Dirichlet auf allen Elektroden. Fernfeld entweder
Abklingbedingung (BEM, natürlich erfüllt) oder $\phi=U_\mathrm{Gehäuse}$ auf dem
Rechteckrand der Skizze (Volumenverfahren).

**Annahmen.** (a) Alle Elektroden sind ideale Leiter. (b) Keine Volumenladung.
(c) Stationär, keine Magnetfelder. (d) Der Meniskus ist Teil der
Emitterelektrode (nur zulässig unter der Bedingung in 2.2).

**Gültigkeitsbereich.** Vor Emissionsbeginn exakt. Im emittierenden Betrieb nur,
solange das Eigenfeld des Strahls klein gegen das angelegte Feld ist. Als
prüfbares Kriterium: für einen Strahl mit Strom $I$, Radius $R_b$ und
Geschwindigkeit $v$ ist das Raumladungsfeld am Strahlrand

$$E_\mathrm{sc} \approx \frac{I}{2\pi\varepsilon_0 R_b v},$$

und Stufe 1 ist brauchbar, solange $E_\mathrm{sc}/E_\mathrm{angelegt} \ll 1$.
Dieses Verhältnis ist bei jedem Lauf auszugeben, nicht anzunehmen.

**Quelle.** Lehrbuchstoff (Jackson). Für die axialsymmetrische
Randintegralformulierung: die einmal über den Azimut integrierte
Freiraum-Greensfunktion

$$G(r,z;r',z') = \frac{r'\,K(m)}{\pi\varepsilon_0 S},\qquad
S^2=(r+r')^2+(z-z')^2,\qquad m=\frac{4rr'}{S^2}$$

mit $K$ dem vollständigen elliptischen Integral erster Art (Parameter-
konvention). Für einen geschlossenen Leiter gilt an der Oberfläche exakt
$E_n = \sigma/\varepsilon_0$.

**Status.** Implementiert. Gegen Kugel- und Rotationsellipsoid-Analytik geprüft
(Kapazität 2,5·10⁻⁶ relativ, Spitzenfeld bei Aspektverhältnis 20 auf 0,18 %).
Dies ist der einzige Modellteil, für den derzeit ein unabhängiger Nachweis
vorliegt.

### 1.2 Was ausgegeben werden muss

Potential, Feld, Flächenladung, Feldüberhöhung $E_\mathrm{Spitze}/(U/d)$,
Kapazitätsmatrix der Elektrodenanordnung, sowie das Kriterium aus 1.1.

---

## Stufe 2 — Statischer Meniskus, nicht emittierend

### 2.1 Feldgleichung auf der freien Oberfläche

Normalspannungsbilanz mit Maxwell-Zug an einer perfekt leitenden Oberfläche:

$$\gamma\,\kappa(s) \;=\; \Delta p - \rho g\,z \;+\; \frac{\varepsilon_0}{2}E_n^2
\qquad\text{auf }\Sigma$$

mit der Summe der Hauptkrümmungen in Bogenlängenparametrisierung

$$\kappa = \frac{\mathrm{d}\varphi}{\mathrm{d}s} + \frac{\sin\varphi}{r},
\qquad \frac{\mathrm{d}r}{\mathrm{d}s}=\cos\varphi,\qquad
\frac{\mathrm{d}z}{\mathrm{d}s}=-\sin\varphi .$$

**Unbekannte:** die Kurve $\Sigma$, parametrisiert durch $(r(s),z(s),\varphi(s))$,
und — je nach Formulierung — einer der Werte $\Delta p$, $U$ oder die Apexhöhe
$h$. Genau zwei der drei sind vorzugeben.

**Randbedingungen.**

* Achse: $r(0)=0$, $\varphi(0)=0$, Regularität $\mathrm{d}\varphi/\mathrm{d}s|_0
  = \tfrac12\kappa(0)$ (beide Hauptkrümmungen gleich).
* Kontaktlinie, **Fall A (gepinnt, scharfe Kante)**: $r(s_\mathrm{end}) = r_c$.
  Gültig, solange der scheinbare Kontaktwinkel zwischen dem Young-Winkel
  $\theta_Y$ und $\theta_Y + \beta$ liegt, wobei $\beta$ der Kantenwinkel ist
  (Gibbs'sche Pinning-Bedingung / Canthotaxis). Diese Bedingung ist zu
  **prüfen**, nicht anzunehmen: verlässt der Winkel das Intervall, bewegt sich
  die Kontaktlinie und Fall A ist ungültig.
* Kontaktlinie, **Fall B (Kontaktwinkel)**: $\varphi(s_\mathrm{end})$ folgt aus
  dem vorgegebenen Kontaktwinkel $\theta_c$ zur Wandnormalen; $r_c$ ist dann
  unbekannt.
* Elektrisch: $\phi=U_\mathrm{emitter}$ auf $\Sigma$.

Der Prototyp implementiert ausschließlich Fall A und prüft die Pinning-Bedingung
nicht. Die Skizze zeigt eine Stirnfläche, deren Kantenausführung offen ist
(offene Frage 4).

**Annahmen.** (a) Statisch, $\mathbf{u}=0$. (b) Perfekter Leiter: $\Sigma$ ist
Äquipotentialfläche, kein Feld in der Flüssigkeit. (c) Keine Massen- oder
Ladungsabfuhr. (d) Keine viskosen Spannungen. (e) $\gamma$ konstant (keine
Marangoni-Effekte, keine Temperaturgradienten entlang $\Sigma$).

**Gültigkeitsbereich.** Annahme (b) ist über die Ladungsrelaxationszeit
$\tau_e = \varepsilon_0\varepsilon_r/K$ zu rechtfertigen — für EMI-BF4
$\approx 8\cdot10^{-11}$ s. Der zu vergleichende Prozesszeitmaßstab ist im
statischen Fall formal unendlich, weshalb (b) hier zulässig ist. **Im
emittierenden Betrieb ist das nicht mehr der Fall** (Stufe 3). Annahme (a)
verlangt $\mathrm{Ca}=\mu u/\gamma \ll 1$ und einen vernachlässigbaren
Speisefluss. Annahme mit Schwerkraft: Bond-Zahl
$\mathrm{Bo}=\rho g r_c^2/\gamma \approx 3\cdot10^{-5}$ bei $r_c=10\,\mu$m,
also entbehrlich, aber implementiert.

**Quellen.** Young-Laplace mit Maxwell-Zug: Taylor (1964); Übersichtsdarstellung
Fernández de la Mora (2007). Gibbs'sche Kontaktlinienfixierung: Standardstoff
der Benetzungstheorie. Zahlenwert des Taylor-Winkels 49,29° als Nullstelle von
$P_{1/2}(\cos\theta)$: Taylor (1964).

### 2.2 Wann darf der Meniskus als Äquipotentialfläche behandelt werden

Nicht pauschal. Das Kriterium ist nicht $\tau_e$ allein, sondern das Verhältnis
des ohmschen Nachschubstroms zur abgeführten Ladung. Solange keine Ladung
abgeführt wird (Stufe 2), fällt keine Spannung in der Flüssigkeit ab und die
Annahme ist exakt. Sobald ein Strom $I$ fließt, ist der Spannungsabfall in der
Flüssigkeit von der Größenordnung

$$\Delta\phi_\ell \sim \frac{I}{K \, \ell}$$

mit $\ell$ einer charakteristischen Länge des stromführenden Bereichs. Higuera
(2008) ✓ kommt für den rein ionischen Betrieb zu dem Ergebnis, dass genau
dieser Widerstand den Emissionsstrom **kontrolliert**. Das Perfect-Conductor-
Modell ist dort also nicht zulässig.

### 2.3 Numerische Formulierung

Die Kopplung ist: Form $\to$ Feld (Stufe 1 mit $\Sigma$ als Teilrand) $\to$
Maxwell-Zug $\to$ Form. Der Prototyp löst das als Fixpunktiteration mit
Unterrelaxation. Für die Neufassung ist stattdessen ein **Newton-Verfahren auf
dem gekoppelten Residuum** vorzusehen, mit der Formableitung des Randintegral-
operators (shape derivative). Begründung: die Fixpunktiteration konvergiert in
der Nähe des Umkehrpunkts nicht mehr (im Prototyp beobachtet: 34 bis über 60
Iterationen, dann Abbruch), und ohne Jacobi-Matrix gibt es auch kein
Stabilitätskriterium (2.4).

### 2.4 Stabilität — was der Umkehrpunkt bedeutet und was nicht

Der Prototyp bezeichnet das Maximum von $U(h)$ als „Onset". Das ist in dieser
Form nicht haltbar. Korrekt zu trennen sind:

**(a) Verlust der statischen Stabilität.** Für ein konservatives System mit
einem Kontrollparameter wechselt die Stabilität an Umkehrpunkten des
Lösungsastes (Poincaré'scher Austausch der Stabilität). Für elektrifizierte
hängende und liegende Tropfen ist das explizit durchgerechnet von Wohlhuter &
Basaran (1992, JFM 235, 481) ✓, die Formen *und* Stabilität berechnen. Zwei
Einschränkungen, die der Prototyp nicht berücksichtigt:

* Das Kriterium gilt für den betrachteten Störungsmodus. Axialsymmetrische
  Umkehrpunkte sagen nichts über nicht-axialsymmetrische Instabilitäten.
* Es hängt davon ab, **welcher** Parameter festgehalten wird. Ein Ast bei
  festem $\Delta p$ und variabler $U$ hat andere Stabilitätsgrenzen als einer
  bei fester Ladung oder festem Volumen. Der Prototyp hält $\Delta p$ fest und
  parametrisiert über $h$; ob das gelieferte Maximum von $U(h)$ die
  Stabilitätsgrenze bei fester Spannung ist, ist zu zeigen, nicht anzunehmen.

Der belastbare Weg ist die Auswertung des Vorzeichens der kleinsten
Eigenwerte der Jacobi-Matrix des Newton-Systems entlang des Astes. Das ist ein
weiterer Grund für 2.3.

**(b) Ein literaturbelegtes Kriterium für den emittierenden Fall.** Gallud &
Lozano (2022, JFM 933, A43) ✓ berechnen stationäre emittierende IL-Menisken und
geben an, dass die statische Stabilität verloren geht, wenn der elektrische
Druck am Feld unmittelbar an der haltenden Elektrode den doppelten
Oberflächenspannungsdruck einer Kugel gleicher Größe übersteigt. Sie finden
zudem, dass statisch stabile emittierende Lösungen an eine minimale hydraulische
Impedanz, einen maximalen Strom und ein enges Fensterr von Hintergrundfeldern
gebunden sind, mit einem Optimum bei Meniskusradien von 0,5–3 µm. Dieses
Kriterium ist zu implementieren und gegen die eigene Rechnung zu prüfen.

**(c) Emissionsbeginn ist keine Bifurkation.** Im PIR ist die
Verdampfungsrate eine glatte Exponentialfunktion des Feldes ohne Schwellwert.
Was experimentell „Onset" heißt, ist der Punkt, an dem der Strom die
Nachweisgrenze des Messaufbaus überschreitet — eine Geräteeigenschaft. Der Code
muss deshalb den Strom als Funktion der Spannung ausgeben und die
Nachweisgrenze als *Eingabeparameter* führen, statt eine Bifurkation als
„Onset" zu etikettieren.

**(d) Der Übergang in den Cone-Jet** ist eine dritte, davon verschiedene
Erscheinung, die die Strömung einschließt und mit einem statischen Modell
grundsätzlich nicht bestimmbar ist.

---

## Stufe 3 — Emittierender Betrieb

Hier ist das Perfect-Conductor-Modell aufzugeben. Grundlage ist das
Taylor-Melcher-Leaky-Dielectric-Modell.

### 3.1 Bulk-Gleichungen in $\Omega_\ell$

$$\nabla\cdot\mathbf{u}=0$$
$$\rho\,(\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p + \mu\nabla^2\mathbf{u} + \rho\mathbf{g}$$
$$\nabla\cdot(K\,\mathbf{E}_\ell)=0 \;\Rightarrow\;
\nabla^2\phi_\ell = 0 \ \text{für konstantes } K$$

**Annahme.** Im Leaky-Dielectric-Grenzfall ist der Bulk elektroneutral; alle
freie Ladung sitzt auf $\Sigma$. Das setzt voraus, dass die Debye-Länge klein
gegen alle geometrischen Längen ist. Für ionische Flüssigkeiten liegt sie im
Sub-Nanometerbereich, für die *Emissionsstruktur* mit Krümmungsradien um 10 nm
ist das Verhältnis nicht mehr komfortabel — als Gültigkeitsgrenze zu prüfen und
zu dokumentieren.

**Quellen.** Melcher & Taylor (1969); Saville (1997).

### 3.2 Sprungbedingungen auf $\Sigma$

**Flächenladung.**
$$\sigma_s = \varepsilon_0 E_n^\mathrm{v} - \varepsilon_0\varepsilon_r E_n^\ell$$

**Ladungserhaltung, stationär.** Ohmscher Antransport aus dem Bulk plus
Konvektion entlang der Oberfläche gleich abgeführter Emissionsstromdichte:
$$K\,E_n^\ell + \nabla_s\cdot(\sigma_s\mathbf{u}_t) = j_\mathrm{evap}(E_n^\mathrm{v})$$

**Normalspannung.**
$$\gamma\kappa = (p_\ell - p_\mathrm{v}) - 2\mu\frac{\partial u_n}{\partial n}
+ \frac{\varepsilon_0}{2}\big(E_n^\mathrm{v}\big)^2
- \frac{\varepsilon_0\varepsilon_r}{2}\big(E_n^\ell\big)^2
- \frac{\varepsilon_0}{2}(1-\varepsilon_r)E_t^2$$

**Tangentialspannung.** Der tangentiale Maxwell-Zug treibt die Strömung — das
ist das definierende Merkmal des Leaky-Dielectric-Modells und im Prototyp gar
nicht vorhanden:
$$\sigma_s E_t = \mu\left(\frac{\partial u_t}{\partial n}
+ \frac{\partial u_n}{\partial t}\right)$$

**Kinematik mit Massenabfuhr.**
$$\rho\,(\mathbf{u}\cdot\hat{n}) = \dot{m}_\mathrm{evap}
= \frac{j_\mathrm{evap}}{(q/m)_\mathrm{emittiert}}$$

$E_t$ ist über $\Sigma$ stetig.

### 3.3 Ionenverdampfung

$$j_\mathrm{evap}(E) = \sigma_s\,\frac{k_BT}{h}\,
\exp\!\left[-\frac{\Delta G - \sqrt{e^3E/(4\pi\varepsilon_0)}}{k_BT}\right]$$

**Unbekannte:** $j_\mathrm{evap}$ als Funktion des lokalen Feldes; koppelt in
3.2 zurück.

**Annahmen.** (a) Thermisch aktiviertes Überwinden einer feldabsenkten Barriere.
(b) Eine Spezies, ein $\Delta G$. (c) Kein Rückfluss. (d) Der Vorfaktor
$k_BT/h$ ist eine nominale Versuchsfrequenz; die tatsächliche Zustandsdichte an
der Oberfläche steckt unbestimmt darin.

**Gültigkeitsbereich.** Felder um 1 V/nm. Die Schottky-Absenkung beträgt dort
1,20 eV, vergleichbar mit $\Delta G$ — daher die Steilheit. Für mehrere Spezies
(Monomer, Dimer, Trimer) ist die Summe über Spezies mit je eigenem $\Delta G$
und eigener Masse zu bilden; die Aufteilung ist ein Modelleingang, kein
Ergebnis.

**Der kritische Punkt für die Parameteranpassung.** $\Delta G$ geht exponentiell
ein. Publizierte Werte für EMI-BF4 streuen etwa 1,0–1,4 eV je nach
Auswertemethode. Der Code darf $\Delta G$ nicht als Naturkonstante führen,
sondern als anzupassenden Parameter mit anzugebendem Vertrauensbereich, und muss
die Empfindlichkeit $\partial\ln I/\partial\Delta G$ mit ausgeben.

**Quellen.** Iribarne & Thomson (1976); Loscertales & Fernández de la Mora
(1995); Anwendung auf IL-Menisken: Coffman et al. (2016, APL 109, 231602) ✓;
Gallud & Lozano (2022) ✓.

### 3.4 Hydraulische Speisung

Kopplung von Speisedruck und Volumenstrom über eine Impedanz:
$$\Delta p_\mathrm{feed} = Z_h\,Q,\qquad Q = \frac{\dot{m}}{\rho}$$

Für einen zylindrischen Kanal der Länge $L_c$ und des Radius $a$ im laminaren
Grenzfall $Z_h = 8\mu L_c/(\pi a^4)$ (Hagen-Poiseuille). Für poröse Emitter
stattdessen Darcy mit Permeabilität. Welcher Fall zutrifft, ist offene Frage 5.

Gallud & Lozano (2022) ✓ finden eine **minimale** hydraulische Impedanz als
notwendige Bedingung für statisch stabilen emittierenden Betrieb — die Speisung
ist also kein Randdetail, sondern mitbestimmend für die Existenz des
Betriebspunkts.

### 3.5 Regimeabgrenzung — drei getrennte Modelle

| Regime | Modell | Gültigkeitsbereich | Status |
|---|---|---|---|
| **PIR** (rein ionisch) | Stufe 3 vollständig, stationär, kein Jet | hohe Leitfähigkeit, kleiner Meniskus (0,5–3 µm nach Gallud & Lozano ✓), ausreichende hydraulische Impedanz | Zielmodell, physikalisch geschlossen |
| **Cone-Jet** | empirische/asymptotische Skalengesetze | $Q>Q_\mathrm{min}$; die Stromkorrelation ist an $\varepsilon_r\gtrsim40$ etabliert | **keine** selbstkonsistente Rechnung; getrennt zu kennzeichnen |
| **Mischbetrieb** | keines der beiden | — | nicht vorhersagbar; nur mit gemessener Aufteilung als Eingang |

Für den Cone-Jet:
$$I = f(\varepsilon_r)\sqrt{\frac{\gamma K Q}{\varepsilon_r}},\qquad
Q_\mathrm{min}\sim\frac{\gamma\varepsilon_0\varepsilon_r}{\rho K}$$
(Fernández de la Mora & Loscertales 1994). $f\approx18$ für
$\varepsilon_r\gtrsim40$; für ionische Flüssigkeiten mit $\varepsilon_r\approx12$
liegt das außerhalb des Etablierungsbereichs. Diese Formeln dürfen im Ausgabe-
format **nicht** gleichrangig neben der feldgekoppelten Rechnung stehen.

**Anforderung an die Implementierung:** getrennte Module, getrennte
Ausgabeblöcke, und bei jeder empirischen Korrelation ein maschinenlesbares Feld
`empirical=true` samt Gültigkeitsbereich, das in jede Ausgabedatei übernommen
wird.

---

## Stufe 4 — Strahltransport und Raumladung

### 4.1 Feldgleichung

$$\nabla^2\phi = -\frac{\rho_f}{\varepsilon_0}\quad\text{in }\Omega_\mathrm{v}$$

mit $\rho_f$ aus der Teilchenverteilung. **Das ist der Grund, warum die reine
BEM hier nicht ausreicht:** sie stellt Lösungen der homogenen Gleichung dar.

### 4.2 Verfahrenswahl

Drei Optionen wurden gegeneinander abgewogen:

| Option | Bewertung |
|---|---|
| Volumen-BEM (Randintegral + Volumenintegral über $\rho_f$) | Möglich, aber $O(N_\mathrm{vol}\!\cdot\!N_\mathrm{bnd})$ und erfordert dennoch eine Volumendiskretisierung samt regularisiertem Kern. Der Vorteil der BEM entfällt damit weitgehend. |
| FEM/FVM-Poisson mit BEM-Kopplung am Außenrand | Sauber, aber zwei Diskretisierungen und eine nichttriviale Kopplung. |
| **FEM- oder FVM-Poisson auf dem Rechteckgebiet der Skizze + PIC** | **Empfehlung.** Das Gebiet ist in der Skizze bereits vorgesehen („mesh"). Standardverfahren, gut validierbar (Child-Langmuir), Teilchenformfunktionen liefern die fehlende Regularisierung. |

Die vorhandene BEM bleibt für die raumladungsfreie Elektrostatik (Stufe 1) und
als **unabhängige Referenz** zur Verifikation des Volumenlösers erhalten — bei
$\rho_f=0$ müssen beide dasselbe liefern. Das ist ein wertvoller Testfall, den
man sonst nicht hätte.

### 4.3 Warum das Ring-Makroteilchenmodell ersetzt werden muss

Nicht als Stilfrage, sondern weil es kein wohlgestelltes Modell ist:

1. Das Eigenfeld eines unendlich dünnen geladenen Rings divergiert in der
   Meridianebene wie $1/d$, das Potential wie $\ln(1/d)$. Numerisch gemessen:
   `nan`/`inf` exakt am Ring, $|E| = 5{,}7\cdot10^{12}$ V/m bei 10⁻¹⁴ m Abstand.
2. Auch ohne Selbstwechselwirkung ist das Feld eines *benachbarten* Rings in
   einem Abstand von der Größenordnung des Teilchenabstands willkürlich, weil
   das Teilchen keine Ausdehnung hat. Ein Abschneideradius wäre ein freier
   Parameter ohne physikalische Festlegung.
3. Physikalisch ist ein Makroteilchen eine geglättete Verteilung. Die
   Formfunktion (cloud-in-cell o. ä.) eines PIC-Verfahrens ist genau die
   fehlende Regularisierung und zugleich der Grund, warum die Gitterlösung
   selbstkonsistent bleibt.

### 4.4 Spezies

Getrennte Populationen mit je eigenem $q/m$, eigener Startverteilung und
eigenem Anteil:

| Spezies | $q/m$ | Herkunft der Gewichtung |
|---|---|---|
| Monomer, Dimer, Trimer (Kation oder Anion) | $e/(m_\mathrm{ion} + n\,m_\mathrm{Paar})$ | Stufe 3, $j_\mathrm{evap}$ je Spezies |
| Tropfen | $I_\mathrm{cj}/(\rho Q)$ | Cone-Jet-Modul, **nicht** aus $j_\mathrm{evap}$ |

Der Prototyp gewichtet auch die Tropfen mit $j_\mathrm{evap}$ (Befund 6,
Faktor 8·10⁴). Das ist zu trennen.

### 4.5 Polarität

Das Vorzeichen der angelegten Spannung entscheidet, ob Kationen oder Anionen
emittiert werden. Zu unterscheiden sind mindestens:

* die emittierte Ionenmasse ($m_\mathrm{EMI^+}=111{,}17$ vs.
  $m_\mathrm{BF_4^-}=86{,}81$ g/mol; bei EMI-Im $280{,}15$ g/mol für das Anion),
* die Solvatationsenergie $\Delta G$ (im Allgemeinen polaritätsabhängig),
* das Vorzeichen von $E_n$ in $j_\mathrm{evap}$ — nur ein nach außen ziehendes
  Feld treibt Emission; der Prototyp benutzt $|E_n|$ und emittiert deshalb in
  beiden Polaritäten identisch.

### 4.6 Selbstkonsistenz und ihre Grenzen

Die Kopplung Raumladung $\to$ Feld $\to$ Emission $\to$ Raumladung ist zu
iterieren. Higuera (2008) ✓ kommt für den PIR zu dem Ergebnis, dass die
Raumladung *um die verdampfende Oberfläche herum* vernachlässigbar ist. Das ist
eine belegte Vereinfachung für den Einzelemitter, keine allgemeine — für
Emitterarrays und hohe Ströme ist sie erneut zu prüfen.

### 4.7 Kenngrößen

$$F = \sum_i \dot m_i v_i,\quad v_i=\sqrt{2(q/m)_i U_\mathrm{acc}},\qquad
I=\sum_i \dot m_i (q/m)_i,\qquad I_\mathrm{sp}=\frac{F}{\dot m\,g_0}$$
$$\eta_\mathrm{pol}=\frac{F^2}{2\dot m\,P},\qquad P = I\,U_\mathrm{acc}$$

Dazu Divergenzwinkel (stromgewichtete Quantile), Extraktorinterzeption je
Elektrode, Energieverteilung.

---

## 2.5 Kopplungsübersicht

```
Geometrie + Potentiale
        |
        v
  [1] Laplace  ------------------------------- E_n auf Sigma
        |                                            |
        |                                            v
        |                                   [2] Young-Laplace + Maxwell
        |                                            |
        |                                     Form Sigma (Fixpunkt/Newton)
        |                                            |
        +-------- gemeinsame Newton-Unbekannte ------+
                                |
                                v
              [3] Leaky dielectric + Ionenverdampfung
                 (u, p, phi_l, sigma_s, j_evap, Q, dp_feed)
                                |
                     emittierter Strom je Spezies
                                |
                                v
              [4] Poisson + PIC  ---> rho_f ---> zurück nach [1]/[3]
```

Stufen 1 und 2 sind eng gekoppelt (gemeinsames Newton-System). Stufe 3 fügt
weitere Unbekannte auf $\Sigma$ hinzu. Stufe 4 koppelt schwach zurück; für den
Einzelemitter im PIR nach Higuera (2008) ✓ vernachlässigbar, für Arrays nicht.
