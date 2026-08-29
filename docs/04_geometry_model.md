# 4. Parametrisches Geometriemodell

## 4.1 Lesart der Skizze

Die Skizze zeigt einen Längsschnitt mit horizontaler Achse. In der Notation
dieses Projekts liegt die Symmetrieachse bei $r=0$, die Achsrichtung ist $z$ und
zeigt vom Emitter zur Extraktionselektrode. Ursprung $z=0$ ist die
**Stirnebene des Emitters**.

Gelesen wurden:

| Beschriftung | Interpretation | Symbol |
|---|---|---|
| ∅₃ | Außendurchmesser am Fuß der Verjüngung | `d3` |
| ∅₁ | Außendurchmesser an der Stirnfläche | `d1` |
| ∅₂ | Bohrungsdurchmesser am Austritt = Kontaktlinienradius ×2 | `d2` |
| „height" | axiale Länge der verjüngten Struktur | `H` |
| „distance" | Stirnebene bis emitterseitige Extraktorfläche | `L` |
| „diameter" (rot) | Aperturdurchmesser der Extraktionselektrode | `D_a` |
| „thickness" (rot) | Dicke der Extraktionselektrode | `t_e` |
| grün, „ionische Flüssigkeit" | Speisekanal und Meniskus | — |
| „mesh" (äußeres Rechteck) | Berechnungsgebiet bzw. Gehäuse | `R_dom`, `z_min`, `z_max` |

Der linke, senkrecht schraffierte Balken ist als **Basisplatte** gelesen, aus
der die verjüngte Emitterstruktur herausragt und die radial bis zum Gehäuse
reicht. Sie liegt auf Emitterpotential.

```
   r
   ^
   |                                              [Gehaeuse / Domaenenrand]
R_dom +==================================================================+
   |  #                                                                  |
   |  #  Basisplatte (U_emitter)                                         |
   |  #                                                                  |
   |  #\                                                                 |
d3/2  #  \    Aussenkontur: Kegelstumpf (oder gekruemmtes Profil)        |
   |  #    \                                                             |
   |  #      \                                                     +####+|
   |  #        \______________________                             |####||
d1/2  #        |                      |   <- Stirnflaeche (Land)   |####||
   |  #        |######################|      Breite w=(d1-d2)/2    |####||
d2/2  #        +----------------------+                            +####+|
   |  #        :  ionische Fluessigkeit ) <- Meniskus              D_a/2 |
   |  #        :                       :\                                |
  0 +--#--------:-----------------------:--------------------------------+--> z
      z_min    -H                       0                L        L+t_e  z_max
                <--------- H ---------> <------ L ------>
```

## 4.2 Parametersatz

### Emitter

| Schlüssel | Bedeutung | Einheit | Zwang |
|---|---|---|---|
| `emitter.d1` | Außendurchmesser an der Stirnfläche | m | `d2 < d1 <= d3` |
| `emitter.d2` | Bohrungsdurchmesser am Austritt | m | `> 0` |
| `emitter.d3` | Außendurchmesser am Fuß | m | |
| `emitter.height` | axiale Länge der Verjüngung, `H` | m | `> 0` |
| `emitter.outer_profile` | `cone` \| `concave` \| `convex` | — | |
| `emitter.profile_exponent` | Formexponent für gekrümmte Profile | — | `> 0` |
| `emitter.bore_profile` | `cylindrical` \| `tapered` | — | offene Frage 1 |
| `emitter.d2_base` | Bohrungsdurchmesser am Fuß (nur `tapered`) | m | `>= d2` |
| `emitter.edge_radius_outer` | Verrundung Außenkante der Stirnfläche | m | `<= w` |
| `emitter.edge_radius_inner` | Verrundung Innenkante (Kontaktlinie) | m | `<= w` |
| `emitter.base_plate_thickness` | Dicke der Basisplatte | m | |
| `emitter.channel_length` | Speisekanallänge stromauf, für `Z_h` | m | |

Abgeleitet und mit auszugeben:

$$\alpha = \arctan\frac{d_3-d_1}{2H}\quad\text{(Kegelhalbwinkel)},\qquad
w = \frac{d_1-d_2}{2}\quad\text{(Stirnflächenbreite)},\qquad
r_c = \frac{d_2}{2}$$

Das Außenprofil als Funktion von $\zeta = (z+H)/H \in [0,1]$:

$$r_\mathrm{aussen}(\zeta) = \frac{d_3}{2} - \left(\frac{d_3-d_1}{2}\right)\cdot
\begin{cases}
\zeta & \text{cone}\\
\zeta^{\,p} & \text{convex}\\
1-(1-\zeta)^{\,p} & \text{concave}
\end{cases}$$

### Extraktionselektrode

| Schlüssel | Bedeutung |
|---|---|
| `extractor.distance` | `L`, Stirnebene bis emitterseitige Fläche |
| `extractor.aperture_diameter` | `D_a` |
| `extractor.thickness` | `t_e` |
| `extractor.edge_radius` | Verrundung der Aperturkante, `<= t_e/2` |
| `extractor.outer_radius` | radiale Ausdehnung; `auto` = bis Gehäuse |

### Gebiet und weitere Elektroden

| Schlüssel | Bedeutung |
|---|---|
| `domain.radius`, `domain.z_min`, `domain.z_max` | Berechnungsgebiet |
| `domain.boundary` | `open` (BEM-Fernfeld) \| `dirichlet` (Gehäuse) |
| `collector.enabled`, `collector.z`, `collector.radius` | optional, offene Frage 3 |

### Potentiale und Betrieb

`U_emitter`, `U_extractor`, `U_housing`, `U_collector`; Polarität ergibt sich
aus dem Vorzeichen von `U_emitter - U_extractor` und steuert die Speziesauswahl
(Kationen oder Anionen).

`feed.mode` = `pressure` \| `impedance` \| `flow_rate`, dazu `feed.delta_p`,
`feed.impedance` bzw. `feed.Q`, sowie `feed.model` = `poiseuille` \| `darcy`.

`wetting.mode` = `pinned` \| `contact_angle`, dazu `wetting.theta`.

`fluid.*` und `species.*` wie in
[02_model_specification.md](02_model_specification.md), Abschnitt 3.3/4.4.

## 4.3 Zwangsbedingungen, die vor dem Vernetzen geprüft werden

1. $0 < d_2 < d_1 \le d_3$; bei $d_1 = d_3$ ist die Struktur zylindrisch, das ist
   zulässig.
2. $H > 0$, $L > 0$, $t_e > 0$.
3. Verrundungsradien $\le w$ bzw. $\le t_e/2$.
4. $D_a > d_1$ — sonst liegt die Extraktorapertur radial innerhalb der
   Emitterstirnfläche. Das ist geometrisch darstellbar, aber selten gewollt;
   **Warnung**, kein Abbruch.
5. $L > $ maximal erreichbare Apexhöhe. Andernfalls kann der Meniskus die
   Elektrode berühren; **Abbruch**, weil das Modell diesen Fall nicht abdeckt.
6. Bei `wetting.mode = pinned`: `edge_radius_inner` muss klein gegen $r_c$ sein,
   sonst ist die Kontaktlinie nicht durch die Kante festgelegt und die
   Gibbs-Bedingung nicht anwendbar (siehe Spezifikation 2.1, Fall A).

## 4.4 Automatische Vernetzung

Anforderung: der Benutzer gibt kein Netz vor. Vorgesehen ist eine
zweistufige Erzeugung.

**Stufe A — Basisnetz aus Geometriemerkmalen.** Elementgröße aus lokaler
Krümmung und Merkmalsgröße:

$$h(x) = \min\Big(h_\mathrm{max},\ c_\kappa\,R_\mathrm{Krümmung}(x),\
c_f\,\ell_\mathrm{feature}(x),\ c_g\,\mathrm{dist}(x,\Gamma_\mathrm{nah})\Big)$$

mit einer beschränkten Wachstumsrate (typisch 1,2 je Element), damit keine
Größensprünge entstehen.

Vom Benutzer geforderte Verfeinerungszonen und ihre Merkmalsgröße:

| Zone | $\ell_\mathrm{feature}$ |
|---|---|
| Meniskusspitze | Apexkrümmungsradius (aus der laufenden Lösung, siehe Stufe B) |
| Kontaktlinie | `edge_radius_inner`, ersatzweise $w/10$ |
| Elektrodenkanten (Stirnfläche außen) | `edge_radius_outer` |
| Extraktorapertur | `extractor.edge_radius`, $D_a/40$ |
| Strahlkern | anfänglich $r_c$; ab Stufe 4 aus dem Strahlradius |

**Stufe B — a-posteriori-Verfeinerung.** Nach jeder Lösung wird verfeinert, wo
ein Fehlerindikator anschlägt:

* Rand (BEM): Sprung der Flächenladungsdichte zwischen Nachbarelementen,
  normiert; alternativ Vergleich mit der Lösung auf einem gröberen Netz.
* Volumen (Poisson/PIC): gradientenbasierter Indikator auf $\phi$, zusätzlich
  Teilchen-pro-Zelle-Statistik im Strahlbereich.
* Freie Oberfläche: Verschiebung zwischen zwei Newton-Schritten relativ zur
  lokalen Elementgröße.

Abbruch, wenn eine **Zielgröße** netzunabhängig ist — nicht das Netz selbst.
Als Zielgrößen sind mindestens zu überwachen: Apexfeld, emittierter Strom,
emittierende Fläche, Divergenzwinkel. Genau diese Prüfung fehlte im Prototyp
(Befund 9).

## 4.5 Verhältnis zur vorhandenen Geometrieschicht

`make_capillary_open` bildet den Sonderfall $d_1 = d_3$ (zylindrisches
Röhrchen, keine Verjüngung) mit scharfer Innenkante ab. Der neue Aufbau ist
eine Verallgemeinerung; die vorhandenen Verifikationsgeometrien
(`make_sphere`, `make_prolate_spheroid`) bleiben als Testfälle unverändert
erhalten.

## 4.6 Was aus der Skizze nicht hervorgeht

Zusammenfassung der offenen Fragen 1–5 aus [README.md](README.md), an dieser
Stelle mit der jeweils getroffenen Vorannahme:

| Frage | Vorannahme im Entwurf | Auswirkung, falls falsch |
|---|---|---|
| Bohrung zylindrisch oder verjüngt? | zylindrisch, ∅₂ konstant | ändert Speiseimpedanz und Feld nahe der Kontaktlinie |
| Rechteck = Gehäuse oder nur Domäne? | Domäne, `boundary = open` | bei geerdetem Gehäuse ändert sich das Feld messbar |
| Kollektor stromab? | nicht vorhanden | zusätzliche Elektrode, Strahlaufweitung |
| Stirnfläche scharf oder verrundet? | scharf, Kontaktlinie gepinnt | bei Verrundung gilt Fall B (Kontaktwinkel), andere Lösungsmenge |
| massiv oder porös? | massiv, innen benetzt, Poiseuille | bei porös gilt Darcy, andere `Z_h`-Skalierung |
