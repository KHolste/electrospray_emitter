# 19 — P8: Cone-Jet und Tropfenbetrieb

**Stand: 2026-08-29.** Status: **`blocked`**. Eigener Modellvertrag, eigener
Status, strikt getrennt von der Pure-Ion-Emission (P5) und vom statischen
Meniskus (P3a/P3b).

Code: `include/es/cone_jet_contract.hpp`, `src/cone_jet_contract.cpp`,
`apps/es_cone_jet.cpp`, `tests/test_cone_jet_contract.cpp`,
`python/plot_cone_jet.py`.
Ergebnisse: `results/2026-08-29_p8_cone_jet/`.

---

## 19.1 Warum die Trennung nicht kosmetisch ist

Ein Cone-Jet **hat einen Jet**: ein schlankes geladenes Filament, das in Tropfen
zerfällt. Ein gepinnter statischer Meniskus hat keinen. Die beiden sind
verschiedene Zustände derselben Flüssigkeit und lassen sich nicht auseinander
berechnen. Der Prototyp druckte den Tropfenstrom der Cone-Jet-Korrelation neben
den Ionenstrom, als wären es Alternativen eines Modells — sie sind es nicht.

## 19.2 Der Blocker, zweiteilig und unabhängig

**Erstens: die Physik ist nicht gerechnet.** `cone_jet_requirements()` listet
die nötigen Teilmodelle und misst ihre Verfügbarkeit gegen die anderen Phasen:

| Teilmodell | vorhanden? |
|---|---|
| Zweiphasenströmung mit freier Oberfläche | **nein** — P4 spezifiziert den Vertrag (docs/15) |
| endliche Leitfähigkeit mit Oberflächenladungstransport | **nein** — P3 spezifiziert ihn (docs/14, 14.5) |
| belegte relative Permittivität | **nein** — `MissingMaterialData` (docs/13) |
| belegte Leitfähigkeit | ja — P2 |
| belegte Oberflächenspannung, Dichte, Viskosität | ja — P2 |
| Zerfallsmodell des Jets | **nein** |
| geprüfte empirische Skalierung | **nein** |

**Zweitens: die empirische Skalierung wird ebenfalls nicht übernommen.** Der
Auftrag erlaubt das nur nach Prüfung von Originalgleichung, **Erratum**,
Vorfaktor und Gültigkeitsbereich an der Quelle. Die Skalierung von Gañán-Calvo
steht in *Phys. Rev. Lett.* **79**, 217 (1997), das Erratum in *Phys. Rev.
Lett.* **85**, 4193 (2000). **Keiner der beiden Volltexte war in diesem Lauf
erreichbar.** Die Zahlen, die eine Suche zurückgibt, sind Schnipselwerte; dieses
Projekt akzeptiert sie nicht als Quelle. Der Vorfaktor ist genau das, was das
Erratum betrifft — ihn aus einem Schnipsel zu übernehmen wäre die eine Sache,
die hier ausdrücklich verboten ist.

`cone_jet_current()` wirft deshalb `NotImplementedInThisPhase` und nennt **beide**
Gründe. Der alte `ConeJetModel` in `emission.hpp` bleibt unangetastet und wird
von hier aus nicht benutzt.

## 19.3 Was stattdessen geliefert wird: Kennzahlen, die Definitionen sind

Jede der folgenden Größen ist aus einer Bilanz **hergeleitet** und braucht keine
Literaturquelle. Die Herleitung steht jeweils im Header, und der Test prüft sie
gegen ihre eigene Definition — nicht gegen einen Literaturwert.

| Kennzahl | Definition | drückt aus |
|---|---|---|
| $\tau_e=\varepsilon_0\varepsilon_r/K$ | aus $\dot\rho+(K/\varepsilon)\rho=0$ | Ladungsrelaxation |
| $t_c=\sqrt{\rho a^3/\gamma}$ | Trägheit gegen Kapillardruck | kapillar-inertial |
| $t_v=\mu a/\gamma$ | viskose Spannung gegen Kapillardruck | viskokapillar |
| $\mathrm{Oh}=\mu/\sqrt{\rho\gamma a}$ | $t_v/t_c$ | Viskosität gegen Trägheit |
| $r^*=(\gamma\varepsilon_0^2\varepsilon_r^2/(\rho K^2))^{1/3}$ | $\tau_e=t_c(r^*)$ | Größe, ab der eine Struktur der Ladungsrelaxation nicht mehr davonläuft |
| $\mathrm{Bo}_E=\varepsilon_0E^2a/(2\gamma)$ | Maxwell-Druck gegen $\gamma/a$ | **genau die beiden Terme, die P3b bilanziert** |
| $\mathrm{Re}=2\rho Q/(\pi R\mu)$ | Zulaufströmung | Trägheit gegen Viskosität |
| $\mathrm{Ca}=\mu u/\gamma$ | Zulaufströmung | Viskosität gegen Kapillarität |

**$r^*$ ist ausdrücklich nicht „der Jetradius" und nicht das $Q_0$ irgendeiner
publizierten Skalierung** — jene tragen Vorfaktoren, die dieses Projekt nicht
geprüft hat.

Die elektrische Bondzahl ist die einzige Stelle, an der diese Diagnose etwas
berührt, das dieses Projekt tatsächlich gerechnet hat: sie ist das Verhältnis
der beiden Terme der P3b-Gleichgewichtsgleichung.

## 19.4 Die Diagnose mit den belegten Daten

Mit $a=5\ \mu$m, den gewählten Quellen aus P2 und $E=5\cdot10^7$ V/m:

| Größe | Wert |
|---|---|
| $\mathrm{Oh}$ | 1,96 |
| $\mathrm{Bo}_E$ | 1,02 |
| $\mathrm{Re}$ (bei $Q=10^{-13}$ m³/s) | 4,5·10⁻⁴ |
| $\mathrm{Ca}$ | 8,6·10⁻⁴ |
| $\tau_e$, $r^*$ | **nicht berechenbar** — $\varepsilon_r$ fehlt |

Ablesbar: $\mathrm{Oh}\approx2$ — auf dieser Längenskala dominiert die
Viskosität die Trägheit. Das ist eine **Diagnose**, welche Physik ein Modell
enthalten müsste, und **keine** Aussage darüber, ob das Gerät in einem
Cone-Jet-Modus läuft.

## 19.5 Abhängigkeit von den Teilmodellen, dokumentiert

Ein Cone-Jet-Strom hinge ab von: dem **Volumenstrom** (er setzt Jetradius und
Verweilzeit), der **Viskosität** (sie setzt die Streckungsrate im Jet), der
**Oberflächenspannung** (sie hält den Jet zusammen und treibt seinen Zerfall),
der **Leitfähigkeit** (sie liefert die Ladung an die Oberfläche) und der
**elektrischen Belastung** (sie zieht den Kegel). Von diesen fünf sind vier
belegt und eine — über $\varepsilon_r$ die Ladungsrelaxation — nicht. Aber
selbst mit allen fünf fehlte die Strömung, in der sie wirken.

## 19.6 Was P8 ausdrücklich nicht enthält

Kein Cone-Jet-Strom, kein Jetradius, kein Tropfendurchmesser, keine
Tropfengrößenverteilung, keine Regimekarte. Keine Reaktivierung des alten
Cone-Jet-Blocks. Keine Übernahme einer empirischen Skalierung.
