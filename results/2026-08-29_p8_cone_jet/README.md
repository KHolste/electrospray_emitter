# P8 — Cone-Jet und Tropfenbetrieb — 2026-08-29

Status: **`blocked`**. Eigener Modellvertrag, eigener Status, strikt getrennt
von der Pure-Ion-Emission (P5) und vom statischen Meniskus (P3a/P3b).

```sh
./build/es_cone_jet results/2026-08-29_p8_cone_jet meta.commit=$(git rev-parse HEAD)
python python/plot_cone_jet.py results/2026-08-29_p8_cone_jet
```

Vertrag und Herleitungen: [`docs/19_cone_jet.md`](../../docs/19_cone_jet.md).

## Warum die Trennung nicht kosmetisch ist

Ein Cone-Jet HAT einen Jet: ein schlankes geladenes Filament, das in Tropfen
zerfaellt. Ein gepinnter statischer Meniskus hat keinen. Die beiden sind
verschiedene Zustaende derselben Fluessigkeit und lassen sich nicht auseinander
berechnen.

## Der Blocker, zweiteilig und unabhaengig

**1. Die Physik ist nicht gerechnet.** Von sieben noetigen Teilmodellen fehlen
fuenf: Zweiphasenstroemung mit freier Oberflaeche, endliche Leitfaehigkeit mit
Oberflaechenladungstransport, belegte relative Permittivitaet, Zerfallsmodell
des Jets, geprueffte empirische Skalierung. Die Verfuegbarkeit ist gegen die
anderen Phasen GEMESSEN, nicht behauptet (`requirements.csv`).

**2. Die empirische Skalierung wird nicht uebernommen.** Der Auftrag erlaubt das
nur nach Pruefung von Originalgleichung, ERRATUM, Vorfaktor und
Gueltigkeitsbereich an der Quelle. Ganan-Calvo, Phys. Rev. Lett. 79, 217 (1997),
Erratum Phys. Rev. Lett. 85, 4193 (2000) -- **keiner der beiden Volltexte war
erreichbar**. Zahlen aus Suchtreffer-Schnipseln gelten hier nicht als Quelle,
und der Vorfaktor ist genau das, was das Erratum betrifft.

`cone_jet_current()` wirft und nennt BEIDE Gruende. Der alte `ConeJetModel` in
`emission.hpp` bleibt unangetastet und wird von hier nicht benutzt.

## Was stattdessen geliefert wird

Kennzahlen, die **Definitionen** sind und keine Literaturquelle brauchen. Jede
ist aus einer Bilanz hergeleitet, und der Test prueft sie gegen ihre eigene
Herleitung -- etwa: bei r* sind Ladungsrelaxations- und Kapillarzeit gleich,
geprueft auf 1e-12.

| Kennzahl | Wert bei a = 5 um, Q = 1e-13 m^3/s, E = 5e7 V/m |
|---|---|
| Ohnesorge | 1,96 |
| elektrische Bondzahl | 1,02 |
| Reynolds (Zulauf) | 4,5e-04 |
| Kapillarzahl (Zulauf) | 8,6e-04 |
| tau_e | **6,06e-11 s** -- selbstkonsistent (P3), Band 4,18e-11 ... 1,32e-10 s |
| r* | **5,37e-09 m** -- Band 4,19e-09 ... 9,04e-09 m |

**Zu tau_e und r*.** Eine fruehere Fassung fuehrte beide als *nicht berechenbar,
eps_r fehlt*. Das war vor der P3-Korrektur richtig. Beide haengen von eps_r
jedoch **nur ueber tau** ab:

```
r* = (gamma eps0^2 eps_r^2 / (rho K^2))^(1/3) = (gamma tau^2 / rho)^(1/3)
```

und tau ist seit P3 selbstkonsistent auf der **gemessenen** Dispersionskurve
geloest (eps_r = 10,664 bei f* = 2,627 GHz). Das ist **kein Ersatzwert**. Ein
**einzelner** eps_r-Wert fehlt weiterhin, `eps_r_status` bleibt
`MissingMaterialData`, und der Gesamtstatus der Diagnose ebenso. Die Identitaet
zwischen beiden Rechenwegen wird im Test gerechnet, nicht behauptet: r* aus tau
und r* aus eps_r(f*) stimmen auf 1e-12 ueberein.

Ablesbar: Oh ~ 2, auf dieser Laengenskala dominiert die Viskositaet die
Traegheit. Das ist eine DIAGNOSE, welche Physik ein Modell enthalten muesste --
und KEINE Aussage darueber, ob das Geraet in einem Cone-Jet-Modus laeuft.

Die elektrische Bondzahl ist die einzige Stelle, an der diese Diagnose etwas
beruehrt, das dieses Projekt gerechnet hat: sie ist das Verhaeltnis der beiden
Terme der P3b-Gleichgewichtsgleichung.

## NICHT gerechnet

Kein Cone-Jet-Strom, kein Jetradius, kein Tropfendurchmesser, keine
Tropfengroessenverteilung, keine Regimekarte.

## Abbildung

| Datei | Inhalt |
|---|---|
| `fig1_status.png` | Gueltigkeitskarte der Teilmodelle, die Zeitskalen nebeneinander, und die Kennzahlen, die Definitionen sind -- ausdruecklich als Diagnose bezeichnet |
