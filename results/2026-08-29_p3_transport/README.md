# P3 — endliche Leitfähigkeit und Zulaufströmung — 2026-08-29

Status: **`validated_subset`**.

```sh
./build/es_transport results/2026-08-29_p3_transport meta.commit=$(git rev-parse HEAD)
python python/plot_transport.py results/2026-08-29_p3_transport
```

Modellvertrag und die Spezifikation dessen, was fehlt:
[`docs/14_transport.md`](../../docs/14_transport.md).

## Gerechnet und geprueft

**1. Voll ausgebildete Rohrstroemung.** Die Reduktion auf ein Skalarproblem ist
fuer das GERADE Rohr exakt: die Kontinuitaet erzwingt du_z/dz = 0, der
konvektive Term verschwindet identisch, und die z-Impulsgleichung IST der
achsensymmetrische Laplace-Operator.

| Pruefung | Ergebnis |
|---|---|
| Volumenstrom gegen Hagen-Poiseuille (81 radiale Knoten) | 2,6e-05 relativ |
| Netzkonvergenz | Ordnung 2,00 (gemessen) |
| hydraulischer Widerstand gegen die FORMEL von P1 | 2,6e-05 relativ, kein gemeinsamer Code |
| Linearitaet in dp/dz, Umkehrung in mu | exakt |
| ohne Druckgradient | Geschwindigkeit exakt null |

**2. Stationaerer Leitungsstrom.** div(sigma grad phi) = 0 ist dieselbe
Gleichung wie die Elektrostatik; die Knotenreaktion ist dann ein STROM.

| Pruefung | Ergebnis |
|---|---|
| Strom gegen I = sigma A V / L | 3,7e-11 relativ |
| Potential gegen die exakte lineare Loesung | 3,9e-12 |
| Ohmsches Gesetz Term fuer Term (V, sigma, A, L, Polaritaet) | alles geprueft |
| **kein Strom durch die zero-flux-Flaeche** | < 1e-12 |

Die letzte Zeile ist die physikalische Aussage: OHNE EMISSION darf kein
stationaerer Normalstrom durch die freie Oberflaeche flieszen, weil die Ladung
dort nirgendwohin koennte. Die richtige Bedingung ist j.n = 0, und dass der
Loeser sie einhaelt, wird gemessen.

**3. Ladungsrelaxation.** tau = eps0 eps_r / sigma, die geschlossene
Zerfallsfunktion, und die Pruefung, dass sie die ODE d rho/dt = -rho/tau
tatsaechlich loest (zentrale Differenz der zurueckgegebenen Funktion, nicht eine
zweite Kopie der Exponentialfunktion).

## Der Befund — korrigiert

Eine fruehere Fassung schrieb hier: *„tau ist mit den belegten Stoffdaten dieses
Projekts NICHT berechenbar“*, und begruendete das damit, dass die einzige
mehrpunktige Permittivitaetsquelle „von 1 bis 18 GHz misst -- eine
Dispersionsmessung, kein Gleichstromwert“. **Diese Begruendung war ein
Lesefehler der Formel.**

### Welche Permittivitaet in tau_q gehoert

Fuenf verschiedene Groessen heissen bei einer leitfaehigen ionischen Fluessigkeit
„die Permittivitaet“. Sie werden jetzt getrennt gefuehrt
(`permittivity_concepts.csv`):

| | Groesse | im Datensatz |
|---|---|---|
| (1) | statische/niederfrequente **Schein**permittivitaet — vom Leitungsbeitrag und von der Elektrodenpolarisation beherrscht, eine Eigenschaft der **Messzelle** | **kein einziger Punkt** — geprueft, nicht behauptet |
| (2) | intrinsische statische Permittivitaet eps_s — aus Mikrowellenspektren **extrapoliert**, nicht bei null Hertz gemessen | 3 Punkte: 12,8 / 12,9 / 13,6 |
| (3) | frequenzabhaengige komplexe Permittivitaet eps*(f) | 18 Punkte, 1–18 GHz, mit Unsicherheiten |
| (4) | Elektrodenpolarisation — Messartefakt, macht (1) unbrauchbar | — |
| (5) | DC-Leitfaehigkeit K — eigene Groesse, eigene Quelle | 1,5584 S/m, `measured` |

Aus Ladungserhaltung und Gauss folgt d rho_f/dt + (K/(eps0 eps_r)) rho_f = 0,
wobei eps_r die **gebundene** Ladungsantwort ist — die Polarisation, die dem Feld
folgt, *waehrend* die freie Ladung zerfaellt. Die freie Ladung zerfaellt auf der
Zeitskala tau selbst; ihr Spektrum liegt also bei f\* = 1/(2·pi·tau). Da eps_r
dispersiv ist, ist die Gleichung fuer tau **implizit**:

```
tau = eps0 * eps_r(f*) / K ,     f* = 1 / (2 pi tau)
```

Fuer diese Fluessigkeit liegt f\* bei einigen GHz — **also genau dort, wo
gemessen wurde**. Die 1–18-GHz-Daten sind nicht die falsche Frequenz, sondern die
richtige. Sie umgekehrt als Gleichstromwert zu uebernehmen waere ebenso falsch.

### Die Loesung auf der gemessenen Kurve

Fixpunktiteration auf den **Messpunkten**, logarithmisch zwischen benachbarten
Messfrequenzen interpoliert. Es wird **keine Dispersionsfunktion angepasst** und
ausserhalb des Messbereichs nicht extrapoliert.

| | |
|---|---|
| f\* | **2,627 GHz** — innerhalb des gemessenen Bereichs |
| eps_r(f\*) | 10,664 |
| K | 1,5584 S/m |
| **tau** | **6,059e-11 s** (18 Iterationen, Residuum 0) |
| zum Vergleich mit eps_s = 13,6 | 7,727e-11 s, also **+27,5 %** |

### Ein Einzelwert bleibt fehlend — ein Band ist belegt

Keine der vier Permittivitaetsquellen nennt Reinheit **und** Wassergehalt. Die
Auswahlregel von P2 waehlt weiterhin **keine** aus, und die Einzelwertabfrage
`judge_conductor_limit()` schlaegt weiterhin geschlossen fehl — der Lauf prueft
das ausdruecklich und setzt sonst `exit_code = 2`.

Belegt ist stattdessen ein **Band** und die Empfindlichkeit darueber:
eps_r in [7,7 ; 13,6], K in [0,91 ; 1,63] S/m. tau waechst mit eps_r und faellt
mit K; die fuer die Naeherung schlechteste Ecke ist (eps_hi, K_lo) — gerechnet
und geprueft, nicht behauptet.

| Ecke | tau [s] | t_kap/tau | t_vis/tau |
|---|---|---|---|
| eps_lo, K_hi | 4,183e-11 | 4,12e4 | 8,05e4 |
| eps_lo, K_lo | 7,492e-11 | 2,30e4 | 4,49e4 |
| eps_hi, K_hi | 7,388e-11 | 2,33e4 | 4,56e4 |
| **eps_hi, K_lo** (schlechteste) | **1,323e-10** | **1,30e4** | **2,54e4** |
| selbstkonsistent | 6,059e-11 | 2,84e4 | 5,56e4 |

mit t_kap = 1,722e-6 s, t_vis = 3,367e-6 s und einer geforderten Schranke von
100.

### Was sich damit aendert

**Der Aequipotentialansatz von P3b ist fuer die STATISCHE Form nicht mehr eine
Annahme, sondern ueber das gesamte begruendete Band belegt** — an der
unguenstigsten Ecke noch um mehr als zwei Groessenordnungen ueber der Schranke.
Es wird dafuer **kein einziger unbelegter eps_r-Wert** benutzt: das Ergebnis
haengt nicht davon ab, welchen Wert im Band man naehme.

**Was das nicht sagt.** Nichts ueber emittierenden Betrieb — dort ist die
Prozesszeit die Transitzeit durch die Emissionszone, und die ist nicht gerechnet.
Nichts ueber die tangentiale Traktion q_s·E_t (siehe unten). Und eps_r ist
weiterhin keine belegte Zahl: fuer jede Rechnung, die einen *Wert* statt eines
*Verhaeltnisses* braucht, fehlt er.

## NICHT gerechnet

Kein allgemeiner Stokes-Loeser: keine Einlaufstroemung, keine
Druck-Geschwindigkeits-Kopplung, keine gekruemmte Geometrie, keine freie
Oberflaeche. Keine gekoppelte finite-conductivity-Meniskusrechnung -- ihre
Randbedingungen stehen in docs/14, Abschnitt 14.5, als Spezifikation, und der
Code enthaelt KEINE Naeherung, die so taete, als waere sie geloest. Keine
Emission. Keine Zeitintegration einer Form. Keine Raumladung. Keine
Stabilitaetsaussage.

Besonders: bei endlicher Leitfaehigkeit ist E_t auf der Oberflaeche ungleich
null, und q_s E_t ist eine TANGENTIALE Traktion, die die Fluessigkeit in
Bewegung setzt. Diese Traktion existiert im Perfect-Conductor-Grenzfall gar
nicht. Ein endlich leitender Meniskus ist deshalb keine Korrektur des
P3b-Ergebnisses, sondern ein anderes Problem.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_pipe_flow.png` | geloestes Profil gegen Hagen-Poiseuille, die Differenz, und die Netzkonvergenz von Volumenstrom und hydraulischem Widerstand |
| `fig2_charge.png` | Ladungsrelaxation, stationaerer Leitungsstrom, und die Pruefung, dass kein Strom durch eine zero-flux-Flaeche flieszt |
| `fig3_time_scales.png` | welche Permittivitaet in tau_q gehoert: die gemessene Dispersion mit f\*, die implizite Gleichung als zwei sich schneidende Kurven, und die Zeitskalen mit tau als Band |

## Dateien

`permittivity_concepts.csv` -- die fuenf Groessen, die "die Permittivitaet"
heissen, und welche in tau_q gehoert. `permittivity_points.csv` -- jeder
Messpunkt mit Frequenz, Konzept und Zulaessigkeit.
`self_consistency.csv` -- die beiden Kurven, deren Schnittpunkt tau IST.
`relaxation.csv` -- die selbstkonsistente Loesung und alle vier Bandecken je
Prozesszeit. `time_scales.csv` -- tau mit lo/hi als Band.
