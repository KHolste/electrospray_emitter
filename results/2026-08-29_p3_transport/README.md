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

## Der Befund

**tau ist mit den belegten Stoffdaten dieses Projekts NICHT berechenbar.**
eps_r ist fuer EMI-BF4 `MissingMaterialData` (P2, docs/13): keine der vier
Quellen nennt Reinheit und Wassergehalt, und die einzige mehrpunktige misst von
1 bis 18 GHz -- eine Dispersionsmessung, kein Gleichstromwert.

`judge_conductor_limit()` liefert deshalb `MissingMaterialData` statt einer
Zahl. **Der Aequipotentialansatz von P3b ist damit eine Annahme und kein
nachgewiesener Grenzfall.**

Mit einem ausdruecklich UNBELEGTEN eps_r = 12,8 waere
tau = 7,27e-11 s gegen t_kap = 1,72e-06 s und t_vis = 3,37e-06 s, also
Verhaeltnisse von 2e4 bis 5e4. Die Annahme ist plausibel -- belegt ist sie
nicht, und das steht in jeder Ausgabe.

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
| `fig3_time_scales.png` | die Zeitskalen nebeneinander, mit dem, was belegt ist, und dem, was nicht |
