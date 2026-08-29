# P4 — zeitabhängige freie Oberfläche — 2026-08-29

Status: **`infrastructure_only`**.

```sh
./build/es_kinematics results/2026-08-29_p4_kinematics meta.commit=$(git rev-parse HEAD)
python python/plot_kinematics.py results/2026-08-29_p4_kinematics
```

Vollstaendiger Vertrag des fehlenden Loesers:
[`docs/15_free_surface_dynamics.md`](../../docs/15_free_surface_dynamics.md).

## Die Vorbedingung ist NICHT erfuellt

Ein dynamischer Meniskus ist nicht implementiert, und zwar strukturell:

1. Die Stroemung aus P3 hat **keine freie Oberflaeche**. Ihre Exaktheit beruht
   auf du_z/dz = 0, was eine sich verformende Oberflaeche zerstoert.
2. Der Ladungstransport aus P3 ist stationaer und kennt keine Flaechenladung.
   Es fehlt die **tangentiale Traktion q_s E_t**, die im
   Perfect-Conductor-Grenzfall gar nicht existiert und die einen endlich
   leitenden Meniskus antreibt.
3. eps_r ist `MissingMaterialData` (P2).

`solve_dynamic_meniscus()` wirft `NotImplementedInThisPhase` und nennt alle
drei Gruende; der Test prueft, dass die Meldung sie wirklich nennt.

## Was implementiert und geprueft ist

Genau ein Baustein: die **kinematische Randbedingung**
dx/dt . n = u . n auf einem VORGESCHRIEBENEN Feld. Das ist Kinematik, keine
Dynamik -- nichts wird fuer u geloest und keine Kraft ausgewertet.

Zwei Felder mit bekannter exakter Lagrange-Abbildung, die zwei verschiedene
Fehlerarten trennen:

| Feld | div u | exakte Abbildung | prueft |
|---|---|---|---|
| u = alpha x | 3 alpha | x -> x e^{alpha t} | das Volumen folgt einer BEKANNTEN Aenderung |
| u = (-alpha r/2, alpha z) | **0 exakt** | (r,z) -> (r e^{-alpha t/2}, z e^{alpha t}) | das Volumen bleibt EXAKT erhalten, waehrend sich die Form stark aendert |

| Pruefung | Ergebnis |
|---|---|
| Nullfeld bewegt nichts | exakt null, bitgenau |
| Dilatation, Volumen gegen V0 e^{3 alpha T} | 1,1e-15 |
| Dilatation, jeder Knoten gegen die exakte Abbildung | 3,9e-15 von R |
| beobachtete Zeitordnung (80 -> 160 Schritte) | **4,0** (RK4) |
| Squeeze, Volumenerhaltung | 2,1e-15 |
| Squeeze, Form gegen die exakte Abbildung | 6,8e-15 von R |
| gepinnte Kontaktlinie | bitgenau unveraendert |
| Apex auf der Achse | exakt r = 0 |

## Drei Verbote, als Abwesenheit durchgesetzt

Keine Mobilitaet. Keine kuenstliche Daempfung. Keine Stabilitaetsaussage. Es
gibt im Code keinen Koeffizienten, der eine Form mit einer gewaehlten Rate
relaxiert.

## Eine gemessene Grenze

Die tangentiale Umverteilung ist Netzbewegung: sie schiebt Knoten entlang der
Flaeche und kann die Flaeche im Kontinuum nicht aendern. Auf einem POLYGONZUG
schneidet das Neuabtasten Ecken ab. Gemessen ueber 900 Schritte:

| Knoten | 41 | 81 | 161 | 321 |
|---|---|---|---|---|
| \|dV\|/V | 4,1e-3 | 2,0e-3 | 9,1e-4 | 3,5e-4 |

**Beobachtete Ordnung 1,18 -- nicht 2.** Ein dynamischer Loeser, der in jedem
Schritt umverteilt, braucht deshalb eine gekruemmte Rekonstruktion oder muss den
Verlust bilanzieren. Das steht hier, statt weggeschaerft zu werden.

## Der Zeitschrittvertrag lehnt ab, statt zu warnen

Kein Knoten ueber die Achse, kein Segment kollabiert, keine Knotenbewegung ueber
0,25 der kuerzesten Segmentlaenge. Beim Squeeze VERSCHAERFT sich die dritte
Schranke im Lauf, weil die Segmente am Aequator schrumpfen: 200 Schritte bleiben
bei Schritt 180 mit `StepTooLarge` stehen. Das ist der Vertrag bei der Arbeit.

Die Schranke ist eine DISKRETISIERUNGS-Schranke. Hier wird keine Kraft
integriert, also gilt keine kapillare oder akustische Zeitschrittgrenze; die
gehoeren zu dem Loeser, den es nicht gibt.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_advection.png` | die advektierte Oberflaeche gegen die exakte Lagrange-Abbildung, fuer beide Felder ueber fuenf Zeitpunkte |
| `fig2_convergence.png` | Zeitkonvergenz gegen die exakte Abbildung, die Volumenbilanz gegen die EXAKTE Aenderung, und die gemessene Ordnung der Netzbewegung |
