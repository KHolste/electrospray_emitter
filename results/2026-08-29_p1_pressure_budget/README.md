# P1 — Druckhaushalt am Austritt — 2026-08-29

Alle Daten und Abbildungen sind in diesem Lauf frisch erzeugt;
`figures_provenance.txt` weist Commit und Zustand des Arbeitsbaums aus.

```sh
./build/es_feed examples/device_p1.cfg examples/feed_p1.cfg \
    results/2026-08-29_p1_pressure_budget meta.commit=$(git rev-parse HEAD)
python python/plot_feed.py results/2026-08-29_p1_pressure_budget
```

Modellvertrag, Vorzeichenkonvention und Gültigkeitsgrenzen:
[`docs/12_feed_model.md`](../../docs/12_feed_model.md).

## Gerechnet

    delta_p_exit = (p_reservoir - p_vacuum) - delta_p_hydrostatisch - delta_p_viskos
    delta_p_hydrostatisch = -rho g_z (z_exit - z_reservoir)
    delta_p_viskos        =  8 mu L Q / (pi R^4)

Der viskose Term ist der laminare Hagen-Poiseuille-Verlust eines geraden,
VOLLSTAENDIG GEFUELLTEN Kreiskanals. Die direkte Eingabe von delta_p_exit bleibt
als kontrollierter Modus erhalten, damit P3a und P3b reproduzierbar bleiben.

## NICHT gerechnet

Kein Stroemungsloeser. Kein kapillarer Aufstieg, keine bewegliche Kontaktlinie.
Kein Young-Winkel im Kanal -- der Kanal ist voll und hat dort keine freie
Oberflaeche; die einzige liegt an der gepinnten Austrittskante und wird von
P3a/P3b gerechnet, nicht hier ein zweites Mal. Der Vorrat ist eine
Randbedingung, kein Volumen: keine Entleerung, kein Gaspolster.
p_reservoir und Q bleiben EINGABEN.

## Kernzahlen (Rechenbeispiel, Stoffstatus illustrative)

| Groesse | Wert |
|---|---|
| hydraulischer Widerstand R_h | 4,535e16 Pa s/m^3 |
| viskoser Verlust bei Q = 1e-13 m^3/s | 4535 Pa = 0,502 gamma/a |
| Hydrostatik ueber 1 mm auf der Erde | 12,5 Pa = 1,39e-3 gamma/a |
| Hydrostatik ueber einen Bohrungsradius | 0,063 Pa = 6,9e-6 gamma/a |
| gamma/a | 9040 Pa |
| Q, bei dem der Haushalt die ganze Spanne \|Pi\| <= 2 durchfaehrt | rund 4e-13 m^3/s |
| Reynoldszahl ueber die ganze Studie | < 1e-2, also weit laminar |
| Einlauflaenge / L ueber die ganze Studie | < 1e-6 |

Ablesbar: der hydraulische Widerstand ist die Groesse, die den Austrittsdruck
auf der Kapillarskala bewegt; die Schwerkraft ist es selbst im Laborfall nicht.

## Pruefungen

`tests/test_feed.cpp`, alle bestanden. Wichtigste: das Simpson-Integral des
parabolischen Profils ueber den Querschnitt ergibt Q auf 1e-9 -- damit ist der
Widerstand gegen das Geschwindigkeitsfeld geprueft, aus dem er hergeleitet ist,
und nicht gegen eine zweite Kopie derselben Formel. Vorzeichen, Referenzhoehe
und Stroemungsrichtung sind Term fuer Term einzeln geprueft.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_budget_vs_flow.png` | Druckhaushalt Term fuer Term gegen den Volumenstrom, dieselbe Groesse auf der Kapillarskala, und die Gueltigkeitsgrenzen |
| `fig2_sensitivity.png` | viskoser Term gegen Kanallaenge und -radius, mit den Exponenten +1 und -4, die die geschlossene Form verlangt |
| `fig3_scales.png` | Druckskalen nebeneinander, und was der unveraenderte P3a-Loeser aus dem Haushalt macht |
