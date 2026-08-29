# P0 — numerische Bereinigung von P3b — 2026-08-29

Alle Daten und Abbildungen sind in diesem Lauf frisch erzeugt;
`figures_provenance.txt` weist Commit und Zustand des Arbeitsbaums aus.

```sh
./build/es_p3b_audit examples/device_p1.cfg examples/electrocapillary_p3b.cfg \
    results/2026-08-29_p0_p3b_audit meta.commit=$(git rev-parse HEAD)
python python/plot_p3b_audit.py results/2026-08-29_p0_p3b_audit
```

Herleitungen und alle Messungen im Wortlaut:
[`docs/11_p0_p3b_audit.md`](../../docs/11_p0_p3b_audit.md).
Der Volltext des Laufs steht in `report.txt`.

## Gerechnet

**Keine neue Physik.** Drei offene Punkte aus P3b werden gemessen statt
behauptet: die Lastprojektion, das Kanten-Gate und die gekoppelte
Netzkonvergenz. Eine gescheiterte Rechnung steht ueberall als `nan`, nie als
null.

## Die vier Kernbefunde

**1. Die uebergebene Last ist tatsaechlich stetig.** Der Sprung ueber die
Binraender faellt bei zehnfach kleinerem Probenabstand auf das 0,100-fache --
fuer eine Treppenfunktion waere er 1. Die integrierte Maxwell-Kraft ueberlebt
die Rekonstruktion exakt (0,0e+00 relativ gegen das rekonstruierte Flaechenmass;
4e-6 gegen 2 pi r ds -- beide Zahlen stehen da, nicht nur die guenstigere).
Die Segmentprojektion trifft die geschlossene Form der glatten Last mit der
Ordnung 1,98 und die der Singularitaet mit der Ordnung 1+beta. Kein leeres Bin,
keine Ausschlusszone, keine Kappung.

**2. `kTolExclusion` ist widerlegt, nicht getauscht.** Fuer p = C d^beta hat die
geprueffte Groesse einen geschlossenen Grenzwert, der die Netzweite gar nicht
enthaelt. Messung gegen Formel: hoechstens 2,3e-04 Abweichung. Ebene Form:
beta = -0,4456, gemessen 0,1110, Formel 0,1133. Die 11,3 % aus P3b sind also der
von beta festgelegte Wert und kein numerischer Zufall. Gewoelbte Form:
beta = +0,041, gemessen 0,0364 -- unter der Schranke, WEIL beta groesser ist.
Ein Kriterium, dessen Erfuellung so von der Form abhaengt, prueft keine
Netzkonvergenz. Die Schranke bleibt im Code deklariert und wird berichtet.

**3. Die Randsingularitaet setzt die Konvergenzrate der ganzen Rechnung.** Die
integrierte Maxwell-Kraft der ebenen Form konvergiert mit der Ordnung 0,541 --
und 1+beta = 0,554. Zweite Ordnung ist mit einer unverrundeten Kante nicht zu
haben.

**4. Das 1-%-Ziel ist NIRGENDS erreicht.** Vorab festgelegt: geschaetzter
Diskretisierungsfehler der feinsten Stufe unter 1 %.

| Punkt | Groesse | Ordnung | geschaetzter Fehler | Urteil |
|---|---|---|---|---|
| 1000 V | h/a | -- | -- | `NotInAsymptoticRange` |
| 1000 V | E_n kantenfern | -- | -- | `NotInAsymptoticRange` |
| 1000 V | Gesamtkraft | 0,581 | 6,07 % | `DiscretizationNotConverged` |
| 1400 V | h/a | 0,953 | 5,67 % | `DiscretizationNotConverged` |
| 1400 V | E_n kantenfern | 0,231 | 4,83 % | `DiscretizationNotConverged` |
| 1400 V | Gesamtkraft | 0,910 | 4,59 % | `DiscretizationNotConverged` |

Bei 1000 V sind h/a und E_n nicht einmal im asymptotischen Bereich: keine
Ordnung beobachtbar, also kein Fehler schaetzbar. Bekannt ist dort nur die
blosse Aenderung zwischen den beiden feinsten Stufen (0,32 %), und das ist keine
Fehlerschranke.

**Konsequenz:** die gekoppelten Ergebnisse von P3b sind QUALITATIV. Keine
Apexhoehe und keine Kraft aus diesem Ast traegt drei Stellen.

## NICHT gerechnet

Keine Emission, keine endliche Leitfaehigkeit, keine Stroemung, keine
Raumladung, keine Zeitabhaengigkeit, keine Stabilitaetsaussage, kein
Taylor-Kegel. delta_p_exit ist in diesem Lauf eine Eingabe (P1 aendert das
getrennt). Die Stoffwerte bleiben illustrative.

## Netzstufen

Das Kanten-Gate laeuft ueber die Stufen 1 bis 5. Stufe 5 (395 000 Knoten,
3,01 GiB Bandfaktorisierung) ist die feinste rechenbare; dafuer wurde
`DielectricSetup::memory_cap_bytes` konfigurierbar gemacht, Vorgabewert
unveraendert. Eine nicht rechenbare Stufe wird gemeldet, nicht still
weggelassen. Die gekoppelte Studie laeuft ueber die Stufen 1 bis 4 -- Stufe 5
gekoppelt waere rund 50 Faktorisierungen zu je 3 GiB.

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_load_projection.png` | Rohlast, Segmenttreppe und uebergebene Last getrennt, ein Ausschnitt, der die Stetigkeit sichtbar macht, und die kumulierte Kraft G(tau), aus der die uebergebene Last gebaut ist |
| `fig2_projection_error.png` | Krafterhaltung durch die ganze Kette und die gemessene Ordnung gegen die geschlossenen Formen |
| `fig3_exclusion.png` | `kTolExclusion`: Messung gegen geschlossene Form ueber sechs Verfeinerungen, daneben der Diskretisierungsfehler derselben Last |
| `fig4_gate.png` | Kanten-Gate: direkte Kraftfolge und die beiden Extrapolationen desselben Grenzwertes nebeneinander |
| `fig5_coupled.png` | gekoppelte Netzstudie -- ausdruecklich NICHT „Konvergenz" genannt, weil das Ziel verfehlt ist |
