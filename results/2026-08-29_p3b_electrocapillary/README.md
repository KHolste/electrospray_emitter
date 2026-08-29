# P3b — selbstkonsistentes statisches Elektro-Kapillargleichgewicht — 2026-08-29

Alle Daten und Abbildungen sind in diesem Lauf frisch erzeugt; `figures_provenance.txt`
weist `releasable=yes`, `working_tree_dirty=no` und den Code-Commit aus.

```sh
./build/es_electrocapillary examples/device_p1.cfg examples/electrocapillary_p3b.cfg \
    results/2026-08-29_p3b_electrocapillary meta.commit=$(git rev-parse HEAD)
python python/plot_electrocapillary.py results/2026-08-29_p3b_electrocapillary
```

Modellvertrag, Vorzeichenherleitung, Gate und die vier dokumentierten Fehler:
[`docs/10_electrocapillary_model.md`](../../docs/10_electrocapillary_model.md).

## Gerechnet

    gamma kappa(s) = delta_p_exit + eps0 E_n(s)^2 / 2

E_n ist das vakuumseitige Normalfeld auf der freien Oberflaeche der Fluessigkeit,
die ein idealer Leiter auf V_emitter bleibt; die Elektrostatik ist das
unveraenderte dielektrische P2c-Problem.

## NICHT gerechnet

Keine Emission, keine endliche Leitfaehigkeit, keine Stroemung, keine
Viskositaet, keine Raumladung, keine Zeitabhaengigkeit, keine dynamische
Stabilitaet, kein Taylor-Kegel-Onset, kein Cone-Jet, keine Schwerkraft.
delta_p_exit bleibt eine Eingabe. Wo ein Ast endet, ist die Stelle, an der
DIESER Loeser stehen bleibt.

## Kernzahlen

| Groesse | Wert |
|---|---|
| Kanten-Gate | bestanden fuer Kontaktwinkel 0 bis 71,8 Grad |
| Singularitaetsexponent beta | +0,308 / +0,054 / -0,292 / -0,444 / -0,611 / -1,095 |
| Nullfeld gegen P3a | 0,000e+00 (bitgleich) |
| Polaritaet, Formunterschied | 0,000e+00 |
| Ladung sigma gegen FEM-Reaktionen | 1,4e-03 |
| Aeste | 10 Punkte bis 1439,5 V; 10 bis -1439,5 V; 8 bis 1236,3 V; einer abgelehnt |

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_moving_mesh.png` | konformes bewegliches Volumennetz, ebene/positive/negative Form, Detail von Apex und Austrittskante |
| `fig2_field.png` | Potential und Feldstaerke derselben Formen, identische Farb- und Laengenskalen |
| `fig3_edge_field.png` | E_n(d) und p_M(d) an der Kante ueber vier Netzstufen, mit Singularitaetsfit |
| `fig4_gate.png` | integrierte Maxwell-Kraft und schwache Projektion gegen Netzstufe und Ausschlussdistanz |
| `fig5_shapes_vs_voltage.png` | selbstkonsistente Meniskusformen ueber der Spannung |
| `fig6_branch.png` | Apexhoehe, kantenfernes Feld, Kruemmung und mechanisches Residuum ueber der Spannung |
| `fig7_convergence.png` | Netz- und Kopplungskonvergenz an drei Betriebspunkten |
