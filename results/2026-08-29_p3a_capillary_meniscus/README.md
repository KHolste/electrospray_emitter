# P3a — statischer Kapillarmeniskus ohne elektrisches Feld — 2026-08-29

Alle Daten und Abbildungen in diesem Ordner sind in diesem Lauf frisch erzeugt
und jede Abbildung nennt den tatsächlichen HEAD-Commit; die Herkunft steht in
`figures_provenance.txt`.

Reproduktion:

```sh
./build/es_capillary examples/device_p1.cfg examples/capillary_p3a.cfg \
    results/2026-08-29_p3a_capillary_meniscus meta.commit=$(git rev-parse HEAD)
python python/plot_capillary.py results/2026-08-29_p3a_capillary_meniscus
```

Die beiden benutzten Konfigurationsdateien liegen als `device_p1.cfg` und
`capillary_p3a.cfg` mit im Ordner.

Modellvertrag, Herleitung und Prüfungen im Einzelnen:
[`docs/09_capillary_model.md`](../../docs/09_capillary_model.md).

---

## Was gerechnet wurde

Die Form der freien Flüssigkeitsoberfläche an der Austrittskante, gehalten
**allein** von der Oberflächenspannung gegen einen vorgegebenen
Flüssigkeitsdruck:

```
gamma (dpsi/ds + sin(psi)/r) = delta_p_exit,   delta_p_exit = p_fluessig - p_vakuum
```

achsensymmetrisch, in Bogenlängenparametrisierung `(r(s), z(s), psi(s))` vom
Apex zur gepinnten Kante. Der Grenzwert `sin(psi)/r -> dpsi/ds` auf der Achse
ist analytisch eingesetzt (`dpsi/ds|Apex = kappa/2`); es wird nirgends durch
`r = 0` geteilt.

**Vorzeichen:** `delta_p_exit > 0` wölbt nach `+z` (zum Extraktor),
`= 0` ergibt eine **exakt** ebene Fläche, `< 0` zieht die Oberfläche in die
Bohrung. Die Kontaktlinie ist an der scharfen Austrittskante
`r = phi_2/2 = 5,00 µm`, `z = 0` **gepinnt**; ein Kontaktwinkel wird nicht
zusätzlich vorgeschrieben und in Kombination mit dem Pinning ausdrücklich
abgelehnt.

## Was NICHT gerechnet wurde

Kein elektrisches Feld, kein Maxwell-Druck, keine Kopplung an den
Elektrostatiksolver, keine Betriebsspannung. Keine Emission, keine Raumladung.
Keine Strömung. Keine Zeitabhängigkeit und **keine Stabilitätsaussage**. Kein
Taylor-Kegel, kein Cone-Jet. Keine Schwerkraft — die Bond-Zahl, die das
rechtfertigt, steht unten. Die P2c-Vorratsgeometrie ist unverändert und wird
mechanisch nicht gelöst; der gefüllte Zulauf ist eine **Voraussetzung**, und der
Zustand stromauf geht ausschließlich über `delta_p_exit` ein.

---

## Kernzahlen

| Größe | Wert |
|---|---|
| Pinningradius `a = phi_2/2` | 5,000·10⁻⁶ m (aus `DeviceParameters`) |
| Oberflächenspannung γ | 4,520·10⁻² N/m — **Status `illustrative`** |
| Kapillardruckskala γ/a | 9040 Pa |
| Darstellbarer Bereich | \|Π\| ≤ 2, also \|Δp\| ≤ 18080 Pa |
| Bond-Zahl ρ g a²/γ | 6,94·10⁻⁶ |
| Größter Profilfehler gegen die Kugelkappe | 3,3·10⁻¹³·a |
| Größter relativer Fehler von Apexhöhe / Fläche / Volumen | 2,6·10⁻¹² / 4,5·10⁻¹⁴ / 1,3·10⁻¹² |
| Größtes Young-Laplace-Residuum (unabhängig gebildet) | 2,2·10⁻⁸ |
| Beobachtete Konvergenzordnung Profil / Residuum | 3,8 / 2,00 |
| Im Sweep abgelehnte Punkte | 22 von 421, größter gelöster \|Π\| = 1,99 |

---

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_geometry_and_bc.png` | Geometrie und Randbedingungen: Symmetrieachse, Bohrungsradius, gepinnte Kontaktlinie, Druckdefinition mit Vorzeichen, Richtung der äußeren Normale |
| `fig2_profiles.png` | Numerische Meniskusprofile für Π = −1,98 … +1,98 mit überlagerter analytischer Kugelkappe; rechts die punktweise Abweichung |
| `fig3_convergence.png` | Netzkonvergenz von Profilfehler, Apexhöhe, Fläche, Volumen und Young-Laplace-Residuum über fünf Auflösungen |
| `fig4_apex_and_curvature.png` | Apexhöhe und Krümmung über Π einschließlich der Grenze des darstellbaren Bereichs |
| `fig5_liquid_example.png` | Gesondert gekennzeichnetes **Stoffbeispiel** EMI-BF4 |

Jede Abbildung nennt P3a, „kein elektrisches Feld“, „keine Emission“,
„statisches Kapillargleichgewicht“, den Stoffstatus, Geometrie- und
Druckparameter sowie Commit und Konfiguration.

## Datendateien

| Datei | Inhalt |
|---|---|
| `profiles.csv` | Knoten aller gerechneten Profile, dazu die analytische Kappe und die Differenzen `dz_m`, `dn_m` in voller doppelter Genauigkeit |
| `profile_summary.csv` | je Druck: Status, gewählte Intervallzahl, Schätzfehler, numerisch und analytisch Apexhöhe, Bogenlänge, Fläche, Volumen, Krümmungsspanne, Residuum |
| `residuals.csv` | Young-Laplace-Residuum entlang der Oberfläche, allein aus den Knotenkoordinaten gebildet |
| `convergence.csv` | Netzstudie über 16 … 256 Intervalle für vier Drücke |
| `sweep.csv` | Π von −2,10 bis +2,10 in Schritten von 0,01, mit Status je Punkt |
| `range.csv` | Verhalten an und jenseits der Grenze, mit Statustext und der Angabe, ob eine Form zurückgegeben wurde |
| `liquid.csv` | Stoffdatensatz mit Status, Fundstelle und Vorbehalt |
| `liquid_example.csv`, `liquid_example_profiles.csv` | das gekennzeichnete Stoffbeispiel |
| `capillary_parameters.csv` | Größen des Kapillarproblems, SI |
| `parameters.csv`, `regions.csv`, `boundaries.csv`, `features.csv` | die unveränderte P1-Gerätegeometrie |
| `report.txt`, `meta.txt` | Bericht und Laufkennung |

---

## Stoffdaten — was belegt ist und was nicht

**Belegt ist die Stoffidentität.** `KunzeFynn-2024-12-10.pdf`, Abschnitt 2.3.2,
gedruckte Seite 28 (PDF-Seite 36): „EMI-BF4 and EMI-Im were also selected for
this project.“ Die Tabelle „List of Publications“, gedruckte Seite 30
(PDF-Seite 38), ordnet EMI-BF4 den Publikationen I–IV zu, also genau den geraden
10-µm-Kapillaren aus SU-8 bzw. IP-Q.

**Nicht belegt sind die Zahlenwerte.** Der Volltext enthält keine
Stoffwerttabelle für die Treibstoffe und keinen numerischen Wert für
Oberflächenspannung oder Dichte; Abschnitt 2.3.2 vergleicht die beiden
Flüssigkeiten ausschließlich qualitativ. γ und ρ sind deshalb **unverändert**
aus der quellenlosen Tabelle in `src/fluid.cpp` übernommen und tragen den Status
`illustrative`. Nichts wurde aus dem Gedächtnis ergänzt.

Die Prüfung des Lösers läuft davon unabhängig **dimensionslos** gegen die
analytische Kugelkappe. Viskosität, Leitfähigkeit und relative Permittivität
sind nur dokumentarisch vorgemerkt und gehen in P3a nicht ein.
