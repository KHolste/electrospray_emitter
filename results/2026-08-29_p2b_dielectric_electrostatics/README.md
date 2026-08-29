# P2b — dielektrische achsensymmetrische Elektrostatik — 2026-08-29

Alle Daten und Abbildungen in diesem Ordner sind in diesem Lauf frisch erzeugt
und jede Abbildung nennt den tatsächlichen HEAD-Commit; die Herkunft steht in
`figures_provenance.txt`.

Reproduktion:

```sh
./build/es_dielectric examples/device_p1.cfg examples/dielectric_p2b.cfg \
    results/2026-08-29_p2b_dielectric_electrostatics meta.commit=$(git rev-parse HEAD)
python python/plot_dielectric.py results/2026-08-29_p2b_dielectric_electrostatics
```

Die beiden benutzten Konfigurationsdateien liegen als `device_p1.cfg` und
`dielectric_p2b.cfg` mit im Ordner.

---

## Was gerechnet wurde — und was nicht

Genau ein Problem: ∇·(ε(x)∇φ) = 0, achsensymmetrisch, statisch, ρ_f = 0, mit

* der **ionischen Flüssigkeit** als idealem Äquipotentialleiter auf
  `V_emitter`, abgeschnitten an der `liquid_feed_boundary`;
* dem **3D-gedruckten Emitterkörper als Dielektrikum** — er ist *keine*
  Elektrode, keine seiner Flächen liegt auf einem Potential;
* dem **Extraktor als Polymerträger mit metallisierter Fläche**; nur die
  Metallisierung liegt auf `V_extractor`;
* einer **asymptotischen Monopolbedingung** am offenen Fernrand.

**Nicht** modelliert: Raumladung, Emission, Meniskusbewegung, Strömung, endliche
Leitfähigkeit der Flüssigkeit, Zeitabhängigkeit, 3D. Die Fläche bei z = 0 ist
die *anfängliche ebene Flüssigkeitsoberfläche*, kein berechneter Meniskus.

**Dieser Ordner ersetzt `../2026-08-29_p2a_vacuum_electrostatics/`**, in dem der
Emitterkörper als Metall behandelt wurde. Audit und Begründung:
[`docs/08_dielectric_model.md`](../../docs/08_dielectric_model.md).

---

## Zwei Dinge, die dieser Lauf ausdrücklich NICHT belegt

> **KORREKTUR (P2c, siehe `../2026-08-29_p2c_reservoir_decoupling/` und
> `docs/08_dielectric_model.md`, 8.9).** Punkt 1 unten ist in seinen Zahlen
> richtig, in seiner Deutung nicht. `liquid_feed_z` verschob nicht eine
> Randbedingung, sondern gleichzeitig die Länge der leitfähigen
> Flüssigkeitssäule, die Länge des dielektrischen Rückteils und die rückwärtige
> Gerätegeometrie — variiert wurde also die **Geometrie der
> Hochspannungselektrode**. Die dort genannte Abhilfe „Basisplatte auf
> Emitterpotential“ ist zurückgenommen; der fehlende Modellteil war der
> **Flüssigkeitsvorrat als dielektrisch umschlossener Körper**.
>
> Dieser Ordner ist damit die **Diagnose des überholten Säulenmodells**. Er wird
> nicht gelöscht, weil die Messung stimmt und der Vergleich gebraucht wird.
> Hinweis: der Parameter heißt im heutigen Code `base_plate_thickness`, und
> `python/plot_dielectric.py` liest die CSVs dieses Ordners nicht mehr.

**1. Die Lage der Zulaufgrenze ist nicht auskonvergiert.** Verlangt war der
Nachweis, dass eine Rückverlagerung der `liquid_feed_boundary` das Feld am
Meniskus nicht mehr ändert. Sie ändert es: eine Verdopplung der modellierten
Säulenlänge von 400 µm auf 800 µm verschiebt das Potential um 2,2·10⁻² der
angelegten Spannweite und das Feld um 1,4·10⁻¹ relativ — die *vor* der Messung
festgelegten Grenzen von 10⁻³ werden um mehr als eine Größenordnung verfehlt.
Ursache ist die Selbstkapazität der Säule: ein Leiter, den man verlängert, hält
mehr Ladung. Eine geerdete Hülle bei 25 mm ändert daran nichts (< 0,1 %), es
liegt also nicht am offenen Rand. `liquid_feed_z = −200 µm` ist deshalb ein
**BEISPIELWERT**, und **jede Zahl hier gilt für diesen Wert**. Siehe
`convergence_feed.csv` und `docs/08_dielectric_model.md`, 8.7.

**2. Der SU-8-Wert ist vorläufig.** ε_r = 3,3 ist aus 1-GHz-Angaben
(Herstellerdatenblätter und Literatur, 3,24…4,1) gewählt; ein statischer Wert
für unseren Aushärtungszustand existiert nicht. Aus einer Rechnung damit folgt
**keine Validierungsaussage**. Berichtet wird die Sensitivität über 2,8…4,5:
|E| an zwei Bohrungsradien bewegt sich darin um rund 1 %, die Emitterladung um
rund 7 %. Siehe `materials.csv` und `sensitivity_permittivity.csv`.

Zusätzlich unverändert gültig: Austrittskante, äußere Stirnkante und beide
Aperturkanten sind **unverrundet**. Das Feld einer scharfen Kante divergiert und
folgt der Elementgröße; es wird dort kein Spitzenfeld berichtet.

---

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_material_regions.png` | Materialgebiete und Randbedingungen. Die Randbedingungspunkte sind aus `node_roles.csv` gezeichnet, also aus dem Zustand, den der Löser tatsächlich aufgebaut hat. |
| `fig2_volume_mesh.png` | Automatisches Volumennetz: Übersicht, Kapillare mit Kegelflanke, Austrittskante. |
| `fig3_potential_field.png` | Potential und Feldstärke im SU-8-Referenzfall, je Übersicht und Spitzendetail. |
| `fig4_convergence.png` | Netzkonvergenz, Lage der Zulaufgrenze, ε_r-Sensitivität, Fernrandbehandlung. |

## Daten

| Datei | Inhalt |
|---|---|
| `report.txt` | vollständiger Bericht des Laufs, einschließlich des Befunds zur Zulaufgrenze |
| `meta.txt` | Kennzahlen des Laufs, maschinenlesbar |
| `parameters.csv` | selbstbeschreibender Parametersatz |
| `materials.csv`, `materials_library.csv` | Materialbelegung und -bibliothek mit Wert, Status, Quelle, Bedingung, Vorbehalt |
| `device_outline.csv` | geschlossene Meridianumrisse der Gebiete und benannte Randkurven |
| `mesh_nodes.csv`, `mesh_cells.csv`, `mesh_checks.csv` | Volumennetz und seine Prüfliste |
| `node_roles.csv`, `boundary_audit.csv` | jeder festgehaltene Knoten mit seiner Rolle; das Audit gegen versehentlich leitende Polymerflächen |
| `probes.csv` | feste Auswertepunkte mit Abstand zur nächsten unverrundeten Kante |
| `reference_surface_field.csv` | einseitiges E_z unmittelbar über der ebenen Flüssigkeitsreferenz |
| `convergence_mesh.csv` | fünf Netzstufen |
| `convergence_feed.csv` | vier Lagen der Zulaufgrenze — **der Befund oben** |
| `sensitivity_permittivity.csv` | neun Permittivitäten von 1 bis 4,5 |
| `farfield_study.csv` | vier Domänengrößen, asymptotisch gegen geerdet |
| `fem_vs_bem.csv` | Querprüfung gegen die unabhängige BEM bei ε_r = 1 |
| `linearity_polarity.csv` | Linearität und Polaritätsumkehr |
| `field_uebersicht.csv`, `field_spitze.csv` | abgetastete Felder für die Abbildungen |

---

## Kennzahlen des Referenzlaufs

| Größe | Wert |
|---|---|
| Netzstufe / Knoten | 3 / 74 807 |
| Residuum der freien Gleichungen | 1,7·10⁻²⁴ C |
| Q_emitter / Q_extractor / Summe | 6,948·10⁻¹² / −5,671·10⁻¹² / 1,277·10⁻¹² C |
| Netzkonvergenz, letzte Verfeinerung | Δφ = 1,1·10⁻⁴ der Spannweite, Δ\|E\| = 1,7·10⁻³ |
| FEM gegen BEM bei ε_r = 1 | Δφ = 2,8·10⁻⁴ der Spannweite, Δ\|E\| = 4,4·10⁻³ |
| D_n-Sprung an der Dielektrikumsgrenze | 1,5·10⁻² (fällt mit der Verfeinerung) |
| Fernrand bei 12 mm, asymptotisch gegen geerdet | 1,1·10⁻⁴ der Spannweite |
| festgehaltene Knoten auf Polymerflächen | **0** |
| Zulaufgrenze | **nicht konvergiert**, siehe oben |

Die Summe der Ladungen ist nicht null: bei φ → 0 im Unendlichen trägt das System
Nettoladung, weil es keine Rückführelektrode gibt.
