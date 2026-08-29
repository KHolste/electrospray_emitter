# P2c — Gerätegeometrie und Flüssigkeitsvorrat entkoppelt — 2026-08-29

Alle Daten und Abbildungen in diesem Ordner sind in diesem Lauf frisch erzeugt
und jede Abbildung nennt den tatsächlichen HEAD-Commit; die Herkunft steht in
`figures_provenance.txt` (`releasable=yes`, Arbeitsbaum sauber).

Reproduktion:

```sh
./build/es_reservoir examples/device_p1.cfg examples/reservoir_p2c.cfg \
    results/2026-08-29_p2c_reservoir_decoupling meta.commit=$(git rev-parse HEAD)
python python/plot_reservoir.py results/2026-08-29_p2c_reservoir_decoupling
```

Die beiden benutzten Konfigurationsdateien liegen als `device_p1.cfg` und
`reservoir_p2c.cfg` mit im Ordner.

---

## Was korrigiert wurde

Bis P2b verschob **ein einziger Wert** (`liquid_feed_z`) gleichzeitig

* die Länge der **leitfähigen Flüssigkeitssäule**,
* die Länge des **dielektrischen Rückteils**,
* und damit die **rückwärtige Gerätegeometrie**.

Die daraus berichtete „Konvergenz gegen die Lage der Zulaufgrenze“ war deshalb
**irreführend**: verändert wurde die **Geometrie der Hochspannungselektrode**.
Die Lage eines vollständig eingetauchten elektrischen Kontakts ist im Modell des
ideal leitfähigen Flüssigkeitskörpers dagegen irrelevant und wird gar nicht
geometrisch dargestellt.

`liquid_feed_z` existiert nicht mehr. An seine Stelle treten **vier getrennte
Parametergruppen**: vorderer Emitter, dielektrischer Grundkörper
(`base_plate_thickness`), fester Zulaufkanal (`feed_channel_length`,
`feed_channel_radius`) und Vorratsgeometrie (`plenum_*`).

Begründung, Audit und die vollständige Korrektur von 8.7:
[`docs/08_dielectric_model.md`](../../docs/08_dielectric_model.md), Abschnitt 8.9.

---

## Was gerechnet wurde — und was nicht

Genau ein Problem: ∇·(ε(x)∇φ) = 0, achsensymmetrisch, statisch, ρ_f = 0, mit

* der **gesamten zusammenhängenden Flüssigkeit** — Bohrung, Zulaufkanal und
  Plenum — als **einem** idealen Äquipotentialleiter auf `V_emitter`;
* dem gedruckten **Emitterkörper und Grundkörper als Dielektrikum**;
* dem **Vorratskörper als Dielektrikum** (im Gerät PEEK);
* dem **Extraktor als Polymerträger mit metallisierter Fläche**;
* einer **asymptotischen Monopolbedingung** am offenen Fernrand.

**Nicht** modelliert: Oberflächenspannung, Meniskusberechnung, Strömung,
Emission, Raumladung, endliche Leitfähigkeit der Flüssigkeit,
Zeitabhängigkeit, 3D. Die Fläche bei z = 0 ist die *anfängliche ebene
Flüssigkeitsoberfläche*. Schwerkraft und Benetzung im Vorratsraum sind nicht
modelliert; der Füllstand ist eine ebene Fläche und ein Parameter.

**Es wurde keine leitfähige Halterung, keine rückwärtige Metallscheibe und keine
Basisplatte auf Emitterpotential eingeführt.** `emitter_back_length ≠ 0` bleibt
abgelehnt.

---

## Vier Dinge, die dieser Lauf ausdrücklich NICHT belegt

**1. Das lokale Feld erreicht die Toleranz nur am Meniskus, nicht im Maximum.**
Die vor der Messung festgelegte Grenze von 10⁻³ (`es::reservoir_convergence`)
ist **nicht gelockert** worden. Beim letzten Vergrößerungsschritt (Volumen
×3,1) ändert sich

* E_z auf der Achse unmittelbar über der Oberfläche und |E| bei zwei
  Bohrungsradien um **1,1·10⁻⁴** — *innerhalb* der Grenze und unter dem
  Diskretisierungsfehler von 6,5·10⁻⁴;
* das **Maximum über alle kantenfernen Sondenpunkte** um **2,0·10⁻³** —
  *außerhalb* der Grenze, getrieben vom Punkt `axis_gap_three_quarter`.

Das Urteil ist der zweite Wert: **nicht auskonvergiert**. Es wurde keine
künstliche Elektrode eingeführt und keine „Referenzgröße" stillschweigend
ausgewählt. Siehe `reservoir_convergence.csv`, das die Änderung *je Sondenpunkt*
mitschreibt.

**2. Die Vorratsabmessungen sind eine Ersatzgeometrie.** Die Dissertation zeigt
in Abb. A.5 ein **nicht** rotationssymmetrisches Reservoir von etwa 25 mm Höhe
sowie 10 mm Außen- und 8 mm Innenbreite. Diese Angaben sind **ausschließlich zur
Wahl der Größenordnung** benutzt worden. Das achsensymmetrische Plenum bildet
Abb. A.5 nicht nach.

**3. `base_plate_thickness`, `feed_channel_length` und `wall_thickness` sind
VORLÄUFIGE Beispielwerte**, keine belegten Abmessungen. Sie sind in jeder
Variante identisch und werden nirgends als Kunze-Geometrie ausgegeben.

**4. Der Wert für den Vorratskörper ist ein Stellvertreter.** Im Gerät ist er
PEEK; `peek` ist in der Materialbibliothek registriert, trägt aber **absichtlich
keine Zahl**, weil kein geprüftes statisches ε_r vorliegt. Gerechnet wird mit dem
Emitterharz, und `materials.csv` sagt das im Vorbehalt. Ein Wert wird über
`material.peek.relative_permittivity` **und** `material.peek.source` ohne
Codeänderung nachgereicht.

Zusätzlich unverändert gültig: Austrittskante, äußere Stirnkante und beide
Aperturkanten sind **unverrundet**; dort wird kein Spitzenfeld berichtet.

---

## Nachweis der festen Frontgeometrie

Verglichen wird **bitweise** auf den Knotenkoordinaten und **exakt** auf dem
Material jeder Zelle, für alle Gitterzeilen von der Rückfläche des Grundkörpers
bis zum oberen Domänenrand und alle Radien bis zum Extraktoraußenrand
(`front_identity.csv`):

| Größe | Wert |
|---|---|
| verglichene Zeilen je Variante | 260 |
| verglichene Zellen je Variante | 45 066 |
| größte Knotenabweichung radial | **0 m** |
| größte Knotenabweichung axial | **0 m** |
| Materialabweichungen | **0** |

Das gilt auch für die abgeschnittene Säule, die als Diagnose mitläuft.

---

## Abbildungen

| Datei | Inhalt |
|---|---|
| `fig1_reservoir_geometries.png` | Maßstäbliche Geometrie: die feste Frontgeometrie in Mikrometern, daneben das abgeschnittene Säulenmodell und alle Plenumgrößen in Millimetern und untereinander in derselben Skala. |
| `fig2_full_mesh.png` | Das **vollständige** Volumennetz der tatsächlichen Rechendomäne (r bis 12 mm, z von −12 mm bis 9 mm), dazu vier Detailansichten: Vorratskörper, Emitter mit Kapillare, Extraktionsstrecke, Austrittskante. |
| `fig3_near_field.png` | Potential und Feldstärke im **identischen** Nahfeldausschnitt für alle Reservoirvarianten, mit **identischen** Farbskalen, gemeinsam über alle Varianten bestimmt. |
| `fig4_reservoir_convergence.png` | Konvergenz von lokalem Feld und globaler Ladung gegen die Vorratsgröße, mit der vorab festgelegten Toleranz und dem Diskretisierungsfehler als Vergleichslinien. |

## Daten

| Datei | Inhalt |
|---|---|
| `report.txt` | vollständiger Bericht des Laufs, einschließlich des offenen Befunds |
| `meta.txt` | Kennzahlen des Laufs, maschinenlesbar |
| `variants.csv` | alle verglichenen Varianten mit Geometrie, lokalen Feldern und Ladungen |
| `reservoir_convergence.csv` | Änderung je Vergrößerungsschritt, **inklusive der Änderung je einzelnem Sondenpunkt** |
| `front_identity.csv` | der bitweise Nachweis oben |
| `convergence_mesh.csv` | fünf Netzstufen bei festem Plenum — trennt Diskretisierungsfehler von Vorratswirkung |
| `variant_outlines.csv` | maßstäbliche Meridianumrisse aller Varianten |
| `reference_surface_fields.csv` | einseitiges E_z über der ebenen Flüssigkeitsoberfläche, je Variante |
| `field_near_*.csv` | abgetastete Felder im identischen Nahfeldfenster, je Variante |
| `field_window.csv` | die Grenzen dieses Fensters |
| `parameters.csv` | selbstbeschreibender Parametersatz der Referenzvariante |
| `mesh_nodes.csv`, `mesh_cells.csv`, `mesh_checks.csv` | Volumennetz der Referenzvariante und seine Prüfliste (11/11) |
| `node_roles.csv`, `boundary_audit.csv` | jeder festgehaltene Knoten mit seiner Rolle; das strukturelle Audit |
| `probes.csv`, `device_outline.csv`, `reference_surface_field.csv` | Auswertepunkte und Umrisse der Referenzvariante |
| `materials.csv`, `materials_library.csv` | Materialbelegung mit Wert, Status, Quelle, Bedingung, Vorbehalt |

---

## Kennzahlen des Referenzlaufs (`plenum_c`, Netzstufe 3)

| Größe | Wert |
|---|---|
| Knoten | 107 460 (270 × 398) |
| Netzprüfungen | 11 / 11 bestanden |
| E_z(r = 0, z = 0⁺) | 3,4019·10⁷ V/m |
| \|E\| bei zwei Bohrungsradien | 1,5193·10⁷ V/m |
| Q der Flüssigkeit / Gegenladung des Extraktors | 1,1206·10⁻⁹ C / −4,1092·10⁻¹⁰ C |
| festgehaltene Knoten ohne Flüssigkeitskontakt | **0** |
| festgehaltene Knoten auf benannten Polymerflächen | **0** |
| Frontgeometrie über alle Varianten | **bitgenau identisch** |
| lokales Feld am Meniskus, letzte Vergrößerungsstufe | 1,1·10⁻⁴ — Grenze eingehalten |
| Maximum über alle Sondenpunkte, letzte Stufe | 2,0·10⁻³ — **Grenze NICHT eingehalten** |

Zum Vergleich das überholte abgeschnittene Säulenmodell bei sonst gleicher
Frontgeometrie: E_z(0, 0⁺) = 7,6563·10⁷ V/m, Q = 6,94·10⁻¹² C. Der Unterschied
zum Plenum ist ein **Faktor 2,25** im lokalen Feld — nicht eine Konvergenzfrage,
sondern der Unterschied zwischen einem Leiter von 1,6·10⁴ µm³ und einem Vorrat
von Kubikmillimetern.

Die Summe der Ladungen ist nicht null: bei φ → 0 im Unendlichen trägt das System
Nettoladung, weil es keine Rückführelektrode gibt. Die globalen Ladungsgrößen
bleiben deshalb vorratsabhängig, und auf sie wird keine Toleranz gelegt.

---

## Wovon die absoluten Zahlen abhängen

Der Vernetzer verlangt, dass das Plenum **radial außerhalb** des modellierten
Extraktors liegt — nur dann kann seine Größe keinen einzigen Nahfeldknoten
verschieben. `extractor_outer_radius = 2 mm` ist selbst ein **Beispielwert**;
mit dieser Wahl ragt der Vorrat radial über die Elektrode hinaus und wird von
ihr nicht abgeschirmt. Ein **belegter Extraktoraußenradius** und **belegte
Vorratsabmessungen** sind damit die nächsten fehlenden Eingaben — nicht eine
feinere Rechnung.
