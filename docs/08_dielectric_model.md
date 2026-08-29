# 8. P2b — dielektrisches Modell des kapillaren Kunze-Emitters

**Stand: 2026-08-29. Status: gebaut und geprüft. Ein Befund offen, benannt.**

Dieses Dokument hält den **korrigierten physikalischen Modellvertrag** fest, das
Audit der Annahmen aus P2a, das gewählte numerische Verfahren, den verwendeten
SU-8-Wert samt Einschränkung und den einen Punkt, an dem die geforderte
Konvergenz nicht erreicht wird.

Ergebnisse:
[`results/2026-08-29_p2b_dielectric_electrostatics/`](../results/2026-08-29_p2b_dielectric_electrostatics/).

---

## 8.1 Audit: was in P2a falsch war und wo es stand

P2a löste die Elektrostatik so, als wäre der Emitter aus Metall. Das ist für
dieses Gerät falsch: der Emitter ist ein 3D-gedruckter, **nichtleitender**
Photopolymer, auf Hochspannung liegt die ionische Flüssigkeit beziehungsweise
ihr metallischer Kontakt.

| Stelle | falsche Annahme in P2a | Korrektur in P2b |
|---|---|---|
| `device_geometry.hpp`, `Region::EmitterSolid` | ein Festkörper ohne Materialeigenschaft, in P2a implizit Leiter | eigenes Materialmodell, `materials.hpp`; der Emitterkörper trägt `relative_permittivity` und **keine** Randbedingung |
| `vacuum_bem.hpp`, `vacuum_bem_mesh()` | akzeptiert alle Elemente, die Vakuum von Emittermetall, Extraktor **oder** Flüssigkeit trennen, und legt sie auf **einen** Leiter `Tag::Emitter` | in P2b ist nur die Flüssigkeit Leiter; der Polymerkörper ist Feldgebiet. Die BEM bleibt unverändert und dient als unabhängige Vergleichsrechnung für ε_r = 1 |
| `DeviceParameters::emitter_back_length` und `BoundaryId::NumericalEmitterBackClosure` | rückwärtige **leitende** Abschlussscheibe, damit die Einfachschicht auf einem geschlossenen Leiter sitzt | ersatzlos. `build_volume_mesh()` **lehnt** einen Wert ungleich 0 ab; die Flüssigkeitssäule endet an der Zulaufgrenze |
| `extractor_surface.*` (Apertur, Vorder-, Rückfläche, Rand) | **alle** Flächen der Elektrode auf `V_extractor` | nur die **metallisierte** Fläche. `extractor.metallisation` = `front_only` / `front_and_aperture` (Referenz) / `all_surfaces` |
| `free_surface_reference` bei z = 0 | Perfect-Conductor-Referenzfläche des Emitterleiters | unverändert richtig, aber jetzt aus dem richtigen Grund: sie ist die Oberfläche der **Flüssigkeit**, nicht des Emitters |
| Rückwärtiger Modellschnitt | in P2a durch die leitende Scheibe geschlossen | benannte `liquid_feed_boundary`. Dirichlet **ausschließlich** auf dem Flüssigkeitsquerschnitt r ≤ φ₂/2; der Rest der Ebene ist die Rückfläche des Polymers und keine Elektrode |
| `domain.radius = 3 mm` | für die BEM belanglos (freiraum-Greensfunktion, kein Trunkierungsrand) | für eine Volumenrechnung **viel zu klein**: dort unterscheiden sich asymptotische und geerdete Fernrandbedingung noch um 13 % der angelegten Spannweite. P2b rechnet mit 12 mm |

**Alle P2a-Zahlen, in die der Emitterkörper eingeht, sind damit überholt** —
nicht ungenauer, sondern für ein anderes Gerät gültig. Der Ergebnisordner
`results/2026-08-29_p2a_vacuum_electrostatics/` ist entsprechend gekennzeichnet.
Was aus P2a bestehen bleibt: die BEM selbst, ihre analytischen Prüfungen, und
die Erkenntnis, dass die Länge eines energetisierten Leiters eine Abmessung ist
und kein Konvergenzparameter — dieselbe Erkenntnis trifft P2b an der
Zulaufgrenze wieder (Abschnitt 8.7).

---

## 8.2 Korrigierter Modellvertrag

Gelöst wird, achsensymmetrisch in der Meridianhalbebene,

$$\nabla\cdot\big(\varepsilon(\mathbf{x})\,\nabla\varphi\big) = 0,
\qquad \rho_f = 0 .$$

| Gebiet / Fläche | Behandlung |
|---|---|
| Vakuum | ε_r = 1 |
| Emitterkörper | **Dielektrikum**, ε_r aus `materials.hpp`. Keine freie Ladung, keine Randbedingung, keine Fläche auf Potential |
| ionische Flüssigkeit | **idealer Äquipotentialleiter** auf `V_emitter`. Ihr gesamter Abschluss — Bohrungswand, ebene Referenzfläche bei z = 0 und Zulaufquerschnitt — trägt die Dirichlet-Bedingung |
| Extraktorträger | **Dielektrikum**, ε_r aus `materials.hpp` |
| Metallisierung | idealer Leiter auf `V_extractor`, **ohne Dicke** |
| `liquid_feed_boundary` | φ = `V_emitter` **nur** auf r ≤ φ₂/2 bei z = `liquid_feed_z` |
| Symmetrieachse r = 0 | natürliche Bedingung; das Gewicht 2πr verschwindet dort von selbst |
| Fernrand | asymptotische Monopolbedingung (Abschnitt 8.4) |

Ausdrücklich **nicht** modelliert: Raumladung, Emission, Meniskusbewegung,
Strömung, endliche Leitfähigkeit der Flüssigkeit, Zeitabhängigkeit, 3D.

### Warum die Flüssigkeit hier ein idealer Leiter sein darf — und wann nicht mehr

Die Ladungsrelaxationszeit ε/σ einer ionischen Flüssigkeit (σ ≈ 1 S/m,
ε_r ≈ 10) liegt bei etwa 10⁻¹⁰ s und ist gegen jeden in P2b beschriebenen
Vorgang vernachlässigbar. Sie hört auf zu gelten, **sobald Strom durch den
Meniskus fließt**: dann bestimmen der ohmsche Abfall entlang des Kegels und die
endliche Nachlieferung von Oberflächenladung das Apexfeld — das
Leaky-Dielectric-Regime. P3/P4 müssen die Dirichlet-Bedingung auf der
Flüssigkeit durch ein Leitungsmodell ersetzen. Eine P2b-Feldaussage ist eine
Aussage über **vollständige Abschirmung**, keine Vorhersage für den emittierenden
Betrieb.

### Die Zulaufgrenze ist auch der hydraulische Zulauf

`liquid_feed_boundary` ist semantisch bereits als **hydraulischer Zulauf**
vorgesehen: dort wird in P3 Reservoirdruck, Volumenstrom oder eine hydraulische
Impedanz vorgegeben (`feed.mode` aus 04_geometry_model.md, 4.2). In diesem Lauf
ist davon **nichts** implementiert; die Grenze ist benannt, lokalisiert und
parametriert, damit das Hinzufügen später eine neue Bedingung auf einer
vorhandenen Entität ist und kein Umbau.

---

## 8.3 Material- und Gebietsmodell

`include/es/materials.hpp` führt Material = **Wert + Herkunft + Status**.
Registriert sind `vacuum`, `su8`, `ip-q`, `ipx-q`, `emibf4`, `metal`.

### Der verwendete SU-8-Wert

**Was verfügbar ist.** Die Kunze-Dissertation nennt SU-8, IP-Q und IPx-Q als
verwendete Harze und die kleine Geometriegeneration (Kapillaren um 8–10 µm) als
in SU-8 gefertigt, liefert aber **keine** dielektrischen Daten — der Volltext
enthält weder „permittivity“ noch „dielectric constant“. Belastbar sind:

| Quelle | Wert | Bedingungen |
|---|---|---|
| MicroChem/Kayaku, Datenblatt SU-8 3000, Eigenschaftstabelle | ε_r = 3,28; Durchgangswiderstand 7,8·10¹⁴ Ω·cm | **1 GHz** |
| MicroChem/Kayaku, Datenblatt SU-8 2000 | ε_r = 4,1 | **1 GHz**, 50 % r.F. |
| Ghalichechian & Sertel, *IEEE Antennas Wirel. Propag. Lett.* **14** (2015) 723, doi:10.1109/LAWP.2014.2377695 | ε_r = 3,24 / 3,23 / 2,92 | 1 GHz / 200 GHz / 1 THz, vollständig vernetzter 430-µm-Film |
| Melai, Salm, Smits, Visschers, Schmitz, *J. Micromech. Microeng.* **19** (2009) 065012, doi:10.1088/0960-1317/19/6/065012 | Durchschlagfestigkeit 4,4 MV/cm, Leckstrom thermionisch | — |

Die letzte Quelle trägt die Prämisse des ganzen Modells: SU-8 ist bei den
auftretenden Feldern ein Isolator.

**Was nicht verfügbar ist.** Kein Wert im **statischen** Grenzfall für ein
zweiphotonengedrucktes, hartgebackenes Teil im Vakuum. Alle Zahlen oben sind bei
1 GHz oder darüber gemessen. Für ein passives Dielektrikum ist ε′(ω) nicht
steigend in ω, der statische Wert liegt also **mindestens** beim 1-GHz-Wert.

**Die Wahl.** Nominal **ε_r = 3,3**, ausdrücklich **VORLÄUFIG**
(`MaterialStatus::Provisional`). Sensitivitätsbereich **2,8 … 4,5**: unten der
niedrigste berichtete Messwert für SU-8, oben oberhalb der höchsten
Herstellerangabe, um Formulierungsstreuung und den Anstieg zu tiefen Frequenzen
zu decken. Wasseraufnahme erhöht ε_r, das Gerät läuft im Vakuum — der trockene
Rand ist der relevantere, aber „relevanter“ ist auch keine Messung, deshalb wird
der ganze Bereich berichtet.

**Aus einer Rechnung mit diesem Wert darf keine Validierungsaussage abgeleitet
werden.** Gemessene Wirkung über den ganzen Bereich: |E| an zwei
Bohrungsradien ändert sich um **rund 1 %**, die Emitterladung um **rund 7 %**.
Die Feldgrößen an der Referenzfläche sind also gegen die Unsicherheit von ε_r
weitgehend robust; die Ladung ist es nicht.

Es wird ausdrücklich **nicht** n² aus einem optischen Brechungsindex verwendet:
optisches n erfasst nur die elektronische Polarisierbarkeit und verfehlt den
Orientierungsanteil, der den statischen Wert eines Epoxids dominiert — für SU-8
gäbe das etwa 2,6 und läge um ein Viertel daneben.

### IP-Q und IPx-Q

Registriert, mit `MaterialStatus::Unknown` und **ohne Zahl**. Sie abzufragen
bricht mit einer Erklärung ab. Ein Wert wird über
`material.<name>.relative_permittivity` **und** `material.<name>.source` in der
Konfiguration nachgereicht; ohne Quellenangabe wird der Wert abgelehnt. Ein
Codeeingriff ist dafür nicht nötig — geprüft in
`tests/test_dielectric_device.cpp`.

---

## 8.4 Numerisches Verfahren

**Achsensymmetrische Q1-Finite-Elemente** auf einem blockstrukturierten
Viereckgitter, `include/es/axisym_fem.hpp`.

* Schwache Form $\int \varepsilon\,\nabla\varphi\cdot\nabla v\;2\pi r\,\mathrm{d}r\,\mathrm{d}z$,
  isoparametrisch, 2×2-Gauß. Die 2πr-Gewichtung ist gegen einen
  Koaxialkondensator geprüft, dessen Kapazität eine geschlossene Funktion der
  Radien ist und die keine ebene Formulierung reproduziert.
* **Potentialstetigkeit** an Materialgrenzen gilt durch Konstruktion: die
  Grenzfläche trägt **einen** Knotenfreiheitsgrad. Gemessen wird sie trotzdem
  (exakt 0).
* **Stetigkeit der normalen Flussdichte** ohne freie Grenzflächenladung ist die
  Prüfung mit Inhalt: sie gilt nur, wenn die Assemblierung stimmt. Der
  gemessene einseitige Sprung ist ein O(h)-Diskretisierungsfehler und fällt
  monoton mit der Verfeinerung.
* **Symmetrie bei r = 0**: natürliche Bedingung, das Gewicht 2πr verschwindet.
  Kein 1/r-Term, keine Sonderbehandlung.
* **Löser**: symmetrische Band-Cholesky-Zerlegung mit interner Nummerierung in
  der schmalen Richtung. Direkt, deterministisch, keine Iterationstoleranz.
  Gegen eine dichte LU-Zerlegung desselben Systems geprüft.
* **Ladungen** aus den **Knotenreaktionen** des unreduzierten Operators,
  $Q_S=\sum_{i\in S}(K\varphi)_i$ — exakt und deutlich genauer, als φ zu
  differenzieren und das Ergebnis zu integrieren.
* **Feldauswertung** über knotenweise Rückgewinnung aus den
  Zellmittelpunktsgradienten (dem superkonvergenten Q1-Punkt), gewichtet mit dem
  Rotationsvolumen und **nur über Zellen gleicher Permittivität** — ein roher
  Q1-Gradient ist erster Ordnung und springt von Zelle zu Zelle, eine Netzstudie
  darauf misst das Zittern statt des Fehlers. Für den **einseitigen** Wert an
  einer Grenzfläche wird weiterhin `field_in_cell()` benutzt; die beiden
  beantworten verschiedene Fragen.

### Fernrandbehandlung

Das Gerät ist nicht eingehaust, und das System trägt **Nettoladung**: der
Emitter liegt auf Spannung, der Extraktor auf Masse, und nichts führt die
Emitterladung zurück. Verwendet wird deshalb die **asymptotische
Monopolbedingung** als Robin-Term,

$$\frac{\partial\varphi}{\partial n}
= -\,\varphi\,\frac{\mathbf{n}\cdot(\mathbf{x}-\mathbf{x}_0)}{|\mathbf{x}-\mathbf{x}_0|^{2}},$$

die für jedes beschränkte Ladungssystem die führende Ordnung exakt trifft; der
Rest ist Dipolordnung. α > 0 auf dem ganzen Rand einer Box, die x₀ enthält, das
System bleibt positiv definit.

Gemessen (`farfield_study.csv`), Potential in der Mitte der Extraktionsstrecke:

| Domänenradius | asymptotisch | geerdete Hülle | Differenz / Spannweite |
|---|---|---|---|
| 3 mm (P1-Box) | 112,72 V | 110,32 V | 1,6·10⁻³ |
| 12 mm | 112,432 V | 112,271 V | 1,1·10⁻⁴ |
| 24 mm | 112,434 V | 112,375 V | 3,9·10⁻⁵ |
| 48 mm | 112,430 V | 112,404 V | 1,7·10⁻⁵ |

Die Differenz der beiden Bedingungen ist ein direktes Maß des
Trunkierungsfehlers. Die asymptotische Lösung ist bei 12 mm ausgereizt: eine
Vervierfachung auf 48 mm bewegt sie um 1,1·10⁻⁵ der Spannweite. **12 mm ist
deshalb der Referenzwert.** Die 3-mm-Box aus P1 wäre für eine Volumenrechnung um
Größenordnungen zu klein gewesen; die BEM brauchte das nicht, weil ihr Kern
V → 0 im Unendlichen bereits trägt.

---

## 8.5 Volumenvernetzung — die in 07 offen gelassene Entscheidung

[07_mesher_decision.md](07_mesher_decision.md), 7.3, war bindend: **kein
allgemeiner Volumenvernetzer im Eigenbau, keine beiläufige schwere externe
Abhängigkeit.** Beides ist eingehalten.

Gebaut ist ein **blockstrukturierter, radial gewarpter** Vernetzer,
`include/es/volume_mesh.hpp`, der eine reine Funktion der Geräteparameter ist.
Die Meridianhalbebene dieses Geräts ist ein Stapel von z-Streifen aus radialen
Blöcken, und die einzige nicht achsenparallele Grenze ist die gerade
Kegelflanke. Ein Tensorproduktgitter in (r_ref, z), radial so verzerrt, dass
eine Gitterlinie **exakt** auf der Flanke liegt,

$$r(i,j) = W\big(r_\mathrm{ref}[i],\, z[j]\big),$$

stellt damit jede Materialgrenze exakt dar — ohne hängende Knoten, ohne
Treppenstufen, ohne Delaunay-Maschinerie. W ist stetig, stückweise linear,
streng monoton und hält r_Bohrung und den ersten Schlüsselradius außerhalb des
Emitters fest. Da ρ(z) stückweise linear in z ist und die Gitterzeilen die
Knickstellen sind, sind die gewarpten Gitterlinien innerhalb jeder Zelle gerade
— das Netz reproduziert deshalb **jedes Gebietsvolumen exakt** (Prüfung 5).

Größenfunktion wie beim Randvernetzer: Minimum über Kegel auf endlich vielen
Quellen, Lipschitz-beschränkt, Knoten bei gleichen Zuwächsen von ∫ds/h.
Konstanten in `namespace es::volume_mesher`. Kein Elementmaß, keine
Verfeinerungszone in der Benutzerschnittstelle; einziger Steuerwert ist
`mesh.reference_level`.

**Zwei Dinge, die beim Bau falsch waren und gemessen wurden:**

1. Die Stufe skalierte anfangs nur die **Quellgrößen** und h_max, nicht den
   Gradationskegel. Fern der Merkmale ist h ≈ G·d und bewegt sich dann mit der
   Stufe überhaupt nicht — die Elementgröße in der Extraktionsstrecke blieb über
   sechs Stufen bei etwa 60 µm, und die „Netzstudie“ maß das
   Interpolationszittern. Die Stufe skaliert jetzt das **fertige** Feld.
2. Die Integration von ∫dx/h lief über eine feste Zahl gleichmäßiger Stützstellen.
   Das Größenfeld umfasst vier Dekaden innerhalb eines Schlüsselintervalls; die
   ersten Elemente kamen um einen Faktor zwei zu groß heraus, genau an der
   Stirnebene. Die Stützstellen werden jetzt lokal verfeinert.

Was der Vernetzer **nicht** kann: adaptieren, a posteriori verfeinern, eine
bewegte freie Oberfläche führen. P3 braucht alle drei; die Entscheidung aus 7.3
ist dann richtig zu treffen. Dieser Vernetzer ist die Antwort für die statische
Geometrie und billig wegzuwerfen.

**Keine externe Abhängigkeit eingeführt.** Windows/MinGW- und WSL-Build,
reproduzierbare Einrichtung und Lizenz einer Netzbibliothek mussten deshalb
nicht geprüft werden. Das Projekt hat weiterhin außer CMake und einem
C++-Compiler keine Abhängigkeit.

---

## 8.6 Pflichtprüfungen

`tests/test_axisym_fem.cpp` (6 Gruppen) und `tests/test_dielectric_device.cpp`
(10 Gruppen). Beide laufen im Rahmen von `ctest`.

| # | geforderte Prüfung | wo | Ergebnis |
|---|---|---|---|
| 1 | geschichtetes Dielektrikum, analytisch | `test_axisym_fem` 1 und 2 | Koaxial: Potentialfehler 4,9·10⁻⁷; axial geschichtet exakt (Q1 reproduziert die stückweise lineare Lösung, 3·10⁻¹⁵) |
| 2 | 2πr-Gewichtung gegen bekannte achsensymmetrische Geometrie | `test_axisym_fem` 1, 2, 3 | Koaxialkapazität 6,3·10⁻⁶ relativ; die ebene Formel läge um mehr als Faktor 2 daneben; φ = z²−r²/2 fällt wie h² |
| 3 | Stetigkeit von φ und D_n an der Materialgrenze | `test_axisym_fem` 1/2, `test_dielectric_device` 3 | φ exakt stetig; D_n-Sprung koaxial 8,3·10⁻³ → 4,2·10⁻³ bei Verfeinerung (O(h)); am Gerät 5,0·10⁻² → 2,3·10⁻² über drei Stufen |
| 4 | Linearität und Polaritätsumkehr | `test_axisym_fem` 5, `test_dielectric_device` 4 | exakt bis auf Rundung (≤ 4·10⁻¹⁵); Superposition der beiden Einheitslösungen 3,6·10⁻¹⁵ |
| 5 | FEM gegen BEM für ε_r = 1, kantenfern | `test_dielectric_device` 5 | φ ≤ 2,8·10⁻⁴ der Spannweite, \|E\| ≤ 4,4·10⁻³ relativ, an sieben Punkten |
| 6 | Volumennetz-Konvergenz an kantenfernen Punkten | `test_dielectric_device` 6 | fünf Stufen; letzte Verfeinerung Δφ = 1,1·10⁻⁴ der Spannweite, Δ\|E\| = 1,7·10⁻³, ΔQ_E = 1,2·10⁻³ |
| 7 | Konvergenz gegen die Lage der `liquid_feed_boundary` | `test_dielectric_device` 7 | **NICHT ERREICHT** — siehe 8.7 |
| 8 | Sensitivität gegenüber ε_r | `test_dielectric_device` 8 | Q_E monoton steigend; \|E\| ändert sich über 2,8…4,5 um rund 1 % |
| 9 | keine Polymerfläche wird als Leiter behandelt | `test_dielectric_device` 2 | 0 festgehaltene Knoten auf allen benannten Polymerflächen; ein absichtlich verfälschter Rollenvektor wird gefunden |
| 10 | vollständiger Testlauf ohne Regression | `ctest` | 13/13 bestanden |

Zusätzlich: Determinismus des Vernetzers (bitgenau), Unabhängigkeit des Netzes
ab dem Kegelfuß von der Zulaufposition, Ablehnung der P2a-Abschlussscheibe,
Ablehnung eines Materials ohne Wert, Ablehnung der metallischen Referenz mit
ε_r ≠ 1, und der Nachweis, dass das korrigierte Modell **messbar andere**
Ergebnisse liefert als das überholte (Q_E um 20 % verschieden).

**Unverrundete Kanten bleiben ausgeschlossen.** Austrittskante, äußere
Stirnkante und beide Aperturkanten sind scharf; das Feld einer scharfen Kante
divergiert und folgt der Elementgröße. Kein Spitzenfeld wird dort berichtet, und
die kantennahen Zellen des Profils auf der Referenzfläche sind gezählt und
ausgeschlossen. Endliche Kantenradien sind ein P3-Parameter.

---

## 8.7 Offener Befund: die Zulaufgrenze konvergiert nicht

Gefordert war der Nachweis, dass eine weitere Rückverlagerung der Zulaufgrenze
Potential und Feld am Meniskus nicht mehr ändert. **Sie ändert sie.**

Grenzen, **vor** der Messung festgelegt (`es::feed_truncation`, dieselben wie in
P2a und aus denselben Gründen): eine Verdopplung der modellierten Säulenlänge
darf φ um weniger als 10⁻³ der Spannweite und |E| um weniger als 10⁻³ relativ
bewegen.

Gemessen (`convergence_feed.csv`, Netzstufe 2):

| \|z_feed\| | φ(Achse, halbe Strecke) | E_z(2 Bohrungsradien) | Q_emitter |
|---|---|---|---|
| 100 µm | 85,71 V | 3,495·10⁷ V/m | 4,44·10⁻¹² C |
| 200 µm | 112,76 V | 3,285·10⁷ V/m | 6,99·10⁻¹² C |
| 400 µm | 140,15 V | 3,124·10⁷ V/m | 1,141·10⁻¹¹ C |
| 800 µm | 163,03 V | 3,012·10⁷ V/m | 1,923·10⁻¹¹ C |

Letzte Verdopplung: **2,2·10⁻² der Spannweite** und **1,4·10⁻¹** relativ im
Feld — die Grenzen werden um mehr als eine Größenordnung verfehlt. **Die Grenze
wird nicht verschoben.**

**Ursache, und warum es kein numerischer Mangel ist.** Die Flüssigkeitssäule ist
ein Leiter auf `V_emitter`. Eine dünne Säule der Länge L und des Radius a hat
eine Selbstkapazität von etwa $2\pi\varepsilon_0 L/(\ln(2L/a)-1)$; eine längere
Säule trägt proportional mehr Ladung, und diese Ladung wird an der Spitze
gespürt. Die gemessene Emitterladung folgt dieser Formel auf ein über den ganzen
Bereich konstantes Verhältnis von 1,37…1,43. Die Änderungen klingen ab — etwa
wie 1/ln(L/a) —, aber nirgends schnell genug.

**Am offenen Rand liegt es nicht.** Dieselbe Studie in einer **geerdeten Hülle**
bei 25 mm liefert jeden Wert auf besser als 0,1 % gleich. Es ist also keine
Eigenschaft der Fernrandbehandlung, sondern des Geräts, wie es hier modelliert
ist.

**Es ist derselbe Mechanismus wie in P2a.** Dort war es `emitter_back_length`:
ein Leiter, den man verlängert, hält mehr Ladung, Punkt.

**Was es beheben würde.** Der Emitterhalter.
[04_geometry_model.md](04_geometry_model.md), 4.1, sieht bereits eine
**Basisplatte auf Emitterpotential** vor, aus der die verjüngte Struktur
herausragt. Ein echter Leiter mit angegebenen Abmessungen, der die Säule
abschließt, würde die Zulaufposition wieder zu einem numerischen Parameter
machen. Das ist eine **Geometrieentscheidung** und gehört in eine spätere Phase.
Die hintere Schnittebene *als solche* zur Elektrode zu erklären wäre keine
Geometrieentscheidung, sondern eine erfundene Scheibe, damit eine
Konvergenzstudie gelingt — und ist in diesem Auftrag ausdrücklich ausgeschlossen.

**Konsequenz.** `liquid_feed_z = −200 µm` ist ein **BEISPIELWERT**, keine
gemessene Abmessung — genau wie `extractor_outer_radius`. Jede in P2b berichtete
Zahl gilt für diesen Wert. **Nichts in P2b ist bezüglich der Zulauftrunkierung
konvergiert**, und der Ergebnisordner sagt das an jeder Stelle, an der eine Zahl
steht.

---

## 8.8 Was P2b nicht ist

* Keine Meniskuslösung. Die Fläche bei z = 0 heißt weiterhin **anfängliche ebene
  Flüssigkeitsoberfläche**, ist der geometrische Ausgangszustand und keine
  physikalische Meniskuslösung.
* Keine Emissionsaussage. Das Feld gilt für vollständige Abschirmung durch eine
  ideal leitende Flüssigkeit.
* Keine 3D-Validierung. 3D ist eine spätere ergänzende Stufe.
* Keine Validierung des SU-8-Werts. Der Nominalwert ist vorläufig; berichtet wird
  die Sensitivität über einen aus Quellen begründeten Bereich.
