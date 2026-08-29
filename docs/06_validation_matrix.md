# 6. Validierungsmatrix

**Was hier nicht als Validierung zählt:** interne Konsistenzprüfungen, die
Reproduktion einer im Code selbst einprogrammierten Skalierung, der Vergleich
zweier Varianten desselben Verfahrens, und Netzkonvergenz. Netzkonvergenz zeigt,
dass die Diskretisierung gegen *etwas* konvergiert, nicht dass es das Richtige
ist. Sie ist notwendige Voraussetzung, kein Nachweis.

Kategorien:

* **A — analytisch:** geschlossene Lösung, unabhängig vom Code.
* **B — numerisch:** publizierte Ergebnisse eines anderen Verfahrens oder einer
  anderen Arbeitsgruppe.
* **C — experimentell:** Messdaten.

Statusspalte: `offen` / `geprüft (Abweichung)` / `fehlgeschlagen`.

---

## A — Analytische Referenzen

| # | Zielgröße | Referenzlösung | Toleranz | Phase | Status |
|---|---|---|---|---|---|
| A1 | $K(m)$, $E(m)$ | exakte Werte, unabhängige Quadratur | 10⁻¹² | P0 | geprüft (10⁻¹⁴) |
| A2 | Kapazität isolierte Kugel | $4\pi\varepsilon_0R$ | 10⁻⁴ | P2 | geprüft (2,5·10⁻⁶) |
| A3 | Spitzenfeld Rotationsellipsoid | prolat-sphäroidale Lösung | 10⁻³ | P2 | geprüft (1,8·10⁻³, Aspekt 20) |
| A4 | Feldverteilung Rotationsellipsoid, elementweise | dito, $E_n(\eta)$ | 5·10⁻³ | P2 | geprüft (1,8·10⁻³) |
| A5 | Kapazität koaxiale Zylinder | $2\pi\varepsilon_0 L/\ln(b/a)$ | 10⁻³ | P2 | geprüft (6,3·10⁻⁶, P2b-FEM) |
| A5b | Koaxial mit **zwei Dielektrikumsschichten** | $2\pi L/\big[\ln(c/a)/\varepsilon_1+\ln(b/c)/\varepsilon_2\big]$ | 10⁻³ | P2b | geprüft (6,3·10⁻⁶) |
| A5c | Axial geschichtetes Dielektrikum, Q1 exakt | $\varepsilon_0\pi R^2/\big[c/\varepsilon_1+(h-c)/\varepsilon_2\big]$ | 10⁻¹² | P2b | geprüft (1,6·10⁻¹⁴) |
| A5d | Stetigkeit von $D_n$ an der Materialgrenze | $\varepsilon_1E_1=\varepsilon_2E_2$ | O(h), fallend | P2b | geprüft (8,3·10⁻³ → 4,2·10⁻³) |
| A5e | Achsensymmetrisch harmonisch, $z^2-r^2/2$ | exakte Lösung; prüft die 2πr-Gewichtung | O(h²) | P2b | geprüft (Ordnung 1,9) |
| A6 | Kugel im homogenen Feld | $E_\mathrm{max}=3E_0$ | 10⁻³ | P2 | offen |
| A7 | Randbedingungsresiduum | Dirichletdaten an Kollokationspunkten | 10⁻¹⁰ | P2 | geprüft (1,5·10⁻¹⁵) |
| A8 | Fernfeld-Monopol | $Q/4\pi\varepsilon_0R$ | 10⁻² | P2 | geprüft (10⁻³) |
| A8b | FEM gegen BEM, $\varepsilon_r=1$, kantenfern | unabhängiges Verfahren ohne Trunkierungsrand | 10⁻² | P2b | geprüft (2,8·10⁻⁴ in φ, 4,4·10⁻³ in \|E\|) |
| A9 | Meniskus ohne Feld | Kugelkappe $R=2\gamma/\Delta p$ | 10⁻⁵ | P3 | geprüft (4·10⁻⁶) |
| A10 | Hängender Tropfen ohne Feld | Bashforth-Adams-Tafeln | 10⁻³ | P3 | offen |
| A11 | Taylor-Kegelwinkel | 49,29°, Nullstelle von $P_{1/2}(\cos\theta)$ | 0,2° | P3 | geprüft (49,16°) — **siehe Anmerkung** |
| A12 | Rayleigh-Grenzladung | $q_R=8\pi\sqrt{\varepsilon_0\gamma R^3}$ | 10⁻¹² | P6 | geprüft (Formelidentität) |
| A13 | Ringpotential im Fernfeld | Punktladung $Q/4\pi\varepsilon_0d$ | 10⁻⁴ | P4 | geprüft (10⁻⁵) |
| A14 | Energieerhaltung Trajektorie | $E_\mathrm{kin}/q = \Delta\phi$ | 10⁻⁴ | P4 | geprüft (7·10⁻⁶) |
| A15 | **Child-Langmuir** | $J=\frac{4}{9}\varepsilon_0\sqrt{2q/m}\,V^{3/2}/d^2$ | 2·10⁻² | P4 | offen — **zentraler PIC-Test** |
| A16 | Schottky-Absenkung | $\sqrt{e^3E/4\pi\varepsilon_0}$, 1,20 eV bei 1 V/nm | 10⁻³ | P5 | geprüft (Formelidentität) |
| A17 | Hagen-Poiseuille-Impedanz | $Z_h=8\mu L/\pi a^4$ | 10⁻³ | P5 | offen |
| A18 | Ladungsrelaxation, ebene Grenzfläche | $\tau_e=\varepsilon_0\varepsilon_r/K$ | 10⁻³ | P5 | offen |

**Anmerkung zu A11.** Der Wert 49,16° zeigt, dass der Young-Laplace-Löser mit
Maxwell-Term das asymptotische Kegelverhalten reproduziert. Das ist ein
Diskretisierungstest, kein Nachweis der Emissionsphysik, und er liegt auf dem
**instabilen** Ast. Er ist hier als A11 geführt, weil die Referenz analytisch
und codeunabhängig ist — er darf aber nicht als Beleg für den Betriebspunkt
zitiert werden. Im bisherigen Haupt-README war er in dieser Rolle
überbewertet.

**Anmerkung zu A12, A16.** Diese beiden prüfen nur, dass eine Formel richtig
programmiert ist. Sie sind der schwächste Eintrag in dieser Tabelle und zählen
im Sinne der Eingangsbemerkung nur eingeschränkt.

---

## B — Numerische Referenzen aus der Literatur

| # | Zielgröße | Referenz | Was verglichen wird | Phase | Status |
|---|---|---|---|---|---|
| B1 | Stabilitätsgrenze elektrifizierter Tropfen | Wohlhuter & Basaran (1992) ✓ | Deformation und Stabilitätsgrenze über der Feldstärke, für dielektrische Tropfen bei festem Kontaktwinkel | P3 | offen |
| B2 | Meniskusfamilien bei Ionenverdampfung | Coffman et al. (2016) ✓ | Gleichgewichtsformen, die messbare Ladung abgeben, als Funktion der Meniskusgröße | P5 | offen |
| B3 | Emissionseigenschaften und Stabilität | Gallud & Lozano (2022) ✓ | Strom über Hintergrundfeld; minimale hydraulische Impedanz; Maximalstrom; Optimum bei 0,5–3 µm Meniskusradius; Stabilitätskriterium (elektrischer Druck > 2× Oberflächenspannungsdruck der gleich großen Kugel) | P3, P5 | offen |
| B4 | PIR-Skalierungen | Higuera (2008) ✓ | Viskositätsdominanz der Meniskusströmung; Kontrolle des Verdampfungsstroms durch die endliche Leitfähigkeit; Vernachlässigbarkeit der Raumladung an der Oberfläche | P5, P7 | offen |
| B5 | Cone-Jet-Struktur | Higuera (2003) | Strom und Flussrate des Taylor-Kegels aus der vollen numerischen Lösung | P6 | offen |
| B6 | Cone-Jet-Numerik | Herrada et al. (2012) | Kegel-Jet-Form und Übergangsbereich | P6 | offen |
| B7 | Kegelspitzen-Bildung | Collins et al. (2008) | Tip-Streaming, Bildung konischer Spitzen | P3 | offen |

Für B1–B4 gilt: die Vergleichsgrößen sind aus den Arbeiten zu entnehmen und als
digitalisierte Datensätze im Repositorium abzulegen, damit der Vergleich
reproduzierbar bleibt.

---

## C — Experimentelle Referenzen

| # | Zielgröße | Referenz | Phase | Status |
|---|---|---|---|---|
| C1 | I–U-Kennlinie im PIR, EMI-BF4 / EMI-Im | Romero-Sanz et al. (2003) | P5 | offen |
| C2 | Extern benetzte Emitter, Strom und Spektrum | Lozano & Martínez-Sánchez (2005) | P5 | offen |
| C3 | Strahlstruktur im Vakuum | Gamero-Castaño (2008) | P7 | offen |
| C4 | Poröse Emitter, I–U und Lebensdauer | Legge & Lozano (2011) | P5 | offen |
| C5 | TOF-Spektren, Speziesaufteilung Monomer/Dimer/Trimer | Literatur zur IL-Massenspektrometrie | P5 | offen |
| C6 | Kolloidtriebwerk im Flug | ST7-DRS / LISA Pathfinder | P8 | offen |
| C7 | **Eigene Messdaten** | zunächst nicht verfügbar; Import strukturell vorsehen | P8 | zurückgestellt |

C7 wäre der wertvollste Eintrag der Tabelle: nur damit ließe sich $\Delta G$ für
die tatsächlich verwendete Flüssigkeitscharge anpassen, statt einen
Literaturwert mit einem Streubereich von 1,0–1,4 eV zu übernehmen. Bis dahin
gilt: $\Delta G$ bleibt ein deklarierter Anpassungsparameter aus der Literatur,
und jede Stromangabe trägt die zugehörige Unsicherheit.

---

## Anforderungen an jede Netzkonvergenzprüfung

Nicht als Validierung, sondern als deren Voraussetzung. Zu überwachen sind
mindestens (Befund 9):

| Größe | Warum |
|---|---|
| Apexfeld $E_\mathrm{apex}$ | geht exponentiell in den Strom ein |
| Apexkrümmungsradius | bestimmt die Emissionsstruktur |
| Emissionsstrom $I$ | Zielgröße |
| Emittierende Fläche | quantisiert, siehe unten |
| Divergenzwinkel | Zielgröße der Strahlrechnung |
| Stabilitätsgrenze | darf nicht netzabhängig sein |

Die „99-%-Stromfläche" ist als Summe über diskrete Elemente quantisiert und
konvergiert nicht glatt (im Prototyp gemessen: 3,198 / 3,150 / 3,103 /
3,150 ·10⁻¹⁰ m² bei 61/81/121/161 Knoten). Sie ist durch ein stetiges Funktional
zu ersetzen, etwa den flächengewichteten Median der Stromdichteverteilung.

---

## Regeln für die Berichterstattung

1. Jede ausgegebene Zahl trägt die Modellstufe, aus der sie stammt.
2. Empirische Korrelationen tragen `empirical = true` und ihren
   Gültigkeitsbereich, in der Konsolenausgabe wie in jeder Datei.
3. Angepasste Parameter ($\Delta G$, Vorfaktoren) werden mit Wert, Quelle und
   Empfindlichkeit der Zielgröße ausgegeben.
4. Für keine Größe wird „validiert" geschrieben, solange der zugehörige Eintrag
   in dieser Tabelle nicht auf `geprüft` steht.
