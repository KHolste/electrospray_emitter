# electrospray — Simulation von Kolloid-/Electrospray-Emittern

Axialsymmetrischer Simulationscode für Emitter ionischer Flüssigkeiten:
Elektrostatik, elektrifizierter Meniskus (Taylor-Cone), Emissionsmodell und
Strahltransport. C++20, keine externen Abhängigkeiten außer CMake und einem
C++-Compiler.

```
Stufe 1  Elektrostatik        axialsymmetrische BEM (Randelementmethode)
Stufe 2  Meniskus             Young-Laplace + Maxwell-Stress, freie Oberfläche
Stufe 3  Emission             Iribarne-Thomson (PIR) + Cone-Jet-Korrelation
Stufe 4  Strahltransport      Ring-Makroteilchen mit selbstkonsistenter Raumladung
```

---

## 1. Warum BEM und nicht FEM

Die Unbekannte ist die Flächenladungsdichte σ auf dem Rand; das Potential ist

    V(r,z) = ∫ σ(s') G(r,z; r',z') ds'
    G      = r' K(m) / (π ε₀ S),   S² = (r+r')² + (z−z')²,   m = 4rr'/S²

mit K dem vollständigen elliptischen Integral erster Art. Das ist die einmal
über den Azimut integrierte Freiraum-Greensfunktion. Vier Konsequenzen, die
genau auf dieses Problem passen:

* **Nur der Rand wird vernetzt.** Der Meniskus-Solver bewegt die freie
  Oberfläche in jeder Iteration — Neuvernetzung kostet nichts.
* **Kein künstlicher Fernfeldrand.** Das Außenraumproblem geht bis unendlich.
  Bei Nadel-Extraktor-Geometrien ist das der entscheidende Vorteil.
* **E_n = σ/ε₀ exakt.** Das Normalfeld an der Oberfläche ist die *primäre*
  Unbekannte, nicht eine numerische Ableitung des Potentials. Da sowohl der
  Maxwell-Stress als auch die Feldverdampfungsrate von genau dieser Größe
  getrieben werden, ist das der Unterschied zwischen brauchbar und nutzlos.
* **N ~ 10³** in Axialsymmetrie — die dichte Matrix ist irrelevant.

Die logarithmische Kernsingularität wird durch geometrisch verfeinerte Panels
aufgelöst; das komplementäre Argument `mc = 1 − m` wird direkt als `d²/S²`
gebildet, ohne Auslöschung bei m → 1.

## 2. Der Meniskus-Solver — und warum h die Fortsetzungsvariable ist

Auf der freien Oberfläche gilt

    γ ( dφ/ds + sin(φ)/r ) = Δp − ρ g z + (ε₀/2) E_n²

Physikalische Eingangsgrößen sind der Speisedruck Δp und die Spannung U; die
Apex-Höhe h ist Ergebnis. **Aber h(U) hat am Emissionsbeginn eine senkrechte
Tangente** — das *ist* der Onset: eine Sattel-Knoten-Bifurkation, jenseits derer
keine statische Lösung mehr existiert. Ein Marsch in U kann diesen Punkt
prinzipiell nie erreichen; die Iteration hört irgendwo davor auf zu
konvergieren, und das Ergebnis sieht aus wie ein numerisches Versagen statt wie
die Physik, die es ist.

Der Solver dreht deshalb die Rollen um: **h wird vorgegeben, U gesucht.** Diese
Parametrisierung ist durch den Umkehrpunkt hindurch regulär, der ganze Ast
inklusive Fold wird verfolgt, und die Onset-Spannung ist das Maximum von U(h).
Alles danach ist der instabile Ast — reale Lösungen, aber bei fester Spannung
nicht realisierbar.

Ein Nebenprodukt, das als Validierung taugt: **der Taylor-Winkel 49.3° steht
nirgends im Code.** Er entsteht aus der gekoppelten Rechnung, wenn man dem Ast
folgt (siehe unten).

## 3. Bauen

```sh
cmake -S . -B build -G Ninja
cmake --build build
cd build && ctest --output-on-failure
```

Getestet mit MinGW-w64 GCC 16.1 unter Windows 11, CMake 4.3, Ninja. OpenMP wird
automatisch genutzt, wenn vorhanden.

**Zu `-march=native`:** bewusst nicht Standard. Auf MinGW-w64 GCC 16.1 stürzt
der Code beim Programmende ab, sobald `-mavx2` und libgomp zusammenkommen —
keins von beiden allein, und `-mstackrealign` hilft nicht. Das ist ein
Toolchain-Defekt, kein Codefehler (die Zahlen stimmten in allen geprüften
Fällen mit dem skalaren Build überein). Wer es trotzdem will:
`cmake -S . -B build -DES_NATIVE_ARCH=ON`.

## 4. Anwendungen

Alle Programme lesen eine Config-Datei und akzeptieren `key=value`-Overrides
auf der Kommandozeile (Kommandozeile schlägt Datei). Werte tragen Einheiten:
`10um`, `1.5kV`, `500Pa`, `0.5nL/s`, `25C`, `15deg`.

| Programm      | Zweck |
|---------------|-------|
| `es_field`    | Laplace-Lösung, Spitzenfeld, Feldüberhöhung, Spannungssweep, Feldgitter (CSV/VTK) |
| `es_meniscus` | Fortsetzung des Gleichgewichtsastes → Onset-Spannung, Meniskusform |
| `es_operate`  | Vollständiger Betriebspunkt: Meniskus bei U, Ionenstrom, Cone-Jet-Vergleich, Schub/Isp |
| `es_beam`     | Strahltransport, Divergenz, Extraktor-Interzeption, Raumladung |

```sh
./build/es_meniscus examples/capillary_emibf4.cfg
./build/es_meniscus examples/capillary_emibf4.cfg emitter.r_bore=5um extractor.gap=300um
./build/es_field    examples/capillary_emibf4.cfg output.grid=true
./build/es_field --help          # vollständige Schlüsselreferenz
```

Nicht gelesene Config-Schlüssel werden am Ende gemeldet — ein stillschweigend
ignorierter Tippfehler ist in einer Parameterstudie der teuerste Fehler.

Plots aus den CSV-Dateien:

```sh
python python/plot.py out_capillary        # Meniskusast, Form, Oberflächenfeld
python python/plot_beam.py out_capillary   # Trajektorien und Strahlprofil
```

## 5. Verifikation

`ctest` prüft gegen analytische Lösungen, nicht nur gegen sich selbst.

| Test | Referenz | Ergebnis |
|------|----------|----------|
| Elliptische Integrale K, E | AGM gegen Simpson-Quadratur und exakte Werte | 1e-14 |
| Kugelkapazität | 4πε₀R | 2.5e-6 (640 Elemente) |
| Rotationsellipsoid, Aspekt 20 | Prolat-sphäroidale Analytik, Spitzenfeld | 0.18 % (600 Elemente) |
| Randbedingungs-Residuum | Potentialintegral an allen Kollokationspunkten | 1.5e-15 |
| Fernfeld-Monopol | Q/(4πε₀R) bei 0.5 m | 1e-3 |
| Meniskus ohne Feld | Kugelkappe R = 2γ/Δp | 4e-6 |
| Onset, Netzkonvergenz | 81 vs. 141 Knoten | 2.5e-5 |
| Onset vs. Taylor/Smith | geschlossene Formel | Faktor 0.86 |
| **Kegelhalbwinkel entlang des Astes** | Taylor 49.29° | **49.16°** |
| Energieerhaltung im Strahl | E_kin/q = Potentialdifferenz | 7e-6 |
| Ring-Kern | E = −∇V, Punktladungslimes | 1e-9 |

Der Halbwinkel läuft entlang des Astes monoton
81° → 75° → 68° → 62° → 56° → 51° → **49.16°** gegen den Taylor-Winkel.

## 6. Was der Code kann — und was nicht

**Gültig:**
* Onset-Spannung und Meniskusform für kapillare Emitter mit gepinnter
  Kontaktlinie, statisch oder bei kleinen Flüssen (PIR-Bereich).
* Feldüberhöhung beliebiger axialsymmetrischer Elektrodengeometrien.
* Strahloptik und Extraktor-Interzeption, mit und ohne Raumladung.
* Betriebspunkt-Sweeps und Geometrieoptimierung.

**Bekannte Grenzen — bitte lesen:**

1. **Der statische Meniskus löst die Emissionsstruktur nicht auf.** Bei einer
   10-µm-Kapillare beträgt der Apex-Krümmungsradius am Onset noch ~9 µm und das
   Feld nur ~0.05 V/nm. Feldverdampfung braucht ~1 V/nm, also einen
   Krümmungsradius im 10-nm-Bereich. Diese Struktur ist der *Jet* bzw. die
   Kegelspitze, die das statische Modell nicht enthält. Mit 200-nm-Bohrung
   erreicht die Rechnung 0.38 V/nm am Onset und 1.3 V/nm auf dem instabilen Ast
   — das reproduziert korrekt, dass PIR sub-µm-Emitter erfordert, ersetzt aber
   kein transientes Jet-Modell.
2. **Keine Strömung.** Der viskose Druckabfall entlang des Kegels fehlt. Bei
   Cone-Jet-Flussraten überschätzt das Modell die Apex-Höhe.
3. **Kontaktlinie gepinnt und r-monoton.** Korrekt für eine scharfkantige
   Kapillare. Für extern benetzte oder poröse Emitter ist die Behandlung der
   Spitze als „gepinnt bei r_tip" nur eine Näherung erster Ordnung.
   Überhängende Menisken (über die Halbkugel hinaus) werden bewusst nicht
   verfolgt.
4. **Die Cone-Jet-Vorfaktoren sind empirisch.** `I = f(ε_r)·√(γKQ/ε_r)` wurde an
   Flüssigkeiten mit ε_r ≳ 40 etabliert. Ionische Flüssigkeiten liegen bei
   ε_r ≈ 12, weit außerhalb. Der Code extrapoliert **nicht** heimlich: er hält
   `f_current` als expliziten Eingabewert, setzt ihn auf den
   Hochpermittivitätswert und **warnt**, wenn er außerhalb des Gültigkeitsbereichs
   benutzt wird. Vor Vertrauen in Absolutströme an eigenen I(Q)-Daten fitten.
5. **ΔG_solvation ist der dominierende Fitparameter.** Literaturwerte für
   EMI-BF4 streuen etwa 1.0–1.4 eV, und der Strom hängt exponentiell davon ab.
   Als Fitgröße behandeln, nicht als Konstante.
6. **Stoffdaten streuen.** Die Tabelle in `src/fluid.cpp` enthält
   literaturtypische Raumtemperaturwerte als Startpunkt, keine Messwerte der
   eigenen Charge. Leitfähigkeit ist besonders wasserempfindlich. Alles
   überschreibbar über `[fluid]`.
7. **Die Temperaturmodelle (VFT) sind generisch**, nicht flüssigkeitsspezifisch
   gefittet.

## 7. Struktur

```
include/es/     types constants elliptic linalg geometry bem fluid
                emission meniscus beam config io
src/            Implementierungen
apps/           es_field es_meniscus es_operate es_beam
tests/          test_elliptic test_bem test_meniscus test_emission test_beam
examples/       Config-Beispiele
python/         Plot-Skripte (nur Auswertung, nicht Rechnung)
```

Die Bibliothek `es_core` ist ohne die Anwendungen nutzbar; jedes Modul hat einen
Header mit den Annahmen und Gültigkeitsgrenzen im Kopfkommentar.

## 8. Naheliegende Erweiterungen

* Viskoser Druckabfall entlang des Kegels → Cone-Jet-Betriebspunkte
* Blockweises Matrix-Update: bei bewegtem Meniskus ändert sich nur ein Teil der
  Einflussmatrix; der Elektrodenblock ließe sich cachen (~3× im Meniskus-Loop)
* Emitter-Arrays: das Green'sche Funktionsgerüst ist axialsymmetrisch, ein Array
  bräuchte 3D-BEM oder ein Einheitszellenmodell mit Spiegelladungen
* Transientes Jet-Modell für die nm-Skala am Apex
* Time-of-Flight-Spektren aus der Clusterverteilung
