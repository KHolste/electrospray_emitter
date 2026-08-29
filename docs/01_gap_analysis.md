# 1. Gap-Analyse

## 1.1 Verifikation der gemeldeten Fehler

Alle elf Punkte wurden am Code nachgeprüft, nicht übernommen. Das
Prüfprogramm liegt außerhalb des Repositoriums und benutzt nur die öffentliche
API; `src/` wurde nicht verändert.

Geometrie aller Läufe, soweit nicht anders vermerkt: Kapillare ∅ 20 µm
Bohrung, 40 µm außen, Extraktor 500 µm entfernt, Apertur ∅ 400 µm, EMI-BF4.

| # | Befund | Status | Belegte Messung |
|---|---|---|---|
| 1 | `solve_at_voltage` meldet Konvergenz bei falscher Spannung | **bestätigt** | Anforderung 500 V → gelieferte Form gehört zu 871,8 V (+74,4 %), `converged=1`. Bei 800 V → 871,8 V (+9,0 %). Ab ca. 1000 V korrekt. |
| 2 | `find_onset` meldet Umkehrpunkt bei einem einzigen Punkt | **bestätigt, weiter reichend** | Ast mit 1 konvergiertem Punkt → `found=1`. Zusätzlich: streng **monoton steigender** Ast ohne jedes Maximum → `found=1`, gemeldet wird der letzte Punkt. |
| 3 | Ausgabedateien gehören zum letzten nicht konvergierten Zustand | **bestätigt** | `chk_surface.csv`: Spitzenfeld 1,503·10⁸ V/m = exakt der letzte, nicht konvergierte Astpunkt. Ausgewiesener Onset: 4,797·10⁷ V/m. Faktor 3,13. Apexhöhe in `chk_mesh.csv` 7,667 µm statt 4,143 µm. |
| 4 | Eigenfeld des Ring-Makroteilchens singulär | **bestätigt** | Exakt am Ring: `E_r = nan`, `E_z = nan`, `V = inf`. Annäherung: \|E\| = 5,7·10⁶ / 5,7·10⁸ / 5,7·10¹⁰ / 5,7·10¹² V/m bei 10⁻⁸ / 10⁻¹⁰ / 10⁻¹² / 10⁻¹⁴ m — Divergenz wie 1/d, Potential wie ln(1/d). |
| 5 | Keine physikalisch definierte Selbstwechselwirkung | **bestätigt** | Folgt aus 4. Das Modell hat keine Teilchenausdehnung, also auch keinen definierten Nahfeld-Abschneideradius. |
| 6 | Tropfen werden mit der Ionen-Emissionsrate gewichtet | **bestätigt** | Startgewichte in `trace_beam` = Σ `ion_current_density`·A = 4,744·10⁻¹² A. `cone_jet(Q = 0,1 nL/s)` = 3,945·10⁻⁷ A. Verhältnis 8,3·10⁴. Der Cone-Jet-Strom geht nirgends in den Strahltransport ein. |
| 7 | 1100 V liefert nur ~4,2·10⁻¹⁶ A | **bestätigt** | Gemessen 4,215·10⁻¹⁶ A bei Apexfeld 0,032 V/nm. Das sind rund 2600 Ionen/s. Kein Betriebspunkt. |
| 8 | Negative Emission nicht umgesetzt | **bestätigt** | U = +1500 V und U = −1500 V liefern identisches Spitzenfeld (1,193·10⁸ V/m) und identischen Strom (4,744·10⁻¹² A). `qm_cluster()` benutzt ausschließlich `M_cation`. Eine Anion-Clusterreihe hätte hier q/m = 5,81·10⁵ statt 5,07·10⁵ C/kg. |
| 9 | Netzkonvergenz nur für die Onset-Spannung geprüft | **bestätigt, mit Einschränkung** | Der Test prüft nur U_onset. Nachgeholt für 61/81/121/161 Knoten: U_onset 1179,82 → 1179,76 V; E_apex 5,1124 → 5,1056·10⁷ V/m (0,13 %); I_ion 1,9122 → 1,9132·10⁻¹⁵ A (0,05 %); A_emit 3,198 → 3,150·10⁻¹⁰ m² **nicht monoton** (3,198 / 3,150 / 3,103 / 3,150). Die Kritik am Test ist berechtigt; die Größen sind an diesem Punkt zufällig konvergiert, das war nie gezeigt. A_emit ist als diskrete Elementsumme ohnehin quantisiert und kein zur Netzverfeinerung glatt konvergierender Funktionalwert. |
| 10 | OpenMP + `-march=native` nicht als Compilerfehler abschließbar | **bestätigt — meine frühere Bewertung war falsch** | siehe 1.2 |
| 11 | Empirische Skalengesetze ≠ selbstkonsistente Betriebspunktrechnung | **bestätigt als Konstruktionsfehler** | `cone_jet()` ist eine reine Formelauswertung ohne Kopplung an Geometrie, Feld oder Meniskus. Sie steht im selben Ausgabeblock wie die feldgekoppelte Ionenrechnung, was Gleichrangigkeit suggeriert. |

### Zusatzbefund, nicht in der Liste

**Ausgabepräfixe kollidieren.** `es_meniscus` und `es_beam` schreiben bei
gleichem `output.prefix` in dieselben Dateinamen und überschreiben einander
kommentarlos. Das hat bei der Prüfung von Befund 3 zunächst zu einer falschen
Zwischenbewertung geführt (die untersuchte Datei stammte vom Strahllauf, nicht
vom Meniskuslauf). Für eine Parameterstudie ist das eine Fehlerquelle
derselben Klasse wie ein stillschweigend ignorierter Konfigurationsschlüssel.

### 1.2 Zum OpenMP-Absturz: revidierte Bewertung

Ich hatte im Haupt-README geschrieben, das sei „ein Toolchain-Defekt, kein
Codefehler". **Das war nicht belegt und ist nach den jetzigen Messungen eher
unwahrscheinlich.**

Bisektion über die Übersetzungseinheiten (GCC 16.1 MinGW-w64, OpenMP aktiv,
Testfall `test_meniscus`):

| Konfiguration | Ergebnis |
|---|---|
| alle Quellen + Test mit `-mavx2` | SIGSEGV |
| alle Quellen mit `-mavx2`, Test ohne | läuft durch |
| nur Test mit `-mavx2`, Quellen ohne | **SIGSEGV** |
| je genau eine Quelldatei mit `-mavx2` (7 Varianten) | alle laufen durch |
| alle außer `linalg` / `bem` / `meniscus`, Test mit `-mavx2` | läuft durch |
| alle außer `elliptic` / `geometry` / `fluid` / `emission`, Test mit `-mavx2` | SIGSEGV |

Das Bild ist nicht das eines lokal falsch übersetzten Codestücks. Ein Absturz,
der beim Umschalten von Optimierungsflags an *unbeteiligten* Übersetzungs-
einheiten verschwindet und wiederkommt, ist die typische Signatur von latentem
Undefined Behaviour im Programm — ein Zugriff außerhalb der Grenzen, der bei
einem Speicherlayout harmlos landet und bei einem anderen nicht.

Der Absturz tritt nach der letzten Ausgabezeile auf, also beim Abbau statischer
Objekte. Das passt zu einer Heap-Korruption, die erst beim Freigeben auffällt.

Sanitizer-Lage in dieser Umgebung:

* `-fsanitize=address,undefined`: `libasan`/`libubsan` sind in diesem
  MinGW-Build nicht vorhanden.
* `-D_GLIBCXX_DEBUG`: übersetzt, startet aber nicht (Rückgabewert 127 ohne
  Ausgabe, keine fehlenden DLLs) — bekannte ODR-Inkompatibilität mit der
  vorgebauten libstdc++.
* `-D_GLIBCXX_ASSERTIONS` bei `-O0` ohne OpenMP: läuft sauber durch, schlägt
  nicht an. Das schließt einen Bereichsfehler in einem `std::vector::operator[]`
  auf diesem Pfad **nicht** aus, weil dieser Schalter nur einen Teil der
  Prüfungen aktiviert.
* **WSL2 mit Ubuntu ist auf dem Rechner vorhanden.** Das ist der belastbare Weg
  zu ASan/UBSan.

**Status: offen.** Kein Compilerfehler, kein Codefehler — unbelegt. Der
Workaround (`-march=native` nicht als Voreinstellung) bleibt, ist aber als
Notbehelf zu kennzeichnen und nicht als Lösung. Auflösung in Phase P0 des
Stufenplans.

---

## 1.3 Systematische Lücken gegenüber der Zielsetzung

Die Einzelfehler oben sind reparierbar. Die folgenden Punkte sind
Konstruktionslücken und der eigentliche Grund für die Neuausrichtung.

### G1 — Der emittierende Betrieb ist gar nicht modelliert

Der Prototyp hat zwei Zustände: einen statischen Meniskus **ohne** Emission und
eine Ionenrate, die *nachträglich* auf das Feld dieses nicht-emittierenden
Meniskus angewandt wird. Der emittierende Betrieb ist etwas anderes: dort
fließt Flüssigkeit, es gibt einen viskosen Druckabfall, Ladung wird ohmsch
antransportiert und an der Oberfläche abgeführt, und die Oberfläche ist deshalb
**keine** Äquipotentialfläche mehr.

Higuera (2008, Phys. Rev. E 77, 026308) zeigt für genau diesen Fall, dass die
Strömung im Meniskus viskositätsdominiert ist und der Verdampfungsstrom von der
**endlichen Leitfähigkeit** kontrolliert wird. Das Perfect-Conductor-Modell ist
im emittierenden Betrieb also nicht nur ungenau, sondern trifft den
kontrollierenden Mechanismus nicht.

Konsequenz: die Ionenströme des Prototyps sind keine Vorhersage. Sie sind die
Auswertung einer stark nichtlinearen Funktion an einem Feld, das aus einem
Modell stammt, dessen Gültigkeitsbereich diesen Betriebszustand ausschließt.

### G2 — Der Umkehrpunkt wird ohne Nachweis als Onset ausgegeben

Drei verschiedene Dinge werden im Prototyp gleichgesetzt:

1. der Verlust der statischen Stabilität des Meniskus (eine Bifurkation),
2. der Emissionsbeginn im Sinne eines messbaren Stroms (eine Gerätegrenze),
3. der Übergang in den Cone-Jet-/Tröpfchenbetrieb (eine andere Bifurkation,
   die die Strömung einschließt).

Im PIR gibt es zu (2) gar keinen Schwellwert: die Ionenverdampfungsrate ist eine
glatte Exponentialfunktion des Feldes. Was gemessen wird, ist der Punkt, an dem
der Strom die Nachweisgrenze überschreitet. Der Umkehrpunkt eines statischen
Astes ist ein Kandidat für (1) und nicht mehr — und auch das nur unter
Zusatzannahmen über den festgehaltenen Kontrollparameter und den betrachteten
Störungsmodus. Details und ein literaturbelegtes Kriterium in
[02_model_specification.md](02_model_specification.md), Abschnitt 3.4.

### G3 — Kein Volumenverfahren für Raumladung

BEM löst die *homogene* Laplace-Gleichung. Raumladung ist ein Quellterm im
Volumen und damit außerhalb dessen, was eine reine Randintegralformulierung
darstellen kann. Der Prototyp behilft sich mit direkt aufsummierten
Ring-Makroteilchen ohne Regularisierung; das ist aus den unter Befund 4 und 5
genannten Gründen nicht haltbar und zusätzlich O(N²).

### G4 — Keine Adaptivität, Geometrie passt nicht zur Skizze

Die Elementgrößen werden im Konfigurationsfile von Hand gesetzt (`h_tip`,
`h_far`). Es gibt keinen Fehlerschätzer und keine automatische Verfeinerung.
Die vorhandenen Bausteine (`make_capillary`, `make_needle`) bilden die
Skizzengeometrie — konisch verjüngter Emitter mit ∅₁, ∅₂, ∅₃ und Höhe — nicht
ab.

### G5 — Stoffdaten ohne Quellenangabe und ohne Unsicherheit

`src/fluid.cpp` enthält Zahlenwerte ohne Referenz. Der Kopfkommentar warnt
zwar vor der Streuung, aber ein Kommentar ist keine Quellenangabe. Für eine
wissenschaftlich zitierbare Rechnung muss jeder Wert eine Quelle und eine
Unsicherheit tragen.

### G6 — Polarität, Speziesverteilung, Solvatationsverteilung fehlen

Modelliert ist genau eine Spezies mit einer mittleren Solvatationszahl. Reale
IL-Emission liefert eine Verteilung über Monomer, Dimer, Trimer, plus im
Mischbetrieb Tropfen; jede Spezies hat ihr eigenes q/m und einen eigenen
Anteil, und die Anionenseite unterscheidet sich in Masse und
Solvatationsenergie von der Kationenseite.

### G7 — Validierung ausschließlich gegen analytische Elektrostatik

Was tatsächlich gegen eine unabhängige Referenz geprüft ist: die elliptischen
Integrale, die Kugelkapazität, das Spitzenfeld des Rotationsellipsoids, die
Kugelkappe im feldfreien Grenzfall, die Energieerhaltung entlang einer
Trajektorie. Das ist die Elektrostatik und die Bahnintegration — also genau die
Teile, die *nicht* die Physik des Emitters ausmachen.

Der im Haupt-README hervorgehobene Halbwinkel von 49,16° ist **kein**
unabhängiger Nachweis in dem Sinne, wie er dort dargestellt wird: er zeigt, dass
der Young-Laplace-Löser mit Maxwell-Term das asymptotische Verhalten eines
Kegels reproduziert, wenn man dem Ast weit genug folgt. Das ist ein sinnvoller
Konsistenztest der Diskretisierung, aber er sagt nichts über die
Emissionsphysik, den Betriebspunkt oder die Stabilität aus, und er liegt zudem
auf dem instabilen Ast.
