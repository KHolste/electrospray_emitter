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
| 10 | OpenMP + `-march=native` nicht als Compilerfehler abschließbar | **bestätigt; inzwischen aufgeklärt** | Die Kritik traf zu: der Punkt war unbelegt geschlossen. Inzwischen mit Minimalfall, Faktorexperiment (je 30 Läufe) und Disassemblat aufgeklärt — OpenMP ist unbeteiligt, Auslöser ist AVX-Stack-Ausrichtung. Siehe 1.2. |
| 11 | Empirische Skalengesetze ≠ selbstkonsistente Betriebspunktrechnung | **bestätigt als Konstruktionsfehler** | `cone_jet()` ist eine reine Formelauswertung ohne Kopplung an Geometrie, Feld oder Meniskus. Sie steht im selben Ausgabeblock wie die feldgekoppelte Ionenrechnung, was Gleichrangigkeit suggeriert. |

## 1.1b Stand nach Phase P0

Alle elf Befunde sind testgetrieben bearbeitet. Die Tests liegen in
`tests/test_regressions.cpp` und laufen dauerhaft mit.

| # | Korrektur | Regressionstest |
|---|---|---|
| 1 | `solve_at_voltage` liefert `SolveStatus::Converged` nur, wenn die gelieferte Spannung die angeforderte innerhalb `voltage_tol` (rel. 1e-3) trifft. Nicht erreichbare Spannungen: `VoltageNotBracketed`. Vorher wurde bei 500 V eine Form zu 871,8 V als konvergiert gemeldet; jetzt wird 500 V tatsächlich getroffen, weil unterhalb des Scout-Bereichs weitergesucht wird. | B1 |
| 2 | `find_onset` → `find_static_fold`. Verlangt mindestens drei konvergierte Punkte, ein **inneres** Maximum und strengen Anstieg davor und Abfall danach. Sonst `FoldStatus::TooFewPoints` / `Monotone` / `MaximumAtBoundary`. Der Test enthält den alten Algorithmus als Gegenprobe und zeigt, dass dieser genau die jetzt abgelehnten Fälle akzeptiert hätte. | B2 |
| 3 | Dateinamen tragen Anwendung, Zustand und Spannung (`chk_meniscus_fold_U1179p1V_surface.csv`), dazu einen `#`-Provenienzkopf. `MeniscusSolver::realize()` versetzt den Löser vor jedem Schreiben nachweislich in den ausgewiesenen Zustand. | B3 |
| 4, 5 | Raumladung schlägt geschlossen fehl (`NotImplementedInThisPhase`, Phase P4). Keine improvisierte Regularisierung; die Ringkerne bleiben als Testfall erhalten. | B4 |
| 6 | Tropfen als eigene `SpeciesKind::Droplet`; jede Anforderung wird abgelehnt (Phase P6), auch im Mischbetrieb. | B6 |
| 7 | `es_operate` behauptet keinen Betriebspunkt mehr. `print_operating_point` → `print_diagnostic_estimate` mit ausdrücklicher Kennzeichnung als nicht gekoppelte Abschätzung, Begründung über Higuera (2008) und Angabe der dG-Empfindlichkeit. | B7 |
| 8 | Negative Polarität schlägt geschlossen fehl, mit eigener Begründung und vor der Meniskusrechnung. | B8 |
| 9 | „99-%-Fläche" ersetzt durch das stetige Funktional `A_eff = (∫j dA)² / ∫j² dA`. Exakttest: für uniformes j (Kugel) ist `A_eff` die Gesamtfläche. Netzkonvergenz prüft jetzt Faltenspannung, Apexfeld, Apexradius, Ionenstrom und `A_eff`. | B9, test_meniscus |
| 10 | `-march=native` bleibt abgeschaltet; Diagnose und Reproduktionsfall dokumentiert (1.2). Keine weitere Arbeit am Bugreport. | — |
| 11 | Cone-Jet in eigenem Ausgabeblock, `ConeJetState::empirical`, Warnung bei ε_r < 40, ausdrücklicher Hinweis, dass die Werte nicht in den Strahltransport eingehen. | B11 |

Ein überladenes `converged`-Flag gibt es nicht mehr: `MeniscusSolution::status`
ist ein `SolveStatus`, und `ok()` ist der einzige zulässige Erfolgstest.

### Neu gefunden während P0

**`BemSolver::set_mesh` invalidierte die Basislösung nicht.** Ein anschließendes
`solve()` fand `basis_` nicht leer und superponierte die Basisvektoren des
**alten** Netzes auf die neue Geometrie. Aufgefallen ist es erst durch den
B3-Test: `realize()` lieferte 1,49·10⁹ V/m statt 5,04·10⁷ V/m, also Faktor 29.

Im Prototyp blieb der Fehler unbemerkt, weil die abschließende Feldauswertung in
`solve_at_height` ein Netz benutzte, das sich vom vorhergehenden nur um einen
Relaxationsschritt unterschied — der Fehler war klein und wanderte still in
`apex_field` und `peak_field` jeder Lösung. `set_mesh` löscht jetzt Basis,
Lösung und angelegte Potentiale.

### Zusatzbefund, nicht in der Liste

**Ausgabepräfixe kollidieren.** `es_meniscus` und `es_beam` schreiben bei
gleichem `output.prefix` in dieselben Dateinamen und überschreiben einander
kommentarlos. Das hat bei der Prüfung von Befund 3 zunächst zu einer falschen
Zwischenbewertung geführt (die untersuchte Datei stammte vom Strahllauf, nicht
vom Meniskuslauf). Für eine Parameterstudie ist das eine Fehlerquelle
derselben Klasse wie ein stillschweigend ignorierter Konfigurationsschlüssel.

### 1.2 Zum Absturz mit `-march=native`: aufgeklärt

**Zwei meiner früheren Aussagen waren beide unbelegt.** Zuerst hatte ich das im
Haupt-README als „Toolchain-Defekt, kein Codefehler" abgeschlossen. Danach hatte
ich das revidiert und auf latentes Undefined Behaviour im Projektcode getippt.
Beide Male fehlte die Grundlage. Erst die folgende Untersuchung klärt es.

**Fehler in meiner vorigen Bisektion.** Der Absturz ist **nichtdeterministisch**
— je nach Konfiguration 0 bis 19 von 30 Läufen. Meine Leave-one-out-Bisektion
über die Übersetzungseinheiten hatte jede Konfiguration **genau einmal**
ausgeführt. Bei einer Ausfallrate um 50 % ist das reines Rauschen; die daraus
gezogene Folgerung („verschwindet und kehrt zurück, wenn Flags an unbeteiligten
Übersetzungseinheiten umgeschaltet werden") war nicht haltbar.

**Minimaler Reproduktionsfall.** 15 Zeilen, ein einziger `solve_at_height`-Aufruf,
0,25 s Laufzeit, 16 von 30 Läufen mit SIGSEGV. Liegt als
[`docs/repro/avx_stack_alignment.cpp`](repro/avx_stack_alignment.cpp) im
Repositorium; für einen Bugreport genügt diese Datei plus die Bibliothek.

**Faktorexperiment, je 30 Läufe:**

| Konfiguration | SIGSEGV |
|---|---:|
| `-march=native` + OpenMP-Pragmas | 6 / 30 |
| `-march=native`, `-fopenmp` gelinkt, keine Pragmas | 12 / 30 |
| `-march=native`, **kein OpenMP** | **14 / 30** |
| ohne `-march=native`, OpenMP-Pragmas | 0 / 30 |
| weder noch | 0 / 30 |
| `-mavx2` | 16 / 30 |
| `-mavx` | 19 / 30 |
| `-msse4.2` | 0 / 30 |
| `-mavx2 -mstackrealign` | 18 / 30 |
| `-mavx2 -mprefer-vector-width=128` | **0 / 30** |
| `-mavx2 -fno-tree-vectorize` | 11 / 30 |
| `-mavx2 -O0` | 15 / 30 |

**OpenMP ist unbeteiligt** — ohne jedes OpenMP stürzt es genauso ab. Auslöser
ist 256-Bit-Vektorcodegen (AVX/AVX2); SSE4.2 ist unauffällig. Es tritt auch bei
`-O0` und mit abgeschalteter Auto-Vektorisierung auf, ist also keine Frage der
Schleifenoptimierung.

**Ursache, im Disassemblat belegt.** Unter gdb gefangenes Signal:

```
Thread 1 received signal SIGSEGV, Segmentation fault.
0x...49b in main () at min3.cpp:9
      MeniscusSolver s(make_capillary_open(cp), mp);
=> 0x...49b <main()+203>:  vmovdqa %ymm0,0x20(%rsp)
rsp   0x5ffaf0
```

`vmovdqa` ist ein alignment-pflichtiger 256-Bit-Store. Zieladresse
`0x5ffaf0 + 0x20 = 0x5ffb10`, und `0x5ffb10 mod 32 = 16` — die Adresse ist nur
16-Byte-ausgerichtet.

Der Prolog von `main()` bestätigt es:

```
push %rsi ; push %rbx ; sub $0x378,%rsp ; call __main
```

Keine Stack-Nachrichtung (`and $-32,%rsp`: **0 Treffer** im gesamten `main`),
aber **5** alignment-pflichtige 256-Bit-Zugriffe im selben Frame. Die Win64-ABI
garantiert nur 16 Byte Stack-Ausrichtung. GCC nimmt hier 32 Byte an, ohne sie
herzustellen.

Das erklärt auch den Nichtdeterminismus: ob der Stack beim Prozessstart zufällig
auf einer 32-Byte-Grenze liegt, entscheidet die Adressraum-Randomisierung.

**Gegenprüfungen:**

* ASan und UBSan unter Linux (GCC 13.3, Ubuntu 24.04): **keine Befunde**, weder
  im Standardbau noch mit `-march=native`, weder in den Tests noch in den
  Anwendungen.
* Der Absturz **reproduziert unter Linux nicht** (GCC 13.3, gleiche Flags, je
  drei Läufe von `test_meniscus`, `test_beam`, `test_bem` — alle sauber).
* Unter gdb läuft das große Testprogramm durch; erst der schnelle
  Minimalfall ließ sich unter gdb fangen. Klassischer Heisenbug: der Debugger
  verändert das Stack-Layout.
* `-mstackrealign` behebt es **nicht** — das ist Teil des Befunds, denn genau
  dieses Flag sollte die Nachrichtung erzwingen.

**Bewertung.** Codegenerierungsfehler in diesem MinGW-w64-GCC-16.1-Bau, kein
Undefined Behaviour im Projektcode. Der Nachweis stützt sich auf das
Disassemblat und die Adressarithmetik, nicht auf ein Ausschlussverfahren.

**Was noch aussteht:** ein Reproduktionsfall **ohne Projektabhängigkeit**. Ein
erster eigenständiger Kandidat (Struct mit acht `double`, `-mavx2`) erzeugt zwar
ebenfalls zwei alignment-pflichtige 256-Bit-Zugriffe ohne Nachrichtung, stürzt
aber in 0 von 30 Läufen ab — bei ihm fallen die Frame-Offsets zufällig auf
32-Byte-Grenzen. Für einen Bugreport gegen GCC ist das noch zu erarbeiten; die
Diagnose hängt nicht davon ab.

**Konsequenz für den Code.** `-march=native` bleibt abgeschaltet. Das ist kein
fachfremder Workaround, sondern die Vermeidung einer nachgewiesen defekten
Codegenerierung; die Begründung steht jetzt im Haupt-README statt einer
Vermutung. `-mprefer-vector-width=128` wäre die Alternative, falls
Host-CPU-Codegen später gebraucht wird.

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
