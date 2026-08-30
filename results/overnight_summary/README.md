# Nachtlauf 2026-08-29/30 und die Nacharbeit vom 2026-08-30

Vor diesem Bereinigungscommit stand `origin/main` auf `a932ff1` und HEAD
**46 Commits voraus, 0 zurück**, Arbeitsbaum sauber. Dieser Commit ist der 47.
und der Stand, der nach `origin/main` geht.

### Abnahme auf dem endgültigen Code-HEAD

Beide Läufe wurden **nach** allen fachlichen Korrekturen auf `ea2d5a1`
wiederholt; die Logs liegen unter [`logs/`](logs/) mit Befehl, Toolchain,
Commit und Exit-Code.

| Lauf | Toolchain | Ergebnis | Log |
|---|---|---|---|
| Clean Configure + Build + `ctest` | Windows, MinGW-W64 GCC 16.1.0, CMake 4.3.3, Ninja 1.13.2 | configure/build/ctest je Exit-Code 0, **26/26 grün**, 119,46 s | [`logs/build_ctest_final.log`](logs/build_ctest_final.log) |
| ASan + UBSan, `detect_leaks=1` | WSL2 Ubuntu 24.04.4, GCC 13.3.0, CMake 3.28.3 | configure/build/ctest je Exit-Code 0, **26/26 grün**, 606,44 s, **null** ASan-, UBSan- oder Leak-Meldungen | [`logs/sanitizer_final.log`](logs/sanitizer_final.log) |

**Was der frühere Sanitizersatz zu weit ging.** Er las sich, als sei der Lauf
über den endgültigen HEAD selbst gefahren. Tatsächlich geprüft wurde `cc27a46`.
Der Unterschied ist mechanisch nachweisbar unerheblich, aber er gehört benannt:
`git diff cc27a46 HEAD -- src include apps tests tools CMakeLists.txt` ist
**leer**, der kompilierte Inhalt beider Stände ist also identisch; zwischen
`cc27a46` und HEAD liegen nur ein Paneltitel in `python/plot_transport.py` und
Artefakte. Der oben aufgeführte Lauf schließt die Lücke trotzdem, indem er
tatsächlich auf `ea2d5a1` fährt. Beide Logs tragen dazu einen Nachtrag, der
offenlegt, warum die Zeile `Arbeitsbaum` dort `DIRTY` zeigt (unversioniertes
Logverzeichnis und die parallel laufende Textbereinigung) und dass der
Vergleich der kompilierten Dateien gegen HEAD dabei leer war.

### Abbildungen — die Zahlen einzeln

Die frühere Angabe „27 Abbildungen insgesamt" gehört zu keiner nachprüfbaren
Menge und ist ersetzt. Ausgezählt auf dem Dateisystem:

| Menge | Ordner | PNGs |
|---|---|---|
| der **P0–P9-Lauf** (die zehn Ordner aus Abschnitt 6) | 10 | **25** |
| Ordner mit `figures_provenance.txt` | **15** | 49 |
| **alle** Ergebnisordner mit Abbildungen | **18** | **56** |

Die drei Ordner ohne Provenienzstempel — `branch_ambiguity`,
`p1_boundary_mesh`, `p1_device_geometry` — stammen aus der Zeit vor
`provenance.py` und tragen die verbleibenden 7 Abbildungen. Sie sind der Grund,
warum sich 56 **nicht** auf 15 Ordner verteilt; die beiden Zahlen gehören zu
verschiedenen Mengen und dürfen nicht in einem Satz stehen.

Geprüft wurde getrennt, und beides ist nötig:

* **mechanisch**, alle 56: PNG-Signatur, `IEND`-Abschluss, Dateigröße > 0,
  Abmessungen > 0 — 56 von 56 in Ordnung, kleinste Abbildung 1216×832,
  kleinste Datei 90 452 Bytes;
* **inhaltlich**, die neun Abbildungen aus P3, P4, P5, P7 und P8, die nach dem
  P2-Nachzug neu erzeugt und noch nicht im Endzustand vorgelegt worden waren:
  Kurven gegen die zugrunde liegenden CSV nachgerechnet, Achsen, Einheiten,
  Legenden und Bildunterschriften gelesen, Statusaussagen gegen den aktuellen
  Code gehalten. Einzelheiten in Abschnitt 7.

Ein mechanisch gültiges PNG ist keine geprüfte Abbildung; die erste Prüfung
sagt, dass eine Datei lesbar ist, die zweite, dass sie stimmt.

Jeder Provenienzstempel liegt in dieser Historie (`git merge-base
--is-ancestor`), und jedes der 15 gestempelten Ergebnisverzeichnisse trägt
`working_tree_dirty=no` und `releasable=yes`.

Sicherung des Standes vor der Nacharbeit: Branch
`backup/nachtlauf-p0-p9-20260830` und Tag `backup-nachtlauf-20260830-ede1508`,
beide auf `ede1508`. Es wurde nichts verworfen.

---

## 1. Was die Nacharbeit geändert hat

Vier fachliche Korrekturen, eine Reparatur der lokalen Historie, der Nachzug
dieser Korrekturen in acht weitere Abbildungen — in **zwei** Runden, weil die
erste nur die Bildunterschriften traf und ein Paneltitel in P3 die überholte
Aussage noch weitertrug —, und fünf Artefaktsätze, deren Provenienz durch die
Reparatur ungültig geworden war. Jede Korrektur ist ein eigener Code-Commit mit
einem eigenen Artefakt-Commit dahinter, und jede ist testgeprüft.

| Punkt | Befund | Code | Artefakte |
|---|---|---|---|
| **P2** | ν galt pauschal als `MissingMaterialData` — obwohl µ und ρ aus **derselben Publikation und derselben Probe** stammen | `2b95a94` | `572f1e8` |
| **P2** | `maxwell_force_for_fixed_shape (linear in gamma)` ist physikalisch falsch | `2b95a94` | `572f1e8` |
| **P3** | „nicht DC, also unbrauchbar" verwarf genau die Messung, die die Formel verlangt | `bfa2c61` | `fa0e090` |
| **P6** | „keine Divergenz beim Annähern" galt nur bei **festem Netz** | `4c3b83b` | `593cecd` |
| **P9** | eine einachsige Statuskarte ließ „vergleichbar" wie „validiert" aussehen | `3bef623` | `16cdd09` |
| **P9** | eine fehlende Unsicherheit war ein harter Importfehler und warf Messungen weg | `3bef623` | `16cdd09` |
| **Historie** | `713ab86` ließ sich nicht eigenständig konfigurieren | — | — |
| **P0/P1/P4/P5/P7** | Abbildungen stempelten Commits außerhalb der eigenen Historie | — | `8040125`, `e46c941`, `91f254c`, `bed39d7`, `9d1b92c` |
| **P2/P3/P8** | drei Bildunterschriften trugen Aussagen weiter, die die Korrekturen überholt haben | `cc27a46` | `5d8caba`, `7b683d4`, `ee6b743` |
| **P3** | der Paneltitel von Abb. 2 nannte weiterhin einen unbelegten ε_r — die Bildunterschrift allein zu korrigieren hatte nicht gereicht | `42f4d73` | `8f2d9e5` |

---

## 2. Die Korrekturen im Einzelnen

### P2a — ν = µ/ρ ist ableitbar, und zwar belegt

ν ist keine unabhängige Stoffgröße. Sie gar nicht zu bilden war zu streng; sie
aus *irgendeinem* µ und *irgendeinem* ρ zu bilden wäre schlimmer — der Quotient
wäre die Viskosität einer Probe geteilt durch die Dichte einer **anderen**.

Die Ableitung ist deshalb genau dann erlaubt, wenn beide Elternwerte
nachweislich **dieselbe Flüssigkeit im selben Zustand** beschreiben. Vier
mechanisch geprüfte Bedingungen (`include/es/material_data.hpp`), die einzeln
benannt fehlschlagen:

| | Bedingung | bei 298,15 K |
|---|---|---|
| C1 | µ und ρ haben je eine **ausgewählte** Quelle | erfüllt |
| C2 | beide decken T **ohne Extrapolation** ab | erfüllt (283,15 … 373,15 K) |
| C3 | Umgebungsdruck, nicht frequenzaufgelöst | erfüllt |
| C4 | Reinheit, Wassergehalt, Probenherkunft **wörtlich** gleich | erfüllt |

Hier sind beide Eltern sogar dieselbe Publikation und dieselbe Probe
(Miranda & Santos 2025, *Fluid Phase Equilib.* **597**, 114458; 99 mass %,
0,0100 water mass %, dried by vacuum heating):

```
µ  = (36,370 ± 0,760) mPa·s      Concentric cylinders viscometry
ρ  = (1280,9  ± 1,1)  kg/m³      Vibrating tube method
ν  = (28,3941 ± 0,5938) mm²/s    ABGELEITET, nicht gemessen
```

Die Unsicherheit ist **fortgepflanzt** (quadratische Summe, 2,09 %, von µ
dominiert: 2,09 % gegen 0,086 %); die konservativere lineare Addition
(± 0,6177 mm²/s) wird mitgeliefert, weil beide aus demselben Messgang stammen.

Der Status ist `derived`, nicht `measured`: `is_direct_measurement()` ist dafür
falsch, `carries_quantitative_claim()` wahr. Quervergleich, der berichtet und
nicht weggerechnet wird: die drei direkt gemessenen Quellen geben
27,1 … 27,5 mm²/s — der abgeleitete Wert liegt **3,3 % darüber**, klein gegen
die Literaturstreuung von µ (54 %).

### P2b — die Maxwell-Traktion skaliert nicht mit γ

Bei festgehaltener Geometrie, festgehaltener angelegter Spannung und
festgehaltener Permittivitätsverteilung folgt das Feld aus der Laplace- bzw.
Poisson-Gleichung, in der γ **nirgends** vorkommt. Die Maxwell-Traktion
ε₀E²/2 ist dann ein Funktional allein dieses Feldes und **invariant**.

γ entscheidet, *welche* Form im Gleichgewicht steht — nicht, welche Kraft ein
gegebenes Feld auf eine gegebene Form ausübt.

| Größe | festgehalten | Gesetz | Exponent |
|---|---|---|---|
| Kapillardruckskala γ/a | dimensionslose Form, a | γ¹ | +1 |
| mechanische Last für dieselbe dimensionslose Form | dimensionslose Form | γ¹ | +1 |
| Spannung für dieselbe dimensionslose Form | dimensionslose Form, Geometrie | √γ | +½ |
| elektrische Bondzahl bei festem **Feld** | E, Geometrie | 1/γ | −1 |
| **Maxwell-Traktion bei festem Geometrie/U/ε** | Geometrie, U, ε-Verteilung | γ⁰ | **0** |

Die Zeile wird **nicht gelöscht**, sondern als Invariante mit Exponent 0
geführt, damit die Korrektur sichtbar bleibt. Jede Zeile trägt `recomputed=no`:
keine ist eine neu gerechnete gekoppelte Simulation. Die Tabelle liegt als
`es::gamma_scaling_rows()` in der Bibliothek, wo ein Test sie prüft — statt in
einem `printf`.

### P3 — welche Permittivität in τ_q gehört

Fünf Größen heißen bei einer leitfähigen ionischen Flüssigkeit „die
Permittivität", und sie sind keine Varianten einer Zahl: die niederfrequente
**Schein**permittivität (eine Eigenschaft der Messzelle), die intrinsische
statische ε_s (aus Mikrowellenspektren extrapoliert), die frequenzabhängige
ε\*(f), die Elektrodenpolarisation und die DC-Leitfähigkeit K.

Welche in τ_q gehört, folgt aus der Herleitung: die freie Ladung zerfällt auf
der Zeitskala τ *selbst*, ihr Spektrum liegt also bei f\* = 1/(2πτ), und dort
ist ε_r abzulesen. Da ε_r dispersiv ist, ist die Gleichung **implizit**:

```
τ = ε₀·ε_r(f*) / K ,   f* = 1/(2πτ)
```

**f\* = 2,627 GHz liegt innerhalb genau des Messbereichs, der als „nicht DC"
verworfen worden war.** Die 1–18-GHz-Daten sind nicht die falsche Frequenz,
sondern die richtige — sie umgekehrt als Gleichstromwert zu übernehmen wäre
ebenso falsch.

Fixpunktiteration auf den **Messpunkten**, keine angepasste Dispersionsfunktion,
keine Extrapolation:

| | |
|---|---|
| f\* | **2,627 GHz** (innerhalb 1–18 GHz) |
| ε_r(f\*) | 10,664 |
| τ | **6,059·10⁻¹¹ s** (18 Iterationen, Residuum 0) |
| zum Vergleich mit ε_s = 13,6 | 7,727·10⁻¹¹ s, also **+27,5 %** |

**Ein Einzelwert bleibt fehlend** — keine der vier Quellen nennt Reinheit *und*
Wassergehalt, `judge_conductor_limit()` schlägt weiterhin geschlossen fehl.
Belegt ist ein **Band** mit Empfindlichkeitsrechnung: ε_r ∈ [7,7; 13,6],
K ∈ [0,91; 1,63] S/m. Die für die Näherung schlechteste Ecke ist
(ε_hi, K_lo) — gerechnet und geprüft, nicht behauptet:

| Ecke | τ [s] | t_kap/τ | t_vis/τ |
|---|---|---|---|
| ε_lo, K_hi | 4,183·10⁻¹¹ | 4,12·10⁴ | 8,05·10⁴ |
| ε_lo, K_lo | 7,492·10⁻¹¹ | 2,30·10⁴ | 4,49·10⁴ |
| ε_hi, K_hi | 7,388·10⁻¹¹ | 2,33·10⁴ | 4,56·10⁴ |
| **ε_hi, K_lo** (schlechteste) | **1,323·10⁻¹⁰** | **1,30·10⁴** | **2,54·10⁴** |
| selbstkonsistent | 6,059·10⁻¹¹ | 2,84·10⁴ | 5,56·10⁴ |

Geforderte Schranke: 100.

**Damit ändert sich der Befund von P3.** Der Äquipotentialansatz von P3b ist für
die **statische** Form nicht mehr eine Annahme, sondern über das gesamte
begründete Band belegt — an der ungünstigsten Ecke noch um mehr als zwei
Größenordnungen über der Schranke, **ohne einen einzigen unbelegten ε_r-Wert**.
Über emittierenden Betrieb sagt das nichts.

### P6 — die Deposition regularisiert nichts

„Keine Divergenz beim Annähern" war eine Aussage über **ein festes Netz**, und
der Zusatz fehlte. Fünf Fragen, getrennt beantwortet:

| | Frage | bei festem Netz | bei h → 0 |
|---|---|---|---|
| (a) | Ladungserhaltung der Deposition | exakt (2,4·10⁻¹⁶) | bleibt exakt auf **jeder** Stufe |
| (b) | Fremdfeld bei **festem** Abstand | endlich | **konvergiert**, Ordnung 2,05 / 1,91 |
| (c) | Selbstpotential am Teilchen | endlich | **divergiert** |
| (d) | Breite der deponierten Wolke | eine Zelle | fällt **wie das Netz** gegen null |
| (e) | Anteil des Selbstfeldes | — | fällt wie N^−0,82 mit der Makropartikelzahl |

Zu (c), und das Gesetz hängt vom Ort ab:

| Lage | über fünf Stufen (h: 10⁻⁶ → 6,25·10⁻⁸ m) | angepasstes Gesetz |
|---|---|---|
| abseits der Achse, r = 0,35 R | 1,938 → 3,073 V, Faktor **1,59** | 0,409·ln(1/h), Rest 1,2·10⁻⁴ |
| auf der Achse, r = 0 | 32,4 → 524,6 V, Faktor **16,18** | h^−1,004, Rest 3,0·10⁻⁴ |

Genau das physikalisch Erwartete: abseits der Achse ist ein Makropartikel ein
**Ring** mit logarithmisch singulärem Eigenpotential; auf der Achse entartet er
zur **Punktladung** mit 1/d. Der frühere Header-Kommentar behauptete „wie 1/h" —
das ist nur das Achsengesetz.

**Die drei Kandidaten, mit Urteil und mit der Messung dahinter**
(`es::pic_options()`, in der Bibliothek und im Test geprüft):

| Variante | Urteil |
|---|---|
| **Selbstfeldabzug** | `implemented` — eine Ladung übt auf sich selbst keine Kraft aus; die diskrete Aufgabe ist **linear**, also exakt abziehbar statt dämpfbar |
| feste physikalische Formbreite | `rejected_free_parameter` — Messung (d) zeigt, dass in der Wolke nur h steckt; sie wäre ein frei gewählter Glättungsparameter |
| skalierte Makropartikelzahl | `measured_not_implemented` — Eigenschaft einer Schleife, die es hier nicht gibt |

Gemessen: ein *einzelnes* Teilchen im sonst leeren, geerdeten Kasten spürt ohne
Abzug 1,13·10⁵ … 1,42·10⁵ V/m, obwohl nichts anderes darin ist — der Betrag
hängt allein von seiner Lage in der Zelle ab. Mit Abzug spürt es **exakt null**.
Überlagerungsfehler 1,6·10⁻¹⁵. **Preis: eine zusätzliche Lösung je Teilchen.**

Die PIC-Schleife bleibt **blockiert, aus zwei unabhängigen Gründen**: P5 hat
keine Quelle, und der Abzug skaliert nicht.

### P9 — Vergleichbarkeit ist nicht Validierung

Die frühere Matrix trug *eine* Vergleichbarkeit je Zeile, und die Abbildung
färbte die Zeile danach. Der **Gesamtstrom** ist direkt vergleichbar und von
diesem Projekt überhaupt nicht rechenbar — er erschien grün.

Jede Zeile trägt jetzt sechs unabhängige Urteile (`yes`/`partial`/`no`/`n/a`):
`comparable_geometry`, `implemented`, `converged`, `comparable_with_data`,
`validated`, `blocked` + Grund. **Grün heißt nur, dass genau diese eine Achse
erfüllt ist.**

Die Invariante steht im Code (`es::inconsistency()`): `validated=yes` verlangt
`implemented=yes` UND `converged` ∈ {yes, n/a} UND `comparable_with_data=yes`
UND nicht blockiert. Der Test baut zusätzlich Zeilen, die sie verletzen, und
prüft, dass sie abgefangen werden.

| Achse | von 13 erreicht |
|---|---|
| vergleichbar (ganz oder nach Reduktion) | **9** |
| implementiert | **8** |
| numerisch konvergiert | **0** |
| mit Messdaten vergleichbar | **6** |
| **tatsächlich validiert** | **0** |
| blockiert | **5** |

Der Abstand zwischen der ersten und der vorletzten Zahl ist der ganze Punkt der
Tabelle. Dass nichts validiert ist, ist eine **geprüfte** Aussage.

**Der Importvertrag** unterscheidet jetzt harte von weichen Bedingungen. Fehlende
Einheit, Fundstelle oder Geometrieart lehnen den Satz als GANZES ab. Eine in der
Publikation **nicht angegebene Unsicherheit** ist dagegen ein eigener Zustand
(`NotReported` / `OkUncertaintyNotReported`): der Punkt wird importiert und
archiviert, darf qualitativ dargestellt werden, und `usable_quantitatively()`
ist für ihn falsch. Ein gemischter Satz bleibt ganz und wird getrennt gezählt.

---

## 3. Die Historienreparatur

`713ab86` (P0) nannte in `CMakeLists.txt` drei Einträge — `src/feed.cpp`,
`es_feed`, `test_feed` —, deren Dateien erst mit `1cb6108` (P1) hinzukamen. Der
Commit ließ sich damit **nicht eigenständig konfigurieren**.

Genau zwei Commits wurden inhaltlich angefasst: P0 verliert die drei Einträge,
P1 bekommt sie. Alle übrigen sind **Byte für Byte** wiedergespielt; der Baum am
Ende des wiedergespielten Bereichs ist identisch mit dem vorherigen. Nur die
Nachrichten der Artefakt-Commits nennen jetzt den neuen Hash und führen den
alten mit.

**Geprüft, nicht behauptet:** jeder der elf Code-Commits wurde einzeln in einem
eigenen Worktree ausgecheckt, konfiguriert, gebaut und mit seinem Phasentest
geprüft.

| Commit | Punkt | configure | build | Phasentest | Tests am Commit |
|---|---|---|---|---|---|
| `0677e3b` | P3b | ok | ok | `test_electrocapillary` | 16 |
| `117dc28` | P0 | ok | ok | `test_load_projection` | 17 |
| `205ab7a` | P1 | ok | ok | `test_feed` | 18 |
| `beebc38` | P2 | ok | ok | `test_material_data` | 19 |
| `790be1b` | P3 | ok | ok | `test_transport` | 20 |
| `60e71ab` | P4 | ok | ok | `test_surface_kinematics` | 21 |
| `f8d5237` | P5 | ok | ok | `test_emission_contract` | 22 |
| `d1b59d8` | P6 | ok | ok | `test_space_charge` | 23 |
| `a81696d` | P7 | ok | ok | `test_particle_transport` | 24 |
| `03034fa` | P8 | ok | ok | `test_cone_jet_contract` | 25 |
| `b0731f8` | P9 | ok | ok | `test_validation` | 26 |

### Alt → Neu

| alt | neu | | alt | neu |
|---|---|---|---|---|
| `713ab86` | `117dc28` | | `baaf6fd` | `f8d5237` |
| `85e7096` | `e8b329f` | | `114768e` | `8a3eeec` |
| `1cb6108` | `205ab7a` | | `8de4558` | `d1b59d8` |
| `3ffd040` | `2e03b40` | | `93b3be6` | `5030d80` |
| `4c4bbfe` | `beebc38` | | `cc8f52f` | `a81696d` |
| `42e24d4` | `ee94fed` | | `3e97de0` | `2376455` |
| `2a7205d` | `790be1b` | | `4274453` | `03034fa` |
| `26e8a27` | `b973b38` | | `be884ab` | `1decc56` |
| `c4ecbb7` | `60e71ab` | | `e8713ee` | `b0731f8` |
| `5067cea` | `b3fe7d9` | | `f6fc6f2` | `fcef5e6` |
| | | | `134b620` | `1746c83` |
| | | | `ede1508` | `f2c5270` |

`0677e3b` und `6dc833b` sind unverändert; die Reparatur beginnt erst danach.

**Kein Rest.** Ein Provenienzstempel, der auf einen Commit außerhalb der eigenen
Historie zeigt, ist keine Provenienz. Nach der Reparatur traf das auf P0, P1,
P4, P5 und P7 zu; alle fünf sind aus sauberem Arbeitsbaum neu erzeugt. Geprüft
wird das mechanisch: für **jedes** Ergebnisverzeichnis liegt `head_commit` mit
`git merge-base --is-ancestor` in der Historie, und jedes trägt
`working_tree_dirty=no` und `releasable=yes`.

Die **Zahlen** sind dabei unverändert — der Code dieser Punkte ist von den
Korrekturen unberührt. Bei P1 ist **jede CSV byteweise identisch**; bei P0
reproduziert der einstündige Lauf die Befunde (1000 V `NotInAsymptoticRange`,
`total_force` 6,07 %; 1400 V 4,59 … 5,67 %, alles `DiscretizationNotConverged`).

---

## 4. Status je Punkt

| Punkt | Status | Code | Artefakte | Blocker |
|---|---|---|---|---|
| 0 — P3b numerisch bereinigen | `validated_subset` + `DiscretizationNotConverged` | `117dc28` | `8040125` | 1-%-Ziel **verfehlt** (4,6–6,1 %); bei 1000 V nicht einmal asymptotisch |
| 1 — Druck- und Zulaufmodell | `implemented` | `205ab7a` | `e46c941` | keiner; `p_reservoir` und `Q` bleiben Eingaben |
| 2 — Reale Stoffdaten | `implemented` (γ, ρ, µ, K) + **`derived`** (ν) + `MissingMaterialData` (ε_r) | `beebc38`, `2b95a94` | `572f1e8` | ε_r ohne Quelle mit Methode + Reinheit + Wassergehalt |
| 3 — Endliche Leitfähigkeit | `validated_subset` | `790be1b`, `bfa2c61` | `fa0e090` | ε_r-**Einzelwert** fehlt weiterhin; keine gekoppelte finite-conductivity-Meniskusrechnung |
| 4 — Zeitabhängige freie Oberfläche | `infrastructure_only` | `60e71ab` | `b3fe7d9` | kein Feld mit freier Oberfläche, keine tangentiale Traktion q_s·E_t |
| 5 — Ionenemission | `blocked` | `f8d5237` | `8a3eeec` | Ratengleichung an keiner Primärquelle geprüft; kein belegtes ΔG |
| 6 — Raumladung und PIC-Grundlage | `validated_subset` | `d1b59d8`, `4c3b83b` | `593cecd` | PIC-Schleife **blockiert**, zweifach: P5 hat keine Quelle, der Abzug skaliert nicht |
| 7 — Ionenbahnen | `validated_subset` | `a81696d` | `2376455` | keine physikalische Teilchenquelle → Transportantwort, keine Stromvorhersage |
| 8 — Cone-Jet | `blocked` | `03034fa` | `1decc56` | fünf von sieben Teilmodellen fehlen; Gañán-Calvo nicht im Volltext erreichbar |
| 9 — 3D und Validierung | `infrastructure_only` | `b0731f8`, `3bef623` | `16cdd09` | kein 3D-Netz, kein 3D-Löser; **keine Messdaten importiert** |

---

## 5. Was weiterhin NICHT gerechnet wird

Keine Emission. Keine selbstkonsistente Emissions-PIC-Schleife. Keine
Zweiphasenströmung, kein Jet, kein Tropfenzerfall. Kein dynamischer Meniskus.
Keine gekoppelte finite-conductivity-Meniskusrechnung. Keine Stabilitätsaussage,
kein Taylor-Kegel-Onset. Kein 3D-Netz, kein 3D-Löser, kein 3D-Ergebnis. Keine
importierten Messdaten. Keine Stromvorhersage. **Keine netzunabhängige
Regularisierung des Selbstfeldes** — der exakte Abzug entfernt die scheinbare
Selbstkraft, macht das Selbstpotential aber nicht netzunabhängig.

**Nichts ist validiert**, und das ist eine geprüfte Aussage: 0 von 13 Größen.

---

## 6. Reproduktion

```sh
cmake -S . -B build -G Ninja && cmake --build build && (cd build && ctest)

./build/es_p3b_audit     examples/device_p1.cfg examples/electrocapillary_p3b.cfg \
                         results/2026-08-29_p0_p3b_audit    meta.commit=$(git rev-parse HEAD)
./build/es_feed          examples/device_p1.cfg examples/feed_p1.cfg \
                         results/2026-08-29_p1_pressure_budget meta.commit=$(git rev-parse HEAD)
./build/es_material      results/2026-08-29_p2_material_data  meta.commit=$(git rev-parse HEAD)
./build/es_transport     results/2026-08-29_p3_transport      meta.commit=$(git rev-parse HEAD)
./build/es_kinematics    results/2026-08-29_p4_kinematics     meta.commit=$(git rev-parse HEAD)
./build/es_emission_audit results/2026-08-29_p5_emission_audit meta.commit=$(git rev-parse HEAD)
./build/es_space_charge  results/2026-08-29_p6_space_charge   meta.commit=$(git rev-parse HEAD)
./build/es_trajectories  results/2026-08-29_p7_trajectories   meta.commit=$(git rev-parse HEAD)
./build/es_cone_jet      results/2026-08-29_p8_cone_jet       meta.commit=$(git rev-parse HEAD)
./build/es_validation    results/2026-08-29_p9_validation     meta.commit=$(git rev-parse HEAD)
```

Danach je Ordner das zugehörige `python/plot_*.py`. Der Punkt-0-Lauf braucht
rund eine Stunde (gekoppelte Netzstufe 4 bei 1000 V und 1400 V, Kanten-Gate bis
Stufe 5, 3 GiB Bandfaktorisierung); der P6-Lauf einige Minuten (die
Selbstfeldstudie löst über fünf Netzstufen bis 321×641); alle anderen laufen in
Sekunden.

**Ein Ordner nach dem anderen**, mit sonst sauberem Arbeitsbaum: `provenance.py`
ignoriert nur das Zielverzeichnis und stempelt sonst `working_tree_dirty=yes`.

Sanitizerlauf:

```sh
cmake -S . -B build_asan -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer -g -O1'
cmake --build build_asan -j2
(cd build_asan && ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1 ctest)
```

Die Stoffdatentabelle wird von `python tools/fetch_material_data.py` neu erzeugt;
ein leerer `git diff` danach ist die Reproduzierbarkeitsprüfung. Dieser Lauf hat
sie **nicht** ausgeführt — sie braucht die ILThermo-Abfrage.

---

## 7. Die inhaltliche Prüfung der neun Abbildungen

Neun Abbildungen waren nach dem P2-Nachzug (`5d8caba`) neu erzeugt und noch
nicht im Endzustand vorgelegt worden: P3 Abb. 1–3, P4 Abb. 1–2, P5 Abb. 1,
P7 Abb. 1–2 und P8 Abb. 1. Sie sind einzeln nachgerechnet, nicht nur geöffnet.

### Nachgerechnet gegen die CSV

| Abbildung | geprüfte Größe | gerechnet | in der CSV |
|---|---|---|---|
| P3 Abb. 2/3 | τ = ε₀·ε_r(f\*)/K | 6,058912579·10⁻¹¹ s | 6,058912582·10⁻¹¹ s |
| P3 Abb. 3 | f\* = 1/(2πτ) | 2,626790549·10⁹ Hz | 2,627·10⁹ Hz |
| P3 Abb. 3 | Bandecke (ε_lo, K_hi) | 4,182653139·10⁻¹¹ s | 4,182653139·10⁻¹¹ s |
| P3 Abb. 3 | schlechteste Ecke, t_kap/τ | 1,301159·10⁴ | 1,301159355·10⁴ |
| P4 Abb. 2 | exakte Volumenänderung e^{3αT}−1 | 3,4816890703 | 3,481689070 |
| P5 Abb. 1 | E\* = 4πε₀ΔG²/e³ bei 1,09 eV | 8,2521·10⁸ V/m | 8,250897578·10⁸ V/m |
| P5 Abb. 1 | b = ΔG/kT bei 0,80 eV | 31,1374 | 3,113739560·10¹ |
| P7 Abb. 2 | Transmission 29/41 | 0,7073170732 | 7,073170732·10⁻¹ |
| P7 Abb. 2 | Energiegewinn 999 V · e | 1,600574457·10⁻¹⁶ J | 1,600574457·10⁻¹⁶ J |
| P8 Abb. 1 | r\*_lo/r\* = (τ_lo/τ)^{2/3} | 0,7810966516 | 0,7810966516 |
| P8 Abb. 1 | r\*_hi/r\* = (τ_hi/τ)^{2/3} | 1,6833239649 | 1,6833239646 |

**Querschnitt über die Punkte:** P8 trägt exakt dasselbe τ, ε_r(f\*) und f\*
wie P3 — die P3-Korrektur ist also tatsächlich durchgezogen und nicht nur in
P3 vermerkt. Im Datenteil von `diagnosis.csv` steht **kein** `nan` mehr.

### Achsen, Einheiten, Legenden, Bildunterschriften

Gelesen, nicht überflogen. Jede Achse trägt Größe und Einheit oder ist
ausdrücklich als dimensionslos beschriftet; die Legenden nennen die Kurven, die
gezeichnet sind. Die beiden Stellen, an denen eine überholte Aussage saß, sind
im Bild verschwunden: der Paneltitel von P3 Abb. 2 nennt jetzt
`τ = ε₀ε_r(f*)/K` mit dem selbstkonsistenten τ statt eines unbelegten ε_r, und
P8 Abb. 1 zeichnet die Ladungsrelaxation als Band statt als „nicht berechenbar".

### Statusaussagen gegen den aktuellen Code

| Punkt | im Bild und in `meta.txt` | im Code |
|---|---|---|
| P3 | `validated_subset` | `apps/es_transport.cpp:444` |
| P4 | `infrastructure_only` | `apps/es_kinematics.cpp:175` |
| P5 | `blocked` | `apps/es_emission_audit.cpp:147` |
| P7 | `validated_subset` | `apps/es_trajectories.cpp:259` |
| P8 | `blocked` | `apps/es_cone_jet.cpp:180` |

P5 Abb. 1 zeigt außerdem jeden Pfad einzeln mit `Disabled` bzw.
`MissingEmissionParameters` — genau die sechs Zeilen aus `contract.csv`.
P8 Abb. 1 zeigt fünf von sieben Teilmodellen als fehlend, was der Zeile in
Abschnitt 4 entspricht.

### Mechanisch

PNG-Signatur, `IEND`-Abschluss, Dateigröße und Abmessungen aller neun geprüft
(2030×957 bis 2523×928, 191 861 bis 297 170 Bytes). Provenienz: `42f4d732add3`
(P3), `ee6b743d2ab1` (P4), `91f254c5c768` (P5), `bed39d795ee3` (P7),
`7b683d4fb444` (P8) — alle fünf per `git merge-base --is-ancestor` in dieser
Historie, alle mit `working_tree_dirty=no` und `releasable=yes`.

---

## 8. Zum beweglichen Netz — was geplant und was nicht implementiert ist

Neu: [`docs/21_moving_mesh_plan.md`](../../docs/21_moving_mesh_plan.md).
Es trennt drei Dinge, die bisher nebeneinanderstanden, und legt fest, was vor
einer Taylor-Cone-Aussage erfüllt sein muss:

* Das P3b-Netz hat eine **feste Topologie** mit **beweglichen
  Knotenkoordinaten** — keine Knoteneinfügung, keine Neuvernetzung.
* Seine Bewegung sucht ein **statisches Gleichgewicht**; die Iterationsschritte
  sind keine Zeitschritte und die Zwischenformen haben keine physikalische
  Bedeutung.
* P4 enthält **nur die Kinematik** für ein vorgeschriebenes Geschwindigkeitsfeld
  und ist mit P3b nicht gekoppelt.
* Stark zugespitzte beziehungsweise Taylor-Cone-nahe Formen brauchen **adaptive
  Verfeinerung** und voraussichtlich **Neuvernetzung**.
* Vor einer Taylor-Cone-Aussage müssen **Apexkrümmungsradius, Apexhöhe,
  integrierte Maxwell-Kraft und Netzqualität** unter Verfeinerung *und* unter
  Neuvernetzung konvergieren.
* Ein Abbruch der Fortsetzung oder des Netzverfahrens ist **kein
  Taylor-Cone-Onset**, sondern eine Eigenschaft dieses Lösers und dieses Netzes.

**In diesem Lauf wurde keine adaptive Neuvernetzung implementiert**; der Code
ist gegenüber `ea2d5a1` unverändert.
