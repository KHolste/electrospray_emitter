# 16 — P5: Ionenemission — Vertrag und Blocker

**Stand: 2026-08-29.** Status: **`blocked`**. Das Emissionsmodell ist
implementiert als *Vertrag*, standardmäßig **abgeschaltet**, und es liefert auf
keinem Pfad eine Zahl.

Code: `include/es/emission_contract.hpp`, `src/emission_contract.cpp`,
`apps/es_emission_audit.cpp`, `tests/test_emission_contract.cpp`,
`python/plot_emission_audit.py`.
Ergebnisse: `results/2026-08-29_p5_emission_audit/`.

---

## 16.1 Literaturaudit — was gesucht und was gefunden wurde

### Die Gleichung

Die Ionenverdampfungsrate von Iribarne und Thomson wird durchgängig zitiert als

$$
j = \frac{k_BT}{h}\,\sigma_s\,
\exp\!\left[-\frac{\Delta G - \sqrt{e^3E/(4\pi\varepsilon_0)}}{k_BT}\right],
\qquad \sigma_s = \varepsilon_0 E,
$$

und genau diese Form rechnet `src/emission.cpp` bereits.

**Primärquellen:** Iribarne & Thomson (1976), *J. Chem. Phys.* **64**, 2287–2294;
Thomson & Iribarne (1979), *J. Chem. Phys.* **71**, 4451–4463.

**Keiner der beiden Volltexte war in diesem Lauf erreichbar.** Damit sind der
**Vorfaktor**, die genaue Definition von $\sigma_s$ und der angegebene
**Gültigkeitsbereich** nicht an der Quelle geprüft.

**Was gelesen wurde:** Wiley, „A Brief Overview of the Mechanisms Involved in
Electrospray Mass Spectrometry", Kapitel 1 (Leseprobe, PDF geöffnet und
Textinhalt extrahiert). Der Abschnitt zum Ion Evaporation Model zitiert beide
Primärarbeiten mit vollständiger Fundstelle, beschreibt das Modell qualitativ,
**druckt aber keine Gleichung** und stellt fest:

> „the Ion Evaporation Model is experimentally well supported for small ions of
> the kind that one encounters in inorganic and organic chemistry. However, when
> the ions beco[me] …"

— also ein **ausdrücklicher Vorbehalt für größere Ionen**, und genau darum
handelt es sich bei den Clustern einer ionischen Flüssigkeit. Zusätzlich fand
die Recherche die Aussage, dass die für die theoretische Rate nötigen
thermochemischen Daten „of insufficient accuracy" seien.

**Urteil: `EquationNotValidated`.**

### Die Parameter

**Massen: belegt.** Die Molmasse der Verbindung, 197,97 g/mol, steht im
IoLiTec-Datenblatt IL-0006 (Seite 1/4) und im ILThermo-Komponentendatensatz;
die Aufteilung auf C₆H₁₁N₂⁺ (111,17 g/mol) und BF₄⁻ (86,81 g/mol) ist die
Formelzusammensetzung.

**Aktivierungsbarriere: NICHT belegt.** Für EMI-BF4 wurde keine Quelle mit
einer Zahl gefunden. Der Wert 1,09 eV in `src/fluid.cpp` hat keine Quelle; der
Header dieser Datei nennt ihn selbst „the least certain quantity here by far"
mit einer Spanne von etwa 1,0 bis 1,4 eV je nach Auswertungsmethode.

**Was das bedeutet, als Zahl.** $j$ hängt **exponentiell** von $\Delta G$ ab.
Bei 298 K und festem Feld ändert sich die Rate über die Spanne 1,0 → 1,4 eV um
den Faktor **5,8·10⁶**. Ein Modell mit unbelegter Barriere sagt keinen Strom
voraus; es berichtet den Parameter.

**Urteil: `MissingEmissionParameters`.**

## 16.2 Der Vertrag

Jede physikalische Eingabe ist explizit:

| Größe | Feld | fail-closed? |
|---|---|---|
| Spezies und Name | `EmittedSpecies::name` | — |
| **Ladungsvorzeichen** | `charge_number` (signiert) | `PolarityMismatch` |
| Masse | `mass` + `mass_source` | `MissingEmissionParameters` |
| Aktivierungsbarriere | `activation_barrier` + `barrier_source` | `MissingEmissionParameters` |
| Temperatur | `EmissionModel::temperature` | `MissingEmissionParameters` |
| Feldrichtung | Vorzeichen von `E_n` | `PolarityMismatch` |
| Gleichung geprüft | `equation_validated` | `EquationNotValidated` |
| Modell aktiv | `enabled` (Vorgabe **false**) | `Disabled` |

**Positive und negative Polarität getrennt.** Eine Spezies trägt ihr
Vorzeichen; ein Kation wird nur emittiert, wenn das äußere Normalfeld positiv
ist, ein Anion nur, wenn es negativ ist. Der Prototyp benutzte $|E_n|$ und
Kationenmassen für beide Vorzeichen und berichtete deshalb identische Ströme
für physikalisch verschiedene Fälle.

**Kein Rückgriff auf die alten quellenlosen Werte.** Dieses Modul liest
`es::Fluid` nicht und hat keinen Vorgabewert für irgendeinen Parameter.

## 16.3 Kern und Modell sind getrennt

Die *mathematische Form* ist als **reine Funktion** ihrer Argumente
implementiert (`iribarne_thomson_rate`, `schottky_barrier_lowering`,
`barrier_free_field`, `iribarne_thomson_dimensionless`). Sie lässt sich ohne
jede Stoffangabe prüfen — Dimensionen, Grenzwerte, Monotonie in $E$ und in
$\Delta G$, das $\sqrt E$-Gesetz der Barrierensenkung —, und genau das tut der
Test.

**Die Form zu prüfen ist nicht dasselbe wie das Modell zu validieren.** Ob
diese Funktion die Emission einer ionischen Flüssigkeit beschreibt, ist die
offene Frage, und der Kern liefert deshalb nie ein Ergebnis über den
Modell-Einstiegspunkt.

Gemessene Eigenschaften des Kerns:

| Prüfung | Ergebnis |
|---|---|
| $G(0)=0$ | exakt |
| $G(E)=\sqrt{e^3E/(4\pi\varepsilon_0)}$ | exakt |
| $G(4E)=2G(E)$ | exakt |
| $G(E^*)=\Delta G$ mit $E^*=4\pi\varepsilon_0\Delta G^2/e^3$ | $<10^{-12}$ |
| $j(E^*)$ ist genau der Vorfaktor $(k_BT/h)\varepsilon_0E^*$ | $<10^{-12}$ |
| Monotonie in $E$ über zehn Dekaden | streng |
| Monotonie in $\Delta G$ | streng fallend, Faktor $5{,}8\cdot10^6$ über 1,0→1,4 eV |
| $T=0$, $\Delta G=0$ | NaN, nicht null |
| $E=0$ | **exakt** null ($\sigma_s=\varepsilon_0E$) |
| dimensionslose Form = skalierte Form | $<10^{-12}$ |

## 16.4 Ein im Test gefundener Fehler

Die erste Fassung von `synthetic_complete_model()` gab ein `EmissionModel`
zurück, das auf einen **gemeinsamen statischen Puffer** zeigte. Das Anlegen
eines zweiten Modells überschrieb damit die Spezies des ersten, und der Test,
der beide Polaritäten vergleichen sollte, verglich ein Modell mit sich selbst.

**Gefunden hat es der Test, weil er beide Ströme ausdruckt** und einer davon
null war. Behoben durch Besitz: `SyntheticEmissionModel` hält seine Spezies
selbst. Der Test prüft jetzt ausdrücklich, dass beide Polaritäten einen Strom
liefern — was mit dem geteilten Puffer nicht der Fall war.

## 16.5 Drei weitere Verbote

* **Keine Cone-Jet-Formel als Ionenemission.** Das ist andere Physik; der alte
  Cone-Jet-Block bleibt deaktiviert (P8).
* **Kein Strom aus einem Perfect-Conductor-Feld als Vorhersage.** Ohne die
  finite-conductivity-Rückwirkung (P3/P4 zeigen, dass sie fehlt) ist das
  Oberflächenfeld am Emissionsort nicht das Feld, das die Emission sähe.
* **Keine Betriebsprognose in Abbildungen.** Die einzige Abbildung dieses
  Punktes ist eine Statuskarte plus eine **dimensionslose** Sensitivität; ein
  absoluter Strom kommt nirgends vor.

## 16.6 Was zu tun wäre, um den Blocker zu lösen

1. Iribarne & Thomson (1976) und Thomson & Iribarne (1979) im Volltext lesen und
   Vorfaktor, $\sigma_s$-Definition und Gültigkeitsbereich festhalten;
   `equation_validated` erst danach setzen.
2. Eine Primärquelle für $\Delta G$ von EMI-BF4 mit Methode, Spezies (Monomer,
   Dimer, Cluster) und Unsicherheit beschaffen. Ohne Unsicherheit ist der Wert
   wegen der exponentiellen Empfindlichkeit wenig wert.
3. Die Speziesverteilung festlegen: Monomer, Dimer und höhere Cluster haben
   verschiedene Massen **und** verschiedene Barrieren, und die Verteilung setzt
   $q/m$.
4. Erst dann `enabled = true`, und auch dann nur mit einer Feldquelle, die die
   endliche Leitfähigkeit enthält.
