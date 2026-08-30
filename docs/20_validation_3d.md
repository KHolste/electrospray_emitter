# 20 — P9: Achsensymmetrie, 3D-Erweiterung und Validierung

**Stand: 2026-08-29.** Status: **`infrastructure_only`**. Es gibt **kein
3D-Netz, keinen 3D-Löser und kein 3D-Ergebnis**. Geliefert sind die Trennung,
die Rotationsreferenz, der Importvertrag und die Validierungsmatrix.

Code: `include/es/validation.hpp`, `src/validation.cpp`,
`apps/es_validation.cpp`, `tests/test_validation.cpp`,
`python/plot_validation.py`.
Ergebnisse: `results/2026-08-29_p9_validation/`.

---

## 20.1 Die Geometrieart reist mit dem Ergebnis

`LabelledResult` trägt Wert, Größe, Einheit und `GeometryKind`. Die Art wird im
Konstruktor gesetzt und es gibt **keinen Setter**: ein Ergebnis kann
nachträglich nicht umetikettiert werden. `value_as_three_dimensional()` ist das
einzige Gatter, durch das eine 3D-Behauptung muss, und es wirft für alles außer
einer echten 3D-Rechnung.

Das ist keine Pedanterie. Ein achsensymmetrisches Modell sieht genau das nicht,
worin sich eine reale Anordnung von ihm unterscheidet; die Fehlerart ist eine
Zahl, die richtig aussieht und eine andere Frage beantwortet.

Drei Arten, und die dritte ist die wichtige:

| Art | Bedeutung |
|---|---|
| `Axisymmetric` | die Meridianhalbebene. **Alles**, was dieses Projekt rechnet |
| `ThreeDimensional` | eine echte 3D-Geometrie. **Nichts hier erzeugt eine** |
| `RevolvedAxisymmetric` | ein achsensymmetrisches Ergebnis, als 3D-Feld ausgewertet. **Immer noch achsensymmetrische Physik** |

Die dritte existiert, damit eine Rotationsreferenz nicht als 3D-Rechnung
ausgegeben werden kann.

## 20.2 Die Rotationsreferenz

Ein achsensymmetrisches Feld $\varphi(r,z)$ **ist** das 3D-Feld
$\varphi(\sqrt{x^2+y^2},z)$, und ein achsensymmetrischer Vektor $(E_r,E_z)$
**ist** $(E_r\cos\theta,\,E_r\sin\theta,\,E_z)$. Beides ist trivial wahr — und
beides ist im Code prüfenswert, weil die Stelle, an der man es gewöhnlich falsch
macht, das **Integral** ist.

Geprüft wird deshalb dieselbe Größe zweimal:

* achsensymmetrisch als $\int f\,2\pi r\,\mathrm ds$;
* als **explizite 3D-Quadratur** $\iint f\,r\,\mathrm d\theta\,\mathrm ds$ über
  $n$ Azimute.

| Größe | Unterschied | geschlossene Form |
|---|---|---|
| Fläche der Halbkugel | $\le10^{-12}$ | $2\pi R^2$, getroffen auf $8\cdot10^{-8}$ |
| $\int (z/R)\,\mathrm dA$ | $\le10^{-12}$ | $\pi R^2$, getroffen auf $1{,}3\cdot10^{-7}$ |
| Rotationsvolumen | $\le10^{-12}$ | $\tfrac23\pi R^3$, getroffen auf $1{,}5\cdot10^{-7}$ |

**Der verbleibende Unterschied ist die Summationsrundung**, nicht ein
Modellunterschied: die 3D-Quadratur addiert 720 × 2000 Terme, wo die
achsensymmetrische Form 2000 addiert. Und die **Azimutzahl ändert nichts** — 3
und 997 Azimute geben dasselbe, weil die Mittelpunktregel für einen
azimutunabhängigen Integranden exakt ist. Wäre sie es nicht, wäre der Integrand
nicht azimutunabhängig, also das Feld nicht achsensymmetrisch.

**Damit ist die $2\pi r$-Wichtung, die jedes Integral dieses Projekts trägt, als
das dreidimensionale Integral nachgewiesen** — und nicht bloß behauptet.

Das ist **kein 3D-Löser**: es wird nichts Neues gelöst, es wird eine Wichtung
geprüft.

## 20.3 Der Importvertrag für Messdaten

Ein Messwert ohne Einheit ist eine Zahl; ohne Fundstelle ein Gerücht; und ohne
die Geometrie, zu der er gehört, mit nichts vergleichbar. **Diese drei sind
harte Bedingungen**, und ein Satz, der eine davon bricht, wird als GANZES
abgelehnt — einen Punkt stillschweigend wegzulassen wäre ein Vergleich mit einem
anderen Datensatz.

| Feld | fehlt → | Härte |
|---|---|---|
| `unit` | `MissingUnit` | **hart**, Satz abgelehnt |
| `provenance` | `MissingProvenance` | **hart**, Satz abgelehnt |
| `geometry_stated` | `MissingGeometryKind` | **hart**, Satz abgelehnt |
| zwei Einheiten für dieselbe Größe | `UnitMismatch` | **hart**, Satz abgelehnt |
| `uncertainty` bzw. `uncertainty_type` | `OkUncertaintyNotReported` | **nicht hart** — siehe 20.3a |

### 20.3a Eine nicht angegebene Unsicherheit ist kein Defekt des Datensatzes

Eine frühere Fassung dieses Vertrags behandelte die fehlende Unsicherheit wie
eine fehlende Einheit und lehnte den Punkt hart ab. **Das war falsch und es warf
echte Messungen weg.**

Eine Publikation, die einen Strom von 210 nA ohne Fehlerbalken berichtet, hat
keinen *kaputten* Datensatz erzeugt, sondern einen *unvollständigen*. Die beiden
Fälle sind nicht dasselbe:

* eine fehlende **Einheit** macht die Zahl unlesbar — man weiß nicht, was da
  steht;
* eine fehlende **Unsicherheit** macht die Zahl unbrauchbar für einen
  *quantitativen* Vergleich — aber die Messung selbst existiert, ist zitierbar
  und trägt eine Aussage.

Der Vertrag führt deshalb den ausdrücklichen Zustand
`UncertaintyType::NotReported` und den Importstatus
`ImportStatus::OkUncertaintyNotReported`:

* der Punkt wird **importiert und archiviert**;
* er darf **qualitativ dargestellt** werden — mit seinem Status sichtbar;
* `usable_quantitatively()` ist für ihn **falsch**. Jede Abweichung, jedes
  Chi-Quadrat und jedes Bestanden/Durchgefallen muss genau diese Funktion
  abfragen, nicht `is_usable()`.

Ein **gemischter Satz** bleibt ganz: `ImportResult` führt `n_quantitative` und
`n_qualitative_only` getrennt, der Satzstatus verschweigt nicht, dass ein Punkt
keine Unsicherheit trägt, und **jeder Punkt trägt seinen eigenen Status** — der
Satz ist nicht homogen, und so zu tun wäre derselbe Fehler an einer neuen
Stelle.

Zusätzlich aufgenommen: `coverage_factor` $k$ und `conditions` (Temperatur,
Polarität, Volumenstrom) — ohne die Bedingungen ist auch ein vollständiger
Messwert nicht zuordenbar.

## 20.4 Die Validierungsmatrix: sechs Fragen, nicht eine

**Was daran falsch war.** Die frühere Matrix trug *eine* Vergleichbarkeit je
Zeile, und die Abbildung färbte die Zeile danach. Eine Größe, die im Prinzip
vergleichbar, aber blockiert oder nicht konvergiert ist, erschien damit **grün**
— und grün liest sich als Erfolg. Der klarste Fall ist der **Gesamtstrom**: er
ist zwischen einer achsensymmetrischen Rechnung und einem realen Gerät *direkt*
vergleichbar, und dieses Projekt kann ihn **überhaupt nicht rechnen**, weil P5
blockiert ist.

Vergleichbarkeit und Validierung sind verschiedene Fragen. Jede Zeile trägt
jetzt **sechs unabhängige Urteile**:

| Achse | Frage |
|---|---|
| `comparable_geometry` | Lässt sich die Größe zwischen achsensymmetrisch und 3D überhaupt vergleichen? |
| `implemented` | Rechnet dieses Projekt sie? |
| `converged` | Ist das numerische Ergebnis nach einem **vorab** festgelegten Kriterium konvergiert? |
| `comparable_with_data` | Ließe sie sich mit einer Messung vergleichen? |
| `validated` | Ist sie **tatsächlich** mit Messdaten verglichen worden und hat innerhalb der angegebenen Unsicherheiten übereingestimmt? |
| `blocked` + Grund | Ist sie blockiert, und wodurch? |

Die Werte je Achse sind `yes` / `partial` / `no` / `n/a`. **Grün heißt in der
Abbildung nur, dass genau diese eine Achse erfüllt ist** — nie, dass die Größe
validiert wäre.

### Die Invariante steht im Code, nicht in der Bildunterschrift

`es::inconsistency()` prüft je Zeile:

$$
\texttt{validated}=\text{yes}\;\Longrightarrow\;
\texttt{implemented}=\text{yes}\;\wedge\;
\texttt{converged}\in\{\text{yes},\text{n/a}\}\;\wedge\;
\texttt{comparable\_with\_data}=\text{yes}\;\wedge\;
\neg\,\texttt{blocked}
$$

sowie: blockiert genau dann, wenn ein Grund genannt ist. Der Test baut
zusätzlich von Hand Zeilen, die die Invariante verletzen, und prüft, dass sie
abgefangen werden — die Invariante ist also nicht leer.

### Was dabei herauskommt

Bei 13 Größen:

| Achse | erreicht |
|---|---|
| vergleichbar (ganz oder nach Reduktion) | **9** |
| implementiert | **8** |
| numerisch konvergiert | **0** |
| mit Messdaten vergleichbar | **6** |
| **tatsächlich validiert** | **0** |
| blockiert | **5** |

**Der Abstand zwischen der ersten und der vorletzten Zahl ist der ganze Punkt
dieser Tabelle.** Dass nichts validiert ist, ist keine Bescheidenheit, sondern
eine geprüfte Aussage: es sind überhaupt keine Messdaten importiert, und der
Test schlägt fehl, sobald eine Zeile etwas anderes behauptet.

Die **vier grundsätzlich nicht vergleichbaren** Größen — azimutale Asymmetrie
des Strahls, Versatz von Emitter und Blende, Neigung des Emitters,
Emitter-zu-Emitter-Übersprechen im Array — sind genau das, was ein
achsensymmetrisches Modell nicht hat. Sie stehen benannt in der Matrix, statt zu
fehlen, und sie sind die Antwort auf die Frage, wofür ein 3D-Löser gebraucht
würde.

Jede Zeile trägt außerdem **einen Satz zur Konvergenz und einen zur
Validierung**, statt einer Farbe: warum eine Größe nicht konvergiert ist (bei
der integrierten Maxwell-Kraft etwa das von P0 gemessene und verfehlte
1-%-Ziel, 4,6 bis 6,1 %) und warum sie nicht validiert ist.

**Eine Zeile hat sich durch die P3-Korrektur geändert:** die
Ladungsrelaxationszeit stand als `MissingMaterialData` und steht jetzt als
„Band belegt, Einzelwert fehlt" — implementiert, geschlossen lösbar, mit
Messdaten vergleichbar, und weiterhin **nicht validiert**.

## 20.5 Kunze-Geometrien und Messdaten

Die Gerätegeometrie dieses Projekts (P1) modelliert die geraden 10-µm-Kapillaren
aus SU-8 und IP-Q, für die die Kunze-Dissertation EMI-BF4 als Treibstoff
ausweist (Fundstelle in `include/es/liquid.hpp`). **Messdaten aus dieser Arbeit
sind hier nicht importiert**: der Importvertrag verlangt zu jedem Punkt Einheit,
Unsicherheit mit Typ, Fundstelle und Geometrieart, und diese Zuordnung ist eine
Arbeit am Dokument, die dieser Lauf nicht geleistet hat.

Das Schema steht bereit; die Daten fehlen, und das ist der ehrliche Zustand.

## 20.6 Was P9 ausdrücklich nicht enthält

Kein 3D-Netz, kein 3D-Löser, kein 3D-Ergebnis, keine 3D-Vernetzung, keine neue
Abhängigkeit. Keine importierten Messdaten. Die Rotationsreferenz ist eine
Prüfung der Wichtung und keine Rechnung.
