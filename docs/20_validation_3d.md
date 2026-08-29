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

Ein Messwert ohne Einheit ist eine Zahl; ohne Unsicherheit eine Anekdote; ohne
Fundstelle ein Gerücht; und ohne die Geometrie, zu der er gehört, mit nichts
vergleichbar. Alle vier sind Pflicht:

| Feld | fehlt → |
|---|---|
| `unit` | `MissingUnit` |
| `uncertainty` **und** `uncertainty_type` (GUM Typ A oder B) | `MissingUncertainty` |
| `provenance` | `MissingProvenance` |
| `geometry_stated` | `MissingGeometryKind` |
| zwei Einheiten für dieselbe Größe | `UnitMismatch` |

**Ein Satz mit einem unvollständigen Punkt wird als GANZES abgelehnt.** Einen
Punkt stillschweigend wegzulassen wäre ein Vergleich mit einem anderen
Datensatz.

Zusätzlich aufgenommen: `coverage_factor` $k$ und `conditions` (Temperatur,
Polarität, Volumenstrom) — ohne die Bedingungen ist auch ein vollständiger
Messwert nicht zuordenbar.

## 20.4 Die Validierungsmatrix

13 Größen, jede mit ihrer Vergleichbarkeit, der Bedingung dafür, der Phase, die
sie hier rechnet, und deren Status. Die Zusammenfassung:

* **6 direkt vergleichbar** — Spannung, Gesamtstrom, integrierte Maxwell-Kraft,
  Apexhöhe, Apexfeld, Ladungsrelaxationszeit;
* **3 erst nach ausgesprochener Reduktion** — Transmission, Divergenz,
  Auftreffverteilung: achsensymmetrisch sind sie Funktionen von $r$ allein, in
  3D Verteilungen über den Azimut;
* **4 grundsätzlich nicht vergleichbar** — azimutale Asymmetrie des Strahls,
  Versatz von Emitter und Blende, Neigung des Emitters,
  Emitter-zu-Emitter-Übersprechen im Array.

Die letzten vier sind **genau das, was ein achsensymmetrisches Modell nicht
hat**. Sie stehen benannt in der Matrix, statt zu fehlen, und sie sind die
Antwort auf die Frage, wofür ein 3D-Löser gebraucht würde.

Von den 13 rechnet dieses Projekt 8 — und von diesen 8 trägt keine einzige den
Status „validiert": sie sind `geprueft`, `validated_subset`, `qualitativ`,
`DiscretizationNotConverged`, `MissingMaterialData` oder `blocked`.

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
