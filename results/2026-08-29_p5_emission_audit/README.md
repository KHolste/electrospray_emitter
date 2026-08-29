# P5 — Ionenemission — 2026-08-29

Status: **`blocked`**. Das Modell ist implementiert als Vertrag, standardmaessig
**abgeschaltet**, und liefert auf keinem Pfad eine Zahl.

```sh
./build/es_emission_audit results/2026-08-29_p5_emission_audit meta.commit=$(git rev-parse HEAD)
python python/plot_emission_audit.py results/2026-08-29_p5_emission_audit
```

Literaturaudit und Vertrag:
[`docs/16_emission_contract.md`](../../docs/16_emission_contract.md).

## Der Blocker, zweiteilig

**1. Die Gleichung ist nicht geprueft.** Die Ratengleichung von Iribarne und
Thomson wird durchgaengig zitiert als j = (kT/h) sigma_s exp[-(dG - G(E))/kT]
mit sigma_s = eps0 E, und genau das rechnet `src/emission.cpp` bereits. Die
Primaerquellen -- Iribarne & Thomson (1976) J. Chem. Phys. 64, 2287-2294 und
Thomson & Iribarne (1979) J. Chem. Phys. 71, 4451-4463 -- waren in diesem Lauf
**nicht erreichbar**. Vorfaktor, genaue Definition von sigma_s und
Gueltigkeitsbereich sind damit unbelegt.

Gelesen wurde eine Sekundaerquelle (Wiley, Leseprobe zu "A Brief Overview of the
Mechanisms Involved in Electrospray Mass Spectrometry"). Sie zitiert beide
Arbeiten vollstaendig, **druckt aber keine Gleichung** und stellt fest, das
Modell sei "experimentally well supported for small ions" -- mit ausdruecklichem
Vorbehalt fuer groessere Ionen, also genau fuer den Fall ionischer
Fluessigkeiten.

**2. Die Barriere ist nicht belegt.** Fuer EMI-BF4 wurde keine Quelle mit einer
Zahl gefunden. Der Wert 1,09 eV in `src/fluid.cpp` hat keine Quelle; sein
eigener Header nennt ihn "the least certain quantity here by far" mit einer
Spanne von 1,0 bis 1,4 eV.

**Als Zahl:** ueber diese Spanne aendert sich die Rate bei festem Feld um den
Faktor **5,8e6**. Ein Modell mit unbelegter Barriere sagt keinen Strom voraus;
es berichtet den Parameter.

## Was implementiert ist

Ein Vertrag, in dem jede physikalische Eingabe explizit ist -- Spezies,
LADUNGSVORZEICHEN, Masse, Aktivierungsbarriere, Temperatur, Feldrichtung -- und
jede einen eigenen fail-closed-Status hat. Sechs Pfade koennten eine Zahl
liefern; `contract.csv` zeigt, welchen Status jeder stattdessen meldet:

| Pfad | Status |
|---|---|
| ausgeliefert, positiv/negativ | `Disabled` |
| eingeschaltet, positiv/negativ | `MissingEmissionParameters` |
| Gleichung freigegeben, Barriere fehlt | `MissingEmissionParameters` |

**Positive und negative Polaritaet werden nie aus denselben Speziesdaten
gerechnet.** Ein Kation emittiert nur bei nach aussen zeigendem Feld, ein Anion
nur bei umgekehrtem. Der Prototyp benutzte |E_n| und Kationenmassen fuer beide
Vorzeichen.

**Kein Rueckgriff auf die alten quellenlosen Diagnosewerte:** dieses Modul liest
`es::Fluid` nicht und hat fuer keinen Parameter einen Vorgabewert.

## Geprueft, ohne etwas vorherzusagen

Kern und Modell sind getrennt. Die mathematische Form ist eine reine Funktion
ihrer Argumente und wird ohne jede Stoffangabe geprueft: Dimensionen (die Rate
bei E* ist genau der Vorfaktor (kT/h) eps0 E*), Grenzwerte (E = 0 gibt exakt
null; T = 0 und dG = 0 geben nan), Monotonie in E ueber zehn Dekaden, Monotonie
in dG, das sqrt(E)-Gesetz der Barrierensenkung, und dass die dimensionslose Form
dieselbe Funktion ist.

**Die Form zu pruefen ist nicht, das Modell zu validieren.**

## Ein im Test gefundener Fehler

Die erste Fassung von `synthetic_complete_model()` gab ein Modell zurueck, das
auf einen gemeinsamen statischen Puffer zeigte; das Anlegen eines zweiten
Modells ueberschrieb die Spezies des ersten. Gefunden hat es der Test, weil er
beide Polaritaetsstroeme ausdruckt und einer null war. Behoben durch Besitz.

## Abbildung

| Datei | Inhalt |
|---|---|
| `fig1_blocker.png` | Statuskarte jedes Pfades, die dimensionslose Ratenform, und der Blocker als Zahl: die Hebelwirkung der Barriere ueber die Literaturspanne |

Alle Achsen sind dimensionslos oder Verhaeltnisse. **Ein absoluter Strom kommt
nirgends vor.**
