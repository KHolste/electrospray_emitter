# Literatur

## Prüfstatus und dessen Grenzen

Alle Einträge wurden gegen Verlagsseiten bzw. Datenbankeinträge abgeglichen.
Bitte den Unterschied beachten:

| Marke | Bedeutung |
|---|---|
| **B** | **Bibliographisch geprüft**: Autoren, Titel, Jahr, Band, Seiten/Artikelnummer und DOI gegen die Verlagsseite abgeglichen. |
| **I** | Zusätzlich **inhaltlich** insoweit bestätigt, als die vom Verlag veröffentlichte Zusammenfassung die hier zugeschriebene Aussage trägt. |
| **—** | Standardwerk oder Qualifikationsarbeit; nicht abgeglichen, wird auch nicht als Grundlage einer Modellentscheidung geführt. |

**Was hier nicht geleistet ist:** die Originalarbeiten wurden nicht im Volltext
gelesen. Ein **I** bedeutet, dass Titel und Zusammenfassung die zugeschriebene
Aussage stützen — nicht, dass Gleichungen, Gültigkeitsbereiche oder Zahlenwerte
aus dem Volltext nachgeprüft sind. Vor der Übernahme einer Gleichung in den Code
ist die betreffende Stelle im Volltext zu lesen; das ist in den jeweiligen
Phasen des [Stufenplans](05_implementation_plan.md) vorgesehen.

**Stand:** 27 von 27 zitierten Fachartikeln bibliographisch geprüft; davon 11
zusätzlich inhaltlich über die Zusammenfassung abgeglichen. Nicht abgeglichen
sind vier Lehrbücher, eine Dissertation und die noch zu erschließenden
Stoffdatenquellen.

---

## Meniskus im emittierenden Betrieb — trägt die Modellstufen 2 und 3

**B I** **Gallud, X. & Lozano, P. C. (2022).** *The emission properties,
structure and stability of ionic liquid menisci undergoing electrically assisted
ion evaporation.* Journal of Fluid Mechanics **933**, A43.
DOI: [10.1017/jfm.2021.988](https://doi.org/10.1017/jfm.2021.988).
Preprint: arXiv:2109.12274.
→ Stationäres EHD-Modell mit feldverstärkter thermionischer Emission,
axialsymmetrisch. Zusammenfassung bestätigt: statisch stabile Lösungen sind an
eine minimale hydraulische Impedanz, einen maximalen Strom und ein enges Fenster
von Hintergrundfeldern gebunden, mit Optimum bei Meniskusradien 0,5–3 µm; die
statische Stabilität geht verloren, wenn der elektrische Druck am Feld
unmittelbar an der haltenden Elektrode den doppelten
Oberflächenspannungsdruck einer gleich großen Kugel übersteigt.

**B I** **Higuera, F. J. (2008).** *Model of the meniscus of an ionic-liquid ion
source.* Physical Review E **77**, 026308.
DOI: [10.1103/PhysRevE.77.026308](https://doi.org/10.1103/PhysRevE.77.026308).
→ Zusammenfassung bestätigt: im rein ionischen Regime ist die Strömung im
Meniskus viskositätsdominiert, der Massenfluss durch Verdampfung beeinflusst sie
kaum, die Raumladung um die verdampfende Oberfläche ist vernachlässigbar, und
der Verdampfungsstrom wird von der **endlichen elektrischen Leitfähigkeit**
kontrolliert. Zentrale Begründung dafür, das Perfect-Conductor-Modell im
emittierenden Betrieb aufzugeben.

**B I** **Coffman, C., Martínez-Sánchez, M., Higuera, F. J. & Lozano, P. C.
(2016).** *Structure of the menisci of leaky dielectric liquids during
electrically-assisted evaporation of ions.* Applied Physics Letters **109**,
231602. DOI: [10.1063/1.4971778](https://doi.org/10.1063/1.4971778).
→ Familie von Gleichgewichtsformen, die messbare Ladung abgeben, wenn der
Meniskus groß gegen eine charakteristische Emissionslänge ist.

**B I** **Wohlhuter, F. K. & Basaran, O. A. (1992).** *Shapes and stability of
pendant and sessile dielectric drops in an electric field.* Journal of Fluid
Mechanics **235**, 481.
→ Formen **und** Stabilität elektrifizierter Tropfen; konische Spitzen bei
wachsender Deformation, unabhängig vom Permittivitätsverhältnis. Referenz für
die Stabilitätsanalyse in Phase P3.

---

## Taylor-Kegel und Elektrohydrodynamik

**B** **Zeleny, J. (1917).** *Instability of electrified liquid surfaces.*
Physical Review **10** (1), 1–6.
DOI: [10.1103/PhysRev.10.1](https://doi.org/10.1103/PhysRev.10.1).

**B I** **Taylor, G. I. (1964).** *Disintegration of water drops in an electric
field.* Proceedings of the Royal Society A **280** (1382), 383–397.
DOI: [10.1098/rspa.1964.0151](https://doi.org/10.1098/rspa.1964.0151).
→ Gleichgewichtskegel, Halbwinkel 49,29° als Nullstelle von
$P_{1/2}(\cos\theta)$. Widerlegt ausdrücklich Zelenys Instabilitätskriterium.

**B** **Melcher, J. R. & Taylor, G. I. (1969).** *Electrohydrodynamics: a review
of the role of interfacial shear stresses.* Annual Review of Fluid Mechanics
**1**, 111–146.
DOI: [10.1146/annurev.fl.01.010169.000551](https://doi.org/10.1146/annurev.fl.01.010169.000551).
→ Leaky-Dielectric-Modell; Grundlage der Sprungbedingungen in Spezifikation 3.2.

**B I** **Saville, D. A. (1997).** *Electrohydrodynamics: the Taylor–Melcher
leaky dielectric model.* Annual Review of Fluid Mechanics **29**, 27–64.
DOI: [10.1146/annurev.fluid.29.1.27](https://doi.org/10.1146/annurev.fluid.29.1.27).
→ Zusammenfassung weist ausdrücklich darauf hin, dass frühe experimentelle
Prüfungen die qualitativen Züge stützten, die **quantitative** Übereinstimmung
aber schlecht war. Beim Vergleich eigener Rechnungen zu berücksichtigen.

**B** **Fernández de la Mora, J. (2007).** *The fluid dynamics of Taylor cones.*
Annual Review of Fluid Mechanics **39**, 217–243.
DOI: [10.1146/annurev.fluid.39.050905.110159](https://doi.org/10.1146/annurev.fluid.39.050905.110159).
→ Übersicht; deckt Leitfähigkeiten > 10⁻⁴ S/m ab.

**B** **Pantano, C., Gañán-Calvo, A. M. & Barrero, A. (1994).**
*Zeroth-order, electrohydrostatic solution for electrospraying in cone-jet
mode.* Journal of Aerosol Science **25** (6), 1065–1077.
DOI: [10.1016/0021-8502(94)90202-X](https://doi.org/10.1016/0021-8502(94)90202-X).

**B I** **Collins, R. T., Jones, J. J., Harris, M. T. & Basaran, O. A. (2008).**
*Electrohydrodynamic tip streaming and emission of charged drops from liquid
cones.* Nature Physics **4**, 149–154.
DOI: [10.1038/nphys807](https://doi.org/10.1038/nphys807).
→ Zusammenfassung bestätigt einen für dieses Projekt wichtigen Punkt: Tip
Streaming tritt **nicht** auf, wenn die Flüssigkeit perfekt leitend *oder*
perfekt isolierend ist. Ein unabhängiges Argument gegen das
Perfect-Conductor-Modell im emittierenden Betrieb.

---

## Cone-Jet und Skalengesetze

**B I** **Fernández de la Mora, J. & Loscertales, I. G. (1994).** *The current
emitted by highly conducting Taylor cones.* Journal of Fluid Mechanics **260**,
155–184.
→ Zusammenfassung bestätigt die im Code verwendete Form
$I=f(\varepsilon)\,(\gamma Q K/\varepsilon)^{1/2}$ und die Bedingung, dass der
Jetdurchmesser klein gegen den Kapillardurchmesser sein muss.

**B** **Gañán-Calvo, A. M. (1997).** *Cone-jet analytical extension of Taylor's
electrostatic solution and the asymptotic universal scaling laws in
electrospraying.* Physical Review Letters **79** (2), 217–220.
DOI: [10.1103/PhysRevLett.79.217](https://doi.org/10.1103/PhysRevLett.79.217).
→ **Achtung: Erratum in Physical Review Letters 85, 4193 (2000).** Beim
Übernehmen von Vorfaktoren zu berücksichtigen.

**B** **Gañán-Calvo, A. M., López-Herrera, J. M., Herrada, M. A., Ramos, A. &
Montanero, J. M. (2018).** *Review on the physics of electrospray: from
electrokinetics to the operating conditions of single and coaxial Taylor
cone-jets, and AC electrospray.* Journal of Aerosol Science **125**, 32–56.

**B I** **Higuera, F. J. (2003).** *Flow rate and electric current emitted by a
Taylor cone.* Journal of Fluid Mechanics **484**, 303–327.
→ Zusammenfassung bestätigt: numerische Lösung des Kegel-Jet-Übergangsbereichs,
lokales Problem mit drei dimensionslosen Parametern (zwei Stoffgrößen, einer für
den Fluss). Referenz B5 der Validierungsmatrix.

**B I** **Herrada, M. A., López-Herrera, J. M., Gañán-Calvo, A. M., Vega, E. J.,
Montanero, J. M. & Popinet, S. (2012).** *Numerical simulation of electrospray
in the cone-jet mode.* Physical Review E **86**, 026305.
DOI: [10.1103/PhysRevE.86.026305](https://doi.org/10.1103/PhysRevE.86.026305).
→ Zusammenfassung nennt als Hauptvereinfachung, dass alle freie Ladung auf der
Grenzfläche sitzt — dieselbe Annahme wie in Spezifikation 3.1.

**B I** **Smith, D. P. H. (1986).** *The electrohydrodynamic atomization of
liquids.* IEEE Transactions on Industry Applications **IA-22** (3), 527–535.
DOI: [10.1109/TIA.1986.4504754](https://doi.org/10.1109/TIA.1986.4504754).
→ Experimentelle Arbeit zu Onset-Potential, Kapillarradius, Leitfähigkeit,
Viskosität; berichtet ausdrücklich **Hysterese** in der I–U-Kennlinie. Die
daraus gebräuchliche Onset-Formel wird im Code nur als Vergleichswert geführt.

**B** **Rayleigh, Lord (1882).** *On the equilibrium of liquid conducting masses
charged with electricity.* Philosophical Magazine **14** (87), 184–186.
DOI: [10.1080/14786448208628425](https://doi.org/10.1080/14786448208628425).

---

## Ionenverdampfung und rein ionischer Betrieb

**B** **Iribarne, J. V. & Thomson, B. A. (1976).** *On the evaporation of small
ions from charged droplets.* Journal of Chemical Physics **64** (6), 2287–2294.
DOI: [10.1063/1.432536](https://doi.org/10.1063/1.432536).
→ Ursprung der im Code verwendeten Ratengleichung.

**B** **Loscertales, I. G. & Fernández de la Mora, J. (1995).** *Experiments on
the kinetics of field evaporation of small ions from droplets.* Journal of
Chemical Physics **103** (12), 5041–5060.

**B I** **Romero-Sanz, I., Bocanegra, R., Fernández de la Mora, J. &
Gamero-Castaño, M. (2003).** *Source of heavy molecular ions based on Taylor
cones of ionic liquids operating in the pure ion evaporation regime.* Journal of
Applied Physics **94** (5), 3599–3605.
DOI: [10.1063/1.1598281](https://doi.org/10.1063/1.1598281).
→ Zusammenfassung bestätigt: EMI-BF4 bei Raumtemperatur, Flugzeit-
Massenspektrometrie des vollständigen Sprays. Referenz C1 der
Validierungsmatrix.

**B I** **Lozano, P. & Martínez-Sánchez, M. (2005).** *Ionic liquid ion sources:
characterization of externally wetted emitters.* Journal of Colloid and
Interface Science **282** (2), 415–421.
DOI: [10.1016/j.jcis.2004.08.132](https://doi.org/10.1016/j.jcis.2004.08.132).
→ Zusammenfassung bestätigt: elektrochemisch geschärfte Wolframspitzen,
EMI-BF4, Hochvakuum. Referenz C2.

**—** Coffman, C. S. (2016). *Electrically-Assisted Evaporation of Charged
Fluids: Fundamental Modeling and Studies on Ionic Liquids.* Dissertation, MIT.
→ Nicht abgeglichen. Wird nicht als alleinige Grundlage einer
Modellentscheidung geführt; die begutachtete Kurzfassung ist Coffman et al.
(2016) oben.

---

## Kolloidtriebwerke und Strahl

**B** **Gamero-Castaño, M. & Hruby, V. (2001).** *Electrospray as a source of
nanoparticles for efficient colloid thrusters.* Journal of Propulsion and Power
**17** (5), 977–987. DOI: [10.2514/2.5858](https://doi.org/10.2514/2.5858).

**B** **Gamero-Castaño, M. (2008).** *The structure of electrospray beams in
vacuum.* Journal of Fluid Mechanics **604**, 339–368.
DOI: [10.1017/S0022112008001316](https://doi.org/10.1017/S0022112008001316).
→ Referenz C3.

**B** **Legge, R. S. & Lozano, P. C. (2011).** *Electrospray propulsion based on
emitters microfabricated in porous metals.* Journal of Propulsion and Power
**27** (2), 485–495. DOI: [10.2514/1.50037](https://doi.org/10.2514/1.50037).
→ Referenz C4. Relevant erst für das spätere Darcy-Modell poröser Emitter.

**B** **Anderson, G. et al. (2018).** *Experimental results from the ST7 mission
on LISA Pathfinder.* arXiv:1809.08969 (Physical Review D).
**B** **Ziemer, J. K. et al. (2008).** *ST7-DRS Colloid Thruster System
Development and Performance Summary.* AIAA 2008-4824.
DOI: [10.2514/6.2008-4824](https://doi.org/10.2514/6.2008-4824).
→ Referenz C6. Flugdaten der Busek-CMNT; 5–30 µN je Triebwerk mit 0,1 µN
Auflösung, rund 2400 Betriebsstunden im Flug. Autorenlisten beider Arbeiten sind
nur teilweise abgeglichen und vor Zitation zu prüfen.

---

## Numerik

**B** **Telles, J. C. F. (1987).** *A self-adaptive co-ordinate transformation
for efficient numerical evaluation of general boundary element integrals.*
International Journal for Numerical Methods in Engineering **24**, 959–973.
→ Alternative zur hier verwendeten geometrischen Panelverfeinerung.

**B** **Bashforth, F. & Adams, J. C. (1883).** *An Attempt to Test the Theories
of Capillary Action.* Cambridge University Press.
→ Tafeln für hängende Tropfen; Referenz A10 der Validierungsmatrix. Volltext im
Internet Archive verfügbar.

**—** Brebbia, C. A., Telles, J. C. F. & Wrobel, L. C. *Boundary Element
Techniques.* Springer. — Lehrbuch.

**—** Birdsall, C. K. & Langdon, A. B. *Plasma Physics via Computer
Simulation.* — Lehrbuch; PIC-Formfunktionen, Child-Langmuir als Testfall.

**—** Jackson, J. D. *Classical Electrodynamics.* — Lehrbuch;
Maxwell-Spannungstensor, Randbedingungen an Leitern.

---

## Dielektrische Eigenschaften von SU-8 — trägt das P2b-Materialmodell

Diese vier Quellen tragen den in `include/es/materials.hpp` eingetragenen Wert
und, wichtiger, seine **Einschränkung**. Keine von ihnen gibt eine statische
Permittivität für einen zweiphotonengedruckten, hartgebackenen Teil im Vakuum;
alle messen bei 1 GHz oder darüber. Deshalb ist der verwendete Wert
ausdrücklich als **vorläufig** geführt und wird von einer Sensitivitätsrechnung
begleitet (siehe [08_dielectric_model.md](08_dielectric_model.md), 8.3).

**B I** **Ghalichechian, N. & Sertel, K. (2015).** *Permittivity and loss
characterization of SU-8 films for mmW and terahertz applications.* IEEE
Antennas and Wireless Propagation Letters **14**, 723–726.
doi:10.1109/LAWP.2014.2377695 — vollständig vernetzter 430-µm-Film,
ε_r = 3,24 / 3,23 / 2,92 bei 1 GHz / 200 GHz / 1 THz.

**B I** **Melai, J., Salm, C., Smits, S., Visschers, J. & Schmitz, J. (2009).**
*The electrical conduction and dielectric strength of SU-8.* Journal of
Micromechanics and Microengineering **19**, 065012.
doi:10.1088/0960-1317/19/6/065012 — Durchschlagfestigkeit 4,4 MV/cm, Leckstrom
thermionisch dominiert. Trägt die Prämisse, dass SU-8 bei den auftretenden
Feldern ein Isolator ist.

**B** MicroChem / Kayaku Advanced Materials, **Datenblatt *SU-8 3000 Permanent
Epoxy Negative Photoresist***, Eigenschaftstabelle: Dielektrizitätszahl 3,28
bei **1 GHz**, Durchgangswiderstand 7,8·10¹⁴ Ω·cm. Herstellerdokument, keine
Messbedingungen jenseits der Frequenz angegeben.

**B** MicroChem / Kayaku Advanced Materials, **Datenblatt *SU-8 2000***:
Dielektrizitätszahl 4,1 bei **1 GHz**, 50 % relative Feuchte. Herstellerdokument;
der Wert wurde aus Sekundärwiedergaben der Eigenschaftstabelle übernommen, das
Datenblatt selbst war nicht direkt abrufbar — deshalb **B** und nicht **B I**.

**Kunze, F. (2024).** Dissertation, Justus-Liebig-Universität Gießen. Nennt
SU-8, IP-Q und IPx-Q als verwendete Harze und die kleine Geometriegeneration
(Kapillaren um 8–10 µm) als in SU-8 gefertigt. **Enthält keine dielektrischen
Daten** — der Volltext führt weder „permittivity" noch „dielectric constant".
Für IP-Q und IPx-Q liegt damit überhaupt kein Wert vor; sie sind im Code
registriert, tragen bewusst keine Zahl und brechen bei Verwendung ab.

---

## Stoffdaten ionischer Flüssigkeiten — noch zu erschließen

Für den in Phase P1 vorgesehenen Umbau von `src/fluid.cpp` auf Werte mit Quelle
und Unsicherheit sind heranzuziehen:

* NIST **ILThermo** (Ionic Liquids Database) als primäre Sammelquelle.
* Arbeiten zu Dichte, Viskosität, Leitfähigkeit und Oberflächenspannung
  imidazoliumbasierter Flüssigkeiten (u. a. Tokuda et al.; Jacquemin et al.;
  Freire et al.).

**Diese Quellen sind noch nicht abgeglichen.** Die derzeit im Code stehenden
Zahlenwerte sind literaturtypische Raumtemperaturwerte ohne Einzelnachweis und
in dieser Form nicht zitierfähig. Bis zum Abschluss von P1 darf keine
Emissionsvorhersage aus diesem Code als stoffdatenbelegt dargestellt werden.
