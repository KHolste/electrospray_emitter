# Literatur

**Prüfstatus.** Einträge mit **✓** wurden in dieser Arbeitssitzung gegen die
Verlagsseite bzw. eine Preprint-Quelle geprüft; Autoren, Titel, Jahr, Band und
Artikelnummer sind bestätigt. Alle übrigen Einträge sind aus dem Gedächtnis
wiedergegeben. Sie sind inhaltlich zutreffend zugeordnet, aber
**bibliographische Angaben (Band, Seite, teils Jahr) sind vor einer Zitation in
einer Publikation zu prüfen.** Ich kenne die Arbeiten, kann die genauen
Seitenzahlen aber nicht garantieren.

---

## Geprüft

* **✓ Gallud, X. & Lozano, P. C. (2022).** *The emission properties, structure
  and stability of ionic liquid menisci undergoing electrically assisted ion
  evaporation.* Journal of Fluid Mechanics **933**, A43.
  DOI: 10.1017/jfm.2021.988. Preprint: arXiv:2109.12274.
  → Stationäres EHD-Modell mit feldverstärkter thermionischer Emission,
  axialsymmetrisch. Liefert das in der Spezifikation verwendete
  Stabilitätskriterium und die Existenzbedingungen (minimale hydraulische
  Impedanz, maximaler Strom, Optimum bei Meniskusradien 0,5–3 µm).

* **✓ Higuera, F. J. (2008).** *Model of the meniscus of an ionic-liquid ion
  source.* Physical Review E **77**, 026308.
  DOI: 10.1103/PhysRevE.77.026308.
  → Zeigt für den rein ionischen Betrieb: viskositätsdominierte Strömung im
  Meniskus, geringer Einfluss des Massenflusses durch Verdampfung,
  vernachlässigbare Raumladung an der verdampfenden Oberfläche, und den
  **von der endlichen Leitfähigkeit kontrollierten** Verdampfungsstrom. Zentrale
  Begründung dafür, das Perfect-Conductor-Modell im emittierenden Betrieb
  aufzugeben.

* **✓ Coffman, C., Martínez-Sánchez, M., Higuera, F. J. & Lozano, P. C. (2016).**
  *Structure of the menisci of leaky dielectric liquids during
  electrically-assisted evaporation of ions.* Applied Physics Letters **109**,
  231602. DOI: 10.1063/1.4971778.
  → Familie von Gleichgewichtsformen, die messbare Ladung abgeben, wenn der
  Meniskus groß gegen eine charakteristische Emissionslänge ist.

* **✓ Wohlhuter, F. K. & Basaran, O. A. (1992).** *Shapes and stability of
  pendant and sessile dielectric drops in an electric field.* Journal of Fluid
  Mechanics **235**, 481.
  → Formen **und** Stabilität elektrifizierter Tropfen; konische Spitzen bei
  wachsender Deformation. Referenz für die Stabilitätsanalyse in P3.

---

## Ungeprüft — vor Zitation verifizieren

### Taylor-Kegel, Elektrohydrodynamik

* Zeleny, J. (1917). *Instability of electrified liquid surfaces.* Physical
  Review **10**, 1–6.
* Taylor, G. I. (1964). *Disintegration of water drops in an electric field.*
  Proceedings of the Royal Society A **280**, 383–397.
  → Gleichgewichtskegel, Halbwinkel 49,29° als Nullstelle von
  $P_{1/2}(\cos\theta)$.
* Melcher, J. R. & Taylor, G. I. (1969). *Electrohydrodynamics: a review of the
  role of interfacial shear stresses.* Annual Review of Fluid Mechanics **1**,
  111–146. → Leaky-Dielectric-Modell.
* Saville, D. A. (1997). *Electrohydrodynamics: the Taylor–Melcher leaky
  dielectric model.* Annual Review of Fluid Mechanics **29**, 27–64.
* Fernández de la Mora, J. (2007). *The fluid dynamics of Taylor cones.* Annual
  Review of Fluid Mechanics **39**, 217–243. → Übersicht.
* Pantano, C., Gañán-Calvo, A. M. & Barrero, A. (1994). *Zeroth-order
  electrohydrostatic solution for electrospraying in cone-jet mode.* Journal of
  Aerosol Science **25**, 1065–1077.

### Cone-Jet, Skalengesetze

* Fernández de la Mora, J. & Loscertales, I. G. (1994). *The current emitted by
  highly conducting Taylor cones.* Journal of Fluid Mechanics **260**, 155–184.
  → $I=f(\varepsilon_r)\sqrt{\gamma K Q/\varepsilon_r}$; $f\approx18$ für
  $\varepsilon_r\gtrsim40$.
* Gañán-Calvo, A. M. (1997). *Cone-jet analytical extension of Taylor's
  electrostatic solution and the asymptotic universal scaling laws in
  electrospraying.* Physical Review Letters **79**, 217.
* Gañán-Calvo, A. M., López-Herrera, J. M., Herrada, M. A., Ramos, A. &
  Montanero, J. M. (2018). *Review on the physics of electrospray.* Journal of
  Aerosol Science **125**, 32–56.
* Higuera, F. J. (2003). *Flow rate and electric current emitted by a Taylor
  cone.* Journal of Fluid Mechanics **484**, 303–327.
* Herrada, M. A., López-Herrera, J. M., Gañán-Calvo, A. M., Vega, E. J.,
  Montanero, J. M. & Popinet, S. (2012). *Numerical simulation of electrospray
  in the cone-jet mode.* Physical Review E **86**, 026305.
* Collins, R. T., Jones, J. J., Harris, M. T. & Basaran, O. A. (2008).
  *Electrohydrodynamic tip streaming and emission of charged drops from liquid
  cones.* Nature Physics **4**, 149–154.
* Smith, D. P. H. (1986). *The electrohydrodynamic atomization of liquids.* IEEE
  Transactions on Industry Applications **IA-22**, 527–535. → Onset-Formel für
  Kapillaren; im Code nur als Vergleichswert.
* Rayleigh, Lord (1882). *On the equilibrium of liquid conducting masses charged
  with electricity.* Philosophical Magazine **14**, 184–186.

### Ionenverdampfung, rein ionischer Betrieb

* Iribarne, J. V. & Thomson, B. A. (1976). *On the evaporation of small ions
  from charged droplets.* Journal of Chemical Physics **64**, 2287–2294.
  → Ursprung der im Code verwendeten Ratengleichung.
* Loscertales, I. G. & Fernández de la Mora, J. (1995). *Experiments on the
  kinetics of field evaporation of small ions from droplets.* Journal of
  Chemical Physics **103**, 5041.
* Romero-Sanz, I., Bocanegra, R., Fernández de la Mora, J. &
  Gamero-Castaño, M. (2003). *Source of heavy molecular ions based on Taylor
  cones of ionic liquids operating in the pure ion evaporation regime.* Journal
  of Applied Physics **94**, 3599.
* Lozano, P. & Martínez-Sánchez, M. (2005). *Ionic liquid ion sources:
  characterization of externally wetted emitters.* Journal of Colloid and
  Interface Science **282**, 415–421.
* Coffman, C. S. (2016). *Electrically-Assisted Evaporation of Charged Fluids.*
  Dissertation, MIT.

### Kolloidtriebwerke, Strahl

* Gamero-Castaño, M. & Hruby, V. (2001). *Electrospray as a source of
  nanoparticles for efficient colloid thrusters.* Journal of Propulsion and
  Power **17**, 977–987.
* Gamero-Castaño, M. (2008). *The structure of electrospray beams in vacuum.*
  Journal of Fluid Mechanics **604**, 339–368.
* Legge, R. S. & Lozano, P. C. (2011). *Electrospray propulsion based on
  emitters microfabricated in porous metals.* Journal of Propulsion and Power
  **27**, 485–495.
* Ziemer, J. K. et al., Veröffentlichungen zu ST7-DRS / LISA Pathfinder
  (Busek-Kolloidtriebwerke). → Genaue Stelle noch zu bestimmen.

### Numerik

* Brebbia, C. A., Telles, J. C. F. & Wrobel, L. C. *Boundary Element
  Techniques.* Springer. → Randintegralmethoden allgemein.
* Telles, J. C. F. (1987). *A self-adaptive co-ordinate transformation for
  efficient numerical evaluation of general boundary element integrals.*
  International Journal for Numerical Methods in Engineering **24**, 959–973.
  → Alternative zur hier verwendeten geometrischen Panelverfeinerung.
* Birdsall, C. K. & Langdon, A. B. *Plasma Physics via Computer Simulation.*
  → PIC-Grundlagen, Formfunktionen, Child-Langmuir als Testfall.
* Jackson, J. D. *Classical Electrodynamics.* → Maxwell-Spannungstensor,
  Randbedingungen an Leitern.
* Bashforth, F. & Adams, J. C. (1883). *An Attempt to Test the Theories of
  Capillary Action.* → Tafeln für hängende Tropfen, Referenz A10.

### Stoffdaten ionischer Flüssigkeiten

Für den in P0/P1 vorgesehenen Umbau von `src/fluid.cpp` auf Werte mit Quelle und
Unsicherheit sind heranzuziehen:

* NIST **ILThermo** (Ionic Liquids Database) als primäre Sammelquelle.
* Tokuda, H. et al., Arbeiten zu physikochemischen Eigenschaften
  imidazoliumbasierter ionischer Flüssigkeiten, Journal of Physical Chemistry B.
* Jacquemin, J. et al., Dichte und Viskosität.
* Freire, M. G. et al., Oberflächenspannung.

Die derzeit im Code stehenden Werte sind literaturtypische Raumtemperaturwerte
**ohne Einzelnachweis** und in dieser Form nicht zitierfähig.
