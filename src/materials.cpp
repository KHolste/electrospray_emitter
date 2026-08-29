#include "es/materials.hpp"

#include <cstdio>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {

const char* to_string(MaterialStatus s) {
  switch (s) {
    case MaterialStatus::Exact: return "exact";
    case MaterialStatus::IdealConductor: return "ideal_conductor";
    case MaterialStatus::Measured: return "measured";
    case MaterialStatus::ManufacturerSpec: return "manufacturer_spec";
    case MaterialStatus::Literature: return "literature";
    case MaterialStatus::Provisional: return "provisional";
    case MaterialStatus::Unknown: return "unknown";
  }
  return "unknown";
}

Real Material::permittivity_or_throw() const {
  if (ideal_conductor)
    throw std::runtime_error("Material '" + name +
                             "' ist ein idealer Leiter; eine Permittivitaet ist dort nicht "
                             "definiert und darf nicht abgefragt werden.");
  if (status == MaterialStatus::Unknown || relative_permittivity <= 0.0)
    throw std::runtime_error(
        "Material '" + name +
        "' ist registriert, traegt aber keinen Permittivitaetswert.\n"
        "  Es wurde bewusst KEINE Zahl eingetragen: fuer dieses Harz liegt keine "
        "belastbare\n"
        "  elektrostatische Materialkonstante vor, und eine erfundene waere schlimmer "
        "als ein Abbruch.\n"
        "  Abhilfe: material.<name>.relative_permittivity und "
        "material.<name>.source in der\n"
        "  Konfiguration angeben. Ein Codeeingriff ist dafuer nicht noetig.");
  return relative_permittivity;
}

void Material::print(std::FILE* out) const {
  if (ideal_conductor) {
    std::fprintf(out, "  %-20s idealer Leiter (Aequipotential), eps_r geht nicht ein\n",
                 name.c_str());
  } else if (status == MaterialStatus::Unknown) {
    std::fprintf(out, "  %-20s KEIN WERT (%s)\n", name.c_str(), to_string(status));
  } else {
    std::fprintf(out, "  %-20s eps_r = %.4g   [%s]\n", name.c_str(), relative_permittivity,
                 to_string(status));
    if (has_range())
      std::fprintf(out, "  %-20s   Sensitivitaetsbereich %.3g ... %.3g\n", "", eps_r_low,
                   eps_r_high);
  }
  if (!source.empty()) std::fprintf(out, "  %-20s   Quelle    : %s\n", "", source.c_str());
  if (!conditions.empty()) std::fprintf(out, "  %-20s   Bedingung : %s\n", "", conditions.c_str());
  if (!caveat.empty()) std::fprintf(out, "  %-20s   Vorbehalt : %s\n", "", caveat.c_str());
}

// ---------------------------------------------------------------------------

MaterialLibrary::MaterialLibrary() {
  {
    Material m;
    m.name = "vacuum";
    m.relative_permittivity = 1.0;
    m.status = MaterialStatus::Exact;
    m.source = "Definition";
    m.conditions = "-";
    m_.push_back(m);
  }
  {
    // ------------------------------------------------------------------
    // SU-8 -- the first dielectric, and the one the small Kunze geometry
    // generation (capillaries of 8-10 um) was printed in.
    //
    // WHAT IS AVAILABLE
    //   * MicroChem / Kayaku "SU-8 3000 Permanent Epoxy Negative Photoresist"
    //     data sheet, properties table: "Dielectric Constant @ 1GHz" = 3.28,
    //     bulk resistivity 7.8e14 Ohm cm.
    //   * MicroChem / Kayaku "SU-8 2000" data sheet: dielectric constant 4.1
    //     at 1 GHz, 50 % relative humidity.
    //   * N. Ghalichechian and K. Sertel, "Permittivity and Loss
    //     Characterization of SU-8 Films for mmW and Terahertz Applications",
    //     IEEE Antennas Wirel. Propag. Lett. 14 (2015) 723-726,
    //     doi:10.1109/LAWP.2014.2377695: fully cross-linked 430 um film,
    //     eps_r = 3.24 (1 GHz), 3.23 (200 GHz), 2.92 (1 THz).
    //   * J. Melai, C. Salm, S. Smits, J. Visschers, J. Schmitz, "The electrical
    //     conduction and dielectric strength of SU-8", J. Micromech. Microeng.
    //     19 (2009) 065012, doi:10.1088/0960-1317/19/6/065012: dielectric
    //     strength 4.4 MV/cm, leakage dominated by thermionic emission -- i.e.
    //     SU-8 is an insulator at the fields of interest, which is the premise
    //     of the whole dielectric model.
    //
    // WHAT IS NOT AVAILABLE
    //   No source gives eps_r in the STATIC limit for a two-photon-printed,
    //   hard-baked part in vacuum.  Every number above is at 1 GHz or higher.
    //   For a passive dielectric eps'(omega) is non-increasing in omega, so the
    //   static value is at least the 1 GHz value; the orientational
    //   contribution that the GHz measurement partly misses can only raise it.
    //
    // THE CHOICE
    //   Nominal 3.3 -- the 1 GHz manufacturer/literature cluster (3.24 ... 3.3),
    //   rounded up marginally to acknowledge that the static value cannot be
    //   lower.  Marked PROVISIONAL.  Sensitivity range 2.8 ... 4.5: the low end
    //   is the lowest reported measurement on SU-8 (microstrip characterisation,
    //   about 2.85), the high end sits above the highest manufacturer figure
    //   (4.1 at 1 GHz and 50 % RH) to cover both formulation spread and the
    //   low-frequency rise.  Water uptake raises eps_r, and the device runs in
    //   vacuum, so the dry end of the range is the more relevant one -- but the
    //   study reports the whole range, because "more relevant" is not a
    //   measurement either.
    // ------------------------------------------------------------------
    Material m;
    m.name = "su8";
    m.relative_permittivity = 3.3;
    m.status = MaterialStatus::Provisional;
    m.eps_r_low = 2.8;
    m.eps_r_high = 4.5;
    m.source =
        "MicroChem/Kayaku SU-8 3000 Datenblatt (eps_r = 3.28 @ 1 GHz); SU-8 2000 Datenblatt "
        "(4.1 @ 1 GHz, 50 % r.F.); Ghalichechian & Sertel, IEEE AWPL 14 (2015) 723, "
        "doi:10.1109/LAWP.2014.2377695 (3.24 @ 1 GHz an vollstaendig vernetztem 430-um-Film)";
    m.conditions =
        "alle Quellen bei 1 GHz oder darueber; Feuchte teils 50 % r.F.; Vernetzungsgrad "
        "der Quellen nicht der unserer zweiphotonengedruckten, hartgebackenen Teile";
    m.caveat =
        "VORLAEUFIG. Kein statischer Messwert fuer unseren Aushaertungszustand bekannt. "
        "Es darf keine Validierungsaussage aus einer Rechnung mit diesem Wert abgeleitet "
        "werden; berichtet wird die Sensitivitaet ueber 2.8 ... 4.5.";
    m_.push_back(m);
  }
  {
    Material m;
    m.name = "ip-q";
    m.status = MaterialStatus::Unknown;
    m.source = "keine";
    m.conditions = "-";
    m.caveat =
        "In der Kunze-Dissertation als verwendetes Harz genannt, ohne dielektrische Daten. "
        "Struktur vorgesehen, Wert absichtlich leer.";
    m_.push_back(m);
  }
  {
    Material m;
    m.name = "ipx-q";
    m.status = MaterialStatus::Unknown;
    m.source = "keine";
    m.conditions = "-";
    m.caveat =
        "In der Kunze-Dissertation als verwendetes Harz genannt, ohne dielektrische Daten. "
        "Struktur vorgesehen, Wert absichtlich leer.";
    m_.push_back(m);
  }
  {
    Material m;
    m.name = "emibf4";
    m.status = MaterialStatus::IdealConductor;
    m.ideal_conductor = true;
    m.source = "Modellannahme P2b";
    m.conditions = "statisch";
    m.caveat = liquid_model::why_ideal_conductor_is_admissible_in_p2b();
    m_.push_back(m);
  }
  {
    Material m;
    m.name = "metal";
    m.status = MaterialStatus::IdealConductor;
    m.ideal_conductor = true;
    m.source = "Modellannahme P2b";
    m.conditions = "statisch";
    m.caveat =
        "Metallisierung ohne Dicke: eine ideal leitende Flaeche auf dem Polymertraeger. "
        "Die reale Schichtdicke (typisch 10^-7 m) ist gegen jede Geraeteabmessung "
        "vernachlaessigt; das ist zulaessig, solange keine Groesse innerhalb der Schicht "
        "gefragt wird.";
    m_.push_back(m);
  }
}

Material* MaterialLibrary::find(const std::string& name) {
  for (Material& m : m_)
    if (m.name == name) return &m;
  return nullptr;
}

bool MaterialLibrary::has(const std::string& name) const {
  for (const Material& m : m_)
    if (m.name == name) return true;
  return false;
}

const Material& MaterialLibrary::get(const std::string& name) const {
  for (const Material& m : m_)
    if (m.name == name) return m;
  std::string known;
  for (const Material& m : m_) known += (known.empty() ? "" : ", ") + m.name;
  throw std::runtime_error("Unbekanntes Material '" + name + "'. Bekannt: " + known +
                           ". Neue Materialien werden ueber material.<name>.* in der "
                           "Konfiguration eingetragen.");
}

void MaterialLibrary::set(const Material& m) {
  if (Material* p = find(m.name))
    *p = m;
  else
    m_.push_back(m);
}

void MaterialLibrary::override_permittivity(const std::string& name, Real eps_r,
                                            MaterialStatus status, const std::string& source) {
  Material* p = find(name);
  if (!p) {
    Material m;
    m.name = name;
    m_.push_back(m);
    p = find(name);
  }
  if (p->ideal_conductor)
    throw std::runtime_error("Material '" + name +
                             "' ist als idealer Leiter gefuehrt; eine Permittivitaet dafuer "
                             "zu setzen waere widerspruechlich.");
  if (!(eps_r >= 1.0))
    throw std::runtime_error("relative_permittivity fuer '" + name +
                             "' muss >= 1 sein (angegeben: " + std::to_string(eps_r) + ").");
  p->relative_permittivity = eps_r;
  p->status = status;
  p->source = source;
  // A hand-supplied value carries no sensitivity range unless one is supplied
  // with it; leaving the old range would attribute it to the new number.
  p->eps_r_low = p->eps_r_high = 0.0;
}

void MaterialLibrary::print(std::FILE* out) const {
  std::fprintf(out, "Materialbibliothek (statische relative Permittivitaeten)\n");
  for (const Material& m : m_) m.print(out);
}

void MaterialLibrary::write_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# registrierte Materialien; status=unknown bedeutet absichtlich kein Wert\n");
  std::fprintf(f, "name,relative_permittivity,status,ideal_conductor,eps_r_low,eps_r_high,"
                  "source,conditions,caveat\n");
  for (const Material& m : m_)
    std::fprintf(f, "%s,%.9g,%s,%d,%.9g,%.9g,\"%s\",\"%s\",\"%s\"\n", m.name.c_str(),
                 m.relative_permittivity, to_string(m.status), m.ideal_conductor ? 1 : 0,
                 m.eps_r_low, m.eps_r_high, m.source.c_str(), m.conditions.c_str(),
                 m.caveat.c_str());
  std::fclose(f);
}

// ---------------------------------------------------------------------------

DielectricMaterials DielectricMaterials::reference(const MaterialLibrary& lib) {
  DielectricMaterials d;
  d.vacuum = lib.get("vacuum");
  d.emitter_dielectric = lib.get("su8");
  d.liquid = lib.get("emibf4");
  d.extractor_carrier = lib.get("su8");
  d.metallisation = lib.get("metal");
  return d;
}

DielectricMaterials DielectricMaterials::all_vacuum(const MaterialLibrary& lib) {
  DielectricMaterials d = reference(lib);
  const Material v = lib.get("vacuum");
  d.emitter_dielectric = v;
  d.emitter_dielectric.name = "emitter_dielectric_as_vacuum";
  d.emitter_dielectric.caveat =
      "eps_r = 1 erzwungen. KEIN physikalischer Zustand des Emitters, sondern der Sonderfall, "
      "in dem die unabhaengige BEM dasselbe Problem loesen kann.";
  d.extractor_carrier = v;
  d.extractor_carrier.name = "extractor_carrier_as_vacuum";
  d.extractor_carrier.caveat = d.emitter_dielectric.caveat;
  return d;
}

void DielectricMaterials::check_usable() const {
  (void)emitter_dielectric.permittivity_or_throw();
  (void)extractor_carrier.permittivity_or_throw();
  (void)vacuum.permittivity_or_throw();
  if (!liquid.ideal_conductor)
    throw std::runtime_error("P2b behandelt die ionische Fluessigkeit als idealen Leiter; "
                             "ein Dielektrikum an dieser Stelle ist nicht implementiert.");
  if (!metallisation.ideal_conductor)
    throw std::runtime_error("Die Extraktormetallisierung muss ein idealer Leiter sein.");
}

void DielectricMaterials::print(std::FILE* out) const {
  std::fprintf(out, "Materialbelegung der P2b-Gebiete\n");
  std::fprintf(out, "  Vakuum              -> ");
  vacuum.print(out);
  std::fprintf(out, "  Emitterkoerper      -> ");
  emitter_dielectric.print(out);
  std::fprintf(out, "  ionische Fluessigk. -> ");
  liquid.print(out);
  std::fprintf(out, "  Extraktortraeger    -> ");
  extractor_carrier.print(out);
  std::fprintf(out, "  Metallisierung      -> ");
  metallisation.print(out);
}

void DielectricMaterials::write_csv(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::fprintf(f, "# Materialbelegung der Gebiete des P2b-Modells\n");
  std::fprintf(f, "role,material,relative_permittivity,status,ideal_conductor,"
                  "eps_r_low,eps_r_high,source,conditions,caveat\n");
  auto row = [&](const char* role, const Material& m) {
    std::fprintf(f, "%s,%s,%.9g,%s,%d,%.9g,%.9g,\"%s\",\"%s\",\"%s\"\n", role, m.name.c_str(),
                 m.relative_permittivity, to_string(m.status), m.ideal_conductor ? 1 : 0,
                 m.eps_r_low, m.eps_r_high, m.source.c_str(), m.conditions.c_str(),
                 m.caveat.c_str());
  };
  row("vacuum", vacuum);
  row("emitter_dielectric", emitter_dielectric);
  row("liquid", liquid);
  row("extractor_carrier", extractor_carrier);
  row("metallisation", metallisation);
  std::fclose(f);
}

// ---------------------------------------------------------------------------

namespace liquid_model {

Real charge_relaxation_time(Real sigma, Real eps_r) {
  if (!(sigma > 0.0)) throw std::runtime_error("Leitfaehigkeit muss positiv sein");
  return constants::eps0 * eps_r / sigma;
}

const char* why_ideal_conductor_is_admissible_in_p2b() {
  return "Statisch zulaessig: die Ladungsrelaxationszeit eps/sigma der ionischen "
         "Fluessigkeit liegt bei etwa 1e-10 s und ist gegen jeden in P2b beschriebenen "
         "Vorgang vernachlaessigbar. NICHT mehr zulaessig, sobald Strom durch den Meniskus "
         "fliesst: dann bestimmen ohmscher Abfall und endliche Nachlieferung von "
         "Oberflaechenladung das Apexfeld (Leaky-Dielectric-Regime). P3/P4 muessen die "
         "Dirichlet-Bedingung auf der Fluessigkeit durch ein Leitungsmodell ersetzen.";
}

}  // namespace liquid_model

}  // namespace es
