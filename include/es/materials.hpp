#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "es/status.hpp"
#include "es/types.hpp"

namespace es {

// ===========================================================================
// P2b -- materials of the capillary Kunze emitter
// ===========================================================================
//
// WHAT CHANGED AGAINST P2a, AND WHY
//
// P2a treated the emitter body as metal.  The device it is supposed to model
// has a 3D-printed emitter made of a NON-CONDUCTING photopolymer; high voltage
// sits on the ionic liquid (or its metallic contact), not on the polymer, and
// the extractor is a polymer carrier with a metallised electrode face.  A
// perfect conductor where a dielectric belongs is not a small error: it forces
// an equipotential onto a surface that in reality only polarises, so it changes
// the field at the exit edge by more than any mesh effect.  Every P2a number
// that involved the emitter body is therefore superseded, not merely refined.
//
// WHAT A PERMITTIVITY IN THIS FILE IS AND IS NOT
//
// The electrostatic problem needs eps_r in the STATIC limit, for the material
// in the cure state the part was actually printed and baked in.  Nothing in the
// Kunze dissertation supplies one -- it names SU-8, IP-Q and IPx-Q as the
// resins and gives no dielectric data at all (checked: the text contains no
// occurrence of "permittivity" or "dielectric constant").  The values below are
// therefore either manufacturer specifications at a stated frequency, or
// peer-reviewed measurements at a stated frequency, or an explicitly
// PROVISIONAL nominal value -- never a number invented to fill a field, and
// never n^2 from an optical refractive index (optical n probes the electronic
// polarisability alone and misses the orientational contribution that dominates
// the static value of an epoxy; for SU-8 that route would give about 2.6 and be
// wrong by roughly a quarter).
//
// A provisional value is not a licence to claim validation.  It is a licence to
// run a SENSITIVITY study over a range that the sources justify, and to report
// how much the answer moves inside it.  That is what P2b does.
//
// ADDING IP-Q AND IPx-Q LATER
//
// A material is a value plus its provenance, held in one struct, looked up by
// name.  IP-Q and IPx-Q are already registered, with no number and status
// Unknown, so that asking for one fails loudly and can be satisfied by
// supplying a value in the configuration -- no code change, no new type.

/// How much weight the number carries.  Ordered from strongest to weakest.
enum class MaterialStatus {
  Exact = 0,        ///< eps_r = 1 for vacuum; not a measurement, a definition
  IdealConductor,   ///< permittivity does not enter; the region is equipotential
  Measured,         ///< measured on THIS material in THIS process state
  ManufacturerSpec, ///< manufacturer data sheet, at the stated conditions
  Literature,       ///< peer-reviewed measurement on nominally the same material
  Provisional,      ///< nominal value chosen because no usable one exists
  Unknown,          ///< registered but carries no number; using it is an error
};
const char* to_string(MaterialStatus s);

/// True where the value may back a quantitative claim without a caveat.
inline bool is_defensible(MaterialStatus s) {
  return s == MaterialStatus::Exact || s == MaterialStatus::IdealConductor ||
         s == MaterialStatus::Measured;
}

// ---------------------------------------------------------------------------

struct Material {
  std::string name;
  /// Static relative permittivity.  Meaningless (and left at 0) when
  /// ideal_conductor is true or the status is Unknown.
  Real relative_permittivity{0.0};
  MaterialStatus status{MaterialStatus::Unknown};
  /// The region is an equipotential perfect conductor; eps_r is not used.
  bool ideal_conductor{false};

  std::string source;      ///< where the number comes from, precisely enough to find
  std::string conditions;  ///< frequency, humidity, cure state of THAT measurement
  std::string caveat;      ///< what must not be concluded from it
  /// Range over which a sensitivity study is justified by the sources.
  /// Equal endpoints mean no study is called for.
  Real eps_r_low{0.0}, eps_r_high{0.0};

  bool has_range() const { return eps_r_high > eps_r_low; }
  /// Throws unless the material can be used as a dielectric right now.
  Real permittivity_or_throw() const;
  void print(std::FILE* out) const;
};

// ---------------------------------------------------------------------------

/// Built-in materials, by name.  Case-sensitive, lowercase, no spaces.
///
///   vacuum      eps_r = 1, exact
///   su8         cured SU-8 photo-epoxy -- PROVISIONAL static value, see below
///   ip-q        registered, no value  -- supply one to use it
///   ipx-q       registered, no value  -- supply one to use it
///   peek        registered, no value  -- supply one to use it.  The built
///               reservoir body is PEEK; no static permittivity for THAT part
///               in THIS state is in hand, so it carries no number either.
///   emibf4      ionic liquid, treated as an ideal conductor in P2b
///   metal       the extractor metallisation; ideal conductor, zero thickness
class MaterialLibrary {
 public:
  MaterialLibrary();                       ///< builds the table above
  const Material& get(const std::string& name) const;   ///< throws if unknown
  bool has(const std::string& name) const;
  /// Insert or replace.  This is how a configuration file supplies a value for
  /// IP-Q or overrides the provisional SU-8 number.
  void set(const Material& m);
  /// Replace only the permittivity and its provenance, keeping the name.
  void override_permittivity(const std::string& name, Real eps_r,
                             MaterialStatus status, const std::string& source);
  const std::vector<Material>& all() const { return m_; }
  void print(std::FILE* out) const;
  void write_csv(const std::string& path) const;

 private:
  std::vector<Material> m_;
  Material* find(const std::string& name);
};

// ---------------------------------------------------------------------------

/// The four field regions of the P2b model plus the metallisation, resolved to
/// concrete materials.  Every solver input is here; nothing is looked up by
/// name inside the solver.
struct DielectricMaterials {
  Material vacuum;             ///< the surrounding vacuum
  Material emitter_dielectric; ///< SU-8 (or IP-Q / IPx-Q later)
  Material liquid;             ///< ionic liquid; ideal conductor in P2b
  Material extractor_carrier;  ///< polymer carrier of the extraction electrode
  Material metallisation;      ///< the electrode film; ideal conductor
  /// Body of the liquid reservoir.  A DIELECTRIC, never an electrode.  The
  /// built part is PEEK; because no sourced static value for it is in hand, the
  /// reference assignment falls back to the emitter resin and SAYS SO in the
  /// caveat.  Supplying material.peek.relative_permittivity together with
  /// material.peek.source in the configuration replaces it without a code
  /// change -- the same mechanism IP-Q and IPx-Q use.
  Material reservoir_body;

  /// SU-8 emitter, SU-8 extractor carrier, EMI-BF4 liquid.  The P2b reference.
  static DielectricMaterials reference(const MaterialLibrary& lib);
  /// Every dielectric set to eps_r = 1.  Used ONLY to reproduce the vacuum
  /// problem the BEM can solve independently; not a physical configuration.
  static DielectricMaterials all_vacuum(const MaterialLibrary& lib);

  void print(std::FILE* out) const;
  void write_csv(const std::string& path) const;
  /// Throws unless every dielectric carries a usable number.
  void check_usable() const;
};

// ---------------------------------------------------------------------------
// The finite conductivity of the ionic liquid -- deliberately NOT modelled
// ---------------------------------------------------------------------------
//
// P2b treats the ionic liquid as an ideal equipotential conductor.  That is an
// approximation with a known validity condition and a known failure mode, and
// it is written down here so that no later phase has to rediscover it.
//
// An ionic liquid such as EMI-BF4 has a conductivity of order 1 S/m and a
// static permittivity of order 10.  Its charge relaxation time is
// tau_e = eps / sigma ~ 1e-10 s.  Every process P2b describes is static, so
// tau_e is negligible against it and the equipotential assumption holds: the
// liquid screens its interior completely and its surface is an equipotential.
//
// It stops holding as soon as current flows through the meniscus.  In emission
// the ohmic drop along the cone, and the surface charge that a finite
// relaxation time cannot replenish instantly at the apex, both matter -- that
// is the leaky-dielectric regime, and it is what sets the apex field once ions
// leave.  P3/P4 must replace the Dirichlet condition on the liquid by a
// conduction model.  Any statement about emission taken from a P2b field is
// therefore a statement about perfect screening, not a prediction.
namespace liquid_model {
inline constexpr Real kTypicalConductivity = 1.0;      ///< [S/m], order of magnitude
inline constexpr Real kTypicalPermittivity = 10.0;     ///< [-], order of magnitude
/// eps0 * eps_r / sigma [s].  About 1e-10 s; static problems do not resolve it.
Real charge_relaxation_time(Real sigma = kTypicalConductivity,
                            Real eps_r = kTypicalPermittivity);
const char* why_ideal_conductor_is_admissible_in_p2b();
}  // namespace liquid_model

}  // namespace es
