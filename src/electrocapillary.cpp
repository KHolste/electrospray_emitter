#include "es/electrocapillary.hpp"

#include "es/load_projection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "es/constants.hpp"

namespace es {

using constants::eps0;
using constants::pi;

// ===========================================================================
// Status strings
// ===========================================================================

const char* to_string(GateVerdict v) {
  switch (v) {
    case GateVerdict::NotAttempted: return "NotAttempted";
    case GateVerdict::Passed: return "Passed";
    case GateVerdict::FailedTotalForce: return "FailedTotalForce";
    case GateVerdict::FailedEdgeFarLoad: return "FailedEdgeFarLoad";
    case GateVerdict::FailedExclusion: return "FailedExclusion";
    case GateVerdict::FailedNotIntegrable: return "FailedNotIntegrable";
    case GateVerdict::FailedNoData: return "FailedNoData";
  }
  return "?";
}

const char* explain(GateVerdict v) {
  switch (v) {
    case GateVerdict::NotAttempted:
      return "Das Kanten-Gate wurde nicht gerechnet.";
    case GateVerdict::Passed:
      return "Die elektrische Flaechenlast ist an der unverrundeten Kante integrierbar und die "
             "schwache, flaechengemittelte Projektion konvergiert. Der punktweise Kantenwert "
             "bleibt ausdruecklich nicht konvergiert und wird nirgends verwendet.";
    case GateVerdict::FailedTotalForce:
      return "Die integrierte Maxwell-Kraft konvergiert mit der Netzverfeinerung nicht innerhalb "
             "der vorab festgelegten Grenze.";
    case GateVerdict::FailedEdgeFarLoad:
      return "Die flaechengemittelte Last ausserhalb der Kantenzone konvergiert nicht innerhalb "
             "der vorab festgelegten Grenze.";
    case GateVerdict::FailedExclusion:
      return "Die beiden unabhaengigen Extrapolationen der Grenzkraft -- ueber die Netzstufen "
             "und ueber die Ausschlussdistanz mit dem gefitteten Exponenten -- stimmen nicht "
             "ueberein. Der Beitrag der Kantenzone ist damit nicht als summierbar nachgewiesen.";
    case GateVerdict::FailedNotIntegrable:
      return "Der gemessene Singularitaetsexponent ist nicht integrierbar (beta <= -1).";
    case GateVerdict::FailedNoData:
      return "Zu wenige Netzstufen oder Segmente fuer ein Urteil.";
  }
  return "?";
}

const char* to_string(CouplingStatus s) {
  switch (s) {
    case CouplingStatus::NotAttempted: return "NotAttempted";
    case CouplingStatus::Converged: return "Converged";
    case CouplingStatus::ElectrostaticFailure: return "ElectrostaticFailure";
    case CouplingStatus::MeshInvalid: return "MeshInvalid";
    case CouplingStatus::MechanicalResidualNotConverged:
      return "MechanicalResidualNotConverged";
    case CouplingStatus::ContinuationStepTooSmall: return "ContinuationStepTooSmall";
    case CouplingStatus::CapillaryRangeExceeded: return "CapillaryRangeExceeded";
    case CouplingStatus::EdgeLoadNotWellPosed: return "EdgeLoadNotWellPosed";
    case CouplingStatus::MultipleSolutionsDetected: return "MultipleSolutionsDetected";
  }
  return "?";
}

const char* explain(CouplingStatus s) {
  switch (s) {
    case CouplingStatus::NotAttempted:
      return "Nicht gerechnet.";
    case CouplingStatus::Converged:
      return "Form und Feldlast sind gemeinsam auf die vorab festgelegten Grenzen konvergiert.";
    case CouplingStatus::ElectrostaticFailure:
      return "Der Elektrostatiksolver ist gescheitert.";
    case CouplingStatus::MeshInvalid:
      return "Zu dieser Meniskusform liess sich kein gueltiges konformes Volumennetz bauen.";
    case CouplingStatus::MechanicalResidualNotConverged:
      return "Die Fixpunktiteration aus Form und Last hat ihre Grenzen nicht erreicht.";
    case CouplingStatus::ContinuationStepTooSmall:
      return "Die Fortsetzung liess sich nicht weiterfuehren: die Schrittweite ist unter ihre "
             "Untergrenze gefallen. Das ist die Stelle, an der DIESER Loeser stehen bleibt -- "
             "kein Emissionsbeginn, kein Taylor-Kegel, keine Stabilitaetsaussage.";
    case CouplingStatus::CapillaryRangeExceeded:
      return "Zu dieser Last existiert keine an der Kante gepinnte glatte Oberflaeche mehr: die "
             "Meridiankurve wird senkrecht, bevor sie den Kontaktradius erreicht.";
    case CouplingStatus::EdgeLoadNotWellPosed:
      return "Das Kanten-Gate ist nicht bestanden; es wurde nichts gekoppelt.";
    case CouplingStatus::MultipleSolutionsDetected:
      return "Mehrere statische Formen erfuellen dieselben Daten. Keine davon wird "
             "stillschweigend ausgewaehlt und keine als stabil bezeichnet.";
  }
  return "?";
}

// ===========================================================================
// FreeSurface
// ===========================================================================

FreeSurface FreeSurface::flat_surface(Real a, Real z_contact) {
  FreeSurface s;
  s.flat = true;
  s.contact_radius = a;
  s.contact_z = z_contact;
  s.apex_height = 0.0;
  s.arclength = a;
  s.nodes = {{0.0, z_contact}, {a, z_contact}};
  s.psi = {0.0, 0.0};
  return s;
}

FreeSurface FreeSurface::from(const CapillaryMeniscus& m) {
  if (!is_usable(m.status))
    throw std::runtime_error("FreeSurface::from: die Meniskusloesung traegt den Status " +
                             std::string(to_string(m.status)) + " und ist keine Form.");
  FreeSurface s;
  s.nodes = m.nodes;
  s.psi = m.psi;
  s.contact_radius = m.contact_radius;
  s.contact_z = m.contact_z;
  s.apex_height = m.apex_height;
  s.arclength = m.arclength;
  s.flat = (m.delta_p_exit == 0.0 && m.load.empty());
  return s;
}

namespace {

/// Cubic Hermite value and derivative on [0,1].
inline void hermite(Real y0, Real y1, Real d0, Real d1, Real t, Real* y, Real* dy) {
  const Real t2 = t * t, t3 = t2 * t;
  const Real h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + t;
  const Real h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
  *y = h00 * y0 + h10 * d0 + h01 * y1 + h11 * d1;
  const Real g00 = 6 * t2 - 6 * t, g10 = 3 * t2 - 4 * t + 1;
  const Real g01 = -6 * t2 + 6 * t, g11 = 3 * t2 - 2 * t;
  *dy = g00 * y0 + g10 * d0 + g01 * y1 + g11 * d1;
}

}  // namespace

/// Interval index and local parameter t for a given radius.  r(s) is strictly
/// increasing while |psi| < 90 degrees, which is the whole representable range.
static void locate_on_surface(const FreeSurface& s, Real r, std::size_t* k, Real* t) {
  const std::size_t n = s.nodes.size();
  const Real rr = std::min(s.contact_radius, std::max(0.0, r));
  std::size_t lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const std::size_t mid = (lo + hi) / 2;
    (s.nodes[mid].r <= rr ? lo : hi) = mid;
  }
  const Real ds = s.arclength / static_cast<Real>(n - 1);
  const Real r0 = s.nodes[lo].r, r1 = s.nodes[lo + 1].r;
  const Real d0 = ds * std::cos(s.psi[lo]), d1 = ds * std::cos(s.psi[lo + 1]);
  // Newton on the Hermite for r, safeguarded by the bracket [0,1].
  Real tt = (r1 > r0) ? (rr - r0) / (r1 - r0) : 0.0;
  tt = std::min(1.0, std::max(0.0, tt));
  for (int it = 0; it < 30; ++it) {
    Real y, dy;
    hermite(r0, r1, d0, d1, tt, &y, &dy);
    const Real f = y - rr;
    if (std::abs(f) <= 1e-15 * s.contact_radius || dy <= 0.0) break;
    Real tn = tt - f / dy;
    if (!(tn >= 0.0 && tn <= 1.0)) tn = std::min(1.0, std::max(0.0, tt - f / dy));
    if (std::abs(tn - tt) <= 1e-16) { tt = tn; break; }
    tt = tn;
  }
  *k = lo;
  *t = tt;
}

Real FreeSurface::z_at_radius(Real r) const {
  if (flat) return contact_z;
  std::size_t k;
  Real t;
  locate_on_surface(*this, r, &k, &t);
  const Real ds = arclength / static_cast<Real>(nodes.size() - 1);
  Real z, dz;
  hermite(nodes[k].z, nodes[k + 1].z, -ds * std::sin(psi[k]), -ds * std::sin(psi[k + 1]), t, &z,
          &dz);
  return z;
}

Real FreeSurface::psi_at_radius(Real r) const {
  if (flat) return 0.0;
  std::size_t k;
  Real t;
  locate_on_surface(*this, r, &k, &t);
  const Real ds = arclength / static_cast<Real>(nodes.size() - 1);
  Real z, dz, rr, dr;
  hermite(nodes[k].z, nodes[k + 1].z, -ds * std::sin(psi[k]), -ds * std::sin(psi[k + 1]), t, &z,
          &dz);
  hermite(nodes[k].r, nodes[k + 1].r, ds * std::cos(psi[k]), ds * std::cos(psi[k + 1]), t, &rr,
          &dr);
  return std::atan2(-dz, dr);
}

Real FreeSurface::s_at_radius(Real r) const {
  if (flat) return std::min(contact_radius, std::max(0.0, r));
  std::size_t k;
  Real t;
  locate_on_surface(*this, r, &k, &t);
  const Real ds = arclength / static_cast<Real>(nodes.size() - 1);
  return (static_cast<Real>(k) + t) * ds;
}

Real FreeSurface::revolved_volume() const {
  if (flat) return 0.0;
  std::vector<Vec2> loop = nodes;
  loop.push_back({0.0, contact_z});
  return -es::revolved_volume(loop);
}

Real FreeSurface::revolved_area() const {
  if (flat) return pi * contact_radius * contact_radius;
  return es::revolved_area(nodes);
}

// ===========================================================================
// The moving mesh
// ===========================================================================

bool MeniscusMeshQuality::ok() const {
  return inverted_cells == 0 && min_jacobian > 0.0;
}

void MeniscusMeshQuality::print(std::FILE* out) const {
  std::fprintf(out, "Netzqualitaet des beweglichen Netzes\n");
  std::fprintf(out, "  Verformungsband       : z = %.4g .. %.4g m\n", band_z_lo, band_z_hi);
  std::fprintf(out, "  verformte Zellen      : %lld\n", static_cast<long long>(deformed_cells));
  std::fprintf(out, "  min. Jacobi-Det.      : %.6e (muss > 0 sein)\n", min_jacobian);
  std::fprintf(out, "  invertierte Zellen    : %lld (Sollwert 0)\n",
               static_cast<long long>(inverted_cells));
  std::fprintf(out, "  groesstes Seitenverh. : %.4g\n", max_cell_aspect);
  std::fprintf(out, "  groesste Scherung     : %.4g\n", max_shear);
  std::fprintf(out, "  Kontaktlinie          : dr/a = %.3e, dz/a = %.3e\n", contact_radius_error,
               contact_z_error);
  std::fprintf(out, "  Apex                  : dz/a = %.3e\n", apex_error);
  std::fprintf(out, "  Oberflaeche           : max dz/a = %.3e\n", surface_error);
  std::fprintf(out, "  Fluessigkeitsvolumen  : Netz %.9e m^3, Referenz %.9e m^3, rel. %.3e\n",
               liquid_volume_mesh, liquid_volume_reference, liquid_volume_error);
}

MeniscusMesh build_meniscus_mesh(const DielectricDeviceParameters& p, const FreeSurface& surface) {
  MeniscusMesh out;
  out.device = build_volume_mesh(p);
  out.surface = surface;
  DeviceVolumeMesh& m = out.device;
  QuadMesh& g = m.grid;

  const Real a = m.r_bore;
  if (std::abs(surface.contact_radius - a) > 1e-12 * a)
    throw std::runtime_error(
        "build_meniscus_mesh: der Pinningradius der Form und der Bohrungsradius des Netzes "
        "stimmen nicht ueberein. Es gibt keine zweite Geometriebeschreibung.");

  out.j_surface = m.j_tip;
  out.i_contact = m.i_bore;
  out.quality.band_z_lo = 0.0;
  out.quality.band_z_hi = 0.0;

  if (surface.flat) {
    // Nothing moves.  The mesh is bitwise the P2c one, which is what makes the
    // zero-field cross-check a check rather than a comparison of two meshes.
    out.quality.min_jacobian = g.min_jacobian();
    out.quality.liquid_volume_mesh = m.revolved_volume_of(Region::Liquid);
    out.quality.liquid_volume_reference =
        m.analytic_volume_of(Region::Liquid) + surface.revolved_volume();
    out.quality.liquid_volume_error =
        std::abs(out.quality.liquid_volume_mesh - out.quality.liquid_volume_reference) /
        std::max(out.quality.liquid_volume_reference, 1e-300);
    return out;
  }

  // --- the band edges: existing rows, no new levels -------------------------
  const Real want = meniscus_mesh::kBandFactor * a;
  Index j_lo = -1, j_hi = -1;
  for (Index j = m.j_tip; j >= 0; --j)
    if (g.z_of_row(j) <= -want) { j_lo = j; break; }
  for (Index j = m.j_tip; j < g.nz; ++j)
    if (g.z_of_row(j) >= want) { j_hi = j; break; }
  if (j_lo < 0 || j_hi < 0)
    throw std::runtime_error(
        "build_meniscus_mesh: im Netz gibt es keine Zeile mindestens " +
        std::to_string(meniscus_mesh::kBandFactor) +
        " Bohrungsradien ueber bzw. unter der Austrittsebene; das Verformungsband hat keine "
        "Raender.");
  const Real z_lo = g.z_of_row(j_lo), z_hi = g.z_of_row(j_hi);
  out.quality.band_z_lo = z_lo;
  out.quality.band_z_hi = z_hi;

  // --- the axial warp on the bore columns -----------------------------------
  std::vector<Real> z_ref(static_cast<std::size_t>(g.nz));
  for (Index j = 0; j < g.nz; ++j) z_ref[static_cast<std::size_t>(j)] = g.z_of_row(j);

  Real worst_surface = 0.0;
  for (Index i = 0; i <= m.i_bore; ++i) {
    const Real r = g.at(i, m.j_tip).r;
    // The contact node is pinned exactly; nowhere else does an interpolated
    // value decide where a boundary condition sits.
    const Real zm = (i == m.i_bore) ? surface.contact_z : surface.z_at_radius(r);
    if (!(zm > z_lo && zm < z_hi))
      throw std::runtime_error(
          "build_meniscus_mesh: die Meniskusform verlaesst das Verformungsband (z = " +
          std::to_string(zm) + " m, Band " + std::to_string(z_lo) + " .. " +
          std::to_string(z_hi) + " m).");
    worst_surface = std::max(worst_surface, std::abs(zm - surface.z_at_radius(r)));
    // The contact column is pinned: its meniscus height IS the contact plane,
    // so the warp is the identity there.  Applied as arithmetic it is the
    // identity only up to round-off, and the promise that nothing outside the
    // bore moves would then be true to 1e-21 instead of exactly.  It is skipped
    // rather than computed.
    if (zm == surface.contact_z) continue;
    for (Index j = j_lo + 1; j < j_hi; ++j) {
      const Real zr = z_ref[static_cast<std::size_t>(j)];
      const Real zn = (zr <= 0.0) ? z_lo + (zr - z_lo) * (zm - z_lo) / (0.0 - z_lo)
                                  : zm + zr * (z_hi - zm) / z_hi;
      g.nodes[static_cast<std::size_t>(g.node(i, j))].z = zn;
    }
    g.nodes[static_cast<std::size_t>(g.node(i, m.j_tip))].z = zm;
  }

  // --- quality, measured ----------------------------------------------------
  MeniscusMeshQuality& q = out.quality;
  q.surface_error = worst_surface / a;
  q.min_jacobian = g.min_jacobian();
  q.deformed_cells = (m.i_bore) * (j_hi - j_lo);
  for (Index j = j_lo; j < j_hi; ++j)
    for (Index i = 0; i < m.i_bore; ++i) {
      const auto c = g.cell_nodes(i, j);
      const Vec2 p00 = g.nodes[static_cast<std::size_t>(c[0])];
      const Vec2 p10 = g.nodes[static_cast<std::size_t>(c[1])];
      const Vec2 p11 = g.nodes[static_cast<std::size_t>(c[2])];
      const Vec2 p01 = g.nodes[static_cast<std::size_t>(c[3])];
      const Real w = norm(p10 - p00);
      const Real h0 = norm(p01 - p00), h1 = norm(p11 - p10);
      const Real hmin = std::min(h0, h1), hmax = std::max(h0, h1);
      if (hmin > 0.0) q.max_cell_aspect = std::max(q.max_cell_aspect, std::max(w / hmin, hmin / w));
      const Real shear = std::abs((p10.z - p00.z));
      if (hmax > 0.0) q.max_shear = std::max(q.max_shear, shear / hmax);
      if (g.cell_meridian_area(i, j) <= 0.0) ++q.inverted_cells;
    }
  if (q.min_jacobian <= 0.0 || q.inverted_cells > 0)
    throw std::runtime_error(
        "build_meniscus_mesh: das verformte Netz enthaelt entartete Zellen (min. Jacobi-Det. " +
        std::to_string(q.min_jacobian) + ", " + std::to_string(q.inverted_cells) +
        " invertierte Zellen).");

  q.contact_radius_error = std::abs(g.at(m.i_bore, m.j_tip).r - a) / a;
  q.contact_z_error = std::abs(g.at(m.i_bore, m.j_tip).z - surface.contact_z) / a;
  q.apex_error = std::abs(g.at(0, m.j_tip).z - (surface.contact_z + surface.apex_height)) / a;
  q.liquid_volume_mesh = m.revolved_volume_of(Region::Liquid);
  q.liquid_volume_reference = m.analytic_volume_of(Region::Liquid) + surface.revolved_volume();
  q.liquid_volume_error = std::abs(q.liquid_volume_mesh - q.liquid_volume_reference) /
                          std::max(std::abs(q.liquid_volume_reference), 1e-300);
  return out;
}

// ---------------------------------------------------------------------------
// Point location on the deformed mesh
// ---------------------------------------------------------------------------

bool locate_meniscus(const MeniscusMesh& mm, Vec2 x, Index* io, Index* jo, Real* xio,
                     Real* etao) {
  const DeviceVolumeMesh& m = mm.device;
  const QuadMesh& g = m.grid;
  if (g.nr < 2 || g.nz < 2) return false;
  const Index ib = m.i_bore;

  if (x.r >= g.at(ib, m.j_tip).r) {
    // Outside the bore the rows are level: bisect on the z of the column at the
    // bore wall, which is untouched by the deformation.
    auto row_z = [&](Index j) { return g.at(ib, j).z; };
    if (x.z < row_z(0) || x.z > row_z(g.nz - 1)) return false;
    Index lo = 0, hi = g.nz - 1;
    while (hi - lo > 1) {
      const Index mid = (lo + hi) / 2;
      (row_z(mid) <= x.z ? lo : hi) = mid;
    }
    const Index j = lo;
    const Real za = row_z(j), zb = row_z(j + 1);
    Real eta = (zb > za) ? (x.z - za) / (zb - za) : 0.0;
    eta = std::min(1.0, std::max(0.0, eta));
    auto radius = [&](Index i) { return (1.0 - eta) * g.at(i, j).r + eta * g.at(i, j + 1).r; };
    if (x.r < radius(0) - 1e-15 || x.r > radius(g.nr - 1) + 1e-15) return false;
    Index rlo = 0, rhi = g.nr - 1;
    while (rhi - rlo > 1) {
      const Index mid = (rlo + rhi) / 2;
      (radius(mid) <= x.r ? rlo : rhi) = mid;
    }
    const Index i = rlo;
    const Real ra = radius(i), rb = radius(i + 1);
    Real xi = (rb > ra) ? (x.r - ra) / (rb - ra) : 0.0;
    *io = i;
    *jo = j;
    *xio = std::min(1.0, std::max(0.0, xi));
    *etao = eta;
    return true;
  }

  // Inside the bore the columns are vertical: the radius depends on i alone, so
  // the column comes first and the row follows from the interpolated height of
  // the cell edges at that radius.  Exact for bilinear cells.
  Index rlo = 0, rhi = ib;
  while (rhi - rlo > 1) {
    const Index mid = (rlo + rhi) / 2;
    (g.at(mid, m.j_tip).r <= x.r ? rlo : rhi) = mid;
  }
  const Index i = rlo;
  const Real ra = g.at(i, m.j_tip).r, rb = g.at(i + 1, m.j_tip).r;
  const Real xi = (rb > ra) ? std::min(1.0, std::max(0.0, (x.r - ra) / (rb - ra))) : 0.0;
  auto row_z_at = [&](Index j) { return (1.0 - xi) * g.at(i, j).z + xi * g.at(i + 1, j).z; };
  if (x.z < row_z_at(0) || x.z > row_z_at(g.nz - 1)) return false;
  Index lo = 0, hi = g.nz - 1;
  while (hi - lo > 1) {
    const Index mid = (lo + hi) / 2;
    (row_z_at(mid) <= x.z ? lo : hi) = mid;
  }
  const Index j = lo;
  const Real za = row_z_at(j), zb = row_z_at(j + 1);
  Real eta = (zb > za) ? (x.z - za) / (zb - za) : 0.0;
  *io = i;
  *jo = j;
  *xio = xi;
  *etao = std::min(1.0, std::max(0.0, eta));
  return true;
}

Real potential_at_meniscus(const MeniscusMesh& mm, const std::vector<Real>& phi, Vec2 x) {
  Index i, j;
  Real xi, eta;
  if (!locate_meniscus(mm, x, &i, &j, &xi, &eta))
    throw std::runtime_error("potential_at_meniscus: Punkt liegt ausserhalb des Netzes");
  return potential_in_cell(mm.device.grid, phi, i, j, xi, eta);
}

Vec2 field_recovered_at_meniscus(const MeniscusMesh& mm, const std::vector<Real>& phi,
                                 const std::vector<Real>& eps_r, const std::vector<char>& active,
                                 Vec2 x) {
  Index i, j;
  Real xi, eta;
  if (!locate_meniscus(mm, x, &i, &j, &xi, &eta))
    throw std::runtime_error("field_recovered_at_meniscus: Punkt liegt ausserhalb des Netzes");
  const QuadMesh& g = mm.device.grid;
  const Index c = g.cell(i, j);
  if (!active[static_cast<std::size_t>(c)])
    throw std::runtime_error("field_recovered_at_meniscus: Punkt liegt in einer inaktiven Zelle");
  const Real eps_select = eps_r[static_cast<std::size_t>(c)];
  const Vec2 e00 = field_recovered_at_node(g, phi, eps_r, active, i, j, eps_select);
  const Vec2 e10 = field_recovered_at_node(g, phi, eps_r, active, i + 1, j, eps_select);
  const Vec2 e11 = field_recovered_at_node(g, phi, eps_r, active, i + 1, j + 1, eps_select);
  const Vec2 e01 = field_recovered_at_node(g, phi, eps_r, active, i, j + 1, eps_select);
  const Real n0 = (1 - xi) * (1 - eta), n1 = xi * (1 - eta), n2 = xi * eta, n3 = (1 - xi) * eta;
  return n0 * e00 + n1 * e10 + n2 * e11 + n3 * e01;
}

// ===========================================================================
// The Maxwell load
// ===========================================================================

Real MaxwellLoad::force_beyond(Real d_exclude) const {
  // DEFECT FOUND AND FIXED: counting whole segments by their midpoint quantises
  // this measure by the force content of one segment, which near the edge is
  // several per cent of the total -- so the quantity moved with the mesh for a
  // reason that had nothing to do with the field.  The segment that straddles
  // the exclusion distance is now split by the fraction of its arclength that
  // lies beyond it.  tests/test_electrocapillary.cpp pins the split.
  Real f = 0.0;
  for (std::size_t k = 0; k < seg_force.size(); ++k) {
    const Real d_far = std::max(node_d_edge[k], node_d_edge[k + 1]);
    const Real d_near = std::min(node_d_edge[k], node_d_edge[k + 1]);
    Real w = 1.0;
    if (d_exclude >= d_far)
      w = 0.0;
    else if (d_exclude > d_near)
      w = (d_far - d_exclude) / std::max(d_far - d_near, 1e-300);
    f += w * seg_force[k];
  }
  return f;
}

Real MaxwellLoad::pressure_beyond(Real d_exclude) const {
  Real p = 0.0;
  for (std::size_t k = 0; k < seg_pressure.size(); ++k)
    if (seg_d_mid[k] >= d_exclude) p = std::max(p, seg_pressure[k]);
  return p;
}

Real MaxwellLoad::pressure_at_tau(Real tau) const {
  if (seg_pressure.empty()) return 0.0;
  const Real t = std::min(1.0, std::max(0.0, tau));
  // Piecewise constant per surface segment: this is exactly the conservative
  // projection, so the integral of the returned load over the surface is the
  // integrated Maxwell force and not an interpolation of it.
  std::size_t lo = 0, hi = seg_tau1.size() - 1;
  while (lo < hi && seg_tau1[lo] < t) ++lo;
  return seg_pressure[lo];
}

MaxwellLoad maxwell_load(const MeniscusMesh& mm, const DielectricSolution& sol,
                         Real gamma_over_a) {
  MaxwellLoad out;
  out.gamma_over_a = gamma_over_a;
  const DeviceVolumeMesh& m = mm.device;
  const QuadMesh& g = m.grid;
  const Index j = mm.j_surface;
  const Index ib = mm.i_contact;
  const FreeSurface& fs = mm.surface;
  const Real L = fs.arclength;

  const Real eps_vac = 1.0;
  for (Index i = 0; i <= ib; ++i) {
    const Vec2 x = g.at(i, j);
    // One-sided limit from the VACUUM: the recovery averages only cells whose
    // permittivity matches, and the liquid cells below are inactive, so nothing
    // from inside the conductor or from the polymer enters.
    const Vec2 E =
        field_recovered_at_node(g, sol.fem.phi, sol.cell_eps_r, sol.cell_active, i, j, eps_vac);
    const Real ps = fs.psi_at_radius(x.r);
    const Vec2 n{std::sin(ps), std::cos(ps)};
    const Real En = dot(E, n);
    const Real Emag = norm(E);
    const Real s = fs.s_at_radius(x.r);
    out.node_r.push_back(x.r);
    out.node_z.push_back(x.z);
    out.node_s.push_back(s);
    out.node_tau.push_back(L > 0.0 ? s / L : 0.0);
    out.node_d_edge.push_back(std::max(0.0, L - s));
    out.node_En.push_back(En);
    out.node_pM.push_back(0.5 * eps0 * En * En);
    out.node_tangential_fraction.push_back(
        Emag > 0.0 ? std::sqrt(std::max(0.0, Emag * Emag - En * En)) / Emag : 0.0);
  }

  assemble_load_segments(out, fs);
  return out;
}

// --- segments: the conservative projection ----------------------------------
//
// Trapezoidal integration of the nodal load against the revolved area element.
// Summing the segments returns the total by construction, so the projection
// preserves the integrated Maxwell force exactly; the segment mean of an
// integrable singularity stays finite where the pointwise value does not.
//
// Lifted out of maxwell_load() in P0 so that the manufactured loads of
// load_projection.hpp run THIS quadrature and not a copy of it.  Nothing about
// it changed; it does not know where the nodal values came from.
void assemble_load_segments(MaxwellLoad& out, const FreeSurface& fs) {
  out.seg_tau0.clear();
  out.seg_tau1.clear();
  out.seg_tau_mid.clear();
  out.seg_d_mid.clear();
  out.seg_area.clear();
  out.seg_force.clear();
  out.seg_pressure.clear();
  out.total_force = 0.0;
  out.axial_force = 0.0;
  for (std::size_t k = 0; k + 1 < out.node_r.size(); ++k) {
    const Real r0 = out.node_r[k], r1 = out.node_r[k + 1];
    const Real z0 = out.node_z[k], z1 = out.node_z[k + 1];
    const Real p0 = out.node_pM[k], p1 = out.node_pM[k + 1];
    const Real ds = std::hypot(r1 - r0, z1 - z0);
    // int 2 pi r ds with r linear, and int p 2 pi r ds with p linear too.
    const Real area = 2.0 * pi * 0.5 * (r0 + r1) * ds;
    const Real force =
        2.0 * pi * ds * ((2.0 * r0 * p0 + r0 * p1 + r1 * p0 + 2.0 * r1 * p1) / 6.0);
    const Real nz = std::cos(0.5 * (fs.psi_at_radius(r0) + fs.psi_at_radius(r1)));
    out.seg_tau0.push_back(out.node_tau[k]);
    out.seg_tau1.push_back(out.node_tau[k + 1]);
    out.seg_tau_mid.push_back(0.5 * (out.node_tau[k] + out.node_tau[k + 1]));
    out.seg_d_mid.push_back(0.5 * (out.node_d_edge[k] + out.node_d_edge[k + 1]));
    out.seg_area.push_back(area);
    out.seg_force.push_back(force);
    out.seg_pressure.push_back(area > 0.0 ? force / area : 0.0);
    out.total_force += force;
    out.axial_force += force * nz;
  }
}

Real MaxwellLoad::pressure_at_distance(Real d) const {
  const std::size_t n = seg_pressure.size();
  if (n == 0) return 0.0;
  // seg_d_mid decreases from the apex towards the edge, so scan for the bracket
  // instead of assuming a direction.
  for (std::size_t k = 0; k + 1 < n; ++k) {
    const Real d0 = seg_d_mid[k], d1 = seg_d_mid[k + 1];
    const Real hi = std::max(d0, d1), lo = std::min(d0, d1);
    if (d <= hi && d >= lo) {
      const Real w = (hi > lo) ? (d - lo) / (hi - lo) : 0.0;
      const Real p_lo = (d0 < d1) ? seg_pressure[k] : seg_pressure[k + 1];
      const Real p_hi = (d0 < d1) ? seg_pressure[k + 1] : seg_pressure[k];
      return (1.0 - w) * p_lo + w * p_hi;
    }
  }
  return 0.0;
}

// ===========================================================================
// The edge gate
// ===========================================================================

namespace {

/// Least-squares fit of log p = A + beta log d over the segments whose midpoint
/// lies in [d_lo, d_hi].  Returns false if fewer than four points qualify.
bool fit_exponent(const MaxwellLoad& L, Real d_lo, Real d_hi, Real* beta, Real* r2,
                  Index* used) {
  std::vector<Real> x, y;
  for (std::size_t k = 0; k < L.seg_pressure.size(); ++k) {
    const Real d = L.seg_d_mid[k], p = L.seg_pressure[k];
    if (d < d_lo || d > d_hi || !(p > 0.0)) continue;
    x.push_back(std::log(d));
    y.push_back(std::log(p));
  }
  *used = static_cast<Index>(x.size());
  if (x.size() < 4) return false;
  const Real n = static_cast<Real>(x.size());
  Real sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    sx += x[k];
    sy += y[k];
    sxx += x[k] * x[k];
    sxy += x[k] * y[k];
    syy += y[k] * y[k];
  }
  const Real den = n * sxx - sx * sx;
  if (std::abs(den) < 1e-300) return false;
  *beta = (n * sxy - sx * sy) / den;
  const Real num = n * sxy - sx * sy;
  const Real d2 = (n * sxx - sx * sx) * (n * syy - sy * sy);
  *r2 = (d2 > 0.0) ? (num * num) / d2 : 0.0;
  return true;
}

}  // namespace

EdgeGateResult run_edge_gate(const DielectricDeviceParameters& base,
                             const DielectricMaterials& materials, Real V_emitter,
                             Real V_extractor, Metallisation metallisation, FarField far_field,
                             const FreeSurface& surface, const std::string& tag, Real Pi,
                             const std::vector<int>& levels, Real gamma_over_a,
                             std::size_t memory_cap_bytes) {
  EdgeGateResult gate;
  gate.shape_tag = tag;
  gate.Pi = Pi;
  const Real a = surface.contact_radius;

  for (int lvl : levels) {
    DielectricDeviceParameters p = base;
    p.mesh_level = lvl;
    MeniscusMesh mesh = build_meniscus_mesh(p, surface);
    DielectricSetup s;
    s.geometry = p;
    s.materials = materials;
    s.conductor_model = ConductorModel::Dielectric;
    s.metallisation = metallisation;
    s.far_field = far_field;
    s.V_emitter = V_emitter;
    s.V_extractor = V_extractor;
    s.memory_cap_bytes = memory_cap_bytes;
    DielectricSolution sol =
        solve_dielectric_on(mesh.device, s, DielectricDiagnostics::FieldOnly);
    MaxwellLoad L = maxwell_load(mesh, sol, gamma_over_a);

    EdgeStudyPoint pt;
    pt.mesh_level = lvl;
    pt.n_nodes = sol.fem.n_nodes;
    pt.n_surface_segments = static_cast<Index>(L.seg_force.size());
    pt.smallest_d = L.seg_d_mid.empty() ? 0.0 : L.seg_d_mid.back();
    pt.total_force = L.total_force;
    pt.force_coarse = L.force_beyond(edge_gate::kExclusionCoarse * a);
    pt.force_mid = L.force_beyond(edge_gate::kExclusionMid * a);
    pt.force_fine = L.force_beyond(edge_gate::kExclusionFine * a);
    for (Real v : L.node_pM) pt.peak_node_pM = std::max(pt.peak_node_pM, v);
    // Only away from the contact line.  AT the edge the recovery necessarily
    // averages vacuum cells that sit above the emitter land as well, where the
    // surface normal of the meniscus is not the normal of the boundary, so a
    // large tangential fraction there is the edge, not an error.
    for (std::size_t k = 0; k < L.node_tangential_fraction.size(); ++k)
      if (L.node_d_edge[k] >= edge_gate::kExclusionCoarse * a)
        pt.max_tangential_fraction =
            std::max(pt.max_tangential_fraction, L.node_tangential_fraction[k]);
    pt.fit_d_lo = std::max(2.0 * pt.smallest_d, 5.0e-3 * a);
    pt.fit_d_hi = 0.3 * a;
    fit_exponent(L, pt.fit_d_lo, pt.fit_d_hi, &pt.fit_exponent, &pt.fit_r2, &pt.fit_points);
    gate.levels.push_back(pt);
    gate.loads.push_back(std::move(L));
  }

  if (gate.levels.size() < 2) {
    gate.verdict = GateVerdict::FailedNoData;
    gate.note = explain(gate.verdict);
    return gate;
  }

  const EdgeStudyPoint& fine = gate.levels.back();
  const EdgeStudyPoint& coarse = gate.levels[gate.levels.size() - 2];
  const MaxwellLoad& Lf = gate.loads.back();
  const MaxwellLoad& Lc = gate.loads[gate.loads.size() - 2];

  gate.fitted_exponent = fine.fit_exponent;
  // Reference value only: for a conducting wedge with vacuum opening angle
  // alpha the field goes as d^(pi/alpha - 1), so the pressure goes as twice
  // that.  It ignores the dielectric on the other side of the edge and is
  // therefore NOT a prediction of the number that was fitted.
  {
    const Real psi_c = surface.flat ? 0.0 : surface.psi.back();
    const Real alpha = pi - psi_c;
    gate.wedge_reference_exponent = (alpha > 0.0) ? 2.0 * (pi / alpha - 1.0) : 0.0;
  }

  gate.measured_total_force_change =
      (std::abs(fine.total_force) > 0.0)
          ? std::abs(fine.total_force - coarse.total_force) / std::abs(fine.total_force)
          : 0.0;
  {
    Real worst = 0.0;
    for (int k = 0; k <= 16; ++k) {
      const Real frac = static_cast<Real>(k) / 16.0;
      const Real d = (edge_gate::kExclusionCoarse +
                      (0.9 - edge_gate::kExclusionCoarse) * frac) * a;
      const Real pf = Lf.pressure_at_distance(d);
      const Real pc = Lc.pressure_at_distance(d);
      if (pf > 0.0) worst = std::max(worst, std::abs(pf - pc) / pf);
    }
    gate.measured_edge_far_change = worst;
  }
  gate.measured_exclusion_change =
      (std::abs(fine.force_fine) > 0.0)
          ? std::abs(fine.force_fine - fine.force_mid) / std::abs(fine.force_fine)
          : 0.0;

  // --- two extrapolations of the same limit force ---------------------------
  //
  // Over the MESH: Aitken on the total force of the three finest levels.  The
  // total force approaches its limit from below as h^(1+beta), so the sequence
  // is monotone and the extrapolation is the honest way to name the limit the
  // finite meshes are heading for.
  if (gate.levels.size() >= 3) {
    const Real f1 = gate.levels[gate.levels.size() - 3].total_force;
    const Real f2 = coarse.total_force;
    const Real f3 = fine.total_force;
    const Real d1 = f2 - f1, d2 = f3 - f2;
    gate.limit_force_mesh = (std::abs(d2 - d1) > 0.0) ? f3 - d2 * d2 / (d2 - d1) : f3;
  } else {
    gate.limit_force_mesh = fine.total_force;
  }
  // Over the EXCLUSION DISTANCE: least squares of F(d0) against d0^(1+beta)
  // with the independently fitted beta, on the finest mesh.
  {
    const Real q = 1.0 + gate.fitted_exponent;
    const Real d[3] = {edge_gate::kExclusionFine * a, edge_gate::kExclusionMid * a,
                       edge_gate::kExclusionCoarse * a};
    const Real F[3] = {fine.force_fine, fine.force_mid, fine.force_coarse};
    Real sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int k = 0; k < 3; ++k) {
      const Real x = std::pow(d[k], q);
      sx += x;
      sy += F[k];
      sxx += x * x;
      sxy += x * F[k];
    }
    const Real n = 3.0;
    const Real den = n * sxx - sx * sx;
    gate.limit_force_exclusion =
        (std::abs(den) > 0.0) ? (sy * sxx - sx * sxy) / den : fine.total_force;
  }
  gate.measured_limit_agreement =
      (std::abs(gate.limit_force_exclusion) > 0.0)
          ? std::abs(gate.limit_force_mesh - gate.limit_force_exclusion) /
                std::abs(gate.limit_force_exclusion)
          : 1.0;

  if (!(gate.fitted_exponent > edge_gate::kMinExponent))
    gate.verdict = GateVerdict::FailedNotIntegrable;
  else if (gate.measured_edge_far_change > edge_gate::kTolEdgeFarLoad)
    gate.verdict = GateVerdict::FailedEdgeFarLoad;
  else if (gate.measured_total_force_change > edge_gate::kTolTotalForce)
    gate.verdict = GateVerdict::FailedTotalForce;
  else if (gate.measured_limit_agreement > edge_gate::kTolLimitAgreement)
    gate.verdict = GateVerdict::FailedExclusion;
  else
    gate.verdict = GateVerdict::Passed;
  gate.note = explain(gate.verdict);
  return gate;
}

// ===========================================================================
// The load actually used by the coupling
// ===========================================================================
//
// The type that carries it -- ProjectedLoad -- used to live here in an
// anonymous namespace, which is exactly why its central claim (the load handed
// to the capillary solver is CONTINUOUS and carries the integrated Maxwell
// force) could not be checked from a test.  P0 moved it to
// include/es/load_projection.hpp unchanged in substance and put an audit and a
// test around it.  What is left here is the use of it.

namespace {

DielectricSetup setup_from(const CoupledRequest& q) {
  DielectricSetup s;
  s.geometry = q.geometry;
  s.materials = q.materials;
  s.conductor_model = ConductorModel::Dielectric;
  s.metallisation = q.metallisation;
  s.far_field = q.far_field;
  s.V_emitter = q.V_emitter;
  s.V_extractor = q.V_extractor;
  s.memory_cap_bytes = q.memory_cap_bytes;
  return s;
}

}  // namespace

// ===========================================================================
// The coupled fixed point
// ===========================================================================

CoupledPoint solve_coupled(const CoupledRequest& q) {
  CoupledPoint out;
  out.V_emitter = q.V_emitter;
  out.V_extractor = q.V_extractor;
  out.delta_p_exit = q.delta_p_exit;
  out.mesh_level = q.geometry.mesh_level;

  const Real a = 0.5 * q.geometry.device.phi_2;
  const Real gamma = q.liquid.gamma;
  const Real gamma_over_a = gamma / a;

  // --- the field-free member of this branch ---------------------------------
  CapillaryMeniscus shape;
  if (q.initial_shape != nullptr && is_usable(q.initial_shape->status)) {
    shape = *q.initial_shape;
  } else {
    CapillaryRequest cr;
    cr.delta_p_exit = q.delta_p_exit;
    cr.target_relative_accuracy = 1.0e-10;
    shape = solve_capillary_meniscus(a, 0.0, q.liquid, cr);
    if (!is_usable(shape.status)) {
      out.status = CouplingStatus::CapillaryRangeExceeded;
      out.message = std::string("Startform: ") + explain(shape.status);
      return out;
    }
  }

  ProjectedLoad load;             // starts at zero: the field-free problem
  MaxwellLoad measured;
  Real shape_change = 0.0, load_change = 0.0;

  const DielectricSetup setup = setup_from(q);

  for (int it = 1; it <= coupling::kMaxIterations; ++it) {
    out.iterations = it;

    // 1. conformal volume mesh for the current shape
    MeniscusMesh mesh;
    try {
      mesh = build_meniscus_mesh(q.geometry, FreeSurface::from(shape));
    } catch (const std::exception& e) {
      out.status = CouplingStatus::MeshInvalid;
      out.message = e.what();
      return out;
    }
    out.min_jacobian = mesh.quality.min_jacobian;

    // 2. dielectric electrostatics on it
    DielectricSolution sol;
    try {
      sol = solve_dielectric_on(mesh.device, setup, DielectricDiagnostics::FieldOnly);
    } catch (const std::exception& e) {
      out.status = CouplingStatus::ElectrostaticFailure;
      out.message = e.what();
      return out;
    }
    out.fem_residual = sol.fem.residual_inf;
    out.fem_residual_relative =
        std::abs(sol.Q_emitter) > 0.0 ? sol.fem.residual_inf / std::abs(sol.Q_emitter) : 0.0;

    // 3./4./5. one-sided normal field, Maxwell pressure, conservative projection
    measured = maxwell_load(mesh, sol, gamma_over_a);
    const ProjectedLoad fresh = ProjectedLoad::from(measured);
    load_change = ProjectedLoad::difference(load, fresh);
    load = (it == 1) ? fresh : ProjectedLoad::blend(load, fresh, q.relaxation);

    // 6. Young-Laplace with the space-dependent load
    CapillaryRequest cr;
    cr.delta_p_exit = q.delta_p_exit;
    // The load is known to the surface discretisation of the volume mesh, a few
    // tens of segments, and its reconstruction is continuous but not smooth.
    // Asking the integrator for 1e-10 against it would be asking for an
    // accuracy the right-hand side does not have; 1e-8 is still two decades
    // below the shape tolerance of the fixed point.
    cr.target_relative_accuracy = 1.0e-8;
    if (!load.is_zero())
      cr.extra_normal_load = [&load](Real tau) { return load.at(tau); };
    else
      cr.target_relative_accuracy = 1.0e-10;   // the P3a problem, solved as it
    const CapillaryMeniscus fresh_shape = solve_capillary_meniscus(a, 0.0, q.liquid, cr);
    if (!is_usable(fresh_shape.status)) {
      out.status = (fresh_shape.status == CapillaryStatus::PressureOutsideCapillaryRange)
                       ? CouplingStatus::CapillaryRangeExceeded
                       : CouplingStatus::MechanicalResidualNotConverged;
      out.message = explain(fresh_shape.status);
      out.load = measured;
      return out;
    }
    if (fresh_shape.crossings > 1) {
      out.status = CouplingStatus::MultipleSolutionsDetected;
      out.message = std::string(explain(out.status)) + "  Der Schiessvorgang trifft den "
                    "Kontaktradius " + std::to_string(fresh_shape.crossings) + " mal.";
      out.crossings = fresh_shape.crossings;
      out.load = measured;
      return out;
    }

    // shape change, on the common parametrisation
    shape_change = 0.0;
    {
      const std::size_t n = std::min(shape.nodes.size(), fresh_shape.nodes.size());
      for (std::size_t k = 0; k < n; ++k) {
        const Real t = static_cast<Real>(k) / static_cast<Real>(n - 1);
        const std::size_t ka =
            static_cast<std::size_t>(t * static_cast<Real>(shape.nodes.size() - 1) + 0.5);
        const std::size_t kb =
            static_cast<std::size_t>(t * static_cast<Real>(fresh_shape.nodes.size() - 1) + 0.5);
        shape_change = std::max(shape_change, norm(fresh_shape.nodes[kb] - shape.nodes[ka]) / a);
      }
    }
    shape = fresh_shape;
    out.crossings = fresh_shape.crossings;

    if (q.verbose)
      std::printf("    Iteration %2d: dForm/a = %.3e, dp_M/(gamma/a) = %.3e, h/a = %+.6f\n", it,
                  shape_change, load_change / gamma_over_a, shape.apex_height / a);

    if (shape_change <= coupling::kTolShape &&
        load_change / gamma_over_a <= coupling::kTolLoad && it >= 2)
      break;
    if (it == coupling::kMaxIterations) {
      out.status = CouplingStatus::MechanicalResidualNotConverged;
      out.message = "Die Fixpunktiteration hat nach " +
                    std::to_string(coupling::kMaxIterations) +
                    " Schritten die Grenzen nicht erreicht: dForm/a = " +
                    std::to_string(shape_change) + ", dp_M/(gamma/a) = " +
                    std::to_string(load_change / gamma_over_a) + ".";
      out.shape = shape;
      out.load = measured;
      out.final_shape_change = shape_change;
      out.final_load_change = load_change / gamma_over_a;
      return out;
    }
  }

  // --- what came out ---------------------------------------------------------
  out.shape = shape;
  out.load = measured;
  out.final_shape_change = shape_change;
  out.final_load_change = load_change / gamma_over_a;
  out.apex_height = shape.apex_height;
  out.contact_error = std::abs(shape.contact().r - a) / a;
  out.total_force = measured.total_force;
  out.liquid_volume = shape.revolved_volume;
  out.surface_area = shape.revolved_area;

  const ResidualProfile res = young_laplace_residual(shape);
  out.mechanical_residual = res.max_abs;
  for (std::size_t k = 0; k < res.residual.size(); ++k) {
    const Real kappa = res.residual[k] / a +
                       (q.delta_p_exit + (shape.load.empty() ? 0.0 : shape.load[k])) / gamma;
    out.max_curvature = std::max(out.max_curvature, std::abs(kappa));
    // Away from the contact line, where the load is a converged quantity.
    if (shape.arclength - res.s[k] >= edge_gate::kExclusionCoarse * a)
      out.mechanical_residual_edge_far =
          std::max(out.mechanical_residual_edge_far, std::abs(res.residual[k]));
  }

  // Field quantities the report quotes: the apex value, and one at a fixed
  // distance from the edge that is deliberately OUTSIDE the singular zone.
  if (!measured.node_En.empty()) {
    out.E_apex = measured.node_En.front();
    Real best = std::numeric_limits<Real>::max();
    for (std::size_t k = 0; k < measured.node_d_edge.size(); ++k) {
      const Real dd = std::abs(measured.node_d_edge[k] - 0.25 * a);
      if (dd < best) {
        best = dd;
        out.E_edge_far = measured.node_En[k];
      }
    }
  }

  if (out.contact_error > coupling::kTolContact ||
      out.mechanical_residual_edge_far > coupling::kTolMechanical) {
    out.status = CouplingStatus::MechanicalResidualNotConverged;
    out.message = "Kontaktlinie oder mechanisches Residuum ausserhalb der Grenzen: "
                  "|dr|/a = " + std::to_string(out.contact_error) +
                  ", Residuum kantenfern = " + std::to_string(out.mechanical_residual_edge_far) +
                  " (ueber die ganze Oberflaeche " + std::to_string(out.mechanical_residual) +
                  ").";
    return out;
  }
  out.status = CouplingStatus::Converged;
  return out;
}

// ===========================================================================
// Continuation
// ===========================================================================

ContinuationResult continue_over_voltage(CoupledRequest q, Real V_max) {
  ContinuationResult out;
  const Real sign = (V_max >= 0.0) ? 1.0 : -1.0;
  const Real target = std::abs(V_max);

  // The field-free member of the branch, computed once and used as the start.
  CoupledRequest zero = q;
  zero.V_emitter = 0.0;
  zero.V_extractor = 0.0;
  CoupledPoint p0 = solve_coupled(zero);
  if (!is_usable(p0.status)) {
    out.end_status = p0.status;
    out.end_message = "Der feldfreie Ast konnte nicht gerechnet werden: " + p0.message;
    return out;
  }
  out.points.push_back(p0);
  out.last_converged_voltage = 0.0;

  Real V = 0.0;
  Real step = continuation::kFirstStep;
  CapillaryMeniscus last_shape = p0.shape;

  while (V < target) {
    const Real V_try = std::min(target, V + step);
    ++out.steps_attempted;
    CoupledRequest qi = q;
    qi.V_emitter = sign * V_try;
    qi.V_extractor = 0.0;
    qi.initial_shape = &last_shape;
    CoupledPoint p = solve_coupled(qi);
    if (is_usable(p.status)) {
      out.points.push_back(p);
      last_shape = p.shape;
      V = V_try;
      out.last_converged_voltage = sign * V;
      step = std::min(step * continuation::kGrowth, continuation::kFirstStep);
      continue;
    }
    ++out.steps_rejected;
    out.first_failed_voltage = sign * V_try;
    if (step / 2.0 < continuation::kMinStep) {
      out.end_status = (p.status == CouplingStatus::CapillaryRangeExceeded)
                           ? CouplingStatus::CapillaryRangeExceeded
                           : CouplingStatus::ContinuationStepTooSmall;
      out.end_message = "Bei " + std::to_string(sign * V_try) + " V: " +
                        std::string(to_string(p.status)) + " -- " + p.message +
                        "  Die Schrittweite ist unter " + std::to_string(continuation::kMinStep) +
                        " V gefallen; der Ast wird hier nicht weiter verfolgt. Das ist die "
                        "Stelle, an der dieser Loeser stehen bleibt, und keine Aussage ueber "
                        "Stabilitaet, Emissionsbeginn oder einen Taylor-Kegel.";
      return out;
    }
    step /= 2.0;
  }
  out.end_status = CouplingStatus::Converged;
  out.end_message = "Die Fortsetzung hat die vorgegebene Hoechstspannung erreicht.";
  return out;
}

}  // namespace es
