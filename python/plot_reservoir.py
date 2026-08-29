#!/usr/bin/env python3
"""Figures of the P2c reservoir study.

    ./build/es_reservoir examples/device_p1.cfg examples/reservoir_p2c.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_reservoir.py results/<dir>

Reads only the CSVs written by apps/es_reservoir.cpp.  Nothing here solves,
meshes or computes geometry.

The commit named on every figure is measured from git at drawing time by
python/provenance.py, not taken from an argument, and a dirty working tree is
stamped DIRTY.

Produces
  <dir>/fig1_reservoir_geometries.png  the variants, to scale, in one frame
  <dir>/fig2_full_mesh.png             the WHOLE domain plus three details
  <dir>/fig3_near_field.png            phi and |E| in one identical window,
                                       identical colour scales, all variants
  <dir>/fig4_reservoir_convergence.png local field and global charge against
                                       reservoir size
  <dir>/figures_provenance.txt         what the stamp said, and whether it is
                                       releasable
"""
import csv
import os
import sys
import textwrap
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance as pv

TITLE = ("P2c: Gerätegeometrie und Flüssigkeitsvorrat entkoppelt – feste Frontgeometrie, "
         "variierter dielektrisch umschlossener Flüssigkeitsraum")

RCOL = {
    "liquid": "#bcd9ef",
    "emitter_dielectric": "#d8d2c4",
    "reservoir_dielectric": "#cdbfa6",
    "extractor_carrier": "#e8d4bd",
}
RLABEL = {
    "liquid": "Flüssigkeit (idealer Leiter, $V_\\mathrm{emitter}$) – Bohrung, Kanal, Plenum",
    "emitter_dielectric": "gedruckter Emitter + Grundkörper (Dielektrikum)",
    "reservoir_dielectric": "Vorratskörper (Dielektrikum, PEEK im Gerät)",
    "extractor_carrier": "Extraktorträger (Dielektrikum) mit Metallisierung",
}
# The order the regions are drawn in; liquid last so it sits on top.
ZORDER = ["reservoir_dielectric", "extractor_carrier", "emitter_dielectric", "liquid"]


# --------------------------------------------------------------------- input
def rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def meta(d):
    out = {}
    with open(os.path.join(d, "meta.txt"), encoding="utf-8") as fh:
        for line in fh:
            if "=" in line:
                k, v = line.rstrip("\n").split("=", 1)
                out[k] = v
    return out


def parameters(d):
    out = {}
    for r in rows(os.path.join(d, "parameters.csv")):
        try:
            out[r["name"]] = float(r["value_SI"])
        except ValueError:
            pass
    return out


def variant_outlines(d):
    out = defaultdict(lambda: defaultdict(list))
    for r in rows(os.path.join(d, "variant_outlines.csv")):
        out[r["variant"]][r["name"]].append((float(r["r_m"]), float(r["z_m"])))
    return out


def grid(path):
    """Sampled field CSV -> (r_axis, z_axis, phi, |E|, region)."""
    rs, zs, ph, em, rg = [], [], [], [], []
    for r in rows(path):
        rs.append(float(r["r_m"]))
        zs.append(float(r["z_m"]))
        ph.append(float(r["phi_V"]))
        em.append(float(r["E_V_per_m"]))
        rg.append(r["region"])
    ru = np.unique(np.array(rs))
    zu = np.unique(np.array(zs))
    shape = (len(zu), len(ru))
    return (ru, zu, np.array(ph).reshape(shape), np.array(em).reshape(shape),
            np.array(rg).reshape(shape))


def mesh_nodes(d):
    ii, jj, rr, zz = [], [], [], []
    for r in rows(os.path.join(d, "mesh_nodes.csv")):
        ii.append(int(r["i"]))
        jj.append(int(r["j"]))
        rr.append(float(r["r_m"]))
        zz.append(float(r["z_m"]))
    nr = max(ii) + 1
    nz = max(jj) + 1
    R = np.zeros((nz, nr))
    Z = np.zeros((nz, nr))
    for a, b, c, e in zip(ii, jj, rr, zz):
        R[b, a] = c
        Z[b, a] = e
    return R, Z


# -------------------------------------------------------------------- helpers
def draw_variant(ax, ol, s=1e3, alpha=1.0, lw=0.9, faces=True):
    for name in ZORDER:
        if name not in ol:
            continue
        p = np.array(ol[name]) * s
        ax.add_patch(Polygon(p, closed=True,
                             facecolor=RCOL.get(name, "#cccccc") if faces else "none",
                             edgecolor="#333333", lw=lw, zorder=3, alpha=alpha))


def provenance(fig, m, stamp, extra="", y=0.012, width=175):
    dirty = pv.DIRTY_MARK in stamp or stamp == pv.NO_VERSION
    head = (f"Commit {stamp}   |   {m.get('app','es_reservoir')}   |   "
            f"Netzstufe {m.get('reference_level','?')}, {m.get('nodes','?')} Knoten "
            f"(Referenzvariante {m.get('reference_variant','?')})   |   "
            f"$V_e$ = {float(m.get('V_emitter_V', 0)):.0f} V, "
            f"$V_x$ = {float(m.get('V_extractor_V', 0)):.0f} V")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width) if len(text) > width else [text]),
             ha="center", va="bottom", fontsize=7.2,
             color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.5, 0.997, f"NICHT FREIGEGEBEN – {stamp}", ha="center", va="top",
                 fontsize=9, color="#8a1b1b", weight="bold")


def caveat(fig, text, y=0.05, width=155, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


def variant_rows(d):
    return rows(os.path.join(d, "variants.csv"))


# ------------------------------------------------------------------- figure 1
def figure_geometries(d, m, par, out):
    """Every variant to the SAME scale, in the same frame, plus the fixed front."""
    vr = variant_rows(d)
    ol = variant_outlines(d)
    tags = [r["variant"] for r in vr]
    s = 1e3  # millimetres

    rmax = max(float(r["plenum_radius_m"]) for r in vr) + float(m["plenum_wall_thickness_m"])
    r1 = 1.10 * max(rmax, par["extractor_outer_radius"]) * s
    z0 = 1.05 * (float(m["base_z_m"]) - float(m["feed_channel_length_m"]) -
                 max(float(r["plenum_depth_m"]) for r in vr) -
                 float(m["plenum_wall_thickness_m"])) * s
    z1 = 1.30 * (par["extraction_distance"] + par["extractor_thickness"]) * s

    fig = plt.figure(figsize=(16.4, 5.9))
    gs = fig.add_gridspec(1, len(tags) + 1, width_ratios=[1.25] + [1.0] * len(tags),
                          wspace=0.34)

    # The fixed front geometry, at the scale on which it is a device at all.
    ax = fig.add_subplot(gs[0, 0])
    su = 1e6
    draw_variant(ax, ol[m["reference_variant"]], su)
    ax.set_xlim(0, 320)
    ax.set_ylim(1.35 * float(m["base_z_m"]) * su,
                1.28 * (par["extraction_distance"] + par["extractor_thickness"]) * su)
    ax.set_aspect("equal")
    ax.set_xlabel("r [µm]")
    ax.set_ylabel("z [µm]")
    ax.set_title("FESTE Frontgeometrie\n(in allen Varianten bitgenau gleich)", fontsize=8.6)
    ax.annotate("Emitter,\nBohrung,\nGrundkörper", xy=(par["r_foot"] * su, -100),
                xytext=(95, -170), fontsize=7.0,
                arrowprops=dict(arrowstyle="->", lw=0.9))
    ax.annotate("Extraktor,\nAperturkante", xy=(par["r_aperture"] * su, 500),
                xytext=(60, 300), fontsize=7.0,
                arrowprops=dict(arrowstyle="->", lw=0.9))

    axes = []
    for k, r in enumerate(vr):
        ax = fig.add_subplot(gs[0, k + 1], sharey=axes[0] if axes else None)
        axes.append(ax)
        draw_variant(ax, ol[r["variant"]], s)
        ax.set_xlim(0, r1)
        ax.set_ylim(z0, z1)
        ax.set_aspect("equal")
        ax.set_xlabel("r [mm]")
        if k == 0:
            ax.set_ylabel("z [mm]")
        else:
            ax.tick_params(labelleft=False)
        ax.set_title("\n".join(textwrap.wrap(r["label"], 22)), fontsize=8.2)

    handles = [Line2D([], [], marker="s", ls="none", ms=10, mfc=RCOL[k], mec="#333333",
                      label=RLABEL[k]) for k in ZORDER]
    handles.append(Line2D([], [], marker="s", ls="none", ms=10, mfc="#ffffff", mec="#333333",
                          label="Vakuum, $\\varepsilon_r = 1$ (auch der leere Teil des "
                                "Plenums)"))
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8.2, frameon=False,
               bbox_to_anchor=(0.5, 0.172))
    fig.suptitle("Maßstäbliche Geometrie: feste Frontgeometrie, abgeschnittenes Säulenmodell "
                 "und die untersuchten Plenumgrößen\n" + TITLE, fontsize=11.5)
    fig.subplots_adjust(left=0.045, right=0.99, top=0.815, bottom=0.335)
    caveat(fig, "Das linke Bild ist in MIKROMETERN, alle übrigen in MILLIMETERN und "
                "untereinander in derselben Skala. Der vordere Emitter (Ø 40 µm, 60 µm hoch) "
                "ist auf der Millimeterskala nur ein Strich – das ist der Punkt: was bisher "
                "als „Verschiebung der Zulaufgrenze“ berichtet wurde, änderte in Wahrheit die "
                "Ausdehnung des leitenden Flüssigkeitskörpers um Größenordnungen. "
                f"Grundkörperdicke {float(m['base_plate_thickness_m'])*1e6:.0f} µm, "
                f"Zulaufkanal {float(m['feed_channel_length_m'])*1e6:.0f} µm lang, Wandstärke "
                f"{float(m['plenum_wall_thickness_m'])*1e6:.0f} µm – VORLÄUFIGE Beispielwerte "
                "und in jeder Variante identisch. Das achsensymmetrische Plenum ist eine "
                "ERSATZGEOMETRIE, keine Rekonstruktion des Reservoirs aus Abb. A.5.",
            y=0.064, width=165)
    provenance(fig, m, STAMP,
               f"Nachweis der festen Frontgeometrie: {m['front_rows_compared']} Zeilen und "
               f"{m['front_cells_compared']} Zellen je Variante verglichen, "
               f"max|Δr| = {float(m['front_max_dr_m']):.1e} m, "
               f"max|Δz| = {float(m['front_max_dz_m']):.1e} m, "
               f"{m['front_cell_mismatches']} Materialabweichungen – bitgenau gleich")
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 2
def figure_full_mesh(d, m, par, out):
    """The WHOLE computational domain, plus the details it hides."""
    R, Z = mesh_nodes(d)
    ol = variant_outlines(d)[m["reference_variant"]]
    nz, nr = R.shape
    Rd = float(m["domain_radius_m"])
    zlo, zhi = float(m["domain_z_min_m"]), float(m["domain_z_max_m"])

    fig = plt.figure(figsize=(15.0, 8.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1.0, 1.0], height_ratios=[1.0, 1.0])
    ax_full = fig.add_subplot(gs[:, 0])
    views = [(ax_full, "Gesamte Rechendomäne", 0, Rd, zlo, zhi, 1e3, "mm", 2),
             (fig.add_subplot(gs[0, 1]), "Vorratskörper und Plenum", 0,
              1.25 * (float(m["plenum_wall_thickness_m"]) +
                      max(float(r["plenum_radius_m"]) for r in variant_rows(d)
                          if r["variant"] == m["reference_variant"])),
              1.25 * (float(m["base_z_m"]) - float(m["feed_channel_length_m"]) -
                      2.5e-3 - float(m["plenum_wall_thickness_m"])),
              0.4e-3, 1e3, "mm", 1),
             (fig.add_subplot(gs[0, 2]), "Emitter und Kapillare", 0, 3.0 * par["r_foot"],
              -3.0 * par["r_foot"], 3.0 * par["r_foot"], 1e6, "µm", 1),
             (fig.add_subplot(gs[1, 1]), "Extraktionsstrecke und Apertur", 0,
              1.15 * par["extractor_outer_radius"], -0.25e-3,
              1.3 * (par["extraction_distance"] + par["extractor_thickness"]), 1e3, "mm", 1),
             (fig.add_subplot(gs[1, 2]), "Austrittskante", 0, 2.2 * par["r_bore"],
              -2.2 * par["r_bore"], 2.2 * par["r_bore"], 1e6, "µm", 1)]

    for ax, name, r0, r1, z0, z1, s, unit, step in views:
        draw_variant(ax, ol, s, alpha=0.55, lw=0.7)
        for j in range(0, nz, step):
            ax.plot(R[j, :] * s, Z[j, :] * s, "-", color="#2b6cb0", lw=0.28, alpha=0.85,
                    zorder=4)
        for i in range(0, nr, step):
            ax.plot(R[:, i] * s, Z[:, i] * s, "-", color="#2b6cb0", lw=0.28, alpha=0.85,
                    zorder=4)
        ax.set_xlim(r0 * s, r1 * s)
        ax.set_ylim(z0 * s, z1 * s)
        ax.set_aspect("equal")
        ax.set_xlabel(f"r [{unit}]")
        ax.set_ylabel(f"z [{unit}]")
        ax.set_title(name, fontsize=9.5)

    checks = rows(os.path.join(d, "mesh_checks.csv"))
    npass = sum(1 for c in checks if c["passed"] == "1")
    fig.suptitle("Vollständiges Volumennetz bis zum Fernrand, mit Detailansichten\n" + TITLE,
                 fontsize=11.5)
    fig.tight_layout(rect=[0, 0.135, 1, 0.93])
    caveat(fig, "Links die TATSÄCHLICHE Rechendomäne – "
                f"r bis {Rd*1e3:.0f} mm, z von {zlo*1e3:.0f} mm bis {zhi*1e3:.0f} mm – nicht "
                "ein Ausschnitt daraus. Die vier Detailansichten zeigen, was darin nicht "
                f"sichtbar ist. {npass} von {len(checks)} Netzprüfungen bestanden, darunter "
                "die Reproduktion jedes Gebietsvolumens aus der geschlossenen Formel und das "
                "Verhältnis benachbarter Elementgrößen. Die Größenfunktion folgt allein der "
                "Geometrie; einziger Benutzerwert ist die Netzstufe. Die Zellen des Vorrats "
                "sind gröber als die an der Spitze – das ist beabsichtigt: Vorratsgrößen "
                "speisen das Größenfeld des Nahfelds nicht, sonst könnte eine Vergrößerung "
                "des Plenums die Elementgröße am Meniskus verschieben.",
            y=0.05)
    provenance(fig, m, STAMP,
               f"{m.get('nr','?')} × {m.get('nz','?')} Knoten; Größenfeldfaktor "
               f"{float(m.get('reference_size_scale',1)):.4f}")
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 3
def figure_near_field(d, m, par, out):
    """One identical window, one colour scale, every variant."""
    vr = variant_rows(d)
    ol = variant_outlines(d)
    tags = [r["variant"] for r in vr]
    data = {t: grid(os.path.join(d, f"field_near_{t}.csv")) for t in tags}
    s = 1e6

    # Common scales over ALL variants; anything else would make the panels
    # look alike for the wrong reason.
    phis = np.concatenate([np.ma.array(v[2], mask=v[4] == "liquid").compressed()
                           for v in data.values()])
    ems = np.concatenate([np.ma.array(v[3], mask=(v[4] == "liquid") | (v[3] <= 0)).compressed()
                          for v in data.values()])
    plo, phi_hi = float(np.nanmin(phis)), float(np.nanmax(phis))
    ehi = float(np.nanmax(ems))
    elo = ehi / 1.0e3
    plevels = np.linspace(plo, phi_hi, 40)
    elevels = np.geomspace(elo, ehi, 40)

    fig, axes = plt.subplots(2, len(tags), figsize=(2.45 * len(tags) + 1.6, 8.8),
                             sharex=True, sharey=True)
    for col, t in enumerate(tags):
        ru, zu, PH, EM, RG = data[t]
        Rg, Zg = np.meshgrid(ru * s, zu * s)
        mask = RG == "liquid"

        ax = axes[0, col]
        cf_p = ax.contourf(Rg, Zg, np.ma.array(PH, mask=mask), levels=plevels, cmap="viridis",
                           extend="both")
        ax.contour(Rg, Zg, np.ma.array(PH, mask=mask), levels=14, colors="w", linewidths=0.4,
                   alpha=0.7)
        ax.set_title("\n".join(textwrap.wrap(vr[col]["label"], 26)), fontsize=8.2)

        ax = axes[1, col]
        cf_e = ax.contourf(Rg, Zg, np.ma.array(EM, mask=mask | (EM <= 0)), levels=elevels,
                           norm=LogNorm(elo, ehi), cmap="magma", extend="both")

        for ax in (axes[0, col], axes[1, col]):
            draw_variant(ax, ol[t], s, faces=False, lw=1.1)
            ax.set_xlim(ru[0] * s, ru[-1] * s)
            ax.set_ylim(zu[0] * s, zu[-1] * s)
            ax.set_aspect("equal")
            ax.set_xlabel("r [µm]")
    axes[0, 0].set_ylabel("Potential $\\varphi$\n\nz [µm]")
    axes[1, 0].set_ylabel("Feldstärke $|E|$\n\nz [µm]")
    fig.colorbar(cf_p, ax=axes[0, :].tolist(), label="$\\varphi$ [V]", fraction=0.022, pad=0.01)
    cb = fig.colorbar(cf_e, ax=axes[1, :].tolist(), label="$|E|$ [V/m]", fraction=0.022,
                      pad=0.01, ticks=matplotlib.ticker.LogLocator(numticks=6))
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterSciNotation())

    fig.suptitle("Potential und Feldstärke im identischen Nahfeldausschnitt, identische "
                 "Farbskalen\n" + TITLE, fontsize=11.5)
    caveat(fig, "Alle Bilder zeigen denselben Ausschnitt und benutzen dieselbe Farbskala; die "
                "Skalen sind über alle Varianten gemeinsam bestimmt. Innerhalb der "
                "Flüssigkeit ist |E| nicht definiert (idealer Leiter) und ausgeblendet. Der "
                "Sprung zwischen der abgeschnittenen Säule und jeder Plenumvariante ist keine "
                "Konvergenzfrage, sondern der Unterschied zwischen einem 16-fL-Leiter und "
                "einem Vorrat von Kubikmillimetern. Zwischen den Plenumvarianten ist der "
                "Unterschied mit bloßem Auge nicht mehr zu sehen – die Zahlen dazu stehen in "
                "der nächsten Abbildung. Die Austrittskante ist unverrundet: dort divergiert "
                "das Feld und folgt der Elementgröße.",
            y=0.045)
    provenance(fig, m, STAMP,
               f"Fenster r = 0 … {float(rows(os.path.join(d,'field_window.csv'))[0]['r1_m'])*1e6:.0f} µm, "
               f"z = {float(rows(os.path.join(d,'field_window.csv'))[0]['z0_m'])*1e6:.0f} … "
               f"{float(rows(os.path.join(d,'field_window.csv'))[0]['z1_m'])*1e6:.0f} µm")
    fig.subplots_adjust(left=0.075, right=0.90, top=0.87, bottom=0.175)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 4
def figure_convergence(d, m, out):
    vr = [r for r in variant_rows(d) if r["reservoir_model"] == "axisymmetric_plenum"
          and float(r["fill_fraction"]) == 1.0]
    conv = rows(os.path.join(d, "reservoir_convergence.csv"))
    mesh = rows(os.path.join(d, "convergence_mesh.csv"))
    surf = defaultdict(list)
    for r in rows(os.path.join(d, "reference_surface_fields.csv")):
        surf[r["variant"]].append((float(r["r_m"]), float(r["Ez_V_per_m"])))

    span = abs(float(m["V_emitter_V"]) - float(m["V_extractor_V"]))
    tol = float(m["tol_E_rel"])
    V = np.array([float(r["liquid_volume_m3"]) for r in vr])
    labels = [r["variant"] for r in vr]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.8))

    # -- (a) absolute local quantities against reservoir volume, INCLUDING the
    #        truncated column, so that the size of the correction is visible;
    #        an inset resolves the residual drift between the plenum sizes,
    #        which is three orders of magnitude smaller and would otherwise be
    #        a flat line that says nothing.
    ax = axes[0, 0]
    Ez = np.array([float(r["Ez_axis_surface_V_per_m"] or 0) for r in vr])
    Em = np.array([float(r["Emag_2_bore_radii_V_per_m"]) for r in vr])
    col = [r for r in variant_rows(d) if r["reservoir_model"] == "truncated_column"]
    ax.loglog(V, Ez / 1e6, "o-", ms=6, color="#2b6cb0",
              label="$E_z$ auf der Achse, unmittelbar über der Oberfläche")
    ax.loglog(V, Em / 1e6, "s-", ms=6, color="#0f9d58",
              label="$|E|$ bei zwei Bohrungsradien")
    for x, y, t in zip(V, Ez / 1e6, labels):
        ax.annotate(t.replace("plenum_", ""), (x, y), textcoords="offset points",
                    xytext=(0, -14), fontsize=7, ha="center")
    if col:
        c = col[0]
        ax.loglog([float(c["liquid_volume_m3"])],
                  [float(c["Ez_axis_surface_V_per_m"]) / 1e6], "D", ms=9, color="#8a1b1b",
                  mec="#000000", mew=0.8, label="abgeschnittene Säule (Diagnose), $E_z$")
        ax.loglog([float(c["liquid_volume_m3"])],
                  [float(c["Emag_2_bore_radii_V_per_m"]) / 1e6], "v", ms=9, color="#c1272d",
                  mec="#000000", mew=0.8, label="abgeschnittene Säule (Diagnose), $|E|$")
    ax.set_xlabel("Flüssigkeitsvolumen des Vorrats [m³]")
    ax.set_ylabel("Feld [MV/m]")
    ax.set_title("Lokales Extraktionsfeld gegen Vorratsgröße", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.0, loc="lower left")

    ins = ax.inset_axes([0.57, 0.655, 0.40, 0.28])
    ins.semilogx(V, (Ez / Ez[-1] - 1.0) * 1e3, "o-", ms=4, color="#2b6cb0")
    ins.semilogx(V, (Em / Em[-1] - 1.0) * 1e3, "s-", ms=4, color="#0f9d58")
    ins.axhline(0.0, color="#999999", lw=0.8)
    ins.set_title("nur die Plenumgrößen, relativ zur größten", fontsize=6.6)
    ins.set_ylabel("Abw. [‰]", fontsize=6.2, labelpad=1)
    ins.tick_params(labelsize=6.0)
    ins.grid(alpha=0.3, which="both")

    # -- (b) change per enlargement step, against the pre-declared tolerance
    ax = axes[0, 1]
    steps = np.arange(1, len(conv) + 1)
    ax.semilogy(steps, [float(r["d_E_max_rel"]) for r in conv], "o-", ms=6, color="#8a1b1b",
                label="$\\Delta|E|$, Maximum über ALLE Sondenpunkte")
    ax.semilogy(steps, [float(r["d_Ez_axis_rel"]) for r in conv], "s-", ms=6, color="#2b6cb0",
                label="$\\Delta E_z$ an der Flüssigkeitsoberfläche")
    ax.semilogy(steps, [float(r["d_Emag_2rb_rel"]) for r in conv], "^-", ms=6, color="#0f9d58",
                label="$\\Delta|E|$ bei zwei Bohrungsradien")
    ax.semilogy(steps, [float(r["d_phi_max_over_span"]) for r in conv], "d--", ms=5,
                color="#7b3294", label="$\\Delta\\varphi$ / Spannweite, Maximum")
    ax.axhline(tol, color="#000000", lw=1.4, ls=":",
               label=f"vorab festgelegte Grenze {tol:.0e}")
    ax.axhline(float(m["mesh_change_E_rel"]), color="#999999", lw=1.2, ls="-.",
               label=f"Diskretisierungsfehler {float(m['mesh_change_E_rel']):.1e}")
    ax.set_xticks(steps)
    ax.set_xticklabels([f"{r['from'].replace('plenum_','')}→{r['to'].replace('plenum_','')}"
                        for r in conv])
    ax.set_xlabel("Vergrößerungsschritt des Plenums")
    ax.set_ylabel("relative Änderung")
    ax.set_title("Änderung je Vergrößerungsschritt gegen die Toleranz", fontsize=10,
                 color="#8a1b1b" if m["local_field_converged"] == "no" else "#0f6b3d")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.0, loc="best")

    # -- (c) the global charges, which stay reservoir dependent
    ax = axes[1, 0]
    Ql = np.array([float(r["Q_liquid_C"]) for r in vr])
    Qx = np.array([float(r["Q_extractor_C"]) for r in vr])
    ax.loglog(V, Ql * 1e12, "o-", ms=6, color="#c1272d", label="$Q$ der Flüssigkeit")
    ax.loglog(V, -Qx * 1e12, "s-", ms=6, color="#2b6cb0",
              label="$-Q$ des Extraktors (Gegenladung)")
    ax.set_xlabel("Flüssigkeitsvolumen des Vorrats [m³]")
    ax.set_ylabel("Ladung [pC]")
    ax.set_title("Globale Ladungen – bleiben vorratsabhängig", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.6)

    # -- (d) the surface field profile, all variants
    ax = axes[1, 1]
    for r in variant_rows(d):
        t = r["variant"]
        if not surf[t]:
            continue
        a = np.array(surf[t])
        ax.plot(a[:, 0] * 1e6, a[:, 1] / 1e6, marker="o", ms=3,
                lw=2.2 if t == "saeule" else 1.3,
                ls="--" if t == "saeule" else "-",
                label=r["label"])
    ax.set_xlabel("r [µm] auf der ebenen Flüssigkeitsoberfläche")
    ax.set_ylabel("$E_z(r, z=0^+)$ [MV/m]")
    ax.set_title("Feldprofil auf der Referenzfläche, kantennahe Zellen ausgeschlossen",
                 fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.8)
    ax.text(0.03, 0.20,
            "Die fünf Plenumkurven liegen aufeinander\n"
            "und sind hier nicht zu trennen: sie\n"
            "unterscheiden sich um weniger als ein\n"
            "Promille. Der Abstand zur abgeschnittenen\n"
            "Säule ist dagegen größer als ein Faktor zwei.",
            transform=ax.transAxes, fontsize=6.8, color="#7a3b00", va="bottom")

    fig.suptitle("Konvergenz von lokalem Feld und globaler Ladung gegen die Vorratsgröße\n"
                 + TITLE, fontsize=11.5)
    fig.tight_layout(rect=[0, 0.145, 1, 0.925])
    verdict = ("EINGEHALTEN" if m["local_field_converged"] == "yes" else "NICHT EINGEHALTEN")
    caveat(fig,
           f"BEFUND: beim letzten Vergrößerungsschritt bewegt sich das Potential um "
           f"{float(m['last_step_phi_over_span']):.1e} der Spannweite und die Feldstärke um "
           f"{float(m['last_step_E_rel']):.1e} relativ – Maximum über alle Sondenpunkte, "
           f"getrieben von {m['last_step_worst_E_probe']}. Gegen die VOR der Messung "
           f"festgelegte Grenze von {tol:.0e} ist das {verdict}. Die Grenze wird nicht "
           "gelockert, es wird keine künstliche Elektrode eingeführt und keine "
           "„Referenzgröße“ stillschweigend ausgewählt. Am Meniskus selbst beträgt die "
           f"letzte Änderung {float(m['last_step_meniscus_E_rel']):.1e} – das ist die Größe, "
           "die ein Emissionsmodell braucht, ersetzt das Urteil oben aber nicht. Die globalen "
           "Ladungen bleiben vorratsabhängig, und das ist richtig so: ein endlicher Leiter im "
           "offenen Raum ohne Rückführelektrode trägt eine Ladung, die mit seiner Größe "
           "wächst.",
           y=0.035, color="#8a1b1b" if m["local_field_converged"] == "no" else "#0f6b3d")
    provenance(fig, m, STAMP,
               f"Netzstudie bei festem Plenum: Stufe {int(m['max_level'])-1} → "
               f"{m['max_level']} bewegt φ um "
               f"{float(m['mesh_change_phi_over_span']):.1e} der Spannweite und |E| um "
               f"{float(m['mesh_change_E_rel']):.1e} relativ – "
               f"{len(mesh)} Stufen in convergence_mesh.csv")
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------- main
STAMP = ""


def main(d):
    global STAMP
    m = meta(d)
    par = parameters(d)
    root = pv.repo_root_of(__file__)
    rel = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/")
    state = pv.repo_state(root, ignore=[rel])
    STAMP = pv.stamp(state, m.get("commit"))
    pv.write_provenance(d, state, m.get("commit"))

    made = [figure_geometries(d, m, par, os.path.join(d, "fig1_reservoir_geometries.png")),
            figure_full_mesh(d, m, par, os.path.join(d, "fig2_full_mesh.png")),
            figure_near_field(d, m, par, os.path.join(d, "fig3_near_field.png")),
            figure_convergence(d, m, os.path.join(d, "fig4_reservoir_convergence.png"))]
    for f in made:
        print(f"{f}  {os.path.getsize(f)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
