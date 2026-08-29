#!/usr/bin/env python3
"""Figures of the P2b dielectric axisymmetric electrostatics.

    ./build/es_dielectric examples/device_p1.cfg examples/dielectric_p2b.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_dielectric.py results/<dir>

Reads only the CSVs written by apps/es_dielectric.cpp.  Nothing here solves,
meshes or computes geometry.

The commit named on every figure is measured from git at drawing time by
python/provenance.py, not taken from an argument, and a dirty working tree is
stamped DIRTY.

Produces
  <dir>/fig1_material_regions.png   materials and boundary conditions
  <dir>/fig2_volume_mesh.png        automatic volume mesh, overview and detail
  <dir>/fig3_potential_field.png    potential and field, overview and tip
  <dir>/fig4_convergence.png        mesh, feed position, eps_r, far field
  <dir>/figures_provenance.txt      what the stamp said and whether it is releasable
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

TITLE = ("P2b/Diagnose: dielektrische Elektrostatik am ABGESCHNITTENEN Säulenmodell – "
         "der Rückschnitt ändert Elektrode UND Dielektrikum, siehe docs/08, 8.9")

RCOL = {
    "liquid": "#bcd9ef",
    "emitter_dielectric": "#d8d2c4",
    "extractor_carrier": "#e8d4bd",
}
RLABEL = {
    "liquid": "ionische Flüssigkeit (idealer Leiter, $V_\\mathrm{emitter}$)",
    "emitter_dielectric": "Emitterkörper SU-8 (Dielektrikum, $\\varepsilon_r$)",
    "extractor_carrier": "Extraktorträger SU-8 (Dielektrikum, $\\varepsilon_r$)",
}
ROLECOL = {
    "liquid_conductor": ("#0b6fb4", "Flüssigkeitsoberfläche: $\\varphi=V_\\mathrm{emitter}$"),
    "liquid_feed_boundary": ("#0f9d58",
                             "Schnitt durch die Flüssigkeitssäule: "
                             "$\\varphi=V_\\mathrm{emitter}$ NUR auf dem Querschnitt"),
    "extractor_metallisation": ("#c1272d",
                                "Metallisierung: $\\varphi=V_\\mathrm{extractor}$"),
    "far_field_dirichlet": ("#888888", "geerdete Hülle (nur far_field=grounded)"),
    "emitter_metal_reference": ("#7b3294", "metallischer Referenzemitter (KEIN P2b-Modell)"),
}
FREE_LABEL = ("Polymerflächen: KEINE Randbedingung – Dielektrikum, kein Leiter\n"
              "(Kegelflanke, Stirnfläche, Rückfläche, unbeschichtete Extraktorflächen)")


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
    return {r["name"]: float(r["value_SI"]) for r in rows(os.path.join(d, "parameters.csv"))}


def outline(d):
    out = defaultdict(lambda: defaultdict(list))
    for r in rows(os.path.join(d, "device_outline.csv")):
        out[r["kind"]][r["name"]].append((float(r["r_m"]), float(r["z_m"])))
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


def node_roles(d):
    out = defaultdict(list)
    for r in rows(os.path.join(d, "node_roles.csv")):
        out[r["role"]].append((float(r["r_m"]), float(r["z_m"])))
    return out


# -------------------------------------------------------------------- helpers
def draw_regions(ax, ol, s=1e6, alpha=1.0, lw=1.0):
    for name, pts in ol["region"].items():
        p = np.array(pts) * s
        ax.add_patch(Polygon(p, closed=True, facecolor=RCOL.get(name, "#cccccc"),
                             edgecolor="#333333", lw=lw, zorder=2, alpha=alpha))


def provenance(fig, m, stamp, extra="", y=0.012, width=170):
    dirty = pv.DIRTY_MARK in stamp or stamp == pv.NO_VERSION
    head = (f"Commit {stamp}   |   {m.get('app','es_dielectric')}   |   "
            f"Netzstufe {m.get('reference_level','?')}, {m.get('nodes','?')} Knoten   |   "
            f"$V_e$ = {float(m.get('V_emitter_V', 0)):.0f} V, "
            f"$V_x$ = {float(m.get('V_extractor_V', 0)):.0f} V")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width) if len(text) > width else [text]),
             ha="center", va="bottom", fontsize=7.2,
             color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.5, 0.997, f"NICHT FREIGEGEBEN – {stamp}", ha="center", va="top",
                 fontsize=9, color="#8a1b1b", weight="bold")


def caveat(fig, text, y=0.055, width=150):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color="#7a3b00")


# ------------------------------------------------------------------- figure 1
def figure_materials(d, m, par, ol, stamp, out):
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 6.6),
                             gridspec_kw={"width_ratios": [2.5, 1.0, 1.15]})
    roles = node_roles(d)
    s = 1e6
    views = [("Gesamtansicht mit Extraktionselektrode", 0,
              1.18 * par["extractor_outer_radius"] * s,
              1.7 * par["base_z"] * s,
              1.35 * (par["extraction_distance"] + par["extractor_thickness"]) * s),
             ("Emitterkörper und Säulenschnitt", 0, 2.6 * par["r_foot"] * s,
              1.08 * par["base_z"] * s, 0.35 * par["emitter_height"] * s),
             ("Austrittskante", 0, 3.2 * par["r_bore"] * s,
              -3.0 * par["r_bore"] * s, 2.4 * par["r_bore"] * s)]
    for ax, (name, r0, r1, z0, z1) in zip(axes, views):
        ax.add_patch(Polygon([[r0, z0], [r1, z0], [r1, z1], [r0, z1]], closed=True,
                             facecolor="#fbfbfb", edgecolor="none", zorder=1))
        draw_regions(ax, ol, s)
        for role, (col, _) in ROLECOL.items():
            pts = np.array(roles.get(role, []))
            if len(pts) == 0:
                continue
            ax.plot(pts[:, 0] * s, pts[:, 1] * s, ".", color=col, ms=3.0, zorder=5)
        # The feed boundary is eleven nodes on a five-micron disc; at any scale
        # that shows the device it would be invisible as dots alone, and it is
        # the one boundary this phase is about.
        p = np.array(ol["boundary"]["liquid_column_cut"]) * s
        ax.plot(p[:, 0], p[:, 1], "-", color=ROLECOL["liquid_feed_boundary"][0], lw=3.0,
                solid_capstyle="butt", zorder=8)
        p = np.array(ol["boundary"]["free_surface_reference"]) * s
        ax.plot(p[:, 0], p[:, 1], "-", color=ROLECOL["liquid_conductor"][0], lw=2.4,
                solid_capstyle="butt", zorder=8)
        # Named boundary curves that carry NO condition, drawn as the dielectric
        # interfaces they are.
        for nm in ("emitter_tip_land", "emitter_outer_surface", "emitter_rear_face",
                   "extractor_back", "extractor_rim"):
            if nm in ol["boundary"]:
                p = np.array(ol["boundary"][nm]) * s
                ax.plot(p[:, 0], p[:, 1], "--", color="#444444", lw=1.6, zorder=4)
        for nm, pts in ol["feature"].items():
            p = np.array(pts) * s
            ax.plot(p[:, 0], p[:, 1], "x", color="#c1272d", ms=7, mew=1.8, zorder=6)
        ax.axvline(0, color="#999999", lw=0.8, ls=":", zorder=3)
        ax.set_xlim(r0, r1)
        ax.set_ylim(z0, z1)
        ax.set_aspect("equal")
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        ax.set_title(name, fontsize=10)
    zf = par["base_z"] * s
    axes[1].annotate("Säulenschnitt:\n$\\varphi = V_\\mathrm{emitter}$ NUR auf\n"
                     "$r \\leq r_\\mathrm{Bohrung}$",
                     xy=(0.4 * par["r_bore"] * s, zf),
                     xytext=(0.95 * par["r_foot"] * s, zf + 0.34 * abs(zf)),
                     fontsize=7.4, color="#0f6b3d",
                     arrowprops=dict(arrowstyle="->", color="#0f6b3d", lw=1.1))
    axes[1].annotate("Rückfläche des Polymers –\nKEINE Elektrode",
                     xy=(1.6 * par["r_bore"] * s, zf),
                     xytext=(0.95 * par["r_foot"] * s, zf + 0.09 * abs(zf)), fontsize=7.4,
                     color="#7a3b00", arrowprops=dict(arrowstyle="->", color="#7a3b00", lw=1.1))
    axes[2].annotate("gepinnte Austrittskante,\nunverrundet",
                     xy=(par["r_bore"] * s, 0.0),
                     xytext=(1.15 * par["r_bore"] * s, 1.35 * par["r_bore"] * s), fontsize=7.4,
                     color="#c1272d", arrowprops=dict(arrowstyle="->", color="#c1272d", lw=1.1))

    handles = [Line2D([], [], marker="s", ls="none", ms=10, mfc=RCOL[k], mec="#333333",
                      label=RLABEL[k]) for k in ("liquid", "emitter_dielectric",
                                                 "extractor_carrier")]
    handles.append(Line2D([], [], marker="s", ls="none", ms=10, mfc="#fbfbfb",
                          mec="#333333", label="Vakuum, $\\varepsilon_r = 1$"))
    for role, (col, lab) in ROLECOL.items():
        if roles.get(role):
            handles.append(Line2D([], [], marker=".", ls="none", ms=9, color=col, label=lab))
    handles.append(Line2D([], [], ls="--", color="#444444", lw=1.6, label=FREE_LABEL))
    handles.append(Line2D([], [], marker="x", ls="none", color="#c1272d", ms=7, mew=1.8,
                          label="unverrundete Kanten – dort ist kein Spitzenfeld konvergent"))
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8.2, frameon=False,
               bbox_to_anchor=(0.5, 0.115))
    fig.suptitle("Materialgebiete und Randbedingungen\n" + TITLE, fontsize=11.5)
    fig.tight_layout(rect=[0, 0.33, 1, 0.94])
    caveat(fig, "Die Randbedingungspunkte sind aus node_roles.csv gezeichnet, also aus dem "
                "Zustand, den der Löser tatsächlich aufgebaut hat, nicht aus der Zeichnung. "
                f"Auf Polymerflächen liegen {m.get('polymer_dirichlet_nodes','?')} "
                "festgehaltene Knoten (Sollwert 0), auf der Schnittebene außerhalb der "
                f"Flüssigkeit {m.get('feed_plane_outside_liquid_nodes','?')} (Sollwert 0).",
            y=0.075)
    provenance(fig, m, stamp,
               f"Säulenschnitt bei z = {float(m['column_cut_z_m'])*1e6:.0f} µm "
               f"(base_plate_thickness = "
               f"{float(m['base_plate_thickness_m'])*1e6:.0f} µm) – VORLÄUFIGER "
               f"BEISPIELWERT; Metallisierung: {m.get('metallisation','?')}")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 2
def figure_mesh(d, m, par, ol, stamp, out):
    R, Z = mesh_nodes(d)
    s = 1e6
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 6.4),
                             gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]})
    views = [("Übersicht des automatischen Volumennetzes", 0,
              1.6 * par["extraction_distance"] * s, 1.4 * par["base_z"] * s,
              1.3 * (par["extraction_distance"] + par["extractor_thickness"]) * s, 1),
             ("Kapillare und Kegelflanke", 0, 3.0 * par["r_foot"] * s,
              -3.0 * par["r_foot"] * s, 3.0 * par["r_foot"] * s, 1),
             ("Austrittskante", 0, 2.2 * par["r_bore"] * s,
              -2.2 * par["r_bore"] * s, 2.2 * par["r_bore"] * s, 1)]
    for ax, (name, r0, r1, z0, z1, step) in zip(axes, views):
        draw_regions(ax, ol, s, alpha=0.55, lw=0.8)
        nz, nr = R.shape
        for j in range(0, nz, step):
            ax.plot(R[j, :] * s, Z[j, :] * s, "-", color="#2b6cb0", lw=0.32, alpha=0.85,
                    zorder=4)
        for i in range(0, nr, step):
            ax.plot(R[:, i] * s, Z[:, i] * s, "-", color="#2b6cb0", lw=0.32, alpha=0.85,
                    zorder=4)
        for nm, pts in ol["feature"].items():
            p = np.array(pts) * s
            ax.plot(p[:, 0], p[:, 1], "x", color="#c1272d", ms=8, mew=2.0, zorder=6)
        ax.set_xlim(r0, r1)
        ax.set_ylim(z0, z1)
        ax.set_aspect("equal")
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        ax.set_title(name, fontsize=10)
    checks = rows(os.path.join(d, "mesh_checks.csv"))
    npass = sum(1 for c in checks if c["passed"] == "1")
    fig.suptitle("Automatisches Volumennetz aus den Geräteparametern\n" + TITLE, fontsize=11.5)
    fig.tight_layout(rect=[0, 0.17, 1, 0.92])
    caveat(fig, "Blockstrukturiert und radial gewarpt: eine Gitterlinie liegt exakt auf der "
                "Kegelflanke, jede Materialgrenze fällt exakt auf Gitterlinien, keine "
                f"hängenden Knoten. {npass} von {len(checks)} Netzprüfungen bestanden, "
                "darunter die Reproduktion jedes Gebietsvolumens aus der geschlossenen "
                "Formel. Die Größenfunktion folgt allein der Geometrie; einziger "
                "Benutzerwert ist die Netzstufe. Die markierten Kanten sind unverrundet: "
                "dort divergiert das Feld und folgt der Elementgröße.", y=0.055)
    provenance(fig, m, stamp,
               f"{m.get('nr','?')} x {m.get('nz','?')} Knoten, Halbbandbreite aus dem "
               f"Bandlöser; Größenfeldfaktor {float(m.get('reference_size_scale',1)):.4f}")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 3
def figure_field(d, m, par, ol, stamp, out):
    s = 1e6
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 10.4))
    surf = rows(os.path.join(d, "reference_surface_field.csv"))
    for col, tag in enumerate(("uebersicht", "spitze")):
        ru, zu, PH, EM, RG = grid(os.path.join(d, f"field_{tag}.csv"))
        Rg, Zg = np.meshgrid(ru * s, zu * s)
        mask = RG == "liquid"

        ax = axes[0, col]
        ph = np.ma.array(PH, mask=mask)
        cf = ax.contourf(Rg, Zg, ph, levels=40, cmap="viridis")
        ax.contour(Rg, Zg, ph, levels=14, colors="w", linewidths=0.45, alpha=0.75)
        fig.colorbar(cf, ax=ax, label="$\\varphi$ [V]", fraction=0.046, pad=0.02)
        ax.set_title(f"Potential – {'Übersicht' if col == 0 else 'Spitzendetail'}",
                     fontsize=10)

        ax = axes[1, col]
        em = np.ma.array(EM, mask=mask | (EM <= 0))
        hi = float(np.nanmax(em.compressed()))
        # Three decades below the maximum.  The field spans five decades between
        # the exit edge and the far field, and a full-range scale shows nothing
        # but the edge; values outside the range are drawn with the end colours
        # (extend="both") rather than left blank.
        lo = hi / 1.0e3
        cf = ax.contourf(Rg, Zg, em, levels=np.geomspace(lo, hi, 40), norm=LogNorm(lo, hi),
                         cmap="magma", extend="both")
        cb = fig.colorbar(cf, ax=ax, label="$|E|$ [V/m]", fraction=0.046, pad=0.02,
                          ticks=matplotlib.ticker.LogLocator(numticks=6))
        cb.ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterSciNotation())
        ax.set_title(f"Feldstärke – {'Übersicht' if col == 0 else 'Spitzendetail'}",
                     fontsize=10)

        for ax in (axes[0, col], axes[1, col]):
            for name, pts in ol["region"].items():
                p = np.array(pts) * s
                ax.add_patch(Polygon(p, closed=True, facecolor="none", edgecolor="#ffffff",
                                     lw=1.3, zorder=5))
                ax.add_patch(Polygon(p, closed=True, facecolor="none", edgecolor="#222222",
                                     lw=0.7, zorder=6))
            p = np.array(ol["boundary"]["liquid_column_cut"]) * s
            ax.plot(p[:, 0], p[:, 1], "-", color="#0f9d58", lw=2.6, zorder=7)
            ax.set_xlim(ru[0] * s, ru[-1] * s)
            ax.set_ylim(zu[0] * s, zu[-1] * s)
            ax.set_aspect("equal")
            ax.set_xlabel("r [µm]")
            ax.set_ylabel("z [µm]")

    fig.suptitle("Potential und Feldstärke, SU-8-Referenzfall "
                 f"($\\varepsilon_r$ = {float(m['emitter_eps_r']):.2f}, "
                 f"{m['emitter_eps_r_status']})\n" + TITLE, fontsize=11.5)
    fig.tight_layout(rect=[0, 0.115, 1, 0.94])
    e0 = float(surf[0]["Ez_V_per_m"]) if surf else float("nan")
    caveat(fig, "Weiß umrandet: Flüssigkeit, Emitter-SU-8, Extraktorträger. Grün: die "
                "Zulaufgrenze, an der die Flüssigkeitssäule abgeschnitten ist. Im Inneren "
                "der Flüssigkeit ist |E| nicht definiert (idealer Leiter) und ausgeblendet. "
                f"Auf der ebenen Flüssigkeitsreferenz bei z = 0+ ist E_z(r=0) = {e0:.3e} V/m. "
                "Die Austrittskante ist unverrundet: das Feld dort divergiert, folgt der "
                "Elementgröße und wird nicht als Spitzenfeld berichtet.", y=0.045)
    provenance(fig, m, stamp,
               f"Q_emitter = {float(m['Q_emitter_C']):.4e} C, "
               f"Q_extractor = {float(m['Q_extractor_C']):.4e} C, "
               f"Summe = {float(m['Q_net_C']):.4e} C (nicht null: das System trägt "
               f"Nettoladung)")
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 4
def figure_convergence(d, m, stamp, out):
    span = abs(float(m["V_emitter_V"]) - float(m["V_extractor_V"]))
    conv = rows(os.path.join(d, "convergence_mesh.csv"))
    feed = rows(os.path.join(d, "convergence_feed.csv"))
    eps = rows(os.path.join(d, "sensitivity_permittivity.csv"))
    far = rows(os.path.join(d, "farfield_study.csv"))
    probes = [k[4:-2] for k in conv[0].keys() if k.startswith("phi_") and k.endswith("_V")]

    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.6))

    ax = axes[0, 0]
    n = np.array([float(r["nodes"]) for r in conv])
    for p in probes:
        v = np.array([float(r[f"phi_{p}_V"]) for r in conv])
        # Successive changes, not the distance to the finest level: the latter is
        # zero at the last point by construction and would draw a cliff that
        # means nothing.
        ax.semilogx(n[:-1], np.abs(np.diff(v)) / span, "o-", ms=4, label=p)
    ax.set_xlabel("Knoten der gröberen der beiden Stufen")
    ax.set_ylabel("$|\\Delta\\varphi|$ je Verfeinerung / $|V_e - V_x|$")
    ax.set_yscale("log")
    ax.set_title("Netzkonvergenz des Potentials an kantenfernen Punkten", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=6.6, ncol=2)
    q = np.array([float(r["Q_emitter_C"]) for r in conv])
    ax2 = ax.twinx()
    ax2.semilogx(n, q * 1e12, "s--", color="#c1272d", ms=4)
    ax2.set_ylabel("$Q_\\mathrm{emitter}$ [pC]", color="#c1272d")
    ax2.tick_params(axis="y", colors="#c1272d")

    ax = axes[0, 1]
    zf = np.array([abs(float(r["column_cut_z_m"])) for r in feed]) * 1e6
    for p in probes:
        v = np.array([float(r[f"phi_{p}_V"]) for r in feed])
        ax.semilogx(zf, v, "o-", ms=4, label=p)
    ax.set_xlabel("modellierte Säulenlänge $|z_\\mathrm{Schnitt}|$ [µm]")
    ax.set_ylabel("$\\varphi$ [V]")
    ax.set_title("Rückwärtige Ausdehnung des Säulenmodells –\n"
                 "KEINE Randverschiebung, sondern eine Geometrieänderung", fontsize=9.5,
                 color="#8a1b1b")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=6.6, ncol=2)

    ax = axes[1, 0]
    e = np.array([float(r["eps_r"]) for r in eps])
    key = probes[0]
    Ev = np.array([float(r[f"E_{key}_V_per_m"]) for r in eps])
    Qv = np.array([float(r["Q_emitter_C"]) for r in eps]) * 1e12
    nom = float(m["emitter_eps_r"])
    lo, hi = float(m["emitter_eps_r_low"]), float(m["emitter_eps_r_high"])
    ax.plot(e, Ev / 1e6, "o-", color="#2b6cb0", ms=4, label=f"$|E|$ bei {key}")
    if hi > lo:
        ax.axvspan(lo, hi, color="#f0c674", alpha=0.35,
                   label=f"begründeter Bereich {lo:g} … {hi:g}")
    ax.axvline(nom, color="#8a1b1b", ls="--", lw=1.4,
               label=f"Nominalwert {nom:g} (VORLÄUFIG)")
    ax.set_xlabel("$\\varepsilon_r$ von SU-8")
    ax.set_ylabel("$|E|$ [MV/m]", color="#2b6cb0")
    ax.set_title("Empfindlichkeit gegenüber der Permittivität", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    ax2 = ax.twinx()
    ax2.plot(e, Qv, "s--", color="#c1272d", ms=4)
    ax2.set_ylabel("$Q_\\mathrm{emitter}$ [pC]", color="#c1272d")
    ax2.tick_params(axis="y", colors="#c1272d")

    ax = axes[1, 1]
    Rd = np.array([float(r["domain_radius_m"]) for r in far]) * 1e3
    pa = np.array([float(r["phi_mid_asymptotic_V"]) for r in far])
    pg = np.array([float(r["phi_mid_grounded_V"]) for r in far])
    # The last point of each series is the reference itself; plotting it would
    # draw a cliff to zero that carries no information.
    ax.loglog(Rd[:-1], np.abs(pa[:-1] - pa[-1]) / span, "o-", ms=5,
              label="asymptotisch (Monopol-Robin)")
    ax.loglog(Rd[:-1], np.abs(pg[:-1] - pa[-1]) / span, "s-", ms=5, label="geerdete Hülle")
    ax.loglog(Rd, np.abs(pa - pg) / span, "^:", ms=5, color="#666666",
              label="Differenz der beiden = Trunkierungsmaß")
    ax.set_xlabel("Domänenradius [mm]")
    ax.set_ylabel("Abweichung / $|V_e - V_x|$")
    ax.set_title("Offene Fernrandbehandlung", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)

    fig.suptitle("Konvergenz und Empfindlichkeit\n" + TITLE, fontsize=11.5)
    fig.tight_layout(rect=[0, 0.135, 1, 0.93])
    caveat(fig,
           "BEFUND oben rechts: hier wurde nicht eine Randbedingung verschoben, sondern die "
           "Geometrie des leitenden Flüssigkeitskörpers geändert – Säule und Dielektrikum "
           "wachsen gemeinsam. Eine "
           f"Verdopplung der Säulenlänge verschiebt das Potential um "
           f"{float(m['column_change_phi_over_span']):.1e} der Spannweite und das Feld um "
           f"{float(m['column_change_E_rel']):.1e} relativ; die vorab festgelegte Grenze war "
           f"{float(m['tol_phi_over_span']):.0e}. Ursache: die Flüssigkeitssäule ist ein "
           "Leiter, eine längere Säule trägt proportional mehr Ladung. Eine geerdete Hülle "
           "ändert daran nichts. Alle Zahlen gelten für "
           f"einen Schnitt bei z = {float(m['column_cut_z_m'])*1e6:.0f} µm. Das "
           "Ersatzmodell mit vernetztem Flüssigkeitsraum rechnet es_reservoir. Der "
           "Nominalwert von "
           "$\\varepsilon_r$ ist vorläufig – daraus folgt keine Validierungsaussage.",
           y=0.038)
    provenance(fig, m, stamp,
               f"Netz: letzte Verfeinerung Δφ = "
               f"{float(m['mesh_change_phi_over_span']):.1e} der Spannweite, "
               f"Δ|E| = {float(m['mesh_change_E_rel']):.1e}; FEM gegen BEM bei "
               f"$\\varepsilon_r$ = 1: {float(m['fem_vs_bem_phi_over_span']):.1e} bzw. "
               f"{float(m['fem_vs_bem_E_rel']):.1e}")
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------- main
def main(d):
    m = meta(d)
    par = parameters(d)
    ol = outline(d)
    root = pv.repo_root_of(__file__)
    rel = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/")
    state = pv.repo_state(root, ignore=[rel])
    stamp = pv.stamp(state, m.get("commit"))
    pv.write_provenance(d, state, m.get("commit"))

    made = [figure_materials(d, m, par, ol, stamp,
                             os.path.join(d, "fig1_material_regions.png")),
            figure_mesh(d, m, par, ol, stamp, os.path.join(d, "fig2_volume_mesh.png")),
            figure_field(d, m, par, ol, stamp, os.path.join(d, "fig3_potential_field.png")),
            figure_convergence(d, m, stamp, os.path.join(d, "fig4_convergence.png"))]
    for f in made:
        print(f"{f}  {os.path.getsize(f)} Byte")
    print(f"Provenienz: {stamp}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
