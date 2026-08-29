#!/usr/bin/env python3
"""Figures of the P2a static vacuum electrostatics.

    ./build/es_vacuum examples/device_p1.cfg examples/vacuum_p2a.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_vacuum.py results/<dir>

Reads only the CSVs written by apps/es_vacuum.cpp.  Nothing here solves,
meshes or computes geometry.

Produces
  <dir>/fig1_potential.png       potential map, to scale, with equipotentials
  <dir>/fig2_field_magnitude.png |E| maps, sharp edges marked and cut out
  <dir>/fig3_convergence.png     mesh convergence and the electrode-radius study
  <dir>/fig4_surface_charge.png  surface charge density along the conductors
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
from matplotlib.patches import Circle, Polygon

TITLE = ("P2a: Vakuum-Elektrostatik, rho = 0, ebene Perfect-Conductor-"
         "Referenzfläche, keine Emission")
EPS0 = 8.8541878128e-12

RSTYLE = {
    "liquid": "#dceaf5",
    "emitter_solid": "#d0d0d0",
    "extractor_solid": "#f0d9c6",
}
PAINT_ORDER = ["emitter_solid", "liquid", "extractor_solid"]

ELECTRODE = {
    "emitter": ("#404040", "Emitter (Metall)"),
    "flat_liquid_surface_reference":
        ("#0f9d58", "anfängliche ebene Flüssigkeitsoberfläche (Perfect-Conductor-Referenz)"),
    "extractor": ("#e6550d", "Extraktor"),
}
EDGE_COLOUR = {"sharp_feature": "#c1272d", "truncation_end": "#7b3294"}


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


def regions(d):
    reg = defaultdict(lambda: defaultdict(list))
    for r in rows(os.path.join(d, "regions.csv")):
        reg[r["region"]][int(r["loop"])].append((float(r["r_m"]), float(r["z_m"])))
    return reg


def grid(path):
    """CSV written by write_grid_csv -> (r_axis, z_axis, V, |E|)."""
    rr, zz, vv, ee = [], [], [], []
    for r in rows(path):
        rr.append(float(r["r"])); zz.append(float(r["z"]))
        vv.append(float(r["V"])); ee.append(float(r["Emag"]))
    ur = np.unique(np.array(rr)); uz = np.unique(np.array(zz))
    shape = (uz.size, ur.size)
    return ur, uz, np.array(vv).reshape(shape), np.array(ee).reshape(shape)


def edge_zones(d):
    return [(r["kind"], r["name"], float(r["r_m"]), float(r["z_m"]), float(r["radius_m"]))
            for r in rows(os.path.join(d, "edge_zones.csv"))]


def convergence_summary(d):
    c = rows(os.path.join(d, "convergence.csv"))
    a, b = c[-2], c[-1]
    rel = lambda k: abs(float(b[k]) - float(a[k])) / abs(float(b[k]))
    return {"rows": c, "d_cEE": rel("c_EE_F"), "d_Cm": rel("C_mutual_F"),
            "d_Eref": rel("Ez_ref_V_per_m")}


# ------------------------------------------------------------------ drawing
def draw_regions(ax, reg, s, zorder=2):
    for name in PAINT_ORDER:
        if name not in reg:
            continue
        pts = [(rr * s, zz * s) for rr, zz in reg[name][0]]
        ax.add_patch(Polygon(pts, closed=True, facecolor=RSTYLE[name], edgecolor="0.3",
                             lw=0.9, zorder=zorder))


def cut_edge_zones(ax, zones, s):
    """White out the places where |E| is not a converged quantity."""
    for kind, _name, r, z, rad in zones:
        ax.add_patch(Circle((r * s, z * s), rad * s, facecolor="white",
                            edgecolor=EDGE_COLOUR[kind], lw=1.5, ls=(0, (3, 2)), zorder=6))


def edge_legend_handles():
    return [plt.Line2D([0], [0], ls=(0, (3, 2)), color=EDGE_COLOUR["sharp_feature"], lw=1.5,
                       label="unverrundete Gerätekante – ausgespart, |E| nicht konvergent"),
            plt.Line2D([0], [0], ls=(0, (3, 2)), color=EDGE_COLOUR["truncation_end"], lw=1.5,
                       label="offenes Bogenende des Modellschnitts – ausgespart, "
                             "|E| nicht konvergent")]


def provenance(fig, m, conv, extra="", y=0.012, width=175):
    head = (f"Commit {m.get('commit', 'unbekannt')[:12]}   |   "
            f"Konfiguration: {m.get('config', '?')}   |   "
            f"V_emitter = {float(m['V_emitter_V']):.6g} V, "
            f"V_extractor = {float(m['V_extractor_V']):.6g} V   |   "
            f"Netzstufe {m['reference_level']} (size_scale {m['reference_size_scale']}), "
            f"{m['n_bem_panels']} BEM-Panels   |   "
            f"Konvergenz zwischen den beiden feinsten Stufen: "
            f"|Δc_EE| = {conv['d_cEE']:.1e}, |ΔC_m| = {conv['d_Cm']:.1e}, "
            f"|ΔE_z(ref)| = {conv['d_Eref']:.1e}")
    txt = "\n".join(textwrap.wrap(head, width))
    if extra:
        txt += "\n" + "\n".join(textwrap.wrap(extra, width))
    fig.text(0.5, y, txt, ha="center", va="bottom", fontsize=7.2, color="0.25")


# ----------------------------------------------------------------- figure 1
def figure_potential(d, m, par, reg, conv, out):
    s, sT = 1e3, 1e6
    ur, uz, V, _ = grid(os.path.join(d, "field_full.csv"))
    tr, tz, tV, _ = grid(os.path.join(d, "field_tip.csv"))
    VE = float(m["V_emitter_V"])

    fig, (ax, az) = plt.subplots(1, 2, figsize=(13.8, 7.0))

    im = ax.pcolormesh(ur * s, uz * s, V, cmap="viridis", shading="nearest",
                       vmin=min(0.0, VE), vmax=max(0.0, VE), zorder=1)
    lv = np.sort(np.array([1, 2, 5, 10, 20, 50, 100, 200, 400, 700, 1000, 1300]) * (VE / 1500.0))
    cs = ax.contour(ur * s, uz * s, V, levels=lv, colors="white", linewidths=0.7, zorder=3)
    ax.clabel(cs, fmt="%g V", fontsize=6.4, inline=True)
    draw_regions(ax, reg, s)
    ax.add_patch(plt.Rectangle((0, tz.min() * s), tr.max() * s, (tz.max() - tz.min()) * s,
                               fill=False, ec="#c1272d", lw=1.4, zorder=7))
    ax.annotate("Emitterspitze\nAusschnitt (b)", (tr.max() * s, tz.max() * s), xytext=(84, 104),
                textcoords="offset points", fontsize=8.6, color="#c1272d", ha="center",
                arrowprops=dict(arrowstyle="->", color="#c1272d", lw=1.1),
                bbox=dict(fc="white", ec="#c1272d", alpha=0.95, pad=2.0), zorder=9)
    ax.annotate("Extraktor, V = %g V" % float(m["V_extractor_V"]),
                (0.6 * par["extractor_outer_radius"] * s,
                 (par["extraction_distance"] + 0.5 * par["extractor_thickness"]) * s),
                xytext=(-4, 62), textcoords="offset points", fontsize=8.6, color="0.1",
                ha="center", arrowprops=dict(arrowstyle="->", color="0.1", lw=1.1),
                bbox=dict(fc="white", ec="0.6", alpha=0.9, pad=1.8), zorder=8)
    ax.set_xlim(0, ur.max() * s)
    ax.set_ylim(uz.min() * s, uz.max() * s)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Radiale Koordinate  r  (mm)")
    ax.set_ylabel("Axiale Koordinate  z  (mm)")
    ax.set_title("(a) Potential in der Meridianhalbebene, maßstäblich  (Maßstab in mm)",
                 fontsize=10.2)
    fig.colorbar(im, ax=ax, fraction=0.043, pad=0.02, label="elektrisches Potential  V  (V)")

    # Own colour range: the zoom spans only the top of the global range.
    lo = float(np.nanmin(tV))
    im2 = az.pcolormesh(tr * sT, tz * sT, tV, cmap="viridis", shading="nearest",
                        vmin=lo, vmax=max(lo, VE), zorder=1)
    cs2 = az.contour(tr * sT, tz * sT, tV, levels=np.linspace(lo, VE, 11)[1:-1], colors="white",
                     linewidths=0.8, zorder=3)
    az.clabel(cs2, fmt="%.0f V", fontsize=6.6, inline=True)
    draw_regions(az, reg, sT)
    az.plot([0, par["contact_radius"] * sT], [0, 0], "-", color="#0f9d58", lw=3.0, zorder=5)
    az.annotate("anfängliche ebene Flüssigkeitsoberfläche\n– noch kein berechneter Meniskus,\n"
                "hier reine Perfect-Conductor-Referenz",
                (0.5 * par["contact_radius"] * sT, 0.0), xytext=(30, 96),
                textcoords="offset points", fontsize=8.4, color="#0b6b3d", ha="left",
                arrowprops=dict(arrowstyle="->", color="#0f9d58", lw=1.4),
                bbox=dict(fc="white", ec="#0f9d58", alpha=0.95, pad=2.4), zorder=9)
    az.set_xlim(0, tr.max() * sT)
    az.set_ylim(tz.min() * sT, tz.max() * sT)
    az.set_aspect("equal", adjustable="box")
    az.set_xlabel("Radiale Koordinate  r  (µm)")
    az.set_ylabel("Axiale Koordinate  z  (µm)")
    az.set_title("(b) Emitterspitze  (Maßstab in µm, eigene Farbskala – beides anders als (a))",
                 fontsize=10.2)
    fig.colorbar(im2, ax=az, fraction=0.043, pad=0.02, label="elektrisches Potential  V  (V)")

    fig.suptitle(TITLE, fontsize=11.8)
    provenance(fig, m, conv,
               "Weiß ist Vakuum; grau, hellblau und beige sind Emitterfestkörper, "
               "Flüssigkeitssäule und Elektrode, in denen die Randintegraldarstellung nicht "
               "ausgewertet wird. Beide Teilbilder haben verschiedene Längen- und Farbskalen.")
    fig.tight_layout(rect=(0, 0.065, 1, 0.945))
    fig.savefig(out, dpi=155)
    print("wrote", out)


# ----------------------------------------------------------------- figure 2
def figure_field(d, m, par, reg, conv, out):
    zones = edge_zones(d)
    panels = [
        ("field_full.csv", 1e3, "mm", "(a) gesamte Domäne"),
        ("field_tip.csv", 1e6, "µm", "(b) Emitterspitze"),
        ("field_aperture.csv", 1e6, "µm", "(c) Extraktoröffnung"),
        ("field_truncation.csv", 1e6, "µm", "(d) offenes Bogenende (Modellschnitt)"),
    ]
    fig, axgrid = plt.subplots(2, 2, figsize=(13.6, 11.2))
    axes = axgrid.ravel()

    for ax, (name, s, unit, title) in zip(axes, panels):
        ur, uz, _, E = grid(os.path.join(d, name))
        fin = E[np.isfinite(E)]
        # The colour range is set from the 99.5th percentile, not the maximum:
        # the corner values diverge with the mesh and must not set the scale.
        vmax = float(np.percentile(fin, 99.5))
        vmin = max(float(np.percentile(fin, 2.0)), vmax / 1e3)
        im = ax.pcolormesh(ur * s, uz * s, np.clip(E, vmin, vmax), cmap="inferno",
                           norm=LogNorm(vmin=vmin, vmax=vmax), shading="nearest", zorder=1)
        draw_regions(ax, reg, s)
        cut_edge_zones(ax, zones, s)
        ax.set_xlim(0, ur.max() * s)
        ax.set_ylim(uz.min() * s, uz.max() * s)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"Radiale Koordinate  r  ({unit})")
        ax.set_ylabel(f"Axiale Koordinate  z  ({unit})")
        ax.set_title(f"{title}\n(Maßstab in {unit}, eigene Farbskala)", fontsize=10.0)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02,
                     label="Feldstärke  |E|  (V/m), logarithmisch")

    axes[1].plot([0, par["contact_radius"] * 1e6], [0, 0], "-", color="#0f9d58", lw=3.0, zorder=5)
    axes[1].annotate("anfängliche ebene\nFlüssigkeitsoberfläche –\nnoch kein berechneter Meniskus",
                     (0.5 * par["contact_radius"] * 1e6, 0.0), xytext=(14, 74),
                     textcoords="offset points", fontsize=8.0, color="#0b6b3d",
                     arrowprops=dict(arrowstyle="->", color="#0f9d58", lw=1.3),
                     bbox=dict(fc="white", ec="#0f9d58", alpha=0.95, pad=2.0), zorder=9)
    axes[3].axhline(par["domain_z_min"] * 1e6, color="#7b3294", lw=2.6, zorder=8)
    axes[3].annotate("Modellschnitt z = z_min:\ndarunter ist nichts modelliert,\n"
                     "der Emitterbogen endet hier offen",
                     xy=(0.5 * par["phi_3"] * 1e6, par["domain_z_min"] * 1e6),
                     xycoords="data", xytext=(0.05, 0.90), textcoords="axes fraction",
                     fontsize=8.2, color="#5a1f70", va="top", ha="left",
                     arrowprops=dict(arrowstyle="->", color="#7b3294", lw=1.3),
                     bbox=dict(fc="white", ec="#7b3294", alpha=0.95, pad=2.4), zorder=9)
    fig.legend(handles=edge_legend_handles(), loc="lower center", ncol=1, fontsize=8.8,
               bbox_to_anchor=(0.5, 0.062), framealpha=0.95)

    fig.suptitle(TITLE, fontsize=11.8)
    provenance(fig, m, conv,
               "Die weiß ausgesparten Kreise sind die drei unverrundeten Gerätekanten und das "
               "offene Bogenende des Modellschnitts bei r = 20 µm, z = −400 µm. Dort divergiert "
               "|E| mit der Elementgröße; kein Wert von dort ist ein Ergebnis. Jedes Teilbild "
               "hat seine eigene Längen- und Farbskala.")
    fig.tight_layout(rect=(0, 0.125, 1, 0.955))
    fig.savefig(out, dpi=150)
    print("wrote", out)


# ----------------------------------------------------------------- figure 3
def figure_convergence(d, m, par, conv, out):
    c = conv["rows"]
    n = np.array([int(r["n_bem_panels"]) for r in c], float)
    cEE = np.array([float(r["c_EE_F"]) for r in c])
    Cm = np.array([float(r["C_mutual_F"]) for r in c])
    Er = np.array([float(r["Ez_ref_V_per_m"]) for r in c])
    res = np.array([float(r["residual_rel_clear"]) for r in c])
    st = rows(os.path.join(d, "extractor_radius_study.csv"))
    rext = np.array([float(r["extractor_outer_radius_m"]) for r in st])
    st_Cm = np.array([float(r["C_mutual_F"]) for r in st])
    st_E = np.array([float(r["Ez_ref_V_per_m"]) for r in st])

    fig, ax = plt.subplots(1, 3, figsize=(16.6, 5.6))

    a = ax[0]
    a.plot(n, cEE * 1e15, "o-", color="#1f6fb4")
    a.set_xscale("log")
    a.set_xlabel("Zahl der BEM-Panels  N  (1)")
    a.set_ylabel(r"$c_{EE}$  (fF)", color="#1f6fb4")
    a.tick_params(axis="y", labelcolor="#1f6fb4")
    a.grid(alpha=0.3, which="both")
    a2 = a.twinx()
    a2.plot(n, Cm * 1e15, "s--", color="#e6550d")
    a2.set_ylabel(r"$C_m = -c_{EX}$  (fF)", color="#e6550d")
    a2.tick_params(axis="y", labelcolor="#e6550d")
    a.set_title("(a) Netzkonvergenz der beiden Kapazitätsgrößen\n"
                r"($c_{EE}$ Maxwell-Selbstkoeffizient, $C_m$ Gegenkapazität)", fontsize=9.8)
    a.margins(y=0.18)
    a2.margins(y=0.18)
    for x, y in zip(n, cEE * 1e15):
        a.annotate(f"{y:.4f}", (x, y), textcoords="offset points", xytext=(0, 8), fontsize=7,
                   ha="center", color="#1f6fb4")
    for x, y in zip(n, Cm * 1e15):
        a2.annotate(f"{y:.4f}", (x, y), textcoords="offset points", xytext=(0, -14), fontsize=7,
                    ha="center", color="#e6550d")

    # Deviation from the finest level; the finest level itself is zero by
    # construction and is therefore not plotted.
    b = ax[1]
    k = slice(0, -1)
    b.plot(n[k], np.abs(cEE[k] - cEE[-1]) / abs(cEE[-1]), "o-", color="#1f6fb4", label=r"$c_{EE}$")
    b.plot(n[k], np.abs(Cm[k] - Cm[-1]) / abs(Cm[-1]), "s-", color="#e6550d", label=r"$C_m$")
    b.plot(n[k], np.abs(Er[k] - Er[-1]) / abs(Er[-1]), "^-", color="#0f9d58",
           label=r"$E_z$ am Referenzpunkt")
    b.plot(n, res, "d--", color="0.45", label="Potentialresiduum, relativ\n(Kantenzonen aus)")
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xlabel("Zahl der BEM-Panels  N  (1)")
    b.set_ylabel("relative Abweichung von der feinsten Stufe  (1)")
    b.set_title("(b) Konvergenz gegen die feinste Stufe\n"
                "(der feinste Punkt ist per Definition null und fehlt)", fontsize=9.8)
    b.grid(alpha=0.3, which="both")
    b.legend(fontsize=7.6, loc="lower left")
    b.text(0.97, 0.96,
           f"feinste Stufe, N = {int(n[-1])}\n"
           f"c_EE      = {cEE[-1]:.6e} F\n"
           f"C_m       = {Cm[-1]:.6e} F\n"
           f"E_z(ref)  = {Er[-1]:.6e} V/m",
           transform=b.transAxes, fontsize=7.8, va="top", ha="right", family="monospace",
           bbox=dict(fc="white", ec="0.7", alpha=0.95, pad=3.0))

    cx = ax[2]
    cx.plot(rext * 1e3, st_Cm * 1e15, "s-", color="#e6550d")
    cx.set_xlabel(r"Außenradius der Elektrode  $r_{ext}$  (mm)")
    cx.set_ylabel(r"Gegenkapazität  $C_m$  (fF)", color="#e6550d")
    cx.tick_params(axis="y", labelcolor="#e6550d")
    cx.grid(alpha=0.3)
    cy = cx.twinx()
    cy.plot(rext * 1e3, st_E * 1e-7, "^-", color="#0f9d58")
    cy.set_ylabel(r"$E_z$ am Referenzpunkt  ($10^{7}$ V/m)", color="#0f9d58")
    cy.tick_params(axis="y", labelcolor="#0f9d58")
    cx.set_xlim(0, par["domain_radius"] * 1e3 * 1.08)
    cx.axvline(par["extractor_outer_radius"] * 1e3, color="0.3", ls=":", lw=1.5)
    cx.axvline(par["domain_radius"] * 1e3, color="#7b3294", ls="--", lw=1.5)
    cx.text(par["extractor_outer_radius"] * 1e3 - 0.07, 0.06,
            f"Beispielwert\n$r_{{ext}}$ = {par['extractor_outer_radius']*1e3:.3g} mm",
            transform=cx.get_xaxis_transform(), fontsize=8.0, color="0.3", ha="right", va="bottom")
    cx.text(par["domain_radius"] * 1e3 - 0.07, 0.55,
            f"domain_radius = {par['domain_radius']*1e3:.3g} mm\nkein Leiter, keine "
            f"Randbedingung",
            transform=cx.get_xaxis_transform(), fontsize=8.0, color="#7b3294", ha="right",
            va="bottom")
    cx.set_title("(c) Einfluss des endlichen Elektrodenaußenradius\n"
                 "(Pflichtangabe, nicht mit dem Domänenrand gleichgesetzt)", fontsize=9.8)

    fig.suptitle(TITLE, fontsize=11.8)
    provenance(fig, m, conv,
               "Konvergenzgrößen: c_EE = Q_Emitter bei V_Emitter = 1 V und V_Extraktor = 0 V "
               "(Maxwell-Selbstkoeffizient), C_m = −c_EX (Gegenkapazität) und das Axialfeld "
               f"E_z am kantenfernen Referenzpunkt r = 0, z = {float(m['z_ref_m'])*1e6:.3g} µm. "
               "Kein unkommentiertes „C“: Q_Emitter/(V_E − V_X) ist bei geerdetem Extraktor "
               "identisch c_EE und nicht C_m, weil das System bei V(unendlich) = 0 Nettoladung "
               "trägt.")
    fig.tight_layout(rect=(0, 0.115, 1, 0.94))
    fig.savefig(out, dpi=150)
    print("wrote", out)


# ----------------------------------------------------------------- figure 4
def figure_surface(d, m, par, conv, out):
    sr = rows(os.path.join(d, "surface.csv"))
    by = defaultdict(lambda: defaultdict(list))
    for r in sr:
        t = by[r["tag"]]
        t["arc"].append(float(r["arclen_m"]) * 1e6)
        t["sig"].append(float(r["sigma_C_per_m2"]))
        t["ez"].append(int(r["in_edge_zone"]))
        t["r"].append(float(r["r_m"]) * 1e6)
        t["z"].append(float(r["z_m"]) * 1e6)

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(16.6, 5.6))

    for ax, tags, title in (
        (a, ["flat_liquid_surface_reference", "emitter"],
         "(a) Emitterseite: ebene Flüssigkeitsoberfläche und Emittermetall"),
        (b, ["extractor"], "(b) Extraktor"),
    ):
        for tag in tags:
            t = by[tag]
            col, lab = ELECTRODE[tag]
            arc = np.array(t["arc"]); sig = np.array(t["sig"]); ez = np.array(t["ez"], bool)
            short = {"flat_liquid_surface_reference":
                     "ebene Flüssigkeitsoberfläche (PC-Referenz)"}.get(tag, lab)
            ax.plot(arc[~ez], sig[~ez] * 1e6, ".", ms=5, color=col, label=short)
            if ez.any():
                ax.plot(arc[ez], sig[ez] * 1e6, "x", ms=6, mew=1.3, color=col)
        ax.plot([], [], "kx", ms=6, mew=1.3, label="in einer markierten Zone – nicht konvergent")
        ax.set_xlabel("Bogenlänge entlang der Leiterkontur  s  (µm)")
        ax.set_ylabel(r"Flächenladungsdichte  $\sigma$  (µC/m²)")
        ax.set_title(title, fontsize=9.8)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.4, loc="best")
    a.set_yscale("symlog", linthresh=1.0)

    t = by["flat_liquid_surface_reference"]
    rr = np.array(t["r"]); sig = np.array(t["sig"]); ez = np.array(t["ez"], bool)
    En = sig / EPS0
    c.plot(rr[~ez], En[~ez] * 1e-7, "o-", ms=5, color="#0f9d58",
           label="ebene Flüssigkeitsoberfläche")
    c.plot(rr[ez], En[ez] * 1e-7, "x", ms=8, mew=1.7, color="#c1272d",
           label="Kantenzone der gepinnten Austrittskante –\nnicht konvergent")
    c.axvline(par["contact_radius"] * 1e6, color="0.3", ls=":", lw=1.5)
    c.text(par["contact_radius"] * 1e6 - 0.08, 0.04,
           "gepinnte Austrittskante\n(unverrundet)", transform=c.get_xaxis_transform(),
           fontsize=8.2, color="0.3", ha="right", va="bottom")
    zr = float(m["z_ref_m"]) * 1e6
    c.set_xlabel("Radiale Koordinate  r  (µm)")
    c.set_ylabel(r"Normalfeld  $E_n = \sigma/\varepsilon_0$  ($10^{7}$ V/m)")
    c.set_title("(c) Normalfeld auf der anfänglichen ebenen Flüssigkeitsoberfläche\n"
                "– noch kein berechneter Meniskus, keine Emission", fontsize=9.8)
    c.grid(alpha=0.3)
    c.legend(fontsize=7.6, loc="upper left")
    c.text(0.04, 0.66, f"Referenzpunkt der Konvergenz:\nr = 0, z = {zr:.3g} µm",
           transform=c.transAxes, fontsize=7.8, family="monospace", va="top",
           bbox=dict(fc="white", ec="0.75", alpha=0.95, pad=2.5))

    fig.suptitle(TITLE, fontsize=11.8)
    provenance(fig, m, conv,
               "sigma ist die Dichte der Einfachschicht. Auf geschlossenen bzw. abgeschirmten "
               "Leiterflächen ist sigma/eps0 das einseitige Vakuum-Normalfeld; die Abschirmung "
               "des offenen Emitterbogens ist in report.txt beziffert. Die Vorzeichen folgen "
               "aus V_emitter > V_extractor.")
    fig.tight_layout(rect=(0, 0.105, 1, 0.94))
    fig.savefig(out, dpi=150)
    print("wrote", out)


def main(d):
    m = meta(d)
    par = parameters(d)
    reg = regions(d)
    conv = convergence_summary(d)
    figure_potential(d, m, par, reg, conv, os.path.join(d, "fig1_potential.png"))
    figure_field(d, m, par, reg, conv, os.path.join(d, "fig2_field_magnitude.png"))
    figure_convergence(d, m, par, conv, os.path.join(d, "fig3_convergence.png"))
    figure_surface(d, m, par, conv, os.path.join(d, "fig4_surface_charge.png"))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: plot_vacuum.py <results-directory>")
    main(sys.argv[1])
