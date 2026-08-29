#!/usr/bin/env python3
"""Figures of the P2a static vacuum electrostatics.

    ./build/es_vacuum examples/device_p1.cfg examples/vacuum_p2a.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_vacuum.py results/<dir>

Reads only the CSVs written by apps/es_vacuum.cpp.  Nothing here solves,
meshes or computes geometry.

The commit named on every figure is measured from git at drawing time by
python/provenance.py, not taken from an argument, and a dirty working tree is
stamped DIRTY.  An earlier set of these figures carried a commit id that
predated the code they showed; that is what the measurement is for.

Produces
  <dir>/fig1_potential.png       potential map, to scale, with equipotentials
  <dir>/fig2_field_magnitude.png |E| maps, including the closed numerical back
  <dir>/fig3_convergence.png     mesh, electrode-radius and truncation studies
  <dir>/fig4_surface_charge.png  surface charge density along the conductors
  <dir>/figures_provenance.txt   what the stamp said and whether it is releasable
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance as pv

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
EDGE_COLOUR = {"sharp_feature": "#c1272d",
               "truncation_end": "#7b3294",
               "numerical_closure": "#1f6fb4"}
CLOSURE_COLOUR = "#1f6fb4"


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


def draw_back_closure(ax, m, par, s, lw=3.2):
    """The conducting disc that closes the emitter conductor."""
    zc = float(m["back_closure_z_m"]) * s
    r3 = 0.5 * par["phi_3"] * s
    ax.plot([0.0, r3], [zc, zc], "-", color=CLOSURE_COLOUR, lw=lw, solid_capstyle="butt",
            zorder=7)


def cut_edge_zones(ax, zones, s):
    """White out the places where no field value may be read off."""
    for kind, _name, r, z, rad in zones:
        ax.add_patch(Circle((r * s, z * s), rad * s, facecolor="white",
                            edgecolor=EDGE_COLOUR.get(kind, "#c1272d"), lw=1.5,
                            ls=(0, (3, 2)), zorder=8))


def edge_legend_handles(zones):
    kinds = {k for k, *_ in zones}
    h = []
    if "sharp_feature" in kinds:
        h.append(plt.Line2D([0], [0], ls=(0, (3, 2)), color=EDGE_COLOUR["sharp_feature"], lw=1.5,
                            label="unverrundete Gerätekante – ausgespart, |E| nicht konvergent"))
    if "numerical_closure" in kinds:
        h.append(plt.Line2D([0], [0], ls=(0, (3, 2)), color=EDGE_COLOUR["numerical_closure"],
                            lw=1.5,
                            label="Kante der numerischen Rückschließung – ausgespart, "
                                  "keine Gerätefläche"))
    if "truncation_end" in kinds:
        h.append(plt.Line2D([0], [0], ls=(0, (3, 2)), color=EDGE_COLOUR["truncation_end"], lw=1.5,
                            label="offenes Bogenende – darf nicht mehr auftreten"))
    h.append(plt.Line2D([0], [0], ls="-", color=CLOSURE_COLOUR, lw=3.0,
                        label="numerische Rückschließung des Leiters (voll leitende Scheibe "
                              "auf V_emitter, kein Bauteil)"))
    return h


def provenance(fig, m, conv, stamp, extra="", y=0.012, width=175):
    dirty = pv.DIRTY_MARK in stamp or stamp == pv.NO_VERSION
    head = (f"Commit {stamp}   |   "
            f"Konfiguration: {m.get('config', '?')}   |   "
            f"V_emitter = {float(m['V_emitter_V']):.6g} V, "
            f"V_extractor = {float(m['V_extractor_V']):.6g} V   |   "
            f"emitter_back_length = {float(m['emitter_back_length_m'])*1e6:.4g} µm "
            f"(Beispielwert, Pflichtangabe)   |   "
            f"Netzstufe {m['reference_level']} (size_scale {m['reference_size_scale']}), "
            f"{m['n_bem_panels']} BEM-Panels, davon {m['n_closure_panels']} Rückschließung   |   "
            f"Netzkonvergenz zwischen den beiden feinsten Stufen: "
            f"|Δc_EE| = {conv['d_cEE']:.1e}, |ΔC_m| = {conv['d_Cm']:.1e}, "
            f"|ΔE_z(ref)| = {conv['d_Eref']:.1e}")
    txt = "\n".join(textwrap.wrap(head, width))
    if extra:
        txt += "\n" + "\n".join(textwrap.wrap(extra, width))
    fig.text(0.5, y, txt, ha="center", va="bottom", fontsize=7.2,
             color="#a01010" if dirty else "0.25")
    if dirty:
        fig.text(0.5, 0.997, f"NICHT FREIGEGEBEN – {stamp}", ha="center", va="top",
                 fontsize=10.5, color="white", weight="bold",
                 bbox=dict(fc="#a01010", ec="#a01010", pad=3.0))


def truncation_note(m):
    ok = m.get("truncation_converged", "no") == "yes"
    if ok:
        return ("Trunkierungskonvergenz erreicht: E_z(ref), c_EX und die Sondenpotentiale "
                "ändern sich bei Verdopplung von emitter_back_length um weniger als die "
                "vorab festgelegte Toleranz.")
    return ("Trunkierungskonvergenz NICHT erreicht. Bei Verdopplung von emitter_back_length "
            f"ändert sich E_z(ref) um {float(m['truncation_change_Ez_ref']):.1e}, c_EX um "
            f"{float(m['truncation_change_c_EX']):.1e} und die Sondenpotentiale um "
            f"{float(m['truncation_change_V_probe']):.1e} – gegen eine vorab festgelegte "
            f"Toleranz von {float(m['truncation_tol_Ez_ref']):.0e}. Kein Wert dieser "
            "Abbildungen ist von der rückwärtigen Länge unabhängig; sie ist eine "
            "anzugebende Abmessung, kein Konvergenzparameter.")


# ----------------------------------------------------------------- figure 1
def figure_potential(d, m, par, reg, conv, stamp, out):
    s, sT = 1e3, 1e6
    ur, uz, V, _ = grid(os.path.join(d, "field_full.csv"))
    tr, tz, tV, _ = grid(os.path.join(d, "field_tip.csv"))
    VE = float(m["V_emitter_V"])

    fig, (ax, az) = plt.subplots(1, 2, figsize=(13.8, 7.4))

    im = ax.pcolormesh(ur * s, uz * s, V, cmap="viridis", shading="nearest",
                       vmin=min(0.0, VE), vmax=max(0.0, VE), zorder=1)
    lv = np.sort(np.array([1, 2, 5, 10, 20, 50, 100, 200, 400, 700, 1000, 1300]) * (VE / 1500.0))
    cs = ax.contour(ur * s, uz * s, V, levels=lv, colors="white", linewidths=0.7, zorder=3)
    ax.clabel(cs, fmt="%g V", fontsize=6.4, inline=True)
    draw_regions(ax, reg, s)
    draw_back_closure(ax, m, par, s, lw=2.6)
    ax.add_patch(plt.Rectangle((0, tz.min() * s), tr.max() * s, (tz.max() - tz.min()) * s,
                               fill=False, ec="#c1272d", lw=1.4, zorder=9))
    ax.annotate("Emitterspitze\nAusschnitt (b)", (tr.max() * s, tz.max() * s), xytext=(86, 96),
                textcoords="offset points", fontsize=8.6, color="#c1272d", ha="center",
                arrowprops=dict(arrowstyle="->", color="#c1272d", lw=1.1),
                bbox=dict(fc="white", ec="#c1272d", alpha=0.95, pad=2.0), zorder=10)
    ax.annotate("numerische Rückschließung:\nvolle leitende Scheibe auf V_emitter,\n"
                "kein Bauteil und kein Domänenrand",
                (0.5 * par["phi_3"] * 0.5 * s, float(m["back_closure_z_m"]) * s),
                xytext=(78, -34), textcoords="offset points", fontsize=8.4,
                color=CLOSURE_COLOUR, ha="left",
                arrowprops=dict(arrowstyle="->", color=CLOSURE_COLOUR, lw=1.3),
                bbox=dict(fc="white", ec=CLOSURE_COLOUR, alpha=0.95, pad=2.4), zorder=10)
    ax.annotate("Extraktor, V = %g V" % float(m["V_extractor_V"]),
                (0.6 * par["extractor_outer_radius"] * s,
                 (par["extraction_distance"] + 0.5 * par["extractor_thickness"]) * s),
                xytext=(26, 40), textcoords="offset points", fontsize=8.6, color="0.1",
                ha="left", arrowprops=dict(arrowstyle="->", color="0.1", lw=1.1),
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
    az.annotate("anfängliche ebene\nFlüssigkeitsoberfläche –\nnoch kein berechneter\n"
                "Meniskus, hier reine\nPerfect-Conductor-Referenz",
                (0.5 * par["contact_radius"] * sT, 0.0), xytext=(22, 78),
                textcoords="offset points", fontsize=8.2, color="#0b6b3d", ha="left",
                arrowprops=dict(arrowstyle="->", color="#0f9d58", lw=1.4),
                bbox=dict(fc="white", ec="#0f9d58", alpha=0.95, pad=2.4), zorder=9)
    az.set_xlim(0, tr.max() * sT)
    az.set_ylim(tz.min() * sT, tz.max() * sT)
    az.set_aspect("equal", adjustable="box")
    az.set_xlabel("Radiale Koordinate  r  (µm)")
    az.set_ylabel("Axiale Koordinate  z  (µm)")
    az.set_title("(b) Emitterspitze  (Maßstab in µm, eigene Farbskala –\n"
                 "Längen- und Farbskala anders als in (a))", fontsize=10.2)
    fig.colorbar(im2, ax=az, fraction=0.043, pad=0.02, label="elektrisches Potential  V  (V)")

    fig.suptitle(TITLE, fontsize=11.8, y=0.935)
    provenance(fig, m, conv, stamp,
               "Weiß ist Vakuum; grau, hellblau und beige sind Emitterfestkörper, "
               "Flüssigkeitssäule und Elektrode, in denen die Randintegraldarstellung nicht "
               "ausgewertet wird. Der Emitterleiter ist geschlossen: er endet bei "
               f"z = {float(m['back_closure_z_m'])*1e3:.4g} mm an der blau gezeichneten "
               "Scheibe. " + truncation_note(m))
    fig.tight_layout(rect=(0, 0.165, 1, 0.905))
    fig.savefig(out, dpi=155)
    print("wrote", out)


# ----------------------------------------------------------------- figure 2
def figure_field(d, m, par, reg, conv, stamp, out):
    zones = edge_zones(d)
    panels = [
        ("field_full.csv", 1e3, "mm", "(a) gesamte Domäne"),
        ("field_tip.csv", 1e6, "µm", "(b) Emitterspitze"),
        ("field_aperture.csv", 1e6, "µm", "(c) Extraktoröffnung"),
        ("field_back_closure.csv", 1e6, "µm",
         "(d) numerische Rückschließung – der Leiter ist hier zu"),
    ]
    fig, axgrid = plt.subplots(2, 2, figsize=(13.6, 11.6))
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
        draw_back_closure(ax, m, par, s, lw=3.2 if name.endswith("closure.csv") else 2.2)
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

    zc = float(m["back_closure_z_m"]) * 1e6
    axes[3].annotate("volle leitende Scheibe\nauf V_emitter:\nkein offenes Bogenende,\n"
                     "Inneres feldfrei,\nFläche NICHT auswertbar",
                     xy=(0.25 * par["phi_3"] * 1e6, zc), xycoords="data",
                     xytext=(0.04, 0.04), textcoords="axes fraction", fontsize=7.8,
                     color=CLOSURE_COLOUR, va="bottom", ha="left",
                     arrowprops=dict(arrowstyle="->", color=CLOSURE_COLOUR, lw=1.3),
                     bbox=dict(fc="white", ec=CLOSURE_COLOUR, alpha=0.95, pad=2.2), zorder=10)
    axes[3].annotate("Abstand zur ausgewerteten\nRegion (z ≥ %.4g µm): %.4g µm"
                     % (float(m["evaluation_z_min_m"]) * 1e6,
                        float(m["back_closure_clearance_m"]) * 1e6),
                     xy=(0.04, 0.965), xycoords="axes fraction", fontsize=7.8, color="0.2",
                     va="top", ha="left",
                     bbox=dict(fc="white", ec="0.7", alpha=0.95, pad=2.2), zorder=10)

    fig.legend(handles=edge_legend_handles(zones), loc="lower center", ncol=1, fontsize=8.6,
               bbox_to_anchor=(0.5, 0.068), framealpha=0.95)

    fig.suptitle(TITLE, fontsize=11.8, y=0.962)
    provenance(fig, m, conv, stamp,
               "Die weiß ausgesparten Kreise sind die vier unverrundeten Gerätekanten und die "
               "Kante der numerischen Rückschließung. Dort ist |E| kein Ergebnis: an den "
               "Gerätekanten divergiert es mit der Elementgröße, an der Rückschließung "
               "gehört die Fläche nicht zum Gerät. Jedes Teilbild hat seine eigene Längen- "
               "und Farbskala.")
    fig.tight_layout(rect=(0, 0.135, 1, 0.942))
    fig.savefig(out, dpi=150)
    print("wrote", out)


# ----------------------------------------------------------------- figure 3
def figure_convergence(d, m, par, conv, stamp, out):
    c = conv["rows"]
    n = np.array([int(r["n_bem_panels"]) for r in c], float)
    cEE = np.array([float(r["c_EE_F"]) for r in c])
    Cm = np.array([float(r["C_mutual_F"]) for r in c])
    Er = np.array([float(r["Ez_ref_V_per_m"]) for r in c])
    res = np.array([float(r["residual_rms_physical"]) for r in c])
    st = rows(os.path.join(d, "extractor_radius_study.csv"))
    rext = np.array([float(r["extractor_outer_radius_m"]) for r in st])
    st_Cm = np.array([float(r["C_mutual_F"]) for r in st])
    st_E = np.array([float(r["Ez_ref_V_per_m"]) for r in st])
    tr = rows(os.path.join(d, "truncation.csv"))
    tL = np.array([float(r["emitter_back_length_m"]) for r in tr])
    t_E = np.array([float(r["Ez_ref_V_per_m"]) for r in tr])
    t_cEX = np.array([float(r["c_EX_F"]) for r in tr])
    t_cEE = np.array([float(r["c_EE_F"]) for r in tr])
    tol = float(m["truncation_tol_Ez_ref"])

    fig, axg = plt.subplots(2, 2, figsize=(15.2, 11.4))
    a, b, cx, tx = axg[0, 0], axg[0, 1], axg[1, 0], axg[1, 1]

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
    a.set_title("(a) NETZkonvergenz der beiden Kapazitätsgrößen\n"
                r"($c_{EE}$ Maxwell-Selbstkoeffizient, $C_m$ Gegenkapazität)"
                "\nbei fester rückwärtiger Länge", fontsize=9.8)
    # No offset notation: a flat line plus the annotated digits says "converged"
    # far more directly than an axis labelled "+1.425e1".
    for axis in (a.yaxis, a2.yaxis):
        axis.get_major_formatter().set_useOffset(False)
    a.margins(y=0.30)
    a2.margins(y=0.30)
    for x, y in zip(n, cEE * 1e15):
        a.annotate(f"{y:.4f}", (x, y), textcoords="offset points", xytext=(0, 8), fontsize=7,
                   ha="center", color="#1f6fb4")
    for x, y in zip(n, Cm * 1e15):
        a2.annotate(f"{y:.4f}", (x, y), textcoords="offset points", xytext=(0, -14), fontsize=7,
                    ha="center", color="#e6550d")

    # Deviation from the finest level; the finest level itself is zero by
    # construction and is therefore not plotted.
    k = slice(0, -1)
    b.plot(n[k], np.abs(cEE[k] - cEE[-1]) / abs(cEE[-1]), "o-", color="#1f6fb4", label=r"$c_{EE}$")
    b.plot(n[k], np.abs(Cm[k] - Cm[-1]) / abs(Cm[-1]), "s-", color="#e6550d", label=r"$C_m$")
    b.plot(n[k], np.abs(Er[k] - Er[-1]) / abs(Er[-1]), "^-", color="#0f9d58",
           label=r"$E_z$ am Referenzpunkt")
    b.plot(n, res, "d--", color="0.45",
           label="RMS-Potentialresiduum, relativ\n(auswertbare Geräteflächen)")
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xlabel("Zahl der BEM-Panels  N  (1)")
    b.set_ylabel("relative Abweichung von der feinsten Stufe  (1)")
    b.set_title("(b) NETZkonvergenz gegen die feinste Stufe\n"
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
    cx.set_title("(c) GEOMETRIEabhängigkeit: endlicher Elektrodenaußenradius\n"
                 "(Pflichtangabe, nicht mit dem Domänenrand gleichgesetzt)", fontsize=9.8)

    # (d) the truncation study -- the panel this correction exists for.
    tx.plot(tL * 1e6, t_E * 1e-7, "^-", color="#0f9d58", label=r"$E_z$ am Referenzpunkt")
    tx.set_xscale("log")
    tx.set_xlabel(r"rückwärtige Länge  emitter_back_length  (µm)")
    tx.set_ylabel(r"$E_z$ am Referenzpunkt  ($10^{7}$ V/m)", color="#0f9d58")
    tx.tick_params(axis="y", labelcolor="#0f9d58")
    tx.grid(alpha=0.3, which="both")
    ty = tx.twinx()
    ty.plot(tL * 1e6, np.abs(t_cEX) * 1e15, "s--", color="#e6550d", label=r"$|c_{EX}|$")
    ty.plot(tL * 1e6, t_cEE * 1e15, "o:", color="#1f6fb4", label=r"$c_{EE}$")
    ty.set_yscale("log")
    ty.set_ylabel(r"$|c_{EX}|$ und $c_{EE}$  (fF, logarithmisch)")
    L0 = float(m["emitter_back_length_m"])
    tx.axvline(L0 * 1e6, color="0.3", ls=":", lw=1.5)
    tx.text(L0 * 1e6 * 0.95, 0.97, f"benutzter Wert\n{L0*1e6:.4g} µm\n(Beispielwert)",
            transform=tx.get_xaxis_transform(), fontsize=8.0, color="0.3", ha="right", va="top")
    # Headroom at the bottom so the verdict box does not sit on the curves.
    tx.margins(y=0.10)
    lo, hi = tx.get_ylim()
    tx.set_ylim(lo - 0.55 * (hi - lo), hi)
    lo2, hi2 = ty.get_ylim()
    ty.set_ylim(lo2 * (lo2 / hi2) ** 0.55, hi2)
    hl, la = tx.get_legend_handles_labels()
    h2, l2 = ty.get_legend_handles_labels()
    tx.legend(hl + h2, la + l2, fontsize=7.8, loc="upper left")
    ok = m.get("truncation_converged", "no") == "yes"
    tx.set_title("(d) TRUNKIERUNGSstudie gegen die rückwärtige Länge\n"
                 "(Netz an Spitze und Extraktor über alle Längen bitweise identisch: %s)"
                 % ("ja" if m.get("front_mesh_identical") == "yes" else "NEIN"), fontsize=9.8)
    tx.text(0.97, 0.04,
            "Toleranz vorab: %.0e\n"
            "gemessen bei Verdopplung:\n"
            "  E_z(ref) %.2e\n"
            "  c_EX     %.2e\n"
            "  V(Sonden) %.2e\n"
            "  c_EE     %.2e (ohne Toleranz)\n"
            "%s"
            % (tol, float(m["truncation_change_Ez_ref"]), float(m["truncation_change_c_EX"]),
               float(m["truncation_change_V_probe"]), float(m["truncation_change_c_EE"]),
               "KONVERGIERT" if ok else "NICHT truncation-konvergiert"),
            transform=tx.transAxes, fontsize=7.8, va="bottom", ha="right", family="monospace",
            bbox=dict(fc="white", ec="#a01010" if not ok else "0.7", lw=1.6 if not ok else 1.0,
                      alpha=0.97, pad=3.0))

    fig.suptitle(TITLE, fontsize=11.8, y=0.962)
    provenance(fig, m, conv, stamp,
               "(a) und (b) zeigen NETZkonvergenz bei fester Geometrie – dort konvergiert "
               "alles. (c) und (d) zeigen Abhängigkeiten von der Gerätegeometrie: vom "
               "Elektrodenaußenradius und von der rückwärtigen Emitterlänge. Diese sind "
               "keine Netzfehler und verschwinden nicht mit Verfeinerung. "
               "Kein unkommentiertes „C“: Q_Emitter/(V_E − V_X) ist bei geerdetem Extraktor "
               "identisch c_EE und nicht C_m, weil das System bei V(unendlich) = 0 "
               "Nettoladung trägt. " + truncation_note(m))
    fig.tight_layout(rect=(0, 0.125, 1, 0.942))
    fig.savefig(out, dpi=150)
    print("wrote", out)


# ----------------------------------------------------------------- figure 4
def figure_surface(d, m, par, conv, stamp, out):
    sr = rows(os.path.join(d, "surface.csv"))
    by = defaultdict(lambda: defaultdict(list))
    for r in sr:
        key = "numerical_closure" if r["numerical"] == "1" else r["tag"]
        t = by[key]
        t["arc"].append(float(r["arclen_m"]) * 1e6)
        t["sig"].append(float(r["sigma_C_per_m2"]))
        t["ez"].append(int(r["in_edge_zone"]))
        t["ev"].append(int(r["evaluable"]))
        t["r"].append(float(r["r_m"]) * 1e6)
        t["z"].append(float(r["z_m"]) * 1e6)

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(16.6, 6.6))

    for ax, tags, title in (
        (a, ["flat_liquid_surface_reference", "emitter", "numerical_closure"],
         "(a) Emitterseite: ebene Flüssigkeitsoberfläche, Emittermetall\nund die numerische "
         "Rückschließung"),
        (b, ["extractor"], "(b) Extraktor"),
    ):
        for tag in tags:
            if tag not in by:
                continue
            t = by[tag]
            if tag == "numerical_closure":
                col, short = CLOSURE_COLOUR, "numerische Rückschließung (kein Bauteil)"
            else:
                col, lab = ELECTRODE[tag]
                short = {"flat_liquid_surface_reference":
                         "ebene Flüssigkeitsoberfläche (PC-Referenz)"}.get(tag, lab)
            arc = np.array(t["arc"]); sig = np.array(t["sig"])
            ev = np.array(t["ev"], bool)
            ax.plot(arc[ev], sig[ev] * 1e6, ".", ms=5, color=col, label=short)
            if (~ev).any():
                ax.plot(arc[~ev], sig[~ev] * 1e6, "x", ms=6, mew=1.3, color=col)
        ax.plot([], [], "kx", ms=6, mew=1.3,
                label="nicht auswertbar: Kantenzone, Rückschließung\noder Schaft hinter dem "
                      "Kegelfuß")
        # Logarithmic arc length: the flat liquid surface and the tip land are five
        # micrometres each and would otherwise vanish next to a shank of hundreds.
        ax.set_xscale("log")
        ax.set_xlabel("Bogenlänge entlang der Leiterkontur  s  (µm, logarithmisch)")
        ax.set_ylabel(r"Flächenladungsdichte  $\sigma$  (µC/m²)")
        ax.set_title(title, fontsize=9.6)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7.0, loc="lower left", framealpha=0.95)
    a.set_yscale("symlog", linthresh=1.0)
    a.margins(y=0.25)

    t = by["flat_liquid_surface_reference"]
    rr = np.array(t["r"]); sig = np.array(t["sig"]); ev = np.array(t["ev"], bool)
    En = sig / EPS0
    c.plot(rr[ev], En[ev] * 1e-7, "o-", ms=5, color="#0f9d58",
           label="ebene Flüssigkeitsoberfläche")
    c.plot(rr[~ev], En[~ev] * 1e-7, "x", ms=8, mew=1.7, color="#c1272d",
           label="Kantenzone der gepinnten Austrittskante –\nnicht konvergent")
    c.axvline(par["contact_radius"] * 1e6, color="0.3", ls=":", lw=1.5)
    c.text(par["contact_radius"] * 1e6 - 0.08, 0.04,
           "gepinnte Austrittskante\n(unverrundet)", transform=c.get_xaxis_transform(),
           fontsize=8.2, color="0.3", ha="right", va="bottom")
    zr = float(m["z_ref_m"]) * 1e6
    c.set_xlabel("Radiale Koordinate  r  (µm)")
    c.set_ylabel(r"Normalfeld  $E_n = \sigma/\varepsilon_0$  ($10^{7}$ V/m)")
    c.set_title("(c) Normalfeld auf der anfänglichen ebenen Flüssigkeitsoberfläche\n"
                "– noch kein berechneter Meniskus, keine Emission", fontsize=9.6)
    c.grid(alpha=0.3)
    c.legend(fontsize=7.6, loc="upper left")
    c.text(0.04, 0.66,
           f"Referenzpunkt der Konvergenz:\nr = 0, z = {zr:.3g} µm\n"
           f"gilt für emitter_back_length\n= {float(m['emitter_back_length_m'])*1e6:.4g} µm",
           transform=c.transAxes, fontsize=7.8, family="monospace", va="top",
           bbox=dict(fc="white", ec="0.75", alpha=0.95, pad=2.5))

    fig.suptitle(TITLE, fontsize=11.8, y=0.930)
    provenance(fig, m, conv, stamp,
               "sigma ist die Dichte der Einfachschicht. Der Leiter ist geschlossen, also ist "
               "sigma/eps0 das einseitige Vakuum-Normalfeld – aber nur dort, wo die Fläche das "
               "Gerät begrenzt. Die numerische Rückschließung und der Schaft hinter dem "
               "Kegelfuß sind mitgelöst und mitgezeichnet, liefern aber keinen berichteten "
               "Feldwert. Die Vorzeichen folgen aus V_emitter > V_extractor.")
    fig.tight_layout(rect=(0, 0.155, 1, 0.900))
    fig.savefig(out, dpi=150)
    print("wrote", out)


def main(d):
    m = meta(d)
    par = parameters(d)
    reg = regions(d)
    conv = convergence_summary(d)
    root = pv.repo_root_of(__file__)
    # The one exclusion: this run's own output directory.  Its files are what the
    # run produces, so demanding that they be committed first is circular.
    here = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/") + "/"
    state = pv.repo_state(root, ignore=(here,))
    stamp = pv.stamp(state, m.get("commit"))
    text = pv.write_provenance(d, state, m.get("commit"))
    print("Provenienz:", text,
          "-- freigebbar" if pv.clean_release(state, m.get("commit")) else "-- NICHT freigebbar")
    figure_potential(d, m, par, reg, conv, stamp, os.path.join(d, "fig1_potential.png"))
    figure_field(d, m, par, reg, conv, stamp, os.path.join(d, "fig2_field_magnitude.png"))
    figure_convergence(d, m, par, conv, stamp, os.path.join(d, "fig3_convergence.png"))
    figure_surface(d, m, par, conv, stamp, os.path.join(d, "fig4_surface_charge.png"))
    if not pv.clean_release(state, m.get("commit")):
        print("ACHTUNG: Arbeitsbaum schmutzig oder Daten nicht aus HEAD -- "
              "die Abbildungen sind als NICHT FREIGEGEBEN gekennzeichnet.")
        return 1
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: plot_vacuum.py <results-directory>")
    sys.exit(main(sys.argv[1]))
