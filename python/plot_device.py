#!/usr/bin/env python3
"""Figures of the parametric P1 device geometry.

    ./build/device_figure examples/device_p1.cfg results/<dir>
    python python/plot_device.py results/<dir>

Reads only regions.csv, boundaries.csv, features.csv and parameters.csv written
by tools/device_figure.cpp, so the figures are reproducible from source.
Nothing here computes geometry or physics.

Produces
  <dir>/fig1_device_overview.png   full r-z domain, to scale
  <dir>/fig2_device_detail.png     bore and exit edge; extractor aperture
"""
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Region -> (face colour, edge colour, label)
STYLE = {
    "vacuum":          ("#f7f7f7", "#bdbdbd", "Vakuum / offene Domäne"),
    "liquid":          ("#9ecae1", "#3182bd", "ionische Flüssigkeit"),
    "emitter_solid":   ("#969696", "#4d4d4d", "Emitter (massiv)"),
    "extractor_solid": ("#fdae6b", "#e6550d", "Extraktionselektrode"),
}
PAINT_ORDER = ["vacuum", "emitter_solid", "liquid", "extractor_solid"]


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def load(d):
    reg = defaultdict(lambda: defaultdict(list))
    for r in rows(os.path.join(d, "regions.csv")):
        reg[r["region"]][int(r["loop"])].append((float(r["r_m"]), float(r["z_m"])))
    bnd = defaultdict(list)
    for r in rows(os.path.join(d, "boundaries.csv")):
        bnd[r["name"]].append((float(r["r_m"]), float(r["z_m"]), r["id"]))
    feat = {r["name"]: (float(r["r_m"]), float(r["z_m"]))
            for r in rows(os.path.join(d, "features.csv"))}
    par = {r["name"]: float(r["value_SI"]) for r in rows(os.path.join(d, "parameters.csv"))}
    return reg, bnd, feat, par


def draw_regions(ax, reg, s):
    for name in PAINT_ORDER:
        if name not in reg:
            continue
        face, edge, _ = STYLE[name]
        pts = [(r * s, z * s) for r, z in reg[name][0]]
        ax.add_patch(Polygon(pts, closed=True, facecolor=face, edgecolor=edge,
                             lw=1.2, zorder=1 if name == "vacuum" else 2))


def draw_axis_and_open(ax, bnd, s, lim, label_axis=True):
    for name, pts in bnd.items():
        kind = pts[0][2]
        xy = [(r * s, z * s) for r, z, _ in pts]
        if kind == "symmetry_axis":
            ax.plot([p[0] for p in xy], [p[1] for p in xy], ls=(0, (7, 2, 1, 2)),
                    lw=1.6, color="crimson", zorder=4)
        elif kind == "open_boundary":
            ax.plot([p[0] for p in xy], [p[1] for p in xy], ls=(0, (5, 3)), lw=1.4,
                    color="#7b3294", zorder=3)
    if label_axis:
        ax.text(lim[0] * 0.01, lim[3] * 0.92, " Symmetrieachse r = 0", color="crimson",
                fontsize=8.5, rotation=90, va="top", ha="left")


def dim(ax, p0, p1, text, offset=(0, 0), color="0.2", fs=8.4):
    """Double-headed dimension arrow with a label."""
    ax.annotate("", xy=p1, xytext=p0,
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.1))
    mx, my = 0.5 * (p0[0] + p1[0]) + offset[0], 0.5 * (p0[1] + p1[1]) + offset[1]
    ax.text(mx, my, text, fontsize=fs, color=color, ha="center", va="center",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))


# ------------------------------------------------------------------ figure 1
def figure_overview(d, reg, bnd, feat, par, out):
    s = 1e3  # m -> mm
    R = par["domain_radius"] * s
    z0, z1 = par["domain_z_min"] * s, par["domain_z_max"] * s

    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    draw_regions(ax, reg, s)
    draw_axis_and_open(ax, bnd, s, (0, R, z0, z1))

    ze = par["extraction_distance"] * s
    t = par["extractor_thickness"] * s
    ra = 0.5 * par["extractor_aperture_diameter"] * s

    dim(ax, (0.35 * R, 0.0), (0.35 * R, ze), f"extraction_distance\n{ze * 1e3:.0f} µm",
        offset=(0.18 * R, 0))
    dim(ax, (0.0, ze), (ra, ze), f"Apertur ⌀ {2 * ra * 1e3:.0f} µm", offset=(0, -0.055 * z1))
    dim(ax, (0.72 * R, ze), (0.72 * R, ze + t), f"t = {t * 1e3:.0f} µm", offset=(0.14 * R, 0))
    ax.plot([0, ra], [ze, ze], ls=":", lw=1.0, color="0.35", zorder=5)

    ax.annotate("Emitter (Detail siehe Abb. 2)", xy=(0.03, 0.0), xytext=(0.42 * R, -0.55 * abs(z0)),
                fontsize=8.6, color="0.2",
                arrowprops=dict(arrowstyle="->", color="0.2", lw=1.1))

    handles = [Polygon([(0, 0)], facecolor=STYLE[k][0], edgecolor=STYLE[k][1], label=STYLE[k][2])
               for k in PAINT_ORDER]
    ax.legend(handles=handles, fontsize=8.4, loc="upper right", framealpha=0.95)

    ax.set_xlim(-0.03 * R, 1.03 * R)
    ax.set_ylim(z0 - 0.03 * (z1 - z0), z1 + 0.03 * (z1 - z0))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Radiale Koordinate  r  (mm)")
    ax.set_ylabel("Axiale Koordinate  z  (mm)")
    ax.set_title("Parametrische P1-Geometrie — Übersicht der r-z-Domäne (maßstäblich)",
                 fontsize=11)
    ax.grid(alpha=0.25, zorder=0)
    fig.text(0.5, 0.012,
             "Achsensymmetrisch; das äußere Rechteck ist die offene Rechendomäne, "
             "kein Leiter. Keine Randbedingungen in dieser Stufe.",
             ha="center", fontsize=8, color="0.35")
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(out, dpi=160)
    print("wrote", out)


# ------------------------------------------------------------------ figure 2
def figure_detail(d, reg, bnd, feat, par, out):
    s = 1e6  # m -> um
    p1 = par["phi_1"] * s
    p2 = par["phi_2"] * s
    p3 = par["phi_3"] * s
    H = par["emitter_height"] * s
    ze = par["extraction_distance"] * s
    t = par["extractor_thickness"] * s
    ra = 0.5 * par["extractor_aperture_diameter"] * s

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 6.6),
                                 gridspec_kw={"width_ratios": [1.0, 1.5]})

    # -- (a) bore, tip land, exit edge -------------------------------------
    draw_regions(a1, reg, s)
    draw_axis_and_open(a1, bnd, s, (0, 34, -H - 22, 34), label_axis=False)
    a1.plot([feat["pinned_contact_edge"][0] * s], [feat["pinned_contact_edge"][1] * s],
            "o", ms=11, mfc="none", mew=2.4, color="crimson", zorder=6)
    a1.annotate("gepinnte\nAustrittskante", (feat["pinned_contact_edge"][0] * s, 0),
                xytext=(11, -30), textcoords="offset points", fontsize=8.6, color="crimson",
                arrowprops=dict(arrowstyle="->", color="crimson", lw=1.0))
    a1.plot([feat["emitter_outer_edge"][0] * s], [0], "s", ms=8, mfc="none", mew=2.0,
            color="0.2", zorder=6)

    dim(a1, (0, 12), (0.5 * p2, 12), f"\u23002 = {p2:.0f} \u00b5m", offset=(0, 3.6),
        color="#3182bd")
    dim(a1, (0, 24), (0.5 * p1, 24), f"\u23001 = {p1:.0f} \u00b5m", offset=(0, 3.6))
    dim(a1, (0, -H - 11), (0.5 * p3, -H - 11), f"\u23003 = {p3:.0f} \u00b5m", offset=(0, -4.2))
    dim(a1, (28.5, -H), (28.5, 0), f"emitter_height\n{H:.0f} \u00b5m", offset=(0, 0))
    a1.plot([0, 0.5 * p3], [-H, -H], ls=":", lw=1.0, color="0.35", zorder=5)

    a1.set_xlim(-3, 34)
    a1.set_ylim(-H - 24, 34)
    a1.set_aspect("equal", adjustable="box")
    a1.set_xlabel("Radiale Koordinate  r  (\u00b5m)")
    a1.set_ylabel("Axiale Koordinate  z  (\u00b5m)")
    a1.set_title("(a) Bohrung, Stirnfl\u00e4che und Austrittskante", fontsize=10.5)
    a1.grid(alpha=0.25, zorder=0)

    # -- (b) extractor aperture --------------------------------------------
    draw_regions(a2, reg, s)
    draw_axis_and_open(a2, bnd, s, (0, 330, ze - 95, ze + t + 95), label_axis=False)
    for key in ("extractor_aperture_edge_front", "extractor_aperture_edge_back"):
        a2.plot([feat[key][0] * s], [feat[key][1] * s], "o", ms=10, mfc="none", mew=2.2,
                color="#e6550d", zorder=6)
    a2.annotate("Aperturkanten", (ra, ze), xytext=(16, -30), textcoords="offset points",
                fontsize=8.8, color="#e6550d",
                arrowprops=dict(arrowstyle="->", color="#e6550d", lw=1.0))
    dim(a2, (0, ze - 55), (ra, ze - 55), f"Apertur \u2300 {2 * ra:.0f} \u00b5m",
        offset=(0, -16))
    dim(a2, (ra + 85, ze), (ra + 85, ze + t), f"t = {t:.0f} \u00b5m", offset=(40, 0))

    a2.set_xlim(-8, 335)
    a2.set_ylim(ze - 100, ze + t + 100)
    a2.set_aspect("equal", adjustable="box")
    a2.set_xlabel("Radiale Koordinate  r  (\u00b5m)")
    a2.set_ylabel("Axiale Koordinate  z  (\u00b5m)")
    a2.set_title("(b) Extraktor\u00f6ffnung", fontsize=10.5)
    a2.grid(alpha=0.25, zorder=0)

    # -- shared legend, below both panels ----------------------------------
    handles = [Polygon([(0, 0)], facecolor=STYLE[k][0], edgecolor=STYLE[k][1], label=STYLE[k][2])
               for k in PAINT_ORDER]
    handles.append(plt.Line2D([0], [0], ls=(0, (7, 2, 1, 2)), color="crimson", lw=1.6,
                              label="Symmetrieachse r = 0"))
    handles.append(plt.Line2D([0], [0], ls=(0, (5, 3)), color="#7b3294", lw=1.4,
                              label="offener Domaenenrand"))
    fig.legend(handles=handles, fontsize=8.6, ncol=6, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), framealpha=0.95)

    fig.suptitle("Parametrische P1-Geometrie \u2014 Detail "
                 "(die beiden Ausschnitte haben verschiedene Ma\u00dfst\u00e4be)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.075, 1, 0.955))
    fig.savefig(out, dpi=160)
    print("wrote", out)


def main(d):
    reg, bnd, feat, par = load(d)
    figure_overview(d, reg, bnd, feat, par, os.path.join(d, "fig1_device_overview.png"))
    figure_detail(d, reg, bnd, feat, par, os.path.join(d, "fig2_device_detail.png"))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: plot_device.py <results-directory>")
    main(sys.argv[1])
