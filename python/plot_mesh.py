#!/usr/bin/env python3
"""Figures of the automatic axisymmetric P1 boundary mesh.

    ./build/mesh_figure examples/device_p1.cfg results/<dir>
    python python/plot_mesh.py results/<dir>

Reads only the CSVs written by tools/mesh_figure.cpp -- regions.csv and
features.csv for the background, mesh_nodes.csv, mesh_elements.csv and
mesh_boundaries.csv for the mesh itself.  Nothing here meshes, computes
geometry, or solves anything.

Produces
  <dir>/fig1_mesh_overview.png    whole boundary mesh, to scale, plus statistics
  <dir>/fig2_mesh_tip_detail.png  bore, initial flat liquid surface, exit edge
  <dir>/fig3_mesh_aperture.png    extractor aperture
"""
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

# Boundary identifier -> (colour, short label for the legend)
BSTYLE = {
    "symmetry_axis":          ("#d6221f", "Symmetrieachse r = 0"),
    "emitter_outer_surface":  ("#404040", "Emitteraußenfläche"),
    "emitter_tip_land":       ("#8c6d31", "Stirnfläche"),
    "bore_wall":              ("#1f6fb4", "Bohrungswand"),
    "free_surface_reference": ("#0f9d58", "anfängliche ebene Flüssigkeitsoberfläche"),
    "liquid_inlet":           ("#66c2a5", "Zulaufschnitt"),
    "extractor_surface":      ("#e6550d", "Extraktorflächen / Apertur"),
    "numerical_emitter_back_closure":
                              ("#1f6fb4", "numerische Rückschließung (kein Bauteil)"),
    "open_boundary":          ("#7b3294", "offener Domänenrand"),
}
ORDER = list(BSTYLE)

RSTYLE = {
    "vacuum":          "#fbfbfb",
    "liquid":          "#dceaf5",
    "emitter_solid":   "#e3e3e3",
    "extractor_solid": "#fdeadb",
}
PAINT_ORDER = ["vacuum", "emitter_solid", "liquid", "extractor_solid"]


def rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def load(d):
    reg = defaultdict(lambda: defaultdict(list))
    for r in rows(os.path.join(d, "regions.csv")):
        reg[r["region"]][int(r["loop"])].append((float(r["r_m"]), float(r["z_m"])))
    nodes = [(float(r["r_m"]), float(r["z_m"]), r["feature"])
             for r in rows(os.path.join(d, "mesh_nodes.csv"))]
    elems = [(float(r["r_a_m"]), float(r["z_a_m"]), float(r["r_b_m"]), float(r["z_b_m"]),
              r["boundary_id"], r["kind"])
             for r in rows(os.path.join(d, "mesh_elements.csv"))]
    feat = {r["name"]: (float(r["r_m"]), float(r["z_m"]))
            for r in rows(os.path.join(d, "features.csv"))}
    par = {r["name"]: float(r["value_SI"]) for r in rows(os.path.join(d, "parameters.csv"))}
    stats = [r for r in rows(os.path.join(d, "mesh_boundaries.csv"))]
    return reg, nodes, elems, feat, par, stats


def draw_regions(ax, reg, s):
    for name in PAINT_ORDER:
        if name not in reg:
            continue
        pts = [(r * s, z * s) for r, z in reg[name][0]]
        ax.add_patch(Polygon(pts, closed=True, facecolor=RSTYLE[name], edgecolor="none",
                             zorder=0 if name == "vacuum" else 1))


def draw_mesh(ax, nodes, elems, s, lw=1.6, ms=2.4, node_edge=False):
    """Elements as coloured segments, nodes as dots on top."""
    by_id = defaultdict(list)
    for ra, za, rb, zb, bid, _kind in elems:
        by_id[bid].append([(ra * s, za * s), (rb * s, zb * s)])
    for bid in ORDER:
        if bid not in by_id:
            continue
        colour = BSTYLE[bid][0]
        ax.add_collection(LineCollection(by_id[bid], colors=colour, linewidths=lw, zorder=3))
    ax.plot([n[0] * s for n in nodes], [n[1] * s for n in nodes], ls="none", marker="o",
            ms=ms, mfc="white" if node_edge else "0.1",
            mec="0.1", mew=0.5 if node_edge else 0.0, zorder=4)


def legend_handles(stats_by_id, with_counts=True):
    out = []
    for bid in ORDER:
        label = BSTYLE[bid][1]
        if with_counts and bid in stats_by_id:
            label += f"  (n = {stats_by_id[bid]['n']})"
        out.append(plt.Line2D([0], [0], color=BSTYLE[bid][0], lw=2.4, label=label))
    out.append(plt.Line2D([0], [0], ls="none", marker="o", ms=4, mfc="0.1", mec="0.1",
                          label="Netzknoten"))
    return out


def id_stats(stats):
    out = {}
    for r in stats:
        if r["scope"] == "boundary_id":
            out[r["name"]] = {
                "n": int(r["n_elements"]),
                "min": float(r["len_min_m"]),
                "med": float(r["len_median_m"]),
                "max": float(r["len_max_m"]),
            }
    return out


# ------------------------------------------------------------------ figure 1
def figure_overview(reg, nodes, elems, feat, par, stats, out):
    s = 1e3  # m -> mm
    R = par["domain_radius"] * s
    z0, z1 = par["domain_z_min"] * s, par["domain_z_max"] * s
    st = id_stats(stats)
    total = next(r for r in stats if r["scope"] == "total")

    fig, (ax, tx) = plt.subplots(1, 2, figsize=(13.4, 7.0),
                                 gridspec_kw={"width_ratios": [1.18, 1.0]})

    draw_regions(ax, reg, s)
    draw_mesh(ax, nodes, elems, s, lw=1.3, ms=1.7)
    for key, colour in (("pinned_contact_edge", "#0f9d58"),
                        ("emitter_outer_edge", "#8c6d31"),
                        ("extractor_aperture_edge_front", "#e6550d"),
                        ("extractor_aperture_edge_back", "#e6550d")):
        ax.plot([feat[key][0] * s], [feat[key][1] * s], marker="o", ms=9, mfc="none", mew=1.8,
                color=colour, ls="none", zorder=6)

    ax.annotate("Emitterspitze\n(Abb. 2)", xy=(0.005, 0.0), xytext=(0.30 * R, -0.62 * abs(z0)),
                fontsize=8.6, color="0.15",
                arrowprops=dict(arrowstyle="->", color="0.15", lw=1.1))
    ax.annotate("Extraktoröffnung\n(Abb. 3)", xy=(0.5 * par["extractor_aperture_diameter"] * s,
                                                  par["extraction_distance"] * s),
                xytext=(0.34 * R, 1.28 * par["extraction_distance"] * s),
                fontsize=8.6, color="0.15",
                arrowprops=dict(arrowstyle="->", color="0.15", lw=1.1))

    ax.set_xlim(-0.04 * R, 1.04 * R)
    ax.set_ylim(z0 - 0.04 * (z1 - z0), z1 + 0.04 * (z1 - z0))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Radiale Koordinate  r  (mm)")
    ax.set_ylabel("Axiale Koordinate  z  (mm)")
    ax.set_title("(a) vollständiges Randnetz, maßstäblich  (Maßstab in mm)", fontsize=10.5)
    ax.grid(alpha=0.22, zorder=0)
    fig.legend(handles=legend_handles(st), fontsize=7.6, ncol=3, loc="lower left",
               bbox_to_anchor=(0.035, 0.008), framealpha=0.95)

    # -- statistics panel ---------------------------------------------------
    tx.axis("off")
    lines = ["Elementstatistik je Rand   (Längen in µm)", ""]
    lines.append(f"{'Rand':<38}{'n':>5}{'min':>10}{'median':>10}{'max':>10}")
    lines.append("-" * 73)
    for bid in ORDER:
        if bid not in st:
            continue
        d = st[bid]
        lines.append(f"{bid:<38}{d['n']:>5}{d['min']*1e6:>10.3f}"
                     f"{d['med']*1e6:>10.3f}{d['max']*1e6:>10.3f}")
    lines.append("-" * 73)
    lines.append(f"{'gesamt':<38}{int(total['n_elements']):>5}"
                 f"{float(total['len_min_m'])*1e6:>10.3f}"
                 f"{float(total['len_median_m'])*1e6:>10.3f}"
                 f"{float(total['len_max_m'])*1e6:>10.3f}")
    lines += [
        "",
        "Größenfunktion  h(x) = min( min_s [ h_s + G·|x − x_s| ], h_max )",
        "  Quellen s: benannte Merkmale   h = lfs / 32",
        "             übrige Geometrieecken  h = lfs / 8",
        "  lfs = lokale Merkmalsgröße (kürzeste anliegende Kante bzw.",
        "        Abstand zur nächsten nicht anliegenden Kante)",
        "  G = 0.25 (Lipschitz-Schranke)   h_max = Diagonale / 40",
        "  keine Benutzereingabe: kein h_tip, kein h_far.",
        "",
        "Die Fläche bei z = 0 ist die anfängliche ebene Flüssigkeits-",
        "oberfläche – noch kein berechneter Meniskus.",
        "",
        "Elemente auf r = 0 sind Symmetrieelemente, keine Ringelemente;",
        "ihre Rotationsfläche ist exakt null.",
        "",
        "Scharfe Kanten sind bewusst verfeinert. Ein dort später",
        "berechnetes maximales elektrisches Feld ist deshalb KEINE",
        "netzkonvergente Größe.",
    ]
    tx.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=7.9, va="top", ha="left",
            transform=tx.transAxes)

    fig.suptitle("Automatischer achsensymmetrischer Randvernetzer — P1-Beispielparametersatz",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0.135, 1, 0.945))
    fig.savefig(out, dpi=160)
    print("wrote", out)


# ------------------------------------------------------------------ figure 2
def figure_tip(reg, nodes, elems, feat, par, stats, out):
    s = 1e6  # m -> um
    r1 = 0.5 * par["phi_1"] * s
    r2 = 0.5 * par["phi_2"] * s
    st = id_stats(stats)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.8, 6.4),
                                 gridspec_kw={"width_ratios": [1.0, 1.0]})

    for ax, (rlo, rhi, zlo, zhi), title in (
        (a1, (-1.8, 16.5, -13.5, 9.8), "(a) Bohrung, ebene Flüssigkeitsoberfläche und "
                                       "Austrittskante  (Maßstab in µm)"),
        (a2, (4.15, 6.15, -1.2, 1.2), "(b) Ausschnitt um die gepinnte Austrittskante "
                                      "(≈ 9× vergrößert gegenüber (a))"),
    ):
        draw_regions(ax, reg, s)
        big = ax is a1
        draw_mesh(ax, nodes, elems, s, lw=2.0 if big else 2.6,
                  ms=4.2 if big else 6.5, node_edge=not big)
        ax.plot([feat["pinned_contact_edge"][0] * s], [feat["pinned_contact_edge"][1] * s],
                marker="o", ms=15, mfc="none", mew=2.2, color="#0f9d58", ls="none", zorder=6)
        if big:
            ax.plot([feat["emitter_outer_edge"][0] * s], [0], marker="s", ms=11, mfc="none",
                    mew=2.0, color="#8c6d31", ls="none", zorder=6)
            ax.annotate("gepinnte Austrittskante\n(Merkmal, exakt ein Knoten)",
                        (r2 + 0.2, -0.8), xytext=(46, -60), textcoords="offset points",
                        fontsize=8.6, color="#0f9d58",
                        arrowprops=dict(arrowstyle="->", color="#0f9d58", lw=1.1))
            ax.annotate("äußere Stirnkante", (r1 + 0.35, 0.25), xytext=(12, 42),
                        textcoords="offset points", fontsize=8.6, color="#8c6d31",
                        arrowprops=dict(arrowstyle="->", color="#8c6d31", lw=1.1))
            ax.annotate("anfängliche ebene\nFlüssigkeitsoberfläche —\n"
                        "noch kein berechneter Meniskus",
                        (1.2, 0.15), xytext=(-26, 60), textcoords="offset points",
                        fontsize=8.6, color="#0f9d58", ha="left", va="bottom",
                        arrowprops=dict(arrowstyle="->", color="#0f9d58", lw=1.1))
            ax.add_patch(plt.Rectangle((4.15, -1.2), 2.0, 2.4, fill=False, ls=(0, (4, 2)),
                                       lw=1.2, ec="0.25", zorder=7))
            ax.text(5.15, 1.7, "Ausschnitt (b)", ha="center", va="bottom", fontsize=8.2,
                    color="0.25", zorder=7)
        ax.set_xlim(rlo, rhi)
        ax.set_ylim(zlo, zhi)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Radiale Koordinate  r  (µm)")
        ax.set_ylabel("Axiale Koordinate  z  (µm)")
        ax.set_title(title, fontsize=10.0)
        ax.grid(alpha=0.22, zorder=0)

    sub = ", ".join(f"{BSTYLE[b][1]}: n = {st[b]['n']}, "
                    f"min {st[b]['min']*1e6:.3f} / median {st[b]['med']*1e6:.3f} / "
                    f"max {st[b]['max']*1e6:.3f} µm"
                    for b in ("free_surface_reference", "emitter_tip_land", "bore_wall"))
    fig.legend(handles=legend_handles(st, with_counts=False), fontsize=8.2, ncol=5,
               loc="lower center", bbox_to_anchor=(0.5, 0.005), framealpha=0.95)
    fig.suptitle("Randnetz an der Emitterspitze — die beiden Ausschnitte haben "
                 "verschiedene Maßstäbe", fontsize=11.5)
    fig.text(0.5, 0.105, sub, ha="center", fontsize=7.9, color="0.25")
    fig.tight_layout(rect=(0, 0.135, 1, 0.945))
    fig.savefig(out, dpi=160)
    print("wrote", out)


# ------------------------------------------------------------------ figure 3
def figure_aperture(reg, nodes, elems, feat, par, stats, out):
    s = 1e6  # m -> um
    ra = 0.5 * par["extractor_aperture_diameter"] * s
    ze = par["extraction_distance"] * s
    t = par["extractor_thickness"] * s
    st = id_stats(stats)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.8, 6.2),
                                 gridspec_kw={"width_ratios": [1.0, 1.0]})

    windows = (
        (a1, (ra - 130, ra + 210, ze - 130, ze + t + 130),
         "(a) Extraktoröffnung  (Maßstab in µm)"),
        (a2, (ra - 14, ra + 26, ze - 16, ze + 22),
         "(b) vordere Aperturkante (≈ 9× vergrößert gegenüber (a))"),
    )
    for ax, (rlo, rhi, zlo, zhi), title in windows:
        draw_regions(ax, reg, s)
        big = ax is a1
        draw_mesh(ax, nodes, elems, s, lw=2.0 if big else 2.6,
                  ms=4.2 if big else 6.5, node_edge=not big)
        for key in ("extractor_aperture_edge_front", "extractor_aperture_edge_back"):
            ax.plot([feat[key][0] * s], [feat[key][1] * s], marker="o", ms=14, mfc="none",
                    mew=2.2, color="#e6550d", ls="none", zorder=6)
        if big:
            ax.annotate("Aperturkante vorn\n(Merkmal)", (ra, ze), xytext=(30, -44),
                        textcoords="offset points", fontsize=8.6, color="#e6550d",
                        arrowprops=dict(arrowstyle="->", color="#e6550d", lw=1.1))
            ax.annotate("Aperturkante hinten", (ra, ze + t), xytext=(30, 40),
                        textcoords="offset points", fontsize=8.6, color="#e6550d",
                        arrowprops=dict(arrowstyle="->", color="#e6550d", lw=1.1))
            ax.add_patch(plt.Rectangle((ra - 14, ze - 16), 40, 38, fill=False, ls=(0, (4, 2)),
                                       lw=1.2, ec="0.25", zorder=7))
            ax.annotate("Ausschnitt (b)", (ra + 26, ze + 22), xytext=(26, 22),
                        textcoords="offset points", fontsize=8.2, color="0.25",
                        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.0))
        ax.set_xlim(rlo, rhi)
        ax.set_ylim(zlo, zhi)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Radiale Koordinate  r  (µm)")
        ax.set_ylabel("Axiale Koordinate  z  (µm)")
        ax.set_title(title, fontsize=10.0)
        ax.grid(alpha=0.22, zorder=0)

    d = st["extractor_surface"]
    sub = (f"{BSTYLE['extractor_surface'][1]}: n = {d['n']}, min {d['min']*1e6:.3f} / "
           f"median {d['med']*1e6:.3f} / max {d['max']*1e6:.3f} µm   —   "
           "die Verfeinerung folgt der Elektrodendicke, nicht dem Aperturdurchmesser")
    fig.legend(handles=legend_handles(st, with_counts=False), fontsize=8.2, ncol=5,
               loc="lower center", bbox_to_anchor=(0.5, 0.005), framealpha=0.95)
    fig.suptitle("Randnetz an der Extraktoröffnung — die beiden Ausschnitte haben "
                 "verschiedene Maßstäbe", fontsize=11.5)
    fig.text(0.5, 0.105, sub, ha="center", fontsize=7.9, color="0.25")
    fig.tight_layout(rect=(0, 0.135, 1, 0.945))
    fig.savefig(out, dpi=160)
    print("wrote", out)


def main(d):
    reg, nodes, elems, feat, par, stats = load(d)
    figure_overview(reg, nodes, elems, feat, par, stats,
                    os.path.join(d, "fig1_mesh_overview.png"))
    figure_tip(reg, nodes, elems, feat, par, stats,
               os.path.join(d, "fig2_mesh_tip_detail.png"))
    figure_aperture(reg, nodes, elems, feat, par, stats,
                    os.path.join(d, "fig3_mesh_aperture.png"))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: plot_mesh.py <results-directory>")
    main(sys.argv[1])
