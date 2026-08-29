#!/usr/bin/env python3
"""Figures for the branch-ambiguity check.

    build/branch_figure.exe results/<run>
    python python/plot_branch.py results/<run>

Reads only files written by tools/branch_figure.cpp in that directory, so the
figures are reproducible from source.  Nothing here computes physics.

Two figures are produced:
  <dir>/fig1_branch_U_vs_h.png      U(h) with the target voltage, both
                                    crossings, the fold candidate and the end
                                    of each investigated range
  <dir>/fig2_meniscus_profiles.png  both meniscus profiles at the target
                                    voltage, with emitter and extractor
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def read_csv(path):
    """CSV with a leading '#'-prefixed provenance header."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    cols = {k: [] for k in rows[0]}
    for r in rows:
        for k, v in r.items():
            try:
                cols[k].append(float(v))
            except (TypeError, ValueError):
                cols[k].append(v)
    return cols


def markers(path):
    out = {}
    for row in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        out[row["name"]] = {k: (float(v) if v else None) for k, v in row.items()
                            if k != "name"}
    return out


# ---------------------------------------------------------------- figure 1
def figure_branch(d, mk, out):
    full = read_csv(os.path.join(d, "branch_full.csv"))
    conv = [i for i, st in enumerate(full["status"]) if st == "converged"]
    bad = [i for i, st in enumerate(full["status"]) if st != "converged"]
    h = [full["h_over_rc"][i] for i in conv]
    u = [full["voltage_V"][i] for i in conv]
    u_target = mk["target_voltage"]["voltage_V"]
    f = mk["fold_candidate"]
    h_stop = min([full["h_over_rc"][i] for i in bad], default=max(h))

    ranges = (("range_end_before_fold", (0, (4, 3)), "Bereichsende vor dem Umkehrpunkt"),
              ("range_end_past_fold", (0, (1, 2)), "Bereichsende hinter dem Umkehrpunkt"),
              ("range_end_full", (0, (6, 2, 1, 2)), "angefordertes Bereichsende (h_max)"))

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.8, 8.8),
                                 gridspec_kw={"height_ratios": [1.0, 1.35]})

    # -- (a) the whole requested range -------------------------------------
    a1.plot(h, u, "-o", ms=4, lw=1.5, color="0.25")
    a1.axhline(u_target, ls="--", lw=1.2, color="tab:red")
    for key, style, lab in ranges:
        x = mk[key]["h_over_rc"]
        a1.axvline(x, ls=style, lw=1.3, color="tab:purple", alpha=0.9)
        a1.text(x, 0.03, f" {lab}", rotation=90, va="bottom", ha="right", fontsize=7.4,
                color="tab:purple", transform=a1.get_xaxis_transform())
    a1.axvline(h_stop, ls="-", lw=1.3, color="crimson", alpha=0.8)
    a1.text(h_stop, 0.97, f" Fortsetzung gestoppt bei h/r_c = {h_stop:.2f}", rotation=90,
            va="top", ha="left", fontsize=7.6, color="crimson",
            transform=a1.get_xaxis_transform())
    a1.set_xlim(0, mk["range_end_full"]["h_over_rc"] * 1.05)
    a1.set_ylim(min(u) - 40, max(u) + 55)
    a1.set_xlabel("Bezogene Apexhöhe  h / r_c  (–)")
    a1.set_ylabel("Spannung  U  (V)")
    a1.set_title("(a) angeforderter Bereich gegen tatsächlich verfolgten Bereich", fontsize=10)
    a1.grid(alpha=0.3)

    # -- (b) detail around the crossings -----------------------------------
    a2.plot(h, u, "-o", ms=5.5, lw=1.7, color="0.25", label="statischer Ast (konvergiert)")
    a2.axhline(u_target, ls="--", lw=1.3, color="tab:red")
    a2.text(0.015, u_target, f" Zielspannung {u_target:.1f} V", color="tab:red", va="bottom",
            ha="left", fontsize=9, transform=a2.get_yaxis_transform())

    for key, col, lab, dx, dy, ha in (
            ("crossing_lower_height", "tab:blue", "LowerHeight", 14, -52, "left"),
            ("crossing_upper_height", "tab:green", "UpperHeight", 10, -62, "left")):
        m = mk[key]
        a2.plot([m["h_over_rc"]], [m["voltage_V"]], "o", ms=13, mfc="none", mew=2.4,
                color=col, label=f"Schnittpunkt {lab}")
        a2.annotate(f"{lab}\nh/r_c = {m['h_over_rc']:.3f}\n"
                    f"R_apex/r_c = {m['apex_radius_m'] / 1e-5:.3f}\n"
                    f"E_apex = {m['apex_field_Vpm']:.2e} V/m",
                    (m["h_over_rc"], m["voltage_V"]), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, fontsize=8.4, color=col)

    a2.plot([f["h_over_rc"]], [f["voltage_V"]], "*", ms=18, color="tab:orange",
            label="Kandidat statischer Umkehrpunkt", zorder=5)
    a2.annotate(f"Kandidat Umkehrpunkt {f['voltage_V']:.1f} V",
                (f["h_over_rc"], f["voltage_V"]), textcoords="offset points",
                xytext=(-8, 12), ha="right", fontsize=8.6, color="tab:orange")

    for key, style, lab in ranges[:2]:
        x = mk[key]["h_over_rc"]
        a2.axvline(x, ls=style, lw=1.3, color="tab:purple", alpha=0.9)
        a2.text(x, 0.02, f" {lab}", rotation=90, va="bottom", ha="right", fontsize=7.6,
                color="tab:purple", transform=a2.get_xaxis_transform())

    a2.set_xlim(0.05, max(h) + 0.08)
    a2.set_ylim(min(u) - 25, max(u) + 30)
    a2.set_xlabel("Bezogene Apexhöhe  h / r_c  (–)")
    a2.set_ylabel("Emitter-Extraktor-Spannung  U  (V)")
    a2.set_title("(b) Detail: eine Zielspannung, zwei Lösungen", fontsize=10)
    a2.grid(alpha=0.3)
    a2.legend(fontsize=8.4, loc="lower left")

    fig.suptitle("Statischer Meniskusast — LowerHeight/UpperHeight bezeichnen die "
                 "Apexhöhe, nicht die Stabilität", fontsize=10.5)
    fig.text(0.5, 0.016,
             "Der Umkehrpunkt ist ein Kandidat aus diskreten Astpunkten.\n"
             "Daraus folgt weder dynamische Stabilität noch ein Emissionsbeginn "
             "noch ein Cone-Jet-Übergang.",
             ha="center", va="bottom", fontsize=8.2, color="0.35")
    fig.tight_layout(rect=(0, 0.062, 1, 0.965))
    fig.savefig(out, dpi=160)
    print("wrote", out)


# ---------------------------------------------------------------- figure 2
def _draw_geometry(ax, mesh, scale):
    colours = {"emitter": "0.25", "extractor": "tab:orange", "meniscus": None}
    for r0, z0, r1, z1, tag in zip(mesh["r_a"], mesh["z_a"], mesh["r_b"], mesh["z_b"],
                                   mesh["tag"]):
        c = colours.get(tag, "0.7")
        if c is None:
            continue                       # meniscus is drawn from the shape files
        ax.plot([r0 * scale, r1 * scale], [z0 * scale, z1 * scale], "-", color=c, lw=1.5)


def figure_profiles(d, mk, out):
    mesh = read_csv(os.path.join(d, "mesh_electrodes.csv"))
    lo = read_csv(os.path.join(d, "shape_lower_height.csv"))
    hi = read_csv(os.path.join(d, "shape_upper_height.csv"))
    u_target = mk["target_voltage"]["voltage_V"]
    s = 1e6                                # m -> um

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 5.2))

    # -- overview, emitter and extractor -----------------------------------
    _draw_geometry(a1, mesh, s)
    a1.plot([r * s for r in lo["r"]], [z * s for z in lo["z"]], "-", lw=2, color="tab:blue")
    a1.plot([r * s for r in hi["r"]], [z * s for z in hi["z"]], "-", lw=2, color="tab:green")
    a1.set_aspect("equal", adjustable="box")
    a1.set_xlim(-20, 700)
    a1.set_ylim(-120, 700)
    zoom = 22.0
    a1.add_patch(Rectangle((-2, -8), zoom, zoom, fill=False, ec="crimson", lw=1.4))
    a1.annotate("Ausschnitt rechts", (zoom, 8), xytext=(120, 90),
                textcoords="data", fontsize=8.5, color="crimson",
                arrowprops=dict(arrowstyle="->", color="crimson", lw=1.1))
    a1.set_xlabel("Radiale Koordinate  r  (µm)")
    a1.set_ylabel("Axiale Koordinate  z  (µm)")
    a1.set_title("Übersicht: Emitter (grau) und Extraktor (orange)", fontsize=10.5)
    a1.grid(alpha=0.25)

    # -- zoom on the two menisci -------------------------------------------
    _draw_geometry(a2, mesh, s)
    a2.plot([r * s for r in lo["r"]], [z * s for z in lo["z"]], "-", lw=2.6,
            color="tab:blue", label="LowerHeight")
    a2.plot([r * s for r in hi["r"]], [z * s for z in hi["z"]], "-", lw=2.6,
            color="tab:green", label="UpperHeight")
    a2.plot([0], [mk["crossing_lower_height"]["height_m"] * s], "o", ms=5, color="tab:blue")
    a2.plot([0], [mk["crossing_upper_height"]["height_m"] * s], "o", ms=5, color="tab:green")
    a2.set_aspect("equal", adjustable="box")
    a2.set_xlim(-2, 20)
    a2.set_ylim(-8, 12)
    a2.set_xlabel("Radiale Koordinate  r  (µm)")
    a2.set_ylabel("Axiale Koordinate  z  (µm)")
    a2.set_title(f"Beide Meniskusformen bei U = {u_target:.1f} V", fontsize=10.5)
    a2.legend(fontsize=9, loc="upper right")
    a2.grid(alpha=0.25)

    fig.suptitle("Zwei geometrisch verschiedene Menisken zu derselben Spannung — "
                 "LowerHeight/UpperHeight bezeichnen die Apexhöhe, nicht die Stabilität",
                 fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=160)
    print("wrote", out)


def main(d):
    mk = markers(os.path.join(d, "markers.csv"))
    figure_branch(d, mk, os.path.join(d, "fig1_branch_U_vs_h.png"))
    figure_profiles(d, mk, os.path.join(d, "fig2_meniscus_profiles.png"))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: plot_branch.py <results-directory>")
    main(sys.argv[1])
