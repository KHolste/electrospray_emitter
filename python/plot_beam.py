#!/usr/bin/env python3
"""Plot the CSV output of es_beam.

    python python/plot_beam.py out_capillary

Reads <prefix>_paths.csv, <prefix>_rays.csv and <prefix>_mesh.csv.
"""
import csv
import glob
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt


def read_csv(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as fh:
        # Output files carry a '#'-prefixed provenance header (application,
        # state, voltage); skip it before the CSV starts.
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    if not rows:
        return None
    cols = {k: [] for k in rows[0]}
    for r in rows:
        for k, v in r.items():
            try:
                cols[k].append(float(v))
            except (TypeError, ValueError):
                cols[k].append(v)
    return cols


def first(pattern):
    hits = sorted(glob.glob(pattern))
    return read_csv(hits[0]) if hits else None


def main(prefix):
    paths = first(f"{prefix}_beam_*_paths.csv")
    rays = first(f"{prefix}_beam_*_rays.csv")
    mesh = first(f"{prefix}_beam_*_mesh.csv")
    if not rays:
        sys.exit(f"no {prefix}_beam_*_rays.csv found")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # --- trajectories over the geometry -------------------------------------
    if mesh:
        for r0, z0, r1, z1, tag in zip(mesh["r_a"], mesh["z_a"], mesh["r_b"], mesh["z_b"],
                                       mesh["tag"]):
            c = {"emitter": "0.25", "meniscus": "tab:blue",
                 "extractor": "tab:orange"}.get(tag, "0.7")
            a1.plot([r0 * 1e3, r1 * 1e3], [z0 * 1e3, z1 * 1e3], "-", color=c, lw=1.4)
    if paths:
        tracks = defaultdict(list)
        for i, r, z in zip(paths["ray"], paths["r"], paths["z"]):
            tracks[int(i)].append((r * 1e3, z * 1e3))
        status = rays["status"]
        for i, pts in tracks.items():
            col = {"transmitted": "tab:green", "intercepted": "crimson"}.get(
                status[i] if i < len(status) else "", "0.5")
            a1.plot(*zip(*pts), "-", lw=0.8, alpha=0.8, color=col)
    a1.set_xlabel("r  [mm]")
    a1.set_ylabel("z  [mm]")
    a1.set_title("Trajektorien\ngrün = transmittiert, rot = abgefangen")
    a1.grid(alpha=0.3)

    # --- angular current distribution ---------------------------------------
    ang = [a for a, s in zip(rays["angle_deg"], rays["status"]) if s == "transmitted"]
    cur = [c for c, s in zip(rays["current"], rays["status"]) if s == "transmitted"]
    if ang:
        order = sorted(range(len(ang)), key=lambda i: abs(ang[i]))
        x, acc, tot = [], [], sum(cur)
        run = 0.0
        for i in order:
            run += cur[i]
            x.append(abs(ang[i]))
            acc.append(100.0 * run / tot)
        a2.plot(x, acc, "-o", ms=3)
        for frac, lab in ((50, "50 %"), (95, "95 %")):
            a2.axhline(frac, ls="--", lw=0.8, color="0.5")
            j = next((k for k, v in enumerate(acc) if v >= frac), None)
            if j is not None:
                a2.annotate(f"{lab}: {x[j]:.1f}°", (x[j], frac),
                            textcoords="offset points", xytext=(6, -12), fontsize=8)
        a2.set_xlabel("Halbwinkel  [deg]")
        a2.set_ylabel("eingeschlossener Strom  [%]")
        a2.set_title("Winkelverteilung des transmittierten Stroms")
        a2.grid(alpha=0.3)
    else:
        a2.text(0.5, 0.5, "kein transmittierter Strom", ha="center", va="center")

    fig.tight_layout()
    out = prefix + "_beam.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    if os.environ.get("ES_SHOW"):
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "out")
