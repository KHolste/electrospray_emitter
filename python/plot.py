#!/usr/bin/env python3
"""Plot the CSV output of es_meniscus / es_field.

    python python/plot.py out_capillary

Reads <prefix>_branch.csv, <prefix>_shape_*.csv, <prefix>_surface.csv and
<prefix>_mesh.csv, whichever exist.  Pure post-processing -- nothing here
computes physics.
"""
import csv
import os
import sys

import matplotlib.pyplot as plt


def read_csv(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
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


def main(prefix):
    branch = read_csv(prefix + "_branch.csv")
    surf = read_csv(prefix + "_surface.csv")
    mesh = read_csv(prefix + "_mesh.csv")
    shapes = [(lab, read_csv(f"{prefix}_shape_{lab}.csv")) for lab in ("onset", "last")]
    shapes = [(l, s) for l, s in shapes if s]

    npanel = sum(x is not None for x in (branch, surf)) + (1 if (mesh or shapes) else 0)
    if npanel == 0:
        sys.exit(f"no CSV files found for prefix '{prefix}'")

    fig, axes = plt.subplots(1, npanel, figsize=(5.2 * npanel, 4.4))
    if npanel == 1:
        axes = [axes]
    ax = iter(axes)

    # --- equilibrium branch --------------------------------------------------
    if branch:
        a = next(ax)
        h = [x * 1e6 for x in branch["height"]]
        u = branch["voltage"]
        conv = branch["converged"]
        a.plot(h, u, "-", color="0.6", lw=1)
        a.plot([hi for hi, c in zip(h, conv) if c], [ui for ui, c in zip(u, conv) if c],
               "o", ms=4, label="konvergiert")
        bad = [(hi, ui) for hi, ui, c in zip(h, u, conv) if not c]
        if bad:
            a.plot(*zip(*bad), "x", ms=6, color="crimson", label="nicht konvergiert")
        imax = max(range(len(u)), key=lambda i: u[i] if conv[i] else -1e30)
        a.axhline(u[imax], ls="--", lw=0.8, color="crimson")
        a.annotate(f"Onset {u[imax]:.0f} V", (h[imax], u[imax]),
                   textcoords="offset points", xytext=(6, 8), color="crimson")
        a.set_xlabel("Apex-Höhe h  [µm]")
        a.set_ylabel("Spannung U  [V]")
        a.set_title("Gleichgewichtsast U(h)\nMaximum = Onset, danach instabil")
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

    # --- geometry and meniscus shapes ---------------------------------------
    if mesh or shapes:
        a = next(ax)
        if mesh:
            for r0, z0, r1, z1, tag in zip(mesh["r_a"], mesh["z_a"], mesh["r_b"], mesh["z_b"],
                                           mesh["tag"]):
                c = {"emitter": "0.25", "meniscus": "tab:blue",
                     "extractor": "tab:orange"}.get(tag, "0.7")
                a.plot([r0 * 1e6, r1 * 1e6], [z0 * 1e6, z1 * 1e6], "-", color=c, lw=1.2)
        for lab, s in shapes:
            a.plot([r * 1e6 for r in s["r"]], [z * 1e6 for z in s["z"]], "-", lw=2,
                   label={"onset": "am Onset", "last": "letzter Ast-Punkt"}[lab])
        a.set_xlabel("r  [µm]")
        a.set_ylabel("z  [µm]")
        a.set_title("Geometrie und Meniskusform (Zoom auf den Emitter)")
        # An extractor is typically two orders of magnitude larger than the
        # meniscus, and the emitter shank is long.  Frame the tip, or the part
        # that matters is a single pixel.
        pts = []
        for _lab, sh in shapes:
            pts += list(zip(sh["r"], sh["z"]))
        if mesh:
            pts += [(r, z) for r, z, t in zip(mesh["r_mid"], mesh["z_mid"], mesh["tag"])
                    if t == "meniscus"]
        if not pts and mesh:
            pts = [(r, z) for r, z, t in zip(mesh["r_mid"], mesh["z_mid"], mesh["tag"])
                   if t == "emitter"]
        if pts:
            rr = [p[0] * 1e6 for p in pts]
            zz = [p[1] * 1e6 for p in pts]
            half = 2.2 * max(max(rr), max(zz) - min(zz), 1e-3)
            zc = 0.5 * (max(zz) + min(zz))
            a.set_xlim(-0.08 * half, half)
            a.set_ylim(zc - 0.7 * half, zc + 0.7 * half)
        a.set_aspect("equal", adjustable="box")
        if shapes:
            a.legend(fontsize=8)
        a.grid(alpha=0.3)

    # --- surface field -------------------------------------------------------
    if surf:
        a = next(ax)
        for tag, col in (("meniscus", "tab:blue"), ("emitter", "0.35"),
                         ("extractor", "tab:orange")):
            rr = [r * 1e6 for r, t in zip(surf["r"], surf["tag"]) if t == tag]
            ee = [abs(e) for e, t in zip(surf["En"], surf["tag"]) if t == tag]
            if rr:
                a.plot(rr, ee, ".", ms=3, color=col, label=tag)
        a.set_yscale("log")
        a.set_xlabel("r  [µm]")
        a.set_ylabel(r"$|E_n|$  [V/m]")
        a.set_title("Normalfeld an der Oberfläche")
        a.legend(fontsize=8)
        a.grid(alpha=0.3, which="both")

    fig.tight_layout()
    out = prefix + "_plot.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    if os.environ.get("ES_SHOW"):
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "out_capillary")
