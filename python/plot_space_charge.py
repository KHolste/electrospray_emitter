#!/usr/bin/env python3
"""Figures of the P6 space-charge step.

    ./build/es_space_charge results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_space_charge.py results/<dir>

Produces
  <dir>/fig1_space_charge.png   deposited volume charge, the potential change
                                against rho = 0, and the field change
  <dir>/fig2_checks.png         mesh convergence against the manufactured
                                solution, charge conservation, and the approach
                                to a single macroparticle
  <dir>/figures_provenance.txt
"""
import csv
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance as pv

TITLE = "P6: Poisson mit freier Raumladung – die Ladung ist eine VORGESCHRIEBENE Testquelle"

NOT_MODELLED = ("Keine Emission, keine selbstkonsistente PIC-Schleife, keine Teilchenbewegung "
                "(P7), kein Ringmodell. Aus keiner Zahl darf ein emittierter Strom gelesen "
                "werden: P5 ist blockiert.")

STAMP = ""


def rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines))


def meta(d):
    m = {}
    p = os.path.join(d, "meta.txt")
    if os.path.exists(p):
        for ln in open(p, encoding="utf-8"):
            if "=" in ln:
                k, v = ln.rstrip("\n").split("=", 1)
                m[k] = v
    return m


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def col(rs, name):
    return np.array([f(r.get(name)) for r in rs])


def provenance(fig, m, extra="", y=0.006, width=182):
    dirty = pv.DIRTY_MARK in STAMP or STAMP == pv.NO_VERSION
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_space_charge')}   |   "
            f"Status {m.get('status', '?')}   |   "
            f"{m.get('n_particles', '?')} Makropartikel, "
            f"Q = {f(m.get('total_charge_C')):.3g} C   |   "
            f"Erhaltungsfehler {f(m.get('conservation_error')):.1e}")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.0, color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {STAMP}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=178, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


def grid(rs, key):
    r = col(rs, "r_m") * 1e6
    z = col(rs, "z_m") * 1e6
    v = col(rs, key)
    ru = np.unique(r)
    zu = np.unique(z)
    g = np.full((len(zu), len(ru)), np.nan)
    ri = {x: i for i, x in enumerate(ru)}
    zi = {x: i for i, x in enumerate(zu)}
    for a, b, c in zip(r, z, v):
        g[zi[b], ri[a]] = c
    return ru, zu, g


def figure_fields(d, m, out):
    rs = rows(os.path.join(d, "fields.csv"))
    if not rs:
        return None
    fig, axes = plt.subplots(1, 4, figsize=(18.0, 6.4))
    fig.suptitle(TITLE + "\nAbb. 1 – Volumenladung, Potentialänderung und Feldänderung "
                 "gegenüber ρ = 0", fontsize=12.5, y=0.975)

    for ax, key, lab, cmap in (
            (axes[0], "node_charge_C", "deponierte Knotenladung [C]", "viridis"),
            (axes[1], "phi_V", "Potential mit Ladung [V]", "plasma"),
            (axes[2], "dphi_V", "Δφ = φ(ρ) − φ(0) [V]", "magma"),
            (axes[3], "dEz_V_per_m", "ΔE_z [V/m]", "coolwarm")):
        ru, zu, g = grid(rs, key)
        im = ax.pcolormesh(ru, zu, g, shading="nearest", cmap=cmap)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        ax.set_title(lab, fontsize=9.5)
        ax.set_aspect("auto")

    caveat(fig, "Die Ladungsverteilung ist eine VORGESCHRIEBENE Testquelle. P5 ist blockiert, "
                "es gibt keine physikalische Teilchenquelle, und aus keiner Zahl dieses Laufs "
                "darf ein emittierter Strom gelesen werden. Die dritte und vierte Tafel sind "
                "die eigentliche Aussage: sie zeigen ausschließlich das, was die Raumladung "
                "am Feld ändert, weil dasselbe Problem zusätzlich mit ρ = 0 gelöst und "
                "abgezogen wurde.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.150, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


def figure_checks(d, m, out):
    mf = rows(os.path.join(d, "manufactured.csv"))
    cs = rows(os.path.join(d, "conservation.csv"))
    ap = rows(os.path.join(d, "approach.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.4))
    fig.suptitle(TITLE + "\nAbb. 2 – Netzkonvergenz, Ladungserhaltung und die Annäherung "
                 "an ein Makropartikel", fontsize=12.5, y=0.975)

    ax = axes[0]
    if mf:
        h = col(mf, "h_over_R")
        ax.loglog(h, col(mf, "phi_error"), "o-", label="Potential")
        ax.loglog(h, col(mf, "E_error_far"), "s-", label="Feld, r > R/4")
        ax.loglog(h, col(mf, "E_error_near"), "^-", label="Feld, r < R/4")
        ax.loglog(h, col(mf, "phi_error")[0] * (h / h[0]) ** 2.0, "--", color="#888888",
                  label="2. Ordnung")
        ax.loglog(h, col(mf, "E_error_near")[0] * (h / h[0]) ** 1.0, ":", color="#888888",
                  label="1. Ordnung")
    ax.set_xlabel("h / R")
    ax.set_ylabel("relativer Fehler")
    ax.set_title("gegen die hergestellte Lösung\nφ = φ₀(R²−r²)(L²−z²)/(R²L²)", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    ax = axes[1]
    if cs:
        n = col(cs, "n_nodes")
        ax.loglog(n, np.maximum(col(cs, "conservation_error"), 1e-18), "o-", color="#2ca02c",
                  label="Ladungserhaltung der Deposition")
        ax.loglog(n, np.maximum(col(cs, "partition_of_unity_error"), 1e-18), "s-",
                  color="#1f77b4", label="Zerlegung der Eins")
        ax2 = ax.twinx()
        ax2.plot(n, col(cs, "single_particle_peak_phi_V"), "^-", color="#d62728")
        ax2.set_ylabel("Spitzenpotential EINES Makropartikels [V]", color="#d62728")
        ax2.tick_params(axis="y", colors="#d62728")
    ax.set_xscale("log")
    ax.set_xlabel("Knoten")
    ax.set_ylabel("relativer Fehler")
    ax.set_title("Erhaltung (exakt) und Selbstfeld (Netzgröße)", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2, loc="center left")

    ax = axes[2]
    if ap:
        x = col(ap, "d_over_R")
        ax.semilogx(x, col(ap, "phi_V"), "o-", color="#111111", label="Potential [V]")
        ax2 = ax.twinx()
        ax2.semilogx(x, col(ap, "E_magnitude_V_per_m"), "s-", color="#ff7f0e")
        ax2.set_ylabel("|E| [V/m]", color="#ff7f0e")
        ax2.tick_params(axis="y", colors="#ff7f0e")
    ax.set_xlabel("Abstand d / R vom Makropartikel")
    ax.set_ylabel("φ [V]")
    ax.set_title("keine Divergenz beim Annähern", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.6, loc="upper left")

    caveat(fig, "Links steht ein gemessener Befund über die vorhandene Feldrekonstruktion: "
                "sie ist abseits der Achse zweiter, in Achsennähe nur ERSTER Ordnung. Grund "
                "ist die Volumengewichtung, die in Achsensymmetrie einen Faktor 2πr trägt und "
                "die beiden Nachbarzellen eines Knotens im Verhältnis (r+h/2)/(r−h/2) "
                "gewichtet – eine Unsymmetrie der Ordnung h/r. Rechts: ein Ringmodell würde "
                "beim Annähern logarithmisch divergieren; die FEM-Lösung eines "
                "Knotenlastvektors ist stückweise bilinear und bleibt beschränkt. Ihr Betrag "
                "ist dafür eine NETZGRÖSSE (Mitte, rote Kurve) – eine PIC-Schleife darf ein "
                "Teilchen sein eigenes Feld nicht ungefiltert fühlen lassen.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.170, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


def main(d):
    global STAMP
    m = meta(d)
    root = pv.repo_root_of(__file__)
    rel = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/")
    state = pv.repo_state(root, ignore=[rel])
    STAMP = pv.stamp(state, m.get("commit"))
    pv.write_provenance(d, state, m.get("commit"))
    for x in [figure_fields(d, m, os.path.join(d, "fig1_space_charge.png")),
              figure_checks(d, m, os.path.join(d, "fig2_checks.png"))]:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
