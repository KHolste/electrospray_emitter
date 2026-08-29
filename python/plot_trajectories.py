#!/usr/bin/env python3
"""Figures of the P7 particle-transport step.

    ./build/es_trajectories results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_trajectories.py results/<dir>

Produces
  <dir>/fig1_trajectories.png   trajectories in the r-z cut and the impact sites
  <dir>/fig2_balance.png        energy and current balance, and the comparison
                                without / with prescribed space charge
  <dir>/figures_provenance.txt
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance as pv

TITLE = ("P7: Teilchentransport – TRANSPORTANTWORT auf eine vorgeschriebene "
         "Startverteilung, keine Stromvorhersage")

NOT_MODELLED = ("P5 ist blockiert: es gibt keine physikalische Teilchenquelle, und der "
                "gestartete Strom ist eine EINGABE. Kein Tropfenstrahl. Keine "
                "selbstkonsistente Raumladungsschleife – die Teilchen ändern die Ladung "
                "während des Fluges nicht.")

FATE_COLOUR = {"LeftThroughAperture": "#2ca02c", "HitExtractor": "#d62728",
               "HitPolymer": "#ff7f0e", "HitEmitter": "#9467bd", "Flying": "#7f7f7f",
               "LeftDomain": "#1f77b4"}

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
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_trajectories')}   |   "
            f"Status {m.get('status', '?')}   |   "
            f"U = {f(m.get('voltage_V')):.0f} V, Blende {f(m.get('aperture_radius_m')) * 1e6:.2f} µm, "
            f"dt = {f(m.get('dt_s')):.1e} s   |   "
            f"vorgeschriebene Raumladung {f(m.get('prescribed_space_charge_C')):.2e} C")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.0, color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {STAMP}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=178, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


def figure_traj(d, m, out):
    tr = rows(os.path.join(d, "trajectories.csv"))
    im = rows(os.path.join(d, "impacts.csv"))
    if not tr:
        return None
    cases = []
    for r in tr:
        if r["case"] not in cases:
            cases.append(r["case"])
    R = f(m.get("box_radius_m")) * 1e6
    Z = f(m.get("box_length_m")) * 1e6
    ap = f(m.get("aperture_radius_m")) * 1e6
    fig, axes = plt.subplots(1, len(cases), figsize=(5.6 * len(cases), 7.4), squeeze=False)
    fig.suptitle(TITLE + "\nAbb. 1 – Bahnen im r-z-Schnitt und die Auftrefforte",
                 fontsize=12.0, y=0.975)

    by = defaultdict(lambda: defaultdict(list))
    for r in tr:
        by[r["case"]][r["particle"]].append((f(r["r_m"]) * 1e6, f(r["z_m"]) * 1e6))

    for k, c in enumerate(cases):
        ax = axes[0][k]
        fate_of = {r["particle"]: r["fate"] for r in im if r["case"] == c}
        for pid, pts in by[c].items():
            a = np.array(pts)
            ax.plot(a[:, 0], a[:, 1], lw=1.0, color=FATE_COLOUR.get(fate_of.get(pid, ""), "#333"))
        sel = [r for r in im if r["case"] == c]
        for fate in sorted({r["fate"] for r in sel}):
            s2 = [r for r in sel if r["fate"] == fate]
            ax.plot(col(s2, "r_end_m") * 1e6, col(s2, "z_end_m") * 1e6, "o", ms=4.5,
                    color=FATE_COLOUR.get(fate, "#333"), label=f"{fate} ({len(s2)})")
        # the geometry
        ax.plot([0, R], [0, 0], color="#9467bd", lw=3)
        ax.plot([ap, R], [Z, Z], color="#d62728", lw=3)
        ax.plot([0, ap], [Z, Z], color="#2ca02c", lw=3, ls=":")
        ax.plot([R, R], [0, Z], color="#ff7f0e", lw=3)
        ax.set_xlim(-0.02 * R, 1.05 * R)
        ax.set_ylim(-0.03 * Z, 1.05 * Z)
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        ax.set_title(c.replace("_", " "), fontsize=10)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.0, loc="upper right")

    caveat(fig, "Violett ist der Emitter, orange die Polymerwand, rot der Extraktor, grün "
                "gepunktet die Blende. Die Startverteilung ist VORGESCHRIEBEN und der "
                "gestartete Strom eine Eingabe; die Transmission ist deshalb eine "
                "Transportantwort und keine Stromvorhersage. Die Raumladung ist ebenfalls "
                "vorgeschrieben und wird während des Fluges nicht aktualisiert – genau das "
                "macht den Vergleich der ersten beiden Tafeln zu einer sauberen Differenz.",
            0.052)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.140, 1, 0.945])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


def figure_balance(d, m, out):
    ba = rows(os.path.join(d, "balance.csv"))
    im = rows(os.path.join(d, "impacts.csv"))
    if not ba:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.4))
    fig.suptitle(TITLE + "\nAbb. 2 – Energiebilanz, Strombilanz und die Wirkung der "
                 "vorgeschriebenen Raumladung", fontsize=12.0, y=0.975)

    ax = axes[0]
    for c in sorted({r["case"] for r in im}):
        sel = [r for r in im if r["case"] == c]
        ax.plot(col(sel, "q_dV_J") / 1.602176634e-19, col(sel, "energy_gain_eV"), "o", ms=4,
                label=c.replace("_", " "))
    lim = ax.get_xlim()
    ax.plot(lim, lim, "--", color="#888888", label="Energiegewinn = q ΔV")
    ax.set_xlabel("q ΔV [eV]")
    ax.set_ylabel("gemessener Energiegewinn [eV]")
    ax.set_title("Energiegewinn gegen q ΔV", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.2)

    ax = axes[1]
    tags = [r["case"] for r in ba]
    y = np.arange(len(tags))
    ext = col(ba, "extracted_A") * 1e9
    emi = col(ba, "emitter_A") * 1e9
    pol = col(ba, "polymer_A") * 1e9
    exr = col(ba, "extractor_A") * 1e9
    fly = col(ba, "flying_A") * 1e9
    left = np.zeros(len(tags))
    for v, lab, cbar in ((ext, "durch die Blende", "#2ca02c"),
                         (exr, "auf den Extraktor", "#d62728"),
                         (pol, "auf das Polymer", "#ff7f0e"),
                         (emi, "auf den Emitter", "#9467bd"),
                         (fly, "noch fliegend", "#7f7f7f")):
        ax.barh(y, v, left=left, color=cbar, label=lab)
        left = left + v
    ax.plot(col(ba, "launched_A") * 1e9, y, "k|", ms=18, label="gestartet")
    ax.set_yticks(y)
    ax.set_yticklabels([t.replace("_", "\n") for t in tags], fontsize=8)
    ax.set_xlabel("Strom [nA]  – der gestartete ist eine EINGABE")
    ax.set_title("Strombilanz: die Summe ist exakt der gestartete Strom", fontsize=10)
    ax.grid(alpha=0.25, axis="x")
    ax.legend(fontsize=7.0, loc="lower right")

    ax = axes[2]
    ax.barh(y, col(ba, "transmission"), color="#1f77b4")
    for i, r in enumerate(ba):
        ax.text(f(r["transmission"]) + 0.02, i,
                f"{100 * f(r['transmission']):.1f} %\nSchließfehler {f(r['closure_error']):.1e}",
                va="center", fontsize=7.5)
    ax.set_yticks(y)
    ax.set_yticklabels([t.replace("_", "\n") for t in tags], fontsize=8)
    ax.set_xlim(0, 1.25)
    ax.set_xlabel("Transmission durch die Blende")
    ax.set_title("ohne / mit vorgeschriebener Raumladung", fontsize=10)
    ax.grid(alpha=0.25, axis="x")

    caveat(fig, "Links: der gemessene Energiegewinn liegt auf der Geraden q ΔV – das ist die "
                "Energieerhaltung der Bahnintegration, nicht eine Modellaussage. Mitte: "
                "extrahiert + abgefangen + noch fliegend ergibt exakt den gestarteten Strom; "
                "jedes Teilchen hat genau ein Schicksal und noch fliegende werden getrennt "
                "gezählt statt als verloren verbucht. Rechts: die vorgeschriebene positive "
                "Raumladung defokussiert den positiven Strahl. Das ist eine Transportantwort "
                "auf eine gesetzte Ladung, keine Aussage über einen Betriebspunkt.", 0.052)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.150, 1, 0.930])
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
    for x in [figure_traj(d, m, os.path.join(d, "fig1_trajectories.png")),
              figure_balance(d, m, os.path.join(d, "fig2_balance.png"))]:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
