#!/usr/bin/env python3
"""Figures of the P4 kinematic step.

    ./build/es_kinematics results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_kinematics.py results/<dir>

Produces
  <dir>/fig1_advection.png       the advected surface against the exact map
  <dir>/fig2_convergence.png     time convergence and the volume balance
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

TITLE = ("P4: kinematische Randbedingung auf VORGESCHRIEBENEM Feld – "
         "kein dynamischer Meniskus")

NOT_MODELLED = ("Nicht enthalten: dynamischer Meniskus, Kraftbilanz, Strömungslöser mit "
                "freier Oberfläche, Oberflächenladungstransport, Stabilitätsaussage. Keine "
                "Mobilität, keine künstliche Dämpfung. Das Geschwindigkeitsfeld ist "
                "vorgeschrieben; nichts wird dafür gelöst.")

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
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_kinematics')}   |   "
            f"Status {m.get('status', '?')}   |   "
            f"R = {f(m.get('radius_m')) * 1e6:.2f} µm, α = {f(m.get('alpha_per_s')):.3g} 1/s, "
            f"T = {f(m.get('time_s')):.3g} s   |   "
            f"dynamischer Löser: {m.get('dynamic_solver', '?')}")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.0, color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {STAMP}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=178, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


def figure_advection(d, m, out):
    rs = rows(os.path.join(d, "shapes.csv"))
    if not rs:
        return None
    cases = sorted({r["case"] for r in rs})
    fig, axes = plt.subplots(1, len(cases), figsize=(7.0 * len(cases), 6.6), squeeze=False)
    fig.suptitle(TITLE + "\nAbb. 1 – advektierte Oberfläche gegen die exakte "
                 "Lagrange-Abbildung", fontsize=12.5, y=0.975)
    titles = {"dilatation": "u = α x   (div u = 3α, Volumen wächst)",
              "squeeze": "u = (−αr/2, αz)   (div u = 0, Volumen erhalten)"}
    for k, c in enumerate(cases):
        ax = axes[0][k]
        frames = sorted({f(r["t_over_T"]) for r in rs if r["case"] == c})
        for i, t in enumerate(frames):
            sel = [r for r in rs if r["case"] == c and f(r["t_over_T"]) == t]
            if not sel:
                continue
            sel.sort(key=lambda r: int(r["node"]))
            col_t = plt.get_cmap("viridis")(0.15 + 0.7 * i / max(1, len(frames) - 1))
            ax.plot(col(sel, "r_exact_m") * 1e6, col(sel, "z_exact_m") * 1e6, lw=3.4,
                    color=col_t, alpha=0.35)
            ax.plot(col(sel, "r_m") * 1e6, col(sel, "z_m") * 1e6, lw=1.2, color="#111111",
                    label="advektiert" if i == 0 else None)
            if i == 0:
                ax.plot([], [], lw=3.4, color=col_t, alpha=0.35, label="exakte Abbildung")
        ax.set_aspect("equal")
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        ax.set_title(titles.get(c, c) + "\nt/T = 0 … 1, hell nach dunkel", fontsize=9.5)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    caveat(fig, "Das Geschwindigkeitsfeld ist VORGESCHRIEBEN. Dies ist die kinematische "
                "Randbedingung dx/dt·n = u·n, nicht eine Dynamik: es wird keine Kraft "
                "ausgewertet und nichts für u gelöst. Die beiden Felder trennen zwei "
                "Fehlerarten – bei der Dilatation muss das Volumen einer BEKANNTEN Änderung "
                "folgen, beim Squeeze muss es exakt erhalten bleiben, während sich die Form "
                "stark ändert.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.150, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


def figure_convergence(d, m, out):
    cv = rows(os.path.join(d, "convergence.csv"))
    rd = rows(os.path.join(d, "redistribution.csv"))
    if not cv:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.2))
    fig.suptitle(TITLE + "\nAbb. 2 – Zeitkonvergenz, Volumenbilanz und die Netzbewegung",
                 fontsize=12.5, y=0.975)

    ax = axes[0]
    for c in sorted({r["case"] for r in cv}):
        sel = [r for r in cv if r["case"] == c and r["status"] == "Ok"]
        n = col(sel, "steps")
        e = col(sel, "shape_error_over_R")
        ax.loglog(n, np.maximum(e, 1e-18), "o-", label=c)
    n0 = col([r for r in cv if r["status"] == "Ok"], "steps")
    e0 = col([r for r in cv if r["status"] == "Ok"], "shape_error_over_R")
    if len(n0) > 1 and e0[0] > 0:
        ax.loglog(n0[:5], e0[0] * (n0[:5] / n0[0]) ** -4.0, "--", color="#888888",
                  label="4. Ordnung (RK4)")
    ax.set_xlabel("Zeitschritte")
    ax.set_ylabel("max |x − x_exakt| / R")
    ax.set_title("Zeitkonvergenz gegen die exakte Abbildung", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.6)

    ax = axes[1]
    for c in sorted({r["case"] for r in cv}):
        sel = [r for r in cv if r["case"] == c and r["status"] == "Ok"]
        n = col(sel, "steps")
        got = col(sel, "volume_change")
        want = col(sel, "volume_change_exact")
        ax.loglog(n, np.maximum(np.abs(got - want), 1e-18), "s-", label=c)
    ax.set_xlabel("Zeitschritte")
    ax.set_ylabel("|ΔV/V gerechnet − exakt|")
    ax.set_title("Volumenbilanz gegen die EXAKTE Änderung", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.6)

    ax = axes[2]
    if rd:
        n = col(rd, "n_nodes")
        ax.loglog(n, col(rd, "relative_difference"), "o-", color="#d62728",
                  label="Volumenunterschied der beiden Moden")
        good = np.isfinite(col(rd, "relative_difference"))
        if good.sum() > 1:
            v = col(rd, "relative_difference")
            ax.loglog(n[good], v[good][0] * (n[good] / n[good][0]) ** -1.0, "--",
                      color="#888888", label="1. Ordnung")
            ax.loglog(n[good], v[good][0] * (n[good] / n[good][0]) ** -2.0, ":",
                      color="#888888", label="2. Ordnung")
    ax.set_xlabel("Knoten auf der Oberfläche")
    ax.set_ylabel("|V_NormalOnly − V_Lagrange| / V")
    ax.set_title("die tangentiale Umverteilung ist Netzbewegung –\n"
                 "auf einem Polygonzug nur bis auf die Diskretisierung", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.4)

    caveat(fig, "Rechts steht eine gemessene Grenze, keine Beruhigung: die Umverteilung "
                "schiebt Knoten entlang der Fläche und kann die Fläche deshalb im Kontinuum "
                "nicht ändern – auf einem Polygonzug schneidet das Neuabtasten aber Ecken ab. "
                "Die gemessene Ordnung ist rund 1,2 und nicht 2. Ein dynamischer Löser, der "
                "in jedem Schritt umverteilt, braucht deshalb eine gekrümmte Rekonstruktion "
                "oder muss den Verlust bilanzieren.", 0.058)
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
    for x in [figure_advection(d, m, os.path.join(d, "fig1_advection.png")),
              figure_convergence(d, m, os.path.join(d, "fig2_convergence.png"))]:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
