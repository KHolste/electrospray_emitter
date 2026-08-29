#!/usr/bin/env python3
"""Figures of the P5 emission audit.  The model is BLOCKED.

    ./build/es_emission_audit results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_emission_audit.py results/<dir>

The figure is a DIMENSIONLESS sensitivity study and a status map.  It is
explicitly not an operating prediction: no absolute current appears anywhere,
because the activation barrier that would set it is not sourced.

Produces
  <dir>/fig1_blocker.png   the status of every path that could produce a number,
                           and the dimensionless sensitivity of the quoted form
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

TITLE = "P5: Ionenemission – Vertrag BLOCKIERT, Modell abgeschaltet"

NOT_MODELLED = ("Keine Betriebsprognose. Kein absoluter Strom, für keine Polarität, bei "
                "keinem Feld. Keine Cone-Jet-Formel als Ionenemission. Kein Strom aus einem "
                "Perfect-Conductor-Feld als Vorhersage.")

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
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_emission_audit')}   |   "
            f"Status {m.get('status', '?')}   |   "
            f"Modell aktiviert: {m.get('model_enabled', '?')}, "
            f"Gleichung geprüft: {m.get('equation_validated', '?')}, "
            f"Barriere belegt: {m.get('barrier_sourced', '?')}")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.0, color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {STAMP}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=178, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


def figure(d, m, out):
    ct = rows(os.path.join(d, "contract.csv"))
    se = rows(os.path.join(d, "sensitivity.csv"))
    bl = rows(os.path.join(d, "barrier_leverage.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.4))
    fig.suptitle(TITLE + "\nAbb. 1 – Statuskarte und dimensionslose Sensitivität "
                 "(KEINE Betriebsprognose)", fontsize=12.5, y=0.975)

    ax = axes[0]
    if ct:
        labels = [f"{r['case']}\n{r['polarity']}" for r in ct]
        y = np.arange(len(ct))
        ax.barh(y, np.ones(len(ct)), color="#d62728", alpha=0.75)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7.2)
        for i, r in enumerate(ct):
            ax.text(0.02, i, r["status"], va="center", fontsize=8, color="white",
                    weight="bold")
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_title("jeder Pfad, der eine Zahl liefern könnte –\n"
                     "und der Status, den er stattdessen meldet", fontsize=9.5)

    ax = axes[1]
    if se:
        for dg in sorted({f(r["dG_eV_at_298K"]) for r in se}):
            sel = [r for r in se if f(r["dG_eV_at_298K"]) == dg]
            ax.semilogy(col(sel, "x"), np.maximum(col(sel, "j_over_j0"), 1e-30),
                        label=f"b = ΔG/kT für ΔG = {dg:.2f} eV")
        ax.axvline(1.0, color="#7a3b00", ls="--", lw=1.2)
        ax.text(1.02, 1e-10, "x = 1: die Barriere ist weg\n(und die Herleitung mit ihr)",
                fontsize=7, color="#7a3b00")
    ax.set_xlabel("x = E / E*,   E* = 4πε₀ ΔG² / e³")
    ax.set_ylabel("j / j₀   (dimensionslos)")
    ax.set_title("die zitierte Ratenform, dimensionslos", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.0)

    ax = axes[2]
    if bl:
        dg = col(bl, "dG_eV")
        ax.semilogy(dg, col(bl, "rate_ratio_to_1.09eV_at_1e9_V_per_m"), "o-", lw=2.0,
                    color="#111111")
        ax.axvline(1.09, color="#7a3b00", ls="--", lw=1.2)
        ax.axvspan(1.0, 1.4, color="#d62728", alpha=0.10)
        ax.text(1.02, ax.get_ylim()[1] * 0.02, "in der Literatur genannte Spanne",
                fontsize=7.5, color="#8a1b1b", rotation=90)
    ax.set_xlabel("Aktivierungsbarriere ΔG [eV]")
    ax.set_ylabel("Rate relativ zu ΔG = 1,09 eV  (bei festem Feld)")
    ax.set_title("warum eine unbelegte Barriere keine Vorhersage trägt", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")

    caveat(fig, "Rechts steht der Blocker als Zahl: über die in der Literatur genannte Spanne "
                "1,0 bis 1,4 eV ändert sich die Rate bei festem Feld um mehr als sechs "
                "Größenordnungen. Für EMI-BF4 gibt es keine belegte Barriere, und die "
                "Ratengleichung selbst wurde an keiner Primärquelle geprüft. Deshalb liefert "
                "kein Pfad eine Zahl – auch nicht, wenn man das Modell einschaltet. Alle "
                "Achsen sind dimensionslos oder Verhältnisse; ein absoluter Strom kommt "
                "nirgends vor.", 0.058)
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
    x = figure(d, m, os.path.join(d, "fig1_blocker.png"))
    if x:
        print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
