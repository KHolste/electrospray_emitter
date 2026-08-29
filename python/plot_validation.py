#!/usr/bin/env python3
"""Figure of the P9 validation scheme.

    ./build/es_validation results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_validation.py results/<dir>

Produces
  <dir>/fig1_validation.png   the validation matrix, the revolution reference
                              and the import contract
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

TITLE = ("P9: achsensymmetrisch gegen 3D – Validierungsschema, Rotationsreferenz, "
         "Importvertrag")

NOT_MODELLED = ("Es gibt kein 3D-Netz, keinen 3D-Löser und kein 3D-Ergebnis. Die "
                "Rotationsreferenz löst nichts Neues; sie prüft eine Wichtung. Keiner der "
                "Importpunkte ist eine Messung – sie zeigen die Pflichtfelder.")

COL = {"Direct": "#2ca02c", "AfterStatedReduction": "#ff7f0e", "NotComparable": "#d62728"}
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


def provenance(fig, m, extra="", y=0.008, width=182):
    dirty = pv.DIRTY_MARK in STAMP or STAMP == pv.NO_VERSION
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_validation')}   |   "
            f"Status {m.get('status', '?')}   |   "
            f"3D-Löser: {m.get('three_dimensional_solver', '?')}   |   "
            f"Rotationsreferenz: {m.get('revolution_reference', '?')}")
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
    vm = rows(os.path.join(d, "validation_matrix.csv"))
    rv = rows(os.path.join(d, "revolution.csv"))
    ic = rows(os.path.join(d, "import_contract.csv"))
    fig = plt.figure(figsize=(17.5, 8.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.0, 1.0],
                          hspace=0.42, wspace=0.30)
    fig.suptitle(TITLE + "\nAbb. 1 – Vergleichsschema, Rotationsreferenz und Importvertrag",
                 fontsize=12.0, y=0.975)

    ax = fig.add_subplot(gs[:, 0])
    if vm:
        y = np.arange(len(vm))
        ax.barh(y, np.ones(len(vm)), color=[COL.get(r["comparability"], "#888") for r in vm])
        ax.set_yticks(y)
        ax.set_yticklabels(["\n".join(textwrap.wrap(r["quantity"], 34)) for r in vm],
                           fontsize=7.2)
        for i, r in enumerate(vm):
            ax.text(0.015, i, f"{r['comparability']}   –   {r['status']}", va="center",
                    fontsize=6.8, color="white", weight="bold")
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.invert_yaxis()
        ax.set_title("Validierungsmatrix: was sich zwischen einer achsensymmetrischen\n"
                     "Rechnung und einem 3D-Gerät vergleichen lässt – und was nicht",
                     fontsize=10)
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=COL["Direct"], label="direkt vergleichbar"),
                           Patch(color=COL["AfterStatedReduction"],
                                 label="erst nach ausgesprochener Reduktion"),
                           Patch(color=COL["NotComparable"],
                                 label="grundsätzlich nicht vergleichbar")],
                  fontsize=7.4, loc="lower right")

    ax = fig.add_subplot(gs[0, 1])
    if rv:
        for qn in sorted({r["quantity"] for r in rv}):
            sel = [r for r in rv if r["quantity"] == qn]
            ax.loglog(col(sel, "n_azimuth"), np.maximum(col(sel, "relative_difference"), 1e-18),
                      "o-", label=qn)
        ax.axhline(1e-11, color="#7a3b00", ls="--", lw=1.2, label="Schranke des Tests")
    ax.set_xlabel("Zahl der Azimute in der expliziten 3D-Quadratur")
    ax.set_ylabel("|3D − 2πr-Form| / |2πr-Form|")
    ax.set_title("Rotationsreferenz: die 2πr-Wichtung IST das 3D-Integral", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    ax = fig.add_subplot(gs[1, 1])
    if ic:
        y = np.arange(len(ic))
        ok = np.array([r["status"] == "Ok" for r in ic])
        ax.barh(y, np.ones(len(ic)), color=np.where(ok, "#2ca02c", "#d62728"))
        ax.set_yticks(y)
        ax.set_yticklabels([r["case"].replace("_", " ") for r in ic], fontsize=7.6)
        for i, r in enumerate(ic):
            ax.text(0.02, i, r["status"], va="center", fontsize=7.2, color="white",
                    weight="bold")
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.invert_yaxis()
        ax.set_title("Importvertrag: Einheit, Unsicherheit MIT TYP, Fundstelle\n"
                     "und Geometrieart – jede fehlende schlägt geschlossen fehl",
                     fontsize=9.5)

    caveat(fig, "Links: vier Größen sind grundsätzlich nicht vergleichbar – azimutale "
                "Asymmetrie, Versatz, Neigung und Array-Übersprechen sind genau das, was ein "
                "achsensymmetrisches Modell nicht hat. Sie stehen benannt in der Matrix, "
                "statt zu fehlen. Rechts oben: der Unterschied ist die Summationsrundung, "
                "nicht ein Modellunterschied; die Azimutzahl ändert nichts, weil der "
                "Integrand azimutunabhängig ist. Rechts unten: keiner dieser Punkte ist eine "
                "Messung – ein Satz mit einem unvollständigen Punkt wird als GANZES "
                "abgelehnt.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.savefig(out, dpi=145, bbox_inches=None)
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
    x = figure(d, m, os.path.join(d, "fig1_validation.png"))
    if x:
        print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
