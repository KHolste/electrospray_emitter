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

# One colour per PER-AXIS answer -- never one colour per row.  An earlier
# version coloured a row by its comparability alone, so a quantity that is
# comparable in principle but blocked came out green, and green reads as
# success.  Here green only ever means "this one axis is satisfied".
AX = {"yes": "#2ca02c", "partial": "#ff7f0e", "no": "#d0d0d0", "n/a": "#f5f5f5"}
AXTXT = {"yes": "white", "partial": "white", "no": "#555555", "n/a": "#aaaaaa"}
AXES = [("comparable_geometry", "geometrisch/\nphysikalisch\nvergleichbar"),
        ("implemented", "implementiert"),
        ("converged", "numerisch\nkonvergiert"),
        ("comparable_with_data", "mit Messdaten\nvergleichbar"),
        ("validated", "TATSÄCHLICH\nVALIDIERT")]
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
    fig = plt.figure(figsize=(19.0, 9.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.62, 1.0], height_ratios=[1.0, 1.0],
                          hspace=0.52, wspace=0.24, top=0.830, bottom=0.235,
                          left=0.135, right=0.985)
    fig.suptitle(TITLE + "\nAbb. 1 – Vergleichbarkeit und Validierung sind VERSCHIEDENE "
                         "Fragen; dazu Rotationsreferenz und Importvertrag",
                 fontsize=12.0, y=0.978)

    ax = fig.add_subplot(gs[:, 0])
    if vm:
        n = len(vm)
        y = np.arange(n)
        for j, (key, _) in enumerate(AXES):
            for i, r in enumerate(vm):
                v = r.get(key, "no")
                ax.add_patch(plt.Rectangle((j, i - 0.42), 0.92, 0.84,
                                           facecolor=AX.get(v, "#d0d0d0"),
                                           edgecolor="white", lw=1.2))
                ax.text(j + 0.46, i, v, ha="center", va="center", fontsize=7.0,
                        color=AXTXT.get(v, "#555555"), weight="bold")
        # The blocked column is separate: a quantity can be comparable AND blocked.
        for i, r in enumerate(vm):
            blocked = r.get("blocked") == "yes"
            ax.add_patch(plt.Rectangle((len(AXES), i - 0.42), 0.92, 0.84,
                                       facecolor="#8a1b1b" if blocked else "#f5f5f5",
                                       edgecolor="white", lw=1.2))
            ax.text(len(AXES) + 0.46, i, "BLOCKIERT" if blocked else "–", ha="center",
                    va="center", fontsize=6.4,
                    color="white" if blocked else "#aaaaaa", weight="bold")
        ax.set_xlim(-0.05, len(AXES) + 1.0)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_xticks([j + 0.46 for j in range(len(AXES) + 1)])
        ax.set_xticklabels([lbl for _, lbl in AXES] + ["blockiert"], fontsize=7.6)
        ax.xaxis.set_ticks_position("top")
        ax.set_yticks(y)
        ax.set_yticklabels(["\n".join(textwrap.wrap(r["quantity"], 34)) for r in vm],
                           fontsize=7.2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)
        ax.set_title("Jede Spalte ist eine EIGENE Frage.  Grün heißt nur, dass genau DIESE "
                     "Achse erfüllt ist –\nnie, dass die Größe validiert wäre.",
                     fontsize=9.2, pad=34)

        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=AX["yes"], label="ja"),
                           Patch(color=AX["partial"],
                                 label="nur unter ausgesprochener Einschränkung"),
                           Patch(color=AX["no"], label="nein"),
                           Patch(color=AX["n/a"], label="Achse trifft nicht zu"),
                           Patch(color="#8a1b1b", label="blockiert (Grund in der CSV)")],
                  fontsize=7.4, loc="lower center", ncol=3, frameon=False,
                  bbox_to_anchor=(0.5, -0.090))

        n_val = sum(1 for r in vm if r.get("validated") == "yes")
        n_cmp = sum(1 for r in vm if r.get("comparable_geometry") in ("yes", "partial"))
        n_blk = sum(1 for r in vm if r.get("blocked") == "yes")
        fig.text(0.40, 0.170,
                 f"{n_cmp} von {n} Größen sind vergleichbar.  {n_blk} sind blockiert.  "
                 f"TATSÄCHLICH VALIDIERT: {n_val}.",
                 fontsize=10.0, color="#8a1b1b", weight="bold", ha="center")

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
        # THREE outcomes, not two.  Collapsing the middle one into "abgelehnt"
        # is what the earlier contract did, and it threw real measurements away.
        def col_of(st):
            if st == "Ok":
                return "#2ca02c"
            if st == "OkUncertaintyNotReported":
                return "#ff7f0e"
            return "#d62728"
        y = np.arange(len(ic))
        ax.barh(y, np.ones(len(ic)), color=[col_of(r["status"]) for r in ic])
        ax.set_yticks(y)
        ax.set_yticklabels(["\n".join(textwrap.wrap(r["case"].replace("_", " "), 26))
                            for r in ic], fontsize=7.0)
        for i, r in enumerate(ic):
            ax.text(0.02, i, r["status"], va="center", fontsize=6.8, color="white",
                    weight="bold")
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.invert_yaxis()
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color="#2ca02c", label="vollständig – quantitativ verwendbar"),
                           Patch(color="#ff7f0e",
                                 label="Unsicherheit NICHT berichtet – archiviert,\n"
                                       "qualitativ darstellbar, trägt keine Zahl"),
                           Patch(color="#d62728", label="harter Fehler – Satz abgelehnt")],
                  fontsize=6.6, loc="upper left", frameon=True,
                  bbox_to_anchor=(0.0, -0.06))
        ax.set_title("Importvertrag: Einheit, Fundstelle und Geometrieart sind HARTE\n"
                     "Bedingungen – eine nicht berichtete Unsicherheit ist es nicht",
                     fontsize=9.5)

    caveat(fig, "LINKS: eine frühere Fassung färbte jede Zeile nach ihrer Vergleichbarkeit "
                "allein. Damit erschien der Gesamtstrom GRÜN – er ist direkt vergleichbar "
                "und wird von diesem Projekt überhaupt nicht gerechnet, weil P5 blockiert "
                "ist. Vergleichbarkeit und Validierung sind verschiedene Fragen und werden "
                "jetzt getrennt beantwortet. Die Spalte „TATSÄCHLICH VALIDIERT“ ist überall "
                "nein: es sind keine Messdaten importiert, und der Test prüft, dass die "
                "Matrix das sagt statt es durch eine Farbe zu verdecken. Die Invariante "
                "„validiert verlangt implementiert UND konvergiert UND mit Daten "
                "vergleichbar UND nicht blockiert“ steht im Code, nicht in der Bildunterschrift. "
                "Vier Größen sind grundsätzlich nicht vergleichbar – azimutale Asymmetrie, "
                "Versatz, Neigung, Array-Übersprechen sind genau das, was ein "
                "achsensymmetrisches Modell nicht hat; sie stehen benannt da, statt zu "
                "fehlen. RECHTS OBEN: der Unterschied ist Summationsrundung, kein "
                "Modellunterschied. RECHTS UNTEN: keiner dieser Punkte ist eine Messung. "
                "Fehlende Einheit, Fundstelle oder Geometrieart sind harte Fehler und lehnen "
                "den Satz als GANZES ab; eine in der Publikation nicht angegebene "
                "Unsicherheit ist dagegen ein eigener Zustand – der Punkt wird archiviert "
                "und darf qualitativ gezeigt werden, trägt aber keine quantitative "
                "Validierung.", 0.036)
    provenance(fig, m, NOT_MODELLED, y=0.004)
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
