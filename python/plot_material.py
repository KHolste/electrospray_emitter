#!/usr/bin/env python3
"""Figures of the P2 material-data audit.

    ./build/es_material results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_material.py results/<dir>

Reads only the CSVs written by apps/es_material.cpp.  Nothing here fetches,
fits or averages.

Produces
  <dir>/fig1_properties_vs_T.png   every measured point of every source against
                                   temperature, the selected source highlighted
  <dir>/fig2_scatter.png           the literature scatter at 298,15 K, source by
                                   source, with the provenance each one states
  <dir>/fig3_impact.png            what the sourced gamma changes about the
                                   numbers P3a/P3b reported
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

TITLE = "P2: Stoffdaten für EMI-BF4 – Quellen, Streuung und was fehlt"

NOT_MODELLED = ("Es wird NICHT gemittelt: widersprechende Quellen stehen nebeneinander. "
                "Ein leeres Provenienzfeld heißt „die Quelle hat es nicht angegeben“, nie "
                "null. Frequenzaufgelöste und Nicht-Umgebungsdruck-Punkte tragen keinen "
                "Gleichstrom- bzw. Umgebungsdruckwert.")

STAMP = ""
T_REF = 298.15

LABEL = {
    "surface_tension": ("Oberflächenspannung γ", "N/m", 1e3, "mN/m"),
    "density": ("Dichte ρ", "kg/m³", 1.0, "kg/m³"),
    "dynamic_viscosity": ("dynamische Viskosität µ", "Pa·s", 1e3, "mPa·s"),
    "kinematic_viscosity": ("kinematische Viskosität ν", "m²/s", 1e6, "mm²/s"),
    "electrical_conductivity": ("elektrische Leitfähigkeit K", "S/m", 1.0, "S/m"),
    "relative_permittivity": ("relative Permittivität ε_r", "-", 1.0, "-"),
}


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


def short_label(r):
    """First author and year, or the data-sheet name.  A circle marks a source
    that does NOT state both purity and water content."""
    ref = r.get("reference", "")
    if ref.startswith("IoLiTec"):
        short, yr = "IoLiTec IL-0006", "2012"
    else:
        short = ref.split(";")[0].split(",")[0].strip()
        yr = ""
        for tok in ref.replace("(", " (").split():
            if tok.startswith("(") and len(tok) >= 5 and tok[1:5].isdigit():
                yr = tok[1:5]
                break
    mark = "" if (r.get("states_purity") == "yes" and r.get("states_water") == "yes") else " \u25cb"
    return f"{short} {yr}{mark}"


def provenance(fig, m, extra="", y=0.008, width=182):
    dirty = pv.DIRTY_MARK in STAMP or STAMP == pv.NO_VERSION
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_material')}   |   "
            f"{m.get('substance', '?')}, CAS {m.get('cas', '?')}   |   "
            f"Quellen: {m.get('database', '?')}")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.0, color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {STAMP}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=178, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


# ==========================================================================
def figure_vs_T(d, m, out):
    pts = rows(os.path.join(d, "points.csv"))
    summ = rows(os.path.join(d, "summary.csv"))
    if not pts:
        return None
    props = [p for p in LABEL if any(r["property"] == p for r in pts)]
    n = len(props)
    ncol = 3
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.5 * ncol, 4.9 * nrow), squeeze=False)
    fig.suptitle(TITLE + "\nAbb. 1 – jeder gemessene Punkt jeder Quelle gegen die "
                 "Temperatur", fontsize=12.5, y=0.985)

    for k, prop in enumerate(props):
        ax = axes[k // ncol][k % ncol]
        name, unit, scale, ulabel = LABEL[prop]
        sel = None
        for r in pts:
            if r["property"] != prop:
                continue
            T, v = f(r["T_K"]), f(r["value"]) * scale
            amb = r["ambient"] == "yes"
            fr = f(r["frequency_Hz"]) > 0
            if r["selected"] == "yes":
                ax.plot(T, v, "o", ms=6.5, mfc="#d62728", mec="#111111", mew=0.8, zorder=5)
                sel = r["database_id"]
            elif fr:
                ax.plot(T, v, "v", ms=4.5, color="#9467bd", alpha=0.85, zorder=3)
            elif not amb:
                ax.plot(T, v, "x", ms=4.5, color="#7f7f7f", alpha=0.7, zorder=2)
            else:
                ax.plot(T, v, "o", ms=3.0, mfc="none", mec="#1f77b4", alpha=0.65, zorder=4)
        s = next((r for r in summ if r["property"] == prop), None)
        st = s["status"] if s else "?"
        ax.axvline(T_REF, color="#2ca02c", lw=1.0, ls="--")
        ax.set_title(f"{name}  [{ulabel}]\n{st}" +
                     (f" – gewählt {sel}" if sel else ""), fontsize=9.5)
        ax.set_xlabel("Temperatur [K]")
        ax.set_ylabel(f"{name} [{ulabel}]")
        ax.grid(alpha=0.25)
        if prop in ("dynamic_viscosity", "electrical_conductivity"):
            ax.set_yscale("log")
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].set_visible(False)

    from matplotlib.lines import Line2D
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", ms=6.5, mfc="#d62728", mec="#111111",
               label="gewählte Quelle (Auswahlregel, siehe docs/13)"),
        Line2D([], [], marker="o", ls="", ms=4, mfc="none", mec="#1f77b4",
               label="andere Quelle, Umgebungsdruck"),
        Line2D([], [], marker="x", ls="", ms=5, color="#7f7f7f",
               label="Punkt abseits des Umgebungsdrucks – trägt keinen Umgebungswert"),
        Line2D([], [], marker="v", ls="", ms=5, color="#9467bd",
               label="frequenzaufgelöst (1–18 GHz) – KEIN Gleichstromwert"),
        Line2D([], [], color="#2ca02c", ls="--", label="298,15 K")],
        loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.098))

    caveat(fig, "Zwei Größen haben KEINE gewählte Quelle und melden MissingMaterialData: die "
                "relative Permittivität und die kinematische Viskosität. Keine ihrer Quellen "
                "gibt Methode, Reinheit und Wassergehalt zugleich an; die mehrpunktige "
                "Permittivitätsquelle misst zudem von 1 bis 18 GHz und ist damit keine "
                "Gleichstromgröße. Es wird dort kein Ersatzwert gesetzt.", 0.046)
    provenance(fig, m, NOT_MODELLED, y=0.006)
    fig.tight_layout(rect=[0, 0.165, 1, 0.955])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_scatter(d, m, out):
    src = rows(os.path.join(d, "sources.csv"))
    summ = rows(os.path.join(d, "summary.csv"))
    if not src:
        return None
    props = ["surface_tension", "density", "dynamic_viscosity", "electrical_conductivity"]
    props = [p for p in props if any(r["property"] == p for r in src)]
    fig, axes = plt.subplots(1, len(props), figsize=(6.0 * len(props), 8.0), squeeze=False)
    fig.suptitle(TITLE + "\nAbb. 2 – Literaturstreuung im Fenster 298,15 K ± 2 K, Quelle für Quelle",
                 fontsize=12.5, y=0.985)

    for k, prop in enumerate(props):
        ax = axes[0][k]
        name, unit, scale, ulabel = LABEL[prop]
        # value_near_T, not value_at_T: the same +-2 K window the reported band
        # uses, so that the figure and the band cannot disagree.
        sel = [r for r in src if r["property"] == prop and np.isfinite(f(r["value_near_T"]))]
        sel.sort(key=lambda r: f(r["value_near_T"]))
        y = np.arange(len(sel))
        v = np.array([f(r["value_near_T"]) for r in sel]) * scale
        u = np.array([f(r["uncertainty_at_T"]) for r in sel]) * scale
        u = np.where(np.isfinite(u), u, 0.0)
        colours = ["#d62728" if r["selected"] == "yes" else
                   ("#1f77b4" if r["states_purity"] == "yes" and r["states_water"] == "yes"
                    else "#bbbbbb") for r in sel]
        ax.errorbar(v, y, xerr=u, fmt="none", ecolor="#666666", elinewidth=1.0, capsize=2)
        ax.scatter(v, y, c=colours, s=42, zorder=5, edgecolors="#333333", linewidths=0.5)
        labels = [short_label(r) +
                  ("" if abs(f(r["T_near_K"]) - T_REF) < 0.06
                   else f"  @{f(r['T_near_K']):.1f} K") for r in sel]
        # A sixty-source panel cannot carry sixty labels.  Label the selected
        # source, the two extremes and the manufacturer sheet, and say in the
        # axis title how many are unlabelled -- rather than printing a smear.
        if len(sel) <= 26:
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=6.6)
        else:
            keep = {0, len(sel) - 1}
            for i, r in enumerate(sel):
                if r["selected"] == "yes" or r["database_id"].startswith("IoLiTec"):
                    keep.add(i)
            ax.set_yticks(sorted(keep))
            ax.set_yticklabels([labels[i] for i in sorted(keep)], fontsize=7.0)
            ax.set_ylabel(f"{len(sel)} Quellen, aufsteigend sortiert", fontsize=8)
        s = next((r for r in summ if r["property"] == prop), None)
        spread = f(s["relative_spread"]) if s else np.nan
        illu = f(s["illustrative_value"]) * scale if s else np.nan
        if np.isfinite(illu):
            ax.axvline(illu, color="#7a3b00", lw=1.6, ls="--")
            ax.text(illu, len(sel) - 0.4, " bisher benutzt\n (illustrative)", fontsize=7,
                    color="#7a3b00", va="top")
        ax.set_xlabel(f"{name} bei 298,15 K ± 2 K [{ulabel}]")
        ax.set_title(f"{name}\nStreuung {100 * spread:.1f} %", fontsize=10)
        ax.grid(alpha=0.25, axis="x")

    caveat(fig, "Rot ist die gewählte Quelle, blau eine Quelle, die Reinheit UND Wassergehalt "
                "angibt, grau (○) eine, die das nicht tut. Die braune Linie ist der Wert, den "
                "P3a und P3b bisher benutzt haben. Für γ liegt er unter JEDER Quelle mit "
                "dokumentierter Probe – die untere Bandkante darunter ist eine "
                "Kapillaraufstiegsmessung ohne Reinheits- und Wasserangabe. Die Fehlerbalken "
                "sind die von den Quellen selbst angegebenen; die Streuung zwischen den "
                "Quellen ist um ein Vielfaches größer und ist die ehrliche Unsicherheit.",
            0.052)
    provenance(fig, m, NOT_MODELLED, y=0.006)
    fig.tight_layout(rect=[0, 0.135, 1, 0.950])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_impact(d, m, out):
    imp = rows(os.path.join(d, "impact.csv"))
    if not imp:
        return None
    fig, ax = plt.subplots(figsize=(11.0, 6.2))
    fig.suptitle(TITLE + "\nAbb. 3 – was der belegte γ-Wert an den P3a/P3b-Zahlen ändert",
                 fontsize=12.5, y=0.97)
    y = np.arange(len(imp))
    fsel = col(imp, "factor_selected")
    flo = col(imp, "factor_lo")
    fhi = col(imp, "factor_hi")
    ax.barh(y, fsel - 1.0, left=1.0, color="#1f77b4", height=0.45)
    for i in range(len(imp)):
        ax.plot([flo[i], fhi[i]], [y[i], y[i]], color="#111111", lw=2.0, zorder=5)
        ax.plot([flo[i], fhi[i]], [y[i], y[i]], "|", color="#111111", ms=10, zorder=5)
        ax.text(fsel[i], y[i] + 0.28, f"×{fsel[i]:.3f}", fontsize=8.5, ha="center")
    ax.axvline(1.0, color="#7a3b00", lw=1.6, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['quantity']}\n({r['scaling']})" for r in imp], fontsize=8.5)
    ax.set_xlabel("Faktor gegenüber der Rechnung mit dem bisherigen illustrativen γ")
    ax.grid(alpha=0.25, axis="x")
    ax.set_title("Balken: gewählte Quelle.  Schwarzer Strich: das ganze belegte Band.",
                 fontsize=9.5)

    caveat(fig, "γ geht LINEAR in das Gleichgewicht ein, die Spannung skaliert mit √γ. Das "
                "sind Skalierungen der vorhandenen Rechnung, keine neuen Rechnungen. Der "
                "Faktor 1,19 auf jede Druckskala ist größer als jeder Diskretisierungsfehler, "
                "den P0 gemessen hat (4,6 bis 6,1 %). Die bisherigen P3a/P3b-Zahlen bleiben "
                "gültig als das, was sie sind – Rechnungen mit einem ausdrücklich "
                "illustrativen Stoffwert.", 0.055)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.155, 1, 0.930])
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

    made = [figure_vs_T(d, m, os.path.join(d, "fig1_properties_vs_T.png")),
            figure_scatter(d, m, os.path.join(d, "fig2_scatter.png")),
            figure_impact(d, m, os.path.join(d, "fig3_impact.png"))]
    for x in made:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
