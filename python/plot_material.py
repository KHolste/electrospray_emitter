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
                                   numbers P3a/P3b reported -- and, as its own
                                   row, what it demonstrably does NOT change
  <dir>/fig4_kinematic_viscosity.png
                                   nu = mu/rho: a DERIVED value, with both
                                   parents, the propagated uncertainty and the
                                   directly measured band it is compared against
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
        loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.140))

    caveat(fig, "Zwei Größen haben KEINE gewählte Quelle – aber aus verschiedenen Gründen "
                "und mit verschiedenem Ausgang. Für die KINEMATISCHE VISKOSITÄT nennt keine "
                "der drei Quellen Methode, Reinheit und Wassergehalt zugleich, also wählt die "
                "Regel keine direkt gemessene aus; ein Wert wird aber aus µ und ρ ABGELEITET "
                "(Status derived, Abb. 4), weil beide Elternwerte nachweislich dieselbe Probe "
                "beschreiben. Für die RELATIVE PERMITTIVITÄT bleibt es bei "
                "MissingMaterialData: auch dort nennt keine Quelle Reinheit und Wassergehalt. "
                "Dass die mehrpunktige Quelle von 1 bis 18 GHz misst, ist dabei KEIN "
                "Ablehnungsgrund – für die Ladungsrelaxationszeit ist genau diese Frequenz die "
                "richtige (P3, docs/14). Es wird nirgends ein Ersatzwert gesetzt.", 0.040)
    provenance(fig, m, NOT_MODELLED, y=0.006)
    fig.tight_layout(rect=[0, 0.205, 1, 0.955])
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
    """What a change of gamma does to the P3a/P3b numbers.

    Every row is a SCALING of an already computed dimensionless solution, never
    a new coupled simulation -- impact.csv carries `recomputed` on every row and
    this figure refuses to draw a row that claims otherwise.

    The row that matters most is the one that does NOT move: at fixed geometry,
    fixed applied voltage and fixed permittivity distribution the Maxwell
    traction is a functional of the field alone, and gamma appears nowhere in
    the field problem.  It is drawn as a marker at exactly 1, not as a bar.
    """
    imp = rows(os.path.join(d, "impact.csv"))
    if not imp:
        return None
    if any(r.get("recomputed", "no") != "no" for r in imp):
        raise SystemExit("impact.csv behauptet eine neu gerechnete Loesung -- "
                         "diese Abbildung stellt nur Skalierungen dar.")

    fig, ax = plt.subplots(figsize=(12.4, 6.6))
    fig.suptitle(TITLE + "\nAbb. 3 – wie eine Änderung von γ die P3a/P3b-Zahlen "
                         "verschiebt – und was sie nicht verschiebt",
                 fontsize=12.5, y=0.975)
    y = np.arange(len(imp))
    fsel = col(imp, "factor_selected")
    flo = col(imp, "factor_lo")
    fhi = col(imp, "factor_hi")
    expo = col(imp, "exponent")
    invariant = np.array([r.get("category", "").startswith("invariant") for r in imp])

    for i in range(len(imp)):
        if invariant[i]:
            # No bar: a bar of length zero would still read as "a small effect".
            ax.plot([1.0], [y[i]], marker="D", ms=11, color="#7a3b00", zorder=6)
            ax.annotate("skaliert NICHT mit γ – Faktor exakt 1",
                        xy=(1.0, y[i]), xytext=(1.012, y[i]), va="center",
                        fontsize=9.0, color="#7a3b00", weight="bold")
        else:
            ax.barh(y[i], fsel[i] - 1.0, left=1.0, color="#1f77b4", height=0.42, zorder=3)
            ax.plot([flo[i], fhi[i]], [y[i], y[i]], color="#111111", lw=2.0, zorder=5)
            ax.plot([flo[i], fhi[i]], [y[i], y[i]], "|", color="#111111", ms=10, zorder=5)
            ax.text(fsel[i], y[i] + 0.27, f"×{fsel[i]:.3f}", fontsize=8.5, ha="center")

    ax.axvline(1.0, color="#7a3b00", lw=1.6, ls="--", zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['quantity']}\nγ-Gesetz: {r['law']}  (Exponent {e:+.1f})"
                        f"\nfestgehalten: {r.get('held_fixed', '')}"
                        for r, e in zip(imp, expo)], fontsize=7.6)
    ax.set_xlabel("Faktor gegenüber der Rechnung mit dem bisherigen illustrativen γ")
    ax.grid(alpha=0.25, axis="x")
    ax.set_title("Balken: gewählte Quelle.  Schwarzer Strich: das ganze belegte Band.  "
                 "Raute: Größe, die von γ nicht abhängt.", fontsize=9.5)
    allf = np.concatenate([fsel, flo, fhi])
    ax.set_xlim(float(np.nanmin(allf)) - 0.07, float(np.nanmax(allf)) + 0.12)

    caveat(fig, "KEINE Zeile ist eine neu gerechnete gekoppelte Simulation. Alle Zeilen "
                "skalieren eine bereits vorhandene dimensionslose Lösung, in der γ nur über "
                "die Druckskala γ/a und die elektrische Bondzahl Γ = ε₀E²a/(2γ) auftritt. "
                "Deshalb: Druckskala ~ γ, Spannung bei gleicher dimensionsloser Form ~ √γ, "
                "Bondzahl bei festem Feld ~ 1/γ. Die Maxwell-Traktion ε₀E²/2 bei "
                "festgehaltener Geometrie, Spannung und Permittivitätsverteilung hängt "
                "dagegen gar nicht von γ ab – γ kommt im Feldproblem nicht vor. γ entscheidet, "
                "WELCHE Form im Gleichgewicht steht, nicht welche Kraft ein gegebenes Feld auf "
                "eine gegebene Form ausübt. Eine frühere Fassung führte genau diese Zeile "
                "falsch als „linear in gamma“.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.175, 1, 0.925])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_kinematic(d, m, out):
    """nu = mu/rho -- derived, and labelled as derived everywhere it appears."""
    kv = rows(os.path.join(d, "kinematic_viscosity.csv"))
    if not kv:
        return None
    r = kv[0]
    derived = r.get("status") == "derived"

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.6),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})
    fig.suptitle(TITLE + "\nAbb. 4 – kinematische Viskosität ν: abgeleitet, nicht gemessen",
                 fontsize=12.5, y=0.975)

    # ---------------------------------------------------------- left: the record
    ax = axes[0]
    ax.axis("off")
    if not derived:
        ax.text(0.0, 0.95, "ν ist NICHT ableitbar.", fontsize=12, va="top", weight="bold",
                color="#8a1b1b", transform=ax.transAxes)
        ax.text(0.0, 0.85, textwrap.fill(r.get("blocker", ""), 56), fontsize=9.5, va="top",
                color="#8a1b1b", transform=ax.transAxes)
    else:
        mu, dmu = f(r["mu"]), f(r["mu_uncertainty"])
        rho, drho = f(r["rho"]), f(r["rho_uncertainty"])
        nu, dnu = f(r["nu"]), f(r["uncertainty_quadratic"])
        rel = f(r["relative_uncertainty"])
        lin = f(r["uncertainty_linear"])

        y = 0.985
        def line(txt, size=10.0, colr="#111111", bold=False, gap=0.062, mono=False):
            nonlocal y
            if txt:
                ax.text(0.0, y, txt, fontsize=size, color=colr, va="top",
                        family="monospace" if mono else None,
                        weight="bold" if bold else "normal", transform=ax.transAxes)
            y -= gap

        line(r"$\nu \;=\; \mu \,/\, \rho$", 15, "#111111", True, 0.085)
        line(f"µ  = (36.370 ± 0.760) mPa·s".replace("36.370", f"{mu*1e3:.3f}")
                                          .replace("0.760", f"{dmu*1e3:.3f}"),
             10.5, "#111111", False, 0.050, mono=True)
        line(f"ρ  = ({rho:.1f} ± {drho:.1f}) kg/m³", 10.5, "#111111", False, 0.075, mono=True)
        line(f"ν  = ({nu*1e6:.4f} ± {dnu*1e6:.4f}) mm²/s", 12.0, "#1f77b4", True, 0.050,
             mono=True)
        line(f"     relative Unsicherheit {100*rel:.2f} %, fortgepflanzt",
             9.0, "#444444", False, 0.040, mono=True)
        line(f"     konservativ linear addiert: ± {lin*1e6:.4f} mm²/s",
             9.0, "#444444", False, 0.085, mono=True)
        line("Status: DERIVED – die Größe selbst wurde NICHT gemessen.",
             10.5, "#7a3b00", True, 0.075)

        cond = ("Geprüfte Bedingungen (C1–C4): beide Elternwerte sind die vom Vertrag "
                "AUSGEWÄHLTEN Quellen für µ und ρ, bei derselben Temperatur, bei "
                "Umgebungsdruck, nicht frequenzaufgelöst, mit wörtlich gleicher Reinheit, "
                "gleichem Wassergehalt und gleicher Probenherkunft"
                + (" – hier sogar dieselbe Publikation."
                   if r.get("same_publication") == "yes" else "."))
        wrapped = textwrap.fill(cond, 58)
        ax.text(0.0, y, wrapped, fontsize=8.6, va="top", color="#333333",
                transform=ax.transAxes)
        y -= 0.038 * (wrapped.count("\n") + 1) + 0.050

        w2 = textwrap.fill("Bedingungen wörtlich: " + r.get("conditions", ""), 62)
        ax.text(0.0, y, w2, fontsize=7.6, va="top", color="#666666", transform=ax.transAxes)

    # ------------------------------------------------- right: against the direct band
    ax = axes[1]
    if derived:
        nu, dnu = f(r["nu"]) * 1e6, f(r["uncertainty_quadratic"]) * 1e6
        dlo, dhi = f(r["direct_band_lo"]) * 1e6, f(r["direct_band_hi"]) * 1e6
        lo = min(dlo, nu - dnu)
        hi = max(dhi, nu + dnu)
        pad = 0.22 * (hi - lo)
        ax.set_ylim(lo - 0.45 * pad, hi + 1.5 * pad)

        ax.axhspan(dlo, dhi, color="#bbbbbb", alpha=0.55, zorder=1)
        ax.text(1.42, 0.5 * (dlo + dhi),
                "direkt gemessenes Band\n(erfüllt die Auswahlregel NICHT:\n"
                "keine Reinheits- und Wasserangabe)\ngeht in keine Rechnung ein",
                fontsize=8.2, ha="center", va="center", color="#333333")
        ax.errorbar([0.62], [nu], yerr=[dnu], fmt="D", ms=10, color="#1f77b4",
                    capsize=7, lw=2.0, zorder=5)
        ax.text(0.62, nu + dnu + 0.06 * (hi - lo), f"abgeleitet\n{nu:.4f} mm²/s",
                fontsize=9.5, ha="center", va="bottom", color="#1f77b4", weight="bold")
        dev = f(r["deviation_from_direct_band_hi"])
        ax.set_title(f"Der abgeleitete Wert liegt {100*dev:+.1f} % über der oberen Kante des\n"
                     "direkt gemessenen Bandes – berichtet, nicht weggerechnet",
                     fontsize=9.3, pad=10)
        ax.set_xlim(0.0, 2.3)
        ax.set_xticks([])
        ax.set_ylabel("ν bei 298,15 K  [mm²/s]")
        ax.grid(alpha=0.25, axis="y")
    else:
        ax.axis("off")

    caveat(fig, "ν ist keine unabhängige Stoffgröße, sondern µ/ρ. Der Vertrag wählt keine "
                "direkt gemessene ν aus – keine der drei Quellen nennt Reinheit und "
                "Wassergehalt. Aus µ und ρ DARF abgeleitet werden, sobald gezeigt ist, dass "
                "beide Werte dieselbe Flüssigkeit im selben Zustand beschreiben; die "
                "Bedingungen C1–C4 stehen in include/es/material_data.hpp und schlagen "
                "einzeln benannt fehl, wenn eine verletzt ist. Der Wert trägt den Status "
                "„derived“ und darf nirgends als Messung von ν zitiert werden. Die Abweichung "
                "vom direkt gemessenen Band ist kleiner als die Literaturstreuung von µ "
                "(54 %) und damit durch Probenunterschiede erklärbar.", 0.088)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.225, 1, 0.925])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
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
            figure_impact(d, m, os.path.join(d, "fig3_impact.png")),
            figure_kinematic(d, m, os.path.join(d, "fig4_kinematic_viscosity.png"))]
    for x in made:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
