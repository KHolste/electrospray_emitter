#!/usr/bin/env python3
"""Figure of the P8 cone-jet contract.  The mode is BLOCKED.

    ./build/es_cone_jet results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_cone_jet.py results/<dir>

The figure is a STATUS AND VALIDITY MAP plus a dimensionless diagnosis.  It is
explicitly not a regime prediction: there is no current, no jet radius and no
droplet diameter anywhere.

Produces
  <dir>/fig1_status.png
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

TITLE = "P8: Cone-Jet und Tropfenbetrieb – Status BLOCKIERT, Diagnose statt Vorhersage"

NOT_MODELLED = ("Kein Cone-Jet-Strom, kein Jetradius, kein Tropfendurchmesser. Die "
                "empirische Skalierung wird NICHT übernommen: Originalgleichung, Erratum, "
                "Vorfaktor und Gültigkeitsbereich waren in diesem Lauf an keiner Quelle "
                "prüfbar. Der alte ConeJetModel-Block bleibt unbenutzt.")

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
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_cone_jet')}   |   "
            f"Status {m.get('status', '?')}   |   "
            f"Korrelation übernommen: {m.get('correlation_adopted', '?')}   |   "
            f"a = {f(m.get('contact_radius_m')) * 1e6:.2f} µm")
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
    rq = rows(os.path.join(d, "requirements.csv"))
    dg = rows(os.path.join(d, "diagnosis.csv"))
    eb = rows(os.path.join(d, "electric_bond.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 7.4))
    fig.suptitle(TITLE + "\nAbb. 1 – Gültigkeitskarte und dimensionslose Diagnose",
                 fontsize=12.0, y=0.975)

    ax = axes[0]
    if rq:
        y = np.arange(len(rq))
        ok = np.array([r["available"] == "yes" for r in rq])
        ax.barh(y, np.ones(len(rq)), color=np.where(ok, "#2ca02c", "#d62728"))
        ax.set_yticks(y)
        ax.set_yticklabels(["\n".join(textwrap.wrap(r["requirement"], 34)) for r in rq],
                           fontsize=7.0)
        for i, r in enumerate(rq):
            ax.text(0.02, i, ("vorhanden – " if ok[i] else "FEHLT – ") + r["provided_by"],
                    va="center", fontsize=6.6, color="white", weight="bold")
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_title("welches Teilmodell eine Cone-Jet-Rechnung braucht,\n"
                     "und ob dieses Projekt es hat", fontsize=9.5)

    ax = axes[1]
    if dg:
        Q = col(dg, "Q_m3_per_s") * 1e15
        for key, lab in (("t_capillary_s", "kapillar-inertial  √(ρa³/γ)"),
                         ("t_viscous_s", "viskokapillar  µa/γ"),
                         ("t_residence_s", "Verweilzeit  V/Q"),
                         ("tau_e_s", "Ladungsrelaxation  ε₀ε_r/K")):
            v = col(dg, key)
            if np.all(~np.isfinite(v)):
                ax.plot([], [], "--", color="#8a1b1b",
                        label=lab + "  – nicht berechenbar (ε_r fehlt)")
            elif key == "tau_e_s" and dg[0].get("tau_e_self_consistent") == "yes":
                # NOT from a substitute value: the self-consistent P3 solution,
                # drawn with the justified band it carries.
                lo, hi = col(dg, "tau_e_lo_s"), col(dg, "tau_e_hi_s")
                ax.fill_between(Q, lo, hi, color="#1f77b4", alpha=0.18, zorder=1)
                ax.loglog(Q, v, color="#1f77b4", lw=2.0,
                          label=lab + "  – selbstkonsistent (P3), mit Band")
            else:
                ax.loglog(Q, v, label=lab)
    ax.set_xlabel("Volumenstrom Q [pL/s]  – EINGABE")
    ax.set_ylabel("Zeitskala [s]")
    ax.set_title("die Zeitskalen nebeneinander", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.0)

    ax = axes[2]
    if dg:
        Q = col(dg, "Q_m3_per_s") * 1e15
        for key, lab in (("Oh", "Ohnesorge  µ/√(ργa)"), ("Re", "Reynolds (Zulauf)"),
                         ("Ca", "Kapillarzahl (Zulauf)")):
            ax.loglog(Q, np.maximum(col(dg, key), 1e-30), label=lab)
        ax.axhline(1.0, color="#7a3b00", ls="--", lw=1.2)
        ax.text(Q[1], 1.3, "Kennzahl = 1", fontsize=7.5, color="#7a3b00")
    if eb:
        ax2 = ax.twiny()
        ax2.loglog(col(eb, "E_V_per_m") / 1e6, col(eb, "Bo_E"), color="#9467bd", lw=2.0)
        ax2.set_xlabel("Oberflächenfeld E [MV/m]  →  elektrische Bondzahl (violett)",
                       color="#9467bd")
        ax2.tick_params(axis="x", colors="#9467bd")
    ax.set_xlabel("Volumenstrom Q [pL/s]")
    ax.set_ylabel("dimensionslos")
    ax.set_title("die Kennzahlen, die DEFINITIONEN sind", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.0, loc="lower right")

    caveat(fig, "Dies ist eine DIAGNOSE und keine Regimevorhersage. Jede gezeichnete Kennzahl "
                "ist eine Definition und braucht keine Literaturquelle; die Herleitungen "
                "stehen im Header. Zu ε_r: ein EINZELWERT bleibt MissingMaterialData – keine "
                "Quelle nennt Reinheit und Wassergehalt. Die Ladungsrelaxationszeit und die "
                "elektrohydrodynamische Länge hängen von ε_r aber NUR über τ ab, und τ ist "
                "seit der P3-Korrektur selbstkonsistent auf der gemessenen Dispersionskurve "
                "gelöst; beide tragen deshalb einen Wert und das begründete Band, und keines "
                "davon ist ein Ersatzwert. Eine frühere Fassung führte beide als „nicht "
                "berechenbar“. Die elektrische Bondzahl ist genau das Verhältnis der beiden "
                "Terme, die P3b bilanziert, und damit die einzige Stelle, an der diese "
                "Diagnose etwas berührt, das dieses Projekt gerechnet hat.", 0.058)
    provenance(fig, m, NOT_MODELLED, y=0.008)
    fig.tight_layout(rect=[0, 0.175, 1, 0.930])
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
    x = figure(d, m, os.path.join(d, "fig1_status.png"))
    if x:
        print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
