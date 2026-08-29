#!/usr/bin/env python3
"""Figures of the P3 transport step.

    ./build/es_transport results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_transport.py results/<dir>

Reads only the CSVs written by apps/es_transport.cpp.

Produces
  <dir>/fig1_pipe_flow.png     the solved profile against Hagen-Poiseuille, and
                               the mesh convergence of the flow rate and of the
                               hydraulic resistance P1 asserts
  <dir>/fig2_charge.png        charge relaxation and the steady conduction check
  <dir>/fig3_time_scales.png   the time scales side by side, with what is
                               sourced and what is not
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

TITLE = "P3: endliche Leitfähigkeit und Zulaufströmung – minimale ehrliche Stufe"

NOT_MODELLED = ("Nicht enthalten: allgemeiner Stokes-Löser (keine Einlaufströmung, keine "
                "Druck-Geschwindigkeits-Kopplung, keine gekrümmte Geometrie, keine freie "
                "Oberfläche), gekoppelte finite-conductivity-Meniskuskopplung, Emission, "
                "Zeitintegration einer Form.")

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
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_transport')}   |   "
            f"R = {f(m.get('channel_radius_m')) * 1e6:.2f} µm, "
            f"L = {f(m.get('channel_length_m')) * 1e6:.0f} µm, T = {f(m.get('T_K')):.2f} K   |   "
            f"µ {m.get('mu_status', '?')}, σ {m.get('sigma_status', '?')}, "
            f"ε_r {m.get('eps_r_status', '?')}")
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
def figure_pipe(d, m, out):
    prof = rows(os.path.join(d, "pipe_profile.csv"))
    conv = rows(os.path.join(d, "pipe_convergence.csv"))
    if not prof:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.2))
    fig.suptitle(TITLE + "\nAbb. 1 – gelöste Rohrströmung gegen Hagen-Poiseuille",
                 fontsize=12.5, y=0.975)

    ax = axes[0]
    x = col(prof, "r_over_R")
    ax.plot(x, col(prof, "u_closed_form_m_per_s"), lw=3.0, color="#bbbbbb",
            label="geschlossene Form  u = −(dp/dz)(R²−r²)/(4µ)")
    ax.plot(x, col(prof, "u_solved_m_per_s"), lw=1.4, color="#111111",
            label="FEM, achsensymmetrisch")
    ax.set_xlabel("r / R")
    ax.set_ylabel("axiale Geschwindigkeit [m/s]")
    ax.set_title("Profil in der Mittelebene", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.6)

    ax = axes[1]
    ax.plot(x, col(prof, "difference"), lw=1.4, color="#d62728")
    ax.axhline(0, color="#888888", lw=0.8)
    ax.set_xlabel("r / R")
    ax.set_ylabel("FEM − geschlossene Form [m/s]")
    ax.set_title("Differenz – der Rest ist Diskretisierung", fontsize=10)
    ax.grid(alpha=0.25)

    ax = axes[2]
    if conv:
        n = col(conv, "nr")
        ax.loglog(n, col(conv, "relative_error"), "o-", color="#1f77b4",
                  label="Volumenstrom gegen Hagen-Poiseuille")
        ax.loglog(n, col(conv, "R_h_relative_error"), "s-", color="#2ca02c",
                  label="hydraulischer Widerstand gegen die Formel von P1")
        e = col(conv, "relative_error")
        good = np.isfinite(e) & (e > 0)
        if good.sum() >= 2:
            ax.loglog(n[good], e[good][0] * (n[good] / n[good][0]) ** -2.0, "--",
                      color="#888888", label="2. Ordnung")
    ax.set_xlabel("Knoten in radialer Richtung")
    ax.set_ylabel("relativer Fehler")
    ax.set_title("Netzkonvergenz", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.4)

    caveat(fig, "Die Reduktion auf ein Skalarproblem ist für das GERADE Rohr exakt: die "
                "Kontinuität erzwingt du_z/dz = 0, der konvektive Term verschwindet identisch, "
                "und die z-Impulsgleichung ist der achsensymmetrische Laplace-Operator. Das "
                "ist KEIN Stokes-Löser: keine Einlaufströmung, keine Druck-Geschwindigkeits-"
                "Kopplung, keine gekrümmte Geometrie, keine freie Oberfläche. Rechts wird der "
                "hydraulische Widerstand, den P1 als Formel behauptet, gegen den aus dem "
                "gelösten Feld gewonnenen geprüft – die beiden teilen keinen Code.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.150, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_charge(d, m, out):
    dec = rows(os.path.join(d, "decay.csv"))
    con = rows(os.path.join(d, "conduction.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.2))
    fig.suptitle(TITLE + "\nAbb. 2 – Ladungsrelaxation und stationärer Leitungsstrom",
                 fontsize=12.5, y=0.975)

    ax = axes[0]
    if dec:
        t = col(dec, "t_over_tau")
        ax.semilogy(t, col(dec, "rho_over_rho0"), lw=2.0, color="#111111")
        ax.axvline(1.0, color="#2ca02c", ls="--", lw=1.2)
        ax.text(1.05, 0.5, "t = τ", color="#2ca02c", fontsize=8)
    ax.set_xlabel("t / τ")
    ax.set_ylabel("ρ(t) / ρ₀")
    ax.set_title("ρ(t) = ρ₀ e^(−t/τ),  τ = ε₀ε_r/σ\n"
                 "mit einem AUSDRÜCKLICH UNBELEGTEN ε_r", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")

    ax = axes[1]
    if con:
        n = col(con, "n_nodes")
        ax.loglog(n, np.maximum(col(con, "relative_error"), 1e-18), "o-", color="#1f77b4",
                  label="Strom gegen σAV/L")
        ax.loglog(n, np.maximum(col(con, "potential_error"), 1e-18), "s-", color="#ff7f0e",
                  label="Potential gegen Vz/L")
        ax.loglog(n, np.maximum(col(con, "fem_residual"), 1e-30), "^-", color="#7f7f7f",
                  label="FEM-Residuum [A]")
    ax.set_xlabel("Knoten")
    ax.set_ylabel("relativer Fehler")
    ax.set_title("stationärer Leitungsstrom im Zylinder\ndiv(σ∇φ) = 0", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.4)

    ax = axes[2]
    if con:
        n = col(con, "n_nodes")
        ax.semilogy(n, np.maximum(col(con, "lateral_leakage"), 1e-18), "o-", color="#d62728",
                    lw=2.0)
        ax.axhline(1e-10, color="#7a3b00", ls="--", lw=1.2,
                   label="Schranke des Tests, 1e−10")
    ax.set_xscale("log")
    ax.set_xlabel("Knoten")
    ax.set_ylabel("|j_r| an der Mantelfläche / mittleres |j_z|")
    ax.set_title("kein Strom durch eine zero-flux-Fläche", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.4)

    caveat(fig, "Rechts steht die physikalische Aussage, nicht eine numerische Feinheit: "
                "OHNE EMISSION darf kein stationärer Normalstrom durch die freie Oberfläche "
                "fließen, weil die Ladung dort nirgendwohin könnte. Die richtige Bedingung "
                "ist j·n = 0 – die natürliche Bedingung desselben Operators – und dass der "
                "Löser sie einhält, wird gemessen. Links: ε_r ist für EMI-BF4 "
                "MissingMaterialData (docs/13), also ist τ mit den belegten Daten dieses "
                "Projekts NICHT berechenbar. Die Kurve zeigt die Form des Zerfalls mit einem "
                "ausdrücklich unbelegten Wert und ist keine Stoffaussage.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.155, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_scales(d, m, out):
    ts = rows(os.path.join(d, "time_scales.csv"))
    rel = rows(os.path.join(d, "relaxation.csv"))
    if not ts:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.2))
    fig.suptitle(TITLE + "\nAbb. 3 – die Zeitskalen nebeneinander", fontsize=12.5, y=0.975)

    ax = axes[0]
    names, vals, cols = [], [], []
    for r in ts:
        v = f(r["value_s"])
        names.append(r["scale"].replace("_", " "))
        vals.append(v)
        cols.append("#2ca02c" if r["status"] in ("sourced", "measured")
                    else ("#bbbbbb" if not np.isfinite(v) else "#ff7f0e"))
    y = np.arange(len(names))
    plotted = [v if np.isfinite(v) else np.nan for v in vals]
    ax.barh(y, plotted, color=cols)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8.5)
    ax.set_xlabel("Zeit [s]")
    ax.grid(alpha=0.25, axis="x", which="both")
    for i, (v, r) in enumerate(zip(vals, ts)):
        if np.isfinite(v):
            ax.text(v * 1.25, i, f"{v:.3g} s", va="center", fontsize=7.5)
        else:
            lo = ax.get_xlim()[0]
            ax.text(lo * 3.0, i, "nicht berechenbar – MissingMaterialData", va="center",
                    fontsize=8, color="#8a1b1b")
    ax.set_title("grün: aus belegten Stoffdaten.  orange: mit einem unbelegten ε_r.",
                 fontsize=9.5)

    ax = axes[1]
    if rel:
        lab = [r["case"].replace("_", " ") for r in rel]
        ratio = col(rel, "ratio")
        y = np.arange(len(lab))
        colours = ["#2ca02c" if r["verdict"] == "PerfectConductorJustified"
                   else ("#bbbbbb" if r["verdict"] == "MissingMaterialData" else "#d62728")
                   for r in rel]
        ax.barh(y, np.where(np.isfinite(ratio), ratio, np.nan), color=colours)
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels(lab, fontsize=8.5)
        ax.axvline(f(m.get("perfect_conductor_margin", 100)), color="#7a3b00", ls="--",
                   lw=1.4, label="geforderter Abstand")
        for i, (r, v) in enumerate(zip(rel, ratio)):
            ax.text(v * 1.2 if np.isfinite(v) else 1.0, i, r["verdict"], va="center",
                    fontsize=7.5,
                    color="#8a1b1b" if r["verdict"] == "MissingMaterialData" else "#333333")
        ax.set_xlabel("Prozesszeit / τ")
        ax.set_title("der Perfect-Conductor-Grenzfall ist ein VERHÄLTNIS", fontsize=9.5)
        ax.grid(alpha=0.25, axis="x", which="both")
        ax.legend(fontsize=7.6)

    caveat(fig, "Welche Prozesszeit gilt, ist eine Modellentscheidung und keine Tatsache: für "
                "eine statische Form ist es die Zeit, in der sich die Form einstellt, für "
                "einen emittierenden Betrieb die Transitzeit durch die Emissionszone. Diese "
                "zweite ist hier NICHT gerechnet. Und mit den belegten Stoffdaten dieses "
                "Projekts ist τ überhaupt nicht berechenbar, weil ε_r fehlt – der "
                "Äquipotentialansatz von P3b ist damit eine Annahme und kein nachgewiesener "
                "Grenzfall.", 0.058)
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
    made = [figure_pipe(d, m, os.path.join(d, "fig1_pipe_flow.png")),
            figure_charge(d, m, os.path.join(d, "fig2_charge.png")),
            figure_scales(d, m, os.path.join(d, "fig3_time_scales.png"))]
    for x in made:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
