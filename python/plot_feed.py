#!/usr/bin/env python3
"""Figures of the P1 pressure budget at the exit plane.

    ./build/es_feed examples/device_p1.cfg examples/feed_p1.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_feed.py results/<dir>

Reads only the CSVs written by apps/es_feed.cpp.  Nothing here computes a
pressure, a resistance or a meniscus.

Produces
  <dir>/fig1_budget_vs_flow.png   the budget term by term against the flow rate
  <dir>/fig2_sensitivity.png      the viscous term against channel length and
                                  radius, with the exponents it must show
  <dir>/fig3_scales.png           where the budget sits against gamma/a, and
                                  what the meniscus solver makes of it
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

TITLE = "P1: Druckhaushalt am Austritt – KEIN Strömungslöser"

NOT_MODELLED = ("Nicht enthalten: Strömungslöser, kapillarer Aufstieg, bewegliche "
                "Kontaktlinie, Reservoirentleerung, Emission. p_reservoir und Q sind "
                "Eingaben. Der Kapillardruck des Meniskus steckt in P3a/P3b, nicht hier.")

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


def params(d):
    return {r["name"]: float(r["value_SI"]) for r in rows(os.path.join(d, "parameters.csv"))}


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def col(rs, name):
    return np.array([f(r.get(name)) for r in rs])


def provenance(fig, m, par, extra="", y=0.010, width=180):
    dirty = pv.DIRTY_MARK in STAMP or STAMP == pv.NO_VERSION
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_feed')}   |   "
            f"a = {par.get('contact_radius', 0) * 1e6:.2f} µm, "
            f"γ/a = {par.get('gamma_over_a', 0):.0f} Pa, "
            f"µ = {par.get('mu', 0) * 1e3:.1f} mPa·s, ρ = {par.get('rho', 0):.0f} kg/m³ "
            f"(Stoffstatus {m.get('liquid_status', '?')})   |   "
            f"Kanal R = {par.get('channel_radius', 0) * 1e6:.2f} µm, "
            f"L = {par.get('channel_length', 0) * 1e6:.0f} µm")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.0, color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {STAMP}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=176, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


# ==========================================================================
def figure_flow(d, m, par, out):
    rs = rows(os.path.join(d, "budget_vs_flow.csv"))
    if not rs:
        return None
    Q = col(rs, "Q_m3_per_s")
    ga = par.get("gamma_over_a", 1.0)
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.4))
    fig.suptitle(TITLE + "\nAbb. 1 – Δp_exit = (p_res − p_vak) − Δp_hydrostatisch − "
                 "Δp_viskos, gegen den Volumenstrom", fontsize=12.5, y=0.98)

    ax = axes[0]
    ax.plot(Q * 1e15, col(rs, "driving_Pa"), lw=1.8, color="#2ca02c",
            label="Antrieb p_res − p_vak")
    ax.plot(Q * 1e15, -col(rs, "hydrostatic_Pa"), lw=1.8, color="#1f77b4",
            label="− Δp_hydrostatisch")
    ax.plot(Q * 1e15, -col(rs, "viscous_Pa"), lw=1.8, color="#d62728",
            label="− Δp_viskos (Hagen-Poiseuille)")
    ax.plot(Q * 1e15, col(rs, "delta_p_exit_Pa"), lw=2.4, color="#111111",
            label="Δp_exit (Summe)")
    ax.axhline(0, color="#888888", lw=0.8)
    ax.axhline(2 * ga, color="#7a3b00", ls="--", lw=1.1, label="±2 γ/a: Grenze der\ngepinnten statischen Form")
    ax.axhline(-2 * ga, color="#7a3b00", ls="--", lw=1.1)
    ax.set_xlabel("Volumenstrom Q [pL/s = 10⁻¹⁵ m³/s]  – EINGABE")
    ax.set_ylabel("Druck [Pa]")
    ax.set_title("die Terme des Haushalts", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.2)

    ax = axes[1]
    ax.plot(Q * 1e15, col(rs, "Pi"), lw=2.0, color="#111111")
    ax.axhspan(-2, 2, color="#2ca02c", alpha=0.07)
    ax.axhline(2, color="#7a3b00", ls="--", lw=1.1)
    ax.axhline(-2, color="#7a3b00", ls="--", lw=1.1)
    ax.set_xlabel("Q [pL/s]")
    ax.set_ylabel("Π = Δp_exit / (γ/a)")
    ax.set_title("dieselbe Größe auf der Kapillarskala", fontsize=10)
    ax.grid(alpha=0.25)

    ax = axes[2]
    ax.semilogy(Q * 1e15, np.maximum(col(rs, "reynolds"), 1e-30), lw=1.8, color="#1f77b4",
                label="Reynoldszahl")
    ax.axhline(2300, color="#d62728", lw=1.4, label="laminare Schranke 2300")
    ax.semilogy(Q * 1e15, np.maximum(col(rs, "entrance_fraction"), 1e-30), lw=1.8,
                color="#ff7f0e", label="Einlauflänge / L")
    ax.axhline(0.05, color="#7a3b00", ls="--", lw=1.1, label="zugelassener Anteil 0,05")
    ax.set_xlabel("Q [pL/s]")
    ax.set_ylabel("dimensionslos")
    ax.set_title("Gültigkeitsgrenzen der geschlossenen Form", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    caveat(fig, "Q ist eine EINGABE. Dieser Lauf berechnet keinen Volumenstrom – dafür "
                "bräuchte es die Emission, die es in dieser Phase nicht gibt. Die Kurve ist "
                "eine Parameterstudie und kein Betriebspunkt. Absolute Drücke stützen sich "
                "auf µ und ρ mit Status illustrative.", 0.070)
    provenance(fig, m, par, NOT_MODELLED, y=0.012)
    fig.tight_layout(rect=[0, 0.150, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_sensitivity(d, m, par, out):
    rs = rows(os.path.join(d, "sensitivity.csv"))
    if not rs:
        return None
    ga = par.get("gamma_over_a", 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.4))
    fig.suptitle(TITLE + "\nAbb. 2 – Empfindlichkeit des viskosen Terms: R_h ~ L und "
                 "R_h ~ R⁻⁴", fontsize=12.5, y=0.98)

    for ax, sweep, xlabel, expo in ((axes[0], "length", "Kanallänge L [µm]", 1.0),
                                    (axes[1], "radius", "Kanalradius R [µm]", -4.0)):
        sr = [r for r in rs if r["sweep"] == sweep]
        x = col(sr, "value_m") * 1e6
        y = col(sr, "viscous_Pa")
        ax.loglog(x, np.abs(y), lw=2.0, color="#d62728", label="Δp_viskos, gerechnet")
        good = np.isfinite(y) & (np.abs(y) > 0)
        if good.sum() >= 2:
            ref = np.abs(y[good][0]) * (x[good] / x[good][0]) ** expo
            ax.loglog(x[good], ref, "--", color="#888888",
                      label=f"Steigung {expo:+.0f} (geschlossene Form)")
        ax.axhline(ga, color="#2ca02c", lw=1.4, label="γ/a")
        ax.axhline(2 * ga, color="#7a3b00", ls="--", lw=1.1, label="2 γ/a")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("|Δp_viskos| [Pa]")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=7.6)
        ax.set_title(f"Referenzstrom, Sweep über {'die Länge' if sweep == 'length' else 'den Radius'}",
                     fontsize=10)

    caveat(fig, "Die gestrichelten Geraden sind nicht angepasst: sie sind die Exponenten, "
                "die 8 µL Q/(πR⁴) verlangt. Dass die Rechnung auf ihnen liegt, ist die "
                "Prüfung. Wo der viskose Term die Linie 2 γ/a überschreitet, existiert für "
                "diesen Zulauf keine gepinnte statische Meniskusform mehr – das ist eine "
                "Aussage über den Druckhaushalt, keine über einen Betriebspunkt.", 0.070)
    provenance(fig, m, par, NOT_MODELLED, y=0.012)
    fig.tight_layout(rect=[0, 0.150, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_scales(d, m, par, out):
    sc = rows(os.path.join(d, "pressure_scales.csv"))
    cp = rows(os.path.join(d, "coupled_to_p3a.csv"))
    if not sc:
        return None
    ga = par.get("gamma_over_a", 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.4))
    fig.suptitle(TITLE + "\nAbb. 3 – Einordnung gegen γ/a und was der P3a-Löser aus dem "
                 "Haushalt macht", fontsize=12.5, y=0.98)

    ax = axes[0]
    names = [r["scale"] for r in sc]
    vals = np.abs(col(sc, "value_Pa"))
    colours = ["#2ca02c" if n.startswith("gamma") or n.startswith("capillary") else "#1f77b4"
               for n in names]
    ax.barh(range(len(names)), np.maximum(vals, 1e-12), color=colours)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("|Druck| [Pa]")
    ax.axvline(ga, color="#2ca02c", lw=1.4)
    ax.set_title("Druckskalen nebeneinander", fontsize=10)
    ax.grid(alpha=0.25, axis="x", which="both")
    for i, (v, r) in enumerate(zip(vals, sc)):
        ax.text(max(v, 1e-12) * 1.15, i, f"{f(r['relative_to_gamma_over_a']):+.3g} γ/a",
                va="center", fontsize=7)

    ax = axes[1]
    if cp:
        p = col(cp, "p_reservoir_Pa") / ga
        h = col(cp, "h_over_a")
        ok = np.isfinite(h)
        ax.plot(p[ok], h[ok], "o-", ms=3.5, color="#111111",
                label="P3a-Apexhöhe aus dem Haushalt")
        bad = ~ok
        if bad.any():
            for x, r in zip(p[bad], [r for r, o in zip(cp, ok) if not o]):
                ax.axvline(x, color="#d62728", lw=0.6, alpha=0.5)
            ax.plot([], [], color="#d62728", lw=0.6,
                    label="keine gepinnte statische Form (Lücke, kein Nullwert)")
        ax.axhline(0, color="#888888", lw=0.8)
        ax.set_xlabel("p_reservoir / (γ/a)")
        ax.set_ylabel("h/a")
        ax.set_title("dieselbe Rechnung wie P3a, nur mit Δp_exit aus dem Haushalt",
                     fontsize=10)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.6)

    caveat(fig, "Rechts ist KEINE neue Physik: es ist der unveränderte P3a-Kapillarlöser, "
                "dem statt einer freien Eingabe der Wert aus dem Druckhaushalt übergeben "
                "wird. Wo kein gepinnter statischer Meniskus existiert, steht eine Lücke "
                "mit ihrem Status und kein Nullwert.", 0.070)
    provenance(fig, m, par, NOT_MODELLED, y=0.012)
    fig.tight_layout(rect=[0, 0.150, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


def main(d):
    global STAMP
    m = meta(d)
    par = params(d)
    root = pv.repo_root_of(__file__)
    rel = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/")
    state = pv.repo_state(root, ignore=[rel])
    STAMP = pv.stamp(state, m.get("commit"))
    pv.write_provenance(d, state, m.get("commit"))

    made = [figure_flow(d, m, par, os.path.join(d, "fig1_budget_vs_flow.png")),
            figure_sensitivity(d, m, par, os.path.join(d, "fig2_sensitivity.png")),
            figure_scales(d, m, par, os.path.join(d, "fig3_scales.png"))]
    for x in made:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
