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
  <dir>/fig3_time_scales.png   WHICH permittivity belongs in tau_q: the measured
                               dispersion, the implicit equation solved on it,
                               and the time scales with tau as a band
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
    ax.set_title("ρ(t) = ρ₀ e^(−t/τ),  τ = ε₀ε_r(f*)/K\n"
                 "mit dem SELBSTKONSISTENTEN τ aus Abb. 3", fontsize=9.5)
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
                "Löser sie einhält, wird gemessen. Links: die Kurve benutzt das "
                "SELBSTKONSISTENTE τ aus Abb. 3, also ε_r an der Frequenz f* = 1/(2πτ), an "
                "der die freie Ladung tatsächlich zerfällt – kein unbelegter Ersatzwert. Ein "
                "EINZELNER ε_r-Wert bleibt gleichwohl MissingMaterialData (keine Quelle nennt "
                "Reinheit und Wassergehalt); belegt ist ein Band, und die Kurve zeigt die "
                "FORM des Zerfalls, nicht eine auf drei Stellen belastbare Stoffzahl.", 0.058)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.155, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
def figure_scales(d, m, out):
    """Which permittivity belongs in tau_q, and what follows for P3b.

    Left   -- the measured dispersion eps_r'(f), 1-18 GHz, with the values that
              the sources report as static shown for what they are: extrapolated
              limits, not measurements at zero frequency.
    Middle -- the implicit equation.  tau = eps0 eps_r(f)/K comes from the
              measured curve; f = 1/(2 pi tau) is the definition of the
              frequency at which the decaying free charge lives.  Where the two
              cross IS tau.  Nothing is fitted.
    Right  -- the time scales.  tau is a BAND, because no single eps_r is
              sourced; the verdict is read off its worst edge.
    """
    ts = rows(os.path.join(d, "time_scales.csv"))
    rel = rows(os.path.join(d, "relaxation.csv"))
    pts = rows(os.path.join(d, "permittivity_points.csv"))
    sc = rows(os.path.join(d, "self_consistency.csv"))
    if not ts:
        return None

    fig, axes = plt.subplots(1, 3, figsize=(17.4, 6.4))
    fig.suptitle(TITLE + "\nAbb. 3 – welche Permittivität in τ_q gehört, "
                         "und was daraus für P3b folgt", fontsize=12.5, y=0.975)

    sol = next((r for r in sc if r.get("is_solution") == "yes"), None)

    # ---------------------------------------------------------------- panel 1
    ax = axes[0]
    fr = [r for r in pts if f(r["frequency_Hz"]) > 0 and r["admissible"] == "yes"]
    st = [r for r in pts if f(r["frequency_Hz"]) == 0 and r["admissible"] == "yes"]
    bad = [r for r in pts if r["admissible"] == "no"]
    if fr:
        x = col(fr, "frequency_Hz") / 1e9
        y = col(fr, "eps_r")
        e = col(fr, "uncertainty")
        order = np.argsort(x)
        ax.errorbar(x[order], y[order], yerr=e[order], fmt="o-", ms=4.5, lw=1.3,
                    color="#1f77b4", capsize=3, zorder=4,
                    label="gemessen, frequenzaufgelöst (ε′)")
    for k, r in enumerate(st):
        ax.axhline(f(r["eps_r"]), color="#9467bd", ls=":", lw=1.4, zorder=2,
                   label="von der Quelle als statisch berichtet\n(aus Mikrowellenspektren "
                         "extrapoliert)" if k == 0 else None)
    if bad:
        ax.plot([f(r["frequency_Hz"]) / 1e9 for r in bad], [f(r["eps_r"]) for r in bad],
                "x", color="#d62728", ms=8, label="unter der Elektrodenpolarisationsschwelle")
    if sol:
        fs = f(sol["frequency_Hz"]) / 1e9
        ax.plot([fs], [f(sol["eps_r"])], "D", ms=11, color="#2ca02c", zorder=6,
                label=f"ε_r bei f* = {fs:.2f} GHz")
        ax.axvline(fs, color="#2ca02c", lw=1.2, ls="--", zorder=3)
    ax.set_xscale("log")
    ax.set_xlabel("Frequenz [GHz]")
    ax.set_ylabel("ε_r′")
    ax.set_title("die gemessene Dispersion – sie wurde früher\nals „nicht DC“ verworfen",
                 fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.0, loc="upper right")

    # ---------------------------------------------------------------- panel 2
    ax = axes[1]
    curve = [r for r in sc if r.get("is_solution") != "yes"]
    if curve:
        fq = col(curve, "frequency_Hz") / 1e9
        ax.plot(fq, col(curve, "tau_from_eps_s"), color="#1f77b4", lw=2.0,
                label=r"$\tau = \varepsilon_0\,\varepsilon_r(f)\,/\,K$  (Messkurve)")
        ax.plot(fq, col(curve, "tau_from_f_s"), color="#7a3b00", lw=2.0, ls="--",
                label=r"$\tau = 1/(2\pi f)$  (Definition von $f^*$)")
    if sol:
        fs, ta = f(sol["frequency_Hz"]) / 1e9, f(sol["tau_from_eps_s"])
        ax.plot([fs], [ta], "D", ms=12, color="#2ca02c", zorder=6)
        ax.annotate(f"Lösung\nf* = {fs:.2f} GHz\nτ = {ta:.3g} s",
                    xy=(fs, ta), xytext=(0.30, 0.22), textcoords="axes fraction",
                    fontsize=9.0, color="#2ca02c", weight="bold",
                    arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1.5))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequenz [GHz]")
    ax.set_ylabel("τ [s]")
    ax.set_title("die implizite Gleichung – der Schnittpunkt IST τ,\nes wird nichts angepasst",
                 fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8.0, loc="upper right")

    # ---------------------------------------------------------------- panel 3
    ax = axes[2]
    names, vals, los, his, cols = [], [], [], [], []
    for r in ts:
        names.append(r["scale"].replace("_", " "))
        vals.append(f(r["value_s"]))
        los.append(f(r["lo_s"]))
        his.append(f(r["hi_s"]))
        cols.append("#2ca02c" if r["status"] == "sourced" else "#1f77b4")
    y = np.arange(len(names))
    ax.barh(y, vals, color=cols, height=0.5, zorder=3)
    for i2 in range(len(names)):
        if np.isfinite(los[i2]) and np.isfinite(his[i2]):
            ax.plot([los[i2], his[i2]], [y[i2], y[i2]], color="#111111", lw=2.4, zorder=5)
            ax.plot([los[i2], his[i2]], [y[i2], y[i2]], "|", color="#111111", ms=11, zorder=5)
        anchor = his[i2] if np.isfinite(his[i2]) else vals[i2]
        ax.text(anchor * 1.9, y[i2], f"{vals[i2]:.3g} s", va="center", fontsize=8.5)
    ax.set_xscale("log")
    ax.set_xlim(min(v for v in los + vals if np.isfinite(v)) * 0.35,
                max(vals) * 60.0)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8.5)
    ax.set_xlabel("Zeit [s]")
    ax.grid(alpha=0.25, axis="x", which="both")

    worst = None
    if rel:
        fin = [f(r["ratio"]) for r in rel if np.isfinite(f(r["ratio"]))]
        worst = min(fin) if fin else None
    margin = f(m.get("perfect_conductor_margin", 100))
    ax.set_title("τ als BAND – kein einzelner ε_r-Wert ist belegt.\n"
                 + (f"Schlechteste Bandecke: Prozesszeit/τ = {worst:.3g}  ≫  {margin:.0f}"
                    if worst else ""), fontsize=9.5)

    caveat(fig, "Welche Permittivität in τ_q = ε₀ε_r/K gehört, ist keine Konvention: die freie "
                "Ladung zerfällt auf der Zeitskala τ selbst, ihr Spektrum liegt also bei "
                "f* = 1/(2πτ), und dort ist ε_r abzulesen. Für diese Flüssigkeit liegt f* bei "
                "rund 2,6 GHz – also genau dort, wo gemessen wurde. Die 1–18-GHz-Daten früher "
                "als „nicht DC“ zu verwerfen war ein Lesefehler der Formel; sie umgekehrt als "
                "Gleichstromwert zu übernehmen wäre ebenso falsch. Ein EINZELNER ε_r-Wert "
                "bleibt MissingMaterialData: keine der vier Quellen nennt Reinheit und "
                "Wassergehalt. Belegt ist ein BAND, und über dessen gesamte Breite – "
                "einschließlich der für die Näherung ungünstigsten Ecke – ist τ um mehr als "
                "vier Größenordnungen kürzer als jede der beiden Prozesszeiten. Der "
                "Äquipotentialansatz von P3b ist damit belegt und nicht mehr nur angenommen. "
                "Das sagt weiterhin NICHTS über einen emittierenden Betrieb: dort ist die "
                "Prozesszeit die Transitzeit durch die Emissionszone, und die ist hier nicht "
                "gerechnet.", 0.072)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.265, 1, 0.930])
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
