#!/usr/bin/env python3
"""Figures of the P0 numerical clean-up of P3b.

    ./build/es_p3b_audit examples/device_p1.cfg examples/electrocapillary_p3b.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_p3b_audit.py results/<dir>

Reads only the CSVs written by apps/es_p3b_audit.cpp.  Nothing here solves,
meshes or fits.

Three rules this script follows, because P3b broke all three:
  * a figure is called "Konvergenz" only where the pre-declared one-per-cent
    target is met; otherwise it is called what it is;
  * a computation that failed is a GAP, never a point at zero;
  * the raw load, the conservative segment load and the load actually handed to
    the capillary solver are drawn as three different things with three
    different labels.

Produces
  <dir>/fig1_load_projection.png   raw nodal / segment staircase / handed load,
                                   and the cumulative force they come from
  <dir>/fig2_projection_error.png  force conservation and the measured order of
                                   the projection against the closed forms
  <dir>/fig3_exclusion.png         kTolExclusion: measurement against the closed
                                   form, over four refinements
  <dir>/fig4_gate.png              the edge gate: direct force sequence and the
                                   two extrapolations, side by side
  <dir>/fig5_coupled.png           coupled mesh study with Richardson estimates
  <dir>/figures_provenance.txt
"""
import csv
import os
import sys
import textwrap
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance as pv

TITLE = "P0: numerische Bereinigung von P3b – Lastprojektion, Netzkonvergenz, Kanten-Gate"

NOT_MODELLED = ("Keine neue Physik. Keine Emission, keine endliche Leitfähigkeit, keine "
                "Strömung, keine Raumladung, keine Zeitabhängigkeit, keine "
                "Stabilitätsaussage. Δp_exit ist in diesem Lauf eine Eingabe. Stoffwerte "
                "illustrative.")

STAMP = ""


# ------------------------------------------------------------------ io
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
        v = float(x)
    except (TypeError, ValueError):
        return np.nan
    return v


def col(rs, name):
    return np.array([f(r.get(name)) for r in rs])


# -------------------------------------------------------------- decoration
def provenance(fig, m, extra="", y=0.010, width=180):
    dirty = pv.DIRTY_MARK in STAMP or STAMP == pv.NO_VERSION
    head = (f"Commit {STAMP}   |   {m.get('app', 'es_p3b_audit')}   |   "
            f"Konfiguration {m.get('config', '?')}   |   "
            f"a = {f(m.get('contact_radius_m')) * 1e6:.2f} µm, "
            f"γ/a = {f(m.get('capillary_pressure_scale_Pa')):.0f} Pa "
            f"(Stoffstatus {m.get('liquid_status', '?')})   |   "
            f"Ziel Diskretisierungsfehler {f(m.get('discretization_target', 0.01)) * 100:.0f} %")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.0, color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {STAMP}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=176, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.6, color=color)


def level_colour(lvl, lo=1, hi=5):
    t = 0.20 + 0.72 * (lvl - lo) / max(1, hi - lo)
    return plt.get_cmap("viridis")(t)


# ==========================================================================
# fig 1 -- the three loads, told apart
# ==========================================================================
def figure_loads(d, m, out):
    rs = rows(os.path.join(d, "projection_profiles.csv"))
    if not rs:
        return None
    cases = ["singulaer", "glatt", "feld"]
    titles = {"singulaer": "hergestellt: p = C d$^{\\beta}$, β = −0,44",
              "glatt": "hergestellt: p = p₀(1 + (r/a)²)",
              "feld": "gelöstes Feld, ebene Fläche, Netzstufe 2"}
    fig, axes = plt.subplots(3, 3, figsize=(16.0, 11.4))
    fig.suptitle(TITLE + "\nAbb. 1 – Rohlast, Segmentprojektion und die TATSÄCHLICH "
                 "übergebene Last", fontsize=12.5, y=0.988)

    by = defaultdict(lambda: defaultdict(list))
    for r in rs:
        by[r["case"]][r["kind"]].append((f(r["tau"]), f(r["value"])))

    for k, case in enumerate(cases):
        ax = axes[0][k]
        data = by.get(case)
        if not data:
            ax.set_visible(False)
            axes[1][k].set_visible(False)
            axes[2][k].set_visible(False)
            continue
        nodes = np.array(sorted(data["node"]))
        seg = np.array(data["segment"])
        handed = np.array(sorted(data["handed"]))
        if len(seg):
            ax.plot(seg[:, 0], seg[:, 1], color="#1f77b4", lw=1.4, drawstyle="default",
                    label="Segmentprojektion (Treppe, konservativ)")
        if len(nodes):
            ax.plot(nodes[:, 0], nodes[:, 1], "o", ms=3.4, color="#d62728",
                    label="Rohlast, punktweise am Knoten")
        if len(handed):
            ax.plot(handed[:, 0], handed[:, 1], color="#111111", lw=1.9,
                    label="ÜBERGEBENE Last p(τ) = G′/A′")
        ax.set_title(titles[case], fontsize=10)
        ax.set_xlabel("normierte Bogenlänge τ = s/L")
        ax.set_ylabel("p [Pa]")
        ax.set_yscale("log")
        ax.grid(alpha=0.25)
        ax.axvspan(0.9, 1.0, color="#d62728", alpha=0.07)
        if k == 0:
            ax.legend(fontsize=7.6, loc="upper left")

        # ZOOM.  At full width the handed load hugs the segment means and looks
        # like a staircase; the difference between "continuous" and "piecewise
        # constant" only shows over a few bins.  So it is shown over a few bins.
        axz = axes[1][k]
        lo, hi = 0.55, 0.66
        if len(seg):
            msk = (seg[:, 0] >= lo - 0.02) & (seg[:, 0] <= hi + 0.02)
            axz.plot(seg[msk, 0], seg[msk, 1], color="#1f77b4", lw=1.6,
                     label="Segmentprojektion (springt)")
        if len(handed):
            msk = (handed[:, 0] >= lo) & (handed[:, 0] <= hi)
            axz.plot(handed[msk, 0], handed[msk, 1], color="#111111", lw=2.2,
                     label="übergebene Last (stetig)")
        for b in range(int(lo * 128), int(hi * 128) + 1):
            axz.axvline(b / 128.0, color="#bbbbbb", lw=0.5, zorder=0)
        axz.set_xlim(lo, hi)
        axz.set_xlabel("τ (Ausschnitt; graue Linien sind die Binränder)")
        axz.set_ylabel("p [Pa]")
        axz.grid(alpha=0.2)
        axz.set_title("Ausschnitt: die übergebene Last springt an keinem Binrand",
                      fontsize=9)
        if k == 0:
            axz.legend(fontsize=7.4)

        ax2 = axes[2][k]
        cum = np.array(sorted(data["cumulative"]))
        if len(cum):
            ax2.plot(cum[:, 0], cum[:, 1] * 1e9, color="#2ca02c", lw=1.8)
        ax2.set_xlabel("τ")
        ax2.set_ylabel("kumulierte Kraft G(τ) [nN]")
        ax2.grid(alpha=0.25)
        ax2.set_title("G(τ): daraus wird die übergebene Last gebaut", fontsize=9)

    caveat(fig, "Die rot hinterlegte Zone ist die Kantenzone. Dort ist die ROTE Rohlast "
                "nicht netzkonvergent und wird nirgends verwendet; die schwarze übergebene "
                "Last ist stetig und trägt die integrierte Maxwell-Kraft exakt. Die blaue "
                "Treppe ist die konservative Buchhaltung und war in P3b die Ursache dafür, "
                "dass der Integrator seine Genauigkeit nie erreichte.", 0.048)
    provenance(fig, m, NOT_MODELLED, y=0.008)
    fig.tight_layout(rect=[0, 0.090, 1, 0.950])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
# fig 2 -- force conservation and the measured order
# ==========================================================================
def figure_error(d, m, out):
    rs = rows(os.path.join(d, "projection_audit.csv"))
    if not rs:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 6.4))
    fig.suptitle(TITLE + "\nAbb. 2 – Krafterhaltung der Projektion und die gemessene "
                 "Ordnung gegen die geschlossenen Formen", fontsize=12.5, y=0.98)

    ax = axes[0]
    smooth = [r for r in rs if r["case"] == "glatt"]
    n = col(smooth, "segments")
    e = col(smooth, "err_segment_vs_analytic")
    ax.loglog(n, e, "o-", color="#1f77b4", label="glatte Last, gemessen")
    good = np.isfinite(e) & (e > 0)
    if good.sum() >= 2:
        ref = e[good][0] * (n[good] / n[good][0]) ** -2.0
        ax.loglog(n[good], ref, "--", color="#888888", label="2. Ordnung")
    for beta in sorted({f(r["beta"]) for r in rs if r["case"] == "singulaer"}):
        sr = [r for r in rs if r["case"] == "singulaer" and f(r["beta"]) == beta]
        nn = col(sr, "segments")
        ee = col(sr, "err_segment_vs_analytic")
        ax.loglog(nn, ee, "s-", ms=4, label=f"β = {beta:+.2f}, gemessen")
        g = np.isfinite(ee) & (ee > 0)
        if g.sum() >= 2:
            ax.loglog(nn[g], ee[g][0] * (nn[g] / nn[g][0]) ** -(1.0 + beta), ":",
                      color="#888888")
    ax.set_xlabel("Segmente auf der Oberfläche")
    ax.set_ylabel("|F_Projektion − F_exakt| / F_exakt")
    ax.set_title("Segmentprojektion gegen die geschlossene Form", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    ax = axes[1]
    lbl = []
    for i, r in enumerate(rs):
        tag = r["case"] + ("" if r["case"] != "singulaer" else f" β={f(r['beta']):+.2f}")
        lbl.append(tag)
    x = np.arange(len(rs))
    ax.semilogy(x, np.maximum(col(rs, "err_handed_reconstructed"), 1e-18), "o", ms=3.5,
                color="#111111", label="übergebene Last gegen A′ (exakt erhalten)")
    ax.semilogy(x, np.maximum(col(rs, "err_handed_true"), 1e-18), "s", ms=3.5,
                color="#ff7f0e", label="übergebene Last gegen 2πr ds")
    ax.semilogy(x, np.maximum(col(rs, "err_bin_vs_segment"), 1e-18), "^", ms=3.5,
                color="#2ca02c", label="Bins gegen Segmente")
    ax.axhline(1e-13, color="#888888", ls="--", lw=0.9)
    ax.set_xlabel("Fall (Zeile in projection_audit.csv)")
    ax.set_ylabel("relativer Kraftfehler")
    ax.set_title("Krafterhaltung durch die ganze Kette", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    ax = axes[2]
    ax.semilogy(x, col(rs, "jump_handed_Pa") / np.maximum(col(rs, "load_span_Pa"), 1e-300),
                "o", ms=3.5, color="#111111", label="übergebene Last")
    ax.semilogy(x, col(rs, "jump_stair_Pa") / np.maximum(col(rs, "load_span_Pa"), 1e-300),
                "s", ms=3.5, color="#1f77b4", label="Segmenttreppe")
    ax.set_xlabel("Fall")
    ax.set_ylabel("größter Sprung über einen Binrand / Spannweite")
    ax.set_title("Stetigkeit: der Sprung der übergebenen Last", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    caveat(fig, "Links: die Segmentprojektion ist bei glatter Last zweiter Ordnung und bei "
                "p = C d^β genau von der Ordnung 1+β – das ist die Rate, die eine integrable "
                "Randsingularität zulässt, und keine Schwäche der Projektion. Mitte: die "
                "übergebene Last trägt die integrierte Kraft gegen das rekonstruierte "
                "Flächenmaß auf Maschinengenauigkeit; gegen 2πr ds gilt sie zur "
                "Interpolationsordnung. Beide Zahlen stehen hier, nicht nur die günstigere.",
            0.045)
    provenance(fig, m, NOT_MODELLED)
    fig.tight_layout(rect=[0, 0.105, 1, 0.925])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
# fig 3 -- kTolExclusion, resolved
# ==========================================================================
def figure_exclusion(d, m, out):
    rs = rows(os.path.join(d, "exclusion_criterion.csv"))
    if not rs:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.4))
    fig.suptitle(TITLE + "\nAbb. 3 – Warum edge_gate::kTolExclusion nicht eingehalten "
                 "werden KANN", fontsize=12.5, y=0.98)

    ax = axes[0]
    for beta in sorted({f(r["beta"]) for r in rs}):
        sr = [r for r in rs if f(r["beta"]) == beta]
        n = col(sr, "segments")
        meas = col(sr, "measured_change")
        lim = col(sr, "closed_form_limit")
        p = ax.semilogx(n, meas, "o-", ms=4, label=f"β = {beta:+.2f}, gemessen")
        ax.axhline(lim[0], ls="--", lw=1.0, color=p[0].get_color())
    ax.axhline(0.05, color="#d62728", lw=1.6,
               label="Schranke kTolExclusion = 5 %")
    ax.set_xlabel("Segmente auf der Oberfläche (Verfeinerung)")
    ax.set_ylabel("[F(d₀/2) − F(d₀)] / F(d₀/2)")
    ax.set_title("die geprüfte Größe läuft gegen einen Grenzwert ≠ 0", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    ax = axes[1]
    for beta in sorted({f(r["beta"]) for r in rs}):
        sr = [r for r in rs if f(r["beta"]) == beta]
        ax.loglog(col(sr, "segments"), col(sr, "discretisation_error"), "s-", ms=4,
                  label=f"β = {beta:+.2f}")
    ax.set_xlabel("Segmente auf der Oberfläche")
    ax.set_ylabel("|F_Projektion − F_exakt| / F_exakt")
    ax.set_title("der Diskretisierungsfehler DERSELBEN Last fällt weiter", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2)

    agr = f(m.get("exclusion_worst_agreement", "nan"))
    caveat(fig, "Links steht die Größe, die kTolExclusion prüft: sie enthält die Netzweite "
                "gar nicht und läuft gegen einen geschlossenen Grenzwert, der nur von β und "
                f"d₀ abhängt (gestrichelt; größte Abweichung Messung/Formel {agr:.1e}). "
                "Rechts fällt der echte Diskretisierungsfehler derselben Last weiter. Damit "
                "ist das Kriterium widerlegt und nicht getauscht: es fordert von einer "
                "geometrischen Größe, dass sie verschwindet. Die Schranke bleibt im Code "
                "deklariert und wird als nicht eingehalten berichtet.", 0.062)
    provenance(fig, m, NOT_MODELLED, y=0.008)
    fig.tight_layout(rect=[0, 0.135, 1, 0.925])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
# fig 4 -- the gate, all three views at once
# ==========================================================================
def figure_gate(d, m, out):
    lv = rows(os.path.join(d, "gate_levels.csv"))
    vd = rows(os.path.join(d, "gate_verdict.csv"))
    if not lv:
        return None
    shapes = sorted({r["shape"] for r in lv})
    fig, axes = plt.subplots(1, len(shapes), figsize=(7.2 * len(shapes), 6.6), squeeze=False)
    fig.suptitle(TITLE + "\nAbb. 4 – Kanten-Gate: direkte Kraftfolge und die beiden "
                 "Extrapolationen desselben Grenzwertes", fontsize=12.5, y=0.98)

    for k, sh in enumerate(shapes):
        ax = axes[0][k]
        sr = [r for r in lv if r["shape"] == sh]
        lev = col(sr, "level")
        tot = col(sr, "total_force_N") * 1e9
        ax.plot(lev, tot, "o-", color="#111111", lw=1.8,
                label="Gesamtkraft, direkt je Netzstufe")
        for name, cname, c in (("force_beyond_coarse_N", "d ≥ 0,100 a", "#1f77b4"),
                               ("force_beyond_mid_N", "d ≥ 0,050 a", "#2ca02c"),
                               ("force_beyond_fine_N", "d ≥ 0,025 a", "#ff7f0e")):
            ax.plot(lev, col(sr, name) * 1e9, "s--", ms=4, color=c, lw=1.1,
                    label=f"Kraft jenseits {cname}")
        v = next((r for r in vd if r["shape"] == sh), None)
        if v:
            fm = f(v["limit_force_mesh_N"]) * 1e9
            fe = f(v["limit_force_exclusion_N"]) * 1e9
            ax.axhline(fm, color="#9467bd", lw=1.6,
                       label=f"Extrapolation über die Netzstufen: {fm:.4f} nN")
            ax.axhline(fe, color="#8c564b", lw=1.6, ls="-.",
                       label=f"Extrapolation über die Ausschlussdistanz: {fe:.4f} nN")
            ax.set_title(f"{sh} (Π = {f(v['Pi']):+.2f}) – {v['verdict']}, "
                         f"β = {f(v['fitted_exponent']):+.3f}", fontsize=10)
        else:
            ax.set_title(sh, fontsize=10)
        ax.set_xlabel("Netzstufe (Größenfeldfaktor 2^(−Stufe/2))")
        ax.set_ylabel("integrierte Maxwell-Kraft [nN]")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.0, loc="lower right")

    caveat(fig, "Alle drei Sichten stehen nebeneinander: die direkte Kraftfolge über die "
                "Netzstufen, die Extrapolation daraus, und die aus den Ausschlussdistanzen. "
                "Die beiden Extrapolationen teilen keine Daten – dass sie sich treffen, ist "
                "der Grund, aus dem die Kantenlast als summierbar gilt. Der PUNKTWEISE "
                "Kantenwert wächst dagegen mit jeder Verfeinerung und wird nirgends "
                "verwendet.", 0.062)
    provenance(fig, m, NOT_MODELLED, y=0.008)
    fig.tight_layout(rect=[0, 0.135, 1, 0.925])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ==========================================================================
# fig 5 -- the coupled mesh study
# ==========================================================================
def figure_coupled(d, m, out):
    cs = rows(os.path.join(d, "coupled_convergence.csv"))
    rich = rows(os.path.join(d, "coupled_richardson.csv"))
    if not cs:
        return None
    volts = sorted({f(r["V_emitter_V"]) for r in cs})
    quantities = [("h_over_a", "h/a", 1.0), ("E_edge_far_V_per_m", "E_n kantenfern [MV/m]", 1e-6),
                  ("total_force_N", "integrierte Maxwell-Kraft [nN]", 1e9)]
    fig, axes = plt.subplots(len(quantities), len(volts),
                             figsize=(6.6 * len(volts), 4.0 * len(quantities)), squeeze=False)
    target = f(m.get("discretization_target", "0.01"))

    any_missed = False
    for j, V in enumerate(volts):
        for i, (key, label, scale) in enumerate(quantities):
            ax = axes[i][j]
            sr = [r for r in cs if f(r["V_emitter_V"]) == V]
            lev = col(sr, "level")
            val = col(sr, key) * scale
            ok = np.isfinite(val)
            ax.plot(lev[ok], val[ok], "o-", color="#111111", lw=1.8)
            # A failed point is a GAP, marked where it belongs -- never a zero.
            for r in sr:
                if not np.isfinite(f(r[key])):
                    ax.axvline(f(r["level"]), color="#d62728", lw=1.2, ls=":")
                    ax.annotate(r["status"], (f(r["level"]), np.nanmean(val[ok])),
                                rotation=90, fontsize=6.5, color="#d62728",
                                ha="right", va="center")
            rr = next((x for x in rich if f(x["V_emitter_V"]) == V and x["quantity"] ==
                       key.replace("_V_per_m", "").replace("_N", "") and
                       x.get("window", "fein") == "fein"), None)
            note = ""
            if rr:
                err = f(rr["relative_error_finest"])
                ext = f(rr["extrapolated"])
                if np.isfinite(ext):
                    ax.axhline(ext * scale, color="#9467bd", lw=1.3, ls="--")
                if np.isfinite(err):
                    note = (f"Ordnung {f(rr['observed_order']):.2f}, "
                            f"geschätzter Fehler {err * 100:.2f} %")
                    if err >= target:
                        any_missed = True
                else:
                    note = (f"keine Ordnung beobachtbar; bloße Änderung "
                            f"{f(rr['last_relative_change']) * 100:.2f} %")
                    any_missed = True
                note += f" – {rr['verdict']}"
            ax.set_title(f"{int(V)} V – {label}\n{note}", fontsize=9)
            ax.set_xlabel("Netzstufe")
            ax.set_ylabel(label)
            ax.grid(alpha=0.25)

    head = ("Abb. 5 – gekoppelte Netzstudie mit Richardson-Extrapolation"
            if not any_missed else
            "Abb. 5 – gekoppelte Netzstudie: das 1-%-Ziel ist NICHT erreicht")
    fig.suptitle(TITLE + "\n" + head, fontsize=12.5, y=0.975)
    caveat(fig, "Diese Abbildung heißt ausdrücklich NICHT „Konvergenz“, solange der vorab "
                "festgelegte geschätzte Diskretisierungsfehler von 1 % verfehlt wird. Die "
                "violette Linie ist der extrapolierte Grenzwert, wo eine Ordnung überhaupt "
                "beobachtbar war; wo keine beobachtbar ist, steht nur die bloße Änderung "
                "zwischen den beiden feinsten Stufen, und die ist keine Fehlerschranke. "
                "Eine gescheiterte Rechnung erscheint als Lücke mit ihrem Status, niemals "
                "als Wert null."
            if any_missed else
            "Die violette Linie ist der extrapolierte Grenzwert. Eine gescheiterte Rechnung "
            "erscheint als Lücke mit ihrem Status, niemals als Wert null.", 0.038)
    provenance(fig, m, NOT_MODELLED, y=0.006)
    fig.tight_layout(rect=[0, 0.085, 1, 0.930])
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------- main
def main(d):
    global STAMP
    m = meta(d)
    root = pv.repo_root_of(__file__)
    rel = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/")
    state = pv.repo_state(root, ignore=[rel])
    STAMP = pv.stamp(state, m.get("commit"))
    pv.write_provenance(d, state, m.get("commit"))

    made = [figure_loads(d, m, os.path.join(d, "fig1_load_projection.png")),
            figure_error(d, m, os.path.join(d, "fig2_projection_error.png")),
            figure_exclusion(d, m, os.path.join(d, "fig3_exclusion.png")),
            figure_gate(d, m, os.path.join(d, "fig4_gate.png")),
            figure_coupled(d, m, os.path.join(d, "fig5_coupled.png"))]
    for x in made:
        if x:
            print(f"{x}  {os.path.getsize(x)} Byte")
        else:
            print("Eine Abbildung entfiel: die zugehoerigen Daten fehlen.")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
