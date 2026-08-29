#!/usr/bin/env python3
"""Figures of the P3a static capillary meniscus.

    ./build/es_capillary examples/device_p1.cfg examples/capillary_p3a.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_capillary.py results/<dir>

Reads only the CSVs written by apps/es_capillary.cpp.  Nothing here solves,
integrates or computes a geometry; the analytic cap drawn in figure 2 is the
column the C++ run wrote next to the numerical profile.

The commit named on every figure is measured from git at drawing time by
python/provenance.py, not taken from an argument, and a dirty working tree is
stamped DIRTY.

Produces
  <dir>/fig1_geometry_and_bc.png     geometry, symmetry axis, pinned contact
                                     line, pressure definition, normal direction
  <dir>/fig2_profiles.png            numerical profiles for several positive and
                                     negative Pi with the closed-form cap on top
  <dir>/fig3_convergence.png         mesh convergence of profile, apex height,
                                     area, volume and Young-Laplace residual
  <dir>/fig4_apex_and_curvature.png  apex height and curvature over Pi, with the
                                     boundary of the representable range
  <dir>/fig5_liquid_example.png      the marked material example (EMI-BF4)
  <dir>/figures_provenance.txt       what the stamp said, and whether it is
                                     releasable
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
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance as pv

TITLE = ("P3a: statisches Kapillargleichgewicht an der Austrittskante – "
         "KEIN elektrisches Feld, KEINE Emission")

NOT_MODELLED = ("Nicht enthalten: elektrisches Feld, Maxwell-Druck, Kopplung an den "
                "Elektrostatiksolver, Emission, Raumladung, Strömung, Zeitabhängigkeit, "
                "Stabilitätsaussage, Taylor-Kegel, Cone-Jet, Schwerkraft.")


# --------------------------------------------------------------------- input
def rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def meta(d):
    out = {}
    with open(os.path.join(d, "meta.txt"), encoding="utf-8") as fh:
        for line in fh:
            if "=" in line:
                k, v = line.rstrip("\n").split("=", 1)
                out[k] = v
    return out


def keyvalues(path):
    return {r["key"]: r["value"] for r in rows(path)}


def regions(d):
    out = defaultdict(list)
    for r in rows(os.path.join(d, "regions.csv")):
        out[(r["region"], int(r["loop"]))].append((float(r["r_m"]), float(r["z_m"])))
    return out


def profiles(d, path="profiles.csv", key="variant"):
    out = defaultdict(lambda: {"r": [], "z": [], "za": [], "dz": [], "dn": [],
                               "Pi": 0.0, "dp": 0.0})
    for r in rows(os.path.join(d, path)):
        e = out[r[key]]
        e["r"].append(float(r["r_m"]))
        e["z"].append(float(r["z_m"]))
        if "z_analytic_m" in r:
            e["za"].append(float(r["z_analytic_m"]))
        # The difference is written by the C++ run at full precision: taken as a
        # difference of the printed coordinates it would be swamped by their
        # output precision, and every curve would collapse onto zero.
        if "dz_m" in r:
            e["dz"].append(float(r["dz_m"]))
            e["dn"].append(float(r["dn_m"]))
        e["Pi"] = float(r["Pi"])
        if "delta_p_Pa" in r:
            e["dp"] = float(r["delta_p_Pa"])
    for e in out.values():
        for k in ("r", "z", "za", "dz", "dn"):
            e[k] = np.array(e[k])
    return out


def floats(rs, col):
    out = []
    for r in rs:
        try:
            out.append(float(r[col]))
        except (TypeError, ValueError, KeyError):
            out.append(np.nan)
    return np.array(out)


# -------------------------------------------------------------------- helpers
def pi_colour(pi, pmax=2.0):
    """Blue for a meniscus drawn in, red for one bulging out, black for flat.

    A diverging map with a white centre would make the flat surface invisible,
    which is the one case that must be unmistakable.
    """
    t = 0.42 + 0.53 * min(1.0, abs(pi) / pmax)
    if pi > 0:
        return plt.get_cmap("Reds")(t)
    if pi < 0:
        return plt.get_cmap("Blues")(t)
    return "#111111"


def provenance(fig, m, stamp, extra="", y=0.010, width=176):
    dirty = pv.DIRTY_MARK in stamp or stamp == pv.NO_VERSION
    head = (f"Commit {stamp}   |   {m.get('app','es_capillary')}   |   "
            f"Konfiguration {m.get('config','?')}   |   "
            f"a = $\\phi_2/2$ = {float(m['contact_radius_m'])*1e6:.2f} µm, "
            f"$\\gamma$ = {float(m['surface_tension_N_per_m'])*1e3:.2f} mN/m "
            f"(Status {m.get('liquid_status','?')}), "
            f"$\\gamma/a$ = {float(m['capillary_pressure_scale_Pa']):.0f} Pa, "
            f"$|\\Delta p| \\leq$ {float(m['delta_p_max_Pa']):.0f} Pa, "
            f"Bo = {float(m['bond_number']):.2e}")
    text = head + ("\n" + extra if extra else "")
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width) if len(text) > width else [text]),
             ha="center", va="bottom", fontsize=7.0,
             color="#8a1b1b" if dirty else "#444444")
    if dirty:
        fig.text(0.012, 0.992, f"NICHT FREIGEGEBEN – {stamp}", ha="left", va="top",
                 fontsize=8, color="#8a1b1b", weight="bold")


def caveat(fig, text, y, width=172, color="#7a3b00"):
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="bottom",
             fontsize=7.4, color=color)


def heading(fig, m, subtitle, banner_colour="#1f3d6b", extra_banner=""):
    fig.suptitle(subtitle + "\n" + TITLE, fontsize=11.6, y=0.983, va="top")
    fig.text(0.5, 0.918,
             "statisches Kapillargleichgewicht  ·  kein elektrisches Feld  ·  keine Emission"
             f"  ·  Stoffstatus: {m.get('liquid_status','?')}" + extra_banner,
             ha="center", va="top", fontsize=8.6, color=banner_colour)


# ------------------------------------------------------------------- figure 1
def figure_geometry(d, m, out):
    """Geometry and boundary conditions, at two scales."""
    reg = regions(d)
    a = float(m["contact_radius_m"])
    pr = profiles(d)
    s = 1e6  # micrometres
    au = a * s

    fig = plt.figure(figsize=(15.4, 7.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.45], left=0.05, right=0.985,
                          top=0.830, bottom=0.215, wspace=0.14)

    # --- left: the device the edge belongs to --------------------------------
    ax = fig.add_subplot(gs[0, 0])
    colours = {"liquid": "#bcd9ef", "emitter_solid": "#d8d2c4"}
    for (name, loop), pts in sorted(reg.items()):
        if name not in colours or loop != 0:
            continue
        ax.add_patch(Polygon(np.array(pts) * s, closed=True, facecolor=colours[name],
                             edgecolor="#333333", lw=0.9, zorder=2))
    ax.axvline(0, color="#7a1fa2", lw=1.6, ls=(0, (6, 3)), zorder=4)
    ax.text(-1.2, -34, "Symmetrieachse $r=0$", color="#7a1fa2", fontsize=8.2, rotation=90,
            va="center", ha="center")
    ax.plot([au], [0.0], "o", ms=7, mfc="#d62728", mec="#5c0d0d", zorder=6)
    ax.annotate("gepinnte Kontaktlinie\n$r=\\phi_2/2$, $z=0$", xy=(au, 0),
                xytext=(9, 14), fontsize=8.2, color="#5c0d0d",
                arrowprops=dict(arrowstyle="->", lw=1.0, color="#5c0d0d"))
    ax.add_patch(Rectangle((0, -1.15 * au), 2.05 * au, 2.35 * au, fill=False,
                           ec="#1f3d6b", lw=1.3, ls="--", zorder=7))
    ax.text(2.25 * au, -0.30 * au, "Ausschnitt\nrechts", fontsize=8.0, color="#1f3d6b",
            va="center", ha="left")
    ax.set_xlim(-3.5, 24)
    ax.set_ylim(-66, 22)
    ax.set_aspect("equal")
    ax.set_xlabel("r [µm]")
    ax.set_ylabel("z [µm]")
    ax.set_title("Emitterspitze mit Bohrung\n(P1-Geometrie, unverändert)", fontsize=9.4)
    ax.legend(handles=[
        Line2D([], [], marker="s", ls="", mfc=colours["liquid"], mec="#333", ms=10,
               label="Flüssigkeit (Bohrung als\ngefüllt vorausgesetzt)"),
        Line2D([], [], marker="s", ls="", mfc=colours["emitter_solid"], mec="#333", ms=10,
               label="Emitterkörper"),
    ], loc="lower left", fontsize=7.4, framealpha=0.95)

    # --- right: the boundary-value problem at the edge ------------------------
    ax = fig.add_subplot(gs[0, 1])
    zlo, zhi = -1.48 * au, 1.20 * au
    ax.add_patch(Rectangle((0, zlo), au, -zlo, facecolor="#dceaf6", ec="none", zorder=0))
    ax.add_patch(Rectangle((au, zlo), 1.25 * au, -zlo, facecolor="#d8d2c4", ec="none",
                           zorder=1))
    ax.plot([au, au], [zlo, 0], color="#333", lw=1.5, zorder=3)
    ax.plot([au, 2.25 * au], [0, 0], color="#333", lw=1.5, zorder=3)
    ax.text(1.60 * au, -0.28 * au, "Bohrungswand\n(Emitterkörper)", fontsize=8.2,
            ha="center", color="#4a4034")
    ax.text(0.62 * au, -0.30 * au, "Flüssigkeit", fontsize=8.2, ha="center", color="#2c5f8a")

    ax.axvline(0, color="#7a1fa2", lw=1.7, ls=(0, (6, 3)), zorder=4)
    ax.text(0.08 * au, -0.93 * au,
            "$r=0$: Symmetrie,  $\\psi(0)=0$,\n"
            "$\\mathrm{d}\\psi/\\mathrm{d}s|_\\mathrm{Apex}=\\kappa/2$"
            "  (analytischer Grenzwert)",
            color="#7a1fa2", fontsize=8.0, va="center")

    for tag, e in sorted(pr.items(), key=lambda kv: kv[1]["Pi"]):
        if abs(e["Pi"]) not in (0.0, 1.5):
            continue
        ax.plot(e["r"] * s, e["z"] * s, lw=2.4, color=pi_colour(e["Pi"]), zorder=5,
                label=f"$\\Pi$ = {e['Pi']:+.1f}  ($\\Delta p$ = {e['dp']:+.0f} Pa)")

    bulge = [e for e in pr.values() if abs(e["Pi"] - 1.5) < 1e-9][0]
    for i in np.linspace(0, len(bulge["r"]) - 1, 9).astype(int)[1:-1]:
        j = min(i + 1, len(bulge["r"]) - 1)
        tr, tz = bulge["r"][j] - bulge["r"][i], bulge["z"][j] - bulge["z"][i]
        n = np.hypot(tr, tz)
        nr, nz = -tz / n, tr / n           # tangent rotated by +90 degrees
        if nz < 0:
            nr, nz = -nr, -nz
        ax.arrow(bulge["r"][i] * s, bulge["z"][i] * s, 0.16 * au * nr, 0.16 * au * nz,
                 head_width=0.045 * au, head_length=0.06 * au, fc="#0f6b3d", ec="#0f6b3d",
                 lw=1.0, zorder=6, length_includes_head=True)
    ax.annotate("äußere Normale $\\mathbf{n}=(\\sin\\psi,\\ \\cos\\psi)$,\n"
                "aus der Flüssigkeit ins Vakuum",
                xy=(0.62 * au, 0.98 * au), xytext=(1.12 * au, 1.05 * au), fontsize=8.2,
                color="#0f6b3d", va="center",
                arrowprops=dict(arrowstyle="->", lw=1.0, color="#0f6b3d"))

    ax.plot([au], [0.0], "o", ms=9, mfc="#d62728", mec="#5c0d0d", zorder=8)
    ax.annotate("gepinnte Kontaktlinie – KEIN Kontaktwinkel\n"
                "(Pinning und Young-Winkel schließen einander aus)",
                xy=(au, 0), xytext=(1.10 * au, 0.42 * au), fontsize=8.2, color="#5c0d0d",
                va="center",
                arrowprops=dict(arrowstyle="->", lw=1.0, color="#5c0d0d"),
                bbox=dict(fc="#ffffff", ec="#e7b7b7", lw=0.8, alpha=0.95))

    ax.annotate("", xy=(au, -1.30 * au), xytext=(0, -1.30 * au),
                arrowprops=dict(arrowstyle="<->", lw=1.2, color="#1f3d6b"))
    ax.text(0.5 * au, -1.27 * au, f"$a=\\phi_2/2$ = {au:.2f} µm", fontsize=8.6,
            ha="center", color="#1f3d6b")

    ax.text(1.01 * au, -1.40 * au,
            "$\\Delta p_\\mathrm{exit}=p_\\mathrm{fl}-p_\\mathrm{vak}$ an $z=0$\n"
            "$>0$: Wölbung nach $+z$ (zum Extraktor)\n"
            "$=0$: exakt eben\n"
            "$<0$: in die Bohrung gezogen",
            fontsize=7.9, color="#111111", va="bottom",
            bbox=dict(fc="#ffffff", ec="#999999", lw=0.8, alpha=0.95))

    ax.set_xlim(-0.12 * au, 2.20 * au)
    ax.set_ylim(zlo, zhi)
    ax.set_aspect("equal")
    ax.set_xlabel("r [µm]")
    ax.set_ylabel("z [µm]")
    ax.set_title("Randwertproblem an der Austrittskante:   "
                 "$\\gamma\\,(\\mathrm{d}\\psi/\\mathrm{d}s + \\sin\\psi/r)"
                 " = \\Delta p_\\mathrm{exit}$", fontsize=9.8)
    ax.legend(loc="upper left", fontsize=8.0, framealpha=0.95)

    heading(fig, m, "Geometrie und Randbedingungen")
    caveat(fig, NOT_MODELLED + " Die Bohrung gilt als bis zur Kante gefüllt; der Zustand "
                "stromauf geht ausschließlich über delta_p_exit ein. Die P2c-Vorratsgeometrie "
                "ist unverändert und wird hier nicht gelöst.", y=0.070)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 2
def figure_profiles(d, m, out):
    pr = profiles(d)
    summ = rows(os.path.join(d, "profile_summary.csv"))
    a = float(m["contact_radius_m"])
    s = 1e6
    au = a * s

    fig = plt.figure(figsize=(15.4, 7.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], left=0.05, right=0.985,
                          top=0.830, bottom=0.235, wspace=0.20)

    ax = fig.add_subplot(gs[0, 0])
    for tag, e in sorted(pr.items(), key=lambda kv: kv[1]["Pi"]):
        ax.plot(e["r"] * s, e["z"] * s, lw=2.0, color=pi_colour(e["Pi"]),
                label=f"$\\Pi$={e['Pi']:+.2f},  $\\Delta p$={e['dp']:+7.0f} Pa")
        step = max(1, len(e["r"]) // 24)
        ax.plot(e["r"][::step] * s, e["za"][::step] * s, ls="", marker="o", ms=3.6,
                mfc="none", mec="#111111", mew=0.7, zorder=5)
    ax.plot([au], [0.0], "o", ms=8, mfc="#d62728", mec="#5c0d0d", zorder=8)
    ax.axvline(0, color="#7a1fa2", lw=1.2, ls=(0, (6, 3)))
    ax.set_xlim(-0.3, 2.05 * au)
    ax.set_xlabel("r [µm]")
    ax.set_ylabel("z [µm]")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_title("Numerische Profile (Linien), analytische Kugelkappe (Kreise)",
                 fontsize=9.8)
    ax.legend(loc="center right", fontsize=7.4, framealpha=0.95,
              title="Apex $\\rightarrow$ gepinnte Kante", title_fontsize=7.6)

    ax = fig.add_subplot(gs[0, 1])
    # The differences come from the dz_m / dn_m columns, which the C++ run writes
    # at full double precision.  Taken as a difference of the PRINTED coordinates
    # they would be swamped by their output precision and every curve would
    # collapse onto zero -- an artefact of the file format, not a result.
    floor = 1e-19
    for tag, e in sorted(pr.items(), key=lambda kv: kv[1]["Pi"]):
        if e["Pi"] == 0.0:
            continue
        ax.semilogy(e["r"] / a, np.maximum(np.abs(e["dz"]) / a, floor), lw=1.2,
                    color=pi_colour(e["Pi"]), label=f"$\\Pi$={e['Pi']:+.2f}")
        ax.semilogy(e["r"] / a, np.maximum(np.abs(e["dn"]) / a, floor), lw=0.9, ls=":",
                    color=pi_colour(e["Pi"]))
    flat = [e for e in pr.values() if e["Pi"] == 0.0]
    if flat:
        ax.semilogy(flat[0]["r"] / a, np.maximum(np.abs(flat[0]["dz"]) / a, floor), lw=2.4,
                    color="#111111",
                    label="$\\Pi$=0: identisch 0 (auf $10^{-19}$ gesetzt)")
    ax.set_xlabel("$r/a$")
    ax.set_ylabel("Abweichung von der Kugelkappe / a\n"
                  "(durchgezogen: $|\\Delta z|$, gepunktet: Normalabstand)")
    ax.set_ylim(floor, 1e-11)
    ax.grid(alpha=0.25, which="both")
    ax.set_title("Punktweise Abweichung vom geschlossenen Ergebnis", fontsize=9.8)
    ax.text(0.03, 0.06,
            "Die Kurven zu $+\\Pi$ und $-\\Pi$ liegen exakt übereinander:\n"
            "das Problem ist im Vorzeichen des Drucks spiegelsymmetrisch.",
            transform=ax.transAxes, fontsize=7.4, color="#444444", va="bottom")
    ax.legend(fontsize=7.2, ncol=2, framealpha=0.95, loc="lower right")

    ok = [r for r in summ if r["status"] == "Solved"]
    worst = max(float(r["err_profile_normal"]) for r in ok)
    worst_h = max(float(r["err_h_rel"]) for r in ok)
    worst_v = max(float(r["err_V_rel"]) for r in ok)

    heading(fig, m, "Meniskusprofile gegen die geschlossene Kugelkappe")
    caveat(fig,
           "Ohne Feld und ohne Schwerkraft ist die Lösung eine Kugelkappe mit konstanter "
           "mittlerer Krümmung kappa = delta_p/gamma. Der Löser wertet diese geschlossene Form "
           "nie aus: er integriert die Young-Laplace-Gleichung in Bogenlänge und sucht die "
           "Bogenlänge, bei der die Meridiankurve den Pinningradius erreicht. Größter "
           f"Normalabstand über alle gezeigten Drücke: {worst:.1e}·a; größter relativer Fehler "
           f"der Apexhöhe {worst_h:.1e}, des Rotationsvolumens {worst_v:.1e}. " + NOT_MODELLED,
           y=0.070, color="#0f6b3d")
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 3
def figure_convergence(d, m, out):
    rs = rows(os.path.join(d, "convergence.csv"))
    by_pi = defaultdict(list)
    for r in rs:
        by_pi[float(r["Pi"])].append(r)
    res_by = defaultdict(lambda: ([], []))
    for r in rows(os.path.join(d, "residuals.csv")):
        t, y = res_by[r["variant"]]
        t.append(float(r["s_over_arclength"]))
        y.append(float(r["residual"]))

    fig = plt.figure(figsize=(14.6, 8.6))
    gs = fig.add_gridspec(2, 2, left=0.06, right=0.985, top=0.830, bottom=0.245,
                          wspace=0.20, hspace=0.32)
    axes = [fig.add_subplot(gs[i, j]) for i in (0, 1) for j in (0, 1)]

    panels = [("err_profile_normal", "Profilfehler (Normalabstand) / a", 4),
              ("err_h_rel", "relativer Fehler der Apexhöhe", 4),
              ("err_V_rel", "rel. Fehler: Volumen (o), Fläche (x)", 4),
              ("residual_max", "max. Young-Laplace-Residuum $\\cdot\\,a/\\gamma$", 2)]
    for ax, (col, label, order) in zip(axes, panels):
        for pi in sorted(by_pi):
            g = by_pi[pi]
            n = floats(g, "n_intervals")
            ax.loglog(n, floats(g, col), "o-", ms=4.5, lw=1.4, color=pi_colour(pi),
                      label=f"$\\Pi$={pi:+.1f}")
            if col == "err_V_rel":
                ax.loglog(n, floats(g, "err_A_rel"), "x--", ms=5, lw=1.0,
                          color=pi_colour(pi))
        pi0 = sorted(by_pi)[0]
        n0 = floats(by_pi[pi0], "n_intervals")
        y0 = floats(by_pi[pi0], col)[0]
        ax.loglog(n0, y0 * (n0[0] / n0) ** order, ls=":", color="#666666", lw=1.5,
                  label=f"Steigung $n^{{-{order}}}$")
        ax.set_xlabel("Intervalle auf der Meridiankurve (erzwungen)")
        ax.set_ylabel(label, fontsize=8.6)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=7.0, ncol=2, framealpha=0.95, loc="lower left")

    ins = axes[3].inset_axes([0.55, 0.50, 0.42, 0.34])
    for tag in sorted(res_by):
        t, y = res_by[tag]
        ins.plot(t, y, lw=0.9)
    ins.set_xlabel("$s/s_\\mathrm{Kontakt}$", fontsize=6.4, labelpad=1)
    ins.set_ylabel("Residuum", fontsize=6.4, labelpad=1)
    ins.tick_params(labelsize=5.8)
    ins.set_title("Residuum entlang der Oberfläche,\nautomatische Auflösung", fontsize=6.6)
    ins.grid(alpha=0.25)

    heading(fig, m, "Netzkonvergenz: Profil, Apexhöhe, Fläche, Volumen, YL-Residuum")
    caveat(fig,
           "Die Intervallzahl ist KEINE Benutzereingabe: im Normalbetrieb wählt der Löser sie "
           f"aus der geforderten Genauigkeit ({float(m['requested_accuracy']):.0e}) durch "
           "Verdopplung, bis die Änderung zur halben Auflösung darunter liegt – hier ist sie "
           "für die Studie erzwungen. Das Residuum wird allein aus den Knotenkoordinaten "
           "gebildet (Tangentenwinkel aus den Sehnen, Meridiankrümmung aus deren Drehung) und "
           "teilt keinen Code mit dem Integrator; es konvergiert deshalb zweiter statt vierter "
           "Ordnung. Die Kurven zu +Pi und -Pi liegen exakt übereinander: das Problem ist "
           "im Vorzeichen des Drucks spiegelsymmetrisch. " + NOT_MODELLED, y=0.070)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 4
def figure_apex_and_curvature(d, m, out):
    rs = rows(os.path.join(d, "sweep.csv"))
    pi = floats(rs, "Pi")
    ok = np.array([r["status"] == "Solved" for r in rs])
    h = floats(rs, "h_over_a_num")
    ha = floats(rs, "h_over_a_ana")
    kn = floats(rs, "kappa_num_1_per_m")
    ka = floats(rs, "kappa_ana_1_per_m")
    psi = floats(rs, "psi_contact_deg")
    a = float(m["contact_radius_m"])
    pmax = float(m["pi_max"])
    rng = rows(os.path.join(d, "range.csv"))

    fig = plt.figure(figsize=(15.6, 6.6))
    gs = fig.add_gridspec(1, 3, left=0.055, right=0.985, top=0.820, bottom=0.290,
                          wspace=0.26)
    axes = [fig.add_subplot(gs[0, j]) for j in range(3)]

    def mark_range(ax):
        ax.axvspan(pmax, pi.max(), color="#f6d5d5", zorder=0)
        ax.axvspan(pi.min(), -pmax, color="#f6d5d5", zorder=0)
        for x in (-pmax, pmax):
            ax.axvline(x, color="#8a1b1b", lw=1.2, ls="--", zorder=1)
        ax.set_xlim(pi.min(), pi.max())

    ax = axes[0]
    mark_range(ax)
    ax.plot(pi[ok], ha[ok], lw=3.6, color="#c9c9c9", label="analytisch (Kugelkappe)")
    ax.plot(pi[ok], h[ok], lw=1.4, color="#1f3d6b", label="numerisch")
    ax.plot([-pmax, pmax], [-1, 1], "o", ms=7, mfc="#8a1b1b", mec="#4a0a0a",
            label="Halbkugel $|\\Pi|=2$: $h=\\pm a$")
    ax.set_xlabel("$\\Pi = \\Delta p\\,a/\\gamma$")
    ax.set_ylabel("Apexhöhe $h/a$")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.8, loc="upper left")
    ax.set_title("Apexhöhe über dem dimensionslosen Druck", fontsize=9.6)

    ax = axes[1]
    mark_range(ax)
    ax.plot(pi[ok], ka[ok] * a, lw=3.6, color="#c9c9c9",
            label="analytisch: $\\kappa a=\\Pi$")
    ax.plot(pi[ok], kn[ok] * a, lw=1.4, color="#0f6b3d",
            label="numerisch, aus den Knoten gemessen")
    ax.set_xlabel("$\\Pi = \\Delta p\\,a/\\gamma$")
    ax.set_ylabel("Krümmung $\\kappa\\,a$,  "
                  "$\\kappa=\\kappa_\\mathrm{m}+\\kappa_\\varphi$")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.8, loc="upper left")
    ax.set_title("Mittlere Krümmung: über die Fläche konstant,\nim Druck linear",
                 fontsize=9.6)

    ax = axes[2]
    mark_range(ax)
    ax.plot(pi[ok], psi[ok], lw=1.8, color="#7a1fa2")
    for y in (90, -90):
        ax.axhline(y, color="#8a1b1b", lw=1.0, ls=":")
    ax.text(0.0, 93, "senkrechte Tangente an der Kante", fontsize=7.6, ha="center",
            color="#8a1b1b")
    ax.set_xlabel("$\\Pi = \\Delta p\\,a/\\gamma$")
    ax.set_ylabel("$\\psi$ an der Kontaktlinie [Grad]")
    ax.set_ylim(-118, 118)
    ax.grid(alpha=0.25)
    ax.set_title("Grenze des darstellbaren Bereichs", fontsize=9.6)
    ax.text(0.03, 0.03,
            f"rot: $|\\Pi|>2$ – dort existiert keine glatte,\n"
            f"an $a$ gepinnte Fläche. Der Löser meldet\n"
            f"PressureOutsideCapillaryRange\n"
            f"({int(m['sweep_points_refused'])} von {len(rs)} Sweep-Punkten)\n"
            f"und gibt KEINE Ersatzform zurück.",
            transform=ax.transAxes, fontsize=7.4, color="#8a1b1b", va="bottom",
            bbox=dict(fc="#ffffff", ec="#e7b7b7", lw=0.8, alpha=0.95))

    heading(fig, m, "Apexhöhe, Krümmung und die Grenze des darstellbaren Drucks")
    limit_rows = ";  ".join(f"Pi={float(r['Pi']):.12g} → {r['status']}" for r in rng
                            if abs(float(r["Pi"])) >= 1.999)
    caveat(fig,
           f"Grenze: |delta_p| <= 2 gamma/a = {float(m['delta_p_max_Pa']):.0f} Pa. Bei |Pi|=2 "
           "steht die Meridiankurve an der Kante senkrecht (Halbkugel); darüber wäre der "
           "Kugelradius 2 gamma/delta_p kleiner als der Pinningradius. Der Löser stellt das an "
           "der integrierten Form fest – die Kurve wird senkrecht, bevor sie a erreicht – und "
           "gibt einen eigenen Status zurück:  " + limit_rows + ".  " + NOT_MODELLED,
           y=0.078)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 5
def figure_liquid_example(d, m, out):
    pr = profiles(d, "liquid_example_profiles.csv", key="delta_p_Pa")
    ex = rows(os.path.join(d, "liquid_example.csv"))
    liq = keyvalues(os.path.join(d, "liquid.csv"))
    a = float(m["contact_radius_m"])
    s = 1e6
    au = a * s
    dp = floats(ex, "delta_p_Pa")
    h = floats(ex, "h_m")
    dpmax = float(m["example_liquid_delta_p_max_Pa"])

    fig = plt.figure(figsize=(15.4, 7.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], left=0.05, right=0.985,
                          top=0.830, bottom=0.265, wspace=0.20)

    ax = fig.add_subplot(gs[0, 0])
    for tag, e in sorted(pr.items(), key=lambda kv: kv[1]["Pi"]):
        ax.plot(e["r"] * s, e["z"] * s, lw=1.9, color=pi_colour(e["Pi"]),
                label=f"$\\Delta p$ = {float(tag):+7.0f} Pa   ($\\Pi$={e['Pi']:+.2f})")
    ax.plot([au], [0.0], "o", ms=8, mfc="#d62728", mec="#5c0d0d", zorder=8)
    ax.axvline(0, color="#7a1fa2", lw=1.2, ls=(0, (6, 3)))
    ax.set_xlim(-0.3, 2.15 * au)
    ax.set_xlabel("r [µm]")
    ax.set_ylabel("z [µm]")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.4, loc="center right", framealpha=0.95)
    ax.set_title("Meniskusformen für den Stoffbeispiel-Datensatz", fontsize=9.8)

    ax = fig.add_subplot(gs[0, 1])
    ax.axvspan(dpmax, 1.3 * dpmax, color="#f6d5d5")
    ax.axvspan(-1.3 * dpmax, -dpmax, color="#f6d5d5")
    for x in (-dpmax, dpmax):
        ax.axvline(x, color="#8a1b1b", lw=1.2, ls="--")
    ax.plot(dp, h * s, "o-", lw=1.6, ms=6, color="#1f3d6b")
    ax.plot([-dpmax, dpmax], [-au, au], "o", ms=7, mfc="#8a1b1b", mec="#4a0a0a",
            label="Halbkugel, $h=\\pm a$")
    ax.set_xlim(-1.3 * dpmax, 1.3 * dpmax)
    ax.set_xlabel("$\\Delta p_\\mathrm{exit} = p_\\mathrm{fl}-p_\\mathrm{vak}$ [Pa]")
    ax.set_ylabel("Apexhöhe h [µm]")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.8, loc="lower right")
    ax.set_title(f"Apexhöhe über dem Druck – Grenze $2\\gamma/a$ = {dpmax:.0f} Pa",
                 fontsize=9.8)
    ax.text(0.03, 0.96,
            f"$\\gamma$ = {float(m['example_liquid_gamma'])*1e3:.2f} mN/m  (illustrative)\n"
            f"$\\rho$ = {float(m['example_liquid_rho']):.0f} kg/m³  (illustrative)\n"
            f"T = {float(m['temperature_K']):.2f} K\n"
            f"Bo = {float(m['example_liquid_bond_number']):.2e}",
            transform=ax.transAxes, fontsize=8.2, va="top",
            bbox=dict(fc="#ffffff", ec="#999999", lw=0.8, alpha=0.95))

    heading(fig, m, "STOFFBEISPIEL EMI-BF4 – ausdrücklich als illustrativ gekennzeichnet",
            banner_colour="#8a1b1b", extra_banner=" – KEINE Absolutaussage")
    caveat(fig,
           "BELEGT ist allein die Stoffidentität: die Kunze-Dissertation nennt EMI-BF4 als "
           "Treibstoff der geraden 10-µm-Kapillaren (Abschnitt 2.3.2, gedruckte Seite 28; "
           "Tabelle „List of Publications“, gedruckte Seite 30, Publikationen I–IV). NICHT "
           "belegt sind die Zahlenwerte: die Arbeit enthält keine Stoffwerttabelle und keinen "
           "numerischen Wert für Oberflächenspannung oder Dichte. gamma und rho sind "
           "unverändert aus der quellenlosen Tabelle in src/fluid.cpp übernommen und tragen "
           f"deshalb den Status {liq.get('status','?')}. Druck- und Höhenwerte dieser "
           "Abbildung sind ein Rechenbeispiel, kein Betriebspunkt; die Prüfung des Lösers "
           "läuft unabhängig davon dimensionslos. Viskosität, Leitfähigkeit und Permittivität "
           "sind nur vorgemerkt und gehen hier nicht ein. " + NOT_MODELLED,
           y=0.070, color="#8a1b1b")
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------- main
STAMP = ""


def main(d):
    global STAMP
    m = meta(d)
    root = pv.repo_root_of(__file__)
    rel = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/")
    state = pv.repo_state(root, ignore=[rel])
    STAMP = pv.stamp(state, m.get("commit"))
    pv.write_provenance(d, state, m.get("commit"))

    made = [figure_geometry(d, m, os.path.join(d, "fig1_geometry_and_bc.png")),
            figure_profiles(d, m, os.path.join(d, "fig2_profiles.png")),
            figure_convergence(d, m, os.path.join(d, "fig3_convergence.png")),
            figure_apex_and_curvature(d, m, os.path.join(d, "fig4_apex_and_curvature.png")),
            figure_liquid_example(d, m, os.path.join(d, "fig5_liquid_example.png"))]
    for f in made:
        print(f"{f}  {os.path.getsize(f)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
