#!/usr/bin/env python3
"""Figures of the P3b self-consistent electro-capillary equilibrium.

    ./build/es_electrocapillary examples/device_p1.cfg examples/electrocapillary_p3b.cfg \
        results/<dir> meta.commit=$(git rev-parse HEAD)
    python python/plot_electrocapillary.py results/<dir>

Reads only the CSVs written by apps/es_electrocapillary.cpp.  Nothing here
solves, meshes or fits; the singularity exponent drawn in figure 3 is the one
the C++ run fitted and wrote out.

The commit named on every figure is measured from git at drawing time by
python/provenance.py, not taken from an argument, and a dirty working tree is
stamped DIRTY.

Produces
  <dir>/fig1_moving_mesh.png      conformal volume mesh for a flat, an outward
                                  and an inward meniscus, with apex and edge detail
  <dir>/fig2_field.png            potential and field for the same shapes, on
                                  identical colour and length scales
  <dir>/fig3_edge_field.png       E_n(d) and p_M(d) at the exit edge over four
                                  mesh levels, with the fitted singularity
  <dir>/fig4_gate.png             integrated Maxwell force and weak projection
                                  against mesh level and exclusion distance
  <dir>/fig5_shapes_vs_voltage.png self-consistent meniscus shapes over voltage
  <dir>/fig6_branch.png           apex height, edge-far field, curvature and
                                  mechanical residual over voltage
  <dir>/fig7_convergence.png      mesh and coupling convergence at three points
  <dir>/figures_provenance.txt    what the stamp said, and whether it is releasable

Figures 5 to 7 are produced only where the edge gate passed.  Where it did not,
they are deliberately absent and the reason is EdgeLoadNotWellPosed.
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
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance as pv

TITLE = ("P3b: selbstkonsistentes statisches Elektro-Kapillargleichgewicht – "
         "Perfect Conductor, KEINE Emission")

NOT_MODELLED = ("Nicht enthalten: Emission, endliche Leitfähigkeit, Strömung, Viskosität, "
                "Raumladung, Zeitabhängigkeit, dynamische Stabilität, Taylor-Kegel, Cone-Jet, "
                "Schwerkraft. Δp_exit bleibt eine Eingabe.")

EDGE_NOTE = ("Die Kantenzone ist rot hinterlegt: dort ist der PUNKTWEISE Wert nicht "
             "netzkonvergent und wird nirgends verwendet.")


# --------------------------------------------------------------------- input
def rows(path):
    if not os.path.exists(path):
        return []
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


def parameters(d):
    out = {}
    for r in rows(os.path.join(d, "parameters.csv")):
        try:
            out[r["name"]] = float(r["value_SI"])
        except (ValueError, KeyError):
            pass
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
def shape_colour(pi):
    t = 0.42 + 0.53 * min(1.0, abs(pi) / 2.0)
    if pi > 0:
        return plt.get_cmap("Reds")(t)
    if pi < 0:
        return plt.get_cmap("Blues")(t)
    return "#111111"


def level_colour(lvl, lo=1, hi=4):
    t = 0.25 + 0.7 * (lvl - lo) / max(1, hi - lo)
    return plt.get_cmap("viridis")(t)


def provenance(fig, m, stamp, extra="", y=0.010, width=178):
    dirty = pv.DIRTY_MARK in stamp or stamp == pv.NO_VERSION
    head = (f"Commit {stamp}   |   {m.get('app','es_electrocapillary')}   |   "
            f"Konfiguration {m.get('config','?')}   |   "
            f"a = {float(m['contact_radius_m'])*1e6:.2f} µm, "
            f"γ = {float(m['surface_tension_N_per_m'])*1e3:.2f} mN/m "
            f"(Stoffstatus {m.get('liquid_status','?')}), "
            f"γ/a = {float(m['capillary_pressure_scale_Pa']):.0f} Pa, "
            f"ε_r(Emitter) = {float(m.get('emitter_eps_r',0)):.2f} "
            f"({m.get('emitter_eps_r_status','?')}), "
            f"Referenznetzstufe {m.get('reference_level','?')}")
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


def heading(fig, m, subtitle, banner=""):
    fig.suptitle(subtitle + "\n" + TITLE, fontsize=11.6, y=0.983, va="top")
    fig.text(0.5, 0.918,
             "statisch  ·  ideal leitende Flüssigkeit  ·  keine Emission  ·  "
             "keine Stabilitätsaussage" + banner,
             ha="center", va="top", fontsize=8.6, color="#1f3d6b")


# ------------------------------------------------------------------- figure 1
def figure_mesh(d, m, par, out):
    nodes = defaultdict(dict)
    for r in rows(os.path.join(d, "mesh_nodes.csv")):
        nodes[r["shape"]][(int(r["i"]), int(r["j"]))] = (float(r["r_m"]), float(r["z_m"]))
    surf = defaultdict(list)
    for r in rows(os.path.join(d, "mesh_surface.csv")):
        surf[r["shape"]].append((float(r["r_m"]), float(r["z_m"])))
    a = float(m["contact_radius_m"])
    s = 1e6
    au = a * s
    tags = [t for t in ("pi+0.00", "pi+1.50", "pi-1.50") if t in nodes]
    titles = {"pi+0.00": "eben  (Π = 0)", "pi+1.50": "nach außen  (Π = +1.5)",
              "pi-1.50": "nach innen  (Π = −1.5)"}

    def draw(ax, tag, rlim, zlim, lw=0.45):
        nd = nodes[tag]
        ii = sorted({k[0] for k in nd})
        jj = sorted({k[1] for k in nd})
        for j in jj:
            pts = [nd[(i, j)] for i in ii if (i, j) in nd]
            if len(pts) > 1:
                p = np.array(pts) * s
                ax.plot(p[:, 0], p[:, 1], color="#8a8a8a", lw=lw, zorder=1)
        for i in ii:
            pts = [nd[(i, j)] for j in jj if (i, j) in nd]
            if len(pts) > 1:
                p = np.array(pts) * s
                ax.plot(p[:, 0], p[:, 1], color="#8a8a8a", lw=lw, zorder=1)
        sp = np.array(surf[tag]) * s
        ax.plot(sp[:, 0], sp[:, 1], color="#1f3d6b", lw=2.0, zorder=4,
                label="freie Oberfläche")
        ax.plot([au], [0.0], "o", ms=7, mfc="#d62728", mec="#5c0d0d", zorder=6)
        ax.axvline(0, color="#7a1fa2", lw=1.2, ls=(0, (6, 3)), zorder=3)
        ax.set_xlim(*rlim)
        ax.set_ylim(*zlim)
        ax.set_aspect("equal")
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")

    fig = plt.figure(figsize=(15.2, 8.6))
    gs = fig.add_gridspec(2, 3, left=0.055, right=0.985, top=0.845, bottom=0.215,
                          wspace=0.24, hspace=0.36)
    for k, tag in enumerate(tags):
        ax = fig.add_subplot(gs[0, k])
        draw(ax, tag, (-0.1 * au, 1.7 * au), (-1.7 * au, 1.7 * au))
        ax.set_title(titles.get(tag, tag), fontsize=9.6)
        if k == 0:
            ax.legend(fontsize=7.6, loc="upper right", framealpha=0.95)

    # The apex window follows the apex the mesh actually has, not a guess.
    apex = {t: (np.array(surf[t])[0, 1] * s) for t in surf}
    zooms = []
    if "pi+1.50" in nodes:
        h = apex["pi+1.50"]
        zooms.append(("pi+1.50", (-0.02 * au, 0.45 * au), (h - 0.30 * au, h + 0.12 * au),
                      "Apex, Π = +1.5"))
        zooms.append(("pi+1.50", (0.78 * au, 1.30 * au), (-0.30 * au, 0.40 * au),
                      "Austrittskante, Π = +1.5"))
    if "pi-1.50" in nodes:
        zooms.append(("pi-1.50", (0.78 * au, 1.30 * au), (-0.60 * au, 0.25 * au),
                      "Austrittskante, Π = −1.5"))
    for k, (tag, rl, zl, ttl) in enumerate(zooms):
        ax = fig.add_subplot(gs[1, k])
        # The zone in which no POINTWISE value converges: within 0.1 a of the
        # contact line, and only on the liquid side -- the solid emitter beyond
        # r = a is not part of the free surface at all.
        if k >= 1:
            ax.axvspan(0.9 * au, au, color="#f6d5d5", zorder=0)
        draw(ax, tag, rl, zl, lw=0.7)
        ax.set_title(ttl, fontsize=9.6)

    q = rows(os.path.join(d, "mesh_quality.csv"))
    worst_j = min(floats(q, "min_jacobian")) if q else float("nan")
    worst_v = max(floats(q, "V_error")) if q else float("nan")
    inv = int(max(floats(q, "inverted_cells"))) if q else -1

    heading(fig, m, "Konformes bewegliches Volumennetz")
    caveat(fig,
           f"Die Zeile des P2c-Netzes, die die freie Oberfläche trägt, wird auf die "
           f"vorgeschriebene Meridiankurve gezogen; außerhalb der Bohrung bewegt sich kein "
           f"Knoten, und für die ebene Fläche ist das Netz bitgleich das P2c-Netz. Über alle "
           f"geprüften Formen und Netzstufen: kleinste Jacobi-Determinante {worst_j:.2e}, "
           f"{inv} invertierte Zellen, größter Fehler des Flüssigkeitsvolumens gegen die "
           f"geschlossene Form {worst_v:.1e}. " + EDGE_NOTE + " " + NOT_MODELLED, y=0.070)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 2
def figure_field(d, m, par, out):
    data = defaultdict(lambda: ([], [], [], []))
    for r in rows(os.path.join(d, "field_window.csv")):
        rr, zz, pp, ee = data[r["shape"]]
        rr.append(float(r["r_m"]))
        zz.append(float(r["z_m"]))
        pp.append(float(r["phi_V"]) if r["phi_V"] else np.nan)
        ee.append(float(r["E_V_per_m"]) if r["E_V_per_m"] else np.nan)
    surf = defaultdict(list)
    for r in rows(os.path.join(d, "mesh_surface.csv")):
        surf[r["shape"]].append((float(r["r_m"]), float(r["z_m"])))
    a = float(m["contact_radius_m"])
    s = 1e6
    tags = [t for t in ("pi+0.00", "pi+1.50", "pi-1.50") if t in data]
    titles = {"pi+0.00": "eben  (Π = 0)", "pi+1.50": "nach außen  (Π = +1.5)",
              "pi-1.50": "nach innen  (Π = −1.5)"}

    grids = {}
    for t in tags:
        rr, zz, pp, ee = (np.array(x) for x in data[t])
        ru, zu = np.unique(rr), np.unique(zz)
        shape = (len(zu), len(ru))
        grids[t] = (ru * s, zu * s, pp.reshape(shape), ee.reshape(shape))

    V = float(m["V_probe_V"])
    emax = np.nanmax([np.nanmax(g[3]) for g in grids.values()])
    emin = max(1e4, np.nanmin([np.nanmin(g[3]) for g in grids.values()]))

    fig = plt.figure(figsize=(15.2, 8.8))
    gs = fig.add_gridspec(2, 3, left=0.055, right=0.93, top=0.845, bottom=0.155,
                          wspace=0.22, hspace=0.28)
    for k, t in enumerate(tags):
        ru, zu, ph, em = grids[t]
        ax = fig.add_subplot(gs[0, k])
        cf = ax.contourf(ru, zu, ph, levels=np.linspace(0, V, 21), cmap="viridis")
        ax.contour(ru, zu, ph, levels=np.linspace(0, V, 21), colors="w", linewidths=0.35)
        sp = np.array(surf[t]) * s
        ax.plot(sp[:, 0], sp[:, 1], color="#ffffff", lw=2.0)
        ax.set_title(titles.get(t, t), fontsize=9.6)
        ax.set_aspect("equal")
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        if k == len(tags) - 1:
            cax = fig.add_axes([0.945, 0.50, 0.012, 0.30])
            fig.colorbar(cf, cax=cax, label="φ [V]")

        ax = fig.add_subplot(gs[1, k])
        im = ax.pcolormesh(ru, zu, em, norm=LogNorm(vmin=emin, vmax=emax), cmap="magma",
                           shading="auto")
        ax.plot(sp[:, 0], sp[:, 1], color="#00e5ff", lw=1.8)
        ax.set_aspect("equal")
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        ax.set_title("|E|", fontsize=9.0)
        if k == len(tags) - 1:
            cax = fig.add_axes([0.945, 0.17, 0.012, 0.30])
            fig.colorbar(im, cax=cax, label="|E| [V/m]")

    heading(fig, m, f"Potential und Feldstärke auf vorgeschriebenen Formen, {V:.0f} V")
    caveat(fig,
           "Dieselben Farb- und Längenskalen in jeder Spalte, damit die Formen vergleichbar "
           "sind. Die Flüssigkeit ist ein idealer Leiter auf V_emitter und trägt kein Feld; "
           "diese Zellen sind leer. Die Formen sind hier VORGESCHRIEBEN, nicht "
           "selbstkonsistent – das ist der Zustand, in dem das Kanten-Gate gerechnet wird. "
           + NOT_MODELLED, y=0.050)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 3
def figure_edge(d, m, par, out):
    prof = rows(os.path.join(d, "edge_profiles.csv"))
    segs = rows(os.path.join(d, "edge_segments.csv"))
    gate = rows(os.path.join(d, "gate_verdict.csv"))
    lev = rows(os.path.join(d, "gate_levels.csv"))
    a = float(m["contact_radius_m"])
    lo, hi = int(m["gate_min_level"]), int(m["gate_max_level"])

    by = defaultdict(list)
    for r in prof:
        by[(r["shape"], int(r["level"]))].append(r)
    bys = defaultdict(list)
    for r in segs:
        bys[(r["shape"], int(r["level"]))].append(r)

    fig = plt.figure(figsize=(15.4, 6.8))
    gs = fig.add_gridspec(1, 3, left=0.055, right=0.985, top=0.835, bottom=0.245, wspace=0.26)

    ax = fig.add_subplot(gs[0, 0])
    for lvl in range(lo, hi + 1):
        rs = by.get(("pi+0.00", lvl))
        if not rs:
            continue
        dd = floats(rs, "d_m") / a
        en = np.abs(floats(rs, "En_V_per_m"))
        k = dd > 0
        ax.loglog(dd[k], en[k], "o-", ms=3.5, lw=1.2, color=level_colour(lvl, lo, hi),
                  label=f"Stufe {lvl}")
    ax.set_xlabel("Abstand von der Kontaktlinie  d/a")
    ax.set_ylabel("$|E_n|$ [V/m]  (punktweise)")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.6)
    ax.set_title("Einseitiges Normalfeld an der Kante, ebene Form", fontsize=9.6)
    ax.text(0.03, 0.05, "Die Kurven laufen mit jeder Verfeinerung\nweiter nach oben: der "
                        "punktweise Kantenwert\nkonvergiert NICHT.", transform=ax.transAxes,
            fontsize=7.4, color="#8a1b1b", va="bottom")

    ax = fig.add_subplot(gs[0, 1])
    beta = None
    for r in gate:
        if r["shape"] == "pi+0.00":
            beta = float(r["fit_exponent"])
    for lvl in range(lo, hi + 1):
        rs = bys.get(("pi+0.00", lvl))
        if not rs:
            continue
        dd = floats(rs, "d_mid_m") / a
        pp = floats(rs, "pressure_Pa")
        k = (dd > 0) & (pp > 0)
        ax.loglog(dd[k], pp[k], "o-", ms=3.5, lw=1.2, color=level_colour(lvl, lo, hi),
                  label=f"Stufe {lvl}")
    if beta is not None:
        fit_lo, fit_hi = None, None
        for r in lev:
            if r["shape"] == "pi+0.00" and int(r["level"]) == hi:
                fit_lo, fit_hi = float(r["fit_d_lo_m"]) / a, float(r["fit_d_hi_m"]) / a
        if fit_lo:
            xx = np.logspace(np.log10(fit_lo), np.log10(fit_hi), 20)
            rs = bys.get(("pi+0.00", hi))
            dd = floats(rs, "d_mid_m") / a
            pp = floats(rs, "pressure_Pa")
            k = np.argmin(np.abs(dd - np.sqrt(fit_lo * fit_hi)))
            amp = pp[k] / (dd[k] ** beta)
            ax.loglog(xx, amp * xx ** beta, "k--", lw=1.6,
                      label=f"Fit  $p_M \\propto d^{{{beta:.3f}}}$")
            ax.axvspan(fit_lo, fit_hi, color="#eef3c8", zorder=0)
    ax.set_xlabel("Abstand von der Kontaktlinie  d/a")
    ax.set_ylabel("Segmentdruck $\\bar p_M$ [Pa]")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.4, loc="lower left")
    ax.set_title("Flächengemittelte Last und Singularitätsfit", fontsize=9.6)

    ax = fig.add_subplot(gs[0, 2])
    verdicts = {}
    for r in gate:
        verdicts[r["shape"]] = (float(r["Pi"]), r["verdict"], float(r["fit_exponent"]))
    for tag, (piv, verdict, b) in sorted(verdicts.items(), key=lambda kv: kv[1][0]):
        rs = bys.get((tag, hi))
        if not rs:
            continue
        dd = floats(rs, "d_mid_m") / a
        pp = floats(rs, "pressure_Pa")
        k = (dd > 0) & (pp > 0)
        ok = verdict == "Passed"
        ax.loglog(dd[k], pp[k], "-" if ok else "--", lw=1.8 if ok else 1.4,
                  color=shape_colour(piv),
                  label=f"Π={piv:+.1f}, β={b:+.3f} {'' if ok else '✗'}")
    ax.axhline(float(m["capillary_pressure_scale_Pa"]), color="#0f6b3d", lw=1.0, ls=":")
    ax.text(0.03, 0.94, "γ/a", transform=ax.transAxes, fontsize=7.6, color="#0f6b3d")
    ax.set_xlabel("Abstand von der Kontaktlinie  d/a")
    ax.set_ylabel("Segmentdruck $\\bar p_M$ [Pa]")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.0, loc="lower left")
    ax.set_title(f"Alle geprüften Formen, feinste Stufe {hi}\n"
                 "gestrichelt: Gate nicht bestanden", fontsize=9.6)

    heading(fig, m, "Feld und Maxwell-Druck an der unverrundeten Austrittskante")
    caveat(fig,
           "Die Kante ist eine Leiter/Dielektrikum/Vakuum-Kante ohne Verrundung; dort ist das "
           "Feld singulär und kein punktweiser Wert netzkonvergent. Gemessen wird der Exponent "
           f"β in p_M ~ d^β; integrierbar gegen 2πr ds ist er für β > −1. Für die ebene und "
           "die nach außen gewölbte Form ist er das, für die nach innen gezogene wird die "
           "Leiterkante einspringend und β fällt unter −1: dort wird NICHT gekoppelt. Es wird "
           "nichts abgeschnitten, kein Maximalwert codiert, keine Ausschlusszone frei gewählt "
           "und kein Kantenradius erfunden. " + NOT_MODELLED, y=0.045)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 4
def figure_gate(d, m, par, out):
    lev = rows(os.path.join(d, "gate_levels.csv"))
    gate = rows(os.path.join(d, "gate_verdict.csv"))
    a = float(m["contact_radius_m"])
    by = defaultdict(list)
    for r in lev:
        by[r["shape"]].append(r)

    fig = plt.figure(figsize=(15.4, 7.0))
    gs = fig.add_gridspec(1, 3, left=0.055, right=0.985, top=0.835, bottom=0.265, wspace=0.26)

    ax = fig.add_subplot(gs[0, 0])
    for tag, rs in sorted(by.items(), key=lambda kv: float(kv[1][0]["Pi"])):
        piv = float(rs[0]["Pi"])
        n = floats(rs, "n_nodes")
        ax.semilogx(n, floats(rs, "total_force_N") * 1e6, "o-", ms=5, lw=1.5,
                    color=shape_colour(piv), label=f"Π={piv:+.1f}")
    ax.set_xlabel("Knoten des Volumennetzes")
    ax.set_ylabel("integrierte Maxwell-Kraft [µN]")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.4)
    ax.set_title("Gesamtkraft gegen die Netzstufe", fontsize=9.6)

    ax = fig.add_subplot(gs[0, 1])
    styles = [("force_beyond_0.10a_N", "-", "d ≥ 0.10 a"),
              ("force_beyond_0.05a_N", "--", "d ≥ 0.05 a"),
              ("force_beyond_0.025a_N", ":", "d ≥ 0.025 a")]
    for tag, rs in sorted(by.items(), key=lambda kv: float(kv[1][0]["Pi"])):
        piv = float(rs[0]["Pi"])
        if abs(piv) > 0.01:
            continue
        n = floats(rs, "n_nodes")
        for col, ls, lab in styles:
            ax.semilogx(n, floats(rs, col) * 1e6, ls, marker="o", ms=4, lw=1.5,
                        color=shape_colour(piv), label=lab)
        ax.semilogx(n, floats(rs, "total_force_N") * 1e6, "-", lw=2.4, color="#111111",
                    label="ohne Ausschluss")
    for r in gate:
        if r["shape"] != "pi+0.00":
            continue
        ax.axhline(float(r["limit_force_mesh_N"]) * 1e6, color="#0f6b3d", lw=1.2, ls="-.")
        ax.axhline(float(r["limit_force_exclusion_N"]) * 1e6, color="#7a1fa2", lw=1.2, ls="-.")
        ax.text(0.03, 0.90,
                f"Grenzkraft über die Netzstufen: {float(r['limit_force_mesh_N'])*1e6:.4f} µN\n"
                f"Grenzkraft über die Ausschlussdistanz: "
                f"{float(r['limit_force_exclusion_N'])*1e6:.4f} µN\n"
                f"Abweichung {float(r['limit_agreement']):.2e} "
                f"(Grenze {float(r['tol_limit_agreement']):.0e})",
                transform=ax.transAxes, fontsize=7.4, va="top", color="#333333",
                bbox=dict(fc="#ffffff", ec="#999999", lw=0.8, alpha=0.95))
    ax.set_xlabel("Knoten des Volumennetzes")
    ax.set_ylabel("Kraft jenseits der Ausschlussdistanz [µN]")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.2, loc="lower right")
    ax.set_title("Ebene Form: bei FESTER Ausschlussdistanz ist die Kraft\n"
                 "netzkonvergent; ohne Ausschluss läuft sie langsam hoch", fontsize=9.6)

    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    lines = ["Urteil des Kanten-Gates", ""]
    for r in gate:
        ok = r["verdict"] == "Passed"
        lines.append(f"Π = {float(r['Pi']):+.2f}   {r['verdict']}{'' if ok else '  ✗'}")
        lines.append(f"    β = {float(r['fit_exponent']):+.3f}  "
                     f"(integrierbar: {'ja' if float(r['fit_exponent']) > -1 else 'NEIN'})")
        lines.append(f"    kantenferne Last  {float(r['edge_far_change']):.2e} "
                     f"(≤ {float(r['tol_edge_far']):.0e})")
        lines.append(f"    Gesamtkraft       {float(r['total_force_change']):.2e} "
                     f"(≤ {float(r['tol_total_force']):.0e})")
        lines.append(f"    Grenzkraftabgleich {float(r['limit_agreement']):.2e} "
                     f"(≤ {float(r['tol_limit_agreement']):.0e})")
        lines.append("")
    lines.append(f"Gekoppelt wird nur für Kontaktwinkel")
    lines.append(f"ψ ∈ [{float(m.get('gate_psi_admissible_lo_deg',0)):.1f}°, "
                 f"{float(m.get('gate_psi_admissible_hi_deg',0)):.1f}°].")
    ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, fontsize=8.0, va="top",
            family="monospace")

    heading(fig, m, "Kanten-Gate: integrierte Kraft und schwache Lastprojektion")
    caveat(fig,
           "Die schwache Projektion ist segmentweise konservativ: die Kraft jedes Segments wird "
           "mit dem einseitigen Feld integriert, der Segmentdruck ist Kraft durch "
           "Rotationsfläche, und die Summe über die Segmente IST die integrierte Maxwell-Kraft. "
           "Bei fester Ausschlussdistanz ist diese Kraft netzkonvergent; ohne Ausschluss "
           "nähert sie sich ihrem Grenzwert wie d^(1+β) und wird deshalb über die Netzstufen "
           "extrapoliert. Dass diese Extrapolation mit der über die Ausschlussdistanz "
           "übereinstimmt, ist der Nachweis, dass der Kantenbeitrag summierbar ist. "
           + NOT_MODELLED, y=0.055)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 5
def figure_shapes(d, m, par, out):
    prof = rows(os.path.join(d, "coupled_profiles.csv"))
    if not prof:
        return None
    a = float(m["contact_radius_m"])
    s = 1e6
    au = a * s
    by = defaultdict(lambda: defaultdict(list))
    for r in prof:
        by[r["branch"]][float(r["V_emitter_V"])].append(
            (float(r["r_m"]), float(r["z_m"]), float(r["load_Pa"])))
    branches = [b for b in ("dp0_pos", "dp_plus", "dp_minus") if b in by]
    labels = {"dp0_pos": "Δp_exit = 0", "dp_plus": "Δp_exit = +0.5 γ/a",
              "dp_minus": "Δp_exit = −0.5 γ/a"}

    fig = plt.figure(figsize=(15.2, 7.2))
    gs = fig.add_gridspec(1, len(branches) + 1, left=0.05, right=0.985, top=0.835,
                          bottom=0.215, wspace=0.26)
    for k, b in enumerate(branches):
        ax = fig.add_subplot(gs[0, k])
        vs = sorted(by[b])
        for V in vs:
            p = np.array([(x[0], x[1]) for x in by[b][V]]) * s
            t = (abs(V) / max(1e-9, max(abs(v) for v in vs)))
            ax.plot(p[:, 0], p[:, 1], lw=1.6, color=plt.get_cmap("plasma")(0.12 + 0.75 * t),
                    label=f"{V:+.0f} V" if V in (vs[0], vs[-1]) or len(vs) < 7 else None)
        ax.plot([au], [0.0], "o", ms=7, mfc="#d62728", mec="#5c0d0d", zorder=6)
        ax.axvline(0, color="#7a1fa2", lw=1.1, ls=(0, (6, 3)))
        ax.set_xlim(-0.1 * au, 1.35 * au)
        ax.set_aspect("equal")
        ax.set_xlabel("r [µm]")
        ax.set_ylabel("z [µm]")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.2, loc="lower left")
        ax.set_title(f"{labels.get(b, b)}\n{len(vs)} selbstkonsistente Punkte", fontsize=9.4)

    ax = fig.add_subplot(gs[0, len(branches)])
    for b in branches:
        vs = sorted(by[b])
        for V in vs[::max(1, len(vs) // 5)]:
            arr = np.array(by[b][V])
            ss = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(arr[:, 0]),
                                                           np.diff(arr[:, 1])))])
            t = abs(V) / max(1e-9, max(abs(v) for v in vs))
            ax.plot(ss / ss[-1], arr[:, 2] / float(m["capillary_pressure_scale_Pa"]), lw=1.4,
                    color=plt.get_cmap("plasma")(0.12 + 0.75 * t))
    ax.set_xlabel("normierte Bogenlänge  s/L   (0 = Apex, 1 = Kontaktlinie)")
    ax.set_ylabel("$p_M / (\\gamma/a)$, wie übergeben")
    ax.axvspan(1.0 - 0.1, 1.0, color="#f6d5d5", zorder=0)
    ax.grid(alpha=0.25)
    ax.set_title("Die übergebene Maxwell-Last\n(rot: nicht punktweise konvergente Zone)",
                 fontsize=9.4)

    heading(fig, m, "Selbstkonsistente Meniskusformen über der Spannung")
    caveat(fig,
           "Jede Kurve ist ein Fixpunkt aus Form, Netz, Feld und Last. Die Last zieht die "
           "Oberfläche IMMER nach außen, unabhängig von der Polarität, weil p_M = ε₀E_n²/2 "
           "quadratisch ist. Der Ast ist der, der die feldfreie P3a-Lösung enthält; wo er "
           "endet, ist die Stelle, an der DIESER Löser stehen bleibt – kein Emissionsbeginn, "
           "kein Taylor-Kegel, keine Stabilitätsaussage. " + NOT_MODELLED, y=0.058)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 6
def figure_branch(d, m, par, out):
    pts = rows(os.path.join(d, "coupled_points.csv"))
    if not pts:
        return None
    ends = rows(os.path.join(d, "branch_ends.csv"))
    a = float(m["contact_radius_m"])
    by = defaultdict(list)
    for r in pts:
        by[r["branch"]].append(r)
    labels = {"dp0_pos": "Δp = 0, V > 0", "dp0_neg": "Δp = 0, V < 0",
              "dp_plus": "Δp = +0.5 γ/a", "dp_minus": "Δp = −0.5 γ/a"}
    colours = {"dp0_pos": "#1f3d6b", "dp0_neg": "#7a1fa2", "dp_plus": "#b2182b",
               "dp_minus": "#2166ac"}

    fig = plt.figure(figsize=(14.8, 8.4))
    gs = fig.add_gridspec(2, 2, left=0.06, right=0.985, top=0.835, bottom=0.20,
                          wspace=0.22, hspace=0.30)
    panels = [("h_over_a", "Apexhöhe  h/a", False),
              ("E_edge_far_V_per_m", "$E_n$ bei d = 0.25 a  [V/m]", False),
              ("max_curvature_1_per_m", "maximale Krümmung $|\\kappa|$  [1/m]", False),
              ("mech_residual_edge_far", "mechanisches Residuum (kantenfern)", True)]
    for ax, (col, lab, logy) in zip([fig.add_subplot(gs[i, j]) for i in (0, 1) for j in (0, 1)],
                                    panels):
        for b, rs in by.items():
            V = np.abs(floats(rs, "V_emitter_V"))
            y = floats(rs, col)
            o = np.argsort(V)
            (ax.semilogy if logy else ax.plot)(V[o], np.abs(y[o]) if logy else y[o], "o-",
                                               ms=4.5, lw=1.5, color=colours.get(b, "#333"),
                                               label=labels.get(b, b))
        for r in ends:
            if r["branch"] in by and float(r["first_failed_V"]) != 0.0:
                ax.axvline(abs(float(r["first_failed_V"])), color=colours.get(r["branch"], "#333"),
                           lw=1.0, ls=":")
        if col == "mech_residual_edge_far":
            ax.axhline(1e-3, color="#8a1b1b", lw=1.2, ls="--")
            ax.text(0.03, 0.06, "vorab festgelegte Grenze 1e-3", transform=ax.transAxes,
                    fontsize=7.4, color="#8a1b1b")
        ax.set_xlabel("$|V_\\mathrm{emitter}|$ [V]   ($V_\\mathrm{extraktor}$ = 0)")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.4)

    txt = []
    for r in ends:
        txt.append(f"{r['branch']}: {int(float(r['points']))} Punkte bis "
                   f"{float(r['last_converged_V']):.0f} V, Ende {r['end_status']}")
    heading(fig, m, "Der Ast über der Spannung")
    caveat(fig,
           "Gepunktet: die erste Spannung, bei der ein Schritt scheiterte. " + "  ".join(txt) +
           ".  Das Ende eines Astes ist die Stelle, an der dieser Löser mit dieser "
           "Schrittweite und diesen vorab festgelegten Grenzen stehen bleibt. Es ist KEIN "
           "Emissionsbeginn, KEIN Taylor-Kegel-Onset und KEINE Aussage über dynamische "
           "Stabilität; nichts davon wird in P3b gerechnet. " + NOT_MODELLED, y=0.045)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- figure 7
def figure_convergence(d, m, par, out):
    rs = rows(os.path.join(d, "coupled_convergence.csv"))
    if not rs:
        return None
    a = float(m["contact_radius_m"])
    mesh = [r for r in rs if r["variant"] == "mesh"]
    other = [r for r in rs if r["variant"] != "mesh"]
    by = defaultdict(list)
    for r in mesh:
        by[r["point"]].append(r)
    colours = {"klein": "#2166ac", "mittel": "#0f6b3d", "gross": "#b2182b"}

    fig = plt.figure(figsize=(15.4, 7.0))
    gs = fig.add_gridspec(1, 3, left=0.055, right=0.985, top=0.835, bottom=0.245, wspace=0.26)
    cols = [("h_over_a", "Apexhöhe  h/a"),
            ("total_force_N", "integrierte Maxwell-Kraft [N]"),
            ("E_edge_far_V_per_m", "$E_n$ bei d = 0.25 a [V/m]")]
    for ax, (col, lab) in zip([fig.add_subplot(gs[0, k]) for k in range(3)], cols):
        for pt, g in by.items():
            lv = floats(g, "level")
            y = floats(g, col)
            o = np.argsort(lv)
            V = float(g[0]["V_emitter_V"])
            ax.plot(lv[o], y[o], "o-", ms=5.5, lw=1.6, color=colours.get(pt, "#333"),
                    label=f"{pt} ({V:.0f} V)")
        for r in other:
            ax.plot([float(r["level"])], [float(r[col])],
                    "x" if r["variant"] == "farfield" else "+", ms=9, mew=1.8,
                    color=colours.get(r["point"], "#333"))
        ax.set_xlabel("Netzstufe")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.4)
    fig.axes[0].set_title("x = geerdete Hülle statt asymptotisch,\n+ = andere Unterrelaxation",
                          fontsize=9.2)

    heading(fig, m, "Netz- und Kopplungskonvergenz an drei Betriebspunkten")
    caveat(fig,
           "Drei Verformungsgrade, drei Netzstufen. Zusätzlich zwei Kontrollen am "
           "Referenznetz: der Fernrand (geerdete Hülle statt asymptotischer Monopolbedingung) "
           "und eine andere Unterrelaxation der Fixpunktiteration – die Lösung darf von "
           "keiner der beiden abhängen. Die Intervallzahl des Kapillarlösers ist keine "
           "Eingabe; sie folgt aus der geforderten Genauigkeit. " + NOT_MODELLED, y=0.055)
    provenance(fig, m, STAMP)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------- main
STAMP = ""


def main(d):
    global STAMP
    m = meta(d)
    par = parameters(d)
    root = pv.repo_root_of(__file__)
    rel = os.path.relpath(os.path.abspath(d), root).replace(os.sep, "/")
    state = pv.repo_state(root, ignore=[rel])
    STAMP = pv.stamp(state, m.get("commit"))
    pv.write_provenance(d, state, m.get("commit"))

    made = [figure_mesh(d, m, par, os.path.join(d, "fig1_moving_mesh.png")),
            figure_field(d, m, par, os.path.join(d, "fig2_field.png")),
            figure_edge(d, m, par, os.path.join(d, "fig3_edge_field.png")),
            figure_gate(d, m, par, os.path.join(d, "fig4_gate.png"))]
    if m.get("coupled") == "yes":
        made += [figure_shapes(d, m, par, os.path.join(d, "fig5_shapes_vs_voltage.png")),
                 figure_branch(d, m, par, os.path.join(d, "fig6_branch.png")),
                 figure_convergence(d, m, par, os.path.join(d, "fig7_convergence.png"))]
    else:
        print("Das Kanten-Gate ist nicht bestanden: Abbildungen 5 bis 7 entfallen "
              "(EdgeLoadNotWellPosed).")
    for f in made:
        if f:
            print(f"{f}  {os.path.getsize(f)} Byte")
    print(f"Provenienz: {STAMP}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])
