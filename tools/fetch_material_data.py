#!/usr/bin/env python3
"""Fetch the material data of EMI-BF4 from NIST ILThermo and emit the versioned
C++ table `src/material_data_emibf4.cpp` plus the audit `docs/13_material_data.md`.

    python tools/fetch_material_data.py [--cache <dir>] [--offline]

WHY THIS IS A SCRIPT AND THE RESULT IS COMMITTED

The numbers must be reproducible and they must not depend on a network at build
time.  So the acquisition is a script that is versioned, and its OUTPUT -- a
plain C++ source file with no dependencies -- is versioned too.  Re-running the
script and finding the emitted file unchanged is the reproducibility check.

WHAT IS AND IS NOT DONE HERE

  * Every value is read from the JSON that ILThermo returns.  No number is
    typed in by hand, no number is copied from a search-result snippet, no
    number comes from anybody's memory.
  * The measurement method, the sample source, the purity and the water content
    are copied VERBATIM, including the fact that a source did not state them.
  * Nothing is averaged.  Sources that disagree are all emitted and the spread
    is computed by the C++ side at query time.
  * Values that are frequency resolved keep their frequency, so that a value
    measured at 18 GHz cannot back a DC number by accident.

THE SOURCE

NIST Standard Reference Database #147, ILThermo v2.0,
https://ilthermo.boulder.nist.gov/ .  It is a compilation: every dataset in it
carries the full citation of the PRIMARY publication the data were taken from,
together with the method and the sample description that publication reported.
The compilation is what is queried here; the primary reference is what is
recorded and what must be cited.

The manufacturer data sheet (IoLiTec IL-0006) is added by hand below because it
is not in ILThermo -- but only after the PDF itself was opened and read; the
values, the temperatures and the purity are transcribed from it, and the fields
it does not state are left empty rather than guessed.
"""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

COMPOUND = "1-ethyl-3-methylimidazolium tetrafluoroborate"
BASE = "https://ilthermo.boulder.nist.gov/ILT2"

# ILThermo column header -> (property key, SI unit)
UNITS = {
    "Surface tension liquid-gas, N/m": ("SurfaceTension", "N/m"),
    "Specific density, kg/m<SUP>3</SUP>": ("Density", "kg/m^3"),
    "Viscosity, Pa&#8226;s": ("DynamicViscosity", "Pa s"),
    "Kinematic viscosity, m<SUP>2</SUP>/s": ("KinematicViscosity", "m^2/s"),
    "Electrical conductivity, S/m": ("ElectricalConductivity", "S/m"),
    "Relative permittivity at zero frequency": ("RelativePermittivity", "-"),
}
ORDER = ["SurfaceTension", "Density", "DynamicViscosity", "KinematicViscosity",
         "ElectricalConductivity", "RelativePermittivity"]


def get(url, path, offline):
    if os.path.exists(path) and os.path.getsize(path) > 50:
        return True
    if offline:
        return False
    subprocess.run(["curl", "-s", "-m", "60", "-A", "electrospray-research/0.1",
                    url, "-o", path], check=False)
    time.sleep(0.12)
    return os.path.exists(path) and os.path.getsize(path) > 50


def cesc(s):
    """C string literal body."""
    if s is None:
        return ""
    out = []
    for ch in str(s):
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif ch in "\r\n":
            out.append(" ")
        elif ord(ch) < 32 or ord(ch) > 126:
            out.append("?")   # the table is ASCII; a name with an accent is not a value
        else:
            out.append(ch)
    return "".join(out)


# ---------------------------------------------------------------------------
# The manufacturer data sheet, transcribed from the PDF that was opened and read
# ---------------------------------------------------------------------------
#
# IoLiTec Ionic Liquids Technologies GmbH, "Technical Data Sheet,
# 1-Ethyl-3-methylimidazolium tetrafluoroborate", product code IL-0006,
# CAS 143314-16-3, Revision Date 2/13/2012, Date Issued 11/2/2012, page 2/4,
# section "3 PROPERTIES".  Purity from page 1/4, section 2: ">98%".
#
# The sheet states NO surface tension, NO relative permittivity, NO water
# content and NO measurement method.  Those fields are therefore empty here.
IOLITEC = {
    "reference": ("IoLiTec Ionic Liquids Technologies GmbH, Technical Data Sheet "
                  "1-Ethyl-3-methylimidazolium tetrafluoroborate, Produktcode IL-0006, "
                  "CAS 143314-16-3, Revision Date 2/13/2012, Date Issued 11/2/2012, "
                  "Seite 2/4, Abschnitt 3 PROPERTIES."),
    "title": "Technical Data Sheet IL-0006",
    "purity": ">98% (Datenblatt Seite 1/4, Abschnitt 2)",
    "method": "",           # the sheet does not state one
    "water_content": "",    # the sheet does not state one
    "sample_source": "Hersteller IoLiTec",
    "status": "ManufacturerSpec",
    "values": {
        # property        T [K]     value in SI
        "Density": (297.15, 1282.0),            # 1.282 g/cm^3 at 24 degC
        "DynamicViscosity": (298.15, 25.2e-3),  # 25.2 mPa*s at 25 degC
        "ElectricalConductivity": (298.15, 1.41),   # 14.1 mS/cm at 25 degC
    },
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=os.path.join(HERE, "_ilthermo_cache"))
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.cache, exist_ok=True)

    idx = os.path.join(args.cache, "index.json")
    url = (f"{BASE}/ilsearch?cmp={COMPOUND.replace(' ', '%20')}"
           f"&ncmp=1&year=&auth=&keyw=&prp=")
    if not get(url, idx, args.offline):
        print("Die ILThermo-Übersicht ist nicht erreichbar und nicht im Cache.",
              file=sys.stderr)
        return 2
    index = json.load(open(idx, encoding="utf-8"))
    print(f"ILThermo meldet {index['cnt']} Datensätze für {COMPOUND}")

    sources = {k: [] for k in ORDER}
    for row in index["res"]:
        sid = row[0]
        p = os.path.join(args.cache, sid + ".json")
        if not get(f"{BASE}/ilset?set={sid}", p, args.offline):
            continue
        d = json.load(open(p, encoding="utf-8"))
        heads = [h[0] for h in d.get("dhead", [])]
        it = next((k for k, h in enumerate(heads) if h.startswith("Temperature")), None)
        iv = next((k for k, h in enumerate(heads) if h in UNITS), None)
        if it is None or iv is None:
            continue
        ifr = next((k for k, h in enumerate(heads) if h.startswith("Frequency")), None)
        ipr = next((k for k, h in enumerate(heads) if h.startswith("Pressure")), None)
        key, unit = UNITS[heads[iv]]
        comp = (d.get("components") or [{}])[0]
        sample = {k.rstrip(":"): v for k, v in (comp.get("sample") or [])}
        ref = d.get("ref") or {}
        pts = []
        for pt in d["data"]:
            try:
                T = float(pt[it][0])
                v = float(pt[iv][0])
            except (ValueError, IndexError, TypeError):
                continue
            try:
                u = float(pt[iv][1])
            except (ValueError, IndexError, TypeError):
                u = 0.0
            fr = 0.0
            if ifr is not None:
                try:
                    fr = float(pt[ifr][0]) * 1.0e6   # MHz -> Hz
                except (ValueError, IndexError, TypeError):
                    fr = 0.0
            pr = 0.0
            if ipr is not None:
                try:
                    pr = float(pt[ipr][0]) * 1.0e3   # kPa -> Pa
                except (ValueError, IndexError, TypeError):
                    pr = 0.0
            pts.append((T, v, u, fr, pr))
        if not pts:
            continue
        sources[key].append({
            "database_id": sid,
            "reference": ref.get("full", ""),
            "paper_title": ref.get("title", ""),
            "method": d.get("expmeth", "") or "",
            "sample_source": sample.get("Source", ""),
            "purity": sample.get("Initial purity", ""),
            "water_content": sample.get("Final purity", ""),
            "constraints": "; ".join(d.get("constr") or []),
            "status": "Measured",
            "points": pts,
        })

    # the manufacturer sheet
    for key, (T, v) in IOLITEC["values"].items():
        sources[key].append({
            "database_id": "IoLiTec-IL-0006",
            "reference": IOLITEC["reference"],
            "paper_title": IOLITEC["title"],
            "method": IOLITEC["method"],
            "sample_source": IOLITEC["sample_source"],
            "purity": IOLITEC["purity"],
            "water_content": IOLITEC["water_content"],
            "constraints": "",
            "status": IOLITEC["status"],
            "points": [(T, v, 0.0, 0.0, 0.0)],
        })

    for k in ORDER:
        sources[k].sort(key=lambda s: (s["reference"], s["database_id"]))

    # --- the selection rule, applied mechanically ---------------------------
    P_AMB, P_TOL = 101325.0, 0.05

    def ambient(pt):
        return pt[4] <= 0.0 or abs(pt[4] - P_AMB) < P_TOL * P_AMB

    def select(lst):
        # THE SELECTION RULE, applied mechanically and nowhere tuned:
        #   1. method, purity and water content must all be stated;
        #   2. the source must have AMBIENT points and must not be frequency
        #      resolved -- a density at 60 MPa and a permittivity at 18 GHz are
        #      real measurements of something else;
        #   3. its ambient temperature range must contain 298.15 K, and among
        #      those the one with the most ambient points wins;
        #   4. ties are broken by the smaller reported uncertainty.
        cands = []
        for i, s in enumerate(lst):
            if not (s["method"] and s["purity"] and s["water_content"]):
                continue
            if any(p[3] > 0.0 for p in s["points"]):
                continue
            amb = [p for p in s["points"] if ambient(p)]
            if not amb:
                continue
            Ts = [p[0] for p in amb]
            if not (min(Ts) - 1e-9 <= 298.15 <= max(Ts) + 1e-9 or
                    (len(Ts) == 1 and abs(Ts[0] - 298.15) <= 1.0)):
                continue
            unc = min((p[2] for p in amb if p[2] > 0.0), default=float("inf"))
            cands.append((-len(amb), unc, i))
        if not cands:
            return -1
        cands.sort()
        return cands[0][2]

    chosen = {k: select(sources[k]) for k in ORDER}

    # --- emit the C++ table -------------------------------------------------
    out = []
    A = out.append
    A("// GENERATED by tools/fetch_material_data.py -- do not edit by hand.")
    A("//")
    A("// Material data for EMI-BF4.  Every value is read from the JSON that NIST")
    A("// ILThermo (Standard Reference Database #147) returns, plus one manufacturer")
    A("// data sheet that was opened and read.  Nothing here was typed from memory or")
    A("// copied from a search-result snippet, nothing is averaged, and a field a")
    A("// source did not state is EMPTY rather than guessed.")
    A("//")
    A("// Re-running the generator and finding this file unchanged is the")
    A("// reproducibility check.  The audit is docs/13_material_data.md.")
    A("")
    A('#include "es/material_data.hpp"')
    A("")
    A("namespace es {")
    A("namespace {")
    A("")
    for key in ORDER:
        for i, s in enumerate(sources[key]):
            A(f"const PropertyPoint kPts_{key}_{i}[] = {{")
            for (T, v, u, fr, pr) in s["points"]:
                A(f"    {{{T!r}, {v!r}, {u!r}, {fr!r}, {pr!r}}},")
            A("};")
    A("")
    for key in ORDER:
        A(f"const PropertySource kSrc_{key}[] = {{")
        for i, s in enumerate(sources[key]):
            A("    {PropertyKind::" + key + ",")
            A(f'     "{cesc(s["database_id"])}",')
            A(f'     "{cesc(s["reference"])}",')
            A(f'     "{cesc(s["paper_title"])}",')
            A(f'     "{cesc(s["method"])}",')
            A(f'     "{cesc(s["sample_source"])}",')
            A(f'     "{cesc(s["purity"])}",')
            A(f'     "{cesc(s["water_content"])}",')
            A(f'     "{cesc(s["constraints"])}",')
            A(f"     MaterialDataStatus::{s['status']},")
            A(f"     kPts_{key}_{i}, {len(s['points'])}}},")
        A("};")
        A("")
    A("const MaterialProperty kProps[] = {")
    for key in ORDER:
        A(f"    {{PropertyKind::{key}, kSrc_{key}, "
          f"sizeof(kSrc_{key}) / sizeof(kSrc_{key}[0]), {chosen[key]}}},")
    A("};")
    A("")
    A("const MaterialDataset kEmiBf4 = {")
    A('    "EMI-BF4 (1-Ethyl-3-methylimidazolium-tetrafluoroborat)",')
    A('    "143314-16-3",')
    A("    197.97e-3,")
    A('    "Stoffidentitaet des Geraets: KunzeFynn-2024-12-10.pdf, Abschnitt 2.3.2, '
      'gedruckte Seite 28 (PDF-Seite 36), und Tabelle \'List of Publications\', gedruckte '
      'Seite 30 (PDF-Seite 38), Publikationen I-IV.  Molmasse und CAS aus dem '
      'IoLiTec-Datenblatt IL-0006 sowie aus den ILThermo-Komponentendaten.",')
    A('    "Zahlenwerte aus NIST ILThermo (SRD #147) und dem IoLiTec-Datenblatt IL-0006.  '
      'Es wird NICHT gemittelt: widersprechende Quellen stehen nebeneinander und die '
      'Streuung ist eine abfragbare Groesse.",')
    A("    kProps, sizeof(kProps) / sizeof(kProps[0]),")
    A("};")
    A("")
    A("}  // namespace")
    A("")
    A("const MaterialDataset& emibf4_sourced() { return kEmiBf4; }")
    A("")
    A("}  // namespace es")
    A("")
    dst = os.path.join(ROOT, "src", "material_data_emibf4.cpp")
    open(dst, "w", encoding="utf-8", newline="\n").write("\n".join(out))
    print(f"geschrieben: {dst}")

    for key in ORDER:
        n = len(sources[key])
        sel = chosen[key]
        if sel >= 0:
            s = sources[key][sel]
            amb = [p for p in s["points"] if ambient(p)]
            near = min(amb, key=lambda p: abs(p[0] - 298.15))
            print(f"  {key:24s} {n:3d} Quellen -> gewaehlt {s['database_id']:>16s} "
                  f"{near[1]:.6g} bei {near[0]:.2f} K")
        else:
            print(f"  {key:24s} {n:3d} Quellen -> KEINE erfuellt die Auswahlregel "
                  f"(MissingMaterialData)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
