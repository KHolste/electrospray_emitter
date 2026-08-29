#!/usr/bin/env python3
"""The provenance stamp must detect a dirty working tree and say so.

Checked against real throw-away repositories, not against mocks: the whole point
of the module is that it talks to git, so a test that replaces git tests nothing.
"""
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import provenance as pv  # noqa: E402

failures = 0


def expect(what, ok):
    global failures
    if not ok:
        failures += 1
    print(f"  {what:<62s} {'OK' if ok else 'FAIL'}")


def git(root, *args):
    subprocess.run([shutil.which("git"), "-C", root] + list(args), check=True,
                   capture_output=True, text=True)


def make_repo(root):
    git(root, "init", "--quiet")
    git(root, "config", "user.email", "test@example.invalid")
    git(root, "config", "user.name", "test")
    git(root, "config", "commit.gpgsign", "false")
    with open(os.path.join(root, "a.txt"), "w", encoding="utf-8") as fh:
        fh.write("one\n")
    git(root, "add", "a.txt")
    git(root, "commit", "--quiet", "-m", "erste")


def main():
    if shutil.which("git") is None:
        print("git nicht gefunden -- Test kann nichts pruefen")
        return 1

    print("=== Provenienz: sauberer Baum ===")
    with tempfile.TemporaryDirectory() as root:
        make_repo(root)
        st = pv.repo_state(root)
        expect("Repository erkannt", st["available"])
        expect("Commit-ID ist ein voller SHA", st["commit"] is not None and len(st["commit"]) >= 40)
        expect("sauberer Baum -> nicht dirty", not st["dirty"])
        expect("Stempel ist die kurze Commit-ID", pv.stamp(st) == st["commit"][:12])
        expect("Stempel enthaelt kein DIRTY", pv.DIRTY_MARK not in pv.stamp(st))
        expect("freigebbar, wenn die Daten aus genau diesem Commit stammen",
               pv.clean_release(st, st["commit"]))
        expect("nicht freigebbar ohne Datenherkunft", not pv.clean_release(st, None))
        expect("nicht freigebbar bei fremder Datenherkunft",
               not pv.clean_release(st, "0" * 40))
        expect("fremde Datenherkunft wird im Stempel genannt",
               "NICHT HEAD" in pv.stamp(st, "0" * 40))

        print("\n=== Provenienz: geaenderte verfolgte Datei ===")
        with open(os.path.join(root, "a.txt"), "w", encoding="utf-8") as fh:
            fh.write("zwei\n")
        st2 = pv.repo_state(root)
        expect("geaenderte Datei -> dirty", st2["dirty"])
        expect("Stempel enthaelt DIRTY", pv.DIRTY_MARK in pv.stamp(st2))
        expect("nicht freigebbar", not pv.clean_release(st2, st2["commit"]))
        expect("die geaenderte Datei wird benannt", "a.txt" in st2["dirty_files"])

        print("\n=== Provenienz: unversionierte Datei ===")
        git(root, "checkout", "--quiet", "--", "a.txt")
        expect("nach dem Zuruecksetzen wieder sauber", not pv.repo_state(root)["dirty"])
        with open(os.path.join(root, "neu.txt"), "w", encoding="utf-8") as fh:
            fh.write("noch nicht versioniert\n")
        st3 = pv.repo_state(root)
        expect("unversionierte Datei -> dirty", st3["dirty"])
        expect("Stempel enthaelt DIRTY", pv.DIRTY_MARK in pv.stamp(st3))

        print("\n=== Provenienz: die eine erlaubte Ausnahme ===")
        # A run writes into its own output directory.  That directory, and only
        # that one, may be excluded -- otherwise every regeneration would report
        # itself as dirty.  Everything else must still count.
        import os as _os
        _os.remove(_os.path.join(root, "neu.txt"))
        _os.makedirs(_os.path.join(root, "results", "lauf_a"))
        _os.makedirs(_os.path.join(root, "results", "lauf_b"))
        for sub in ("lauf_a", "lauf_b"):
            with open(_os.path.join(root, "results", sub, "out.csv"), "w",
                      encoding="utf-8") as fh:
                fh.write("x\n")
        st_a = pv.repo_state(root, ignore=("results/lauf_a/",))
        expect("das eigene Ausgabeverzeichnis zaehlt nicht", "results/lauf_a/out.csv"
               not in st_a["dirty_files"])
        expect("es wird als ausgenommen protokolliert",
               "results/lauf_a/out.csv" in st_a["ignored_files"])
        expect("ein ANDERES Ergebnisverzeichnis zaehlt weiterhin", st_a["dirty"])
        expect("und wird benannt", "results/lauf_b/out.csv" in st_a["dirty_files"])
        st_both = pv.repo_state(root, ignore=("results/lauf_a/", "results/lauf_b/"))
        expect("beide ausgenommen -> wieder sauber", not st_both["dirty"])
        expect("die Ausnahme steht im Datensatz",
               "results/lauf_a/" in st_both["ignored_prefixes"])

        print("\n=== Provenienz: die Ausnahme deckt keinen Quelltext ===")
        with open(_os.path.join(root, "a.txt"), "w", encoding="utf-8") as fh:
            fh.write("code geaendert\n")
        st_code = pv.repo_state(root, ignore=("results/lauf_a/", "results/lauf_b/"))
        expect("geaenderte Quelldatei bleibt dirty", st_code["dirty"])
        expect("nicht freigebbar", not pv.clean_release(st_code, st_code["commit"]))
        git(root, "checkout", "--quiet", "--", "a.txt")
        import shutil as _sh
        _sh.rmtree(_os.path.join(root, "results"))

        print("\n=== Provenienz: figures_provenance.txt ===")
        st4 = pv.repo_state(root)
        expect("Baum vor dem Schreiben wieder sauber", not st4["dirty"])
        out = tempfile.mkdtemp()
        text = pv.write_provenance(out, st4, st4["commit"])
        with open(os.path.join(out, "figures_provenance.txt"), encoding="utf-8") as fh:
            body = fh.read()
        expect("Stempel wird mitgeschrieben", f"stamp={text}" in body)
        expect("Freigabestatus wird mitgeschrieben", "releasable=yes" in body)
        expect("der Kopf nennt den vollen Commit", f"head_commit={st4['commit']}" in body)
        shutil.rmtree(out, ignore_errors=True)

    print("\n=== Provenienz: kein Repository ===")
    with tempfile.TemporaryDirectory() as plain:
        st5 = pv.repo_state(plain)
        expect("kein Repository -> nicht verfuegbar", not st5["available"])
        expect("kein Repository gilt als nicht sauber", st5["dirty"])
        expect("Stempel sagt es ausdruecklich", pv.stamp(st5) == pv.NO_VERSION)
        expect("nicht freigebbar", not pv.clean_release(st5, "abc"))

    print(f"\n{'FEHLGESCHLAGEN' if failures else 'bestanden'} ({failures} Fehler)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
