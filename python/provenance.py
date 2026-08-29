#!/usr/bin/env python3
"""Provenance of a figure: which commit it was produced from, and whether the
working tree was clean at the time.

A figure that names a commit is making a claim: "run this commit and you get
this picture".  The claim is only true if the tree was clean when the data and
the figure were made, and if the commit named is the one that was actually
checked out.  Neither can be taken from a command line argument -- an argument
says what the caller believed, and the P2a figures once carried a commit id
that predated the code they showed.  So it is measured here, from git, at the
moment the figure is drawn.

Rules, in order:

  * git missing, or the path is not a repository, or HEAD has no commit
        -> "KEINE VERSIONSANGABE".  Not clean, not publishable.
  * working tree dirty (tracked modifications OR untracked files)
        -> "<short> DIRTY".  Untracked files count: an uncommitted script that
           the figure depends on is exactly the case this must catch.
  * a data commit is supplied and differs from HEAD
        -> "<short> DATEN AUS <short> NICHT HEAD".
  * otherwise
        -> "<short>".

ONE EXCLUSION, and only one: the output directory this very run is writing.
Its files are the run's own product, so requiring them to be committed before
the run has made them is circular -- and it would make every regeneration
self-reporting as DIRTY, which would drain the mark of meaning.  Everything
else -- sources, scripts, configuration, documentation, and every OTHER results
directory -- counts.  The excluded paths are listed in the provenance record so
the exclusion is visible rather than assumed.

A consequence worth stating: a figure cannot name a commit that contains the
figure, because the hash is not known until the commit exists.  The stamp
therefore names the commit of the CODE the figure was produced from; the
artefacts are committed afterwards, in a commit that contains no code.

`clean_release(...)` is the single predicate a release step should ask.
"""
import os
import shutil
import subprocess

NO_VERSION = "KEINE VERSIONSANGABE"
DIRTY_MARK = "DIRTY"


def _git(repo_root, *args):
    """Run git in `repo_root`; return stdout, or None if it cannot be run."""
    exe = shutil.which("git")
    if exe is None:
        return None
    try:
        out = subprocess.run([exe, "-C", str(repo_root)] + list(args),
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    return out.stdout


def repo_state(repo_root, ignore=()):
    """Commit and cleanliness of the working tree at `repo_root`.

    `ignore` is a sequence of repository-relative path prefixes whose changes do
    not count -- in practice the one output directory this run is writing.  See
    the module docstring for why that exclusion exists and why it is the only
    one.  Prefixes are matched against git's own forward-slash paths.

    Returns a dict with keys: commit, short, dirty, dirty_files, ignored_files,
    ignored_prefixes, available.  `dirty` is True whenever cleanliness could not
    be established, so an unusable repository can never be mistaken for a clean
    one.
    """
    pre = tuple(str(x).replace(os.sep, "/").lstrip("./") for x in ignore)
    state = {"commit": None, "short": None, "dirty": True, "dirty_files": [],
             "ignored_files": [], "ignored_prefixes": list(pre), "available": False}
    head = _git(repo_root, "rev-parse", "HEAD")
    if head is None or not head.strip():
        return state
    # --untracked-files=all, not "normal": normal collapses an untracked
    # directory into one entry, and a path prefix can then neither be matched
    # against it nor reported.
    status = _git(repo_root, "status", "--porcelain", "--untracked-files=all")
    if status is None:
        return state
    files = [line[3:].strip().strip('"') for line in status.splitlines() if line.strip()]
    kept = [f for f in files if not any(f.startswith(x) for x in pre)]
    state["commit"] = head.strip()
    state["short"] = head.strip()[:12]
    state["dirty_files"] = kept
    state["ignored_files"] = [f for f in files if f not in kept]
    state["dirty"] = bool(kept)
    state["available"] = True
    return state


def stamp(state, data_commit=None):
    """The string that goes into the figure."""
    if not state["available"]:
        return NO_VERSION
    text = state["short"]
    if state["dirty"]:
        return f"{text} {DIRTY_MARK} ({len(state['dirty_files'])} Datei(en) unversioniert)"
    if data_commit and data_commit.strip() and data_commit.strip() != state["commit"]:
        return f"{text} DATEN AUS {data_commit.strip()[:12]} NICHT HEAD"
    return text


def clean_release(state, data_commit=None):
    """May this run be released?  Clean tree, and data from exactly this commit."""
    if not state["available"] or state["dirty"]:
        return False
    if data_commit is None or not data_commit.strip():
        return False
    return data_commit.strip() == state["commit"]


def repo_root_of(path):
    """Repository root containing `path` -- here, the parent of python/."""
    return os.path.dirname(os.path.dirname(os.path.abspath(path)))


def write_provenance(directory, state, data_commit=None):
    """Record what the figures in `directory` were stamped with."""
    text = stamp(state, data_commit)
    with open(os.path.join(directory, "figures_provenance.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(f"stamp={text}\n")
        fh.write(f"head_commit={state['commit'] or ''}\n")
        fh.write(f"data_commit={(data_commit or '').strip()}\n")
        fh.write(f"working_tree_dirty={'yes' if state['dirty'] else 'no'}\n")
        fh.write(f"releasable={'yes' if clean_release(state, data_commit) else 'no'}\n")
        for x in state.get("ignored_prefixes", []):
            fh.write(f"ignored_prefix={x}\n")
        for f in state["dirty_files"]:
            fh.write(f"dirty_file={f}\n")
        for f in state.get("ignored_files", []):
            fh.write(f"ignored_file={f}\n")
    return text


if __name__ == "__main__":
    import sys
    root = repo_root_of(__file__)
    st = repo_state(root)  # nothing ignored when asked directly
    print(f"repo   : {root}")
    print(f"commit : {st['commit']}")
    print(f"dirty  : {st['dirty']} ({len(st['dirty_files'])} Datei(en))")
    print(f"stamp  : {stamp(st, sys.argv[1] if len(sys.argv) > 1 else None)}")
