#!/usr/bin/env python3
"""Update a maps repository (SysBioChalmers/Human-maps, Yeast-maps) to a newer release of its model.

The README says which release the maps match ("The maps match **Human-GEM 2.1.1**."). For a newer release, this
script downloads the release's model file and annotation tables, and the YAML file of the release the maps match, then
  1. adapts subsystemSVG.tsv and compartmentSVG.tsv to the subsystems and compartments of the release (maplists.py),
  2. fits every map to the release with mapedit.py (rules in RULES.md),
  3. writes the transport maps again from the release (transport_map.py),
  4. writes every listed map as SVG, SBGN-ML and Escher (JSON) into subsystem/ and compartment/ (export_formats.py),
     removing the maps that are no longer listed,
  5. sets the version in the README.
A summary for the pull request is written to --summary.

Usage: update_maps.py --model Human-GEM --repo <maps repository checkout> [--version 2.1.2] [--force]
                      [--work-dir <folder>] [--summary summary.md]
       update_maps.py --model Human-GEM --repo <checkout> --latest-version
Without --version the latest release is used; --latest-version only prints it if it is newer than the one in the
README. --force runs again for the release that the README names.
Exit status: 0 when the maps are up to date, 1 when the update failed and nothing was written, 2 when the maps were
written but a map or the transport maps could not be fitted. Used by the update-maps workflow of this repository.
"""

import argparse
import collections
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "update"))
sys.path.insert(0, HERE)
import maplists  # noqa: E402
import update_model as um  # noqa: E402

FORMATS = ("svg", "sbgn", "escher")
LISTS = ("subsystemSVG.tsv", "compartmentSVG.tsv")
MATCH = re.compile(r"The maps match \*\*([^*\d]*?)\s*(\d[^*\s]*)\*\*")


def readme_version(repo, model):
    """The release that the README of the maps repository names; stops if it names another model."""
    with open(os.path.join(repo, "README.md"), encoding="utf-8") as fh:
        found = MATCH.search(fh.read())
    if found is None:
        um.fail(f'{repo}/README.md has no line "The maps match **<Model> <version>**."', 1)
    if found.group(1) != model:
        um.fail(f"{repo}/README.md says that its maps match {found.group(1)}, not {model}", 1)
    return found.group(2)


def latest_version(model, current):
    """The latest release if it is newer than `current`, otherwise None."""
    from packaging.version import Version
    newer = sorted(Version(r["version"]) for r in um.releases(model) if Version(r["version"]) > Version(current))
    return str(newer[-1]) if newer else None


def download(model, version, folder, tables=True):
    """The release's YAML file, and with `tables` its annotation tables, in folder; returns the YAML path."""
    try:
        tag, files = um.release_files(model, version)
    except SystemExit:  # release_files stops with status 2, which here means that maps were written
        sys.exit(1)
    sources = [n for n in files if n.endswith((".yml", ".yaml"))]
    if len(sources) != 1:
        um.fail(f"expected one YAML file in model/ at {tag}, found {sources}", 1)
    missing = [t for t in um.TABLES if t not in files] if tables else []
    if missing:
        um.fail(f"{model} {version} has no {', '.join(missing)} in model/; the maps cannot be fitted without the "
                "annotation tables", 1)
    os.makedirs(folder, exist_ok=True)
    for name in (sources[0], *(um.TABLES if tables else ())):
        with open(os.path.join(folder, name), "wb") as fh:
            fh.write(um.github_get(files[name]))
    return os.path.join(folder, sources[0])


def run(cmd, what):
    result = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "PYTHONHASHSEED": "0"})
    if result.returncode:
        um.fail(f"{what} failed:\n{(result.stderr or result.stdout)[-3000:]}", 1)
    return result.stdout


def listed_maps(table):
    return sorted({r[3] for r in maplists.read_table(table)[1] if r[3].endswith(".svg")})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="Human-GEM or Yeast-GEM: the model repository in SysBioChalmers")
    ap.add_argument("--repo", required=True, help="checkout of the model's maps repository")
    ap.add_argument("--version", help="release to fit the maps to (default: the latest)")
    ap.add_argument("--latest-version", action="store_true",
                    help="only print the latest release if it is newer than the one in the README")
    ap.add_argument("--force", action="store_true", help="run again for the release that the README names")
    ap.add_argument("--work-dir", help="empty or new folder for the downloads and intermediate maps (default: a "
                                       "temporary folder)")
    ap.add_argument("--kegg-dir", help="KGML cache made by maps/kegg_fetch.py: added reactions that a KEGG map of "
                                       "their subsystem shows keep KEGG's arrangement")
    ap.add_argument("--summary", help="write a Markdown summary for the pull request here")
    a = ap.parse_args()
    repo = os.path.abspath(a.repo)
    current = readme_version(repo, a.model)
    if a.latest_version:
        print(latest_version(a.model, current) or "")
        return
    version = a.version or latest_version(a.model, current)
    if version is None or (version == current and not a.force):
        print(f"The maps match {a.model} {current}; nothing to do.")
        return
    from packaging.version import Version
    if Version(version) < Version(current):
        um.fail(f"{version} is older than the release {current} that the maps match", 1)
    work = os.path.abspath(a.work_dir or tempfile.mkdtemp(prefix="update-maps-"))
    if os.path.exists(work) and os.listdir(work):
        um.fail(f"the work folder {work} is not empty", 1)
    print(f"Fitting the {a.model} maps from {current} to {version}", flush=True)

    new_yml = download(a.model, version, os.path.join(work, "new"))
    old_yml = download(a.model, current, os.path.join(work, "old"), tables=False)
    model_dir = os.path.dirname(new_yml)

    # every map in one folder (file names are unique over both kinds), with its own copies of the lists
    svg_dir = os.path.join(work, "maps")
    os.makedirs(svg_dir)
    for kind in maplists.KINDS:
        for path in glob.glob(os.path.join(repo, kind, "svg", "*.svg")):
            if os.path.exists(os.path.join(svg_dir, os.path.basename(path))):
                um.fail(f"{os.path.basename(path)} is a subsystem map and a compartment map", 1)
            shutil.copy(path, svg_dir)
    lists = {name: shutil.copy(os.path.join(repo, name), os.path.join(work, name)) for name in LISTS}
    report = maplists.adapt(old_yml, new_yml, lists["subsystemSVG.tsv"], lists["compartmentSVG.tsv"], svg_dir)

    maps = sorted(glob.glob(os.path.join(svg_dir, "*.svg")))
    transport = [m for m in maps if b"data-transport=" in open(m, "rb").read(3000)]
    boxed = any(b'class="compartment"' in open(m, "rb").read() for m in maps if m not in transport)
    gene_label = "name" if boxed else "both"  # Human-GEM's maps have compartment boxes, Yeast-GEM's do not
    edited = os.path.join(work, "edited")
    cmd = [sys.executable, os.path.join(HERE, "mapedit.py"), new_yml, model_dir, old_yml, "--maps", *maps,
           "--out", edited, "--map-table", lists["subsystemSVG.tsv"], "--compartment-table", lists["compartmentSVG.tsv"],
           "--gene-label", gene_label] + ([] if boxed else ["--one-area"])
    if a.kegg_dir:
        cmd += ["--kegg-dir", os.path.abspath(a.kegg_dir)]
    run(cmd, "the map editor (mapedit.py)")
    changed = 0
    for m in maps:
        new = os.path.join(edited, os.path.basename(m))
        if os.path.exists(new) and open(new, "rb").read() != open(m, "rb").read():
            shutil.copyfile(new, m)
            changed += 1

    problems = []
    if transport:
        listed = [os.path.join(svg_dir, f) for f in listed_maps(lists["subsystemSVG.tsv"])]
        template = next((f for f in listed if os.path.exists(f) and f not in transport), None)
        if template is None:
            problems.append("the transport maps could not be written: no other subsystem map to use as template")
        else:
            try:
                run([sys.executable, os.path.join(HERE, "transport_map.py"), new_yml, model_dir, template, svg_dir,
                     lists["subsystemSVG.tsv"], "--gene-label", gene_label], "the transport maps (transport_map.py)")
            except SystemExit as err:
                problems.append("the transport maps (transport_map.py) failed; some of them may be incomplete")
                print(err, file=sys.stderr)
    rows = [line.rstrip("\n").split("\t") for line in open(os.path.join(edited, "changes.tsv"), encoding="utf-8")][1:]
    errors = sorted({r[0] for r in rows if len(r) > 1 and r[1] == "error"})
    if errors:
        problems.append(f"{len(errors)} maps were left unchanged after an error: {', '.join(errors)}")
    changes = collections.Counter((r[1], r[2]) for r in rows
                                  if len(r) > 2 and r[1] not in ("review", "error") and r[2] != "connection passes")
    review = collections.Counter(r[2] for r in rows if len(r) > 2 and r[1] == "review")

    # export the listed maps; the repository is changed only after every export has succeeded
    if not os.path.exists(os.path.join(HERE, "schema", "SBGN.xsd")):
        print("Warning: maps/schema/SBGN.xsd is missing; the SBGN-ML files are not validated.", flush=True)
    exported = {}
    for kind in maplists.KINDS:
        files = listed_maps(lists[f"{kind}SVG.tsv"])
        missing = [f for f in files if not os.path.exists(os.path.join(svg_dir, f))]
        if missing:
            um.fail(f"maps listed in {kind}SVG.tsv but missing: {', '.join(missing)}", 1)
        out = os.path.join(work, "export", kind)
        if files:
            output = run([sys.executable, os.path.join(HERE, "export_formats.py"), new_yml, model_dir, out,
                          "--model-name", a.model, "--maps", *[os.path.join(svg_dir, f) for f in files]],
                         f"the export of the {kind} maps")
            last = (output.strip().splitlines() or [""])[-1]
            if not last.endswith("errors: 0"):
                um.fail(f"the export of the {kind} maps reported errors: {last}", 1)
        exported[kind] = (out, files)
    for kind, (out, files) in exported.items():
        for fmt in FORMATS:
            if os.path.isdir(os.path.join(repo, kind, fmt)):
                shutil.rmtree(os.path.join(repo, kind, fmt))
            if files:
                shutil.copytree(os.path.join(out, fmt), os.path.join(repo, kind, fmt))
    for name, path in lists.items():
        shutil.copyfile(path, os.path.join(repo, name))
    if report["renamed"]:
        problems.append("some subsystems may have been renamed; their maps are replaced by blank ones (see the lists "
                        "of maps below)")
    if not problems:  # the README names the release only for maps that are fully fitted to it
        readme = os.path.join(repo, "README.md")
        with open(readme, encoding="utf-8") as fh:
            text = fh.read()
        with open(readme, "w", encoding="utf-8") as fh:  # also the version named in the table of formats
            fh.write(text.replace(f"{a.model} {current}", f"{a.model} {version}"))

    lines = [f"**The maps are NOT fully up to date: {'; '.join(problems)}. The README still names {current}; after "
             f"the maps are corrected, change it to {version}.**", ""] if problems else []
    lines += [f"Maps of **{a.model} {version}**, fitted from the maps of {current} with `maps/update_maps.py` of "
              f"[MetabolicAtlas/data-generation](https://github.com/MetabolicAtlas/data-generation), as SVG, SBGN-ML 0.3 "
              f"and Escher (JSON). {changed} of {len(maps)} maps changed; the rules are in "
              "[maps/RULES.md](https://github.com/MetabolicAtlas/data-generation/blob/main/maps/RULES.md).", "",
              f"{len(exported['subsystem'][1])} subsystem maps and {len(exported['compartment'][1])} compartment maps."]
    described = maplists.describe(report)
    if described:
        lines += ["", "Changes to the lists of maps:", ""] + described
    if changes:
        lines += ["", "| rule | change | count |", "|---|---|---|"]
        lines += [f"| {k[0]} | {k[1]} | {n} |" for k, n in sorted(
            changes.items(), key=lambda x: (x[0][0][0], int(x[0][0][1:]) if x[0][0][1:].isdigit() else 0))]
    if review:
        lines += ["", "Left for review:", "", "| item | count |", "|---|---|"]
        lines += [f"| {k} | {n} |" for k, n in review.most_common()]
    summary = "\n".join(lines)[:60000] + "\n"  # GitHub refuses a description over 65536 characters
    print(summary)
    if a.summary:
        open(a.summary, "w", encoding="utf-8").write(summary)
    sys.exit(2 if problems else 0)


if __name__ == "__main__":
    main()
