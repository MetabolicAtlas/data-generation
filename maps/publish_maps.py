#!/usr/bin/env python3
"""Write a model's maps into a checkout of its maps repository (SysBioChalmers/Human-maps, Yeast-maps).

The subsystem and compartment maps listed in integrated-models/<Model>/subsystemSVG.tsv and
compartmentSVG.tsv are written as SVG, SBGN-ML, SBML, Escher (JSON) and PNG (export_formats.py) into
<repo>/subsystem/<format>/ and <repo>/compartment/<format>/; maps that are no longer listed are removed.
The README's line "The maps match **<Model> <version>**." gets the version in integratedModels.json.
A summary for the pull request is written to --summary.

Usage: publish_maps.py --model Human-GEM --repo <maps repository checkout> [--data-files <checkout>]
                       [--summary summary.md]
Used by the publish-maps workflow of MetabolicAtlas/data-files; see UPDATING_MODELS.md there.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FORMATS = ("svg", "sbgn", "sbml", "png", "escher")


def listed_maps(table):
    out = []
    for line in open(table, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and not line.startswith(("#", "@")) and f[2].endswith(".svg"):
            out.append(f[2])
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, help="folder in integrated-models, e.g. Human-GEM")
    ap.add_argument("--repo", required=True, help="checkout of the model's maps repository")
    ap.add_argument("--data-files", default=os.path.join(os.path.dirname(os.path.dirname(HERE)), "data-files"),
                    help="data-files checkout with the maps (default: data-files next to data-generation)")
    ap.add_argument("--summary", help="write a Markdown summary for the pull request here")
    a = ap.parse_args()
    DATA_FILES = os.path.abspath(a.data_files)
    model_dir = os.path.join(DATA_FILES, "integrated-models", a.model)
    svg_dir = os.path.join(DATA_FILES, "svg", a.model)
    index = json.load(open(os.path.join(DATA_FILES, "integrated-models", "integratedModels.json")))
    version = next(e["version"] for e in index if e["short_name"] == a.model)
    yaml_file = [f for f in os.listdir(model_dir) if f.endswith(".yml")][0]
    work = tempfile.mkdtemp(prefix="maps-")
    tables = model_dir
    reactions = os.path.join(model_dir, "reactions.tsv")
    if not os.path.exists(reactions) or sum(1 for _ in open(reactions)) < 2:  # releases without TSV files
        tables = os.path.join(work, "model-tables")
        subprocess.run([sys.executable, os.path.join(HERE, "yaml_to_tsv.py"), os.path.join(model_dir, yaml_file), tables],
                       check=True, stdout=subprocess.DEVNULL)
    counts, removed = {}, {}
    for kind in ("subsystem", "compartment"):
        maps = listed_maps(os.path.join(model_dir, f"{kind}SVG.tsv"))
        missing = [m for m in maps if not os.path.exists(os.path.join(svg_dir, m))]
        if missing:
            sys.exit(f"maps listed in {kind}SVG.tsv but not in svg/{a.model}: {', '.join(missing)}")
        out = os.path.join(work, kind)
        if maps:
            result = subprocess.run([sys.executable, os.path.join(HERE, "export_formats.py"),
                                     os.path.join(model_dir, yaml_file), tables, out, "--model-name", a.model,
                                     "--maps", *[os.path.join(svg_dir, m) for m in maps]],
                                    capture_output=True, text=True)
            print(result.stdout[-2000:], result.stderr[-2000:], sep="\n")
            last = (result.stdout.strip().splitlines() or [""])[-1]
            if result.returncode or not last.endswith("errors: 0"):
                sys.exit(f"export of the {kind} maps failed or reported errors: {last}")
        names = {m[:-4] for m in maps}
        gone = set()
        for fmt in FORMATS:
            target = os.path.join(a.repo, kind, fmt)
            if os.path.isdir(target):
                gone |= {os.path.splitext(f)[0] for f in os.listdir(target)} - names
                shutil.rmtree(target)
            if maps:
                shutil.copytree(os.path.join(out, fmt), target)
        counts[kind], removed[kind] = len(maps), sorted(gone)
    readme = os.path.join(a.repo, "README.md")
    if os.path.exists(readme):
        text = open(readme, encoding="utf-8").read()
        text = re.sub(r"(The maps match \*\*)[^*]+(\*\*)", rf"\g<1>{a.model} {version}\g<2>", text)
        open(readme, "w", encoding="utf-8").write(text)
    shutil.rmtree(work)
    lines = [f"Maps of **{a.model} {version}**, written from MetabolicAtlas/data-files (`svg/{a.model}`) by "
             "`maps/publish_maps.py` of MetabolicAtlas/data-generation, as SVG, SBGN-ML 0.3, SBML Level 3 with layout and groups, Escher (JSON) "
             "and PNG.", "",
             f"- {counts['subsystem']} subsystem maps and {counts['compartment']} compartment maps"]
    for kind in ("subsystem", "compartment"):
        if removed[kind]:
            lines.append(f"- {kind} maps removed (no longer in the model's map table): {', '.join(removed[kind])}")
    summary = "\n".join(lines) + "\n"
    print(summary)
    if a.summary:
        open(a.summary, "w", encoding="utf-8").write(summary)


if __name__ == "__main__":
    main()
