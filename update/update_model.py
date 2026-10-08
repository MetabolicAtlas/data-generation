#!/usr/bin/env python3
"""Update one integrated model in a data-files checkout to a released version, and check the result.

Steps, in order (see UPDATING_MODELS.md in data-files for the full procedure):
  1. Baseline: run data-generation on the files as they are, and keep a copy of the
     current model files, so the new version can be compared with the old one.
  2. Download the model files of the release (model/*.yml, genes.tsv,
     metabolites.tsv, reactions.tsv) from GitHub into integrated-models/<model>.
     The existing YAML file name is kept (e.g. yeastGEM.yml).
  3. Prepare metaData for data-generation: a plain mapping (RAVEN 3 writes an
     !!omap), short_name set to the folder name, version and date quoted, and fields
     the previous copy had but the release lacks (short_name, full_name, authors,
     github, description, ...) carried over.
  4. Check that every TSV row has as many fields as its header. data-generation reads
     the tables by position, so a short row misplaces values. Problems are listed
     and the script stops: fix the rows by hand and rerun with --keep-files.
  5. Set version and date in integratedModels.json.
  6. Add the model's releases up to this version to its gemRepository.json.
  7. Run data-generation on the updated files, then its check/check_generated_data.py,
     which compares the output with the model files and with the baseline.

Usage, in data-generation, with data-files checked out next to it (or given with --data-files):
    python update/update_model.py --model Human-GEM --version 2.0.0
    python update/update_model.py --model Human-GEM --version 2.0.0 --keep-files
    python update/update_model.py --model Human-GEM --latest-version   # print the latest release, if newer

Exit status: 0 when every hard check passed, 1 when a check failed, 2 when the model
files need a manual fix first.
"""

import argparse
import collections
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_GENERATION = os.path.dirname(HERE)  # this checkout: data-generation, its check script and the map tools
DATA_FILES = os.path.join(os.path.dirname(DATA_GENERATION), "data-files")  # set from --data-files
TABLES = ("genes.tsv", "metabolites.tsv", "reactions.tsv")
OWNER = "SysBioChalmers"


def log(message):
    print(message, flush=True)


def fail(message, status=2):
    print(f"\nSTOP: {message}", file=sys.stderr, flush=True)
    sys.exit(status)


# ----------------------------------------------------------------------------- GitHub

def github_get(url):
    headers = {"Accept": "application/vnd.github+json", "User-Agent": "data-files-update-model"}
    token = os.environ.get("GH_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as response:
        return response.read()


def release_files(model, version):
    """{file name: download url} of the model folder at the release tag."""
    for tag in (f"v{version}", version):
        url = f"https://api.github.com/repos/{OWNER}/{model}/contents/model?ref={tag}"
        try:
            listing = json.loads(github_get(url))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                continue
            raise
        return tag, {item["name"]: item["download_url"] for item in listing if item["type"] == "file"}
    fail(f"no tag v{version} or {version} in {OWNER}/{model}")


def releases(model):
    """Published releases of the model, oldest first, as returned by fetch_release_data."""
    sys.path.insert(0, HERE)
    import fetch_release_data
    from github import Github
    return fetch_release_data.get_release_data(f"{OWNER}/{model}", Github(os.environ.get("GH_TOKEN")),
                                               is_add_version=True)


def integrated_version(model):
    index = json.load(open(os.path.join(DATA_FILES, "integrated-models", "integratedModels.json")))
    return next(e["version"] for e in index if e["short_name"] == model)


def latest_version(model):
    """The latest release if it is newer than the integrated version, otherwise None. Metabolic Atlas
    serves one version per model, so an update goes straight to the latest release."""
    from packaging.version import Version
    current = Version(integrated_version(model))
    newer = sorted(Version(r["version"]) for r in releases(model) if Version(r["version"]) > current)
    return str(newer[-1]) if newer else None


# ----------------------------------------------------------------------------- metaData

def split_metadata(text):
    """(lines before, metaData lines, lines after); metaData is the '- metaData:' item."""
    lines = text.split("\n")
    start = next((i for i, line in enumerate(lines) if re.match(r"^- metaData:", line)), None)
    if start is None:
        fail("no '- metaData:' entry at the top level of the YAML file")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("- "))
    return lines[:start], lines[start:end], lines[end:]


def parse_metadata(block):
    """Ordered (key, raw value) pairs of a metaData block in either RAVEN layout:
    RAVEN 2 '    key: value' under '- metaData:'; RAVEN 3 '  - key: value' under '- metaData: !!omap'."""
    pairs = []
    for line in block[1:]:
        if not line.strip():
            continue
        m = re.match(r"^(?:    |  - )([A-Za-z_][A-Za-z0-9_]*):(?: (.*))?$", line)
        if not m:
            fail(f"cannot read this metaData line: {line!r}")
        pairs.append((m.group(1), m.group(2) or '""'))
    return pairs


def quoted(value):
    value = value.strip()
    if value.startswith(('"', "'")):
        return value
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def merge_metadata(new_pairs, old_pairs, model, date=None):
    """Release metaData plus the fields only the previous copy had. Fields that came
    before 'version' in the previous copy go first, the others after 'date'."""
    new = dict(new_pairs)
    keys = [k for k, _ in new_pairs]
    old_keys = [k for k, _ in old_pairs]
    cut = old_keys.index("version") if "version" in old_keys else len(old_keys)
    before = [(k, v) for k, v in old_pairs[:cut] if k not in new]
    after = [(k, v) for k, v in old_pairs[cut:] if k not in new and k not in ("version", "date")]
    merged = before + [(k, new[k]) for k in keys]
    anchor = max((i for i, (k, _) in enumerate(merged) if k in ("version", "date")), default=len(merged) - 1)
    merged = merged[: anchor + 1] + after + merged[anchor + 1:]
    result = []
    for k, v in merged:
        if k == "short_name":
            v = f'"{model}"'
        elif k in ("version", "date"):
            v = quoted(f'"{date}"' if k == "date" and date else v)
        result.append((k, v))
    if "short_name" not in dict(result):
        result.insert(0, ("short_name", f'"{model}"'))
    return result


def unquote(value):
    value = value.strip()
    return value[1:-1] if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'" else value


# ----------------------------------------------------------------------------- tables

def table_lines(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [(n, row) for n, row in enumerate(csv.reader(fh, delimiter="\t"), 1)
                if row and not row[0].startswith(("#", "@"))]


def has_rows(path):
    return bool(table_lines(path))


def human_gem_columns(path):
    """Whether a model table has rows and Human-GEM's column names (rxns, rxnKEGGID, ...), which the
    map editor reads; yeast-GEM's tables start with an id column named after the annotation keys."""
    lines = table_lines(path)
    return bool(lines) and lines[0][1][0].strip('"') != "id"


def ragged_rows(path):
    """(line number, id, fields, header fields) of rows whose field count differs from the header."""
    lines = table_lines(path)
    if not lines:
        return []
    header = lines[0][1]
    return [(n, row[0], len(row), len(header)) for n, row in lines[1:] if len(row) != len(header)]


# ----------------------------------------------------------------------------- steps

def run_generation(data_generation, run_dir, data_files):
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)
    log(f"  data-generation -> {run_dir}")
    with open(os.path.join(run_dir, "generate.log"), "w") as out:
        status = subprocess.run(["node", os.path.join(data_generation, "index.js"), data_files],
                                cwd=run_dir, stdout=out, stderr=subprocess.STDOUT).returncode
    text = open(os.path.join(run_dir, "generate.log")).read()
    # data-generation can log an error and still exit with status 0
    errors = [line for line in text.splitlines() if re.match(r"^\w*Error\b", line)]
    if status != 0 or errors or not os.path.exists(os.path.join(run_dir, "neo4j", "import.cypher")):
        fail(f"data-generation failed (log in {run_dir}/generate.log):\n" + "\n".join(errors or [text[-3000:]]), 1)


def preflight(args, model_dir):
    if not os.path.isdir(model_dir):
        fail(f"{model_dir} does not exist; this script updates a model that is already integrated")
    if shutil.which("node") is None:
        fail("node is not on PATH (on Vera: module load nodejs/20.13.1-GCCcore-13.3.0)")
    if not os.path.isdir(os.path.join(args.data_generation, "node_modules")):
        fail(f"run 'yarn install --frozen-lockfile' in {args.data_generation} first")
    if not os.path.exists(os.path.join(args.data_generation, "check", "check_generated_data.py")):
        fail(f"{args.data_generation} has no check/check_generated_data.py; use a data-generation version that has it")
    utils_js = open(os.path.join(args.data_generation, "utils.js")).read()
    if "formatEcCodes" not in utils_js:
        log("  WARNING: this data-generation does not normalise EC codes or unquote TSV headers;"
            " use a version with those fixes")
    svgs = glob.glob(os.path.join(DATA_FILES, "svg", args.model, "*.svg"))
    if svgs and open(svgs[0], "rb").read(24) == b"version https://git-lfs":
        fail("the SVG maps are Git LFS pointers; run 'git lfs pull' in data-files")
    if not args.keep_files:
        changed = subprocess.run(["git", "status", "--porcelain", "--", model_dir], cwd=DATA_FILES,
                                 capture_output=True, text=True).stdout.strip()
        if changed:
            fail(f"{model_dir} has uncommitted changes; commit or discard them, or rerun with --keep-files")


def update_maps(args, model_dir, new_yaml, old_yaml, work_maps, summary_path):
    """Fit the model's SVG maps to the new version with maps/mapedit.py (rules in maps/RULES.md).

    The maps in svg/<model> are replaced by the edited ones; the change list (edited/changes.tsv), copies
    with the changes highlighted (edited/*.review.svg) and a summary stay in work_maps."""
    svg_dir = os.path.join(DATA_FILES, "svg", args.model)
    maps = sorted(glob.glob(os.path.join(svg_dir, "*.svg")))
    if not maps:
        log(f"7. Maps: no maps in svg/{args.model}")
        return
    if any(open(m, "rb").read(len(b"version https://git-lfs")) == b"version https://git-lfs" for m in maps):
        fail(f"svg/{args.model} holds Git LFS pointers; run `git lfs pull` first, or pass --skip-maps")
    out = os.path.join(work_maps, "edited")
    shutil.rmtree(work_maps, ignore_errors=True)
    os.makedirs(out)
    tables = model_dir
    if not human_gem_columns(os.path.join(model_dir, "reactions.tsv")):  # the editor's tables, written from the YAML
        tables = os.path.join(work_maps, "model-tables")
        subprocess.run([sys.executable, os.path.join(DATA_GENERATION, "maps", "yaml_to_tsv.py"), new_yaml, tables,
                        model_dir], check=True, stdout=subprocess.DEVNULL)
    boxed = any(b'class="compartment"' in open(m, "rb").read() for m in maps if b"data-transport=" not in open(m, "rb").read(3000))
    gene_label = args.gene_label or ("name" if boxed else "both")
    cmd = [sys.executable, os.path.join(DATA_GENERATION, "maps", "mapedit.py"), new_yaml, tables, old_yaml,
           "--maps", *maps, "--out", out, "--review",
           "--map-table", os.path.join(model_dir, "subsystemSVG.tsv"),
           "--compartment-table", os.path.join(model_dir, "compartmentSVG.tsv"),
           "--custom-table", os.path.join(model_dir, "customSVG.tsv"),
           "--gene-label", gene_label]
    if not boxed:
        cmd.append("--one-area")  # maps without compartment boxes draw every compartment in one area
    if args.kegg_dir:
        cmd += ["--kegg-dir", os.path.abspath(args.kegg_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode:
        log(f"7. Maps: the map editor failed; the maps are left unchanged\n{result.stderr[-3000:]}")
        return
    edited = 0
    for m in maps:
        new = os.path.join(out, os.path.basename(m))
        if os.path.exists(new) and open(new, "rb").read() != open(m, "rb").read():
            shutil.copyfile(new, m)
            edited += 1
    # transport maps are written whole from the new model, with their rows in subsystemSVG.tsv
    transport = [m for m in maps if b"data-transport=" in open(m, "rb").read(3000)]
    transport_note = ""
    if transport:
        listed = [line.split("\t")[2].strip() for line in open(os.path.join(model_dir, "subsystemSVG.tsv"), encoding="utf-8")
                  if line.count("\t") >= 2 and not line.startswith(("#", "@"))]
        template = next(os.path.join(svg_dir, f) for f in listed
                        if os.path.join(svg_dir, f) in maps and os.path.join(svg_dir, f) not in transport)
        result = subprocess.run([sys.executable, os.path.join(DATA_GENERATION, "maps", "transport_map.py"), new_yaml, tables,
                                 template, svg_dir, os.path.join(model_dir, "subsystemSVG.tsv"), "--gene-label", gene_label],
                                capture_output=True, text=True)
        if result.returncode:
            transport_note = "; transport maps left unchanged (transport_map.py failed)"
            log(f"7. Maps: transport_map.py failed\n{result.stderr[-3000:]}")
        else:
            transport_note = f"; {result.stdout.count('.svg:')} transport maps written again"
    rows = [line.rstrip("\n").split("\t") for line in open(os.path.join(out, "changes.tsv"), encoding="utf-8")][1:]
    changes = collections.Counter((r[1], r[2]) for r in rows
                                  if len(r) > 2 and r[1] not in ("review", "error") and r[2] != "connection passes")
    review = collections.Counter(r[2] for r in rows if len(r) > 2 and r[1] == "review")
    errors = [r[0] for r in rows if len(r) > 1 and r[1] == "error"]
    lines = [f"{edited} of {len(maps)} maps changed{transport_note}. The rules are listed in maps/RULES.md of "
             "[data-generation](https://github.com/MetabolicAtlas/data-generation); "
             "the workflow artifact holds the full change list and copies of the maps with the changes highlighted.", "",
             "| rule | change | count |", "|---|---|---|"]
    lines += [f"| {k[0]} | {k[1]} | {n} |" for k, n in sorted(changes.items(), key=lambda x: (x[0][0][0], int(x[0][0][1:]) if x[0][0][1:].isdigit() else 0))]
    if review:
        lines += ["", "Left for review:", "", "| item | count |", "|---|---|"]
        lines += [f"| {k} | {n} |" for k, n in review.most_common()]
    if errors:
        lines += ["", f"Maps left unchanged after an error: {', '.join(sorted(set(errors)))}"]
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"7. Maps: {edited} of {len(maps)} maps changed"
        + (f"; {len(set(errors))} left unchanged after an error" if errors else "") + transport_note
        + f"; summary in {summary_path}")


def main():
    global DATA_FILES
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="folder in integrated-models and repository name, e.g. Human-GEM")
    ap.add_argument("--version", help="release to update to, e.g. 2.0.0")
    ap.add_argument("--latest-version", action="store_true",
                    help="only print the latest release if it is newer than the integrated version (nothing if up to date)")
    ap.add_argument("--data-files", default=DATA_FILES,
                    help="data-files checkout whose model is updated (default: data-files next to data-generation)")
    ap.add_argument("--data-generation", default=DATA_GENERATION,
                    help="data-generation checkout that generates and checks the data (default: this one)")
    ap.add_argument("--work-dir",
                    help="where the baseline and generated data go (default: model-update-work next to data-files)")
    ap.add_argument("--keep-files", action="store_true",
                    help="do not download; use the model files in data-files as they are (after a manual fix)")
    ap.add_argument("--date", help="date to use instead of the one in the release YAML (YYYY-MM-DD)")
    ap.add_argument("--skip-timeline", action="store_true", help="leave gemRepository.json unchanged")
    ap.add_argument("--skip-maps", action="store_true", help="leave the SVG maps unchanged")
    ap.add_argument("--gene-label", choices=["name", "orf", "both"],
                    help="gene boxes on the maps: gene name, gene id, or name over id (default: name, or both for "
                         "maps without compartment boxes, such as Yeast-GEM's)")
    ap.add_argument("--kegg-dir", help="KGML cache made by maps/kegg_fetch.py: added reactions that a KEGG map "
                                       "of their subsystem shows keep KEGG's arrangement")
    ap.add_argument("--baseline-data-files",
                    help="generate the baseline from this data-files checkout (e.g. of main) instead of the "
                         "current files; needed with --keep-files when no earlier baseline exists")
    args = ap.parse_args()
    DATA_FILES = os.path.abspath(args.data_files)
    if not os.path.isfile(os.path.join(DATA_FILES, "integrated-models", "integratedModels.json")):
        fail(f"{DATA_FILES} is not a data-files checkout; give one with --data-files")
    if args.latest_version:
        print(latest_version(args.model) or "")
        return
    if not args.version:
        ap.error("--version is required")
    args.data_generation = os.path.abspath(args.data_generation)
    work = os.path.abspath(args.work_dir or os.path.join(os.path.dirname(DATA_FILES), "model-update-work"))
    model_dir = os.path.join(DATA_FILES, "integrated-models", args.model)

    log(f"Updating {args.model} to {args.version}")
    preflight(args, model_dir)

    # 1. baseline; the marker file ties it to this update, so a rerun with --keep-files finds it
    marker = os.path.join(work, f"{args.model}-{args.version}.baseline")
    if args.baseline_data_files:
        source = os.path.abspath(args.baseline_data_files)
        index = json.load(open(os.path.join(source, "integrated-models", "integratedModels.json")))
        old_version = next(e["version"] for e in index if e["short_name"] == args.model)
        baseline = os.path.join(work, f"{args.model}-{old_version}")
        log(f"1. Baseline: generating data from {source} ({old_version})")
        run_generation(args.data_generation, baseline, source)
        shutil.copytree(os.path.join(source, "integrated-models", args.model), os.path.join(baseline, "model-files"))
        os.makedirs(work, exist_ok=True)
        with open(marker, "w") as fh:
            fh.write(baseline + "\n")
    elif os.path.exists(marker):
        baseline = open(marker).read().strip()
        log(f"1. Baseline: reusing {baseline}")
    else:
        if args.keep_files:
            fail(f"no baseline for this update in {work}; run once without --keep-files on the unmodified files")
        old_version = integrated_version(args.model)
        baseline = os.path.join(work, f"{args.model}-{old_version}")
        log(f"1. Baseline: generating data from the current files ({old_version})")
        run_generation(args.data_generation, baseline, DATA_FILES)
        shutil.copytree(model_dir, os.path.join(baseline, "model-files"))
        with open(marker, "w") as fh:
            fh.write(baseline + "\n")
    old_yaml = glob.glob(os.path.join(baseline, "model-files", "*.yml"))[0]

    # 2. download
    yaml_files = glob.glob(os.path.join(model_dir, "*.yml"))
    if len(yaml_files) != 1:
        fail(f"expected one .yml file in {model_dir}")
    target_yaml = yaml_files[0]
    if args.keep_files:
        log("2. Download: skipped (--keep-files)")
    else:
        tag, files = release_files(args.model, args.version)
        sources = [n for n in files if n.endswith((".yml", ".yaml"))]
        if len(sources) != 1:
            fail(f"expected one YAML file in model/ at {tag}, found {sources}")
        # everything is downloaded before anything is written, so a failure leaves the folder unchanged
        downloads = {target_yaml: github_get(files[sources[0]])}
        kept = []
        for name in TABLES:
            path = os.path.join(model_dir, name)
            if name in files:
                downloads[path] = github_get(files[name])
            elif os.path.exists(path) and not has_rows(path):
                kept.append(name)  # the model has no such table; the empty placeholder stays
            else:
                fail(f"{name} is missing in model/ at {tag}, but integrated-models/{args.model}/{name} has content")
        for path, content in downloads.items():
            with open(path, "wb") as fh:
                fh.write(content)
        log(f"2. Download: {tag} {sources[0]} -> {os.path.basename(target_yaml)}"
            + "".join(f", {n}" for n in TABLES if n not in kept)
            + (f" (no {', '.join(kept)} in the release; empty placeholder kept)" if kept else ""))

    # 3. metaData
    text = open(target_yaml, encoding="utf-8").read()
    head, block, rest = split_metadata(text)
    pairs = parse_metadata(block)
    _, old_block, _ = split_metadata(open(old_yaml, encoding="utf-8").read())
    merged = merge_metadata(pairs, parse_metadata(old_block), args.model, args.date)
    version = unquote(dict(merged)["version"])
    date = unquote(dict(merged).get("date", '""'))
    if version != args.version:
        fail(f"the YAML says version {version}, not {args.version}")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date):
        fail(f"metaData date {date!r} is not YYYY-MM-DD; pass --date")
    new_block = ["- metaData:"] + [f"    {k}: {v}" for k, v in merged]
    with open(target_yaml, "w", encoding="utf-8") as fh:
        fh.write("\n".join(head + new_block + rest))
    log(f"3. metaData: {', '.join(k for k, _ in merged)}")
    if re.search(r"^\s+- ec-code:", text, re.M) and not re.search(r"^\s+- eccodes:", text, re.M):
        if "ec-code" not in open(os.path.join(args.data_generation, "utils.js")).read():
            log("  WARNING: EC codes are under annotation/ec-code, which this data-generation does not read;"
                " the site will show no EC numbers")

    # 4. tables
    problems = []
    for name in TABLES:
        for line, row_id, n, expected in ragged_rows(os.path.join(model_dir, name)):
            problems.append(f"  {name}:{line} {row_id}: {n} fields, header has {expected}")
    if problems:
        fail("rows with the wrong number of fields (fix them in integrated-models/"
             f"{args.model}, check the model's develop branch for a corrected row, then rerun with --keep-files):\n"
             + "\n".join(problems))
    log("4. Tables: every row has as many fields as its header")

    # 5. integratedModels.json
    index_path = os.path.join(DATA_FILES, "integrated-models", "integratedModels.json")
    index_text = open(index_path, encoding="utf-8").read()
    entries = json.loads(index_text)
    position = next((i for i, e in enumerate(entries) if e["short_name"] == args.model), None)
    if position is None:
        fail(f"{args.model} is not in integratedModels.json")
    # edit in place, so the rest of the file keeps its layout
    pattern = re.compile(r'(\{\s*"short_name":\s*"' + re.escape(args.model) + r'".*?\n  \})', re.S)
    entry_text = pattern.search(index_text).group(1)
    updated = re.sub(r'"version":\s*"[^"]*"', f'"version": "{version}"', entry_text, count=1)
    updated = re.sub(r'"date":\s*"[^"]*"', f'"date": "{date}"', updated, count=1)
    with open(index_path, "w", encoding="utf-8") as fh:
        fh.write(index_text.replace(entry_text, updated, 1))
    log(f"5. integratedModels.json: version {version}, date {date}")

    # 6. timeline
    if args.skip_timeline:
        log("6. Timeline: skipped")
    else:
        found = releases(args.model)
        for r in found:
            r.pop("version", None)
        ids = [r["id"] for r in found]
        wanted = f"{args.model}-{version}"
        if wanted not in ids:
            fail(f"release {wanted} not found on GitHub (found up to {ids[-1] if ids else 'none'})")
        timeline = found[: ids.index(wanted) + 1]
        with open(os.path.join(model_dir, "gemRepository.json"), "w", encoding="utf-8") as fh:
            fh.write(json.dumps(timeline, indent=2))
        log(f"6. Timeline: {len(timeline)} releases, last {wanted}")

    # 7. maps
    maps_dir = os.path.join(work, f"{args.model}-{version}-maps")  # kept apart: step 8 empties its own folder
    maps_summary = os.path.join(maps_dir, "maps_summary.md")
    if args.skip_maps:
        log("7. Maps: skipped")
    else:
        update_maps(args, model_dir, target_yaml, old_yaml, maps_dir, maps_summary)

    # 8. generate and check
    log("8. Generating and checking the new version")
    run_dir = os.path.join(work, f"{args.model}-{version}")
    run_generation(args.data_generation, run_dir, DATA_FILES)
    report = os.path.join(run_dir, "check_report.md")
    check = [sys.executable, os.path.join(args.data_generation, "check", "check_generated_data.py"), "--model", args.model,
             "--data-files", DATA_FILES, "--data-generation", args.data_generation, "--new", run_dir,
             "--old", baseline, "--old-model-dir", os.path.join(baseline, "model-files"), "--report", report]
    metabolicatlas = os.path.join(os.path.dirname(DATA_FILES), "MetabolicAtlas")
    if os.path.isdir(metabolicatlas):
        check += ["--metabolicatlas", metabolicatlas]
    status = subprocess.run(check).returncode
    log(f"\nReport: {report}")
    log("Next: read the report, commit the changes in data-files, and run the Docker checks (UPDATING_MODELS.md in "
        "data-files).")
    sys.exit(1 if status else 0)


if __name__ == "__main__":
    main()
