#!/usr/bin/env python3
"""Write the pull-request description for a model update prepared by update_model.py.

Used by the update-model workflow of MetabolicAtlas/data-files. The description says what the update
script did, what (if anything) needs a manual fix, summarises the check report, and
gives the steps to test the update in a local deployment.

Usage:
    python update_pr_body.py --model Human-GEM --version 2.0.0 --status 0 --log update.log \\
        --report check_report.md --maps-summary maps_summary.md --branch auto/update-human-gem-2.0.0 \\
        --data-generation-ref main --run-url <url> > body.md
"""

import argparse
import os
import re

# GitHub refuses descriptions over 65536 characters; the full report is in the run's artifact
MAX_BODY = 60000

STATUS = {
    0: "All hard checks of the generated data passed. Next: test the update in a local deployment (below).",
    1: "**Some hard checks failed or data-generation stopped.** See the script output and the failed checks "
       "below; [UPDATING_MODELS.md]({guide}) lists the usual causes and fixes.",
    2: "**The update stopped: the model files need a manual fix first.** The script output below says what "
       "to fix; [UPDATING_MODELS.md]({guide}) explains how.",
}


def section(report, title):
    """The lines of a '## title...' section of the report, without the heading."""
    match = re.search(rf"^## {re.escape(title)}[^\n]*\n(.*?)(?=^## |\Z)", report, re.S | re.M)
    return match.group(1).strip() if match else ""


def details(summary, body):
    return f"<details><summary>{summary}</summary>\n\n{body}\n\n</details>\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True)
    ap.add_argument("--version", required=True)
    ap.add_argument("--status", type=int, required=True, help="exit status of update_model.py")
    ap.add_argument("--log", required=True, help="output of update_model.py")
    ap.add_argument("--report", help="check_report.md, when the script got that far")
    ap.add_argument("--maps-summary", help="maps_summary.md written by the maps step of update_model.py")
    ap.add_argument("--branch", required=True, help="branch of this pull request")
    ap.add_argument("--base", default="main", help="branch the pull request targets")
    ap.add_argument("--data-generation-ref", default="main")
    ap.add_argument("--run-url", default="")
    ap.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", "MetabolicAtlas/data-files"))
    args = ap.parse_args()

    log = open(args.log, encoding="utf-8").read().strip() if os.path.exists(args.log) else "(no output)"
    report = open(args.report, encoding="utf-8").read() if args.report and os.path.exists(args.report) else ""
    guide = f"https://github.com/{args.repo}/blob/{args.base}/UPDATING_MODELS.md"
    status = STATUS.get(args.status, STATUS[1]).format(guide=guide)
    run = f" ([workflow run]({args.run_url}); the full check report is in its artifacts)" if args.run_url else ""

    out = [f"Updates **{args.model}** to [{args.version}](https://github.com/SysBioChalmers/{args.model}"
           f"/releases/tag/v{args.version}). Prepared by the *Update integrated model* workflow{run}.", "", status, ""]

    if args.status != 0:
        out += ["### What to fix", "",
                f"1. Fix the files on this branch (`{args.branch}`), following the script output below.",
                "2. Run the *Update integrated model* workflow again (Actions tab, *Run workflow*) with "
                f"`continue_branch` set to `{args.branch}`. It keeps your fixes, checks again and updates "
                "this description.", ""]

    out += ["### Script output", "", "```", log[-6000:], "```", ""]

    if report:
        result = next((line for line in report.splitlines() if line.startswith("**Result")), "")
        flagged = [line for line in report.splitlines() if line.startswith(("- FAIL", "- WARN"))]
        out += ["### Check of the generated data", "", result, ""]
        out += flagged + [""] if flagged else []
    optional = []
    if args.maps_summary and os.path.exists(args.maps_summary):
        optional.append(details("Map changes", open(args.maps_summary, encoding="utf-8").read().strip()))
    for title, summary in (("Changes compared with", "Changes compared with the integrated version"),
                           ("Data overlay", "Data overlay coverage"),
                           ("Reactions drawn on the maps", "Reactions drawn on the maps that are not in the model"),
                           ("Identifiers used in MetabolicAtlas tests", "Identifiers used in the MetabolicAtlas tests")):
        text = section(report, title)
        if text:
            optional.append(details(summary, text))

    test = [
        "### Test it locally before merging", "",
        f"This needs Docker; see part 2 of [UPDATING_MODELS.md]({guide}) for details.", "",
        "```bash",
        "git clone https://github.com/MetabolicAtlas/MetabolicAtlas",
        "git clone https://github.com/MetabolicAtlas/data-generation",
        "git clone https://github.com/MetabolicAtlas/data-files",
        "(cd data-files && git lfs install && git lfs pull)",
        "(cd data-generation && yarn install --frozen-lockfile)",
        "cd MetabolicAtlas && cp env-local.env.sample env-local.env   # set NEO4J_PASSWORD and POSTGRES_PASSWORD",
        "source proj.sh",
        "",
        "# 1. baseline: the integrated version",
        "build-stack && start-stack",
        "ma-exec api yarn test 2>&1 | tee ../test-before.log",
        "",
        "# 2. this update",
        f"git -C ../data-files switch {args.branch} && git -C ../data-files lfs pull",
        f"git -C ../data-generation switch {args.data_generation_ref}",
        "stop-stack",
        "build-stack && start-stack   # build-stack also imports the data into Neo4j",
        "ma-exec api yarn test 2>&1 | tee ../test-after.log",
        "xvfb-run npx --yes cypress run --project frontend   # without xvfb-run when a display is available",
        "```", "",
        "Do not run `import-neo4j-db` on a filled database (it times out) or `clean-stack` (it also deletes the "
        "GotEnzymes database). Some API tests fail before the update too, so compare the two logs; only failures new "
        "in `test-after.log` matter. Each one is either an expected change (a count, version or identifier that "
        "changed with this release) or a problem. If every suite fails before any test runs, Jest cannot load "
        "`node-fetch` 3; see part 2 of the guide.", "",
        f"Then, on http://localhost with {args.model} {args.version} selected:", "",
        "- [ ] API tests: no unexplained new failures (expected value changes go to a MetabolicAtlas pull request)",
        "- [ ] Cypress: no failure that depends on the data (some tests are timing-dependent; rerun a failing spec)",
        "- [ ] The model list shows the new version and date",
        "- [ ] A reaction with several EC numbers shows each as its own link, and lists its references",
        "- [ ] A gene page shows its cross-references (Ensembl, UniProt, NCBI Gene, Protein Atlas)",
        "- [ ] A metabolite page shows name, formula and charge",
        "- [ ] A subsystem map and a compartment map open, and clicking a reaction opens its panel",
        "- [ ] The 3D viewer opens for a compartment",
        "- [ ] Each data overlay colours the maps",
        "- [ ] Search finds a gene symbol and an EC number",
        "- [ ] Compare Models works with this model",
        "- [ ] The GEM repository timeline ends at this version", "",
        "Tick the boxes as you go, and note anything unexpected in a comment. Merge when the checks above pass.",
    ]
    # the collapsible report sections are left out, largest last, when the description would get too long
    head, tail = "\n".join(out), "\n".join(test)
    kept, skipped = [], 0
    for block in optional:
        if len(head) + len(tail) + sum(map(len, kept)) + len(block) < MAX_BODY:
            kept.append(block)
        else:
            skipped += 1
    if skipped:
        kept.append(f"*{skipped} report section(s) left out for length; the full report is in the workflow run's "
                    "artifacts.*\n")
    print("\n".join([head] + kept + [tail]))


if __name__ == "__main__":
    main()
