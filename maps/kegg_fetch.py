#!/usr/bin/env python3
"""Find the KEGG reference maps of each subsystem and cache their KGML files, for the map editor (D7, D8).

A subsystem gets up to four KEGG maps: those holding at least 3 of its KEGG reactions and at least 10% of
them (overview maps 011xx/012xx excluded). The result is written to subsystem_pathways.json in the cache
folder (or to $KEGG_MAPPING); subsystems already listed there are skipped, so the script can be rerun.

Usage: kegg_fetch.py <model dir with reactions.tsv> <cache dir> <subsystem>...
The subsystems are read from reactions.tsv (column "subsystem", as written by yaml_to_tsv.py) or, for
Human-GEM, from Human-GEM.yml in the model dir.
"""

import collections
import csv
import json
import os
import re
import sys
import time
import urllib.request


def subsystems_of(model_dir, rx):
    if all("subsystem" in r for r in rx.values()):
        return {i: r["subsystem"].split(";")[0] for i, r in rx.items()}
    sub, cur, want = {}, None, False
    for line in open(os.path.join(model_dir, "Human-GEM.yml")):
        m = re.match(r"\s+- id: (\S+)$", line)
        if m:
            cur, want = m.group(1), False
            continue
        if cur and re.match(r"\s+- subsystem:\s*$", line):
            want = True
            continue
        if want:
            sub[cur] = line.strip()[2:]
            want = False
    return sub


def get(url):
    for _ in range(3):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception:
            time.sleep(2)
    return ""


def main():
    model_dir, cache = sys.argv[1], sys.argv[2]
    os.makedirs(cache, exist_ok=True)
    rx = {r["rxns"]: r for r in csv.DictReader(open(os.path.join(model_dir, "reactions.tsv")), delimiter="\t")}
    sub = subsystems_of(model_dir, rx)
    mapping_file = os.environ.get("KEGG_MAPPING") or os.path.join(cache, "subsystem_pathways.json")
    mapping = json.load(open(mapping_file)) if os.path.exists(mapping_file) else {}
    for s in sys.argv[3:]:
        if s in mapping:
            continue
        kegg = sorted({k for r, x in rx.items() if sub.get(r) == s for k in x["rxnKEGGID"].split(";") if k})
        counts = collections.Counter()
        for i in range(0, len(kegg), 50):
            txt = get("https://rest.kegg.jp/link/pathway/" + "+".join("rn:" + k for k in kegg[i:i + 50]))
            for line in txt.splitlines():
                if line.count("\t") != 1:
                    continue
                m = re.match(r"path:(?:rn|map)(\d{5})$", line.split("\t")[1])
                if m and not m.group(1).startswith("01"):
                    counts["rn" + m.group(1)] += 1
        print(s, len(kegg), "KEGG reactions; pathways:", counts.most_common(6))
        mapping[s] = [p for p, n in counts.most_common(4) if n >= max(3, 0.1 * len(kegg))]
        json.dump(mapping, open(mapping_file, "w"), indent=1)
        for p, _ in counts.most_common(4):
            f = os.path.join(cache, p + ".xml")
            if not os.path.exists(f):
                open(f, "w").write(get(f"https://rest.kegg.jp/get/{p}/kgml"))


if __name__ == "__main__":
    main()
