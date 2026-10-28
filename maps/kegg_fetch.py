#!/usr/bin/env python3
"""Find the KEGG reference maps of each subsystem and cache their KGML files, for the map editor (D7, D8).

A subsystem gets up to four KEGG maps: those holding at least 3 of its KEGG reactions and at least 10% of
them (overview maps 011xx/012xx excluded). The result is written to subsystem_pathways.json in the cache
folder (or to $KEGG_MAPPING); subsystems already listed there are skipped, so the script can be rerun.

Usage: kegg_fetch.py <model dir with the YAML file and reactions.tsv> <cache dir> <subsystem>...
The subsystems are read from the model's YAML file (the first subsystem of a reaction) and the KEGG
reaction ids from reactions.tsv, in Human-GEM's or yeast-GEM's column names.
"""

import collections
import glob
import json
import os
import re
import sys
import time
import urllib.request

import mapedit


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
    model = mapedit.Model(glob.glob(os.path.join(model_dir, "*.yml"))[0], model_dir, [])
    kegg_ids = model.rxn_kegg
    sub = {}
    for rid, r in model.rxns.items():
        subsystem = r.get("subsystem")
        sub[rid] = (subsystem[0] if isinstance(subsystem, list) and subsystem else subsystem) or ""
    mapping_file = os.environ.get("KEGG_MAPPING") or os.path.join(cache, "subsystem_pathways.json")
    mapping = json.load(open(mapping_file)) if os.path.exists(mapping_file) else {}
    for s in sys.argv[3:]:
        if s in mapping:
            continue
        kegg = sorted({k for r, ids in kegg_ids.items() if sub.get(r) == s for k in ids.split(";") if k})
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
