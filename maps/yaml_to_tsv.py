"""Write reactions.tsv, metabolites.tsv and genes.tsv (the Human-GEM column names the map editor reads)
from a model YAML, for models without Human-GEM style TSV files (Yeast-GEM).

Usage: yaml_to_tsv.py <model.yml> <out dir> [<model dir>]

KEGG ids come from the annotations in the YAML; where it has none, from yeast-GEM's annotation tables
(an id column, then columns named after the annotation keys) in <model dir>."""
import collections
import csv
import os
import sys

import yaml

try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader


def as_dict(x):
    return dict(x) if isinstance(x, list) else (x or {})


yml, out = sys.argv[1], sys.argv[2]
model_dir = sys.argv[3] if len(sys.argv) > 3 else None
os.makedirs(out, exist_ok=True)
top = {}
for item in yaml.load(open(yml), Loader=Loader):
    k, v = item if isinstance(item, tuple) else next(iter(item.items()))
    top[k] = v


def ann(x, key):
    v = as_dict(x.get("annotation")).get(key, "")
    return ";".join(v) if isinstance(v, list) else str(v or "")


def table_column(name, column):
    """{id: value} of one column of a yeast-GEM annotation table in model_dir; {} for other tables."""
    path = os.path.join(model_dir, name) if model_dir else ""
    if not os.path.exists(path):
        return {}
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["id"]: r.get(column) or "" for r in rows} if rows and "id" in rows[0] else {}


rxn_kegg = table_column("reactions.tsv", "kegg.reaction")
met_kegg = table_column("metabolites.tsv", "kegg.compound")


degree = collections.Counter()
with open(os.path.join(out, "reactions.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["rxns", "rxnKEGGID", "subsystem"])
    for r in top["reactions"]:
        r = as_dict(r)
        subs = r.get("subsystem") or []
        w.writerow([r["id"], ann(r, "kegg.reaction") or rxn_kegg.get(r["id"], ""), ";".join(subs if isinstance(subs, list) else [subs])])
        for m in as_dict(r.get("metabolites")):
            degree[m] += 1
names = {}
with open(os.path.join(out, "metabolites.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["mets", "metKEGGID", "compartment", "metHMR2ID", "metRetired"])
    for m in top["metabolites"]:
        m = as_dict(m)
        names[m["id"]] = m.get("name", "")
        w.writerow([m["id"], ann(m, "kegg.compound") or met_kegg.get(m["id"], ""), m.get("compartment", ""), "", ""])
with open(os.path.join(out, "genes.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genes", "geneSymbols"])
    for g in top.get("genes", []):
        g = as_dict(g)
        w.writerow([g["id"], g.get("name", "")])
by_name = collections.Counter()
for m, n in degree.items():
    by_name[names.get(m, m)] += n
print("most connected metabolite names:")
for n, c in by_name.most_common(45):
    print(f"  {c:5d} {n}")
