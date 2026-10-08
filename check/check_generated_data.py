#!/usr/bin/env python3
"""Check the output of data-generation for one integrated model against its source files.

The model files in integrated-models/<model> are parsed independently of data-generation
and compared with the Neo4j CSV files, the import script, the SVG map links and the data
overlay output. Hard checks fail the run; informational sections describe what changed
compared with a previous run and which identifiers referenced by the MetabolicAtlas tests
are no longer in the model.

Usage (from the data-generation folder):
    python check/check_generated_data.py --model Human-GEM --data-files <data-files checkout> \\
        --new <run dir> [--old <run dir>] [--old-model-dir <dir>] \\
        [--metabolicatlas <repo dir>] [--report report.md]

A run dir is the directory data-generation was run in (it contains neo4j/, dataOverlay/
and gemRepository/).
"""

import argparse
import collections
import csv
import glob
import json
import os
import re
import sys

import yaml

try:
    Loader = yaml.CSafeLoader
except AttributeError:
    Loader = yaml.SafeLoader

ID_PATTERN = re.compile(r"\b(MA[RM]\d{5}[a-z]?|ENSG\d{11})\b")
XREF_HEADER = {"gene": r"gene.*ID$", "reaction": r"rxn.*ID$", "metabolite": r"met.*ID$"}
LFS_POINTER = b"version https://git-lfs"


# ----------------------------------------------------------------------------- parsing

def omap_to_dict(node):
    """Turn !!omap structures (lists of pairs or of one-key dicts) into plain dicts."""
    if isinstance(node, list):
        if node and all(isinstance(x, tuple) and len(x) == 2 for x in node):
            return {k: omap_to_dict(v) for k, v in node}
        if node and all(isinstance(x, dict) and len(x) == 1 for x in node):
            # top level of the model and RAVEN 2 entries written without !!omap
            keys = [next(iter(x)) for x in node]
            if len(set(keys)) == len(keys):
                return {k: omap_to_dict(x[k]) for x, k in zip(node, keys)}
        return [omap_to_dict(x) for x in node]
    if isinstance(node, dict):
        return {k: omap_to_dict(v) for k, v in node.items()}
    return node


def find_yaml(model_dir):
    """The model YAML in a folder; data-generation takes any *.yml or *.yaml file."""
    found = sorted(glob.glob(os.path.join(model_dir, "*.yml")) + glob.glob(os.path.join(model_dir, "*.yaml")))
    if len(found) != 1:
        raise SystemExit(f"expected one YAML file in {model_dir}, found {len(found)}")
    return found[0]


def load_model(model_dir, model):
    path = find_yaml(model_dir)
    with open(path) as f:
        raw = yaml.load(f, Loader=Loader)
    top = omap_to_dict(raw)
    meta = top.get("metaData") or top.get("metadata")
    return {
        "meta": meta,
        "metabolites": {m["id"]: m for m in top["metabolites"]},
        "reactions": {r["id"]: r for r in top["reactions"]},
        "genes": {g["id"]: g for g in top["genes"]},
        "compartments": top["compartments"],
    }


def as_list(value):
    if value is None or value == "":
        return []
    return value if isinstance(value, list) else [value]


def as_text(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        # YAML 1.1 reads unquoted NO/YES as booleans; js-yaml (YAML 1.2) keeps the text
        return "NO" if value is False else "YES"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def format_ec(value):
    codes = []
    for v in as_list(value):
        codes += [c.strip() for c in str(v).split(";") if c.strip()]
    return "; ".join(codes)


def gpr_genes(rule):
    return {t for t in re.split(r"[\s()]+", rule or "") if t and t not in ("and", "or")}


def references_text(reaction):
    """references, or the PubMed ids under annotation/pubmed (yeast-GEM) as "PMID:1;PMID:2"."""
    if reaction.get("references") is not None:
        return as_text(reaction.get("references"))
    ids = as_list((reaction.get("annotation") or {}).get("pubmed"))
    return ";".join(f"PMID:{as_text(i).strip()}" for i in ids)


def pmids(references):
    out = set()
    for ref in str(references or "").split(";"):
        ref = ref.strip()
        if ref.startswith("PMID"):
            ref = re.sub(r"PMID:*", "", ref).strip()
            if ref.isdigit():
                out.add(ref)
    return out


def idfy(name):
    # utils.idfyString in data-generation
    s = re.sub(r"_+", "_", re.sub(r"[^a-z0-9]", "_", name.lower()))
    return re.sub(r"^_|_$", "", s, count=1)


def label_case(short_name):
    # utils.toLabelCase in data-generation: Human-GEM -> HumanGem
    return "".join(p[:1].upper() + p[1:].lower() for p in re.split(r"[^A-Za-z0-9]+", short_name) if p)


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def read_tsv(path):
    """Rows of a model TSV as dicts; header and values with surrounding quotes removed."""
    with open(path, newline="") as f:
        lines = [line.rstrip("\r\n") for line in f if line.strip() and not line.startswith(("#", "@"))]
    if not lines:
        return [], []
    header = [h.strip().strip('"') for h in lines[0].split("\t")]
    rows = []
    for line in lines[1:]:
        values = [v.strip().strip('"') for v in line.split("\t")]
        row = dict(zip(header, values))
        if len(values) != len(header):
            row[None] = len(values)
        rows.append(row)
    return header, rows


def is_keyed(header):
    """yeast-GEM's tables: an id column, then columns named after the annotation keys."""
    return bool(header) and header[0] == "id"


def read_map_tsv(path):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "@")):
                continue
            parts = line.rstrip("\r\n").split("\t")
            rows.append(parts)
    return rows


def read_identifiers(path):
    """Column header -> database name, from crossReferencesDict in identifiers.js."""
    text = open(path).read()
    mapping, duplicates = {}, []
    for block in re.finditer(r"headers:\s*\[([^\]]*)\]\s*,\s*db:\s*\"([^\"]*)\"", text):
        for header in re.findall(r"\"([^\"]+)\"", block.group(1)):
            if header in mapping:
                # data-generation uses the first entry that lists the header
                duplicates.append(f"{header} ({mapping[header]}; ignored: {block.group(2)})")
            else:
                mapping[header] = block.group(2)
    return mapping, duplicates


# ----------------------------------------------------------------------------- report

class Report:
    def __init__(self):
        self.lines = []
        self.failures = []
        self.warnings = []

    def section(self, title):
        self.lines += ["", f"## {title}", ""]

    def text(self, line=""):
        self.lines.append(line)

    def warn(self, text):
        self.lines.append(f"- WARN: {text}")
        self.warnings.append(text)

    def check(self, name, ok, detail=""):
        self.lines.append(f"- {'PASS' if ok else 'FAIL'}: {name}" + (f" ({detail})" if detail else ""))
        if not ok:
            self.failures.append(name)
        return ok

    def compare_sets(self, name, expected, actual, show=5):
        missing = sorted(map(str, expected - actual))
        extra = sorted(map(str, actual - expected))
        detail = f"{len(expected)} expected"
        if missing:
            detail += f"; missing {len(missing)}: {', '.join(missing[:show])}"
        if extra:
            detail += f"; unexpected {len(extra)}: {', '.join(extra[:show])}"
        return self.check(name, not missing and not extra, detail)

    def table(self, header, rows):
        self.lines.append("| " + " | ".join(header) + " |")
        self.lines.append("|" + "---|" * len(header))
        for r in rows:
            self.lines.append("| " + " | ".join(str(x) for x in r) + " |")


# ----------------------------------------------------------------------------- run output

class Run:
    def __init__(self, run_dir, short_name, version=None):
        self.dir = run_dir
        self.neo4j = os.path.join(run_dir, "neo4j")
        label = label_case(short_name)
        if version:
            self.prefix = f"{label}V{version.replace('.', '_')}"
        else:
            found = sorted(glob.glob(os.path.join(self.neo4j, f"{label}V*.reactions.csv")))
            if len(found) != 1:
                raise SystemExit(f"cannot identify {short_name} output in {self.neo4j}: {found}")
            self.prefix = os.path.basename(found[0])[: -len(".reactions.csv")]
        self.version = self.prefix[len(label) + 1:]

    def path(self, kind):
        return os.path.join(self.neo4j, f"{self.prefix}.{kind}.csv")

    def csv(self, kind):
        return read_csv(self.path(kind))

    def ids(self, kind):
        return {r["id"] for r in self.csv(kind)}

    def xrefs(self, component):
        # external database nodes are shared between models, so read every model's file
        nodes = {}
        for f in glob.glob(os.path.join(self.neo4j, "*.externalDbs.csv")):
            nodes.update({r["id"]: (r["dbName"], r["externalId"]) for r in read_csv(f)})
        key = {"gene": "geneId", "reaction": "reactionId",
               "metabolite": "compartmentalizedMetaboliteId"}[component]
        kind = {"gene": "geneExternalDbs", "reaction": "reactionExternalDbs",
                "metabolite": "compartmentalizedMetaboliteExternalDbs"}[component]
        return {(r[key],) + nodes[r["externalDbId"]] for r in self.csv(kind)}


# ----------------------------------------------------------------------------- checks

def check_metadata(rep, model, m, run, data_files):
    rep.section("Metadata and import script")
    meta = m["meta"] or {}
    rep.check("metaData.short_name equals the folder name", meta.get("short_name") == model,
              f"short_name={meta.get('short_name')!r}")
    version = meta.get("version")
    rep.check("metaData.version is a string", isinstance(version, str), f"version={version!r}")
    integrated = json.load(open(os.path.join(data_files, "integrated-models", "integratedModels.json")))
    entry = next((e for e in integrated if e["short_name"] == model), None)
    rep.check("integratedModels.json has an entry", entry is not None)
    if entry:
        rep.check("integratedModels.json version matches the YAML", entry["version"] == version,
                  f"{entry['version']} vs {version}")
        rep.check("integratedModels.json date matches the YAML", entry["date"] == str(meta.get("date")),
                  f"{entry['date']} vs {meta.get('date')}")
    rep.check("output prefix matches the YAML version", run.version == str(version).replace(".", "_"),
              run.prefix)
    cypher = open(os.path.join(run.neo4j, "import.cypher")).read()
    label = label_case(model)
    prefixes = set(re.findall(rf"file:///({label}V[0-9_]+)\.", cypher))
    rep.check("import.cypher loads only this version of the model", prefixes == {run.prefix},
              ", ".join(sorted(prefixes)))
    blocks = [b for b in cypher.split(";") if f"{run.prefix}." in b and "]->" in b]
    rels = set()
    for b in blocks:
        rels |= set(re.findall(r"\[:(V[0-9_]+)\]", b))
    rep.check("state relationships are typed with the model version", rels == {f"V{run.version}"},
              ", ".join(sorted(rels)))
    for f in sorted(glob.glob(os.path.join(run.neo4j, f"{run.prefix}.*.csv"))):
        if f"file:///{os.path.basename(f)}" not in cypher:
            rows = sum(1 for line in open(f) if line.strip()) - 1
            rep.check(f"{os.path.basename(f)} is loaded by import.cypher", rows == 0,
                      f"not loaded, {rows} rows")


def check_components(rep, m, run):
    rep.section("Components")
    mets, rxns, genes = m["metabolites"], m["reactions"], m["genes"]
    rep.compare_sets("genes", set(genes), run.ids("genes"))
    rep.compare_sets("reactions", set(rxns), run.ids("reactions"))
    rep.compare_sets("compartmentalized metabolites", set(mets), run.ids("compartmentalizedMetabolites"))

    comp_states = run.csv("compartmentStates")
    rep.compare_sets("compartments (letter, name)",
                     {(k, v) for k, v in m["compartments"].items()},
                     {(r["letterCode"], r["name"]) for r in comp_states})
    comp_id = {r["letterCode"]: r["compartmentId"] for r in comp_states}
    rep.compare_sets("metabolite compartments",
                     {(i, comp_id.get(x.get("compartment"))) for i, x in mets.items()},
                     {(r["compartmentalizedMetaboliteId"], r["compartmentId"])
                      for r in run.csv("compartmentalizedMetaboliteCompartments")})

    met_of = {r["compartmentalizedMetaboliteId"]: r["metaboliteId"]
              for r in run.csv("compartmentalizedMetaboliteMetabolites")}
    states = {r["metaboliteId"]: r for r in run.csv("metaboliteStates")}
    # metabolites are one node per name; compartments of a metabolite share one state
    variants = collections.defaultdict(set)
    for i, x in mets.items():
        variants[met_of.get(i)].add((as_text(x.get("formula")), as_text(x.get("charge"))))
    bad, inconsistent = [], set()
    for i, x in mets.items():
        s = states.get(met_of.get(i))
        expected = (as_text(x.get("name")), as_text(x.get("formula")), as_text(x.get("charge")))
        if s is None or (s["name"], s["formula"], s["charge"]) != expected:
            if s is not None and s["name"] == expected[0] and len(variants[met_of[i]]) > 1:
                inconsistent.add(met_of[i])
            else:
                bad.append(i)
    rep.check("metabolite name, formula and charge", not bad,
              f"{len(bad)} differ: {', '.join(bad[:5])}" if bad else f"{len(mets)} checked")
    for met in sorted(inconsistent):
        ids = sorted(i for i in mets if met_of.get(i) == met)
        detail = ", ".join(f"{i} {as_text(mets[i].get('formula'))}/{as_text(mets[i].get('charge'))}" for i in ids)
        rep.warn(f"model: compartments of '{states[met]['name']}' differ in formula/charge, "
                 f"only {states[met]['formula']}/{states[met]['charge']} is shown ({detail})")

    subsystems = {s for r in rxns.values() for s in as_list(r.get("subsystem"))}
    by_id = collections.defaultdict(set)
    for s in subsystems:
        by_id[idfy(s)].add(s)
    sub_states = run.csv("subsystemStates")
    shown = {r["subsystemId"]: r["name"] for r in sub_states}
    for sid, names in sorted(by_id.items()):
        if len(names) > 1:
            users = sorted(i for i, r in rxns.items()
                           if set(as_list(r.get("subsystem"))) & (names - {shown.get(sid)}))
            rep.warn(f"model: subsystems {sorted(names)} get the same id '{sid}', shown as "
                     f"'{shown.get(sid)}' (other spelling used by {', '.join(users[:10])})")
    rep.compare_sets("subsystems", set(by_id), set(shown))
    sub_id = {s: idfy(s) for s in subsystems}
    rep.compare_sets("reaction subsystems",
                     {(i, sub_id.get(s)) for i, r in rxns.items() for s in as_list(r.get("subsystem"))},
                     {(r["reactionId"], r["subsystemId"]) for r in run.csv("reactionSubsystems")})


def check_gene_names(rep, model_dir, m, run):
    """Gene symbol, name and aliases from a Human-GEM style genes.tsv, else the name in the YAML."""
    header, rows = read_tsv(os.path.join(model_dir, "genes.tsv"))
    table = {} if not header or is_keyed(header) else {r[header[0]]: r for r in rows}
    states = {r["geneId"]: r for r in run.csv("geneStates")}
    bad = []
    for g, gene in m["genes"].items():
        if g in table:
            expected = (table[g].get("geneSymbols", ""), table[g].get("geneNames", ""),
                        table[g].get("geneAliases", ""))
        else:
            expected = (as_text(gene.get("name")), "", "")
        s = states.get(g)
        if s is None or (s["name"], s["alternateName"], s["synonyms"]) != expected:
            bad.append(g)
    named = sum(1 for s in states.values() if s["name"])
    rep.check("gene symbol, name and aliases", not bad,
              f"{len(bad)} differ: {', '.join(sorted(bad)[:5])}" if bad
              else f"{len(states)} checked, {named} with a name")


def check_reactions(rep, model_dir, m, run):
    rep.section("Reactions")
    rxns = m["reactions"]
    header, rows = read_tsv(os.path.join(model_dir, "reactions.tsv"))
    # yeast-GEM's reactions.tsv has the EC codes that its YAML no longer has
    table_ec = {r["id"]: r.get("ec-code", "") for r in rows} if is_keyed(header) else {}
    states = {r["reactionId"]: r for r in run.csv("reactionStates")}
    fields = collections.Counter()
    for i, r in rxns.items():
        s = states.get(i)
        if s is None:
            fields["missing"] += 1
            continue
        expected = {
            "name": as_text(r.get("name")),
            "lowerBound": as_text(r.get("lower_bound")),
            "upperBound": as_text(r.get("upper_bound")),
            "geneRule": as_text(r.get("gene_reaction_rule")),
            "reversible": "true" if r.get("lower_bound") == -1000 else "false",
            # RAVEN 3 and raven-toolbox write EC codes as annotation/ec-code
            "ec": format_ec(r.get("eccodes") or (r.get("annotation") or {}).get("ec-code"))
            or format_ec(table_ec.get(i)),
            # yeast-GEM writes PubMed ids as annotation/pubmed
            "references": references_text(r),
        }
        for k, v in expected.items():
            if s[k] != v:
                fields[k] += 1
    rep.check("reaction name, bounds, GPR, reversibility, EC and references", not fields,
              ", ".join(f"{k}: {n}" for k, n in fields.items()) or f"{len(rxns)} checked")
    rep.check("EC codes are separated by '; '", all("," not in s["ec"] for s in states.values()))

    expected = set()
    for i, r in rxns.items():
        for met, coef in (r.get("metabolites") or {}).items():
            # PyYAML (YAML 1.1) reads exponents without a decimal point, e.g. -1e-06, as text
            coef = float(coef)
            expected.add((i, met, "in" if coef < 0 else "out", abs(coef)))
    actual = {(r["reactionId"], r["compartmentalizedMetaboliteId"], "in", float(r["stoichiometry"]))
              for r in run.csv("compartmentalizedMetaboliteReactions")}
    actual |= {(r["reactionId"], r["compartmentalizedMetaboliteId"], "out", float(r["stoichiometry"]))
               for r in run.csv("reactionCompartmentalizedMetabolites")}
    rep.compare_sets("stoichiometry (substrates and products)", expected, actual)

    rep.compare_sets("reaction-gene links from the GPRs",
                     {(i, g) for i, r in rxns.items() for g in gpr_genes(r.get("gene_reaction_rule"))},
                     {(r["reactionId"], r["geneId"]) for r in run.csv("reactionGenes")})

    expected = {(i, p) for i, r in rxns.items() for p in pmids(references_text(r))}
    rep.compare_sets("reaction PubMed links", expected,
                     {(r["reactionId"], r["pubmedReferenceId"]) for r in run.csv("reactionPubmedReferences")})
    # PubMed nodes are shared between models, so read every model's file
    nodes = set()
    for f in glob.glob(os.path.join(run.neo4j, "*.pubmedReferences.csv")):
        nodes |= {r["id"] for r in read_csv(f)}
    missing = {p for _, p in expected} - nodes
    rep.check("PubMed nodes exist for every linked PMID", not missing,
              f"missing {len(missing)}: {', '.join(sorted(missing)[:5])}" if missing else "")


def check_xrefs(rep, model_dir, m, run, identifiers):
    rep.section("Cross-references from the TSV files")
    ids = {"gene": set(m["genes"]), "reaction": set(m["reactions"]), "metabolite": set(m["metabolites"])}
    for component, filename in (("gene", "genes.tsv"), ("reaction", "reactions.tsv"),
                                ("metabolite", "metabolites.tsv")):
        header, rows = read_tsv(os.path.join(model_dir, filename))
        if not header:
            rep.compare_sets(f"{component} cross-reference links (no {filename} for this model)", set(),
                             run.xrefs(component))
            continue
        # data-generation reads columns by position, so a short or long row shifts its values
        ragged = [r[header[0]] for r in rows if len(r) != len(header) or None in r]
        rep.check(f"{filename}: every row has {len(header)} columns", not ragged,
                  f"{len(ragged)} rows differ: {', '.join(ragged[:5])}" if ragged else "")
        keyed = is_keyed(header)
        if keyed:
            # columns named after annotation keys: those identifiers.js knows are cross-references
            columns = [h for h in header[1:] if h in identifiers]
        else:
            columns = [h for h in header[1:] if re.search(XREF_HEADER[component], h)]
            unknown = [h for h in columns if h not in identifiers]
            rep.check(f"{filename}: every cross-reference column is known to identifiers.js", not unknown,
                      ", ".join(unknown))
        skipped = [h for h in header[1:] if h not in columns]
        if skipped:
            rep.text(f"  - {filename}: columns not imported as cross-references: {', '.join(skipped)}")
        if component == "gene" and not keyed:
            # Human-GEM's gene ids are Ensembl and Protein Atlas ids
            columns += ["geneEnsemblID", "geneProteinAtlasID"]
        expected = set()
        for row in rows:
            i = row[header[0]]
            if i not in ids[component]:
                continue
            if component == "gene" and not keyed:
                row = dict(row, geneEnsemblID=i, geneProteinAtlasID=i)
            for h in columns:
                db = identifiers.get(h)
                value = row.get(h, "")
                if db is None or not value:
                    continue
                if db == "ChEBI":
                    value = re.sub(r"^CHEBI:", "", value)
                elif db in ("Rhea", "RheaMaster"):
                    value = re.sub(r"^RHEA:", "", value)
                for v in (x.strip() for x in value.split(";")):
                    if v:
                        expected.add((i, db, "CHEBI:" + v if db == "ChEBI" else v))
        rep.compare_sets(f"{component} cross-reference links", expected, run.xrefs(component))
        if component == "gene":
            header_line = open(os.path.join(model_dir, filename)).readline()
            if '"' in header_line:
                rep.text("  - genes.tsv has a quoted header line; data-generation must unquote it")


def check_maps(rep, model_dir, data_files, model, m, run):
    rep.section("SVG maps")
    svg_dir = os.path.join(data_files, "svg", model)
    svg_nodes = {r["id"]: r for r in run.csv("svgMaps")}
    for component, kind in (("compartment", "compartmentSvgMaps"), ("subsystem", "subsystemSvgMaps"),
                            ("custom", "customSvgMaps")):
        rows = read_map_tsv(os.path.join(model_dir, f"{component}SVG.tsv"))
        if component == "custom":
            # custom maps (map name, filename) stand alone: a map node, linked to no component
            rep.compare_sets("custom maps", {r[1] for r in rows if len(r) >= 2},
                             {n["filename"] for n in svg_nodes.values()} & {r[1] for r in rows if len(r) >= 2})
            continue
        expected = {(idfy(r[0]), r[2][:-4]) for r in rows if len(r) >= 3}
        key = f"{component}Id"
        rep.compare_sets(f"{component} map links", expected,
                         {(r[key], r["svgMapId"]) for r in run.csv(kind)})
    missing, pointers = [], []
    for node in svg_nodes.values():
        path = os.path.join(svg_dir, node["filename"])
        if not os.path.exists(path):
            missing.append(node["filename"])
        else:
            with open(path, "rb") as f:
                if f.read(len(LFS_POINTER)) == LFS_POINTER:
                    pointers.append(node["filename"])
    rep.check("every linked SVG file exists", not missing, ", ".join(missing[:5]))
    rep.check("SVG files are real files, not Git LFS pointers (run git lfs pull)", not pointers,
              f"{len(pointers)} pointers" if pointers else "")

    rep.section("Reactions drawn on the maps (informational)")
    rxns = set(m["reactions"])
    rows = []
    for node in sorted(svg_nodes.values(), key=lambda n: n["filename"]):
        path = os.path.join(svg_dir, node["filename"])
        if not os.path.exists(path) or node["filename"] in pointers:
            continue
        text = open(path, encoding="utf8", errors="replace").read()
        drawn = set(re.findall(r'class="rea"[^>]*?\bid="(MAR\d{5}|r_\d{4})"', text)) | set(
            re.findall(r'\bid="(MAR\d{5}|r_\d{4})"[^>]*?class="rea"', text))
        absent = sorted(drawn - rxns)
        if absent:
            rows.append((node["filename"], len(drawn), len(absent), " ".join(absent[:8])))
    if rows:
        rep.text("Maps with reactions that are not in the model:")
        rep.text()
        rep.table(["map", "reactions drawn", "not in model", "examples"], rows)
    else:
        rep.text("All reactions drawn on the maps are in the model.")


def check_overlays(rep, model_dir, model, m, run, old_run):
    rep.section("Data overlay")
    src = os.path.join(model_dir, "dataOverlay")
    out = os.path.join(run.dir, "dataOverlay", model)
    if not os.path.isdir(src):
        rep.text("No data overlay for this model.")
        return
    rep.check("dataOverlay index.json was generated", os.path.exists(os.path.join(out, "index.json")))
    ids = {"gene": set(m["genes"]), "reaction": set(m["reactions"]),
           "metabolite": set(m["metabolites"])}
    rows = []
    for path in sorted(glob.glob(os.path.join(src, "*", "*.tsv"))):
        component, name = path.split(os.sep)[-2:]
        if name == "index.tsv":
            continue
        generated = os.path.join(out, component, name)
        rep.check(f"{component}/{name} was generated", os.path.exists(generated))
        source_ids = {line.split("\t", 1)[0] for line in list(open(path))[1:]}
        in_model = source_ids & ids.get(component, set())
        new_count = sum(1 for _ in open(generated)) - 1 if os.path.exists(generated) else 0
        old_path = os.path.join(old_run.dir, "dataOverlay", model, component, name) if old_run else None
        old_count = sum(1 for _ in open(old_path)) - 1 if old_path and os.path.exists(old_path) else ""
        no_data = len(ids.get(component, set()) - source_ids)
        rows.append((f"{component}/{name}", len(source_ids), len(in_model), old_count, new_count, no_data))
    rep.text()
    rep.table(["file", "ids in file", "ids in model", "rows generated (old)", "rows generated (new)",
               "model ids without data"], rows)


def compare_runs(rep, m, run, old_run, old_model):
    rep.section(f"Changes compared with {old_run.prefix} (informational)")
    rows = []
    for f in sorted(glob.glob(os.path.join(run.neo4j, f"{run.prefix}.*.csv"))):
        kind = os.path.basename(f)[len(run.prefix) + 1: -4]
        old_path = old_run.path(kind)
        old_n = sum(1 for _ in open(old_path)) - 1 if os.path.exists(old_path) else "-"
        new_n = sum(1 for _ in open(f)) - 1
        rows.append((kind, old_n, new_n, new_n - old_n if old_n != "-" else ""))
    rep.table(["file", old_run.prefix, run.prefix, "change"], rows)

    rep.text()
    rep.text("Cross-reference links per database:")
    rep.text()
    rows = []
    for component in ("gene", "reaction", "metabolite"):
        old = collections.Counter(x[1] for x in old_run.xrefs(component))
        new = collections.Counter(x[1] for x in run.xrefs(component))
        for db in sorted(set(old) | set(new)):
            change = new[db] - old[db]
            flag = " (gone)" if not new[db] else (" (dropped >5%)" if new[db] < 0.95 * old[db] else "")
            rows.append((component, db, old[db], new[db], f"{change:+d}{flag}"))
    rep.table(["component", "database", "old", "new", "change"], rows)

    if old_model:
        rep.text()
        rep.text("Identifiers added and removed:")
        rep.text()
        rows = []
        for kind in ("metabolites", "reactions", "genes"):
            a, b = set(old_model[kind]), set(m[kind])
            added, removed = sorted(b - a), sorted(a - b)
            rows.append((kind, len(removed), " ".join(removed), len(added), " ".join(added)))
        rep.table(["type", "removed", "ids", "added", "ids"], rows)


def check_test_references(rep, metabolicatlas, m):
    rep.section("Identifiers used in MetabolicAtlas tests and frontend (informational)")
    all_ids = set(m["metabolites"]) | set(m["reactions"]) | set(m["genes"])
    all_ids |= {i[:-1] for i in m["metabolites"]}  # metabolite ids without compartment
    rows = []
    for sub in ("api/test", "frontend/cypress/e2e", "frontend/src"):
        for path in sorted(glob.glob(os.path.join(metabolicatlas, sub, "**", "*"), recursive=True)):
            if not path.endswith((".js", ".vue", ".json")) or not os.path.isfile(path):
                continue
            for n, line in enumerate(open(path, encoding="utf8", errors="replace"), 1):
                for i in ID_PATTERN.findall(line):
                    if i not in all_ids:
                        rows.append((os.path.relpath(path, metabolicatlas) + f":{n}", i))
    if rows:
        rep.text("Not in the model (some may refer to other models; Cypress fixtures are stubbed and not listed):")
        rep.text()
        rep.table(["file:line", "id"], rows)
    else:
        rep.text("Every MAR, MAM and ENSG identifier referenced there exists in the model.")


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="folder name in integrated-models, e.g. Human-GEM")
    ap.add_argument("--data-files", required=True, help="data-files checkout used for --new")
    ap.add_argument("--data-generation", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    help="data-generation checkout for identifiers.js (default: the one this script is in)")
    ap.add_argument("--new", required=True, help="run dir to check")
    ap.add_argument("--old", help="earlier run dir to compare with")
    ap.add_argument("--old-model-dir", help="integrated-models/<model> folder of the earlier version")
    ap.add_argument("--metabolicatlas", help="MetabolicAtlas checkout, to list stale test identifiers")
    ap.add_argument("--report", help="write the report to this Markdown file")
    args = ap.parse_args()

    model_dir = os.path.join(args.data_files, "integrated-models", args.model)
    m = load_model(model_dir, args.model)
    version = str((m["meta"] or {}).get("version"))
    run = Run(args.new, args.model, version)
    old_run = Run(args.old, args.model) if args.old else None
    old_model = load_model(args.old_model_dir, args.model) if args.old_model_dir else None
    identifiers, duplicates = read_identifiers(os.path.join(args.data_generation, "identifiers.js"))

    rep = Report()
    rep.text(f"# Generated data check: {args.model} {version}")
    rep.text()
    rep.text(f"Run: `{os.path.abspath(args.new)}` ({run.prefix})")
    provenance = os.path.join(args.new, "provenance.txt")
    if os.path.exists(provenance):
        for line in open(provenance):
            rep.text(f"- {line.strip()}")

    check_metadata(rep, args.model, m, run, args.data_files)
    check_components(rep, m, run)
    check_gene_names(rep, model_dir, m, run)
    check_reactions(rep, model_dir, m, run)
    check_xrefs(rep, model_dir, m, run, identifiers)
    if duplicates:
        rep.warn("identifiers.js lists these headers more than once: " + "; ".join(duplicates))
    check_maps(rep, model_dir, args.data_files, args.model, m, run)
    check_overlays(rep, model_dir, args.model, m, run, old_run)
    if old_run:
        compare_runs(rep, m, run, old_run, old_model)
    if args.metabolicatlas:
        check_test_references(rep, args.metabolicatlas, m)

    summary = (f"**Result: {len(rep.failures)} hard check(s) failed**" if rep.failures
               else "**Result: all hard checks passed**") + f", {len(rep.warnings)} warning(s)"
    rep.lines.insert(1, "")
    rep.lines.insert(2, summary)
    text = "\n".join(rep.lines) + "\n"
    if args.report:
        with open(args.report, "w") as f:
            f.write(text)
    print(text if not args.report else summary)
    for name in rep.failures:
        print(f"FAIL: {name}", file=sys.stderr)
    sys.exit(1 if rep.failures else 0)


if __name__ == "__main__":
    main()
