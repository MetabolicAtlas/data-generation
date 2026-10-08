#!/usr/bin/env python3
"""Fit the Metabolic Atlas SVG maps of a model (Human-GEM, Yeast-GEM) to a model release.

The maps are edited in place, keeping the drawing: ids follow the model, reactions and metabolites
that are gone are removed, swapped metabolites are relabelled, gene boxes follow the gene rules,
missing reactions of a map's subsystem or compartment are added and placed, parts that share a
metabolite are connected, and the background band follows the edges. The rules (R1-R10 for
correctness, D1-D10 for layout) are listed in RULES.md next to this file; changes.tsv in the output
folder lists every change by rule, and --review also writes copies with the changes highlighted.

Usage: mapedit.py <model.yml> <model dir with reactions/metabolites/genes.tsv> [<earlier model.yml>...]
       --maps a.svg ... --out <dir> [--review] --map-table subsystemSVG.tsv --compartment-table
       compartmentSVG.tsv --custom-table customSVG.tsv [--one-area] [--gene-label name|orf|both]
       [--kegg-dir <KGML cache> [--kegg-mapping <json>] [--kegg-relayout <new map>...]] [--parts <json>]
"""

import argparse
import collections
import copy
import csv
import math
import os
import re
import sys

import yaml
from lxml import etree

Loader = yaml.CSafeLoader
SVG = "{http://www.w3.org/2000/svg}"
CURRENCY = {"H+", "H2O", "ATP", "ADP", "AMP", "Pi", "PPi", "NAD+", "NADH", "NADP+", "NADPH", "CO2", "O2", "CoA",
            "FAD", "FADH2", "H2O2", "GTP", "GDP", "UTP", "UDP", "CTP", "CDP", "NH3", "HCO3-", "Na+", "K+", "Cl-",
            "ubiquinone", "ubiquinol",
            # Yeast-GEM names
            "phosphate", "diphosphate", "coenzyme A", "NAD", "NADP(+)", "carbon dioxide", "oxygen", "ammonium",
            "hydrogen peroxide", "bicarbonate", "potassium", "sodium", "chloride", "CMP", "UMP", "GMP",
            "ubiquinone-6", "ubiquinol-6", "ferricytochrome c", "ferrocytochrome c"}
if __name__ == "__main__":  # one module for this file, also when run as a script (layout.py imports it)
    # the same maps on every run: placement goes through sets of ids, whose order depends on the hash seed
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable] + sys.argv)
    sys.modules.setdefault("mapedit", sys.modules[__name__])

SUFFIX = {"p": "x", "s": "e"}
RXN_PREFIXES = ("MAR", "r_")  # reaction ids: Human-GEM, Yeast-GEM
MET_COMP = {}  # metabolite id -> compartment code, from the model


def is_rxn(c):
    return c.startswith(RXN_PREFIXES)


def comp_of(m):
    """Compartment code of a metabolite: from the model, else the last letter of a Human-GEM id."""
    return MET_COMP.get(m) or m[-1]
HIGHLIGHT = {"removed": "#e5484d", "relabelled": "#f59f00", "detached": "#e5484d", "id fixed": "#12b886"}


# ----------------------------------------------------------------------------- model

def text(v):
    return "NO" if v is False else ("" if v is None else str(v))


def load_model(path):
    top = {}
    for item in yaml.load(open(path), Loader=Loader):
        k, v = item if isinstance(item, tuple) else next(iter(item.items()))
        top[k] = v
    td = lambda x: dict(x) if isinstance(x, list) else x  # noqa: E731
    mets = {}
    for m in top["metabolites"]:
        m = td(m)
        m["name"] = text(m.get("name"))
        mets[m["id"]] = m
    rxns = {}
    for r in top["reactions"]:
        r = td(r)
        r["stoich"] = {m: float(c) for m, c in td(r.get("metabolites") or {}).items()}
        rule = text(r.get("gene_reaction_rule"))
        r["genes"] = {t for t in re.split(r"[\s()]+", rule) if t and t not in ("and", "or")}
        rxns[r["id"]] = r
    meta = td(top.get("metaData") or {})
    load_model.compartments = {k: text(v) for k, v in td(top.get("compartments") or {}).items()}
    return mets, rxns, text(meta.get("version"))


def read_tsv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = [h.strip('"') for h in rows[0]]
    return [dict(zip(header, [v.strip('"') for v in r])) for r in rows[1:] if r]


class Model:
    def __init__(self, yml, model_dir, old_ymls):
        self.mets, self.rxns, self.version = load_model(yml)
        import layout
        for k, v in load_model.compartments.items():  # the model's names count; Human-GEM's are the defaults
            layout.COMPARTMENT[k] = v[:1].upper() + v[1:]
        MET_COMP.update({m: x["compartment"] for m, x in self.mets.items() if x.get("compartment")})
        self.old_mets, self.old_rxns, self.old_versions = {}, {}, []
        for p in old_ymls:
            om, orx, _ = load_model(p)
            self.old_mets.update(om)
            self.old_rxns.update(orx)
            self.old_versions.append(orx)
        self.symbol = {g["genes"]: g.get("geneSymbols", "") for g in read_tsv(os.path.join(model_dir, "genes.tsv"))}
        self.rxn_kegg = {r["rxns"]: r.get("rxnKEGGID", "") for r in read_tsv(os.path.join(model_dir, "reactions.tsv"))}
        self.met_kegg = {m["mets"]: m.get("metKEGGID", "") for m in read_tsv(os.path.join(model_dir, "metabolites.tsv"))}
        # legacy metabolite ids -> current ids
        self.legacy = {}
        for row in read_tsv(os.path.join(model_dir, "metabolites.tsv")):
            mid = row["mets"]
            comp = mid[-1]
            old_comp = {v: k for k, v in SUFFIX.items()}.get(comp, comp)
            for col in ("metHMR2ID", "metRetired"):
                for v in row.get(col, "").split(";"):
                    if v:
                        self.legacy.setdefault(v, mid)
            for col in ("metBiGGID", "metRecon3DID"):
                for v in row.get(col, "").split(";"):
                    if v:
                        for c in {comp, old_comp}:
                            self.legacy.setdefault(f"{v}_{c}", mid)
                            self.legacy.setdefault(f"{v}[{c}]", mid)
        self.by_stoich = collections.defaultdict(list)
        self.by_core = collections.defaultdict(list)
        for i, r in self.rxns.items():
            self.by_stoich[self._key(r["stoich"])].append(i)
            self.by_core[self._core(r["stoich"])].append(i)

    def name(self, m):
        return text((self.mets.get(m) or self.old_mets.get(m) or {}).get("name", m))

    @staticmethod
    def _key(st):
        a = tuple(sorted((m, round(c, 6)) for m, c in st.items()))
        b = tuple(sorted((m, round(-c, 6)) for m, c in st.items()))
        return min(a, b)

    def _core(self, st):
        a = tuple(sorted((self.name(m), comp_of(m), c > 0) for m, c in st.items() if self.name(m) not in CURRENCY))
        b = tuple(sorted((n, comp, not d) for n, comp, d in a))
        return min(a, b)

    def replacement(self, rid):
        """(new id, 'identical'|'same main metabolites') for a reaction no longer in the model, or None."""
        old = self.old_rxns.get(rid)
        if not old:
            return None
        st = {self.fix_met_id(m, set()) or m: c for m, c in old["stoich"].items()}
        same = [x for x in self.by_stoich.get(self._key(st), []) if x != rid]
        if len(same) == 1:
            return same[0], "identical"
        core = [x for x in self.by_core.get(self._core(st), []) if x != rid]
        if len(core) == 1:
            return core[0], "same main metabolites"
        return None

    def fix_met_id(self, mid, reaction_ids):
        """Current id for a metabolite id on a map, or None if it cannot be resolved."""
        if mid in self.mets:
            return mid
        if mid.startswith("MAM") and mid[-1] in SUFFIX and mid[:-1] + SUFFIX[mid[-1]] in self.mets:
            return mid[:-1] + SUFFIX[mid[-1]]
        if mid in self.legacy:
            return self.legacy[mid]
        # merged metabolite: same name and compartment among the reaction's metabolites
        name = self.name(mid).lower()
        comp = self.old_mets.get(mid, {}).get("compartment") or SUFFIX.get(mid[-1], mid[-1])
        hits = {m for r in reaction_ids if r in self.rxns for m in self.rxns[r]["stoich"]
                if self.mets[m]["name"].lower() == name and comp_of(m) == comp}
        return hits.pop() if len(hits) == 1 else None


# ----------------------------------------------------------------------------- geometry

TOKEN = re.compile(r"[MmLlHhVvZzCcSsQqTtAa]|-?\d*\.?\d+(?:e-?\d+)?")


def subpaths(d):
    """Split a path into sub-paths of absolute points; None if it uses curves or arcs."""
    out, x, y = [], 0.0, 0.0
    for piece in [p for p in re.split(r"(?=[Mm])", d) if p.strip()]:
        toks, pts, cmd, i = TOKEN.findall(piece), [], None, 0
        while i < len(toks):
            t = toks[i]
            if t.isalpha():
                cmd, i = t, i + 1
                if cmd in "Zz" and pts:
                    x, y = pts[0]
                    pts.append(pts[0])
                elif cmd in "CcSsQqTtAa":
                    return None
                continue
            a = float(t)
            if cmd in "MmLl":
                b = float(toks[i + 1])
                x, y = (x + a, y + b) if cmd in "ml" else (a, b)
                pts.append((x, y))
                cmd = {"m": "l", "M": "L"}.get(cmd, cmd)
                i += 2
            elif cmd in "Hh":
                x = x + a if cmd == "h" else a
                pts.append((x, y))
                i += 1
            elif cmd in "Vv":
                y = y + a if cmd == "v" else a
                pts.append((x, y))
                i += 1
            else:
                return None
        out.append(pts)
    return out


def to_d(paths):
    parts = []
    for pts in paths:
        if not pts:
            continue
        s = "M%.1f %.1f" % pts[0]
        s += "".join("L%.1f %.1f" % p for p in pts[1:])
        parts.append(s)
    return "".join(parts)


def seg_dist(p, a, b):
    """Distance from point p to segment a-b."""
    ax, ay, bx, by = a[0], a[1], b[0], b[1]
    dx, dy = bx - ax, by - ay
    if dx == dy == 0:
        return math.hypot(p[0] - ax, p[1] - ay)
    t = max(0, min(1, ((p[0] - ax) * dx + (p[1] - ay) * dy) / (dx * dx + dy * dy)))
    return math.hypot(p[0] - ax - t * dx, p[1] - ay - t * dy)


def segments(paths):
    return [(pts[i], pts[i + 1]) for pts in paths for i in range(len(pts) - 1) if pts[i] != pts[i + 1]]


def near(p, q, r):
    return math.hypot(p[0] - q[0], p[1] - q[1]) <= r


def translate(el):
    t = el.get("transform", "")
    m = re.search(r"matrix\(([^)]*)\)", t)
    if m:
        a = [float(v) for v in m.group(1).split(",")]
        return a[4], a[5]
    for child in el:
        m = re.search(r"matrix\(([^)]*)\)", child.get("transform", ""))
        if m:
            a = [float(v) for v in m.group(1).split(",")]
            return a[4], a[5]
    return None


# ----------------------------------------------------------------------------- map

class Map:
    def __init__(self, path):
        self.path = path
        self.tree = etree.parse(path, etree.XMLParser(remove_blank_text=False))
        self.root = self.tree.getroot()
        self.layer = {g.get("id"): g for g in self.root.iter(SVG + "g") if g.get("id") in ("fes", "fehs", "ees", "nodes")}
        # maps without genes (or without edges) lack some layers; create them empty, in drawing order
        style = {"fes": {"fill": "none", "stroke": "#666", "stroke-width": "3"}, "fehs": {"fill": "#666"},
                 "ees": {"fill": "none", "stroke": "#666", "stroke-width": "3"}}
        for name in ("fes", "fehs", "ees"):
            if name not in self.layer:
                g = etree.Element(SVG + "g", id=name, **style[name])
                g.text, g.tail = "\n  ", "\n  "
                self.layer["nodes"].addprevious(g)
                self.layer[name] = g

    def classes(self, el):
        return (el.get("class") or "").split()

    def reactions(self):
        return {g.get("id"): g for g in self.layer["nodes"] if g.get("class") == "rea"}

    def nodes(self, kind):
        return [g for g in self.layer["nodes"] if self.classes(g)[:1] == [kind]]

    def edges(self, layer, rid):
        return [p for p in self.layer[layer] if rid in self.classes(p)[1:]]

    @staticmethod
    def remove(el):
        parent = el.getparent()
        prev = el.getprevious()
        if el.tail and prev is not None:
            prev.tail = (prev.tail or "")  # keep the previous element's indentation
        elif el.tail and prev is None:
            parent.text = (parent.text or "")
        parent.remove(el)

    def write(self, path):
        self.tree.write(path, xml_declaration=False, encoding="utf-8")


def node_radius(el):
    shape = el.find(SVG + "ellipse")
    return 45 if shape is not None else 22


def trim_band(mp, log, reach=25):
    """R7: cut the background bands back wherever they no longer run along an edge."""
    import layout
    kept_edges = segments([sp for p in mp.layer["fes"] for sp in (subpaths(p.get("d")) or [])])
    grid = collections.defaultdict(list)
    for a, b in kept_edges:
        for i in range(int(min(a[0], b[0]) // 200) - 1, int(max(a[0], b[0]) // 200) + 2):
            for j in range(int(min(a[1], b[1]) // 200) - 1, int(max(a[1], b[1]) // 200) + 2):
                grid[(i, j)].append((a, b))

    def covered(mid):
        return any(seg_dist(mid, a, b) <= reach for a, b in grid.get((int(mid[0] // 200), int(mid[1] // 200)), ()))

    total = 0
    for g, band in layout.band_groups(mp):
        sps = subpaths(band.get("d") or "")
        if sps is None:
            log.append(("review", "band not trimmed (curves)", g.get("id"), "", ""))
            continue
        out, dropped = [], 0
        for pts in sps:
            runs, run, cut_here = [], [pts[0]] if pts else [], False
            for a, b in zip(pts, pts[1:]):
                mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
                if a != b and not covered(mid):
                    dropped += 1
                    cut_here = True
                    if len(run) > 1:
                        runs.append(run)
                    run = [b]
                else:
                    run.append(b)
            if len(run) > 1:
                runs.append(run)
            if cut_here:
                # fragments left by the cut that are too short to read as band (they render as dots)
                runs = [r for r in runs if sum(math.hypot(q[0] - p[0], q[1] - p[1]) for p, q in zip(r, r[1:])) >= 60]
            out += runs
        if dropped:
            band.set("d", to_d(out))
            band.set("stroke-linecap", "round")
            total += dropped
    if total:
        log.append(("R7", "background band trimmed", "", "", f"{total} segments"))
    return total


def edit_map(path, model, out_path, review_path=None, scope=None, added_path=None, kegg=None, relayout=False,
             hubs=None, bridges="off", one_area=False, map_names=frozenset(), cofactors="added", map_name=None):
    mp = Map(path)
    log = []
    if mp.root.get("data-transport"):  # written whole from the model by transport_map.py
        import shutil
        shutil.copy(path, out_path)
        return []
    marks = collections.defaultdict(set)  # element key -> highlight kind, for the review copy
    rmap = mp.reactions()
    import layout
    before_boxes = [b for b in (layout.node_box(mp, n) for n in mp.layer["nodes"]) if b]

    # R1 metabolite ids
    unresolved = []  # reported after the removals, which may delete these nodes
    for g in mp.nodes("met"):
        cls = mp.classes(g)
        mid, linked = cls[1], [c for c in cls[2:] if is_rxn(c)]
        if mid in model.mets:
            continue
        new = model.fix_met_id(mid, linked)
        if new:
            g.set("class", " ".join(["met", new] + cls[2:]))
            log.append(("R1", "metabolite id", mid, new, model.name(new)))
            marks["id fixed"].add(id(g))
        else:
            unresolved.append(g)

    # R2/R3 reactions no longer in the model
    for rid, g in list(rmap.items()):
        if rid in model.rxns:
            continue
        embedded = [m for m in re.findall(r"MAR\d{5}|r_\d{4}", rid) if m in model.rxns and m != rid]
        rep = (embedded[0], "legacy id") if len(set(embedded)) == 1 else model.replacement(rid)
        if rep and rep[0] not in rmap:
            new, how = rep
            g.set("id", new)
            for t in g.iter(SVG + "text"):
                if t.text and rid in t.text:
                    t.text = t.text.replace(rid, new)
            for layer in ("fes", "fehs", "ees"):
                for p in mp.edges(layer, rid):
                    p.set("class", p.get("class").replace(rid, new))
            for n in mp.nodes("met") + mp.nodes("enz"):
                cls = mp.classes(n)
                if rid in cls:
                    n.set("class", " ".join(new if c == rid else c for c in cls))
            rmap[new] = rmap.pop(rid)
            marks["relabelled"].add(id(g))
            log.append(("R2", f"reaction relabelled ({how})", rid, new, model.rxns[new].get("name", "")))
            continue
        marks["removed"].add(id(g))
        log.append(("R3", "reaction removed", rid, "", text(model.old_rxns.get(rid, {}).get("name", ""))))

    removed = [rid for rid in rmap if rid not in model.rxns]
    review_root = copy.deepcopy(mp.root) if review_path else None

    def detach(node, rid):
        cls = [c for c in mp.classes(node) if c != rid]
        node.set("class", " ".join(cls))
        return [c for c in cls[2:] if is_rxn(c)]

    gone = []  # edge segments removed from the map, for trimming the background band (R7)
    records = []  # geometry of removed reactions, for the layout rules D1/D2
    orphans = []  # main metabolite nodes left without a reaction; removed after the layout rules
    for rid in removed:
        g = rmap[rid]
        fe = [sp for p in mp.edges("fes", rid) for sp in (subpaths(p.get("d")) or [])]
        mains = [(n, translate(n), mp.classes(n)[1]) for n in mp.nodes("met")
                 if rid in mp.classes(n) and n.find(SVG + "ellipse") is not None]
        records.append({"id": rid, "pos": translate(g), "fe": fe, "mains": mains,
                        "name": text(model.old_rxns.get(rid, {}).get("name", ""))})
    for rid in removed:
        g = rmap.pop(rid)
        mp.remove(g)
        for layer in ("fes", "fehs", "ees"):
            for p in mp.edges(layer, rid):
                if layer == "fes":
                    gone += segments(subpaths(p.get("d")) or [])
                mp.remove(p)
        for n in mp.nodes("met") + mp.nodes("enz"):
            if rid in mp.classes(n) and not detach(n, rid):
                if n.find(SVG + "ellipse") is not None:
                    orphans.append(n)  # may still be the end of a restored path (D1/D2)
                else:
                    mp.remove(n)

    # R9 swapped metabolites: a node that serves only one reaction and is no longer part of it takes the id
    # and name of a metabolite of that reaction that is not drawn, in the same role (substrate or product)
    import layout
    for rid in list(mp.reactions()):
        if rid not in model.rxns or rid not in model.old_rxns:
            continue
        st = model.rxns[rid]["stoich"]
        old = {}  # role of each metabolite in any earlier version, the maps' own version first
        for version in model.old_versions:
            for m, c in version.get(rid, {}).get("stoich", {}).items():
                old.setdefault(model.fix_met_id(m, []) or m, c)
        nodes = [n for n in mp.nodes("met") if rid in mp.classes(n)]
        missing = [m for m in st if m not in {mp.classes(n)[1] for n in nodes}]
        for n in nodes:
            cls = mp.classes(n)
            mid = cls[1]
            if mid in st or [c for c in cls[2:] if is_rxn(c)] != [rid] or not old.get(mid):
                continue
            is_main = n.find(SVG + "ellipse") is not None
            cand = [m for m in missing if (st[m] > 0) == (old[mid] > 0)
                    and (model.name(m) not in layout.COFACTORS) == is_main]
            if is_main:
                # a main metabolite that is drawn elsewhere is linked there instead (D5), not duplicated here
                cand = [m for m in cand if not any(mp.classes(x)[1] == m and x.find(SVG + "ellipse") is not None
                                                   for x in mp.nodes("met"))]
            if not cand:
                continue
            new = cand[0]
            missing.remove(new)
            n.set("class", " ".join(["met", new] + cls[2:]))
            fresh = (layout.new_main if is_main else lambda *a: layout.new_side(*a, False))(
                new, model.name(new), [rid], 0, 0)
            for child in list(n):
                if child.tag != SVG + ("ellipse" if is_main else "circle"):
                    n.remove(child)
            for child in list(fresh):
                if child.tag not in (SVG + "ellipse", SVG + "circle"):
                    if not is_main:
                        left = (translate(n) or (0, 0))[0] < (translate(mp.reactions()[rid]) or (0, 0))[0]
                        child.set("text-anchor", "end" if left else "start")
                        child.set("x", "-12" if left else "12")
                    n.append(child)
            log.append(("R9", "metabolite relabelled in place (swapped in the reaction)", mid, rid,
                        f"{model.name(mid)} -> {model.name(new)}"))

    # R4 metabolites no longer part of a drawn reaction
    for n in mp.nodes("met"):
        cls = mp.classes(n)
        mid = cls[1]
        for rid in [c for c in cls[2:] if is_rxn(c)]:
            if rid not in model.rxns or mid in model.rxns[rid]["stoich"]:
                continue
            pos = translate(n)
            cut = 0
            for layer in ("fes", "fehs"):
                for p in mp.edges(layer, rid):
                    sps = subpaths(p.get("d"))
                    if sps is None or pos is None:
                        continue
                    keep = [s for s in sps if not (near(s[0], pos, node_radius(n)) or near(s[-1], pos, node_radius(n)))]
                    if layer == "fes":
                        gone += segments([s for s in sps if s not in keep])
                    if len(keep) != len(sps):
                        cut += len(sps) - len(keep)
                        p.set("d", to_d(keep))
            left = detach(n, rid)
            marks["detached"].add(mid + rid)
            log.append(("R4", "metabolite detached from reaction", mid, rid,
                        f"{model.name(mid)}; {cut} edge pieces removed" + ("; node removed" if not left else "")))
            if not left:
                mp.remove(n)
                break
            cls = mp.classes(n)

    # D1-D3 layout: keep pathway shapes, restore connections, add the subsystem's missing reactions
    import layout
    dr = layout.Drawer(mp, model, log)
    import kegg_layout
    dr.kegg = kegg or []
    if hubs:
        dr.hubs = hubs
    # the map's own compartment: the most drawn one outside the boxes, or on a blank map the most used by
    # its reactions
    base = None
    if scope:
        boxes = layout.compartment_boxes(mp).values()
        letters = collections.Counter(comp_of(mp.classes(n)[1]) for n in mp.nodes("met")
                                      if not any(b[0] <= p[0] <= b[2] and b[1] <= p[1] <= b[3]
                                                 for b in boxes for p in [layout.translate(n) or (-1, -1)]))
        if not letters:
            letters = collections.Counter(comp_of(m) for r in scope.missing(set()) for m in model.rxns[r]["stoich"])
        base = letters.most_common(1)[0][0] if letters else None
    if one_area and scope and scope.subsystems and not layout.compartment_boxes(mp) and base:
        dr.one_area = layout.COMPARTMENT.get(base, base)
    redrawn = bool(relayout and kegg) and kegg_layout.relayout(dr, kegg[0], scope, comp=base or "c")
    if scope is not None:
        if not redrawn:
            layout.restore_paths(dr, records, scope)
            layout.respace_rings(dr)
        if scope:
            comps = layout.compartments_on_map(mp) | ({layout.COMPARTMENT[c] for c in scope.compartments}
                                                      if scope.compartments else set())
            boxed = set(layout.compartment_boxes(mp))
            if not comps - boxed - set(layout.ALIASES) and base:  # no unboxed compartment named: the map's own
                comps.add(layout.COMPARTMENT.get(base, ""))
            layout.add_missing(dr, scope, comps)
            if bridges != "off":
                import bridges as B
                limit, _, steps = bridges.partition(":")
                B.add_bridges(dr, scope, comps, None if limit == "all" else int(limit), int(steps or 1))
        layout.complete_reactions(dr)
        import placement
        placement.place_pending(dr, dr.hubs)
        import connect
        connect.connect(dr)
        if bridges != "off":
            import bridges as B
            B.draw_bridges(dr)
            B.check_bridges(dr)
    gone += dr.gone
    for n in orphans:
        if n.getparent() is not None and not [c for c in mp.classes(n)[2:] if is_rxn(c)]:
            mp.remove(n)

    # a new map: its title centred over the final width, where that place is free
    if mp.root.get("data-new"):
        del mp.root.attrib["data-new"]
        for g in mp.root.iter(SVG + "g"):
            t = g.find(SVG + "text") if g.get("class") == "subsystem" else None
            if t is not None and t.text:
                size = float(t.get("font-size") or 100)
                half = 0.28 * size * len(t.text)
                mid = float(mp.root.get("width")) / 2
                x = mid if t.get("text-anchor") == "middle" else mid - half
                box = (mid - half, size, mid + half, size * 2.3)
                if dr.space.free(box, 20):
                    t.set("transform", f"matrix(1,0,0,1,{x:.1f},{size * 2:.1f})")
                break

    # context reactions (D2 paths and D14 bridges of other subsystems) are drawn without the subsystem's band,
    # as the original drawings do
    context = {rid for rid, _ in dr.band_add if scope is not None and scope.subsystems and not scope.contains(rid)}
    layout.extend_band(mp, model, dr.band_add, skip=context)
    # R7 background band: drop the parts that run along no remaining edge (edges removed, or emptied by D10)
    trim_band(mp, log)
    if scope:
        rename_bands(mp, model, scope, map_names, log)
    if redrawn:  # D8: the band follows the new drawing
        layout.regenerate_band(mp)

    # R5 genes
    for n in mp.nodes("enz"):
        cls = mp.classes(n)
        gene = cls[1]
        for rid in [c for c in cls[2:] if is_rxn(c)]:
            if rid in model.rxns and gene not in model.rxns[rid]["genes"]:
                left = detach(n, rid)
                log.append(("R5", "gene detached from reaction", gene, rid, model.symbol.get(gene, "")))
                if not left:
                    mp.remove(n)
                    break
        else:
            sym = model.symbol.get(gene)
            t = n.find(SVG + "text")
            want = layout.gene_label(gene, sym)
            have = [x.strip() for x in t.itertext() if x.strip()] if t is not None else []
            if t is not None and have and have != want:
                log.append(("R5", "gene label", gene, " / ".join(have), " / ".join(want)))
                layout.set_gene_text(t, gene, sym)

    # report what is left for review
    for g in unresolved:
        if g.getparent() is not None:
            cls = mp.classes(g)
            log.append(("review", "metabolite id not resolved", cls[1], " ".join(cls[2:]), model.name(cls[1])))
    drawn = set(mp.reactions())
    for rid in sorted(drawn):
        r = model.rxns[rid]
        on_map = {mp.classes(n)[1] for n in mp.nodes("met") if rid in mp.classes(n)}
        missing = [m for m in r["stoich"] if m not in on_map]
        main_missing = [m for m in missing if model.name(m) not in CURRENCY]
        if main_missing:
            log.append(("review", "main metabolites not drawn", rid, " ".join(main_missing),
                        ", ".join(model.name(m) for m in main_missing)))
        elif missing:
            log.append(("review", "side metabolites not drawn", rid, " ".join(missing),
                        ", ".join(model.name(m) for m in missing)))
        genes_drawn = {mp.classes(n)[1] for n in mp.nodes("enz") if rid in mp.classes(n)}
        if r["genes"] - genes_drawn:
            log.append(("review", "genes not drawn", rid, " ".join(sorted(r["genes"] - genes_drawn)),
                        ", ".join(model.symbol.get(g, g) for g in sorted(r["genes"] - genes_drawn))))

    if scope is not None:
        layout.pack_genes(mp, log, {e[3] for e in log if e[0] == "R5" and e[1] == "gene detached from reaction"})
        layout.close_space(mp, before_boxes, log)
        import compact
        compact.remove_empty_boxes(mp, log)
        compact.settle_boxes(mp, log)
        # the map's title: the text of a subsystem map's band, or of a band named after the map itself (Yeast-GEM's
        # compartment maps); other bands' texts are labels within the drawing
        title = next((g.find(SVG + "text") for g in mp.root.iter(SVG + "g") if g.get("class") == "subsystem"
                      and g.find(SVG + "text") is not None and (scope.subsystems or g.get("id") == map_name)), None)
        if cofactors != "off":  # before D13, which closes the space and centres the content as it ends up
            import declutter
            declutter.declutter(dr, cofactors)
        compact.compact(mp, log, title)
        if scope:
            settle(dr, scope, comps, bridges, title, cofactors, log)
        # R7 once more: white space closed after the last trim can take a band piece off its edge; then D13 again,
        # as the band is part of the content it centres
        if trim_band(mp, log):
            compact.compact(mp, log, title)

    mp.root.set("data-modelversion", model.version)
    mp.write(out_path)
    if added_path:
        write_added(copy.deepcopy(mp.root), dr.placed, added_path)

    # D16 rounds repeat some notes: one per change, and one review note per item
    seen, rows = set(), []
    for row in log:
        key = tuple(row[:3]) if row[0] == "review" else tuple(row)
        if key not in seen:
            seen.add(key)
            rows.append(tuple(row))
    log[:] = rows
    if review_path:
        write_review(review_root, removed, marks, log, review_path)
    return log


def settle(dr, scope, comps, bridges, title, cofactors, log, rounds=3):
    """Repeat the steps that depend on where things are drawn (D5 metabolites of drawn reactions, D10 links, D14
    bridges, with the band, cofactor dots and white space after them) on the finished drawing until a round
    changes nothing, so that the editor run again on its own maps finds nothing to do. Parts are not moved again
    to be linked: a second run would not move them either (it places nothing)."""
    import compact
    import connect
    import layout
    mp = dr.mp
    for _ in range(rounds):
        before = etree.tostring(mp.root)
        dr.space = layout.Space(mp)
        dr.band_add, dr.gone, dr.pending = [], [], []
        layout.complete_reactions(dr)
        connect.connect(dr, movable=0)
        if bridges != "off" and scope.subsystems:
            import bridges as B
            limit, _, steps = bridges.partition(":")
            dr.bridge_groups, dr.bridge_via = [], {}
            B.add_bridges(dr, scope, comps, None if limit == "all" else int(limit), int(steps or 1))
            B.draw_bridges(dr)
            B.check_bridges(dr)
        if etree.tostring(mp.root) == before:
            return
        context = {rid for rid, _ in dr.band_add if scope.subsystems and not scope.contains(rid)}
        layout.extend_band(mp, dr.model, dr.band_add, skip=context)
        trim_band(mp, log)
        if cofactors != "off":
            import declutter
            declutter.declutter(dr, cofactors)
        compact.compact(mp, log, title)


def rename_bands(mp, model, scope, map_names, log):
    """R12: a band (subsystem group) whose id is not a subsystem of the model gets the current name: on a subsystem
    map the map's subsystem, on a compartment map the subsystem of most (at least 60%) of the reactions on it. A
    band that would take the id of another band joins it. A label or title that repeats the old id gets the new
    name. Remnants (trimmed by R7) with no reaction on them are removed, with their label; with one or two they
    join their subsystem's band without the label. Bands named after their map (Yeast-GEM compartment maps) are
    left alone."""
    import layout
    current = set()
    for r in model.rxns.values():
        subs = r.get("subsystem")
        current |= {x for x in (subs if isinstance(subs, list) else [subs]) if x}
    groups = layout.band_groups(mp)
    by_id = {g.get("id"): (g, p) for g, p in groups}
    rpos = {rid: translate(n) for rid, n in mp.reactions().items()}
    norm = lambda x: " ".join((x or "").split()).lower()  # noqa: E731
    for g, path in groups:
        gid = g.get("id")
        if gid in current or gid in map_names:
            continue
        if scope.subsystems and len(scope.subsystems) == 1:
            target = next(iter(scope.subsystems))
        else:
            segs = [(a, b) for sp in (subpaths(path.get("d") or "") or []) for a, b in zip(sp, sp[1:])]
            count = collections.Counter()
            for rid, q in rpos.items():
                if q and rid in model.rxns and any(seg_dist(q, a, b) <= 60 for a, b in segs):
                    subs = model.rxns[rid].get("subsystem")
                    subs = subs if isinstance(subs, list) else [subs]
                    if subs and subs[0]:
                        count[subs[0]] += 1
            if not count:  # a remnant with no reaction on it (or only a label): it goes
                g.getparent().remove(g)
                log.append(("R12", "band of a subsystem not in the model removed (no reaction on it)", gid, "", ""))
                continue
            target, n = count.most_common(1)[0]
            if sum(count.values()) < 3:  # a remnant with one or two reactions: it joins their band, without its label
                for t in list(g.iter(SVG + "text")):
                    t.getparent().remove(t)
            elif n < 0.6 * sum(count.values()):
                log.append(("review", "band name not current, reactions on it of several subsystems", gid,
                            ", ".join(f"{k} {v}" for k, v in count.most_common(3)), ""))
                continue
        for t in g.iter(SVG + "text"):
            if norm("".join(t.itertext())) == norm(gid):
                for child in list(t):
                    t.remove(child)
                t.text = target
        if target in by_id and by_id[target][0] is not g:
            og, other = by_id[target]
            other.set("d", ((other.get("d") or "") + " " + (path.get("d") or "")).strip())
            for t in list(g.iter(SVG + "text")):  # the label goes with it, where it was
                og.append(t)
            g.getparent().remove(g)
            log.append(("R12", "band joined the band of its current subsystem", gid, target, ""))
        else:
            g.set("id", target)
            by_id[target] = (g, path)
            log.append(("R12", "band renamed to its current subsystem", gid, target, ""))


def write_added(root, placed, path):
    """A copy of the edited map with added and moved reactions coloured."""
    style = etree.SubElement(root, SVG + "style")
    style.text = (".rv-added .shape{fill:#2f9e44!important;stroke:#2f9e44!important}"
                  "path.rv-added{stroke:#2f9e44!important;stroke-width:7!important}"
                  "#fehs path.rv-added{fill:#2f9e44!important}")
    for el in root.iter(SVG + "g", SVG + "path"):
        cls = (el.get("class") or "").split()
        if not cls:
            continue
        rids = [el.get("id")] if cls == ["rea"] else [c for c in cls[1:] if is_rxn(c)]
        if cls[0] in ("rea", "fe", "feh", "ee", "met", "enz") and rids and all(r in placed for r in rids):
            el.set("class", " ".join(cls + ["rv-added"]))
    etree.ElementTree(root).write(path, encoding="utf-8")


def write_review(root, removed, marks, log, path):
    """A copy of the map before removal, with the changes coloured."""
    style = etree.SubElement(root, SVG + "style")
    style.text = "\n".join([
        ".rv-removed .shape{fill:%s!important;stroke:%s!important}" % (HIGHLIGHT["removed"], HIGHLIGHT["removed"]),
        "path.rv-removed{stroke:%s!important;stroke-width:8!important}" % HIGHLIGHT["removed"],
        ".rv-relabelled .shape{fill:%s!important}" % HIGHLIGHT["relabelled"],
        ".rv-idfixed .shape{stroke:%s!important;stroke-width:6!important}" % HIGHLIGHT["id fixed"],
        ".rv-detached .shape{stroke:%s!important;stroke-width:6!important;stroke-dasharray:6 4}" % HIGHLIGHT["detached"],
    ])
    relab = {e[3] for e in log if e[0] == "R2"}
    detached = {(e[2], e[3]) for e in log if e[0] == "R4"}
    for el in root.iter(SVG + "g", SVG + "path"):
        cls = (el.get("class") or "").split()
        if not cls:
            continue
        rid = el.get("id") if cls == ["rea"] else None
        kind = None
        if (rid and rid in removed) or (cls[0] in ("fe", "feh", "ee") and cls[1] in removed):
            kind = "rv-removed"
        elif cls[0] in ("met", "enz") and cls[2:] and all(c in removed for c in cls[2:] if is_rxn(c)):
            kind = "rv-removed"
        elif rid and rid in relab:
            kind = "rv-relabelled"
        elif cls[0] == "met" and any((cls[1], r) in detached for r in cls[2:]):
            kind = "rv-detached"
        elif cls[0] == "met" and cls[1] in {e[3] for e in log if e[0] == "R1"}:
            kind = "rv-idfixed"
        if kind:
            el.set("class", " ".join(cls + [kind]))
    etree.ElementTree(root).write(path, encoding="utf-8")


def read_map_table(path):
    """{svg file name: component name} from a subsystemSVG.tsv or compartmentSVG.tsv."""
    out = {}
    for line in open(path):
        if line.strip() and not line.startswith(("#", "@")):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                out[f[2]] = f[0]
    return out


def share_missing(model, paths, letters):
    """Assign each missing reaction of a compartment drawn on several maps (the cytosol) to one map: the one
    that already draws most of its main metabolites, otherwise the one with the fewest additions so far."""
    import layout
    scope = layout.Scope(model, compartments=letters)
    drawn, mains = {}, {}
    for p in paths:
        mp = Map(p)
        drawn[p] = set(mp.reactions())
        mains[p] = collections.Counter()
        for n in mp.nodes("met"):
            if n.find(SVG + "ellipse") is not None:
                cls = mp.classes(n)
                mains[p][model.fix_met_id(cls[1], cls[2:]) or cls[1]] += 1
    everywhere = set().union(*drawn.values())
    assigned = {p: set() for p in paths}
    for rid in scope.missing(everywhere):
        ms = [m for m in model.rxns[rid]["stoich"] if model.name(m) not in layout.COFACTORS]
        best = max(paths, key=lambda p: (sum(1 for m in ms if mains[p][m]), -len(assigned[p])))
        assigned[best].add(rid)
    return assigned


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model_yml")
    ap.add_argument("model_dir")
    ap.add_argument("old_ymls", nargs="*")
    ap.add_argument("--maps", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--review", action="store_true")
    ap.add_argument("--map-table", help="subsystemSVG.tsv, to know each map's subsystem (enables D1-D6)")
    ap.add_argument("--compartment-table", help="compartmentSVG.tsv, to know each compartment map's compartment")
    ap.add_argument("--custom-table", help="customSVG.tsv, for maps linked as custom maps")
    ap.add_argument("--kegg-dir", help="cached KGML files and subsystem_pathways.json (enables D7)")
    ap.add_argument("--gene-label", choices=["name", "orf", "both"], default="name",
                    help="gene boxes show the gene name, the gene id, or the name over the id")
    ap.add_argument("--parts", help="JSON {map file: [reaction ids]}: the reactions a compartment map part draws")
    ap.add_argument("--one-area", action="store_true",
                    help="subsystem maps without compartment boxes draw all compartments in one area (Yeast-GEM)")
    ap.add_argument("--kegg-mapping", help="subsystem -> KEGG maps JSON (default: subsystem_pathways.json in --kegg-dir)")
    ap.add_argument("--bridges", default="all",
                    help="D14: add one-step bridges of other subsystems between unconnected parts of subsystem maps, "
                         "through metabolites in at most N reactions (a number), 'all' (no limit; default) or 'off'; "
                         "append ':2' to also allow two-step bridges (e.g. 30:2)")
    ap.add_argument("--cofactors", choices=["added", "all", "off"], default="added",
                    help="D15: move cofactor dots that overlap other elements to a free place around their reaction: "
                         "those of reactions placed in this run (default), all, or none")
    ap.add_argument("--hubs", choices=["spread", "mid", "scale"], default="mid",
                    help="crowded KEGG modules: move only the close compounds apart, or scale the module up")
    ap.add_argument("--kegg-relayout", nargs="*", default=[],
                    help="map files whose cytosol is redrawn as their first KEGG map lays it out (D8)")
    a = ap.parse_args()
    import layout
    subsystem_of = read_map_table(a.map_table) if a.map_table else {}
    own_name = {f[2].strip(): f[1] for t in (a.map_table, a.compartment_table) if t
                for f in (line.rstrip("\n").split("\t") for line in open(t))
                if len(f) >= 3 and not f[0].startswith(("#", "@"))}  # file -> map name
    map_names = frozenset(own_name.values())
    compartment_of = read_map_table(a.compartment_table) if a.compartment_table else {}
    layout.GENE_LABEL = a.gene_label
    model = Model(a.model_yml, a.model_dir, a.old_ymls)
    letters_of = {v: {k} for k, v in layout.COMPARTMENT.items()}
    for k, v in load_model.compartments.items():  # map tables use the model's own names
        letters_of.setdefault(v, set()).add(k)
    letters_of["Mitochondria"] = {"m", "i"}
    # compartments drawn on several maps share their missing reactions
    only = {}
    by_comp = collections.defaultdict(list)
    for p in a.maps:
        if os.path.basename(p) in compartment_of:
            by_comp[compartment_of[os.path.basename(p)]].append(p)
    if a.parts:  # parts of a large compartment map, set up front (new maps have nothing drawn to go by)
        import json
        for name, rids in json.load(open(a.parts)).items():
            for p in a.maps:
                if os.path.basename(p) == name:
                    only[p] = set(rids)
    for comp, paths in by_comp.items():
        if len(paths) > 1 and not any(p in only for p in paths):
            for p, rids in share_missing(model, paths, letters_of.get(comp, set())).items():
                only[p] = rids
    os.makedirs(a.out, exist_ok=True)
    allrows = []
    linked = set(subsystem_of) | set(compartment_of)
    if a.custom_table:
        linked |= set(read_map_table(a.custom_table)) | {f[1] for f in (line.rstrip("\n").split("\t") for line in
                                                         open(a.custom_table)) if len(f) >= 2 and f[1].endswith(".svg")}
    for p in a.maps:
        name = os.path.basename(p)
        if a.map_table and name not in set(subsystem_of) | set(compartment_of):
            # custom maps (and unlinked files) show content beyond the model; they are left as they are
            print(name, "custom or unlinked map; left unchanged", flush=True)
            import shutil
            shutil.copy(p, os.path.join(a.out, name))
            allrows.append((name, "review", "custom or unlinked map, left unchanged", "", "", ""))
            continue
        if not a.map_table:
            scope = None
        elif name in subsystem_of:
            scope = layout.Scope(model, subsystems={subsystem_of[name]})
        elif name in compartment_of:
            scope = layout.Scope(model, compartments=letters_of.get(compartment_of[name], set()), only=only.get(p))
        else:
            scope = layout.Scope(model)
        kegg = []
        if a.kegg_dir and name in subsystem_of:
            import json
            import kegg_layout
            pathways = json.load(open(a.kegg_mapping or os.path.join(a.kegg_dir, "subsystem_pathways.json"))).get(
                subsystem_of[name], [])
            kegg = [kegg_layout.Kgml(os.path.join(a.kegg_dir, pw + ".xml")) for pw in pathways
                    if os.path.exists(os.path.join(a.kegg_dir, pw + ".xml"))]
        try:
            log = edit_map(p, model, os.path.join(a.out, name),
                           os.path.join(a.out, name.replace(".svg", ".review.svg")) if a.review else None,
                           scope, os.path.join(a.out, name.replace(".svg", ".added.svg")) if a.review else None,
                           kegg, name in a.kegg_relayout, a.hubs, a.bridges, a.one_area, map_names, a.cofactors,
                           own_name.get(name))
        except Exception as err:  # keep the map unchanged and go on with the others
            import shutil
            import traceback
            traceback.print_exc()
            shutil.copy(p, os.path.join(a.out, name))
            log = [("error", "map left unchanged", type(err).__name__, str(err)[:200], "")]
        allrows += [(name,) + row for row in log]
        counts = collections.Counter(r[1] for r in log)
        print(name, dict(counts), flush=True)
    with open(os.path.join(a.out, "changes.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["map", "rule", "change", "id", "detail", "name"])
        w.writerows(allrows)


if __name__ == "__main__":
    main()
