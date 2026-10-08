#!/usr/bin/env python3
"""Write Metabolic Atlas SVG maps in other formats, for the Human-maps and Yeast-maps repositories.

For each map: SBGN-ML (process description: metabolites as simple chemicals, reactions as processes,
genes as macromolecules catalysing them, compartments), SBML Level 3 with the layout and groups packages
(the reactions with their stoichiometry and gene modifiers from the model, every drawn node and edge in
the layout, one group per subsystem), an Escher map (JSON), and a PNG image. Positions and edge paths are
those of the SVG.

Usage: export_formats.py <model.yml> <model dir with reactions/metabolites/genes.tsv> <out dir>
       --maps a.svg ... [--png-width 3000]
Writes <out>/<format>/<map>.<ext> for format svg, sbgn, sbml, escher (.json) and png. Needs python-libsbml,
lxml, cairosvg and Pillow. With the libsbgn schema at schema/SBGN.xsd next to this file, every SBGN file is
validated against it; with Escher's map schema at schema/escher_1-0-0.json (and jsonschema), every Escher map
too; every SBML file is checked with libsbml, and every Escher map with Escher's own consistency checks.
"""

import argparse
import collections
import io
import math
import os
import re
import shutil
import sys

import libsbml
from lxml import etree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mapedit  # noqa: E402
import layout as L  # noqa: E402

SBGN = "http://sbgn.org/libsbgn/0.3"
SBGN_PD = "http://identifiers.org/combine.specifications/sbgn.pd.level-1.version-1.3"


def sid(s):
    """A valid SBML/SBGN identifier."""
    s = re.sub(r"[^A-Za-z0-9_]", "_", s)
    return s if re.match(r"[A-Za-z_]", s) else "_" + s


class MapData:
    """What an SVG map draws: compartment boxes, metabolite and gene nodes, reactions and their edges."""

    def __init__(self, path, model):
        mp = mapedit.Map(path)
        self.mp, self.model = mp, model
        self.name = os.path.basename(path)[:-4]
        self.width, self.height = float(mp.root.get("width")), float(mp.root.get("height"))
        # a subsystem map's one band carries the map title; the bands of a compartment map carry subsystem names
        bands = [g for g in mp.root.iter(L.SVG + "g") if g.get("class") == "subsystem"]
        title = " ".join(text for text, *_ in text_blocks(bands[0])) if len(bands) == 1 else ""
        self.title = title or self.name[:1].upper() + self.name[1:].replace("_", " ")
        self.boxes = [(g.get("id"), b) for g in mp.root.iter(L.SVG + "g") if g.get("class") == "compartment"
                      for b in [L.compartment_boxes(mp).get(g.get("id"))] if b]
        self.rpos = {r: L.translate(g) for r, g in mp.reactions().items() if L.translate(g) and r in model.rxns}
        self.mets, self.genes = [], []
        for n in mp.nodes("met"):
            cls, p = mp.classes(n), L.translate(n)
            if p and len(cls) > 1:
                main = n.find(L.SVG + "ellipse") is not None
                box = (p[0] - 30, p[1] - 20, 60, 40) if main else (p[0] - 7, p[1] - 7, 14, 14)
                self.mets.append({"met": cls[1], "pos": p, "box": box, "main": main, "label": svg_label(n),
                                  "rids": [c for c in cls[2:] if c in self.rpos]})
        for n in mp.nodes("enz"):
            cls, p = mp.classes(n), L.translate(n)
            if p and len(cls) > 1:
                label = " / ".join(x.strip() for x in n.find(L.SVG + "text").itertext() if x.strip()) \
                    if n.find(L.SVG + "text") is not None else cls[1]
                self.genes.append({"gene": cls[1], "pos": (p[0] + 30, p[1] + 15), "box": (p[0], p[1], 60, 30),
                                   "label": label, "rids": [c for c in cls[2:] if c in self.rpos]})
        self.edges = collections.defaultdict(list)   # reaction -> [(met node index, polyline node->reaction)]
        self.gene_edges = collections.defaultdict(list)  # reaction -> [(gene node index, polyline gene->reaction)]
        for rid, rp in self.rpos.items():
            nodes = [(i, m["pos"]) for i, m in enumerate(self.mets) if rid in m["rids"]]
            for path in mp.edges("fes", rid):
                for sp in L.subpaths(path.get("d") or "") or []:
                    self._assign(self.edges[rid], nodes, sp, rp)
            gnodes = [(i, g["pos"]) for i, g in enumerate(self.genes) if rid in g["rids"]]
            for path in mp.edges("ees", rid):
                for sp in L.subpaths(path.get("d") or "") or []:
                    self._assign(self.gene_edges[rid], gnodes, sp, rp)

    @staticmethod
    def _assign(out, nodes, sp, rp):
        """Give an edge piece to the node at one of its ends, oriented from the node to the reaction."""
        if not nodes or len(sp) < 2:
            return
        best = min(((min(L.length(L.sub(sp[0], p)), L.length(L.sub(sp[-1], p))), i, p) for i, p in nodes))
        _, i, p = best
        pts = list(sp) if L.length(L.sub(sp[0], p)) <= L.length(L.sub(sp[-1], p)) else list(reversed(sp))
        if L.length(L.sub(pts[-1], rp)) > 80:  # a piece that does not reach the reaction: extend it
            pts.append(rp)
        out.append((i, pts))

    def compartment_of_box(self, p):
        return next((name for name, b in self.boxes if b[0] <= p[0] <= b[2] and b[1] <= p[1] <= b[3]), None)


def font_size(t):
    v = t.get("font-size") or (t.getparent().get("font-size") if t.getparent() is not None else None) or "18"
    return float(re.sub(r"[^0-9.]", "", v) or 18)


def svg_label(n):
    """Where a metabolite node draws its label, relative to the node: (anchor, x, baseline, font size), with
    anchor start, end or middle and x the left end, the right end or the centre of its lines."""
    lines = []
    for t in n.iter(L.SVG + "text"):
        content = " ".join(x.strip() for x in t.itertext() if x.strip())
        if not content:
            continue
        size = font_size(t)
        anchor = t.get("text-anchor") or t.getparent().get("text-anchor") or "start"
        x, y = float(t.get("x") or 0), float(t.get("y") or 0)
        w = L.label_width(content, size)
        left = x - w / 2 if anchor == "middle" else x - w if anchor == "end" else x
        lines.append((anchor, left, left + w, y, size))
    if not lines:
        return ("middle", 0.0, 6.0, 18.0)
    anchors = {a for a, *_ in lines}
    base, size = sum(x[3] for x in lines) / len(lines), lines[0][4]
    if anchors == {"start"}:
        return ("start", min(x[1] for x in lines), base, size)
    if anchors == {"end"}:
        return ("end", max(x[2] for x in lines), base, size)
    return ("middle", (min(x[1] for x in lines) + max(x[2] for x in lines)) / 2, base, size)


def text_blocks(group):
    """The labels of a compartment or subsystem group: each text, or each group of text lines read as one, as
    (text, centre x, mean baseline, font size)."""
    out = []
    for child in group:
        texts = [child] if child.tag == L.SVG + "text" else list(child.iter(L.SVG + "text")) \
            if child.tag == L.SVG + "g" else []
        lines = []
        for t in texts:
            content = " ".join(x.strip() for x in t.itertext() if x.strip())
            b = L.text_box(t)
            if content and b:
                size = font_size(t)
                lines.append((content, (b[0] + b[2]) / 2, b[3] - 0.4 * size, size))
        if lines:
            out.append((" ".join(x[0] for x in lines), sum(x[1] for x in lines) / len(lines),
                        sum(x[2] for x in lines) / len(lines), max(x[3] for x in lines)))
    return out


# ----------------------------------------------------------------------------- SBGN-ML

def write_sbgn(d, path):
    root = etree.Element("{%s}sbgn" % SBGN, nsmap={None: SBGN})
    m = etree.SubElement(root, "{%s}map" % SBGN, id=sid(d.name), version=SBGN_PD)

    def glyph(cls, gid, label, box, compartment=None):
        g = etree.SubElement(m, "{%s}glyph" % SBGN, **{"class": cls, "id": gid})
        if compartment:
            g.set("compartmentRef", compartment)
        if label is not None:
            etree.SubElement(g, "{%s}label" % SBGN, text=label)
        etree.SubElement(g, "{%s}bbox" % SBGN, x=f"{box[0]:.1f}", y=f"{box[1]:.1f}", w=f"{box[2]:.1f}", h=f"{box[3]:.1f}")
        return g
    comp_ids = {}
    for i, (name, b) in enumerate(d.boxes):
        comp_ids[name] = f"comp{i + 1}"
        glyph("compartment", comp_ids[name], name, (b[0], b[1], b[2] - b[0], b[3] - b[1]))
    count = collections.Counter(x["met"] for x in d.mets)
    for i, x in enumerate(d.mets):
        g = glyph("simple chemical", f"met{i + 1}", d.model.name(x["met"]), x["box"],
                  comp_ids.get(d.compartment_of_box(x["pos"])))
        if count[x["met"]] > 1:  # the clone marker comes before the bounding box
            g.find("{%s}bbox" % SBGN).addprevious(etree.Element("{%s}clone" % SBGN))
    for i, x in enumerate(d.genes):
        glyph("macromolecule", f"gene{i + 1}", x["label"], x["box"], comp_ids.get(d.compartment_of_box(x["pos"])))
    arcs = []  # written after all glyphs, as the schema wants
    for rid, rp in d.rpos.items():
        st = d.model.rxns[rid]["stoich"]
        ins = [pts[-2] for i, pts in d.edges[rid] if st.get(d.mets[i]["met"], 0) < 0 and len(pts) > 1]
        outs = [pts[-2] for i, pts in d.edges[rid] if st.get(d.mets[i]["met"], 0) > 0 and len(pts) > 1]
        a = (sum(p[0] for p in ins) / len(ins), sum(p[1] for p in ins) / len(ins)) if ins else (rp[0] - 1, rp[1])
        b = (sum(p[0] for p in outs) / len(outs), sum(p[1] for p in outs) / len(outs)) if outs else (rp[0] + 1, rp[1])
        u = L.unit(L.sub(b, a)) if L.length(L.sub(b, a)) > 1 else (1.0, 0.0)
        pid = f"process_{sid(rid)}"
        g = glyph("process", pid, rid, (rp[0] - 15, rp[1] - 15, 30, 30))
        p_in, p_out = L.add(rp, u, -25), L.add(rp, u, 25)
        etree.SubElement(g, "{%s}port" % SBGN, id=pid + ".1", x=f"{p_in[0]:.1f}", y=f"{p_in[1]:.1f}")
        etree.SubElement(g, "{%s}port" % SBGN, id=pid + ".2", x=f"{p_out[0]:.1f}", y=f"{p_out[1]:.1f}")
        for k, (i, pts) in enumerate(d.edges[rid]):
            role = st.get(d.mets[i]["met"], 0)
            if role < 0:
                arc(arcs, "consumption", f"{pid}_c{k}", f"met{i + 1}", pid + ".1", pts[:-1] + [p_in])
            elif role > 0:
                line = list(reversed(pts))
                arc(arcs, "production", f"{pid}_p{k}", pid + ".2", f"met{i + 1}", [p_out] + line[1:])
        for k, (i, pts) in enumerate(d.gene_edges[rid]):
            arc(arcs, "catalysis", f"{pid}_e{k}", f"gene{i + 1}", pid, pts)
    m.extend(arcs)
    etree.ElementTree(root).write(path, xml_declaration=True, encoding="UTF-8", pretty_print=True)


def arc(arcs, cls, aid, source, target, pts):
    a = etree.Element("{%s}arc" % SBGN, **{"class": cls, "id": aid, "source": source, "target": target})
    arcs.append(a)
    etree.SubElement(a, "{%s}start" % SBGN, x=f"{pts[0][0]:.1f}", y=f"{pts[0][1]:.1f}")
    for p in pts[1:-1]:
        etree.SubElement(a, "{%s}next" % SBGN, x=f"{p[0]:.1f}", y=f"{p[1]:.1f}")
    etree.SubElement(a, "{%s}end" % SBGN, x=f"{pts[-1][0]:.1f}", y=f"{pts[-1][1]:.1f}")


# ----------------------------------------------------------------------------- SBML with layout

def check(status, what):
    if status not in (libsbml.LIBSBML_OPERATION_SUCCESS, None):
        raise RuntimeError(f"libsbml: {what}: {libsbml.OperationReturnValue_toString(status)}")


def bbox(obj, gid, box):
    obj.setId(gid)
    b = obj.getBoundingBox()
    b.setX(box[0]); b.setY(box[1]); b.setWidth(box[2]); b.setHeight(box[3])  # noqa: E702


def curve(srg, pts):
    c = srg.getCurve()
    for a, b in zip(pts, pts[1:]):
        seg = c.createLineSegment()
        seg.getStart().setX(a[0]); seg.getStart().setY(a[1])  # noqa: E702
        seg.getEnd().setX(b[0]); seg.getEnd().setY(b[1])  # noqa: E702


def write_sbml(d, path):
    model = d.model
    doc = libsbml.SBMLDocument(3, 1)
    check(doc.enablePackage(libsbml.LayoutExtension.getXmlnsL3V1V1(), "layout", True), "layout")
    check(doc.enablePackage(libsbml.GroupsExtension.getXmlnsL3V1V1(), "groups", True), "groups")
    doc.setPackageRequired("layout", False)
    doc.setPackageRequired("groups", False)
    m = doc.createModel()
    m.setId(sid(d.name))
    m.setName(d.title)
    mets = {x for r in d.rpos for x in model.rxns[r]["stoich"]}
    comps = sorted({L.comp_of(x) for x in mets})
    for c in comps:
        comp = m.createCompartment()
        comp.setId(sid(c)); comp.setName(L.COMPARTMENT.get(c, c)); comp.setConstant(True)  # noqa: E702
        comp.setSpatialDimensions(3); comp.setSize(1)  # noqa: E702
    for x in sorted(mets):
        s = m.createSpecies()
        s.setId(sid(x)); s.setName(model.name(x)); s.setCompartment(sid(L.comp_of(x)))  # noqa: E702
        s.setHasOnlySubstanceUnits(False); s.setBoundaryCondition(False); s.setConstant(False)  # noqa: E702
        s.setSBOTerm(247)
    genes = sorted({g for r in d.rpos for g in model.rxns[r]["genes"]})
    gene_comp = {}
    for r in d.rpos:
        for g in model.rxns[r]["genes"]:
            gene_comp.setdefault(g, sid(L.comp_of(next(iter(model.rxns[r]["stoich"])))))
    for g in genes:
        s = m.createSpecies()
        s.setId(sid("gene_" + g)); s.setName(model.symbol.get(g) or g); s.setCompartment(gene_comp[g])  # noqa: E702
        s.setHasOnlySubstanceUnits(False); s.setBoundaryCondition(False); s.setConstant(False)  # noqa: E702
        s.setSBOTerm(252)
    groups = collections.defaultdict(list)
    for rid in d.rpos:
        r = model.rxns[rid]
        rx = m.createReaction()
        rx.setId(sid(rid)); rx.setName(mapedit.text(r.get("name")) or rid)  # noqa: E702
        rx.setReversible((r.get("lower_bound") or 0) < 0); rx.setFast(False); rx.setSBOTerm(176)  # noqa: E702
        for x, c in sorted(r["stoich"].items()):
            ref = rx.createReactant() if c < 0 else rx.createProduct()
            ref.setId(sid(f"{rid}_{x}")); ref.setSpecies(sid(x)); ref.setStoichiometry(abs(c)); ref.setConstant(True)  # noqa: E702
        for g in sorted(r["genes"]):
            mod = rx.createModifier()
            mod.setId(sid(f"{rid}_gene_{g}")); mod.setSpecies(sid("gene_" + g)); mod.setSBOTerm(460)  # noqa: E702
        subs = r.get("subsystem")
        for s in (subs if isinstance(subs, list) else [subs] if subs else []):
            groups[s].append(rid)
    gplug = m.getPlugin("groups")
    for i, (name, rids) in enumerate(sorted(groups.items())):
        grp = gplug.createGroup()
        grp.setId(f"group{i + 1}"); grp.setName(name); grp.setKind(libsbml.GROUP_KIND_PARTONOMY)  # noqa: E702
        for r in rids:
            grp.createMember().setIdRef(sid(r))
    lay = m.getPlugin("layout").createLayout()
    lay.setId("layout")
    lay.setDimensions(libsbml.Dimensions(libsbml.LayoutPkgNamespaces(3, 1, 1), d.width, d.height))
    code = {v: k for k, v in L.COMPARTMENT.items()}
    for i, (name, b) in enumerate(d.boxes):
        c = code.get(name) or code.get(L.ALIASES.get(name, ""), "")
        if sid(c) in [sid(x) for x in comps]:
            cg = lay.createCompartmentGlyph()
            bbox(cg, f"cg{i + 1}", (b[0], b[1], b[2] - b[0], b[3] - b[1]))
            cg.setCompartmentId(sid(c))
    for i, x in enumerate(d.mets):
        if x["met"] not in mets:
            continue
        sg = lay.createSpeciesGlyph()
        bbox(sg, f"sg{i + 1}", x["box"])
        sg.setSpeciesId(sid(x["met"]))
        tg = lay.createTextGlyph()
        bbox(tg, f"tg_sg{i + 1}", x["box"])
        tg.setGraphicalObjectId(f"sg{i + 1}"); tg.setOriginOfTextId(sid(x["met"]))  # noqa: E702
    gene_ids = set(genes)
    for i, x in enumerate(d.genes):
        if x["gene"] not in gene_ids:
            continue
        sg = lay.createSpeciesGlyph()
        bbox(sg, f"gg{i + 1}", x["box"])
        sg.setSpeciesId(sid("gene_" + x["gene"]))
        tg = lay.createTextGlyph()
        bbox(tg, f"tg_gg{i + 1}", x["box"])
        tg.setGraphicalObjectId(f"gg{i + 1}"); tg.setText(x["label"])  # noqa: E702
    for rid, rp in d.rpos.items():
        st = model.rxns[rid]["stoich"]
        rg = lay.createReactionGlyph()
        bbox(rg, f"rg_{sid(rid)}", (rp[0] - 15, rp[1] - 15, 30, 30))
        rg.setReactionId(sid(rid))
        for k, (i, pts) in enumerate(d.edges[rid]):
            met = d.mets[i]["met"]
            if met not in st:
                continue
            srg = rg.createSpeciesReferenceGlyph()
            srg.setId(f"srg_{sid(rid)}_{k}")
            srg.setSpeciesGlyphId(f"sg{i + 1}")
            srg.setSpeciesReferenceId(sid(f"{rid}_{met}"))
            if st[met] < 0:
                srg.setRole(libsbml.SPECIES_ROLE_SUBSTRATE if d.mets[i]["main"] else libsbml.SPECIES_ROLE_SIDESUBSTRATE)
                curve(srg, pts)
            else:
                srg.setRole(libsbml.SPECIES_ROLE_PRODUCT if d.mets[i]["main"] else libsbml.SPECIES_ROLE_SIDEPRODUCT)
                curve(srg, list(reversed(pts)))
        for k, (i, pts) in enumerate(d.gene_edges[rid]):
            g = d.genes[i]["gene"]
            if g not in model.rxns[rid]["genes"]:
                continue
            srg = rg.createSpeciesReferenceGlyph()
            srg.setId(f"mrg_{sid(rid)}_{k}")
            srg.setSpeciesGlyphId(f"gg{i + 1}")
            srg.setSpeciesReferenceId(sid(f"{rid}_gene_{g}"))
            srg.setRole(libsbml.SPECIES_ROLE_MODIFIER)
            curve(srg, pts)
    libsbml.writeSBMLToFile(doc, path)
    back = libsbml.readSBMLFromFile(path)  # validate what was written
    back.checkConsistency()
    errors = [back.getError(i) for i in range(back.getNumErrors())
              if back.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
    return [f"{e.getErrorId()}: {e.getShortMessage()}" for e in errors]


SCHEMA = None


def validate_sbgn(path):
    """Errors of an SBGN-ML file against the libsbgn schema (schema/SBGN.xsd next to this file)."""
    global SCHEMA
    xsd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "schema", "SBGN.xsd")
    if not os.path.exists(xsd):
        return []
    if SCHEMA is None:
        SCHEMA = etree.XMLSchema(etree.parse(xsd))
    ok = SCHEMA.validate(etree.parse(path))
    return [] if ok else [str(e) for e in SCHEMA.error_log][:3]


# ----------------------------------------------------------------------------- Escher

ESCHER_SCHEMA = "https://escher.github.io/escher/jsonschema/1-0-0#"

# Escher's own styles: labels are bold sans-serif drawn from their left end on their baseline, 20 px for
# metabolites, 30 px for reactions (the gene rule, when shown, follows below in 18 px lines) and 50 px for text
# labels; metabolite circles have a radius of 20 px (primary) or 10 px (secondary)
NODE_SIZE, REACTION_SIZE, TEXT_SIZE = 20, 30, 50
RADIUS = {True: 20, False: 10}
# Helvetica Bold advance widths (1/1000 em) of the printable ASCII characters; any other character counts 600
BOLD = dict(zip(map(chr, range(32, 127)), [
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556, 556, 556, 556,
    556, 556, 556, 556, 333, 333, 584, 584, 584, 611, 975, 722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722,
    611, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584, 556, 333, 556,
    611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, 611, 611, 611, 611, 389, 556, 333, 611, 556, 778,
    556, 556, 500, 389, 280, 389, 584]))
REPEAT = 3500   # px between the repeated names of a large compartment
LARGE = 2500    # px: a compartment box this wide or high gets its name at its corners too


def text_width(s, size):
    return size * sum(BOLD.get(ch, 600) for ch in s) / 1000


def text_label_box(x, y, w, size):
    """Box of a label drawn from (x, y), its left end on the baseline."""
    return (x, y - 0.75 * size, x + w, y + 0.25 * size)


class Taken:
    """What an Escher map draws (nodes, segments, labels) as weighted boxes, to find free places for labels."""

    def __init__(self, cell=100):
        self.cell, self.grid = cell, collections.defaultdict(list)

    def cells(self, b):
        c = self.cell
        return [(i, j) for i in range(int(b[0] // c), int(b[2] // c) + 1)
                for j in range(int(b[1] // c), int(b[3] // c) + 1)]

    def add(self, b, weight):
        item = (b, weight)
        for key in self.cells(b):
            self.grid[key].append(item)

    def cost(self, b):
        seen, total = set(), 0
        for key in self.cells(b):
            for item in self.grid.get(key, ()):
                o, w = item
                if id(item) not in seen and o[0] < b[2] and o[2] > b[0] and o[1] < b[3] and o[3] > b[1]:
                    seen.add(id(item))
                    total += w
        return total


def segment_points(a, b, b1, b2, step=15):
    """Points along a segment as Escher draws it: a cubic Bezier through its control points, or a line."""
    if b1 is None or b2 is None:
        ctrl = [a, b]
    else:
        ctrl = [a, (b1["x"], b1["y"]), (b2["x"], b2["y"]), b]
    k = max(2, int(L.poly_length(ctrl) / step))
    out = []
    for i in range(k + 1):
        t = i / k
        if len(ctrl) == 2:
            out.append(L.add(a, L.sub(b, a), t))
        else:
            c = [(1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3]
            out.append((sum(w * p[0] for w, p in zip(c, ctrl)), sum(w * p[1] for w, p in zip(c, ctrl))))
    return out


def best(taken, cands, w, size):
    """The candidate (preference, x, y) whose label box overlaps least of what is drawn, then the preferred."""
    return min(cands, key=lambda c: (taken.cost(text_label_box(c[1], c[2], w, size)), c[0]))


def escher_text(text, cx, base, size):
    """An SVG label block (centre x, mean baseline, font size) as an Escher text label at the same centre."""
    return cx - text_width(text, TEXT_SIZE) / 2, base - 0.35 * size + 0.35 * TEXT_SIZE


def compartment_spots(g, box, head, w):
    """Places to repeat a compartment's name, with the directions they may slide to find room: the corners and
    along the edges of a large box, or along a separator line on the heading's side."""
    asc, desc, m = 0.75 * TEXT_SIZE, 0.25 * TEXT_SIZE, 70
    out = []
    if box:
        x0, y0, x1, y1 = box
        if max(x1 - x0, y1 - y0) < LARGE:
            return out
        left, right = x0 + m, max(x0 + m, x1 - m - w)
        n = max(1, math.ceil((right - left) / REPEAT))
        for k in range(n + 1):
            x = left + (right - left) * k / n
            out += [((x, y0 + m + asc), (1, 0), (0, 1)), ((x, y1 - m - desc), (1, 0), (0, -1))]
        top, bottom = y0 + m + asc, y1 - m - desc
        n = math.ceil((bottom - top) / REPEAT)
        for k in range(1, n):
            y = top + (bottom - top) * k / n
            out += [((left, y), (0, 1), (1, 0)), ((right, y), (0, 1), (-1, 0))]
        return out
    lines = [sp for path in g.iter(L.SVG + "path") for sp in (L.subpaths(path.get("d") or "") or [])
             if len(sp) >= 2 and abs(sp[0][0] - sp[-1][0]) < 5]
    if not lines:
        return out
    xs = sorted(sp[0][0] for sp in lines)
    hx = head[0] + w / 2
    x = (xs[0] + xs[-1] - w) / 2 if xs[0] < hx < xs[-1] else xs[0] - 60 - w if hx < xs[0] else xs[-1] + 60
    bottom = max(max(p[1] for p in sp) for sp in lines) - 200
    y = head[1] + REPEAT
    while y < bottom:
        out.append(((x, y), (0, 1), (0, 0)))
        y += REPEAT
    return out


def write_escher(d, path, model_name):
    """An Escher map (JSON, schema 1-0-0): metabolite nodes where the SVG draws them (main metabolites as
    primary nodes, cofactors as secondary ones), each reaction as a midmarker with two multimarkers along its
    direction, the SVG's edge paths as Bezier segments (their bends as control points), and the title, band names
    and compartment names as text labels. Escher has no boxes, bands or gene boxes, so a large compartment's name
    is repeated at its corners and along its edges or separator line; genes are in the gene rules.

    Labels are placed for names (Escher's "Descriptive names"). A cofactor's label stands where the SVG has it;
    a main metabolite's label is centred like the SVG's, above or below the node (or beside it) where it
    overlaps least of the drawing; a reaction's label starts where the SVG draws its gene boxes, else beside the
    reaction, wherever its name overlaps least."""
    import json
    model = d.model
    ids = iter(range(1, 10 ** 9))
    nodes, reactions, labels = {}, {}, {}
    met_node = {}
    for i, x in enumerate(d.mets):
        nid = str(next(ids))
        met_node[i] = nid
        px, py = x["pos"]
        nodes[nid] = {"node_type": "metabolite", "x": round(px, 1), "y": round(py, 1), "bigg_id": x["met"],
                      "name": model.name(x["met"]), "label_x": round(px, 1), "label_y": round(py, 1),
                      "node_is_primary": bool(x["main"])}
    axes = {}
    for rid, rp in d.rpos.items():
        r = model.rxns[rid]
        st = r["stoich"]
        ins = [pts[-2] for i, pts in d.edges[rid] if st.get(d.mets[i]["met"], 0) < 0 and len(pts) > 1]
        outs = [pts[-2] for i, pts in d.edges[rid] if st.get(d.mets[i]["met"], 0) > 0 and len(pts) > 1]
        a = (sum(p[0] for p in ins) / len(ins), sum(p[1] for p in ins) / len(ins)) if ins else (rp[0] - 1, rp[1])
        b = (sum(p[0] for p in outs) / len(outs), sum(p[1] for p in outs) / len(outs)) if outs else (rp[0] + 1, rp[1])
        u = L.unit(L.sub(b, a)) if L.length(L.sub(b, a)) > 1 else (1.0, 0.0)
        mid, m_in, m_out = str(next(ids)), str(next(ids)), str(next(ids))
        for nid, q in ((mid, rp), (m_in, L.add(rp, u, -20)), (m_out, L.add(rp, u, 20))):
            nodes[nid] = {"node_type": "midmarker" if nid == mid else "multimarker", "x": round(q[0], 1),
                          "y": round(q[1], 1)}
        segments = {str(next(ids)): {"from_node_id": m_in, "to_node_id": mid, "b1": None, "b2": None},
                    str(next(ids)): {"from_node_id": mid, "to_node_id": m_out, "b1": None, "b2": None}}
        for i, pts in d.edges[rid]:
            role = st.get(d.mets[i]["met"], 0)
            if not role:
                continue
            bends = [{"x": round(q[0], 1), "y": round(q[1], 1)} for q in pts[1:-1]]
            b1, b2 = (bends[0], bends[-1]) if bends else (None, None)
            if role < 0:  # metabolite to the multimarker on the substrate side
                segments[str(next(ids))] = {"from_node_id": met_node[i], "to_node_id": m_in, "b1": b1, "b2": b2}
            else:
                segments[str(next(ids))] = {"from_node_id": m_out, "to_node_id": met_node[i], "b1": b2, "b2": b1}
        rule = mapedit.text(r.get("gene_reaction_rule"))
        genes = [{"bigg_id": g, "name": model.symbol.get(g) or g} for g in sorted(r["genes"])]
        reactions[str(next(ids))] = {
            "name": mapedit.text(r.get("name")) or rid, "bigg_id": rid,
            "reversibility": (r.get("lower_bound") or 0) < 0,
            "label_x": round(rp[0], 1), "label_y": round(rp[1], 1),
            "gene_reaction_rule": rule, "genes": genes,
            "metabolites": [{"bigg_id": m, "coefficient": c} for m, c in st.items()],
            "segments": segments}
        axes[rid] = u

    # what is drawn: circles, segments, then each label as it is placed
    taken = Taken()
    for n in nodes.values():
        rad = RADIUS[n["node_is_primary"]] if n["node_type"] == "metabolite" else 5
        taken.add((n["x"] - rad, n["y"] - rad, n["x"] + rad, n["y"] + rad), 5)
    for r in reactions.values():
        for seg in r["segments"].values():
            fa, fb = nodes[seg["from_node_id"]], nodes[seg["to_node_id"]]
            for q in segment_points((fa["x"], fa["y"]), (fb["x"], fb["y"]), seg["b1"], seg["b2"]):
                taken.add((q[0] - 4, q[1] - 4, q[0] + 4, q[1] + 4), 1)

    def place(node, lx, ly):
        node["label_x"], node["label_y"] = round(lx, 1), round(ly, 1)
        taken.add(text_label_box(lx, ly, text_width(node["name"], NODE_SIZE), NODE_SIZE), 5)

    for primary in (False, True):  # cofactors first: their labels stay where the SVG has them
        for i, x in enumerate(d.mets):
            if bool(x["main"]) != primary:
                continue
            node = nodes[met_node[i]]
            (px, py), (anchor, ax, base, size) = x["pos"], x["label"]
            w = text_width(node["name"], NODE_SIZE)
            if not primary:
                lx = px + ax if anchor == "start" else px + ax - w if anchor == "end" else px + ax - w / 2
                place(node, lx, py + base * NODE_SIZE / size)
                continue
            cx = px + ax if anchor == "middle" else px
            rad = RADIUS[True]
            cands = [(0, cx - w / 2, py - rad - 8), (1, cx - w / 2, py + rad + 6 + 0.75 * NODE_SIZE),
                     (2, px + rad + 7, py + 7), (3, px - rad - 7 - w, py + 7)]
            _, lx, ly = best(taken, cands, w, NODE_SIZE)
            place(node, lx, ly)

    gene_boxes = collections.defaultdict(list)
    for g in d.genes:
        for rid in g["rids"]:
            gene_boxes[rid].append(g["box"])
    for r in reactions.values():  # from their left end, so that an id and a long name both grow away
        rid = r["bigg_id"]
        rp, u = d.rpos[rid], axes[rid]
        w_id = text_width(rid, REACTION_SIZE)
        w = max(w_id, text_width(r["name"], REACTION_SIZE))
        cands = []
        gb = gene_boxes.get(rid)
        if gb:  # where the SVG draws the gene boxes, which Escher does not show
            x0, y0 = min(b[0] for b in gb), min(b[1] for b in gb)
            x1, y1 = max(b[0] + b[2] for b in gb), max(b[1] + b[3] for b in gb)
            if L.length(L.sub(((x0 + x1) / 2, (y0 + y1) / 2), rp)) < 400:
                cands.append((0, x0, y0 + 0.75 * REACTION_SIZE))
        n = L.normal(u)
        sides = [d.mets[i]["pos"] for i, _ in d.edges[rid] if not d.mets[i]["main"]]
        lean = sum((p[0] - rp[0]) * n[0] + (p[1] - rp[1]) * n[1] for p in sides)
        away = (-n[0], -n[1]) if lean > 0 else n  # the side without cofactors
        for j, gap in enumerate((0, 25, 55, 90)):
            spots = {"above": (rp[0] - w_id / 2, rp[1] - 38 - gap),
                     "below": (rp[0] - w_id / 2, rp[1] + 38 + 0.75 * REACTION_SIZE + gap),
                     "right": (rp[0] + 38 + gap, rp[1] + 0.25 * REACTION_SIZE),
                     "left": (rp[0] - 38 - gap - w_id, rp[1] + 0.25 * REACTION_SIZE)}
            facing = {"above": (0, -1), "below": (0, 1), "right": (1, 0), "left": (-1, 0)}
            for side, (x, y) in spots.items():
                f = facing[side]
                pref = 1 + j * 0.3 + (0 if f[0] * away[0] + f[1] * away[1] > 0.5 else 0.5) + (0.3 if side == "left" else 0)
                cands.append((pref, x, y))
        _, lx, ly = best(taken, cands, w, REACTION_SIZE)
        r["label_x"], r["label_y"] = round(lx, 1), round(ly, 1)
        taken.add(text_label_box(lx, ly, w, REACTION_SIZE), 5)

    def text_label(text, x, y):
        labels[str(next(ids))] = {"text": text, "x": round(x, 1), "y": round(y, 1)}
        taken.add(text_label_box(x, y, text_width(text, TEXT_SIZE), TEXT_SIZE), 5)

    groups = [g for g in d.mp.root.iter(L.SVG + "g") if g.get("class") in ("compartment", "subsystem")]
    heads = []
    for g in groups:
        for k, (text, cx, base, size) in enumerate(text_blocks(g)):
            x, y = escher_text(text, cx, base, size)
            text_label(text, x, y)
            if k == 0 and g.get("class") == "compartment":
                heads.append((g, text, (x, y)))
    boxes = L.compartment_boxes(d.mp)
    for g, text, head in heads:  # the name again where a large compartment reaches beyond its heading
        w = text_width(text, TEXT_SIZE)
        placed = [head]
        for (x, y), along, inward in compartment_spots(g, boxes.get(g.get("id")), head, w):
            if any(L.length(L.sub((x, y), p)) < REPEAT * 0.4 for p in placed):
                continue
            for s, t in [(s, t) for t in (0, 150, 300) for s in (0, 150, -150, 300, -300, 600, -600, 900, -900)]:
                q = (x + along[0] * s + inward[0] * t, y + along[1] * s + inward[1] * t)
                if taken.cost(text_label_box(q[0], q[1], w, TEXT_SIZE)) == 0:
                    text_label(text, *q)
                    placed.append(q)
                    break

    header = {"map_name": f"{d.title} ({model_name})", "map_id": d.name,
              "map_description": f"{d.title}: a Metabolic Atlas map of {model_name}, converted from its SVG map. "
                                 "Its labels are placed for descriptive names.",
              "homepage": "https://escher.github.io", "schema": ESCHER_SCHEMA}
    body = {"reactions": reactions, "nodes": nodes, "text_labels": labels,
            "canvas": {"x": 0, "y": 0, "width": d.width, "height": d.height}}
    with open(path, "w") as fh:
        json.dump([header, body], fh, separators=(",", ":"))


ESCHER = None


def validate_escher(path):
    """Errors of an Escher map against its JSON schema (schema/escher_1-0-0.json next to this file) and Escher's
    own checks (escher.validate): segments end on nodes, connected metabolites have a coefficient, genes have
    names."""
    import json
    global ESCHER
    errors = []
    data = json.load(open(path))
    js = os.path.join(os.path.dirname(os.path.abspath(__file__)), "schema", "escher_1-0-0.json")
    if os.path.exists(js):
        import jsonschema
        if ESCHER is None:
            ESCHER = json.load(open(js))
        try:
            jsonschema.validate(data, ESCHER)
        except jsonschema.ValidationError as err:
            errors.append(f"schema: {err.message[:200]}")
    nodes = data[1]["nodes"]
    for r in data[1]["reactions"].values():
        coef = {m["bigg_id"] for m in r["metabolites"] if m["coefficient"]}
        names = {g["bigg_id"] for g in r["genes"] if g.get("name")}
        for sid_, seg in r["segments"].items():
            for k in ("from_node_id", "to_node_id"):
                n = nodes.get(seg[k])
                if n is None:
                    errors.append(f"{r['bigg_id']}: segment {sid_} has no {k} node")
                elif n["node_type"] == "metabolite" and n["bigg_id"] not in coef:
                    errors.append(f"{r['bigg_id']}: no coefficient for {n['bigg_id']}")
        rule = re.sub(r"([()\s])(?:and|or)([)(\s])", r"\1\2", f" {r['gene_reaction_rule']} ")
        for g in re.sub(r"[()]", " ", rule).split():
            if g not in names:
                errors.append(f"{r['bigg_id']}: no name for gene {g}")
    return errors[:3]


def write_png(svg, path, width):
    import cairosvg
    from PIL import Image
    png = cairosvg.svg2png(url=svg, output_width=width, background_color="white")
    Image.open(io.BytesIO(png)).convert("RGB").save(path, optimize=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model_yml")
    ap.add_argument("model_dir")
    ap.add_argument("out")
    ap.add_argument("--maps", nargs="+", required=True)
    ap.add_argument("--png-width", type=int, default=3000)
    ap.add_argument("--model-name", help="model name for the Escher maps (default: the YAML file name)")
    a = ap.parse_args()
    model = mapedit.Model(a.model_yml, a.model_dir, [])
    for fmt in ("svg", "sbgn", "sbml", "png", "escher"):
        os.makedirs(os.path.join(a.out, fmt), exist_ok=True)
    model_name = f"{a.model_name or os.path.basename(a.model_yml).rsplit('.', 1)[0]} {model.version}"
    problems = 0
    for p in a.maps:
        d = MapData(p, model)
        shutil.copyfile(p, os.path.join(a.out, "svg", d.name + ".svg"))
        write_sbgn(d, os.path.join(a.out, "sbgn", d.name + ".sbgn"))
        errors = validate_sbgn(os.path.join(a.out, "sbgn", d.name + ".sbgn"))
        errors += write_sbml(d, os.path.join(a.out, "sbml", d.name + ".sbml"))
        write_png(p, os.path.join(a.out, "png", d.name + ".png"), a.png_width)
        write_escher(d, os.path.join(a.out, "escher", d.name + ".json"), model_name)
        errors += validate_escher(os.path.join(a.out, "escher", d.name + ".json"))
        problems += len(errors)
        edges = sum(len(v) for v in d.edges.values())
        print(f"{d.name}: {len(d.rpos)} reactions, {len(d.mets)} metabolite nodes, {len(d.genes)} gene nodes, "
              f"{edges} edges" + (f"; errors: {errors[:3]}" if errors else ""), flush=True)
    print("SBGN schema, SBML consistency and Escher errors:", problems)


if __name__ == "__main__":
    main()
