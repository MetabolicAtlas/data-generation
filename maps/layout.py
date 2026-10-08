"""Layout rules for the map editor: keep pathway shapes, restore connections, add missing reactions.

D1  keep a pathway's shape: when a removed reaction ran between two main metabolites and the
    model still connects them through reactions that are already drawn elsewhere, those
    reactions and their intermediates are moved onto the removed reaction's line.
D2  restore connections: when the model connects them through reactions that are not drawn
    (preferring reactions of the map's subsystem), those reactions are drawn on that line.
D3  every reaction of the subsystem is on its map: a missing reaction is attached next to a
    drawn main metabolite, or placed as a small module in free space of its compartment.
New and moved elements copy the maps' style: diamond reaction nodes, blue ellipses for main
metabolites, dots for cofactors, yellow gene boxes, arrowheads on products (and substrates
of reversible reactions), and the background band along main edges.
"""

import collections
import heapq
import math
import re

from lxml import etree

from mapedit import CURRENCY, SVG, comp_of, is_rxn, segments, subpaths, to_d, translate

COFACTORS = CURRENCY | {"thiamin-PP", "lipoamide", "dihydrolipoamide", "FAD", "FADH2", "ubiquinone", "ubiquinol",
                        "AMP", "PPi", "NH4+", "glutathione", "GSSG", "SAM", "SAH", "UDP-glucuronate",
                        "PAPS", "PAP", "ferricytochrome C", "ferrocytochrome C", "glycine", "taurine",
                        # Yeast-GEM names
                        "S-adenosyl-L-methionine", "S-adenosyl-L-homocysteine", "thiamine diphosphate",
                        "3'-phosphoadenylyl sulfate", "adenosine 3',5'-bismethylphosphate", "FMN", "reduced FMN"}
COMPARTMENT = {"c": "Cytosol", "m": "Mitochondria", "x": "Peroxisome", "r": "Endoplasmic reticulum",
               "g": "Golgi apparatus", "l": "Lysosome", "n": "Nucleus", "e": "Extracellular", "i": "Inner mitochondria"}
STEP = 150  # distance between a reaction node and a main metabolite node in new drawings
ROW_STEP = 210  # the same in the bottom section, where long names sit side by side


# ----------------------------------------------------------------------------- geometry

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def add(a, b, k=1.0):
    return (a[0] + k * b[0], a[1] + k * b[1])


def length(v):
    return math.hypot(v[0], v[1])


def unit(v):
    n = length(v) or 1.0
    return (v[0] / n, v[1] / n)


def normal(v):
    return (-v[1], v[0])


def poly_length(pts):
    return sum(length(sub(b, a)) for a, b in zip(pts, pts[1:]))


def point_at(pts, d):
    """Point and unit tangent at arc length d along a polyline."""
    for a, b in zip(pts, pts[1:]):
        seg = length(sub(b, a))
        if d <= seg or b is pts[-1]:
            t = 0 if seg == 0 else min(1, d / seg)
            return add(a, sub(b, a), t), unit(sub(b, a))
        d -= seg
    return pts[-1], unit(sub(pts[-1], pts[-2]))


def cut(pts, d0, d1):
    """The part of a polyline between arc lengths d0 and d1."""
    out, acc = [point_at(pts, d0)[0]], 0.0
    for a, b in zip(pts, pts[1:]):
        acc += length(sub(b, a))
        if d0 < acc < d1:
            out.append(b)
    out.append(point_at(pts, d1)[0])
    return out


def boundary(center, toward, kind):
    """Point where an edge from `center` toward `toward` leaves a node of the given kind."""
    u = unit(sub(toward, center))
    if u == (0.0, 0.0):
        return center
    if kind == "rea":
        r = 30 / ((abs(u[0]) + abs(u[1])) or 1)
    elif kind == "main":
        r = 1 / math.sqrt((u[0] / 32) ** 2 + (u[1] / 22) ** 2)
    else:
        r = 10
    return add(center, u, r)


def arrow(tip, frm):
    u = unit(sub(tip, frm))
    base = add(tip, u, -13)
    n = normal(u)
    a, b = add(base, n, 6), add(base, n, -6)
    return [tip, a, b, tip]


def label_width(s, size):
    return 0.62 * size * len(s)


# ----------------------------------------------------------------------------- elements

def fmt(v):
    return ("%.1f" % v).rstrip("0").rstrip(".")


def matrix(x, y):
    return "matrix(1,0,0,1,%s,%s)" % (fmt(x), fmt(y))


def wrap(name, width=16, lines=3):
    """Split a name into at most `lines` lines of about `width` characters, at spaces, hyphens or commas."""
    out = []
    rest = name
    while rest and len(out) < lines - 1 and len(rest) > width:
        cut_at = max((i for i, ch in enumerate(rest[: width + 1]) if ch in " -,"), default=width - 1)
        out.append(rest[: cut_at + 1].rstrip())
        rest = rest[cut_at + 1:].lstrip()
    if rest:
        out.append(rest)
    return out


def element(tag, **attrs):
    el = etree.Element(SVG + tag)
    for k, v in attrs.items():
        el.set(k.replace("_", "-"), v)
    return el


def text_el(content, **attrs):
    t = element("text", **attrs)
    t.text = "\n        %s\n      " % content
    return t


def new_reaction(rid, x, y):
    g = element("g", id=rid, **{"class": "rea"})
    g.append(element("path", **{"class": "shape", "d": "M0-30L30 0 0 30-30 0 0-30", "stroke": "#000", "fill": "#fff",
                                "transform": matrix(x, y)}))
    g.append(text_el(rid, **{"class": "lbl", "font-family": "sans-serif", "font-size": "16", "font-weight": "700",
                             "text-anchor": "middle", "y": "5", "transform": matrix(x, y)}))
    return g


def new_main(mid, name, rids, x, y):
    g = element("g", **{"class": " ".join(["met", mid] + rids), "transform": matrix(x, y)})
    g.append(element("ellipse", **{"class": "shape", "rx": "30", "ry": "20", "fill": "#9df", "stroke": "#333"}))
    lines = wrap(name)
    if len(lines) == 1:
        g.append(text_el(lines[0], **{"class": "lbl", "font-family": "sans-serif", "font-size": "18",
                                      "font-weight": "700", "text-anchor": "middle", "y": "6"}))
    else:
        lb = element("g", **{"class": "lbl", "font-family": "sans-serif", "font-size": "18", "font-weight": "700",
                             "text-anchor": "middle"})
        top = 6 - 10 * (len(lines) - 1)
        for k, line in enumerate(lines):
            lb.append(text_el(line, y=fmt(top + 20 * k)))
        g.append(lb)
    return g


def new_side(mid, name, rids, x, y, left):
    g = element("g", **{"class": " ".join(["met", mid] + rids), "transform": matrix(x, y)})
    g.append(element("circle", **{"class": "shape", "r": "7", "fill": "#222"}))
    g.append(text_el(name, **{"class": "lbl", "font-family": "sans-serif", "font-size": "18", "font-weight": "700",
                              "text-anchor": "end" if left else "start", "x": "-12" if left else "12", "y": "6"}))
    return g


GENE_LABEL = "name"  # gene boxes show the gene name ("name"), the gene id ("orf"), or the name over the id ("both")


def gene_label(gene, symbol):
    """The text a gene box shows, as lines."""
    if GENE_LABEL == "orf" or not symbol:
        return [gene]
    if GENE_LABEL == "both" and symbol != gene:
        return [symbol, gene]
    return [symbol]


def set_gene_text(t, gene, symbol):
    """Write the gene label into a gene box's text element (one line, or the name over the id)."""
    lines = gene_label(gene, symbol)
    for c in list(t):
        t.remove(c)
    if len(lines) == 1:
        t.text = lines[0]
        t.set("font-size", "12")
        t.set("y", "20")
        return
    t.text = None
    t.set("font-size", "10")
    t.set("y", "13")
    a = etree.SubElement(t, SVG + "tspan", x=t.get("x") or "30", y="13")
    a.text = lines[0]
    b = etree.SubElement(t, SVG + "tspan", x=t.get("x") or "30", y="25", **{"font-size": "8", "font-weight": "400"})
    b.text = lines[1]


def new_gene(gene, symbol, rids, x, y):
    g = element("g", **{"class": " ".join(["enz", gene] + rids)})
    g.append(element("path", **{"class": "shape", "fill": "#feb", "stroke": "#000", "d": "M0 0h60v30H0z",
                                "transform": matrix(x, y)}))
    t = text_el("", **{"class": "lbl", "font-family": "sans-serif", "font-size": "12", "font-weight": "700",
                       "text-anchor": "middle", "x": "30", "y": "20", "transform": matrix(x, y)})
    set_gene_text(t, gene, symbol)
    g.append(t)
    return g


def set_position(el, x, y):
    if el.get("transform"):
        el.set("transform", matrix(x, y))
    for child in el.iter():
        if child is not el and child.get("transform") and "matrix" in child.get("transform"):
            child.set("transform", matrix(x, y))


def append(layer, el):
    last = layer[-1] if len(layer) else None
    if last is not None:
        el.tail = last.tail
        last.tail = "\n    "
    else:
        layer.text = "\n    "
        el.tail = "\n  "
    layer.append(el)


# ----------------------------------------------------------------------------- free space

class Space:
    """Bounding boxes of what is drawn, in a grid index, to find free places for new elements."""
    CELL = 200

    def __init__(self, mp):
        self.grid = collections.defaultdict(list)
        self.recent = []  # boxes added since an Occupancy grid last took them in
        for n in mp.layer["nodes"]:
            b = node_box(mp, n)
            if b:
                self.add(b)
        self.keep = list(mp.root.iter(SVG + "text"))  # an element's id() stays the same while it is referred to
        in_nodes = {id(t) for n in mp.layer["nodes"] for t in n.iter(SVG + "text")}
        for t in self.keep:  # headings, titles and any other free text
            if id(t) not in in_nodes:
                b = text_box(t)
                if b:
                    self.add(b)
        for layer in ("fes", "ees"):
            for p in mp.layer[layer]:
                self.add_edges(subpaths(p.get("d") or "") or [])

    def cells(self, box):
        c = self.CELL
        for i in range(int(box[0] // c), int(box[2] // c) + 1):
            for j in range(int(box[1] // c), int(box[3] // c) + 1):
                yield (i, j)

    def add(self, box):
        if box:
            for cell in self.cells(box):
                self.grid[cell].append(box)
            self.recent.append(box)

    def add_edges(self, pts_list):
        for a, b in segments(pts_list):
            k = max(1, int(length(sub(b, a)) / 25))
            for i in range(k + 1):
                q = add(a, sub(b, a), i / k)
                self.add((q[0], q[1], q[0], q[1]))

    def free(self, box, margin=12):
        x0, y0, x1, y1 = box[0] - margin, box[1] - margin, box[2] + margin, box[3] + margin
        for cell in self.cells((x0, y0, x1, y1)):
            for b in self.grid.get(cell, ()):
                if b[0] <= x1 and b[2] >= x0 and b[1] <= y1 and b[3] >= y0:
                    return False
        return True


def text_box(t):
    """Approximate box of a positioned text element (titles and compartment labels)."""
    content = "".join(t.itertext()).strip()
    if not content:
        return None
    size = float(t.get("font-size") or (t.getparent().get("font-size") if t.getparent() is not None else None) or 40)
    pos = None
    for el in (t, t.getparent()):
        if el is not None and el.get("transform"):
            m = re.search(r"matrix\(([^)]*)\)", el.get("transform"))
            if m:
                a = [float(v) for v in m.group(1).split(",")]
                pos = (a[4], a[5])
                break
    if pos is None:
        try:
            pos = (float(t.get("x") or 0), float(t.get("y") or 0))
        except ValueError:
            return None
    y = pos[1] + float(t.get("y") or 0) if t.get("transform") else pos[1]
    w = label_width(content, size)
    anchor = t.get("text-anchor") or (t.getparent().get("text-anchor") if t.getparent() is not None else None)
    x0 = pos[0] - w / 2 if anchor == "middle" else pos[0] - w if anchor == "end" else pos[0]
    return (x0 - 10, y - size, x0 + w + 10, y + size * 0.4)


def node_box(mp, n):
    cls = mp.classes(n)
    if not cls:
        return None
    if cls[0] == "rea":
        p = translate(n)
        return (p[0] - 32, p[1] - 32, p[0] + 32, p[1] + 32) if p else None
    if cls[0] == "enz":
        p = translate(n)
        return (p[0], p[1], p[0] + 60, p[1] + 30) if p else None
    if cls[0] == "met":
        p = translate(n)
        if not p:
            return None
        texts = ["".join(t.itertext()).strip() for t in n.iter(SVG + "text")]
        w = max([label_width(t, 18) for t in texts] + [60])
        if n.find(SVG + "ellipse") is not None:
            extra = 10 * (len(texts) - 1)
            return (p[0] - w / 2, p[1] - 24 - extra, p[0] + w / 2, p[1] + 24 + extra)
        t = n.find(SVG + "text")
        left = t is not None and t.get("text-anchor") == "end"
        return (p[0] - (w + 12 if left else 8), p[1] - 12, p[0] + (8 if left else w + 12), p[1] + 12)
    return None


def compartment_boxes(mp):
    out = {}
    for g in mp.root.iter(SVG + "g"):
        if g.get("class") == "compartment":
            r = g.find(SVG + "rect")
            if r is not None and r.get("width"):
                x, y = float(r.get("x")), float(r.get("y"))
                out[g.get("id")] = (x, y, x + float(r.get("width")), y + float(r.get("height")))
    return out


def canvas(mp):
    return (0.0, 0.0, float(mp.root.get("width")), float(mp.root.get("height")))


ALIASES = {"Inner mitochondria": "Inner membrane"}


def region(mp, comp_name):
    """Rectangle of a compartment on the map. A compartment without its own box (the cytosol on subsystem
    maps, the map's own compartment on compartment maps) is the canvas outside the other boxes."""
    boxes = compartment_boxes(mp)
    name = comp_name if comp_name in boxes else ALIASES.get(comp_name, comp_name)
    if name in boxes:
        return boxes[name], []
    reg = canvas(mp)
    # a compartment without a box starts below its heading
    for g in mp.root.iter(SVG + "g"):
        if g.get("class") == "compartment" and g.get("id") == comp_name and g.find(SVG + "rect") is None:
            heads = [text_box(t) for t in g.iter(SVG + "text")]
            heads = [h for h in heads if h]
            if heads:
                reg = (reg[0], min(h[1] for h in heads), reg[2], reg[3])
    return reg, list(boxes.values())


def content_center(mp, reg, holes):
    pts = []
    for n in mp.layer["nodes"]:
        p = translate(n)
        if p and reg[0] <= p[0] <= reg[2] and reg[1] <= p[1] <= reg[3] and not any(
                h[0] <= p[0] <= h[2] and h[1] <= p[1] <= h[3] for h in holes):
            pts.append(p)
    if not pts:
        return None
    return (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))


def find_spot(space, reg, holes, w, h, near=None, step=40):
    """Top-left corner of a free w x h box inside reg (and outside holes), closest to `near`."""
    x0, y0, x1, y1 = reg
    cands = []
    y = y0 + 80
    while y + h < y1 - 150:
        x = x0 + 60
        while x + w < x1 - 60:
            box = (x, y, x + w, y + h)
            if not any(box[0] < hb[2] + 30 and box[2] > hb[0] - 30 and box[1] < hb[3] + 30 and box[3] > hb[1] - 30
                       for hb in holes):
                d = 0 if near is None else math.hypot(x + w / 2 - near[0], y + h / 2 - near[1])
                cands.append((d, x, y))
            x += step
        y += step
    # modules go in a row along the bottom of the compartment: lowest free row first, then left to right
    for _, x, y in sorted(cands, key=lambda t: (-t[2], t[1])):
        if space.free((x, y, x + w, y + h), margin=20):
            return x, y
    return None


# ----------------------------------------------------------------------------- drawing reactions

class Drawer:
    def __init__(self, mp, model, log):
        self.mp, self.model, self.log = mp, model, log
        self.space = Space(mp)
        self.band_add = []   # polylines to add to the background band
        self.gone = []       # edge segments removed (for band trimming)
        self.placed = set()  # ids of reactions placed or moved in this run
        self.lines = []      # lines that paths were laid on (D1/D2), checked for rings (D4)
        self.pending = []    # (reaction, compartment) for placement (D9)
        self.context = set()  # reactions of other subsystems added as context (D14): drawn without the band
        self.comps_on_map = set()  # compartments the map has an area for, its own compartment included
        self.one_area = ""   # name of the single area all compartments are drawn in (maps without boxes)
        self.kegg = []       # KEGG pathway layouts for this map (D7)
        self.hubs = "spread" # how crowded KEGG modules open up: "spread" (move close compounds) or "scale"

    # nodes -------------------------------------------------------------------
    def met_nodes(self, mid):
        return [n for n in self.mp.nodes("met") if self.mp.classes(n)[1] == mid]

    def is_main(self, n):
        return n.find(SVG + "ellipse") is not None

    def link(self, n, rid):
        cls = self.mp.classes(n)
        if rid not in cls:
            n.set("class", " ".join(cls + [rid]))

    def unlink(self, n, rid):
        cls = [c for c in self.mp.classes(n) if c != rid]
        n.set("class", " ".join(cls))
        if not [c for c in cls[2:] if is_rxn(c)]:
            self.mp.remove(n)

    def side_name(self, m):
        return self.model.name(m)

    # edges -------------------------------------------------------------------
    def remove_edges(self, rid):
        for layer in ("fes", "fehs", "ees"):
            for p in self.mp.edges(layer, rid):
                if layer == "fes":
                    self.gone += segments(subpaths(p.get("d")) or [])
                self.mp.remove(p)

    def draw(self, rid, rpos, links, genes_side):
        """Edges, arrowheads and genes of reaction rid at rpos.

        links: (node position, node kind 'main'|'side', role -1|+1, polyline from node to reaction or None)."""
        rev = (self.model.rxns[rid].get("lower_bound") or 0) < 0
        fe, feh = [], []
        for npos, kind, role, line in links:
            pts = list(line) if line else [npos, rpos]
            pts[0] = boundary(npos, pts[1], kind)
            pts[-1] = boundary(rpos, pts[-2], "rea")
            fe.append(pts)
            if role > 0 or rev:
                feh.append(arrow(pts[0], pts[1]))
            if kind == "main":
                self.band_add.append((rid, [npos] + pts[1:-1] + [rpos]))
        append(self.mp.layer["fes"], element("path", **{"class": "fe " + rid, "d": to_d(fe)}))
        append(self.mp.layer["fehs"], element("path", **{"class": "feh " + rid, "d": to_d(feh)}))
        self.space.add_edges(fe)
        self.space.add((rpos[0] - 32, rpos[1] - 32, rpos[0] + 32, rpos[1] + 32))
        self.place_genes(rid, rpos, genes_side)

    def place_genes(self, rid, rpos, side):
        """Gene boxes of rid in a stack beside the reaction (side: unit vector), with an edge to it."""
        genes = sorted(self.model.rxns[rid]["genes"])
        old = [n for n in self.mp.nodes("enz") if rid in self.mp.classes(n)]
        for n in old:
            self.unlink(n, rid)
        for p in self.mp.edges("ees", rid):
            self.mp.remove(p)
        if not genes:
            return
        cols = min(3, len(genes))
        rows = math.ceil(len(genes) / cols)
        w, h = 60 * cols, 30 * rows
        best = None
        along = normal(side)
        for sd in (side, (-side[0], -side[1])):
            for dist in (55, 80, 110, 150, 200):
                for shift in (0, 60, -60, 120, -120):
                    c = add(add(rpos, sd, dist + (h / 2 if abs(sd[1]) > 0.7 else w / 2)), along, shift)
                    box = (c[0] - w / 2, c[1] - h / 2, c[0] + w / 2, c[1] + h / 2)
                    if self.space.free(box, margin=6):
                        best = box
                        break
                if best:
                    break
            if best:
                break
        if best is None:
            c = add(rpos, side, 55 + (h / 2 if abs(side[1]) > 0.7 else w / 2))
            best = (c[0] - w / 2, c[1] - h / 2, c[0] + w / 2, c[1] + h / 2)
        for i, gene in enumerate(genes):
            x, y = best[0] + 60 * (i % cols), best[1] + 30 * (i // cols)
            append(self.mp.layer["nodes"], new_gene(gene, self.model.symbol.get(gene, gene), [rid], x, y))
        self.space.add(best)
        center = ((best[0] + best[2]) / 2, (best[1] + best[3]) / 2)
        start = boundary(rpos, center, "rea")
        end = (min(max(start[0], best[0]), best[2]), min(max(start[1], best[1]), best[3]))
        append(self.mp.layer["ees"], element("path", **{"class": "ee " + rid, "d": to_d([[start, end]])}))

    def side_slots(self, rpos, tangent, side, subs, prods):
        """Positions for cofactors: substrates before the node, products after it, on the free side."""
        out = []
        for group, sign in ((subs, -1), (prods, 1)):
            for i, m in enumerate(group):
                row, k = divmod(i, 2)
                p = add(add(rpos, side, 62 + 34 * row), tangent, sign * (40 + 46 * k))
                out.append((m, p, -1 if sign < 0 else 1))
        return out

    def choose_side(self, rpos, tangent):
        n = normal(tangent)
        scores = []
        for s in (n, (-n[0], -n[1])):
            probe = add(rpos, s, 80)
            box = (probe[0] - 70, probe[1] - 35, probe[0] + 70, probe[1] + 35)
            scores.append((0 if self.space.free(box, 0) else 1, s))
        scores.sort(key=lambda t: t[0])
        return scores[0][1]

    def add_side_mets(self, rid, rpos, tangent, side, exclude):
        st = self.model.rxns[rid]["stoich"]
        subs = [m for m in st if st[m] < 0 and m not in exclude]
        prods = [m for m in st if st[m] > 0 and m not in exclude]
        links = []
        for m, p, role in self.side_slots(rpos, tangent, side, subs, prods):
            left = (p[0] - rpos[0]) < 0
            el = new_side(m, self.side_name(m), [rid], p[0], p[1], left)
            append(self.mp.layer["nodes"], el)
            b = node_box(self.mp, el)
            if b:
                self.space.add(b)
            links.append((p, "side", role, None))
        return links


# ----------------------------------------------------------------------------- scope

ARTEFACT_SUBSYSTEMS = {"Artificial reactions", "Exchange/demand reactions"}
MAX_METABOLITES = 12  # a reaction with more metabolites is a lumped model reaction


def artefact(model, rid):
    """R11: why a reaction is a model artefact rather than a biochemical step (or None): a reaction of an
    artificial subsystem (biomass, exchange), a reaction with metabolites on one side only (an exchange, demand
    or sink, whatever its subsystem: Yeast-GEM files them under "Exchange reaction"), a pool reaction that forms
    a pool, or a lumped reaction with more than MAX_METABOLITES metabolites."""
    r = model.rxns.get(rid)
    if not r:
        return None
    subs = r.get("subsystem")
    subs = set(subs if isinstance(subs, list) else [subs])
    if subs & ARTEFACT_SUBSYSTEMS:
        return "artificial subsystem"
    coefs = r["stoich"].values()
    if coefs and (all(c < 0 for c in coefs) or all(c > 0 for c in coefs)):
        return "exchange (metabolites on one side only)"
    if "Pool reactions" in subs and "pool" in (r.get("name") or "").lower():
        return "pool reaction"
    if len(r["stoich"]) > MAX_METABOLITES:
        return f"{len(r['stoich'])} metabolites"
    return None


class Scope:
    """What belongs on a map: the reactions of a subsystem, or of a compartment (all metabolites in it).
    `only` limits the missing reactions to add, for compartments drawn on several maps (the cytosol)."""

    def __init__(self, model, subsystems=(), compartments=(), only=None):
        self.model, self.subsystems, self.compartments, self.only = model, set(subsystems), set(compartments), only

    def __bool__(self):
        return bool(self.subsystems or self.compartments)

    def contains(self, rid):
        r = self.model.rxns.get(rid)
        if not r:
            return False
        if self.subsystems:
            subs = r.get("subsystem")
            return bool(self.subsystems & set(subs if isinstance(subs, list) else [subs]))
        if self.compartments:
            return bool(r["stoich"]) and all(comp_of(m) in self.compartments for m in r["stoich"])
        return False

    def missing(self, drawn):
        out = [r for r in self.model.rxns if r not in drawn and self.contains(r) and not artefact(self.model, r)]
        if self.only is not None:
            out = [r for r in out if r in self.only]
        return sorted(out)


# ----------------------------------------------------------------------------- path search

def path_search(model, s, t, drawn, scope, max_steps=4, max_cost=4.0):
    """Cheapest reaction path s -> t over main metabolites in one compartment."""
    adj = path_search.cache.get(id(model))
    if adj is None:
        adj = collections.defaultdict(list)
        for rid, r in model.rxns.items():
            st = r["stoich"]
            mains = [m for m in st if model.name(m) not in COFACTORS]
            subs, prods = [m for m in mains if st[m] < 0], [m for m in mains if st[m] > 0]
            rev = (r.get("lower_bound") or 0) < 0
            for a in subs:
                for b in prods:
                    adj[a].append((rid, b))
                    if rev:
                        adj[b].append((rid, a))
        path_search.cache[id(model)] = adj

    def cost(rid):
        if rid in drawn:
            return 0.1
        if artefact(model, rid):
            return 1e9
        return 1.0 if scope and scope.contains(rid) else 3.0

    heap = [(0.0, 0, s, [])]
    best = {}
    while heap:
        c, n, m, path = heapq.heappop(heap)
        if m == t and path:
            return c, path
        if n >= max_steps or best.get(m, 1e9) <= c:
            continue
        best[m] = c
        for rid, nxt in adj[m]:
            if comp_of(nxt) != comp_of(s) or any(rid == r for r, _ in path):
                continue
            nc = c + cost(rid)
            if nc <= max_cost:
                heapq.heappush(heap, (nc, n + 1, nxt, path + [(rid, nxt)]))
    return None


path_search.cache = {}


# ----------------------------------------------------------------------------- D1 / D2

def removed_line(rec, s_pos, p_pos):
    """Polyline from the substrate node to the product node along the removed reaction's edges."""
    rpos = rec["pos"]
    def piece(npos):
        for sp in rec["fe"]:
            if math.hypot(sp[0][0] - npos[0], sp[0][1] - npos[1]) < 50:
                return sp
            if math.hypot(sp[-1][0] - npos[0], sp[-1][1] - npos[1]) < 50:
                return list(reversed(sp))
        return None
    a, b = piece(s_pos), piece(p_pos)
    pts = [s_pos] + (a[1:-1] if a else []) + [rpos] + (list(reversed(b))[1:-1] if b else []) + [p_pos]
    return [p for i, p in enumerate(pts) if i == 0 or p != pts[i - 1]]


def restore_paths(dr, removed, scope):
    mp, model, log = dr.mp, dr.model, dr.log
    drawn0 = set(mp.reactions())

    def plen(rec):
        old = model.old_rxns.get(rec["id"])
        if not old:
            return 0
        st = {model.fix_met_id(m, []) or m: c for m, c in old["stoich"].items()}
        ends = [e for e in rec["mains"] if e[2] in model.mets and model.name(e[2]) not in COFACTORS]
        best = 0
        for s in [e for e in ends if st.get(e[2], 0) < 0]:
            for p in [e for e in ends if st.get(e[2], 0) > 0]:
                f = path_search(model, s[2], p[2], drawn0, scope)
                if f:
                    best = max(best, len(f[1]))
        return best

    for rec in sorted(removed, key=plen, reverse=True):
        old = model.old_rxns.get(rec["id"])
        if not old:
            continue
        st = {model.fix_met_id(m, []) or m: c for m, c in old["stoich"].items()}
        ends = [(n, pos, mid) for n, pos, mid in rec["mains"] if n.getparent() is not None and mid in model.mets]
        subs = [e for e in ends if st.get(e[2], 0) < 0 and model.name(e[2]) not in COFACTORS]
        prods = [e for e in ends if st.get(e[2], 0) > 0 and model.name(e[2]) not in COFACTORS]
        drawn = set(mp.reactions())
        best = None
        for s in subs:
            for p in prods:
                found = path_search(model, s[2], p[2], drawn, scope)
                if found and (best is None or found[0] < best[0]):
                    best = (found[0], found[1], s, p)
        if not best:
            ends_txt = " / ".join(model.name(e[2]) for e in subs) + " -> " + " / ".join(model.name(e[2]) for e in prods)
            log.append(("review", "connection not restored (no path in the model)", rec["id"], ends_txt, rec["name"]))
            continue
        cost, path, s, p = best
        rids = [r for r, _ in path]
        if any(r in dr.placed for r in rids):
            continue
        all_drawn = all(r in drawn for r in rids)
        if all_drawn and len(rids) == 1 and connects(mp, rids[0], s[0], p[0]):
            # an existing reaction already joins the same two nodes: move it onto the line only when
            # the line is part of a circle (a cycle drawn as a ring) and the reaction is drawn away from it
            line = removed_line(rec, s[1], p[1])
            g = mp.reactions()[rids[0]]
            if not (on_circle(line) and min(length(sub(translate(g), q)) for q in densify(line)) > 100):
                continue
        inter = [m for _, m in path[:-1]]
        if all_drawn and not (rids[0] in mp.classes(s[0]) and rids[-1] in mp.classes(p[0])):
            log.append(("review", "path drawn elsewhere between other copies, not moved", rec["id"], " ".join(rids),
                        rec["name"]))
            continue
        if all_drawn and not movable(mp, rids, inter):
            log.append(("review", "path drawn elsewhere, not moved (shared nodes)", rec["id"], " ".join(rids), rec["name"]))
            continue
        if any(r in drawn for r in rids) and not all_drawn:
            log.append(("review", "path partly drawn, not changed", rec["id"], " ".join(rids), rec["name"]))
            continue
        line = removed_line(rec, s[1], p[1])
        lay_path(dr, rids, [s[2]] + inter + [p[2]], s, p, line)
        dr.lines.append(line)
        rule = "D1" if all_drawn else "D2"
        what = "path moved onto the removed reaction's line" if all_drawn else "missing path drawn on the removed reaction's line"
        log.append((rule, what, rec["id"], " ".join(rids),
                    " -> ".join(model.name(m) for m in [s[2]] + inter + [p[2]])))


def densify(line, step=20):
    out = []
    for a, b in zip(line, line[1:]):
        k = max(1, int(length(sub(b, a)) / step))
        out += [add(a, sub(b, a), i / k) for i in range(k)]
    return out + [line[-1]]


def on_circle(line):
    fit = fit_circle(densify(line))
    return bool(fit) and fit[1] >= 200 and fit[2] <= 12


def connects(mp, rid, a, b):
    return rid in mp.classes(a) and rid in mp.classes(b)


def movable(mp, rids, inter):
    for m in inter:
        for n in [n for n in mp.nodes("met") if mp.classes(n)[1] == m]:
            others = [c for c in mp.classes(n)[2:] if is_rxn(c) and c not in rids]
            if others:
                return False
    return True


def lay_path(dr, rids, mets, s, p, line):
    """Place reactions rids (joining mets[0] -> mets[-1]) evenly along line, moving or creating nodes."""
    mp, model = dr.mp, dr.model
    total = poly_length(line)
    k = len(rids)
    slots = [total * i / (2 * k) for i in range(2 * k + 1)]
    node_pos = {0: s[1], 2 * k: p[1]}
    node_el = {0: s[0], 2 * k: p[0]}
    # intermediate metabolite nodes
    for i, m in enumerate(mets[1:-1], start=1):
        pos, _ = point_at(line, slots[2 * i])
        existing = [n for n in dr.met_nodes(m) if any(r in mp.classes(n) for r in rids)]
        if existing:
            el = existing[0]
            set_position(el, *pos)
            for extra in existing[1:]:
                for r in rids:
                    if r in mp.classes(extra):
                        dr.unlink(extra, r)
        else:
            el = new_main(m, model.name(m), [], pos[0], pos[1])
            append(mp.layer["nodes"], el)
        node_pos[2 * i], node_el[2 * i] = pos, el
        dr.space.add(node_box(mp, el))
    for i, rid in enumerate(rids):
        rpos, tangent = point_at(line, slots[2 * i + 1])
        a_el, b_el = node_el[2 * i], node_el[2 * i + 2]
        a_id, b_id = mets[i], mets[i + 1]
        # detach this reaction from every other drawn copy of anything, then rebuild it
        for n in mp.nodes("met"):
            if rid in mp.classes(n) and n is not a_el and n is not b_el:
                dr.unlink(n, rid)
        dr.link(a_el, rid)
        dr.link(b_el, rid)
        g = mp.reactions().get(rid)
        if g is not None:
            dr.remove_edges(rid)
            set_position(g, *rpos)
        else:
            append(mp.layer["nodes"], new_reaction(rid, *rpos))
        side = dr.choose_side(rpos, tangent)
        st = model.rxns[rid]["stoich"]
        links = [(node_pos[2 * i], "main", 1 if st.get(a_id, 0) > 0 else -1, cut(line, slots[2 * i], slots[2 * i + 1])),
                 (node_pos[2 * i + 2], "main", 1 if st.get(b_id, 0) > 0 else -1,
                  list(reversed(cut(line, slots[2 * i + 1], slots[2 * i + 2]))))]
        links += dr.add_side_mets(rid, rpos, tangent, side, {a_id, b_id})
        dr.draw(rid, rpos, links, (-side[0], -side[1]))
        dr.placed.add(rid)


# ----------------------------------------------------------------------------- D3

def add_missing(dr, scope, comps_on_map):
    """D3: every reaction in scope on the map: attached next to a drawn metabolite where that is clean,
    otherwise collected for placement (D9, placement.py)."""
    mp = dr.mp
    dr.comps_on_map = set(comps_on_map)
    drawn = set(mp.reactions())
    for rid in scope.missing(drawn):
        draw_missing(dr, rid, comps_on_map)
    for rid in sorted(r for r in dr.model.rxns if r not in drawn and scope.contains(r)):
        why = artefact(dr.model, rid)
        if why:
            dr.log.append(("R11", "model artefact, not drawn", rid, why, dr.model.rxns[rid].get("name", "") or ""))


def draw_missing(dr, rid, comps_on_map):
    model, log = dr.model, dr.log
    st = model.rxns[rid]["stoich"]
    mains = [m for m in st if model.name(m) not in COFACTORS] or list(st)
    comp = COMPARTMENT.get(comp_of(mains[0]), "")
    if dr.one_area:  # a map without compartment areas: every reaction goes in the one area
        comp = dr.one_area
    elif len({comp_of(m) for m in mains}) > 1:
        log.append(("review", "missing reaction spans compartments, not drawn", rid, "", model.rxns[rid].get("name", "")))
        return
    subs = [m for m in mains if st[m] < 0]
    prods = [m for m in mains if st[m] > 0]
    anchors = [(m, ns[0]) for m in mains for ns in [[n for n in dr.met_nodes(m) if dr.is_main(n)]] if ns]
    # a reaction joining two drawn metabolites goes between them; everything else is placed by D9
    if comp in comps_on_map and any(m in subs for m, _ in anchors) and any(m in prods for m, _ in anchors) \
            and attach_local(dr, rid, subs, prods, anchors):
        log.append(("D14" if rid in dr.context else "D3", "context reaction drawn" if rid in dr.context else
                    "missing reaction added", rid, "between " + ", ".join(model.name(m) for m, _ in anchors),
                    model.rxns[rid].get("name", "") or ""))
        return
    dr.pending.append((rid, comp))


def module_box(model, center_x, y, subs, prods, rid):
    """Approximate box of a reaction drawn in a row with its main metabolites and cofactors."""
    labels = [label_width(line, 18) for m in subs + prods for line in wrap(model.name(m))] or [60]
    half = STEP + max(labels) / 2 + 20
    n_side = len(model.rxns[rid]["stoich"]) - len(subs) - len(prods)
    up = 70 + 34 * max(0, (n_side - 1) // 4)
    genes = math.ceil(len(model.rxns[rid]["genes"]) / 3) * 30
    down = 50 + 55 * (max(len(subs), len(prods)) - 1) + (genes + 40 if genes else 0)
    return (center_x - half, y - up, center_x + half, y + down)


def attach_local(dr, rid, subs, prods, anchors, near=450):
    """Draw rid next to its drawn main metabolite(s) if there is free room close by; False otherwise."""
    model = dr.model
    ids = [m for m, _ in anchors]
    a_sub = [(m, n) for m, n in anchors if m in subs]
    a_prod = [(m, n) for m, n in anchors if m in prods]
    if a_sub and a_prod:
        (ms, ns), (mp_, np_) = a_sub[0], a_prod[0]
        pa, pb = translate(ns), translate(np_)
        if length(sub(pb, pa)) > near:
            return False
        tangent = unit(sub(pb, pa))
        mid = add(pa, sub(pb, pa), 0.5)
        for off in (0, 120, -120, 240, -240):
            c = add(mid, normal(tangent), off)
            if dr.space.free((c[0] - 45, c[1] - 45, c[0] + 45, c[1] + 45), 4):
                pos = {ms: pa, mp_: pb}
                extra = [m for m in subs + prods if m not in pos]
                if extra:
                    return False
                elbow = {ms: add(pa, normal(tangent), off), mp_: add(pb, normal(tangent), off)} if off else {}
                build(dr, rid, subs, prods, pos, c, tangent, dict(anchors), elbow)
                return True
        return False
    if len(set(ids)) != 1:
        return False
    m0, n0 = anchors[0]
    p0 = translate(n0)
    others = [m for m in (prods if m0 in subs else subs) if m != m0]
    same = [m for m in (subs if m0 in subs else prods) if m != m0]
    for ang in range(0, 360, 45):
        u = (math.cos(math.radians(ang)), math.sin(math.radians(ang)))
        rpos = add(p0, u, STEP)
        far = add(p0, u, 2 * STEP)
        pos = {m0: p0}
        for k, m in enumerate(others):
            pos[m] = add(far, normal(u), 55 * k)
        for k, m in enumerate(same):
            pos[m] = add(p0, normal(u), 55 * (k + 1))
        xs = [q[0] for q in pos.values()] + [rpos[0]]
        ys = [q[1] for q in pos.values()] + [rpos[1]]
        lw = max([label_width(line, 18) for m in others + same for line in wrap(model.name(m))] or [60]) / 2
        box = (min(xs) - lw, min(ys) - 80, max(xs) + lw, max(ys) + 80)
        if dr.space.free(box, 6):
            tangent = u if m0 in subs else (-u[0], -u[1])
            build(dr, rid, subs, prods, pos, rpos, tangent, {m0: n0}, {})
            return True
    return False


def build(dr, rid, subs, prods, pos, rpos, tangent, existing, elbow, cof_side=None, gene_side=None):
    """Draw reaction rid at rpos with its main metabolites at pos (reusing the nodes in existing)."""
    mp, model = dr.mp, dr.model
    links = []
    for m in subs + prods:
        n = existing.get(m)
        if n is None:
            n = new_main(m, model.name(m), [rid], *pos[m])
            append(mp.layer["nodes"], n)
            dr.space.add(node_box(mp, n))
        else:
            dr.link(n, rid)
        bend = elbow.get(m)
        links.append((pos[m], "main", -1 if m in subs else 1, [pos[m], bend, rpos] if bend else None))
    append(mp.layer["nodes"], new_reaction(rid, *rpos))
    side = cof_side or dr.choose_side(rpos, tangent)
    links += dr.add_side_mets(rid, rpos, tangent, side, set(subs + prods))
    dr.draw(rid, rpos, links, gene_side or (-side[0], -side[1]))
    dr.placed.add(rid)


def chains(model, rids):
    """Order reactions into chains where a main product of one is a main substrate of the next."""
    def mains(r, sign):
        st = model.rxns[r]["stoich"]
        return {m for m in st if st[m] * sign > 0 and model.name(m) not in COFACTORS}
    succ = {r: [q for q in rids if q != r and mains(r, 1) & mains(q, -1)] for r in rids}
    pred = collections.Counter(q for r in rids for q in succ[r])
    left, out = list(rids), []
    while left:
        start = next((r for r in left if pred[r] == 0), left[0])
        chain = [start]
        left.remove(start)
        while True:
            nxt = next((q for q in succ[chain[-1]] if q in left), None)
            if nxt is None:
                break
            chain.append(nxt)
            left.remove(nxt)
        out.append(chain)
    return out


def new_compartment_box(dr, comp, height, width=None, at=None):
    """A labelled box for a compartment below the existing content, for its added reactions."""
    mp = dr.mp
    template = next((g for g in mp.root.iter(SVG + "g") if g.get("class") == "compartment"
                     and g.find(SVG + "rect") is not None), None)
    w = float(mp.root.get("width"))
    bw = min(w - 500, width) if width else w - 500
    # next to the boxes added before, in the same row, when it fits; otherwise in a new row below
    added = [g.find(SVG + "rect") for g in mp.root.iter(SVG + "g")
             if g.get("class") == "compartment" and (g.get("id") or "").endswith(" (added)")]
    added = [(float(r.get("x")), float(r.get("y")), float(r.get("x")) + float(r.get("width")),
              float(r.get("y")) + float(r.get("height"))) for r in added if r is not None]
    row = [b for b in added if added and b[1] == max(a[1] for a in added)]
    if at is not None:
        x0, top = at
        need = top + height + 150 + 300 - float(mp.root.get("height"))
        if need > 0:
            grow_canvas(dr, need)
    elif row and max(b[2] for b in row) + 150 + bw <= w - 300:
        x0, top = max(b[2] for b in row) + 150, row[0][1]
        need = top + height + 150 + 300 - float(mp.root.get("height"))
        if need > 0:
            grow_canvas(dr, need)
    else:
        x0, top = 200.0, float(mp.root.get("height")) - 100
        grow_canvas(dr, height + 300)
    box = (x0, top, x0 + bw, top + height + 150)
    g = element("g", id=comp + " (added)", **{"class": "compartment", "data-moved": "1"})
    rect = element("rect", **{"class": "shape", "x": fmt(box[0]), "y": fmt(box[1]), "width": fmt(box[2] - box[0]),
                              "height": fmt(box[3] - box[1]), "rx": "150", "ry": "150", "fill": "none", "stroke": "#000",
                              "stroke-width": "22"})
    g.append(rect)
    size = "100"
    if template is not None:
        t = template.find(SVG + "text")
        if t is not None and t.get("font-size"):
            size = t.get("font-size")
    g.append(text_el(comp, **{"class": "lbl", "font-family": "sans-serif", "font-size": size, "font-weight": "700",
                              "x": fmt(box[0] + 60), "y": fmt(box[1] + float(size) + 30)}))
    mp.layer["fes"].addprevious(g)
    g.tail = "\n  "
    dr.space.add(text_box(g.find(SVG + "text")))
    return (box[0], box[1] + float(size) + 60, box[2], box[3])


def grow_canvas(dr, dh):
    """Make the map taller by dh; the licence badge moves down with the bottom edge."""
    mp = dr.mp
    h = float(mp.root.get("height"))
    mp.root.set("height", fmt(h + dh))
    for g in mp.root.iter(SVG + "g"):
        if g.get("id") == "license":
            m = re.match(r"translate\(([-0-9.e]+),([-0-9.e]+)\)", g.get("transform") or "")
            dx, dy = (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0)
            g.set("transform", "translate(%s,%s)" % (fmt(dx), fmt(dy + dh)))
    dr.log.append(("D3", "map made taller for the bottom row", "", f"{dh:.0f} px", ""))
    return canvas(mp)


def move_box(mp, g, dx, dy):
    """Move a compartment box with everything inside it (nodes, edges, arrowheads, band, heading)."""
    rect = g.find(SVG + "rect")
    x0, y0 = float(rect.get("x")), float(rect.get("y"))
    x1, y1 = x0 + float(rect.get("width")), y0 + float(rect.get("height"))

    def f(p):
        return (p[0] + dx, p[1] + dy) if x0 <= p[0] <= x1 and y0 <= p[1] <= y1 else p
    for n in list(mp.layer["nodes"]):
        p = translate(n)
        cls = mp.classes(n)
        if not p or not cls:
            continue
        ref = (p[0] + 30, p[1] + 15) if cls[0] == "enz" else p
        if f(ref) != ref:
            set_position(n, p[0] + dx, p[1] + dy)
    paths = [path for layer in ("fes", "fehs", "ees") for path in mp.layer[layer]] + [b for _, b in band_groups(mp)]
    for path in paths:
        sps = subpaths(path.get("d") or "")
        if sps:
            path.set("d", to_d([[f(q) for q in sp] for sp in sps]))
    rect.set("x", fmt(x0 + dx))
    rect.set("y", fmt(y0 + dy))
    for t in g.iter(SVG + "text"):
        m = re.match(r"matrix\(1,0,0,1,([-0-9.e]+),([-0-9.e]+)\)", t.get("transform") or "")
        if m:
            t.set("transform", "matrix(1,0,0,1,%s,%s)" % (fmt(float(m.group(1)) + dx), fmt(float(m.group(2)) + dy)))
        elif t.get("x"):
            t.set("x", fmt(float(t.get("x")) + dx))
            t.set("y", fmt(float(t.get("y") or 0) + dy))


def widen_canvas(dr, dw):
    """Make the map wider by dw; the licence badge moves right with the edge."""
    mp = dr.mp
    mp.root.set("width", fmt(float(mp.root.get("width")) + dw))
    for g in mp.root.iter(SVG + "g"):
        if g.get("id") == "license":
            m = re.match(r"translate\(([-0-9.e]+),([-0-9.e]+)\)", g.get("transform") or "")
            dx, dy = (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0)
            g.set("transform", "translate(%s,%s)" % (fmt(dx + dw), fmt(dy)))
    dr.log.append(("D3", "map made wider", "", f"{dw:.0f} px", ""))
    return canvas(mp)


# ----------------------------------------------------------------------------- band

def band_groups(mp):
    """(subsystem group, band path) for every subsystem highlight on the map."""
    out = []
    for g in mp.root.iter(SVG + "g"):
        if g.get("class") == "subsystem":
            path = g.find(SVG + "path")
            if path is not None:
                out.append((g, path))
    return out


def extend_band(mp, model, band_add, skip=frozenset()):
    """Band along new main edges, in the group of the reaction's subsystem (a new grey group if the map has none).
    Reactions in skip (context reactions of other subsystems on a subsystem map) get no band."""
    groups = band_groups(mp)
    if not band_add or not groups:
        return
    by_id = {g.get("id"): (g, p) for g, p in groups}
    for rid, line in band_add:
        if rid in skip:
            continue
        subs = model.rxns.get(rid, {}).get("subsystem")
        subs = subs if isinstance(subs, list) else [subs]
        target = next((by_id[x] for x in subs if x in by_id), None)
        if target is None and len(groups) == 1:
            target = groups[0]
        if target is None and subs and subs[0]:
            g = element("g", id=subs[0], **{"class": "subsystem"})
            path = element("path", **{"class": "shape", "d": "", "fill": "none", "stroke": "#bdbdbd",
                                      "stroke-width": "80", "opacity": ".5", "stroke-linejoin": "round",
                                      "stroke-linecap": "round"})
            g.append(path)
            mp.layer["fes"].addprevious(g)
            g.tail = "\n  "
            target = (g, path)
            by_id[subs[0]] = target
            groups.append(target)
        if target is None:
            continue
        path = target[1]
        path.set("d", (path.get("d") or "") + to_d([line]))
        path.set("stroke-linecap", "round")


def map_scope(mp, subsystem_of_map, name):
    s = subsystem_of_map.get(name)
    return {s} if s else set()


def compartments_on_map(mp):
    names = {g.get("id") for g in mp.root.iter(SVG + "g") if g.get("class") == "compartment"}
    return names | {k for k, v in ALIASES.items() if v in names}



# ----------------------------------------------------------------------------- D4 rings

def fit_circle(pts):
    """Least-squares circle through points: (center, radius, rms residual)."""
    n = len(pts)
    if n < 3:
        return None
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    u = [p[0] - mx for p in pts]
    v = [p[1] - my for p in pts]
    suu = sum(a * a for a in u)
    svv = sum(b * b for b in v)
    suv = sum(a * b for a, b in zip(u, v))
    suuu = sum(a ** 3 for a in u)
    svvv = sum(b ** 3 for b in v)
    suvv = sum(a * b * b for a, b in zip(u, v))
    svuu = sum(b * a * a for a, b in zip(u, v))
    det = suu * svv - suv * suv
    if abs(det) < 1e-9:
        return None
    uc = (svv * (suuu + suvv) - suv * (svvv + svuu)) / (2 * det)
    vc = (suu * (svvv + svuu) - suv * (suuu + suvv)) / (2 * det)
    c = (uc + mx, vc + my)
    r = math.sqrt(uc * uc + vc * vc + (suu + svv) / n)
    rms = math.sqrt(sum((length(sub(p, c)) - r) ** 2 for p in pts) / n)
    return c, r, rms


def rotate(p, c, a):
    dx, dy = p[0] - c[0], p[1] - c[1]
    return (c[0] + dx * math.cos(a) - dy * math.sin(a), c[1] + dx * math.sin(a) + dy * math.cos(a))


def arc(c, r, a0, a1, step=math.radians(6)):
    d = (a1 - a0 + math.pi) % (2 * math.pi) - math.pi
    n = max(1, int(abs(d) / step))
    return [(c[0] + r * math.cos(a0 + d * i / n), c[1] + r * math.sin(a0 + d * i / n)) for i in range(n + 1)]


def respace_rings(dr):
    """D4: when a path was laid on a circular backbone, spread the ring evenly.

    Ring nodes that connect to structures off the ring move together with those structures as rigid
    blocks (rotated about the circle centre); the free nodes between blocks are spaced evenly."""
    mp = dr.mp
    done = []
    for line in dr.lines:
        fit = fit_circle(densify(line))
        if not fit:
            continue
        c, r, rms = fit
        if r < 200 or rms > 12 or any(length(sub(c, d[0])) < 40 and abs(r - d[1]) < 40 for d in done):
            continue
        ring = ring_members(mp, c, r)
        if ring is None:
            continue
        done.append((c, r))
        respace(dr, c, r, *ring)
    if done and len(band_groups(mp)) == 1:
        regenerate_band(mp)


def ring_members(mp, c, r, tol=30):
    """Reaction and main metabolite nodes forming a closed ring on the circle, in angle order."""
    cand_r, cand_m = {}, {}
    for n in mp.layer["nodes"]:
        cls = mp.classes(n)
        p = translate(n) if cls else None
        if not p or abs(length(sub(p, c)) - r) > tol:
            continue
        if cls[0] == "rea":
            cand_r[n.get("id")] = n
        elif cls[0] == "met" and n.find(SVG + "ellipse") is not None:
            cand_m[id(n)] = n
    changed = True
    while changed:
        changed = False
        for k, n in list(cand_m.items()):
            if len([x for x in mp.classes(n)[2:] if x in cand_r]) < 2:
                del cand_m[k]
                changed = True
        for rid in list(cand_r):
            if not any(rid in mp.classes(n) for n in cand_m.values()):
                del cand_r[rid]
                changed = True
    items = [(math.atan2(translate(n)[1] - c[1], translate(n)[0] - c[0]), n) for n in
             list(cand_r.values()) + list(cand_m.values())]
    if len(items) < 6:
        return None
    items.sort(key=lambda t: t[0])
    gaps = [(items[(i + 1) % len(items)][0] - items[i][0]) % (2 * math.pi) for i in range(len(items))]
    if max(gaps) > 3 * 2 * math.pi / len(items):
        return None
    return items, set(cand_r), cand_m


def respace(dr, c, r, items, ring_rids, ring_mets):
    mp, model, log = dr.mp, dr.model, dr.log
    ring_ids = {id(n) for _, n in items}
    # off-ring components hanging from ring metabolites
    comp_of = {}
    for _, n in items:
        if id(n) not in ring_mets:
            continue
        comp = {"reas": set(), "mets": set(), "anchors": {id(n)}, "keep": []}
        stack = [x for x in mp.classes(n)[2:] if is_rxn(x) and x not in ring_rids]
        while stack:
            rid = stack.pop()
            if rid in comp["reas"]:
                continue
            comp["reas"].add(rid)
            for m in mp.nodes("met"):
                if rid in mp.classes(m):
                    if id(m) in ring_mets:
                        comp["anchors"].add(id(m))
                    elif id(m) not in comp["mets"]:
                        comp["keep"].append(m)  # keeps the proxy alive, so id(m) stays this element's
                        comp["mets"].add(id(m))
                        stack += [x for x in mp.classes(m)[2:] if is_rxn(x) and x not in ring_rids]
        comp_of[id(n)] = comp
    # group anchors joined by a component (union-find)
    parent = {id(n): id(n) for _, n in items}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for comp in comp_of.values():
        a = list(comp["anchors"])
        for b in a[1:]:
            parent[find(b)] = find(a[0])
    # blocks: runs of items between the first and last member of a group (in angle order)
    N = len(items)
    order = {id(n): i for i, (_, n) in enumerate(items)}
    groups = collections.defaultdict(list)
    for _, n in items:
        if comp_of.get(id(n), {}).get("reas"):
            groups[find(id(n))].append(order[id(n)])
    block_of = [None] * N
    for gid, idxs in groups.items():
        if len(idxs) < 2:
            block_of[idxs[0]] = gid
            continue
        idxs.sort()
        # the shorter way round between the outermost members
        span_fwd = idxs[-1] - idxs[0]
        best = (span_fwd, idxs[0], idxs[-1])
        for k in range(len(idxs)):
            a, b = idxs[k], idxs[k - 1]
            span = (b - a) % N
            if span < best[0]:
                best = (span, a, b)
        _, a, b = best
        k = a
        while True:
            block_of[k] = gid
            if k == b:
                break
            k = (k + 1) % N
    # units: consecutive items with the same block form one rigid unit
    ang = [t[0] for t in items]
    start = next((i for i in range(N) if block_of[i] is None or block_of[i] != block_of[i - 1]), 0)
    units, i = [], start
    seen = 0
    while seen < N:
        unit = [i]
        seen += 1
        while seen < N and block_of[(i + 1) % N] is not None and block_of[(i + 1) % N] == block_of[i]:
            i = (i + 1) % N
            unit.append(i)
            seen += 1
        units.append(unit)
        i = (i + 1) % N
    spans = [((ang[u[-1]] - ang[u[0]]) % (2 * math.pi)) if len(u) > 1 else 0.0 for u in units]
    gap = (2 * math.pi - sum(spans)) / len(units)
    if gap <= 0:
        return
    new = [0.0] * N
    cur = 0.0
    for u, sp in zip(units, spans):
        for k in u:
            new[k] = cur + ((ang[k] - ang[u[0]]) % (2 * math.pi))
        cur += sp + gap
    diffs = [new[k] - ang[k] for k in range(N)]
    shift = math.atan2(sum(math.sin(-d) for d in diffs), sum(math.cos(-d) for d in diffs))
    delta = {id(items[k][1]): ((new[k] + shift - ang[k] + math.pi) % (2 * math.pi)) - math.pi for k in range(N)}
    new_angle = {id(items[k][1]): ang[k] + delta[id(items[k][1])] for k in range(N)}
    old_pos = {id(n): translate(n) for _, n in items}
    new_pos = {k: (c[0] + r * math.cos(a), c[1] + r * math.sin(a)) for k, a in new_angle.items()}

    def near_node(pt, npos, rad=55):
        return length(sub(pt, npos)) < rad

    rotated = []  # nodes already rotated (kept alive, compared by identity)

    def done_before(n):
        return any(n is x for x in rotated)

    def rotate_reaction_attachments(rid, a, keep_near, fn=None):
        """Move a reaction's own side nodes, genes, gene edges and the edges not ending at keep_near:
        rotated about the circle centre by a, or by the point function fn."""
        if fn is None:
            def fn(q):
                return rotate(q, c, a)
        for n in mp.nodes("met") + mp.nodes("enz"):
            cls = mp.classes(n)
            if rid in cls and id(n) not in ring_ids and [x for x in cls[2:] if is_rxn(x)] == [rid] \
                    and not done_before(n):
                rotated.append(n)
                p = translate(n)
                if p is None:
                    continue
                if cls[0] == "enz":
                    q = fn((p[0] + 30, p[1] + 15))
                    set_position(n, q[0] - 30, q[1] - 15)
                else:
                    set_position(n, *fn(p))
        for p in mp.edges("ees", rid):
            p.set("d", to_d([[fn(q) for q in sp] for sp in (subpaths(p.get("d")) or [])]))
        fe, feh = [], []
        for p in mp.edges("fes", rid):
            fe += [[fn(q) for q in sp] for sp in (subpaths(p.get("d")) or [])
                   if not any(near_node(sp[0], o) or near_node(sp[-1], o) for o in keep_near)]
        for p in mp.edges("fehs", rid):
            feh += [[fn(q) for q in sp] for sp in (subpaths(p.get("d")) or [])
                    if not any(near_node(sp[0], o, 60) for o in keep_near)]
        return fe, feh

    # ring reactions: rotate attachments, redraw ring edges as arcs
    for _, g in items:
        if id(g) in ring_mets:
            continue
        rid = g.get("id")
        a = delta[id(g)]
        nbrs = [n for n in ring_mets.values() if rid in mp.classes(n)]
        fe, feh = rotate_reaction_attachments(rid, a, [old_pos[id(n)] for n in nbrs])
        rev = (model.rxns.get(rid, {}).get("lower_bound") or 0) < 0
        st = model.rxns.get(rid, {}).get("stoich", {})
        rpos = new_pos[id(g)]
        for n in nbrs:
            pts = arc(c, r, new_angle[id(n)], new_angle[id(g)])
            pts[0] = boundary(new_pos[id(n)], pts[1], "main")
            pts[-1] = boundary(rpos, pts[-2], "rea")
            fe.append(pts)
            if st.get(mp.classes(n)[1], 0) > 0 or rev:
                feh.append(arrow(pts[0], pts[1]))
        for k, p in enumerate(mp.edges("fes", rid)):
            p.set("d", to_d(fe) if k == 0 else "")
        for k, p in enumerate(mp.edges("fehs", rid)):
            p.set("d", to_d(feh) if k == 0 else "")
        set_position(g, *rpos)
    # ring metabolites move along the circle; the structures hanging from them are shifted, not turned,
    # by the mean movement of the ring nodes they hang from, and their edges into those nodes follow
    group_shift = {}
    for k in ring_mets:
        root = find(k)
        if root not in group_shift:
            members = [m for m in ring_mets if find(m) == root]
            olds = [old_pos[m] for m in members]
            news = [new_pos[m] for m in members]
            group_shift[root] = ((sum(p[0] for p in news) - sum(p[0] for p in olds)) / len(members),
                                 (sum(p[1] for p in news) - sum(p[1] for p in olds)) / len(members))
    moved_reas = set()
    for k, n in ring_mets.items():
        set_position(n, *new_pos[k])
        comp = comp_of.get(k)
        if not comp:
            continue
        t = group_shift[find(k)]
        anchors = [(add(old_pos[m], t), new_pos[m]) for m in comp["anchors"]]

        def shift_pt(q, t=t):
            return add(q, t)
        for rid in comp["reas"] - moved_reas:
            moved_reas.add(rid)
            g = mp.reactions().get(rid)
            if g is not None and translate(g):
                set_position(g, *shift_pt(translate(g)))
            fe, feh = rotate_reaction_attachments(rid, 0, [], shift_pt)
            for sp in fe:
                for idx in (0, -1):
                    for shifted_old, new in anchors:
                        if near_node(sp[idx], shifted_old, 60):
                            sp[idx] = boundary(new, sp[1] if idx == 0 else sp[-2], "main")
            moved_feh = []
            for tri in feh:
                hit = next(((so, nw) for so, nw in anchors if near_node(tri[0], so, 60)), None)
                moved_feh.append([add(q, sub(hit[1], hit[0])) for q in tri] if hit else tri)
            feh = moved_feh
            for i2, p in enumerate(mp.edges("fes", rid)):
                p.set("d", to_d(fe) if i2 == 0 else "")
            for i2, p in enumerate(mp.edges("fehs", rid)):
                p.set("d", to_d(feh) if i2 == 0 else "")
        for m in comp["keep"]:
            if not done_before(m) and translate(m):
                rotated.append(m)
                set_position(m, *shift_pt(translate(m)))
    log.append(("D4", "ring re-spaced", "", f"{N} nodes, {len(units)} units", f"circle r={r:.0f}"))


def regenerate_band(mp):
    """Background band along every edge between a main metabolite node and a reaction."""
    band = next((g.find(SVG + "path") for g in mp.root.iter(SVG + "g") if g.get("class") == "subsystem"), None)
    if band is None:
        return
    mains = [translate(n) for n in mp.nodes("met") if n.find(SVG + "ellipse") is not None]
    reas = [translate(g) for g in mp.reactions().values()]
    lines = []
    for p in mp.layer["fes"]:
        for sp in subpaths(p.get("d")) or []:
            ends = [sp[0], sp[-1]]
            m = next((q for q in mains if q and any(length(sub(e, q)) < 45 for e in ends)), None)
            rr = next((q for q in reas if q and any(length(sub(e, q)) < 45 for e in ends)), None)
            if m and rr:
                pts = list(sp)
                if length(sub(pts[0], m)) > length(sub(pts[-1], m)):
                    pts.reverse()
                lines.append([m] + pts + [rr])
    band.set("d", to_d(lines))
    band.set("stroke-linecap", "round")


# ----------------------------------------------------------------------------- D5

def erase_reaction(dr, rid):
    mp = dr.mp
    g = mp.reactions().get(rid)
    if g is not None:
        mp.remove(g)
    dr.remove_edges(rid)
    for n in mp.nodes("met") + mp.nodes("enz"):
        if rid in mp.classes(n):
            dr.unlink(n, rid)


def add_cofactor_dots(dr, rid, g, mets):
    """Missing cofactors as dots around the reaction, near the dots of the same role where there are any."""
    if not mets:
        return
    mp, model, log = dr.mp, dr.model, dr.log
    st = model.rxns[rid]["stoich"]
    rpos = translate(g)
    dots = [n for n in mp.nodes("met") if rid in mp.classes(n) and n.find(SVG + "circle") is not None]
    fe, feh = [], []
    rev = (model.rxns[rid].get("lower_bound") or 0) < 0
    for m in mets:
        same = [translate(n) for n in dots if st.get(mp.classes(n)[1], 0) * st[m] > 0 and translate(n)]
        target = (sum(p[0] for p in same) / len(same), sum(p[1] for p in same) / len(same)) if same else None
        best = None
        for rad in (75, 105, 140):
            for k in range(18):
                a = 2 * math.pi * k / 18
                c = add(rpos, (math.cos(a), math.sin(a)), rad)
                left = c[0] < rpos[0]
                w = label_width(model.name(m), 18) + 20
                box = (c[0] - (w if left else 10), c[1] - 12, c[0] + (10 if left else w), c[1] + 12)
                if dr.space.free(box, 4):
                    score = length(sub(c, target)) if target else rad
                    if best is None or score < best[0]:
                        best = (score, c, left)
            if best:
                break
        if best is None:
            log.append(("review", "cofactor could not be placed", rid, m, model.name(m)))
            continue
        _, c, left = best
        el = new_side(m, model.name(m), [rid], c[0], c[1], left)
        append(mp.layer["nodes"], el)
        dr.space.add(node_box(mp, el))
        pts = [boundary(c, rpos, "side"), boundary(rpos, c, "rea")]
        fe.append(pts)
        if st[m] > 0 or rev:
            feh.append(arrow(pts[0], pts[1]))
        log.append(("R9", "missing cofactor added as a dot", m, rid, model.name(m)))
    for p in mp.edges("fes", rid)[:1]:
        p.set("d", (p.get("d") or "") + to_d(fe))
    for p in mp.edges("fehs", rid)[:1]:
        p.set("d", (p.get("d") or "") + to_d(feh))


def complete_reactions(dr, reach=300):
    """D5: link a drawn reaction to the nearest node of each main metabolite it is missing."""
    mp, model, log = dr.mp, dr.model, dr.log
    for rid, g in list(mp.reactions().items()):
        if rid not in model.rxns or g.getparent() is None:
            continue
        st = model.rxns[rid]["stoich"]
        on = {mp.classes(n)[1] for n in mp.nodes("met") if rid in mp.classes(n)}
        if rid not in dr.placed:
            add_cofactor_dots(dr, rid, g, [m for m in st if m not in on and model.name(m) in COFACTORS])
        missing = [m for m in st if m not in on and model.name(m) not in COFACTORS]
        if not missing:
            continue
        mains = [m for m in st if model.name(m) not in COFACTORS]
        if mains and not (set(mains) & on) and rid not in dr.placed:
            # nothing it is drawn with is still part of it: redraw it as a missing reaction
            erase_reaction(dr, rid)
            draw_missing(dr, rid, compartments_on_map(mp))
            log.append(("D5", "reaction redrawn (none of its drawn metabolites is still part of it)", rid, "",
                        model.rxns[rid].get("name", "") or ""))
            continue
        rpos = translate(g)
        fe, feh = [], []
        for m in missing:
            cands = [(length(sub(translate(n), rpos)), n) for n in dr.met_nodes(m) if dr.is_main(n) and translate(n)]
            cands.sort(key=lambda t: t[0])
            if cands and cands[0][0] <= reach:
                n = cands[0][1]
                npos = translate(n)
            else:
                npos = None
                for ang in range(0, 360, 45):
                    u = (math.cos(math.radians(ang)), math.sin(math.radians(ang)))
                    c = add(rpos, u, STEP)
                    if dr.space.free((c[0] - 60, c[1] - 25, c[0] + 60, c[1] + 25), 0):
                        npos = c
                        break
                if npos is None:
                    log.append(("review", "main metabolite could not be placed", rid, m, model.name(m)))
                    continue
                n = new_main(m, model.name(m), [], *npos)
                append(mp.layer["nodes"], n)
                dr.space.add(node_box(mp, n))
            dr.link(n, rid)
            pts = [boundary(npos, rpos, "main"), boundary(rpos, npos, "rea")]
            fe.append(pts)
            rev = (model.rxns[rid].get("lower_bound") or 0) < 0
            if st[m] > 0 or rev:
                feh.append(arrow(pts[0], pts[1]))
            dr.band_add.append((rid, [npos, rpos]))
            log.append(("D5", "main metabolite linked to drawn reaction", rid, m,
                        model.name(m) + (" (existing node)" if cands and cands[0][0] <= reach else " (new node)")))
        for p in mp.edges("fes", rid)[:1]:
            p.set("d", (p.get("d") or "") + to_d(fe))
        for p in mp.edges("fehs", rid)[:1]:
            p.set("d", (p.get("d") or "") + to_d(feh))
        dr.placed.add(rid)


# ----------------------------------------------------------------------------- R8 gene stacks

def pack_genes(mp, log, only):
    """R8: close the gaps in the gene stacks of reactions that lost a gene (only)."""
    by_rid = collections.defaultdict(list)
    for n in mp.nodes("enz"):
        rids = [c for c in mp.classes(n)[2:] if is_rxn(c)]
        if len(rids) == 1 and rids[0] in only and translate(n):
            by_rid[rids[0]].append(n)
    for rid, boxes in by_rid.items():
        pos = [translate(n) for n in boxes]
        rows = collections.defaultdict(list)
        for n, p in zip(boxes, pos):
            rows[round(p[1] / 5)].append((p[0], n))
        cols = max(len(r) for r in rows.values())
        x0, y0 = min(p[0] for p in pos), min(p[1] for p in pos)
        order = [n for _, n in sorted(zip(((round(p[1]), p[0]) for p in pos), boxes), key=lambda t: t[0])]
        target = {id(n): (x0 + 60 * (i % cols), y0 + 30 * (i // cols)) for i, n in enumerate(order)}
        if all(abs(target[id(n)][0] - p[0]) < 1 and abs(target[id(n)][1] - p[1]) < 1 for n, p in zip(boxes, pos)):
            continue
        old_box = (x0, y0, max(p[0] for p in pos) + 60, max(p[1] for p in pos) + 30)
        for n in order:
            set_position(n, *target[id(n)])
        nb = (x0, y0, x0 + 60 * min(cols, len(order)), y0 + 30 * math.ceil(len(order) / cols))
        for p in mp.edges("ees", rid):
            sps = subpaths(p.get("d")) or []
            for sp in sps:
                for idx in (0, -1):
                    q = sp[idx]
                    if old_box[0] - 3 <= q[0] <= old_box[2] + 3 and old_box[1] - 3 <= q[1] <= old_box[3] + 3:
                        sp[idx] = (min(max(q[0], nb[0]), nb[2]), min(max(q[1], nb[1]), nb[3]))
            p.set("d", to_d(sps))
        log.append(("R8", "gene stack packed", rid, f"{len(order)} boxes", ""))


# ----------------------------------------------------------------------------- D6 empty space

def close_space(mp, before_boxes, log, keep=120, min_gap=300):
    """D6: inside each compartment box, close stripes that are empty now but held content before."""
    for comp, (bx0, by0, bx1, by1) in compartment_boxes(mp).items():
        group = next(g for g in mp.root.iter(SVG + "g") if g.get("class") == "compartment" and g.get("id") == comp)
        rect = group.find(SVG + "rect")
        grown = group.attrib.pop("data-grown", None) is not None
        inside = lambda b: bx0 <= (b[0] + b[2]) / 2 <= bx1 and by0 <= (b[1] + b[3]) / 2 <= by1  # noqa: E731
        now = [b for b in (node_box(mp, n) for n in mp.layer["nodes"]) if b and inside(b)]
        before = [b for b in before_boxes if inside(b)]
        # edges, arrowheads and the band count as content too (the band is 80 px wide)
        for layer, half in (("fes", 4), ("fehs", 4), ("ees", 4)):
            for path in mp.layer[layer]:
                for sp in subpaths(path.get("d") or "") or []:
                    now += [b for b in ((q[0] - half, q[1] - half, q[0] + half, q[1] + half) for q in sp) if inside(b)]
        for _, band_el in band_groups(mp):
            for sp in subpaths(band_el.get("d") or "") or []:
                for a, b2 in zip(sp, sp[1:]):
                    k = max(1, int(length(sub(b2, a)) / 40))
                    now += [bb for bb in ((x - 40, y - 40, x + 40, y + 40) for x, y in
                                          (add(a, sub(b2, a), i / k) for i in range(k + 1))) if inside(bb)]
        # the compartment's own heading stays where it is; content must not move under it
        heading = None
        for g in mp.root.iter(SVG + "g"):
            if g.get("class") == "compartment" and g.get("id") == comp:
                hb = [text_box(t) for t in g.iter(SVG + "text")]
                hb = [b for b in hb if b]
                if hb:
                    heading = (min(b[0] for b in hb), min(b[1] for b in hb), max(b[2] for b in hb), max(b[3] for b in hb))
        if not now:
            continue
        cuts = {}
        for axis, lo, hi in ((0, bx0, bx1), (1, by0, by1)):
            spans = sorted((b[axis], b[axis + 2]) for b in now)
            free, cur = [], lo
            for a, b in spans:
                if a > cur:
                    free.append((cur, a))
                cur = max(cur, b)
            free.append((cur, hi))
            reductions = []
            for a, b in free:
                if b - a < min_gap or (b - a) - (100 if a == lo or b == hi else keep) <= 0:
                    continue
                if not grown and not any(x[axis] < b and x[axis + 2] > a for x in before):
                    continue  # this space was already empty
                leading, trailing = a == lo, b == hi
                target = 100 if (leading or trailing) else keep
                if leading and heading:
                    # keep content in the heading's rows clear of the heading
                    other = 1 - axis
                    rows = [x for x in now if x[other] < heading[other + 2] and x[other + 2] > heading[other]]
                    if rows:
                        target = max(target, heading[axis + 2] - lo + 40)
                if (b - a) - target < 100:
                    continue  # not worth moving content for
                reductions.append((a, b, (b - a) - target, trailing))
            if reductions:
                cuts[axis] = reductions
        if not cuts:
            continue

        def shift(v, axis):
            d = 0.0
            for a, b, red, trailing in cuts.get(axis, []):
                if trailing:
                    continue
                if v >= b:
                    d += red
                elif v > a:
                    d += red * (v - a) / (b - a)
            return v - d

        def f(p):
            if not (bx0 <= p[0] <= bx1 and by0 <= p[1] <= by1):
                return p
            return (shift(p[0], 0), shift(p[1], 1))

        for n in list(mp.layer["nodes"]):
            p = translate(n)
            cls = mp.classes(n)
            if not p or not cls:
                continue
            ref = (p[0] + 30, p[1] + 15) if cls[0] == "enz" else p
            q = f(ref)
            if q != ref:
                set_position(n, q[0] - (ref[0] - p[0]), q[1] - (ref[1] - p[1]))
        for layer in ("fes", "fehs", "ees"):
            for path in mp.layer[layer]:
                sps = subpaths(path.get("d") or "")
                if sps:
                    path.set("d", to_d([[f(q) for q in sp] for sp in sps]))
        for _, band in band_groups(mp):
            sps = subpaths(band.get("d") or "")
            if sps:
                band.set("d", to_d([[f(q) for q in sp] for sp in sps]))
        dw = sum(r for a, b, r, t in cuts.get(0, []))
        dh = sum(r for a, b, r, t in cuts.get(1, []))
        rect.set("width", fmt(float(rect.get("width")) - dw))
        rect.set("height", fmt(float(rect.get("height")) - dh))
        log.append(("D6", "empty space closed", comp, f"{dw:.0f} px narrower, {dh:.0f} px lower", ""))
