"""Transport maps: the reactions that move metabolites across a membrane, one map per membrane.

An organelle is a box with the cytosol around it, so the box edge is its membrane; for the plasma membrane the
cell (the cytosol) is the box and the extracellular space is around it. Each reaction is a row across the
membrane: the metabolite outside, the reaction on the membrane, the metabolite inside (an exchanger has two such
lanes, one each way), each labelled on the side away from the membrane. Metabolites that cross along with it or
are used on the way (H+, Na+, ATP) are dots above the row, labelled above, on their compartment's side; the
genes are below the row, outside. Rows are grouped by transporter: the reaction's gene that transports most
across this membrane, when it does at least MIN_GROUP of them (the others are "Other genes", and reactions
without genes come last). Groups are sorted largest first and by metabolite within a group, and fill the box's
left edge, then its right edge. Long lists are shared between several boxes side by side, so that the map is
about ASPECT times as wide as it is tall; a membrane with more than MAX_ROWS reactions gets several maps, cut
between groups where possible. A compartment with fewer than TINY metabolites (Human-GEM's intermembrane space,
for the protons of proton symport) is drawn on the outer side of the membrane it takes part in.

Reactions between two compartments other than the cytosol (vesicular transport, lipid flip-flop) share one map,
grouped by compartment pair; each group has its own membrane line, with the first compartment on the left.

A model that files its transport reactions by membrane, in subsystems named "Transport [a, b]" (Yeast-GEM), gets a
map per such subsystem instead, named after it: a box for an organelle (or for the cell, when the other side is
the extracellular space or the cell envelope), and a membrane line, with a on the left, for a pair without the
cytosol.

The band of the transport subsystem runs along the membrane; reactions of other subsystems that cross the same
membrane are drawn too, as context, with a gap in the band (on maps of a per-membrane subsystem only those that
carry a metabolite across: the model files the others, such as lipid pseudo-reactions, elsewhere). The maps
are written whole from the model on every update, with their rows in subsystemSVG.tsv; mapedit leaves maps
marked data-transport alone.

usage: transport_map.py <model yml> <model dir> <template svg> <svg dir> <subsystemSVG.tsv> [--gene-label name|orf|both]
"""

import argparse
import collections
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import layout as L  # noqa: E402
import mapedit as me  # noqa: E402
import newmap  # noqa: E402

SUBSYSTEM = "Transport reactions"
PAIR = re.compile(r"^Transport\s*\[\s*(\w+)\s*,\s*(\w+)\s*\]$")  # a per-membrane transport subsystem
CYTOSOL, OUTSIDE = "c", "e"
OUTSIDE_WORDS = ("extracellular", "envelope", "wall", "periplasm")  # compartments outside the cell
TINY = 5                 # metabolites: a smaller compartment is drawn on the outer side of its membrane
MAX_ROWS = 600           # reactions on one map
D_OUT, D_IN = 280, 150   # px from the membrane to the metabolite node outside and inside
BOX_W = 800              # narrowest width of the inner compartment's box
WRAP, WRAP_LINES = 18, 6  # metabolite names wrap at about this many characters, in at most this many lines
LANE = 84                # between the lanes of an exchanger
DOT_Y = 64               # dots: this far above the reaction
DOT_RISE = 8             # each further dot this much higher
GENE_GAP = 14            # gene boxes: this far below the lowest lane
ROW_GAP = 18             # between rows
GROUP_GAP = 150          # between groups, with room for the group heading
TOP = 520                # first row
MARGIN = 260             # map edge to the outermost node or label
BETWEEN = 280            # between the reach of two boxes or membrane lines
MIN_GROUP = 3            # reactions a transporter needs for a group of its own
ASPECT = 1.4             # width / height the number of boxes aims at
COLOUR = "#64ff64"       # band and title
MIN_COLUMN = 10          # rows a column needs before the rows are shared with another box or line


def compartment(model, m):
    return model.mets.get(m, {}).get("compartment") or me.comp_of(m)


def name_of(c):
    return L.COMPARTMENT.get(c, c)


def slug(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def norm(subsystem):
    return re.sub(r"\s+", "", subsystem or "").lower()


def subsystems_of(model, rid):
    subs = model.rxns[rid].get("subsystem")
    return {norm(x) for x in (subs if isinstance(subs, list) else [subs]) if x}


def rows_of(model, outer, inner, subsystem=None):
    """The reactions with metabolites in `inner` and in the compartments `outer` (a set) only, each as (rid, lanes,
    dots, inner): lanes are (outside metabolite, inside metabolite) pairs (either can be None), dots the others.
    With a per-membrane subsystem, reactions of other subsystems only when a metabolite crosses."""
    out = []
    for rid, r in model.rxns.items():
        st = r["stoich"]
        if not st or L.artefact(model, rid):
            continue
        comps = {compartment(model, m) for m in st}
        if inner not in comps or not comps - {inner} or not comps - {inner} <= outer:
            continue
        side = {m: "in" if compartment(model, m) == inner else "out" for m in st}
        by_name = collections.defaultdict(dict)
        for m in st:
            by_name[model.name(m)][side[m]] = m
        crossing = [n for n, d in by_name.items() if len(d) == 2 and st[d["out"]] * st[d["in"]] < 0]
        if subsystem and not crossing and norm(subsystem) not in subsystems_of(model, rid):
            continue
        main = [n for n in crossing if n not in L.COFACTORS]
        if main:
            lanes = [(by_name[n]["out"], by_name[n]["in"]) for n in main]
        else:
            others = [m for m in st if model.name(m) not in L.COFACTORS]
            if others:  # converted on the way: the main metabolites on their own side, in pairs
                o = [m for m in others if side[m] == "out"]
                i = [m for m in others if side[m] == "in"]
                lanes = [(o[k] if k < len(o) else None, i[k] if k < len(i) else None) for k in range(max(len(o), len(i)))]
            else:  # only cofactors: those that cross are the main metabolites
                lanes = [(by_name[n]["out"], by_name[n]["in"]) for n in crossing] or \
                        [(m, None) if side[m] == "out" else (None, m) for m in st]
        on_lane = {m for lane in lanes for m in lane if m}
        out.append((rid, lanes, [m for m in st if m not in on_lane], inner))
    return out


def outside_like(c):
    return any(w in name_of(c).lower() for w in OUTSIDE_WORDS)


def per_membrane(model):
    """The maps of a model with per-membrane transport subsystems: one per subsystem."""
    spellings = collections.defaultdict(collections.Counter)  # compartment pair -> names of its subsystem
    for rid, r in model.rxns.items():
        subs = r.get("subsystem")
        for x in (subs if isinstance(subs, list) else [subs]):
            m = PAIR.match(x or "")
            if m:
                spellings[frozenset(m.groups())][x] += 1
    out = []
    for pair, names in spellings.items():
        sub = names.most_common(1)[0][0]
        a, b = PAIR.match(sub).groups()
        other = b if a == CYTOSOL else a
        if CYTOSOL in pair and not outside_like(other):
            kind, outer, inner = "box", {CYTOSOL}, other
        elif CYTOSOL in pair:
            kind, outer, inner = "box", {other}, CYTOSOL
        else:
            kind, outer, inner = "line", {a}, b
        rows = rows_of(model, outer, inner, sub)
        if rows:
            out.append(dict(stem=slug(sub), title=sub, kind=kind, outer=outer, inner=inner, rows=rows, subsystem=sub))
    return sorted(out, key=lambda x: x["stem"])


def membranes(model):
    """The maps to write, as dicts: file stem, title, kind ("box": rows grouped by transporter on the edges of the
    inner compartment's boxes; "line": the same on membrane lines; "pairs": one group per compartment pair), outer
    compartments, inner compartment, rows and the subsystem of the band."""
    if any(PAIR.match(x or "") for r in model.rxns.values()
           for x in (r.get("subsystem") if isinstance(r.get("subsystem"), list) else [r.get("subsystem")])):
        return per_membrane(model)
    count = collections.Counter(compartment(model, m) for m in model.mets)
    tiny = {c for c, n in count.items() if n < TINY}
    sets = collections.Counter()
    for rid, r in model.rxns.items():
        cs = frozenset(compartment(model, m) for m in r["stoich"])
        if len(cs) > 1 and not L.artefact(model, rid):
            sets[cs] += 1
    specs = []  # (inner, outer compartments, file stem, title)
    organelles = sorted({c for cs in sets if CYTOSOL in cs for c in cs} - {CYTOSOL, OUTSIDE} - tiny, key=name_of)
    for x in organelles:
        outer = {CYTOSOL} | {t for t in tiny if any(t in cs and x in cs for cs in sets)}
        specs.append((x, outer, "transport_" + slug(name_of(x)), "Transport: " + name_of(x)))
    if frozenset((CYTOSOL, OUTSIDE)) in sets:
        specs.append((CYTOSOL, {OUTSIDE}, "transport_plasma_membrane", "Transport: Plasma membrane"))

    def covered(cs):
        return any(inner in cs and cs - {inner} and cs - {inner} <= outer for inner, outer, _, _ in specs)
    out = []
    for inner, outer, stem, title in specs:
        rows = rows_of(model, outer, inner)
        if rows:
            out.append(dict(stem=stem, title=title, kind="box", outer=outer, inner=inner, rows=rows, subsystem=SUBSYSTEM))
    pairs = sorted({tuple(sorted(cs, key=lambda c: (c != OUTSIDE, name_of(c)))) for cs in sets
                    if len(cs) == 2 and not covered(cs)})
    rows = [row for a, b in pairs for row in rows_of(model, {a}, b)]
    if rows:
        out.append(dict(stem="transport_between_organelles", title="Transport: between organelles", kind="pairs",
                        outer=None, inner=None, rows=rows, subsystem=SUBSYSTEM))
    return out


def group_keys(model, rows):
    """Reaction -> its transporter group: the gene with most reactions across this membrane among the
    reaction's genes, if it has at least MIN_GROUP; "Other genes"; or "" for no gene."""
    count = collections.Counter(g for row in rows for g in model.rxns[row[0]]["genes"])
    out = {}
    for row in rows:
        genes = model.rxns[row[0]]["genes"]
        if not genes:
            out[row[0]] = ""
            continue
        top = min(genes, key=lambda g: (-count[g], model.symbol.get(g) or g))
        out[row[0]] = (model.symbol.get(top) or top) if count[top] >= MIN_GROUP else "Other genes"
    return out


def label_lines(model, m):
    return L.wrap(model.name(m), WRAP, WRAP_LINES)


def label_width(model, m):
    return max(L.label_width(line, 18) for line in label_lines(model, m))


def node_half(model, m):
    return max(22, 10 * len(label_lines(model, m)) + 4)


def main_node(model, m, rid, x, y, left):
    """A main metabolite node with its name beside it, on the left or the right, in lines of about WRAP
    characters."""
    n = L.new_main(m, model.name(m), [rid], x, y)
    for old in [c for c in n if c.get("class") == "lbl"]:
        n.remove(old)
    lines = label_lines(model, m)
    style = {"class": "lbl", "font-family": "sans-serif", "font-size": "18", "font-weight": "700",
             "text-anchor": "end" if left else "start"}
    x0 = "-38" if left else "38"
    if len(lines) == 1:
        n.append(L.text_el(lines[0], x=x0, y="6", **style))
    else:
        lbl = L.element("g", **style)
        top = 6 - 10 * (len(lines) - 1)
        for k, line in enumerate(lines):
            lbl.append(L.text_el(line, x=x0, y=L.fmt(top + 20 * k)))
        n.append(lbl)
    return n


def row_extent(model, lanes, dots, genes):
    """(above, below): px the row takes above its first lane and below its last lane."""
    above = max(node_half(model, m) for lane in lanes[:1] for m in lane if m)
    if dots:
        per_side = max(collections.Counter(me.comp_of(m) for m in dots).values())
        above = max(above, (len(lanes) - 1) * LANE / 2 + DOT_Y + DOT_RISE * (per_side - 1) + 30)
    below = max(node_half(model, m) for lane in lanes[-1:] for m in lane if m)
    if genes:
        below = max(below, GENE_GAP + 30 * ((len(genes) + 2) // 3))
    return above, below


def row_height(model, row):
    above, below = row_extent(model, row[1], row[2], model.rxns[row[0]]["genes"])
    return above + LANE * (len(row[1]) - 1) + below + ROW_GAP


def heading(text, x, y, size, anchor="middle", fill="#00000"):
    return L.text_el(text, **{"y": "0", "transform": L.matrix(x, y), "fill": fill, "font-family": "sans-serif",
                              "font-size": str(size), "font-weight": "700", "text-anchor": anchor})


def draw_row(dr, row, y0, membrane, out_dir, context):
    """Draw one reaction with its first lane at height y0; out_dir is -1 when outside is to the left."""
    rid, lanes, dots, inner = row
    mp, model = dr.mp, dr.model
    st = model.rxns[rid]["stoich"]
    ys = [y0 + LANE * k for k in range(len(lanes))]
    rpos = (membrane, sum(ys) / len(ys))
    links = []
    for (mo, mi), y in zip(lanes, ys):
        for m, dx in ((mo, out_dir * D_OUT), (mi, -out_dir * D_IN)):
            if not m:
                continue
            p = (membrane + dx, y)
            n = main_node(model, m, rid, p[0], p[1], dx < 0)
            L.append(mp.layer["nodes"], n)
            links.append((p, "main", 1 if st[m] > 0 else -1))
    L.append(mp.layer["nodes"], L.new_reaction(rid, *rpos))
    # dots: one row per side, substrates nearest the membrane, each a little higher than the one before so
    # that their edges fan out; labelled above (the edges come from below)
    for inside, sign in ((False, out_dir), (True, -out_dir)):
        x = membrane + sign * 40
        mine = [m for m in dots if (compartment(model, m) == inner) == inside]
        for k, m in enumerate(sorted(mine, key=lambda m: (st[m] > 0, model.name(m)))):
            name = dr.side_name(m)
            w = max(L.label_width(name, 18), 14)
            x += sign * w / 2
            p = (x, rpos[1] - DOT_Y - DOT_RISE * k - (len(lanes) - 1) * LANE / 2)
            dot = L.new_side(m, name, [rid], p[0], p[1], False)
            t = dot.find(L.SVG + "text")
            t.set("text-anchor", "middle")
            t.set("x", "0")
            t.set("y", "-14")
            L.append(mp.layer["nodes"], dot)
            links.append((p, "side", 1 if st[m] > 0 else -1))
            x += sign * (w / 2 + 24)
    rev = (model.rxns[rid].get("lower_bound") or 0) < 0
    fe, feh = [], []
    for npos, kind, role in links:
        pts = [L.boundary(npos, rpos, kind), L.boundary(rpos, npos, "rea")]
        fe.append(pts)
        if role > 0 or rev:
            feh.append(L.arrow(pts[0], pts[1]))
    L.append(mp.layer["fes"], L.element("path", **{"class": "fe " + rid, "d": L.to_d(fe)}))
    L.append(mp.layer["fehs"], L.element("path", **{"class": "feh " + rid, "d": L.to_d(feh)}))
    # genes below the lowest lane, between the outer node and the membrane, in rows of three
    genes = sorted(model.rxns[rid]["genes"], key=lambda g: model.symbol.get(g) or g)
    if genes:
        cols = min(3, len(genes))
        top = ys[-1] + GENE_GAP
        x0 = membrane - 46 - 60 * cols if out_dir < 0 else membrane + 46
        for k, g in enumerate(genes):
            L.append(mp.layer["nodes"], L.new_gene(g, model.symbol.get(g, g), [rid], x0 + 60 * (k % cols),
                                                   top + 30 * (k // cols)))
        box = (x0, top, x0 + 60 * cols, top + 30 * ((len(genes) + cols - 1) // cols))
        centre = ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2)
        start = L.boundary(rpos, centre, "rea")
        end = (min(max(start[0], box[0]), box[2]), min(max(start[1], box[1]), box[3]))
        L.append(mp.layer["ees"], L.element("path", **{"class": "ee " + rid, "d": L.to_d([[start, end]])}))


def columns(items, n):
    """Cut the sequence of (kind, height, payload) items into n columns of about equal height, in order. Cuts go
    between groups where one is near; inside a group they leave at least MIN_GROUP rows on each side, and the
    column that starts inside a group repeats its heading."""
    total = sum(h for _, h, _ in items)
    pos, acc = [], 0.0
    for kind, h, _ in items:
        pos.append(acc)
        acc += h
    group_at, rows_before, g, k = [], [], -1, 0
    for kind, _, _ in items:
        if kind == "head":
            g, k = g + 1, 0
        group_at.append(g)
        rows_before.append(k)
        k += kind == "row"
    size = collections.Counter(group_at[i] for i, (kind, _, _) in enumerate(items) if kind == "row")
    rows_after = [size[group_at[i]] - rows_before[i] for i in range(len(items))]
    cuts, last = [], 0
    for c in range(1, n):
        target = total * c / n
        best = None
        for i in range(last + 1, len(items)):
            kind = items[i][0]
            if kind == "row" and (rows_before[i] < MIN_GROUP or rows_after[i] < MIN_GROUP):
                continue
            score = abs(pos[i] - target) + (0 if kind == "head" else 0.12 * total / n)
            if best is None or score < best[0]:
                best = (score, i)
        if best is None:
            break
        cuts.append(best[1])
        last = best[1]
    cols = []
    for a, b in zip([0] + cuts, cuts + [len(items)]):
        col = items[a:b]
        if col and col[0][0] == "row":
            head = next(items[i] for i in range(a, -1, -1) if items[i][0] == "head")
            col = [("head", head[1], dict(head[2], label=head[2]["label"] + ", continued"))] + col
        cols.append(col)
    return cols + [[] for _ in range(n - len(cols))]


def items_of(model, rows, kind):
    """The rows as a sequence of group headings and rows: by transporter ("box", "line") or compartment pair
    ("pairs")."""
    if kind != "pairs":
        key = group_keys(model, rows)
        group = {row[0]: key[row[0]] for row in rows}
    else:
        outer = {row[0]: next(compartment(model, m) for m in model.rxns[row[0]]["stoich"]
                              if compartment(model, m) != row[3]) for row in rows}
        group = {row[0]: (outer[row[0]], row[3]) for row in rows}
    groups = collections.defaultdict(list)
    for row in rows:
        groups[group[row[0]]].append(row)
    for g in groups.values():
        g.sort(key=lambda r: (min(model.name(m).lower() for lane in r[1] for m in lane if m), r[0]))
    if kind != "pairs":
        order = sorted(groups, key=lambda k: (k == "", k == "Other genes", -len(groups[k]), k))
    else:
        order = sorted(groups, key=lambda k: (k[0] != OUTSIDE, name_of(k[0]), name_of(k[1])))
    items = []
    for k in order:
        label = (k or "No gene") if kind != "pairs" else f"{name_of(k[0])} | {name_of(k[1])}"
        items.append(("head", GROUP_GAP, {"label": label, "n": len(groups[k]), "pair": k if kind == "pairs" else None}))
        items += [("row", row_height(model, r), r) for r in groups[k]]
    return items


def write_map(model, template, path, title, kind, outer, inner, items, subsystem, colour=None):
    """One transport map from its items; returns (boxes or membrane lines, reactions, context reactions)."""
    rows = [x for k, _, x in items if k == "row"]
    total = sum(h for _, h, _ in items)
    dots_w = max([sum(max(L.label_width(model.name(m), 18), 14) + 24 for m in r[2]) for r in rows] + [0])
    out_w = max([label_width(model, lane[0]) for r in rows for lane in r[1] if lane[0]] + [0])
    in_w = max([label_width(model, lane[1]) for r in rows for lane in r[1] if lane[1]] + [0])
    reach = max(D_OUT + 38 + out_w, 40 + dots_w)
    reach_in = max(D_IN + 38 + in_w, 40 + dots_w)
    box_w = max(BOX_W, 2 * reach_in + 40)
    if kind == "box":  # columns in pairs, on the left and right edge of each box
        unit, per = box_w + 2 * reach + BETWEEN, 2
    else:  # one membrane line per column, outside on its left
        unit, per = reach + reach_in + BETWEEN, 1
    best = None
    for k in range(1, 13):  # as many boxes or lines as bring the map closest to ASPECT
        if k > 1 and len(rows) < MIN_COLUMN * per * k:
            break
        w = k * unit - BETWEEN + 2 * MARGIN
        h = TOP + total / (per * k) + 300
        score = abs(w / h - ASPECT) / ASPECT
        if best is None or score < best[0]:
            best = (score, k, w)
    _, n, width = best
    cols = columns(items, per * n)
    height = TOP + max(sum(h for _, h, _ in c) for c in cols) + 300
    newmap.blank(template, title, path, width, height, colour or COLOUR)
    mp = me.Map(path)
    mp.root.attrib.pop("data-new", None)
    mp.root.set("data-transport", f"{''.join(sorted(outer))} {inner}" if kind != "pairs" else "pairs")
    mp.root.set("data-modelversion", model.version)
    band = None
    for g, p in L.band_groups(mp):
        g.set("id", subsystem)
        band = p
    dr = L.Drawer(mp, model, [])
    context = {row[0] for row in rows if norm(subsystem) not in subsystems_of(model, row[0])}
    notes = L.element("g", id="transport-groups", **{"class": "note"})
    band_lines, lines, places = [], [], []
    for c, col in enumerate(cols):
        if kind == "box":
            x0 = MARGIN + reach + (c // 2) * unit
            membrane, s = (x0, -1) if c % 2 == 0 else (x0 + box_w, 1)
            if c % 2 == 0:
                places.append((x0, x0 + box_w))
        else:
            membrane, s = MARGIN + reach + c * unit, -1
        y, run, line = TOP - GROUP_GAP + 40, [], None
        if kind == "line" and col:  # the two compartments over each membrane line
            (a,), b = tuple(outer), inner
            for c_, x_, anchor in ((a, membrane - 46, "end"), (b, membrane + 46, "start")):
                parts = L.wrap(name_of(c_), 16, 3)
                for k, part in enumerate(parts):
                    notes.append(heading(part, x_, TOP - 150 - 46 * (len(parts) - 1 - k), 40, anchor))
        for what, h, x in col:
            if what == "head":
                if run:
                    band_lines.append(run)
                    run = []
                if line and kind == "pairs":
                    lines.append(line)
                    line = None
                if kind != "pairs":
                    notes.append(heading(f"{x['label']} ({x['n']})", membrane + s * 46, y + GROUP_GAP - 50, 30,
                                         "end" if s < 0 else "start", "#555"))
                else:
                    a, b = x["pair"]
                    cont = ", continued" if x["label"].endswith(", continued") else ""
                    notes.append(heading(f"{name_of(a)} ({x['n']}{cont})", membrane - 46, y + GROUP_GAP - 50, 30,
                                         "end", "#555"))
                    notes.append(heading(name_of(b), membrane + 46, y + GROUP_GAP - 50, 30, "start", "#555"))
                y += GROUP_GAP
                continue
            above, below = row_extent(model, x[1], x[2], model.rxns[x[0]]["genes"])
            y += above
            draw_row(dr, x, y, membrane, s, context)
            mid = y + LANE * (len(x[1]) - 1) / 2
            line = (line or [(membrane, y - above)])[:1] + [(membrane, y + LANE * (len(x[1]) - 1) + below)]
            if x[0] in context:
                if run:
                    band_lines.append(run)
                run = []
            else:
                run = (run or [(membrane, mid - 30)])[:1] + [(membrane, mid + 30)]
            y += LANE * (len(x[1]) - 1) + below + ROW_GAP
        if run:
            band_lines.append(run)
        if line:
            lines.append(line)
    if band is not None:  # along the membrane, where the rows are transport reactions
        band.set("d", L.to_d(band_lines))
    top, bottom = TOP - 260, height - 180
    if kind == "box":
        for b, (x0, x1) in enumerate(places):
            g = L.element("g", id=name_of(inner) + (f" ({b + 1})" if len(places) > 1 else ""), **{"class": "compartment"})
            g.append(L.element("rect", **{"class": "shape", "x": L.fmt(x0), "y": L.fmt(top), "width": L.fmt(box_w),
                                          "height": L.fmt(bottom - top), "rx": "150", "ry": "150", "fill": "none",
                                          "stroke": "#000", "stroke-width": "22"}))
            g.append(heading(name_of(inner), (x0 + x1) / 2, top + 120, 72))
            mp.layer["fes"].addprevious(g)
            g.tail = "\n  "
        main_outer = CYTOSOL if CYTOSOL in outer else sorted(outer)[0]
        h = L.element("g", id=name_of(main_outer), **{"class": "compartment"})
        h.append(heading(name_of(main_outer), MARGIN + reach / 2, top + 120, 72))
        mp.layer["fes"].addprevious(h)
        h.tail = "\n  "
    else:  # the membrane lines: one per column, or per group of a compartment pair
        notes.append(L.element("path", **{"d": L.to_d(lines), "fill": "none", "stroke": "#000", "stroke-width": "22"}))
    mp.layer["fes"].addprevious(notes)
    notes.tail = "\n  "
    for t in mp.root.iter(L.SVG + "text"):  # the title over the boxes
        if t.getparent().get("class") == "subsystem":
            t.set("transform", L.matrix(width / 2, 200))
            t.set("text-anchor", "middle")
    mp.write(path)
    return n, len(rows), len(context)


def part_label(items):
    """What one of several maps of a membrane holds: its groups, or the first two and how many more."""
    labels = list(dict.fromkeys(x["label"].replace(", continued", "") for what, _, x in items if what == "head"))
    labels = [x if x != "No gene" else "no gene" for x in labels]
    if len(labels) <= 3:
        return ", ".join(labels)
    return f"{labels[0]}, {labels[1]} and {len(labels) - 2} more"


def write_all(model, template, svg_dir, table):
    """Every transport map of the model into svg_dir, and their rows in the subsystem map table; transport maps
    written before and not now are removed."""
    written = []
    for spec in membranes(model):
        stem, title, kind, outer, inner, rows = (spec[k] for k in ("stem", "title", "kind", "outer", "inner", "rows"))
        items = items_of(model, rows, kind)
        parts = columns(items, math.ceil(len(rows) / MAX_ROWS)) if len(rows) > MAX_ROWS else [items]
        for k, part in enumerate(parts):
            name = f"{stem}_{k + 1}" if len(parts) > 1 else stem
            ttl = f"{title} {k + 1} ({part_label(part)})" if len(parts) > 1 else title
            n, nrows, nctx = write_map(model, template, os.path.join(svg_dir, name + ".svg"), ttl, kind, outer, inner,
                                       part, spec["subsystem"])
            written.append((name + ".svg", ttl, spec["subsystem"]))
            groups = sum(1 for what, _, _ in part if what == "head")
            print(f"{name}.svg: {nrows} reactions ({nctx} of other subsystems), {groups} groups, "
                  f"{n} {'box(es)' if kind == 'box' else 'membrane line(s)'}", flush=True)
    names = {f for f, _, _ in written}
    for f in sorted(os.listdir(svg_dir)):
        p = os.path.join(svg_dir, f)
        if f.startswith("transport_") and f.endswith(".svg") and f not in names and "data-transport=" in open(p).read(2000):
            os.remove(p)
            print("removed", f)
    lines = open(table).read().rstrip("\n").split("\n")
    keep = [x for x in lines if not re.match(r"[^\t]*\t[^\t]*\ttransport_[a-z0-9_]+\.svg\s*$", x)]
    keep += [f"{sub}\t{t}\t{f}" for f, t, sub in written]
    with open(table, "w") as fh:
        fh.write("\n".join(keep) + "\n")
    return written


def main():
    ap = argparse.ArgumentParser(description="Write every transport map of a model, with their rows in subsystemSVG.tsv.")
    for arg in ("model_yml", "model_dir", "template", "svg_dir", "table"):
        ap.add_argument(arg)
    ap.add_argument("--gene-label", choices=["name", "orf", "both"], default="name",
                    help="gene boxes show the gene name, the gene id, or the name over the id")
    a = ap.parse_args()
    L.GENE_LABEL = a.gene_label
    yml, model_dir, template, svg_dir, table = a.model_yml, a.model_dir, a.template, a.svg_dir, a.table
    model = me.Model(yml, model_dir, [])
    written = write_all(model, template, svg_dir, table)
    print(f"{len(written)} transport maps; rows written to {table}")


if __name__ == "__main__":
    main()
