"""D7: use KEGG pathway layouts (KGML) as a reference for placing reactions.

D7b  KEGG-shaped blocks: added reactions that a KEGG map shows start from the KEGG arrangement
     (scaled), one block per connected module of a KEGG map; a compound shared within a module is one
     node. D9 (placement.py) places the blocks.
D8   KEGG relayout: a map's cytosol redrawn as its KEGG map lays it out (opt-in per map).
"""

import collections
import xml.etree.ElementTree as ET

import layout as L

SCALE = 4.0  # KEGG pixels to map pixels in the bottom section


class Kgml:
    def __init__(self, path):
        import os
        load_anomers(os.path.dirname(path))
        root = ET.parse(path).getroot()
        self.name = root.get("name", "")
        self.title = root.get("title", "")
        self.compounds = {}  # entry id -> (KEGG compound id, (x, y))
        lines = {}
        for e in root.findall("entry"):
            g = e.find("graphics")
            if e.get("type") == "compound" and g is not None:
                cid = e.get("name", "").split()[0].replace("cpd:", "").replace("gl:", "")
                self.compounds[e.get("id")] = (cid, (float(g.get("x")), float(g.get("y"))))
            elif e.get("type") == "reaction" and g is not None and g.get("coords"):
                c = [float(v) for v in g.get("coords").split(",")]
                lines[e.get("id")] = list(zip(c[0::2], c[1::2]))
            elif e.get("type") == "reaction" and g is not None and g.get("x"):
                # the enzyme box on KEGG's reaction line: a two-point line around it, so its midpoint is the box
                x, y = float(g.get("x")), float(g.get("y"))
                lines[e.get("id")] = [(x - 1.0, y), (x + 1.0, y)]
        self.reactions = collections.defaultdict(list)  # KEGG reaction id -> [(subs, prods, line)]
        for r in root.findall("reaction"):
            subs = [s.get("id") for s in r.findall("substrate") if s.get("id") in self.compounds]
            prods = [p.get("id") for p in r.findall("product") if p.get("id") in self.compounds]
            for name in r.get("name", "").split():
                self.reactions[name.replace("rn:", "")].append((subs, prods, lines.get(r.get("id"))))
        self.by_compound = collections.defaultdict(list)
        for eid, (cid, pos) in self.compounds.items():
            self.by_compound[cid].append(pos)


def kegg_ids(model, rid):
    return [k for k in model.rxn_kegg.get(rid, "").split(";") if k]


ANOMERS = {}  # KEGG compound id -> ids of the same compound with or without an alpha-/beta- anomer prefix


def load_anomers(kegg_dir):
    """Read compound_list.tsv (KEGG /list/compound) once; alpha-X and beta-X are taken as X."""
    import os
    import re
    path = os.path.join(kegg_dir, "compound_list.tsv")
    if ANOMERS or not os.path.exists(path):
        return
    by_name, base = {}, {}
    for line in open(path):
        cid, _, names = line.rstrip("\n").partition("\t")
        first = names.split(";")[0].strip()
        by_name.setdefault(first.lower(), cid)
        m = re.match(r"(alpha|beta)-(.+)", first)
        if m:
            base[cid] = m.group(2).lower()
    for cid, name in base.items():
        b = by_name.get(name)
        if b:
            ANOMERS.setdefault(b, set()).add(cid)
            ANOMERS.setdefault(cid, set()).add(b)
    for b in list(ANOMERS):  # alpha and beta forms of one base are equivalent too
        for x in list(ANOMERS[b]):
            ANOMERS[x] |= ANOMERS[b] - {x}


def met_kegg(model, m):
    ks = {k for k in model.met_kegg.get(m, "").split(";") if k}
    return ks | {a for k in ks for a in ANOMERS.get(k, ())}


def match(model, kgmls, rid):
    """(kgml, subs entry ids, prods entry ids, line, orientation flag) for a reaction, or None."""
    for kg in kgmls:
        for k in kegg_ids(model, rid):
            for subs, prods, line in kg.reactions.get(k, []):
                return kg, subs, prods, line
    return None


def assign(model, kg, subs_e, prods_e, rid, mains_sub, mains_prod):
    """Map our main metabolites to KGML compound entries (either orientation); {met: entry id}."""
    def pairs(ours, entries):
        out = {}
        for exact in (True, False):  # exact KEGG ids first, then anomers
            for m in ours:
                if m in out:
                    continue
                ks = {k for k in model.met_kegg.get(m, "").split(";") if k} if exact else met_kegg(model, m)
                e = next((e for e in entries if kg.compounds[e][0] in ks and e not in out.values()), None)
                if e:
                    out[m] = e
        return out
    a = {**pairs(mains_sub, subs_e), **pairs(mains_prod, prods_e)}
    b = {**pairs(mains_sub, prods_e), **pairs(mains_prod, subs_e)}
    return a if len(a) >= len(b) else b


def mid_line(line, a=None, b=None):
    if line and len(line) >= 2:
        return L.point_at(line, L.poly_length(line) / 2)[0]
    return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2) if a and b else None


# ----------------------------------------------------------------------------- D7b

def bottom_blocks(dr, rids):
    """Split rids into KEGG-shaped blocks and the rest. Returns ([(kg, [(rid, amap, se, pe, line)])], rest)."""
    model = dr.model
    groups = collections.defaultdict(list)
    rest = []
    for rid in rids:
        st = model.rxns[rid]["stoich"]
        subs = [m for m in st if st[m] < 0 and model.name(m) not in L.COFACTORS]
        prods = [m for m in st if st[m] > 0 and model.name(m) not in L.COFACTORS]
        m = match(model, dr.kegg, rid)
        if not m:
            rest.append(rid)
            continue
        kg, se, pe, line = m
        amap = assign(model, kg, se, pe, rid, subs, prods)
        if not amap:
            rest.append(rid)
            continue
        groups[kg.name].append((rid, subs, prods, amap, line))
    kgs = {kg.name: kg for kg in dr.kegg}
    blocks = []
    for k, items in groups.items():
        for comp in modules(items):
            blocks.append((kgs[k], comp))
    blocks.sort(key=lambda b: -len(b[1]))
    return blocks, rest


def modules(items):
    """Split a KEGG map's items into connected modules (items sharing a compound entry)."""
    parent = list(range(len(items)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen = {}
    for i, it in enumerate(items):
        for e in it[3].values():
            if e in seen:
                parent[find(i)] = find(seen[e])
            else:
                seen[e] = i
    comps = collections.defaultdict(list)
    for i, it in enumerate(items):
        comps[find(i)].append(it)
    return list(comps.values())


def edge_lengths(kg, items):
    """KEGG (dx, dy) between the substrate and product entries of each item."""
    out = []
    for _, subs, prods, amap, _ in items:
        a = [kg.compounds[amap[x]][1] for x in subs if x in amap]
        b = [kg.compounds[amap[x]][1] for x in prods if x in amap]
        if a and b:
            out.append((abs(b[0][0] - a[0][0]), abs(b[0][1] - a[0][1])))
    return out


def module_scale(kg, items, longest=400):
    """SCALE, reduced so that no reaction of the module is drawn longer than `longest` px."""
    lens = [L.length(e) for e in edge_lengths(kg, items)]
    return min(SCALE, longest / max(lens)) if lens and max(lens) > 0 else SCALE


def block_size(kg, items):
    pts = [kg.compounds[e][1] for _, _, _, amap, _ in items for e in amap.values()]
    sc = module_scale(kg, items)
    w = sc * (max(p[0] for p in pts) - min(p[0] for p in pts)) + 300
    h = sc * (max(p[1] for p in pts) - min(p[1] for p in pts)) + 400
    return w, h


# ----------------------------------------------------------------------------- D8

def grid_spread(P, size, linked, between=200, align=15, reach=1500):
    """Open up crowded places without breaking KEGG's rows and columns: along each axis the distinct
    coordinates (columns, then rows) keep their order and alignment; a column moves further out only as far
    as needed to clear the labels (and, for compounds joined by a reaction, the reaction between them) of
    compounds in earlier columns that share its row."""
    for axis in (0, 1):
        other = 1 - axis
        keys = sorted(P, key=lambda k: P[k][axis])
        cols = []
        for k in keys:
            if cols and P[k][axis] - cols[-1][0] <= align:
                cols[-1][1].append(k)
            else:
                cols.append([P[k][axis], [k]])
        new = []
        for i, (orig, ks) in enumerate(cols):
            pos = orig if not new else new[-1] + (orig - cols[i - 1][0])
            for j in range(i - 1, -1, -1):
                if orig - cols[j][0] > reach:
                    break
                for b in ks:
                    for a in cols[j][1]:
                        half = (size[a][other] + size[b][other]) / 2
                        if abs(P[a][other] - P[b][other]) < half + 20:  # they share a row: no overlap allowed
                            need = (size[a][axis] + size[b][axis]) / 2 + 40
                            if frozenset((a, b)) in linked:
                                need += between
                            pos = max(pos, new[j] + need)
            new.append(pos)
        for (orig, ks), pos in zip(cols, new):
            for k in ks:
                P[k] = (pos, P[k][1]) if axis == 0 else (P[k][0], pos)


def orthogonal(pos, rpos, a, b, straight=40):
    """KEGG-style edges: one bend per edge, the segment at the reaction along the reaction's main direction
    (substrate to product). {metabolite: bend point}; edges already (nearly) straight get none."""
    horizontal = abs(b[0][0] - a[0][0]) >= abs(b[0][1] - a[0][1]) if a and b else True
    out = {}
    for m, p in pos.items():
        if abs(p[0] - rpos[0]) > straight and abs(p[1] - rpos[1]) > straight:
            out[m] = (p[0], rpos[1]) if horizontal else (rpos[0], p[1])
    return out


def items_for(dr, kg, rids):
    """(rid, subs, prods, amap, line) for the reactions kg shows; and the rest."""
    model = dr.model
    items, rest = [], []
    for rid in rids:
        st = model.rxns[rid]["stoich"]
        subs = [m for m in st if st[m] < 0 and model.name(m) not in L.COFACTORS]
        prods = [m for m in st if st[m] > 0 and model.name(m) not in L.COFACTORS]
        m = match(model, [kg], rid)
        amap = assign(model, kg, m[1], m[2], rid, subs, prods) if m else {}
        if m and amap:
            items.append((rid, subs, prods, amap, m[3]))
            continue
        # no shared reaction id: a KEGG reaction joining a substrate and a product of ours
        best = None
        for entries in kg.reactions.values():
            for se, pe, line in entries:
                a = assign(model, kg, se, pe, rid, subs, prods)
                if any(x in a for x in subs) and any(x in a for x in prods):
                    exact = sum(kg.compounds[e][0] in model.met_kegg.get(x, "").split(";") for x, e in a.items())
                    if best is None or (len(a), exact) > (len(best[0]), best[2]):
                        best = (a, line, exact)
        if best:
            items.append((rid, subs, prods, best[0], best[1]))
        else:
            rest.append(rid)
    return items, rest


def relayout(dr, kg, scope, comp="c", margin=150, top=850, min_scale=3.0, max_fill=4.0, lane_gap=200,
             min_share=0.5):
    """D8: redraw the map's reactions of one compartment as the KEGG map lays them out.

    The compartment's area is the canvas outside the compartment boxes (the cytosol on subsystem maps).
    Its reactions, drawn and missing, are erased; those KEGG shows are drawn at the KEGG positions
    (one scale, filling the area, larger where compounds would crowd), one node per KEGG compound entry; the others are
    left to D3, which attaches them next to the new layout or puts them in the bottom section."""
    mp, model, log = dr.mp, dr.model, dr.log
    boxes = list(L.compartment_boxes(mp).values())
    cw, ch = L.canvas(mp)[2:]

    def outside(p):
        return p and not any(b[0] <= p[0] <= b[2] and b[1] <= p[1] <= b[3] for b in boxes)

    def in_comp(rid):  # with one area, only the main metabolites need to be in the compartment
        st = model.rxns.get(rid, {}).get("stoich", {})
        mets = [m for m in st if not dr.one_area or model.name(m) not in L.COFACTORS]
        return st and mets and all(L.comp_of(m) == comp for m in mets)
    drawn = [rid for rid, g in mp.reactions().items() if in_comp(rid) and outside(L.translate(g))]
    missing = [rid for rid in scope.missing(set(mp.reactions())) if in_comp(rid)] if scope else []
    items, rest = items_for(dr, kg, drawn + missing)
    share = len(items) / max(1, len(drawn) + len(missing))
    if len(items) < 5 or share < min_share:
        log.append(("review", "KEGG relayout skipped: too few reactions on the KEGG map", kg.name,
                    f"{len(items)} of {len(drawn) + len(missing)}", ""))
        return False
    for rid in drawn:
        L.erase_reaction(dr, rid)
    for n in list(mp.nodes("met")) + list(mp.nodes("enz")):
        if outside(L.translate(n)) and not [c for c in mp.classes(n)[2:] if L.is_rxn(c)]:
            mp.remove(n)
    dr.space = L.Space(mp)

    # area: left of the boxes that reach into the content band, below the headings (just below the title
    # on a map without compartment headings)
    heads = [L.text_box(t) for g in mp.root.iter(L.SVG + "g") if g.get("class") in ("subsystem", "compartment")
             for t in g.iter(L.SVG + "text")]
    if not [g for g in mp.root.iter(L.SVG + "g") if g.get("class") == "compartment"] and any(heads):
        top = min(top, max(h[3] for h in heads if h) + 250)
    left = margin
    for g in mp.root.iter(L.SVG + "g"):
        if g.get("class") != "compartment" or g.find(L.SVG + "rect") is not None:
            continue
        if g.get("id") == L.COMPARTMENT.get(comp):  # the compartment's own heading: start just below it
            hb = [L.text_box(t) for t in g.iter(L.SVG + "text")]
            if any(hb):
                top = min(top, max(h[3] for h in hb if h) + 200)
        for path in g.iter(L.SVG + "path"):  # a divider line of another compartment (extracellular)
            pts = [q for sp in (L.subpaths(path.get("d") or "") or []) for q in sp]
            if pts and max(q[0] for q in pts) - min(q[0] for q in pts) < 50 and max(q[0] for q in pts) < cw / 2:
                left = max(left, max(q[0] for q in pts) + margin)
    area = (left, top, cw - margin, ch - margin)  # compartment boxes in the way move aside afterwards
    pts = [kg.compounds[e][1] for it in items for e in it[3].values()]
    kx0, ky0 = min(p[0] for p in pts), min(p[1] for p in pts)
    kw = max(p[0] for p in pts) - kx0
    aw = area[2] - area[0] - 300
    # one scale for both axes, so the layout keeps KEGG's proportions: as large as the area allows, up to
    # max_fill; crowded places are opened up afterwards (grid_spread)
    ents = sorted({e for it in items for e in it[3].values()})
    s_fill = aw / max(kw, 1)
    sx = sy = max(min_scale, min(max_fill, s_fill))  # grid_spread opens up the crowded places
    scale = (sx, sy)

    def place(q):
        return (area[0] + 150 + sx * (q[0] - kx0), area[1] + 150 + sy * (q[1] - ky0))
    # crowded places open up; rows and columns stay aligned
    P = {e: place(kg.compounds[e][1]) for e in ents}
    label = {}
    for it in items:
        for met, e in it[3].items():
            lines = L.wrap(model.name(met))
            label.setdefault(e, (max(L.label_width(x, 18) for x in lines), 24 * len(lines) + 20))
    linked = {frozenset((it[3][x], it[3][y])) for it in items for x in it[1] for y in it[2]
              if x in it[3] and y in it[3] and it[3][x] != it[3][y]}
    grid_spread(P, {e: label.get(e, (120, 40)) for e in P}, linked)
    # compartment boxes the layout would cover move to a column right of it
    lay = (min(q[0] for q in P.values()) - 250, area[1], max(q[0] for q in P.values()) + 250,
           max(q[1] for q in P.values()) + 250)
    def rect_of(g):
        r = g.find(L.SVG + "rect")
        x0, y0 = float(r.get("x")), float(r.get("y"))
        return (x0, y0, x0 + float(r.get("width")), y0 + float(r.get("height")))

    def overlaps(a, b, gap=100):
        return a[0] < b[2] + gap and a[2] > b[0] - gap and a[1] < b[3] + gap and a[3] > b[1] - gap
    groups = [g for g in mp.root.iter(L.SVG + "g") if g.get("class") == "compartment" and g.find(L.SVG + "rect") is not None]
    for g in groups:
        b = rect_of(g)
        if not overlaps(b, lay, 0):
            continue
        w, h = b[2] - b[0], b[3] - b[1]
        others = [rect_of(o) for o in groups if o is not g]
        x, y = lay[2] + 200, area[1] - 100
        while any(overlaps((x, y, x + w, y + h), o) for o in others):  # first free place in the column, then
            y += 50                                                     # the next column
            if y > area[1] + 6000:
                x, y = x + 1200, area[1] - 100
        if x + w + margin > L.canvas(mp)[2]:
            L.widen_canvas(dr, x + w + margin - L.canvas(mp)[2])
        if y + h + margin > L.canvas(mp)[3]:
            L.grow_canvas(dr, y + h + margin - L.canvas(mp)[3])
        L.move_box(mp, g, x - b[0], y - b[1])
        log.append(("D8", "compartment box moved beside the KEGG layout", g.get("id"), "", ""))
    dr.space = L.Space(mp)
    need = max(q[1] for q in P.values()) + 300 + margin - L.canvas(mp)[3]
    if need > 0:
        L.grow_canvas(dr, need)
    # reactions between the same two KEGG compounds get parallel lanes
    pair = {}
    for it in items:
        es = [it[3][x] for x in it[1] if x in it[3]][:1] + [it[3][x] for x in it[2] if x in it[3]][:1]
        pair.setdefault(frozenset(es), []).append(it[0])
    nodes = {}
    for rid, subs, prods, amap, line in items:
        es = [amap[x] for x in subs if x in amap][:1] + [amap[x] for x in prods if x in amap][:1]
        lane = pair[frozenset(es)]
        off = (lane.index(rid) - (len(lane) - 1) / 2) * lane_gap
        pos, existing = {}, {}
        for met, eid in amap.items():
            pos[met] = P[eid]
            near = [k for k in nodes if k[1] == met and L.length(L.sub(P[k[0]], pos[met])) < 600]
            if (eid, met) in nodes:
                existing[met] = nodes[(eid, met)]
            elif near:  # KEGG shows this compound twice close by: one node
                existing[met] = nodes[near[0]]
                pos[met] = L.translate(nodes[near[0]])
                amap = {**amap, met: near[0][0]}
            else:
                # another metabolite already on this KEGG entry (an anomer): set this one beside it
                taken = [k for k in nodes if k[0] == eid and k[1] != met]
                if taken:
                    pos[met] = L.add(pos[met], (0.0, 1.0), 110 * len(taken))
        a = [pos[x] for x in subs if x in pos]
        b = [pos[x] for x in prods if x in pos]
        rpos = mid_line([place(q) for q in line] if line else None, a[0] if a else None, b[0] if b else None)
        rpos = rpos or (a[0] if a else b[0])
        if any(k[0] in amap.values() and k[1] not in amap for k in nodes) and a and b:
            rpos = L.add(a[0], L.sub(b[0], a[0]), 0.5)
        if a and b and (L.length(L.sub(rpos, a[0])) < 70 or L.length(L.sub(rpos, b[0])) < 70):
            rpos = L.add(a[0], L.sub(b[0], a[0]), 0.5)
        tangent = L.unit(L.sub(b[0], a[0])) if a and b and L.length(L.sub(b[0], a[0])) > 1 else (1.0, 0.0)
        for k, x in enumerate([x for x in subs + prods if x not in pos]):
            side = a[0] if x in subs and a else b[0] if b else rpos
            pos[x] = L.add(side, L.normal(tangent), 70 * (k + 1))
        if a and b and not (min(a[0][0], b[0][0]) - 40 <= rpos[0] <= max(a[0][0], b[0][0]) + 40 and
                            min(a[0][1], b[0][1]) - 40 <= rpos[1] <= max(a[0][1], b[0][1]) + 40):
            rpos = L.add(a[0], L.sub(b[0], a[0]), 0.5)  # the enzyme box no longer between them (spread)
        if a and b:  # substrate and product on one row or column: the reaction on that line
            d = L.sub(b[0], a[0])
            if abs(d[1]) < 0.15 * abs(d[0]) and min(a[0][0], b[0][0]) < rpos[0] < max(a[0][0], b[0][0]):
                rpos = (rpos[0], a[0][1] + d[1] * (rpos[0] - a[0][0]) / d[0])
            elif abs(d[0]) < 0.15 * abs(d[1]) and min(a[0][1], b[0][1]) < rpos[1] < max(a[0][1], b[0][1]):
                rpos = (a[0][0] + d[0] * (rpos[1] - a[0][1]) / d[1], rpos[1])
        if off:
            rpos = L.add(rpos, L.normal(tangent), off)
        L.build(dr, rid, subs, prods, pos, rpos, tangent, existing, orthogonal(pos, rpos, a, b))
        for met, eid in amap.items():
            node = next((n for n in dr.met_nodes(met) if rid in mp.classes(n) and dr.is_main(n)), None)
            if node is not None:
                nodes[(eid, met)] = node
        log.append(("D8", "reaction drawn as in KEGG" + (" (was missing)" if rid in missing else ""), rid, kg.name,
                    model.rxns[rid].get("name", "") or ""))
    # the compartment's other reactions: D3 places them next to the new layout or in the bottom section
    comps = {L.COMPARTMENT[comp]}
    for rid in rest:
        if rid in drawn and not (scope and scope.contains(rid)):
            L.draw_missing(dr, rid, comps)
    log.append(("D8", "KEGG relayout", kg.name, f"{len(items)} of {len(drawn) + len(missing)} reactions ({share:.0%}), scale {scale[0]:.1f}",
                f"{len([r for r in rest if r in drawn])} drawn reactions not on the KEGG map re-added by D3"))
    return True
