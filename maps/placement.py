"""D9: place added reactions as blocks.

A block is a group of added reactions drawn together: a module of a KEGG map (D7b), a chain of reactions
(the product of one is the substrate of the next), or a single reaction. Each block has its own layout
(local coordinates). It is placed, in this order of preference,
  1. docked: one of its metabolites is put on a node of that metabolite already on the map, in any of
     the eight orientations (turned or mirrored), where the rest of the block is clear of the drawing;
  2. in free space of its compartment, as close as possible to the map's nodes of its metabolites;
  3. below the content (the canvas grows); for a boxed compartment, in its box made larger (D11), and for a
     compartment the map has no box for, in a new box sized to its blocks.
KEGG gives a block its first shape; it is a guide, not a rule: blocks are turned, mirrored and
spread apart as the map needs.
"""

import collections
import math

import layout as L

ORIENT = [(1, 0, 0, 1), (-1, 0, 0, 1), (1, 0, 0, -1), (-1, 0, 0, -1),
          (0, 1, 1, 0), (0, -1, 1, 0), (0, 1, -1, 0), (0, -1, -1, 0)]
NODE_HALF = (85, 40)
MIN_GAP = 250      # px between two compound nodes of a block (hubs)
LANE_GAP = 200     # px between reactions that join the same two compounds
MID_SCALE = 1.6    # "mid" hubs: a crowded module grows at most this much before its close compounds spread
STEP = L.ROW_STEP  # metabolite to reaction in chains
BLOCK_GAP = 30     # px kept free around a placed block


class Block:
    def __init__(self, kind, label):
        self.kind, self.label = kind, label
        self.nodes = {}   # key -> (metabolite, local (x, y))
        self.steps = []   # (rid, subs, prods, {met: key}, local rpos, {met: local elbow})

    def mets(self):
        return {m for m, _ in self.nodes.values()}

    def bbox(self, o=ORIENT[0]):
        pts = [tr(o, p) for _, p in self.nodes.values()] + [tr(o, s[4]) for s in self.steps]
        return (min(p[0] for p in pts), min(p[1] for p in pts), max(p[0] for p in pts), max(p[1] for p in pts))


def tr(o, p, t=(0.0, 0.0)):
    return (o[0] * p[0] + o[1] * p[1] + t[0], o[2] * p[0] + o[3] * p[1] + t[1])


def mains(model, rid):
    st = model.rxns[rid]["stoich"]
    subs = [m for m in st if st[m] < 0 and model.name(m) not in L.COFACTORS] or [m for m in st if st[m] < 0][:1]
    prods = [m for m in st if st[m] > 0 and model.name(m) not in L.COFACTORS] or [m for m in st if st[m] > 0][:1]
    return subs, prods


# ----------------------------------------------------------------------------- building blocks

def kegg_block(dr, kg, items, hubs="mid"):
    """A block from a KEGG module: compounds where KEGG has them, scaled; hubs opened up."""
    import kegg_layout as K
    b = Block("kegg", kg.title or kg.name)
    sc = K.module_scale(kg, items)
    pts = [kg.compounds[e][1] for it in items for e in it[3].values()]
    x0, y0 = min(p[0] for p in pts), min(p[1] for p in pts)
    ents = sorted({e for it in items for e in it[3].values()})
    if hubs in ("scale", "mid"):  # the whole module larger until no two compounds are closer than MIN_GAP
        top = 8 if hubs == "scale" else sc * MID_SCALE  # "mid": at most MID_SCALE x larger, then spread
        for _ in range(30):
            close = [1 for i, a in enumerate(ents) for c in ents[i + 1:]
                     if sc * L.length(L.sub(kg.compounds[a][1], kg.compounds[c][1])) < MIN_GAP]
            if not close or sc >= top:
                break
            sc *= 1.1
    pos = {e: (sc * (kg.compounds[e][1][0] - x0), sc * (kg.compounds[e][1][1] - y0)) for e in ents}
    if hubs in ("spread", "mid"):  # only the compounds that are (still) too close move apart
        spread(pos)
    for rid, subs, prods, amap, line in items:
        keys = {}
        for m in subs + prods:
            if m in amap:
                keys[m] = (amap[m], m)
                b.nodes.setdefault(keys[m], (m, pos[amap[m]]))
        b.steps.append([rid, subs, prods, keys, None, {}])
    finish(b)
    return b


def spread(pos, gap=MIN_GAP, rounds=60, size=None, linked=frozenset(), between=180):
    """Move nodes that are too close apart along the line joining them. With size {key: (w, h)} (label
    sizes), the distance needed depends on the labels along that direction, at least gap; pairs in linked
    (joined by a reaction) need `between` more."""
    keys = list(pos)
    for _ in range(rounds):
        moved = False
        for i, a in enumerate(keys):
            for c in keys[i + 1:]:
                d = L.sub(pos[c], pos[a])
                dist = L.length(d)
                u = L.unit(d) if dist > 1 else (math.cos(i), math.sin(i))
                need = gap
                if size:
                    (wa, ha), (wc, hc) = size[a], size[c]
                    need = max(gap * 0.7, abs(u[0]) * (wa + wc) / 2 + abs(u[1]) * (ha + hc) / 2 + 40)
                if frozenset((a, c)) in linked:  # a reaction (diamond, cofactors, genes) sits between them
                    need += between
                if dist < need:
                    push = (need - dist) / 2
                    pos[a] = L.add(pos[a], u, -push)
                    pos[c] = L.add(pos[c], u, push)
                    moved = True
        if not moved:
            break


def chain_block(model, chain, per_row=6):
    """A block from a chain of reactions: left to right, wrapping as a serpentine after per_row reactions."""
    b = Block("chain", "")
    pos, heading = (0.0, 0.0), (1.0, 0.0)
    prev = {}
    for i, rid in enumerate(chain):
        if i and i % per_row == 0:
            heading = (0.0, 1.0)  # one step down, then back the other way
        elif i and i % per_row == 1 and i > 1:
            heading = (-1.0 if (i // per_row) % 2 else 1.0, 0.0)
        subs, prods = mains(model, rid)
        side = L.normal(heading)
        keys = {}
        # new substrates go beside the previous reaction's products, not on them
        taken = [b.nodes[key][1] for key in prev.values()]
        k = 0
        for m in subs:
            if m in prev:
                keys[m] = prev[m]
            else:
                while any(L.length(L.sub(L.add(pos, side, 90.0 * k), q)) < 45 for q in taken):
                    k += 1
                keys[m] = (rid, m)
                b.nodes[keys[m]] = (m, L.add(pos, side, 90.0 * k))
                taken.append(b.nodes[keys[m]][1])
                k += 1
        end = L.add(pos, heading, 2 * STEP)
        for k2, m in enumerate(prods):
            keys[m] = (rid, m)
            b.nodes[keys[m]] = (m, L.add(end, side, 90.0 * k2))
        b.steps.append([rid, subs, prods, keys, None, {}])
        prev = {m: keys[m] for m in prods}
        pos = end
    finish(b)
    return b


def finish(b):
    """Main metabolites without a place get one beside their reaction; reaction positions, with lanes."""
    for s in b.steps:
        rid, subs, prods, keys = s[:4]
        a = [b.nodes[keys[m]][1] for m in subs if m in keys]
        c = [b.nodes[keys[m]][1] for m in prods if m in keys]
        base = a[0] if a else c[0] if c else (0.0, 0.0)
        tangent = L.unit(L.sub(c[0], a[0])) if a and c and L.length(L.sub(c[0], a[0])) > 1 else (1.0, 0.0)
        for k, m in enumerate([m for m in subs + prods if m not in keys]):
            keys[m] = (rid, m)
            side = a[0] if m in subs and a else c[0] if c else base
            b.nodes[keys[m]] = (m, L.add(L.add(side, L.normal(tangent), 90 * (k + 1)), tangent, -60 if m in subs else 60))
    pair = collections.defaultdict(list)
    for s in b.steps:
        subs, prods, keys = s[1], s[2], s[3]
        pair[frozenset([keys[subs[0]], keys[prods[0]]] if subs and prods else [s[0]])].append(s)
    for lane in pair.values():
        for i, s in enumerate(lane):
            rid, subs, prods, keys = s[:4]
            a = b.nodes[keys[subs[0]]][1] if subs else None
            c = b.nodes[keys[prods[0]]][1] if prods else None
            if a and c:
                mid = L.add(a, L.sub(c, a), 0.5)
                tangent = L.unit(L.sub(c, a)) if L.length(L.sub(c, a)) > 1 else (1.0, 0.0)
            else:
                mid, tangent = L.add(a or c, (1.0, 0.0), STEP), (1.0, 0.0)
            off = (i - (len(lane) - 1) / 2) * LANE_GAP
            s[4] = L.add(mid, L.normal(tangent), off)
            # KEGG-style edges with one bend (lanes come out parallel)
            import kegg_layout as K
            s[5] = K.orthogonal({m: b.nodes[keys[m]][1] for m in subs + prods}, s[4], [a] if a else [], [c] if c else [])


# ----------------------------------------------------------------------------- footprint and fit

class Occupancy:
    """What is drawn, as a grid of 40 px cells with a summed-area table: a rectangle is free in one lookup."""
    CELL = 40

    def __init__(self, dr):
        c = self.CELL
        w, h = L.canvas(dr.mp)[2:]
        self.nx, self.ny = int(w // c) + 2, int(h // c) + 2
        self.occ = bytearray(self.nx * self.ny)
        seen = set()
        for boxes in dr.space.grid.values():
            for b in boxes:
                if id(b) not in seen:
                    seen.add(id(b))
                    self.mark(b)
        dr.space.recent = []
        self.table()

    def mark(self, b):
        c = self.CELL
        i0, i1 = max(0, int(b[0] // c)), min(self.nx - 1, int(b[2] // c))
        j0, j1 = max(0, int(b[1] // c)), min(self.ny - 1, int(b[3] // c))
        for j in range(j0, j1 + 1):
            row = j * self.nx
            self.occ[row + i0:row + i1 + 1] = b"\x01" * (i1 - i0 + 1) if i1 >= i0 else b""

    def table(self):
        n1, occ = self.nx + 1, self.occ
        sat = [0] * (n1 * (self.ny + 1))
        for j in range(self.ny):
            run = 0
            row, up, cur = j * self.nx, j * n1, (j + 1) * n1
            for i in range(self.nx):
                run += occ[row + i]
                sat[cur + i + 1] = sat[up + i + 1] + run
        self.sat, self.n1 = sat, n1

    def update(self, dr):
        """Take in what was drawn since the last update."""
        if dr.space.recent:
            for b in dr.space.recent:
                self.mark(b)
            dr.space.recent = []
            self.table()

    def free(self, box, margin):
        c = self.CELL
        i0, i1 = max(0, int((box[0] - margin) // c)), min(self.nx - 1, int((box[2] + margin) // c))
        j0, j1 = max(0, int((box[1] - margin) // c)), min(self.ny - 1, int((box[3] + margin) // c))
        if i0 > i1 or j0 > j1:
            return True
        s, n1 = self.sat, self.n1
        return s[(j1 + 1) * n1 + i1 + 1] - s[j0 * n1 + i1 + 1] - s[(j1 + 1) * n1 + i0] + s[j0 * n1 + i0] == 0


def footprint(dr, b, o, t, skip=()):
    boxes = []
    for key, (m, p) in b.nodes.items():
        if key in skip:
            continue
        q = tr(o, p, t)
        boxes.append((q[0] - NODE_HALF[0], q[1] - NODE_HALF[1], q[0] + NODE_HALF[0], q[1] + NODE_HALF[1]))
    for s in b.steps:
        q = tr(o, s[4], t)
        n = len(dr.model.rxns[s[0]]["genes"])
        h = 95 + 16 * math.ceil(n / 3)
        boxes.append((q[0] - h, q[1] - h, q[0] + h, q[1] + h))
    return boxes


def fits(dr, boxes, reg, holes):
    edge = 120 if holes or reg[0] == 0.0 else 60  # 60 inside a compartment box
    for bx in boxes:
        if bx[0] < reg[0] + edge or bx[2] > reg[2] - edge or bx[1] < reg[1] + edge or bx[3] > reg[3] - edge:
            return False
        if any(bx[0] < h[2] + 30 and bx[2] > h[0] - 30 and bx[1] < h[3] + 30 and bx[3] > h[1] - 30 for h in holes):
            return False
        occ = getattr(dr, "occ", None)
        if not (occ.free(bx, BLOCK_GAP) if occ else dr.space.free(bx, BLOCK_GAP)):
            return False
    return True


def in_region(p, reg, holes):
    return p and reg[0] <= p[0] <= reg[2] and reg[1] <= p[1] <= reg[3] and not any(
        h[0] <= p[0] <= h[2] and h[1] <= p[1] <= h[3] for h in holes)


def node_index(dr):
    """Metabolite -> [(position, node)] of its main nodes."""
    index = collections.defaultdict(list)
    for n in dr.mp.nodes("met"):
        if dr.is_main(n) and L.translate(n):
            index[dr.mp.classes(n)[1]].append((L.translate(n), n))
    return index


def clash(b, o):
    """Pairs of the block's own nodes whose labels overlap in orientation o: metabolites of one reaction sit 90 px
    apart, one above the other, and a block turned by 90 degrees puts them side by side."""
    pts = [tr(o, p) for _, p in b.nodes.values()]
    return sum(1 for i, p in enumerate(pts) for q in pts[i + 1:]
               if abs(p[0] - q[0]) < 2 * NODE_HALF[0] and abs(p[1] - q[1]) < 2 * NODE_HALF[1])


def dock(dr, b, reg, holes, index):
    """Best docking of b onto drawn nodes: (orientation, translation, {key: node}) or None. The most nodes landing on
    their own drawn node first; among those, the orientation in which the fewest of the block's labels overlap."""
    best = None
    clashes = [clash(b, o) for o in ORIENT]
    for key, (m, p) in b.nodes.items():
        for q, n in index.get(m, ()):
            if not in_region(q, reg, holes):
                continue
            for oi, o in enumerate(ORIENT):
                t = L.sub(q, tr(o, p))
                anchors = {key: n}
                for k2, (m2, p2) in b.nodes.items():  # other block metabolites landing on their own node
                    if k2 != key:
                        q2 = tr(o, p2, t)
                        n2 = next((x for qx, x in index.get(m2, ()) if L.length(L.sub(qx, q2)) < 60), None)
                        if n2 is not None:
                            anchors[k2] = n2
                cost = (-len(anchors), clashes[oi], oi)
                if best is not None and cost >= best[0]:
                    continue
                if fits(dr, footprint(dr, b, o, t, set(anchors)), reg, holes):
                    best = (cost, o, t, anchors)
    return best[1:] if best else None


def free_spot(dr, b, reg, holes, near, step=None):
    """(orientation, translation) putting b in free space of reg, closest to `near` (or top-left first)."""
    # about 4000 positions per orientation, whatever the size of the region
    step = step or max(100.0, math.sqrt((reg[2] - reg[0]) * (reg[3] - reg[1]) / 4000))
    cands = []
    for o in (ORIENT[0], ORIENT[4]):
        bb = b.bbox(o)
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        penalty = (0 if o is ORIENT[0] else 50) + 300 * min(clash(b, o), 3)  # px: overlapping labels count as further away
        y = reg[1] + 200
        while y + h < reg[3] - 200:
            x = reg[0] + 200
            while x + w < reg[2] - 200:
                c = (x + w / 2, y + h / 2)
                d = math.hypot(c[0] - near[0], c[1] - near[1]) if near else y * 10 + x
                cands.append((d + penalty, o, (x - bb[0], y - bb[1])))
                x += step
            y += step
    cands.sort(key=lambda c: c[0])
    for _, o, t in cands:
        if fits(dr, footprint(dr, b, o, t), reg, holes):
            return o, t
    return None


def draw(dr, b, o, t, anchors, how):
    mp, model = dr.mp, dr.model
    created = dict(anchors)
    for rid, subs, prods, keys, rloc, eloc in b.steps:
        pos, existing = {}, {}
        for m in subs + prods:
            k = keys[m]
            if k in created:
                existing[m] = created[k]
                pos[m] = L.translate(created[k])
            else:
                pos[m] = tr(o, b.nodes[k][1], t)
        rpos = tr(o, rloc, t)
        a = [pos[m] for m in subs]
        c = [pos[m] for m in prods]
        tangent = L.unit(L.sub(c[0], a[0])) if a and c and L.length(L.sub(c[0], a[0])) > 1 else (1.0, 0.0)
        elbow = {m: tr(o, e, t) for m, e in eloc.items()}
        L.build(dr, rid, subs, prods, pos, rpos, tangent, existing, elbow)
        tail = mp.layer["nodes"][-(len(subs) + len(prods) + 40):]  # new nodes are appended
        for m in subs + prods:
            n = existing[m] if m in existing else next((x for x in reversed(tail) if x.get("class", "").startswith("met ") and
                                         mp.classes(x)[1] == m and rid in mp.classes(x) and dr.is_main(x)), None)
            if n is not None:
                created.setdefault(keys[m], n)
        ctx = rid in dr.context
        dr.log.append(("D14" if ctx else "D3", "context reaction drawn" if ctx else "missing reaction added", rid, how,
                       model.rxns[rid].get("name", "") or ""))
    return created


# ----------------------------------------------------------------------------- placing

def place_pending(dr, hubs="mid"):
    """D9: place the collected reactions (dr.pending) block by block."""
    mp, model, log = dr.mp, dr.model, dr.log
    by_comp = collections.defaultdict(list)
    for rid, comp in dr.pending:
        if rid not in mp.reactions() and rid not in by_comp[comp]:
            by_comp[comp].append(rid)
    dr.pending = []
    on_map = L.compartments_on_map(mp) | dr.comps_on_map
    for comp, rids in by_comp.items():
        blocks = []
        if dr.kegg:
            import kegg_layout as K
            kblocks, rids = K.bottom_blocks(dr, rids)
            blocks += [kegg_block(dr, kg, items, hubs) for kg, items in kblocks]
        blocks += [chain_block(model, ch) for ch in L.chains(model, sorted(rids))]
        blocks.sort(key=lambda b: -len(b.steps))
        boxed = comp in L.compartment_boxes(mp) or L.ALIASES.get(comp) in L.compartment_boxes(mp)
        leftover = []
        dr.occ = Occupancy(dr)
        index = node_index(dr)

        def drawn(created):
            for n in dict.fromkeys(created.values()):  # drawing order, not memory order
                m, q = mp.classes(n)[1], L.translate(n)
                if q and not any(x is n for _, x in index[m]):
                    index[m].append((q, n))
        for b in blocks:
            if comp not in on_map:
                leftover.append(b)
                continue
            dr.occ.update(dr)
            reg, holes = L.region(mp, comp)
            d = dock(dr, b, reg, holes, index)
            if d:
                o, t, anchors = d
                names = ", ".join(sorted({model.name(b.nodes[k][0]) for k in anchors}))
                drawn(draw(dr, b, o, t, anchors, f"docked onto {names}" + (f" ({b.label})" if b.kind == "kegg" else "")))
                continue
            near = [q for m in b.mets() for q, _ in index.get(m, ())]
            near = [p for p in near if in_region(p, reg, holes)]
            c = (sum(p[0] for p in near) / len(near), sum(p[1] for p in near) / len(near)) if near else \
                L.content_center(mp, reg, holes)
            f = free_spot(dr, b, reg, holes, c)
            if f:
                drawn(draw(dr, b, f[0], f[1], {}, "in free space" + (" near its metabolites" if near else "") +
                           (f" ({b.label})" if b.kind == "kegg" else "")))
                continue
            if boxed:
                leftover.append(b)
                continue
            created = below_drawing(dr, b, comp, c, "below the drawing" + (f" ({b.label})" if b.kind == "kegg" else ""))
            if created is not None:
                drawn(created)
            else:
                log.append(("review", "added reaction could not be placed", b.steps[0][0], comp, ""))
        dr.occ = None
        if leftover:
            import compact
            if not (boxed and compact.grow_box(dr, comp, leftover, draw, ORIENT[0])):
                new_box(dr, comp, leftover)


def below_drawing(dr, b, comp, near, how):
    """Step 3 for a block of an unboxed compartment: the canvas grows (wider, or taller while the map stays
    landscape) and the block goes into the new space, as close to `near` as there is room. The created nodes,
    or None if it still finds no place."""
    mp = dr.mp
    reg, holes = L.region(mp, comp)
    bb = b.bbox()
    cw, ch = L.canvas(mp)[2:]
    wide = (bb[2] - bb[0]) + 800 - (reg[2] - reg[0])
    if wide > 0 or ch > 0.7 * cw:  # keep the map from becoming a tall strip
        L.widen_canvas(dr, max(wide, bb[2] - bb[0] + 800, 0.25 * cw))
    else:
        L.grow_canvas(dr, bb[3] - bb[1] + 500)
    reg, holes = L.region(mp, comp)  # the block must fit the other way too
    if (bb[3] - bb[1]) + 500 > reg[3] - reg[1]:
        L.grow_canvas(dr, (bb[3] - bb[1]) + 500 - (reg[3] - reg[1]))
    if (bb[2] - bb[0]) + 500 > reg[2] - reg[0]:
        L.widen_canvas(dr, (bb[2] - bb[0]) + 500 - (reg[2] - reg[0]))
    dr.occ = Occupancy(dr)
    reg, holes = L.region(mp, comp)
    f = free_spot(dr, b, reg, holes, near)
    return draw(dr, b, f[0], f[1], {}, how) if f else None


PAD = 180  # px around a block in a new box (cofactors, genes)


def new_box(dr, comp, blocks, gap=200):
    """A new box of compartment comp below the drawing, sized to its blocks, laid out in rows."""
    import compact
    mp = dr.mp
    width = min(float(mp.root.get("width")) - 800,
                max(max(b.bbox()[2] - b.bbox()[0] + 2 * PAD for b in blocks), compact.squarish(blocks, PAD)))
    rows, row, x = [], [], 0.0
    for b in blocks:
        bb = b.bbox()
        w = bb[2] - bb[0] + 2 * PAD
        if row and x + w > width:
            rows.append(row)
            row, x = [], 0.0
        row.append((b, x, bb))
        x += w + gap
    if row:
        rows.append(row)
    heights = [max(bb[3] - bb[1] for _, _, bb in r) + 2 * PAD for r in rows]
    inner_w = max(sum(bb[2] - bb[0] + 2 * PAD + gap for _, _, bb in r) - gap for r in rows)
    height = sum(heights) + gap * (len(rows) - 1) - 100
    at = compact.free_box_spot(dr, inner_w + 300, height + 150 + 160)  # D13: in free space next to the boxes
    reg = L.new_compartment_box(dr, comp, height, inner_w + 300, at=at)
    y = reg[1]
    for r, h in zip(rows, heights):
        for b, x, bb in r:
            t = (reg[0] + 150 + x + PAD - bb[0], y + PAD - bb[1])
            draw(dr, b, ORIENT[0], t, {}, f"in a new {comp} box" + (f" ({b.label})" if b.kind == "kegg" else ""))
        y += h + gap
