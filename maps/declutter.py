"""D15: cofactor dots clear of other elements. A cofactor dot whose dot or label overlaps another node, label or
gene box moves to a free place around its reaction, with its label on the side away from the reaction; its edge
and arrowhead move with it. A place counts as free when the dot and label overlap nothing and no edge crosses the
label, the dot's edge passes no node and crosses at most one more edge than before, and it keeps clear of the
reaction's main edges. The nearest such
place is taken; a dot with no free place stays. Only dots joined to one reaction by a straight edge move, and
only those of reactions placed in this run unless every dot may move (scope "all").
"""

import collections
import math

import layout as L
from mapedit import near

RADII = (62, 80, 100, 125, 150, 180, 220)
ANGLES = 24          # directions tried around the reaction
SHRINK, MIN_AREA = 3, 150  # an overlap smaller than this is a touch, as in the overlap measure
MAIN_CLEAR = math.radians(25)  # between a dot's edge and a main edge of its reaction


def overlap(a, b):
    w = min(a[2], b[2]) - max(a[0], b[0]) - 2 * SHRINK
    h = min(a[3], b[3]) - max(a[1], b[1]) - 2 * SHRINK
    return w > 0 and h > 0 and w * h >= MIN_AREA


class Grid:
    def __init__(self, cell=200):
        self.cell, self.items = cell, collections.defaultdict(list)

    def keys(self, b):
        c = self.cell
        return [(i, j) for i in range(int(b[0] // c), int(b[2] // c) + 1) for j in range(int(b[1] // c), int(b[3] // c) + 1)]

    def add(self, b, owner):
        for k in self.keys(b):
            self.items[k].append((b, owner))

    def remove(self, owner):
        for k in list(self.items):
            self.items[k] = [x for x in self.items[k] if x[1] is not owner]

    def near(self, b):
        seen = set()
        for k in self.keys(b):
            for x in self.items.get(k, ()):
                if id(x) not in seen:
                    seen.add(id(x))
                    yield x


def side_box(p, name, left):
    w = max(L.label_width(name, 18), 60)
    return (p[0] - (w + 12 if left else 8), p[1] - 12, p[0] + (8 if left else w + 12), p[1] + 12)


def cross(a, b, c, d):
    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    return o1 * o2 < 0 and o3 * o4 < 0


def declutter(dr, scope="added"):
    """Move overlapping cofactor dots (D15); returns the number moved."""
    mp, log = dr.mp, dr.log
    rpos = {rid: L.translate(g) for rid, g in mp.reactions().items()}
    boxes = Grid()
    owner_rxns = {}
    for n in mp.layer["nodes"]:
        cls = mp.classes(n)
        b = L.node_box(mp, n) if cls else None
        if not b:
            continue
        boxes.add(b, n)
        owner_rxns[id(n)] = {n.get("id")} if cls[0] == "rea" else {c for c in cls[2:] if L.is_rxn(c)}
    for g in mp.root.iter(L.SVG + "g"):
        if g.get("class") in ("compartment", "subsystem"):
            for t in g.iter(L.SVG + "text"):
                b = L.text_box(t)
                if b:
                    boxes.add(b, t)
                    owner_rxns[id(t)] = set()
    segs = Grid(100)
    gone = set()  # edges of dots that moved
    for path in mp.layer["fes"]:
        rids = {c for c in mp.classes(path)[1:] if L.is_rxn(c)}
        for sp in L.subpaths(path.get("d") or "") or []:
            for a, b in zip(sp, sp[1:]):
                segs.add((min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1])), (a, b, frozenset(rids)))

    def clashes(b, n, rid):
        return [o for ob, o in boxes.near(b) if o is not n and overlap(b, ob) and
                not (rid in owner_rxns.get(id(o), ()) and mp.classes(o)[:1] == ["rea"])]

    def label_crossed(nb, q, own):
        """An edge (other than the dot's own) runs through the label (the box beside the dot)."""
        lab = (nb[0] + (0 if nb[0] < q[0] - 9 else 10), nb[1] + 2, nb[2] - (0 if nb[2] > q[0] + 9 else 10), nb[3] - 2)
        lab = (min(lab[0], lab[2]), lab[1], max(lab[0], lab[2]), lab[3])
        corners = [(lab[0], lab[1]), (lab[2], lab[1]), (lab[2], lab[3]), (lab[0], lab[3])]
        for _, (c, d, rids) in segs.near(lab):
            if (c, d) == own or (d, c) == own or (c, d) in gone:
                continue
            inside = any(lab[0] <= x[0] <= lab[2] and lab[1] <= x[1] <= lab[3] for x in (c, d))
            if inside or any(cross(c, d, corners[i], corners[(i + 1) % 4]) for i in range(4)):
                return True
        return False

    def crossings(a, b, rid):
        box = (min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1]))
        return sum(1 for _, (c, d, rids) in segs.near(box) if rid not in rids and (c, d) not in gone and cross(a, b, c, d))

    def through_node(a, b, n, rid):
        k = max(1, int(L.length(L.sub(b, a)) / 8))
        for i in range(1, k):
            q = L.add(a, L.sub(b, a), i / k)
            for ob, o in boxes.near((q[0] - 1, q[1] - 1, q[0] + 1, q[1] + 1)):
                if o is n or (rid in owner_rxns.get(id(o), ()) and mp.classes(o)[:1] == ["rea"]):
                    continue
                if ob[0] + 2 <= q[0] <= ob[2] - 2 and ob[1] + 2 <= q[1] <= ob[3] - 2:
                    return True
        return False

    moved = 0
    sides = [n for n in mp.nodes("met") if not dr.is_main(n)]
    for n in sides:
        cls = mp.classes(n)
        rids = [c for c in cls[2:] if L.is_rxn(c)]
        if len(rids) != 1 or rids[0] not in rpos or not rpos[rids[0]]:
            continue
        rid = rids[0]
        if scope != "all" and rid not in dr.placed:
            continue
        p = L.translate(n)
        b = L.node_box(mp, n)
        if not p or not b or not clashes(b, n, rid):
            continue
        fe = mp.edges("fes", rid)
        sps = L.subpaths(fe[0].get("d") or "") if fe else None
        if not sps:
            continue
        mine = [k for k, sp in enumerate(sps) if len(sp) == 2 and near(sp[0], p, 16)]
        if len(mine) != 1:
            continue
        r = rpos[rid]
        main_dirs = [math.atan2(sp[-1][1] - r[1], sp[-1][0] - r[0]) if near(sp[-1], r, 40) else
                     math.atan2(sp[0][1] - r[1], sp[0][0] - r[0]) for k, sp in enumerate(sps) if k not in mine and
                     any(near(q, r, 40) for q in (sp[0], sp[-1])) and len(sp) >= 2]
        name = "".join(n.find(L.SVG + "text").itertext()).strip() if n.find(L.SVG + "text") is not None else ""
        old_edge = sps[mine[0]]
        before = crossings(old_edge[0], old_edge[1], rid)
        best = None
        for rad in RADII:
            for k in range(ANGLES):
                ang = 2 * math.pi * k / ANGLES
                if any(abs((ang - d + math.pi) % (2 * math.pi) - math.pi) < MAIN_CLEAR for d in main_dirs):
                    continue
                q = (r[0] + rad * math.cos(ang), r[1] + rad * math.sin(ang))
                left = q[0] < r[0] - 1
                nb = side_box(q, name, left)
                if clashes(nb, n, rid):
                    continue
                a, e = L.boundary(q, r, "side"), L.boundary(r, q, "rea")
                if label_crossed(nb, q, (tuple(old_edge[0]), tuple(old_edge[1]))):
                    continue
                if through_node(a, e, n, rid) or crossings(a, e, rid) > before + 1:
                    continue
                d = L.length(L.sub(q, p))
                if best is None or d < best[0]:
                    best = (d, q, left, nb, a, e)
            if best:
                break
        if not best:
            continue
        _, q, left, nb, a, e = best
        L.set_position(n, q[0], q[1])
        t = n.find(L.SVG + "text")
        if t is not None:
            t.set("text-anchor", "end" if left else "start")
            t.set("x", "-12" if left else "12")
        gone.add((tuple(old_edge[0]), tuple(old_edge[1])))
        sps[mine[0]] = [a, e]
        fe[0].set("d", L.to_d(sps))
        for h in mp.edges("fehs", rid):
            hs = L.subpaths(h.get("d") or "")
            if hs:
                hs = [L.arrow(a, e) if near(s[0], old_edge[0], 4) else s for s in hs]
                h.set("d", L.to_d(hs))
        boxes.remove(n)
        boxes.add(nb, n)
        segs.add((min(a[0], e[0]), min(a[1], e[1]), max(a[0], e[0]), max(a[1], e[1])), (a, e, frozenset([rid])))
        moved += 1
    if moved:
        log.append(("D15", "cofactor dots moved clear of other elements", "", f"{moved} dot(s)", ""))
    return moved
