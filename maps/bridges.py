"""D14: bridges. On a subsystem map, a reaction of another subsystem that joins two parts of the map that are
not connected, in one step, is added as a context reaction (drawn without the subsystem's band). With
steps=2, when no single reaction joins two parts, a pair of reactions through one metabolite not yet on the
map may (one-step bridges are always tried first).

Parts are the drawn reactions (and those about to be placed) joined through shared metabolites. A bridge
joins two parts through metabolites already on the map; cofactors never join parts, and metabolites in more
than `limit` reactions of the model (hubs) only when `limit` is None ("all", the default). Model artefacts (R11),
exchange and transport reactions are not used. The bridge joining the largest parts is added first, until
no reaction joins two parts. The subsystems in SKIP get no bridges: Isolated (isolated by definition) and
tRNA charging (Aminoacyl-tRNA biosynthesis, Yeast-GEM's tRNA metabolism), where amino acids would join
everything. After placement, D10 links a bridge's metabolites over longer lines than usual; each bridge is then
drawn between the nodes it joins (draw_bridges), moving the smaller or the larger of its two parts (at most
MOVE_MAX reactions) as a whole next to the other when they are far apart or there is no room between them,
and is removed if, as drawn, it does not join two parts (check_bridges).
"""

import collections

import layout as L
from mapedit import comp_of

# subsystems without bridges: Isolated (isolated by definition) and tRNA charging, whose amino acids are hubs
# that would join every reaction to every other (round 10 answer)
SKIP = {"Isolated", "Aminoacyl-tRNA biosynthesis", "tRNA metabolism"}


def degree(model):
    if id(model) not in degree.cache:
        degree.cache[id(model)] = collections.Counter(m for r in model.rxns.values() for m in r["stoich"])
    return degree.cache[id(model)]


degree.cache = {}


def mains(model, rid):
    return [m for m in model.rxns[rid]["stoich"] if model.name(m) not in L.COFACTORS]


def usable(model, rid):
    r = model.rxns[rid]
    subs = r.get("subsystem")
    subs = set(subs if isinstance(subs, list) else [subs])
    return (r["stoich"] and not L.artefact(model, rid) and "Exchange/demand reactions" not in subs
            and len({comp_of(m) for m in r["stoich"]}) == 1)


def add_bridges(dr, scope, comps_on_map, limit, steps=1):
    """Add bridges (D14) to dr.pending or between their drawn metabolites; limit None means no hub limit."""
    model, log = dr.model, dr.log
    if not scope or not scope.subsystems or scope.subsystems & SKIP:
        return
    deg = degree(model)

    def links(m):
        return model.name(m) not in L.COFACTORS and (limit is None or deg[m] <= limit)

    drawn = set(dr.mp.reactions()) | {rid for rid, _ in dr.pending}
    drawn = {r for r in drawn if r in model.rxns}
    # metabolites that can join: drawn as main nodes, or main metabolites of reactions still to be placed
    on_map = {dr.mp.classes(n)[1] for n in dr.mp.nodes("met") if dr.is_main(n)}
    on_map |= {m for rid, _ in dr.pending if rid in model.rxns for m in mains(model, rid)}
    index = collections.defaultdict(set)  # metabolite -> reactions that could bridge through it
    for rid in model.rxns:
        if rid not in drawn and not scope.contains(rid) and usable(model, rid):
            for m in mains(model, rid):
                if links(m):
                    index[m].add(rid)
    added = 0
    if not hasattr(dr, "bridge_groups"):
        dr.bridge_groups = []
        dr.bridge_via = {}
    while True:
        parent = {r: r for r in drawn}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        first = {}
        for r in sorted(drawn):
            for m in mains(model, r):
                if links(m) and m in on_map:
                    if m in first:
                        parent[find(r)] = find(first[m])
                    else:
                        first[m] = r
        size = collections.Counter(find(r) for r in drawn)
        owner = {m: find(r) for m, r in first.items()}
        best = None
        for m, part in owner.items():
            for rid in index.get(m, ()):
                if rid in drawn:
                    continue
                touched = {owner[x] for x in mains(model, rid) if x in owner and links(x)}
                if len(touched) < 2:
                    continue
                if L.COMPARTMENT.get(comp_of(m), comp_of(m)) not in comps_on_map and not dr.one_area:
                    continue
                joined = sorted((size[p] for p in touched), reverse=True)
                key = (-sum(joined[:2]), max(deg[x] for x in mains(model, rid) if x in owner), rid)
                if best is None or key < best[0]:
                    best = (key, rid, [x for x in mains(model, rid) if x in owner])
        pair = None
        if best is None and steps == 2:
            for m, part in owner.items():
                if L.COMPARTMENT.get(comp_of(m), comp_of(m)) not in comps_on_map and not dr.one_area:
                    continue
                for r1 in index.get(m, ()):
                    if r1 in drawn:
                        continue
                    for x in mains(model, r1):
                        if x in owner or not links(x):
                            continue
                        for r2 in index.get(x, ()):
                            if r2 in drawn or r2 == r1:
                                continue
                            far = {owner[y] for y in mains(model, r2) if y in owner and links(y)} - {part}
                            if not far:
                                continue
                            key = (-(size[part] + max(size[q] for q in far)), deg[m] + deg[x], r1, r2)
                            if pair is None or key < pair[0]:
                                pair = (key, r1, r2, m, x)
        if best is None and pair is None:
            break
        if best is not None:
            _, rid, via = best
            new = [(rid, "bridge added (context reaction)", "joins parts through " +
                    ", ".join(sorted({model.name(x) for x in via})))]
        else:
            _, r1, r2, m, x = pair
            new = [(r, "two-step bridge added (context reaction)",
                    f"joins parts through {model.name(m)} and {model.name(x)} (new on the map)") for r in (r1, r2)]
        dr.bridge_groups.append([rid for rid, _, _ in new])
        if best is not None:  # drawn after placement and D10, between the nodes it joins (draw_bridges)
            dr.bridge_via[rid] = sorted(set(via))
        for rid, what, how in new:
            drawn.add(rid)
            dr.context.add(rid)
            if best is None:
                L.draw_missing(dr, rid, comps_on_map)
            added += 1
            if best is None:  # drawn now; one-step bridges are logged when they are drawn (draw_bridges)
                log.append(("D14", what, rid, how, subsystem_and_name(model, rid)))
    return added


def subsystem_and_name(model, rid):
    subs = model.rxns[rid].get("subsystem")
    subs = subs[0] if isinstance(subs, list) else subs
    return f"{subs}: {model.rxns[rid].get('name', '') or ''}"


MOVE_MAX = 25  # reactions: the largest part moved to join a bridge


def joins(dr, group):
    """True if the reactions in group join, as drawn, two parts of the other drawn reactions (other bridges
    included, so a bridge may join through a node another bridge added)."""
    parent = {}

    def find(k):
        parent.setdefault(k, k)
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k
    touching = set()
    for n in dr.mp.nodes("met"):
        if not dr.is_main(n):
            continue
        rids = [c for c in dr.mp.classes(n)[2:] if L.is_rxn(c)]
        others = [r for r in rids if r not in group]
        for r in others[1:]:
            parent[find(r)] = find(others[0])
        if len(others) < len(rids):
            touching.update(others)
    return len({find(r) for r in touching}) >= 2


def sample(pts_list, step=40):
    out = []
    for sp in pts_list:
        out.append(sp[0])
        for a, b in zip(sp, sp[1:]):
            k = max(1, int(L.length(L.sub(b, a)) / step))
            out += [L.add(a, L.sub(b, a), i / k) for i in range(1, k + 1)]
    return out


def group_elements(mp, rids):
    """Nodes, edge paths and band pieces that belong only to the reactions rids."""
    nodes = []
    for n in mp.layer["nodes"]:
        cls = mp.classes(n)
        if not cls:
            continue
        if cls[0] == "rea":
            if n.get("id") in rids:
                nodes.append(n)
            continue
        rx = [c for c in cls[2:] if L.is_rxn(c)]
        if rx and all(r in rids for r in rx):
            nodes.append(n)
    paths = []
    for layer in ("fes", "fehs", "ees"):
        for path in mp.layer[layer]:
            rx = [c for c in mp.classes(path)[1:] if L.is_rxn(c)]
            if rx and all(r in rids for r in rx):
                paths.append(path)
    pts = [L.translate(n) for n in nodes if L.translate(n)]
    pts += sample([sp for path in paths for sp in (L.subpaths(path.get("d") or "") or [])])
    grid = collections.defaultdict(list)
    for q in pts:
        grid[(int(q[0] // 60), int(q[1] // 60))].append(q)

    def near(q):
        i, j = int(q[0] // 60), int(q[1] // 60)
        return any(L.length(L.sub(q, x)) < 60 for di in (-1, 0, 1) for dj in (-1, 0, 1) for x in grid.get((i + di, j + dj), ()))
    band = []
    for _, bp in L.band_groups(mp):
        for k, sp in enumerate(L.subpaths(bp.get("d") or "") or []):
            if sample([sp], 30) and all(near(q) for q in sample([sp], 30)):
                band.append((bp, k))
    return nodes, paths, band


def footprint(mp, nodes, paths):
    boxes = [b for b in (L.node_box(mp, n) for n in nodes) if b]
    boxes += [(q[0] - 4, q[1] - 4, q[0] + 4, q[1] + 4)
              for q in sample([sp for path in paths for sp in (L.subpaths(path.get("d") or "") or [])])]
    return boxes


def obstacles(mp, skip_nodes, skip_paths):
    skip_n, skip_p = {id(n) for n in skip_nodes}, {id(p) for p in skip_paths}
    others = [n for n in mp.layer["nodes"] if id(n) not in skip_n]
    paths = [p for layer in ("fes", "fehs", "ees") for p in mp.layer[layer] if id(p) not in skip_p]
    boxes = footprint(mp, others, paths)
    for g in mp.root.iter(L.SVG + "g"):
        if g.get("class") in ("compartment", "subsystem"):
            boxes += [b for b in (L.text_box(t) for t in g.iter(L.SVG + "text")) if b]
    grid = collections.defaultdict(list)
    for b in boxes:
        for i in range(int(b[0] // 100), int(b[2] // 100) + 1):
            for j in range(int(b[1] // 100), int(b[3] // 100) + 1):
                grid[(i, j)].append(b)
    return grid


def fits(boxes, off, grid, reg, holes, margin=25):
    for b in boxes:
        x0, y0, x1, y1 = b[0] + off[0], b[1] + off[1], b[2] + off[0], b[3] + off[1]
        if x0 < reg[0] + 80 or x1 > reg[2] - 80 or y0 < reg[1] + 80 or y1 > reg[3] - 80:
            return False
        if any(x0 < h[2] + 40 and x1 > h[0] - 40 and y0 < h[3] + 40 and y1 > h[1] - 40 for h in holes):
            return False
        for i in range(int((x0 - margin) // 100), int((x1 + margin) // 100) + 1):
            for j in range(int((y0 - margin) // 100), int((y1 + margin) // 100) + 1):
                for o in grid.get((i, j), ()):
                    if x0 - margin < o[2] and x1 + margin > o[0] and y0 - margin < o[3] and y1 + margin > o[1]:
                        return False
    return True


def shift(mp, nodes, paths, band, off):
    for n in nodes:
        p = L.translate(n)
        if p:
            L.set_position(n, p[0] + off[0], p[1] + off[1])
    for path in paths:
        sps = L.subpaths(path.get("d") or "")
        if sps:
            path.set("d", L.to_d([[(q[0] + off[0], q[1] + off[1]) for q in sp] for sp in sps]))
    by_path = collections.defaultdict(set)
    for bp, k in band:
        by_path[id(bp)].add(k)
    for _, bp in L.band_groups(mp):
        ks = by_path.get(id(bp))
        if ks:
            sps = L.subpaths(bp.get("d") or "") or []
            bp.set("d", L.to_d([[(q[0] + off[0], q[1] + off[1]) for q in sp] if k in ks else sp
                                for k, sp in enumerate(sps)]))


def area(mp, comp):
    """The compartment's region (L.region), narrowed to its side of the separator lines between compartments
    drawn without boxes (the extracellular line on many subsystem maps)."""
    reg, holes = L.region(mp, comp)
    head = None
    lines = []
    for g in mp.root.iter(L.SVG + "g"):
        if g.get("class") != "compartment" or g.find(L.SVG + "rect") is not None:
            continue
        if g.get("id") == comp:
            hb = [b for b in (L.text_box(t) for t in g.iter(L.SVG + "text")) if b]
            if hb:
                head = (hb[0][0] + hb[0][2]) / 2
        for path in g.iter(L.SVG + "path"):
            for sp in L.subpaths(path.get("d") or "") or []:
                if len(sp) >= 2 and abs(sp[0][0] - sp[-1][0]) < 5:
                    lines.append(sp[0][0])
    x0, y0, x1, y1 = reg
    if head is not None:
        for x in lines:
            if head > x:
                x0 = max(x0, x + 60)
            else:
                x1 = min(x1, x - 60)
    return (x0, y0, x1, y1), holes


def move_group(dr, rids, a_node, t_node, comp, dists=(450, 600, 800, 1000, 1300, 1700, 2200), tries=6):
    """Places for the reactions rids as a whole (nodes, edges, band; the shape is kept) with a_node in free
    space dists px from t_node, inside the compartment's area, nearest first: each is applied in turn and
    yielded as (distance, what moved); the caller keeps one or lets the next undo it."""
    import math
    mp = dr.mp
    nodes, paths, band = group_elements(mp, rids)
    if not nodes:
        return
    boxes = footprint(mp, nodes, paths)
    grid = obstacles(mp, nodes, paths)
    reg, holes = area(mp, comp)
    a, t = L.translate(a_node), L.translate(t_node)
    for dist in dists:
        for k in range(24):
            ang = k * math.pi / 12
            new = (t[0] + dist * math.cos(ang), t[1] + dist * math.sin(ang))
            off = (new[0] - a[0], new[1] - a[1])
            if fits(boxes, off, grid, reg, holes):
                shift(mp, nodes, paths, band, off)
                dr.space = L.Space(mp)
                yield dist, (nodes, paths, band, off)
                shift(mp, nodes, paths, band, (-off[0], -off[1]))  # not kept: back where it was
                dr.space = L.Space(mp)
                tries -= 1
                if tries == 0:
                    return
                break  # the next distance


def compact_boxes(mp):
    import compact
    return compact.boxes(mp).items()


def part_index(dr):
    """Parts of the drawn reactions (bridges included as they are drawn so far): find and part members."""
    mp = dr.mp
    parent = {}

    def find(k):
        parent.setdefault(k, k)
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k
    for n in mp.nodes("met"):
        if dr.is_main(n):
            rids = [c for c in mp.classes(n)[2:] if L.is_rxn(c)]
            for r in rids[1:]:
                parent[find(r)] = find(rids[0])
    members = collections.defaultdict(set)
    for r in mp.reactions():
        members[find(r)].add(r)
    return find, members


def end_nodes(dr, via, find):
    """The closest pair of main nodes of the joining metabolites that lie in different parts, in one area."""
    import connect
    mp = dr.mp
    boxes = list(L.compartment_boxes(mp).values())
    cands = []
    for m in via:
        for n in mp.nodes("met"):
            if mp.classes(n)[1] == m and dr.is_main(n) and L.translate(n):
                rids = [c for c in mp.classes(n)[2:] if L.is_rxn(c)]
                if rids:
                    cands.append((m, n, find(rids[0])))
    best = None
    for i, (ma, na, pa) in enumerate(cands):
        for mb, nb, pb in cands[i + 1:]:
            if ma == mb or pa == pb:
                continue
            qa, qb = L.translate(na), L.translate(nb)
            if connect.area_of(qa, boxes) != connect.area_of(qb, boxes):
                continue
            d = L.length(L.sub(qa, qb))
            if best is None or d < best[0]:
                best = (d, (ma, na, pa), (mb, nb, pb))
    return best


def draw_bridge(dr, rid, ends, why=None, crossings=None):
    """Draw rid in free space between or beside the existing nodes in ends {metabolite: node}, with routed edges
    to them (crossing at most `crossings` edges each, BRIDGE_CROSSINGS by default); its other main metabolites
    get new nodes beside it. Places with the shortest edges are tried first. False when there is no room; why
    (a Counter) then counts what was missing at the places tried."""
    import math
    import connect
    mp, model = dr.mp, dr.model
    st = model.rxns[rid]["stoich"]
    mains_ = mains(model, rid)
    subs, prods = [m for m in mains_ if st[m] < 0], [m for m in mains_ if st[m] > 0]
    (m1, n1), (m2, n2) = list(ends.items())
    p1, p2 = L.translate(n1), L.translate(n2)
    reg, holes = area(mp, dr.one_area or L.COMPARTMENT.get(comp_of(m1), comp_of(m1)))
    inside = lambda q: (reg[0] + 60 <= q[0] <= reg[2] - 60 and reg[1] + 60 <= q[1] <= reg[3] - 60 and  # noqa: E731
                        not any(h[0] - 40 <= q[0] <= h[2] + 40 and h[1] - 40 <= q[1] <= h[3] + 40 for h in holes))
    axis = L.unit(L.sub(p2, p1)) if L.length(L.sub(p2, p1)) > 1 else (1.0, 0.0)
    if crossings is None:
        crossings = connect.BRIDGE_CROSSINGS
    cands = [L.add(L.add(p1, L.sub(p2, p1), f), L.normal(axis), off) for f in (0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8)
             for off in (0, 150, -150, 300, -300, 450, -450, 650, -650, 900, -900)]
    cands += [L.add(q, (math.cos(k * math.pi / 8), math.sin(k * math.pi / 8)), r)  # beside either end
              for q in (p1, p2) for r in (250, 400, 600) for k in range(16)]
    cands.sort(key=lambda c: round(L.length(L.sub(c, p1)) + L.length(L.sub(c, p2))))  # stable: ties keep order
    for c in cands:
        if not inside(c):
            continue
        if not dr.space.free((c[0] - 55, c[1] - 55, c[0] + 55, c[1] + 55), 4):
            if why is not None:
                why["no free spot"] += 1
            continue
        routes = {}
        for m, n, q in ((m1, n1, p1), (m2, n2, p2)):
            pts = connect.route(dr, c, q, [(c, 140), (q, 80)], crossings)
            if pts is None:
                break
            routes[m] = list(reversed(pts))
        if len(routes) != 2:
            if why is not None:
                why["no clear route"] += 1
            continue
        pos = {m1: p1, m2: p2}
        extra = [m for m in mains_ if m not in pos]
        for m in extra:  # beside the reaction, where there is room for the label
            for k in range(16):
                u = (math.cos(k * math.pi / 8), math.sin(k * math.pi / 8))
                q = L.add(c, u, L.STEP * (1 if k % 2 == 0 else 1.4))
                if inside(q) and dr.space.free((q[0] - 90, q[1] - 30, q[0] + 90, q[1] + 30), 6) and \
                        all(L.length(L.sub(q, x)) > 120 for x in pos.values()):
                    pos[m] = q
                    break
            else:
                break
        if any(m not in pos for m in extra):
            if why is not None:
                why["no room for its other metabolites"] += 1
            continue
        links = []
        for m in subs + prods:
            if m in ends:
                dr.link(ends[m], rid)
                line = routes[m]
            else:
                node = L.new_main(m, model.name(m), [rid], *pos[m])
                L.append(mp.layer["nodes"], node)
                dr.space.add(L.node_box(mp, node))
                line = None
            links.append((pos[m], "main", -1 if m in subs else 1, line))
        L.append(mp.layer["nodes"], L.new_reaction(rid, *c))
        sp = [pos[m] for m in subs] or [c]
        pp = [pos[m] for m in prods] or [c]
        tangent = L.unit(L.sub(pp[0], sp[0])) if L.length(L.sub(pp[0], sp[0])) > 1 else axis
        side = dr.choose_side(c, tangent)
        links += dr.add_side_mets(rid, c, tangent, side, set(mains_))
        dr.draw(rid, c, links, (-side[0], -side[1]))
        dr.placed.add(rid)
        return True
    return False


LONG_BRIDGE = 3000     # px: the furthest apart two parts may be for a bridge drawn without moving either
LONG_CROSSINGS = 3     # edges each line of such a bridge may cross


def draw_bridges(dr, near=1500):
    """D14 after placement and D10: each one-step bridge is drawn between the nodes it joins. When they are
    more than `near` px apart, or there is no room between them, the smaller of the two parts, else the larger
    (at most MOVE_MAX reactions), is moved as a whole next to the other first. When neither can move, a bridge
    whose parts are at most LONG_BRIDGE px apart is drawn where they are, over lines that cross at most
    LONG_CROSSINGS edges each. A bridge that cannot be drawn is dropped, with the reason."""
    mp, model, log = dr.mp, dr.model, dr.log
    for rid, via in sorted(getattr(dr, "bridge_via", {}).items()):
        find, members = part_index(dr)
        pair = end_nodes(dr, via, find)
        if pair is None:  # the parts are joined already
            dr.context.discard(rid)
            continue
        d, (ma, na, pa), (mb, nb, pb) = pair
        ends = {ma: na, mb: nb}
        why = collections.Counter()
        done = d <= near and draw_bridge(dr, rid, ends, why)
        sizes = sorted((len(members[pa]), len(members[pb])))
        moves = 0
        parts = sorted(((pa, na, nb, ma), (pb, nb, na, mb)), key=lambda x: len(members[x[0]]))  # smaller first
        for root, a_node, t_node, m in parts:
            if done or len(members[root]) > MOVE_MAX:
                continue
            comp = dr.one_area or L.COMPARTMENT.get(comp_of(m), comp_of(m))
            for dist, (nodes, paths, band, off) in move_group(dr, members[root], a_node, t_node, comp):
                moves += 1
                done = draw_bridge(dr, rid, ends, why)
                if done:
                    # D6 may close the space the part left in its box
                    for g, b in compact_boxes(mp):
                        if any(b[0] <= L.translate(n)[0] <= b[2] and b[1] <= L.translate(n)[1] <= b[3]
                               for n in nodes[:1] if L.translate(n)):
                            g.set("data-grown", "1")
                    log.append(("D14", "part moved to join a bridge", rid,
                                f"{len(members[root])} reaction(s), {dist} px from the other copy", model.name(m)))
                    break  # kept here; otherwise move_group puts it back and tries the next place
        if not done and near < d <= LONG_BRIDGE:
            done = draw_bridge(dr, rid, ends, why, LONG_CROSSINGS)
            if done:
                log.append(("D14", "bridge drawn over long lines", rid, f"{d:.0f} px between the parts", ""))
        if not done:
            if sizes[0] > MOVE_MAX and d > LONG_BRIDGE:
                reason = "parts too large to move and too far apart"
            elif not moves and d > LONG_BRIDGE:
                reason = "no free place to move either part, and too far apart"
            else:
                reason = ", ".join(f"{k} ({v})" for k, v in why.most_common()) or "no free place to move either part"
            log.append(("review", "bridge left out: no room to draw it", rid,
                        f"{d:.0f} px between parts of {sizes[0]} and {sizes[1]} reactions; {reason}; "
                        f"{moves} place(s) tried for a part", ", ".join(model.name(x) for x in via)))
            dr.context.discard(rid)
            continue
        log.append(("D14", "bridge drawn (context reaction)", rid, "between " + model.name(ma) + " and " + model.name(mb),
                    subsystem_and_name(model, rid)))


def check_bridges(dr):
    """After placement, D10 and drawing: a bridge stays only if, as drawn, it joins two parts that are not
    joined without it; otherwise it is removed again (in the order they were added)."""
    mp, log = dr.mp, dr.log
    for g in [g for g in getattr(dr, "bridge_groups", []) if any(r in mp.reactions() for r in g)]:
        if joins(dr, set(g)):
            continue
        for r in g:
            if r in mp.reactions():
                L.erase_reaction(dr, r)
            dr.context.discard(r)
            log.append(("D14", "bridge removed: not joined to two parts on the drawing", r, "", ""))
