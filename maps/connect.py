"""D10: connect the drawing. A metabolite (not a cofactor) drawn as main nodes in parts of the map that are not
connected to each other gets one node: the reactions at one copy are linked to the other copy by a clear
edge (at most two bends, crossing at most LINK_CROSSINGS edges; of the clear ones the shortest, with each
crossed edge counted as CROSSING_COST px), and the emptied copy is removed. Unconnected parts come first,
smallest first. A part that cannot be linked, made only of reactions added in this run and small, is moved
once: its reactions are placed again (D9) next to the other copy, or as close to it as there is room (the map
or its box grows when there is none), so that a later pass can link it. Each link joins two parts, which can
make new pairs possible, so it runs in passes until a pass adds no link (at most `passes`); pairs that still
cannot be linked cleanly are listed for review. Hub metabolites (in more than `hub` reactions of the model)
are left as drawn: maps show them several times on purpose; the metabolites of bridges (D14) are linked all
the same, over longer lines (BRIDGE_REACH, crossing at most BRIDGE_CROSSINGS edges).
"""

import collections

import layout as L


def components(dr):
    """Union-find over drawn reactions and the main nodes linked to them (keyed by the node elements, which the
    table keeps alive: an lxml element's id() changes once nothing refers to it); returns find, size by root."""
    mp = dr.mp
    parent = {}

    def find(k):
        parent.setdefault(k, k)
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    for n in mp.nodes("met"):
        if dr.is_main(n):
            find(n)
            for rid in mp.classes(n)[2:]:
                if L.is_rxn(rid):
                    union(n, rid)
    size = collections.Counter(find(r) for r in mp.reactions())
    return find, union, size


def area_of(p, boxes):
    return next((i for i, b in enumerate(boxes) if b[0] <= p[0] <= b[2] and b[1] <= p[1] <= b[3]), -1)


def crossed(dr, pts, skip, margin=8):
    """Edges the polyline crosses (60 px cells where it meets one), or None if it passes a node or label other
    than those in skip."""
    cross = set()
    for a, b in zip(pts, pts[1:]):
        k = max(1, int(L.length(L.sub(b, a)) / 20))
        for i in range(k + 1):
            q = L.add(a, L.sub(b, a), i / k)
            if any(L.length(L.sub(q, s)) < r for s, r in skip):
                continue
            for cell in dr.space.cells((q[0] - margin, q[1] - margin, q[0] + margin, q[1] + margin)):
                for bx in dr.space.grid.get(cell, ()):
                    if bx[0] - margin <= q[0] <= bx[2] + margin and bx[1] - margin <= q[1] <= bx[3] + margin:
                        if bx[2] - bx[0] > 1 or bx[3] - bx[1] > 1:
                            return None       # a node, label or gene box
                        cross.add((round(q[0] / 60), round(q[1] / 60)))
    return len(cross)


CROSSING_COST = 600  # px of extra line worth one crossed edge less


def route(dr, frm, to, skip, crossings=3):
    """The best clear path from frm to to, of: straight when nearly aligned, else one bend (either way round);
    two bends, across at some height (or along at some x) between or just beyond the two ends. Clear: it passes
    no node or label and crosses at most `crossings` edges; best: the shortest, with each crossed edge counted
    as CROSSING_COST px."""
    shortest = []  # options as short as a path with bends at right angles can be
    if abs(frm[0] - to[0]) < 30 or abs(frm[1] - to[1]) < 30:
        shortest.append([frm, to])
    else:
        shortest += [[frm, (frm[0], to[1]), to], [frm, (to[0], frm[1]), to]]
    longer = []
    for f in (0.5, 0.25, 0.75, -0.2, 1.2):
        my = frm[1] + f * (to[1] - frm[1])
        mx = frm[0] + f * (to[0] - frm[0])
        (shortest if 0 <= f <= 1 else longer).extend([[frm, (frm[0], my), (to[0], my), to],
                                                      [frm, (mx, frm[1]), (mx, to[1]), to]])
    best = None
    for pts in shortest + longer:
        if best is not None and L.poly_length(pts) >= best[0]:
            continue
        n = crossed(dr, pts, skip)
        if n is None or n > crossings:
            continue
        cost = L.poly_length(pts) + CROSSING_COST * n
        if best is None or cost < best[0]:
            best = (cost, pts)
            if n == 0 and pts in shortest:
                break  # nothing is shorter
    return best[1] if best else None


LINK_CROSSINGS = 2       # edges a D10 link may cross
BRIDGE_REACH = 6000      # px a bridge's metabolite may be linked over (D14)
BRIDGE_CROSSINGS = 6     # edges such a link may cross


def connect(dr, reach=1500, passes=8, hub=30, movable=8):
    total, failed, moved = 0, [], set()
    for k in range(passes):
        dr.space = L.Space(dr.mp)
        n, failed = one_pass(dr, reach, hub)
        total += n
        mv = move_parts(dr, failed, moved, movable)
        if not n and not mv:
            break
    biggest = max(components(dr)[2].values() or [0])
    for m, su, d, _, _ in failed:
        if su < biggest:
            dr.log.append(("review", "part not connected: no clear link", m, f"{su} reaction(s), {d:.0f} px away",
                           dr.model.name(m)))
    if total:
        dr.log.append(("D10", "connection passes", "", f"{total} links in {k + 1} pass(es)", ""))
    return total


def one_pass(dr, reach, hub):
    mp, model, log = dr.mp, dr.model, dr.log
    boxes = list(L.compartment_boxes(mp).values())
    find, union, size = components(dr)
    rpos = {rid: L.translate(g) for rid, g in mp.reactions().items()}
    if not hasattr(model, "degree"):  # reactions per metabolite: hubs are drawn several times on purpose
        model.degree = collections.Counter(m for r in model.rxns.values() for m in r["stoich"])
    by_met = collections.defaultdict(list)
    # the metabolites of a bridge (D14) are linked even when they are hubs: joining parts is its purpose
    bridged = {m for r in getattr(dr, "context", ()) if r in model.rxns for m in model.rxns[r]["stoich"]}
    for n in mp.nodes("met"):
        m = mp.classes(n)[1]
        if dr.is_main(n) and L.translate(n) and model.name(m) not in L.COFACTORS and (model.degree[m] <= hub or
                                                                                      m in bridged):
            by_met[mp.classes(n)[1]].append(n)
    biggest = max(size.values()) if size else 0
    cands = []
    for m, ns in by_met.items():
        for i, u in enumerate(ns):
            for v in ns[i + 1:]:
                pu, pv = L.translate(u), L.translate(v)
                if find(u) == find(v) or area_of(pu, boxes) != area_of(pv, boxes):
                    continue
                d = L.length(L.sub(pu, pv))
                bridge = any(r in getattr(dr, "context", ()) for r in mp.classes(u)[2:] + mp.classes(v)[2:])
                if d > (BRIDGE_REACH if bridge else reach):
                    continue
                su, sv = size[find(u)], size[find(v)]
                if su > sv:
                    u, v, su, sv = v, u, sv, su  # u: the copy in the smaller part, which moves to v
                cands.append((su >= biggest, su, d, m, u, v))
    cands.sort(key=lambda c: c[:3])
    linked, failed = 0, []
    for _, su, d, m, u, v in cands:
        if u.getparent() is None or v.getparent() is None or find(u) == find(v):
            continue
        pu, pv = L.translate(u), L.translate(v)
        rids = [r for r in mp.classes(u)[2:] if L.is_rxn(r) and r in rpos and r in model.rxns]
        paths = {}
        bridge = any(r in getattr(dr, "context", ()) for r in mp.classes(u)[2:] + mp.classes(v)[2:])
        for r in rids:
            pts = route(dr, rpos[r], pv, [(rpos[r], 140), (pv, 80), (pu, 80)],
                        BRIDGE_CROSSINGS if bridge else LINK_CROSSINGS)
            if pts is None:
                break
            paths[r] = pts
        if len(paths) != len(rids) or not rids:
            failed.append((m, su, d, u, v))
            continue
        for r, pts in paths.items():
            st = model.rxns[r]["stoich"]
            role = 1 if st.get(m, 0) > 0 else -1
            rev = (model.rxns[r].get("lower_bound") or 0) < 0
            # the edge pieces ending at the emptied copy go
            for layer in ("fes", "fehs"):
                for p in mp.edges(layer, r):
                    sps = L.subpaths(p.get("d") or "") or []
                    rad = 45
                    left = [s for s in sps if not (L.length(L.sub(s[0], pu)) < rad or L.length(L.sub(s[-1], pu)) < rad)]
                    if len(left) != len(sps):
                        if layer == "fes":
                            dr.gone += L.segments([s for s in sps if s not in left])
                        p.set("d", L.to_d(left))
            # the new edge, from the node to the reaction
            line = list(reversed(pts))
            line[0] = L.boundary(pv, line[1], "main")
            line[-1] = L.boundary(rpos[r], line[-2], "rea")
            fe = mp.edges("fes", r)
            if fe:
                fe[0].set("d", L.to_d((L.subpaths(fe[0].get("d") or "") or []) + [line]))
            if role > 0 or rev:
                feh = mp.edges("fehs", r)
                if feh:
                    feh[0].set("d", L.to_d((L.subpaths(feh[0].get("d") or "") or []) + [L.arrow(line[0], line[1])]))
            dr.space.add_edges([line])
            dr.band_add.append((r, [pv] + line[1:-1] + [rpos[r]]))
            dr.link(v, r)
            union(v, r)
        mp.remove(u)
        linked += 1
        log.append(("D10", "metabolite drawn twice, linked to one node", m, f"{len(rids)} reaction(s), {d:.0f} px",
                    model.name(m)))
    return linked, failed


def move_parts(dr, failed, moved, movable):
    """Place small unconnected parts of added reactions again, next to the other copy of their metabolite."""
    import placement as P
    mp, model = dr.mp, dr.model
    find, _, size = components(dr)
    count = 0
    for m, su, d, u, v in failed:
        if u.getparent() is None or v.getparent() is None or find(u) == find(v):
            continue
        root = find(u)
        rids = tuple(sorted(r for r in mp.reactions() if find(r) == root))
        if not rids or len(rids) > movable or rids in moved or not all(r in dr.placed for r in rids):
            continue
        moved.add(rids)  # by its reactions: once moved, the part has new nodes
        rids = list(rids)
        target = L.translate(v)
        old = [L.translate(g) for r, g in mp.reactions().items() if r in rids]
        home = (sum(p[0] for p in old) / len(old), sum(p[1] for p in old) / len(old))
        for r in rids:
            L.erase_reaction(dr, r)
        dr.band_add = [x for x in dr.band_add if x[0] not in rids]  # their band goes with them
        dr.space = L.Space(mp)
        dr.occ = P.Occupancy(dr)
        index = P.node_index(dr)
        comp = dr.one_area or L.COMPARTMENT.get(L.comp_of(m), "")
        reg, holes = L.region(mp, comp)
        for b in sorted((P.chain_block(model, ch) for ch in L.chains(model, rids)), key=lambda b: -len(b.steps)):
            dr.occ.update(dr)
            dk = P.dock(dr, b, reg, holes, index)
            if dk:
                created = P.draw(dr, b, dk[0], dk[1], dk[2], "moved onto its other copy (D10)")
            else:
                f = P.free_spot(dr, b, reg, holes, target) or P.free_spot(dr, b, reg, holes, home)
                if f:
                    created = P.draw(dr, b, f[0], f[1], {}, "moved next to its other copy (D10)")
                elif comp in L.compartment_boxes(mp) or L.ALIASES.get(comp) in L.compartment_boxes(mp):
                    import compact  # no room: its box grows (D11), as for a block D9 finds no room for
                    compact.grow_box(dr, comp, [b], P.draw, P.ORIENT[0])
                    created = {}
                    index = P.node_index(dr)
                    dr.occ = P.Occupancy(dr)
                else:  # no room: the canvas grows (D9 step 3)
                    created = P.below_drawing(dr, b, comp, home, "moved, below the drawing (D10)")
                reg, holes = L.region(mp, comp)  # the box or the canvas grew
                if created is None:
                    dr.log.append(("review", "part could not be placed again", b.steps[0][0], comp, ""))
                    continue
            for n in dict.fromkeys(created.values()):  # drawing order, not memory order
                q = L.translate(n)
                if q:
                    index[mp.classes(n)[1]].append((q, n))
        dr.occ = None
        count += 1
        dr.log.append(("D10", "part moved next to its other copy", m, f"{len(rids)} reaction(s)", model.name(m)))
    return count
