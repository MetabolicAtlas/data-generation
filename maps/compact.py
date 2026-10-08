"""D11-D13: compartment boxes and white space.

D11 a compartment keeps one box: reactions added to a boxed compartment that has no free room go into
    that box, which grows into free space next to it, or moves as a whole, made larger, to free space or
    beside the drawing.
D12 a compartment box with nothing in it is removed.
D13 white space: a stripe across the whole map (top to bottom or left to right) with nothing in it is
    narrowed, the canvas is cut to the content, and the title of a subsystem map is centred again; boxes
    moved or added in the run look again for free space at their final size.
"""

import math
import re

import layout as L
from mapedit import SVG, subpaths, to_d, translate

MARGIN = 200   # px kept between the content and the canvas edge
KEEP = 250     # px an empty stripe is narrowed to
TITLE_GAP = 60  # px between the title and a drawing moved down from under it
MIN_GAP = 400  # narrower empty stripes are left alone


# ----------------------------------------------------------------------------- moving everything

def boxes(mp):
    """Compartment group -> rectangle, for compartments drawn as boxes."""
    out = {}
    for g in mp.root.iter(SVG + "g"):
        if g.get("class") == "compartment":
            r = g.find(SVG + "rect")
            if r is not None and r.get("width"):
                x, y = float(r.get("x")), float(r.get("y"))
                out[g] = (x, y, x + float(r.get("width")), y + float(r.get("height")))
    return out


def text_pos(t):
    m = re.match(r"matrix\(1,0,0,1,([-0-9.e]+),([-0-9.e]+)\)", t.get("transform") or "")
    if m:
        return float(m.group(1)), float(m.group(2))
    if t.get("x"):
        return float(t.get("x")), float(t.get("y") or 0)
    return None


def set_text_pos(t, p):
    if (t.get("transform") or "").startswith("matrix(1,0,0,1,"):
        t.set("transform", "matrix(1,0,0,1,%s,%s)" % (L.fmt(p[0]), L.fmt(p[1])))
    else:
        t.set("x", L.fmt(p[0]))
        t.set("y", L.fmt(p[1]))


def free_texts(mp):
    """Texts outside the nodes: compartment headings, subsystem labels and the title."""
    out = []
    for g in mp.root.iter(SVG + "g"):
        if g.get("class") in ("compartment", "subsystem"):
            out += [t for t in g.iter(SVG + "text")]
    return out


def apply(mp, f, rect_f, skip_texts=()):
    """Move every element: nodes, edges, arrowheads, bands and separators by the point map f, boxes by
    rect_f (rectangle -> rectangle), free texts by f."""
    for n in list(mp.layer["nodes"]):
        p = translate(n)
        cls = mp.classes(n)
        if not p or not cls:
            continue
        off = (30, 15) if cls[0] == "enz" else (0, 0)  # a gene box is placed by its corner
        q = f((p[0] + off[0], p[1] + off[1]))
        if q != (p[0] + off[0], p[1] + off[1]):
            L.set_position(n, q[0] - off[0], q[1] - off[1])
    paths = [path for layer in ("fes", "fehs", "ees") for path in mp.layer[layer]]
    for g in mp.root.iter(SVG + "g"):
        if g.get("class") in ("compartment", "subsystem"):
            paths += list(g.iter(SVG + "path"))
    for path in paths:
        sps = subpaths(path.get("d") or "")
        if sps:
            path.set("d", to_d([[f(q) for q in sp] for sp in sps]))
    for g, b in boxes(mp).items():
        r = g.find(SVG + "rect")
        nb = rect_f(g, b)
        r.set("x", L.fmt(nb[0]))
        r.set("y", L.fmt(nb[1]))
        r.set("width", L.fmt(nb[2] - nb[0]))
        r.set("height", L.fmt(nb[3] - nb[1]))
    for t in free_texts(mp):
        if t in skip_texts:
            continue
        p = text_pos(t)
        if p:
            set_text_pos(t, f(p))


def set_canvas(mp, w, h):
    """Canvas size; the licence badge stays in the bottom right corner."""
    dw, dh = w - float(mp.root.get("width")), h - float(mp.root.get("height"))
    mp.root.set("width", L.fmt(w))
    mp.root.set("height", L.fmt(h))
    for g in mp.root.iter(SVG + "g"):
        if g.get("id") == "license":
            m = re.match(r"translate\(([-0-9.e]+),([-0-9.e]+)\)", g.get("transform") or "")
            dx, dy = (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0)
            g.set("transform", "translate(%s,%s)" % (L.fmt(dx + dw), L.fmt(dy + dh)))


def inside(p, b):
    return b[0] <= p[0] <= b[2] and b[1] <= p[1] <= b[3]


# ----------------------------------------------------------------------------- D11

def rows(blocks, width, pad, gap):
    """Blocks laid out left to right in rows of at most `width`: [(block, x, y, bbox)], total w, total h."""
    out, x, y, row_h, used_w = [], 0.0, 0.0, 0.0, 0.0
    for b in blocks:
        bb = b.bbox()
        w, h = bb[2] - bb[0] + 2 * pad, bb[3] - bb[1] + 2 * pad
        if x and x + w > width:
            x, y, row_h = 0.0, y + row_h + gap, 0.0
        out.append((b, x, y, bb))
        used_w = max(used_w, x + w)
        x += w + gap
        row_h = max(row_h, h)
    return out, used_w, y + row_h


def grow_box(dr, comp, blocks, draw, orient, pad=180, gap=200):
    """D11: put blocks that found no room in compartment comp's box into that box, made larger."""
    mp = dr.mp
    g = next((g for g in boxes(mp) if g.get("id") in (comp, L.ALIASES.get(comp))), None)
    if g is None:
        return False
    x0, y0, x1, y1 = boxes(mp)[g]
    inner = x1 - x0 - 300
    widest = max(b.bbox()[2] - b.bbox()[0] + 2 * pad for b in blocks)
    down, w_down, h_down = rows(blocks, max(inner, widest, squarish(blocks, pad)), pad, gap)
    # to the right: rows as wide as needed to fit the box height, if possible
    height = y1 - y0 - 300
    width = widest
    right, w_right, h_right = rows(blocks, width, pad, gap)
    while h_right > height and width < 20 * widest:
        width *= 1.25
        right, w_right, h_right = rows(blocks, width, pad, gap)
    cw, ch = L.canvas(mp)[2:]
    occ_free = lambda b: dr.space.free(b, 30)  # noqa: E731
    below = (x0, y1, x0 + max(w_down + 300, x1 - x0), y1 + h_down + 100)
    beside = (x1, y0, x1 + w_right + 150, y0 + max(h_right + 300, y1 - y0))
    if occ_free(below) and below[2] <= cw - 100:
        how, layout, start = "down", down, (x0 + 150, y1 - 50)
        dx, dy = max(0.0, below[2] - x1), h_down + 100
        if below[3] > ch - MARGIN:
            set_canvas(mp, cw, below[3] + MARGIN)
        grow_rect(g, dx, dy)
    elif occ_free(beside) and beside[3] <= ch - 100:
        how, layout, start = "right", right, (x1 - 50, y0 + 250)
        dx, dy = w_right + 150, max(0.0, beside[3] - y1)
        if beside[2] > cw - MARGIN:
            set_canvas(mp, beside[2] + MARGIN, ch)
        grow_rect(g, dx, dy)
    else:
        # no room around the box: the box moves as a whole, made larger, to free space or below the drawing
        nw, nh = max(x1 - x0, w_down + 300), (y1 - y0) + h_down + 100
        at = free_box_spot(dr, nw, nh, skip=g)
        how = "moved to free space"
        if at is not None and at[1] + nh + 300 > ch:
            set_canvas(mp, cw, at[1] + nh + 300)
        if at is None:
            # below the drawing or right of it, whichever leaves the smaller canvas closer to landscape
            top = min([b[1] for b in boxes(mp).values()] or [300.0])
            options = [("moved below the drawing", (200.0, ch - 100), (max(cw, 200 + nw + MARGIN), ch + nh + 200)),
                       ("moved right of the drawing", (cw - 100, top), (cw + nw + 200, max(ch, top + nh + 300)))]
            how, at, size = min(options, key=lambda o: o[2][0] * o[2][1] * (1 + abs(math.log(o[2][0] / o[2][1] / 1.5))))
            set_canvas(mp, *size)
        L.move_box(mp, g, at[0] - x0, at[1] - y0)
        grow_rect(g, nw - (x1 - x0), nh - (y1 - y0))
        g.set("data-moved", "1")
        layout, start = down, (at[0] + 150, at[1] + (y1 - y0) - 50)
        dr.space = L.Space(mp)
    for b, x, y, bb in layout:
        t = (start[0] + x + pad - bb[0], start[1] + y + pad - bb[1])
        draw(dr, b, orient, t, {}, f"in the {comp} box, made larger ({how})")
    g.set("data-grown", "1")  # D6 may close any empty stripe in it, not only stripes that held content before
    dr.log.append(("D11", "compartment box made larger", comp, how, f"{len(blocks)} block(s)"))
    return True


def squarish(blocks, pad):
    """Row width that makes a layout of the blocks about as wide as it is high."""
    area = sum((b.bbox()[2] - b.bbox()[0] + 2 * pad) * (b.bbox()[3] - b.bbox()[1] + 2 * pad) for b in blocks)
    return 1.6 * area ** 0.5


def free_box_spot(dr, w, h, step=100, margin=80, skip=None):
    """Top-left corner of a free w x h place on the canvas, next to the compartment boxes if possible."""
    mp = dr.mp
    cw, ch = L.canvas(mp)[2:]
    bx = [b for g, b in boxes(mp).items() if g is not skip]
    best = None
    y = min([b[1] for b in bx] or [200.0])  # not above the boxes (title, headings)
    while y <= ch - 300:  # a place may reach below the canvas (it grows); growth costs
        x = 200.0
        while x + w <= cw - 200:
            cand = (x, y, x + w, y + h)
            clear = free_of(dr, cand, margin, skip) if skip is not None else dr.space.free(cand, margin)
            if clear and not any(cand[0] < b[2] + margin and cand[2] > b[0] - margin and
                                 cand[1] < b[3] + margin and cand[3] > b[1] - margin for b in bx):
                gap = min([max(b[0] - cand[2], cand[0] - b[2], b[1] - cand[3], cand[1] - b[3], 0) for b in bx]
                          or [0])
                score = gap + 0.2 * y + 1.5 * max(0.0, y + h + 300 - ch)
                if best is None or score < best[0]:
                    best = (score, x, y)
            x += step
        y += step
    return best[1:] if best else None


def free_of(dr, cand, margin, skip):
    """cand is clear of the drawing, ignoring what lies inside the box `skip` (which is about to move)."""
    if skip is None:
        return True
    own = boxes(dr.mp)[skip]
    x0, y0, x1, y1 = cand[0] - margin, cand[1] - margin, cand[2] + margin, cand[3] + margin
    for cell in dr.space.cells((x0, y0, x1, y1)):
        for b in dr.space.grid.get(cell, ()):
            if b[0] <= x1 and b[2] >= x0 and b[1] <= y1 and b[3] >= y0 and not (
                    own[0] <= b[0] and b[2] <= own[2] and own[1] <= b[1] and b[3] <= own[3]):
                return False
    return True


def grow_rect(g, dx, dy):
    r = g.find(SVG + "rect")
    r.set("width", L.fmt(float(r.get("width")) + dx))
    r.set("height", L.fmt(float(r.get("height")) + dy))


def box_score(b, others):
    gap = min([max(o[0] - b[2], b[0] - o[2], o[1] - b[3], b[1] - o[3], 0) for o in others] or [0])
    return gap + 0.2 * b[1]


def settle_boxes(mp, log):
    """D13: a box moved or added in this run looks again for a free place at its final size (D6 and D10
    can make it smaller after it was placed), next to the other boxes and as high as possible."""
    class Probe:  # what free_box_spot needs
        pass
    for g in [g for g in boxes(mp) if g.attrib.pop("data-moved", None) is not None]:
        probe = Probe()
        probe.mp, probe.space = mp, L.Space(mp)
        b = boxes(mp)[g]
        w, h = b[2] - b[0], b[3] - b[1]
        others = [o for x, o in boxes(mp).items() if x is not g]
        at = free_box_spot(probe, w, h, skip=g)
        if at is None or at[1] + h + 300 > float(mp.root.get("height")):
            continue
        if box_score((at[0], at[1], at[0] + w, at[1] + h), others) < box_score(b, others) - 100:
            L.move_box(mp, g, at[0] - b[0], at[1] - b[1])
            log.append(("D13", "compartment box moved to free space", g.get("id"), f"to {at[0]:.0f}, {at[1]:.0f}", ""))


# ----------------------------------------------------------------------------- D12

def remove_empty_boxes(mp, log):
    remove_empty_headings(mp, log)
    for g, b in list(boxes(mp).items()):
        if any(translate(n) and inside(translate(n), b) for n in mp.layer["nodes"]):
            continue
        log.append(("D12", "empty compartment box removed", g.get("id"), "", ""))
        g.getparent().remove(g)


def remove_empty_headings(mp, log):
    """A compartment drawn without a box (a heading, perhaps a separator line) whose metabolites are not
    drawn anywhere outside the boxes."""
    from mapedit import comp_of
    bx = list(boxes(mp).values())
    names = set()
    for n in mp.nodes("met"):
        p = translate(n)
        if p and not any(inside(p, b) for b in bx):
            letter = comp_of(mp.classes(n)[1])
            names.add(L.COMPARTMENT.get(letter, letter))
    names |= {k for k, v in L.ALIASES.items() if v in names} | {L.ALIASES.get(k) for k in names if k in L.ALIASES}
    for g in list(mp.root.iter(SVG + "g")):
        if g.get("class") == "compartment" and g.find(SVG + "rect") is None and g.get("id") not in names:
            log.append(("D12", "heading of an empty compartment removed", g.get("id"), "", ""))
            g.getparent().remove(g)


# ----------------------------------------------------------------------------- D13

def content(mp, title):
    """Boxes of everything drawn, except the licence badge (and the title, returned apart)."""
    out = []
    for n in mp.layer["nodes"]:
        b = L.node_box(mp, n)
        if b:
            out.append(b)
    for layer in ("fes", "fehs", "ees"):
        for path in mp.layer[layer]:
            for sp in subpaths(path.get("d") or "") or []:
                out += [(q[0] - 4, q[1] - 4, q[0] + 4, q[1] + 4) for q in sp]
    for g in mp.root.iter(SVG + "g"):
        if g.get("class") in ("compartment", "subsystem"):
            for path in g.iter(SVG + "path"):
                band = g.get("class") == "subsystem"
                if not band:
                    continue  # a separator line between compartments stretches with the content
                for sp in subpaths(path.get("d") or "") or []:
                    half = 45 if band else 12
                    out += [(q[0] - half, q[1] - half, q[0] + half, q[1] + half) for q in sp]
                    for a, c in zip(sp, sp[1:]):  # long straight pieces: their ends are not enough
                        k = max(1, int(L.length(L.sub(c, a)) / 200))
                        for i in range(1, k):
                            q = L.add(a, L.sub(c, a), i / k)
                            out.append((q[0] - half, q[1] - half, q[0] + half, q[1] + half))
            for t in g.iter(SVG + "text"):
                if t is not title:
                    b = L.text_box(t)
                    if b:
                        out.append(b)
    out += list(boxes(mp).values())
    return out


def compact(mp, log, title=None, rounds=3):
    """D13: narrow empty stripes across the whole map, cut the canvas to the content, centre the title; again until
    nothing changes (moving the drawing from under the title can leave a margin to cut)."""
    for _ in range(rounds):
        if not compact_once(mp, log, title):
            return


def compact_once(mp, log, title=None):
    """One D13 round; False when it changes nothing."""
    items = content(mp, title)
    if not items:
        return False
    tb = L.text_box(title) if title is not None else None
    cw, ch = float(mp.root.get("width")), float(mp.root.get("height"))
    cuts = {}
    for axis in (0, 1):
        spans = sorted((b[axis], b[axis + 2]) for b in items + ([tb] if tb and axis == 1 else []))
        free, cur = [], None
        for a, b in spans:
            if cur is not None and a > cur:
                free.append((cur, a))
            cur = b if cur is None else max(cur, b)
        lead = spans[0][0]
        red = [(a, b, (b - a) - KEEP) for a, b in free if b - a >= MIN_GAP]
        if lead > MARGIN + 100:
            red.insert(0, (0.0, lead, lead - MARGIN))
        cuts[axis] = (red, cur)

    def shift(v, axis):
        d = 0.0
        for a, b, r in cuts[axis][0]:
            if v >= b:
                d += r
            elif v > a:
                d += r * (v - a) / (b - a)
        return v - d

    def f(p):
        return (shift(p[0], 0), shift(p[1], 1))

    def rect_f(g, b):
        return (shift(b[0], 0), shift(b[1], 1), shift(b[2], 0), shift(b[3], 1))
    w = shift(cuts[0][1], 0) + MARGIN + 100
    h = shift(cuts[1][1], 1) + MARGIN + 150  # room for the licence badge
    p = text_pos(title) if title is not None else None
    dx = 0.0
    if p:
        need = min(cw, (tb[2] - tb[0]) + 2 * MARGIN)  # the old canvas held the title
        left, right = shift(min(b[0] for b in items), 0), shift(cuts[0][1], 0)
        if need > (right - left) + 2 * MARGIN + 100:
            # content narrower than the title: centred under it, wherever its left edge was before (as if its margin
            # were MARGIN; the licence badge side has 100 more)
            dx = (need - (right - left) - 100) / 2 - left
            w = need
        else:
            w = max(w, need)
    # the title sits above the drawing: content that would lie under it, where it ends up, moves down with the rest
    down = 0.0
    if p:
        tw = tb[2] - tb[0]
        tx0, tx1, ty0, ty1 = w / 2 - tw / 2, w / 2 + tw / 2, shift(tb[1], 1), shift(tb[3], 1)
        for b in items:
            (x0, y0), (x1, y1) = f((b[0], b[1])), f((b[2], b[3]))
            if x0 + dx < tx1 and x1 + dx > tx0 and y0 < ty1 and y1 > ty0:  # under the title, not just near it
                down = max(down, ty1 + TITLE_GAP - y0)
        down = down if down > 1 else 0.0  # rounding is not a change
        h += down
    # nothing to do when the content would move (cuts and centring together) and the canvas change by less than
    # a pixel: rounding is not a change
    if abs(w - cw) <= 1 and abs(h - ch) <= 1 and not down and all(
            abs(f(q)[0] + dx - q[0]) <= 1 and abs(f(q)[1] - q[1]) <= 1 for b in items for q in ((b[0], b[1]), (b[2], b[3]))):
        return False
    apply(mp, f, rect_f, skip_texts=(title,) if title is not None else ())
    if dx or down:
        apply(mp, lambda q: (q[0] + dx, q[1] + down), lambda g, b: (b[0] + dx, b[1] + down, b[2] + dx, b[3] + down),
              skip_texts=(title,))
    if down:
        log.append(("D13", "drawing moved down from under the title", "", f"{down:.0f} px", ""))
    if p:
        half = (tb[2] - tb[0] - 20) / 2
        x = w / 2 if title.get("text-anchor") == "middle" or (title.getparent() is not None and
                                                              title.getparent().get("text-anchor") == "middle") \
            else w / 2 - half
        set_text_pos(title, (x, shift(p[1], 1)))
    nw, nh = w, h
    if abs(dx) > 1 and not any(cuts[a][0] for a in (0, 1)) and abs(nw - cw) <= 1 and abs(nh - ch) <= 1:
        log.append(("D13", "content centred under the title again", "", f"{dx:.0f} px", ""))
    elif abs(nw - cw) > 1 or abs(nh - ch) > 1 or any(cuts[a][0] for a in (0, 1)):  # not for rounding differences
        set_canvas(mp, nw, nh)
        # separator lines between compartments end where the content ends
        last_y = shift(cuts[1][1], 1) + down
        for g in mp.root.iter(SVG + "g"):
            if g.get("class") == "compartment":
                for path in g.iter(SVG + "path"):
                    sps = subpaths(path.get("d") or "")
                    if sps:
                        path.set("d", to_d([[(q[0], min(q[1], last_y)) for q in sp] for sp in sps]))
        nx = sum(r for _, _, r in cuts[0][0])
        ny = sum(r for _, _, r in cuts[1][0])
        log.append(("D13", "white space closed", "", f"{nx:.0f} px narrower inside, {ny:.0f} px lower inside",
                    f"canvas {cw:.0f} x {ch:.0f} -> {nw:.0f} x {nh:.0f}"))
    return True
