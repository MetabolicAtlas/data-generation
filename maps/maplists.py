"""Adapt the lists of maps, subsystemSVG.tsv and compartmentSVG.tsv, to the subsystems and compartments of a new release.

A listed subsystem or compartment that the new model lacks loses its row and its SVG. One that is in the new model but
not in the old one, and that no row lists, gets a row and a blank map (newmap.blank) that the map editor fills. One that
no row lists and that the old model had stays unlisted. Rows of transport maps (SVG marked data-transport) and new
subsystems whose name starts with "Transport" are left to transport_map.py. A removed and an added subsystem that share
at least RENAME_SHARE of the removed one's reactions are reported as a possible rename; the lists are not changed for it.

Usage: maplists.py <old model yml> <new model yml> <subsystemSVG.tsv> <compartmentSVG.tsv> <folder with the SVG maps>
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mapedit  # noqa: E402
import newmap  # noqa: E402

RENAME_SHARE = 0.8
KINDS = ("subsystem", "compartment")


def subsystem_reactions(rxns):
    """{subsystem name: set of reaction ids}."""
    out = {}
    for rid, r in rxns.items():
        subs = r.get("subsystem")
        for s in (subs if isinstance(subs, list) else [subs]):
            if s:
                out.setdefault(str(s), set()).add(rid)
    return out


def read_model(path):
    """(subsystem name -> reaction ids, {compartment name in lower case: the model's spelling}) of a model YAML file."""
    _, rxns, _ = mapedit.load_model(path)
    return subsystem_reactions(rxns), {c.lower(): c for c in mapedit.load_model.compartments.values()}


def read_table(path):
    """(lines, rows): all lines of the table and, for each map row, its line number, name, map name and file."""
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    rows = []
    for i, line in enumerate(lines):
        f = line.split("\t")
        if line.strip() and not line.startswith(("#", "@")) and len(f) >= 3:
            rows.append((i, f[0].strip(), f[1].strip(), f[2].strip()))
    return lines, rows


def is_transport(svg_dir, file):
    path = os.path.join(svg_dir, file)
    return os.path.exists(path) and b"data-transport=" in open(path, "rb").read(3000)


def file_name(title, taken):
    stem = re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_") or "map"
    name, n = f"{stem}.svg", 2
    while name in taken:
        name, n = f"{stem}_{n}.svg", n + 1
    return name


def adapt(old_yml, new_yml, subsystem_table, compartment_table, svg_dir):
    """Change the two tables and the SVGs in svg_dir. Returns {"removed": {kind: [(name, file)]},
    "added": {kind: [(name, file)]}, "renamed": [(old name, new name, share)]}."""
    old_subs, old_comps = read_model(old_yml)
    new_subs, new_comps = read_model(new_yml)
    existing = {"subsystem": {n: n for n in new_subs}, "compartment": new_comps}
    previous = {"subsystem": set(old_subs), "compartment": set(old_comps)}
    tables = {"subsystem": subsystem_table, "compartment": compartment_table}
    report = {"removed": {k: [] for k in KINDS}, "added": {k: [] for k in KINDS}, "renamed": []}
    taken = {f for f in os.listdir(svg_dir) if f.endswith(".svg")}
    palette = 0
    listed_names = {r[1] for r in read_table(tables["subsystem"])[1]}  # before the table is rewritten below
    for kind in KINDS:
        lines, rows = read_table(tables[kind])
        key = (lambda n: n) if kind == "subsystem" else (lambda n: n.lower())
        drop = [r for r in rows if key(r[1]) not in existing[kind] and not is_transport(svg_dir, r[3])]
        for _, name, _, file in drop:
            report["removed"][kind].append((name, file))
            if os.path.exists(os.path.join(svg_dir, file)):
                os.remove(os.path.join(svg_dir, file))
            taken.discard(file)
        kept = [r for r in rows if r not in drop]
        listed = {key(r[1]) for r in rows}
        # compared in lower case for compartments; the row and the map get the model's spelling, which mapedit.py
        # looks up case-sensitively
        new = sorted(existing[kind][n] for n in set(existing[kind]) - previous[kind] - listed
                     if not (kind == "subsystem" and n.lower().startswith("transport")))
        templates = [r[3] for r in kept if not is_transport(svg_dir, r[3])] or \
            [r[3] for o in KINDS for r in read_table(tables[o])[1] if os.path.exists(os.path.join(svg_dir, r[3]))
             and not is_transport(svg_dir, r[3])]
        added = []
        for name in new:
            if not templates:
                raise SystemExit(f"no map to use as the template of the blank map for {name}")
            file = file_name(name, taken)
            taken.add(file)
            newmap.blank(os.path.join(svg_dir, templates[0]), name, os.path.join(svg_dir, file),
                         colour=newmap.PALETTE[(len(kept) + palette) % len(newmap.PALETTE)])
            palette += 1
            added.append((name, file))
            report["added"][kind].append((name, file))
        dropped_lines = {r[0] for r in drop}
        out = [line for i, line in enumerate(lines) if i not in dropped_lines]
        while out and not out[-1].strip():
            out.pop()
        out += [f"{name}\t{name}\t{file}" for name, file in added]
        with open(tables[kind], "w", encoding="utf-8", newline="") as fh:
            fh.write("\n".join(out) + "\n")
    removed = {n for n in old_subs if n not in new_subs and n in listed_names}
    added = {n for n in new_subs if n not in old_subs}
    for gone in sorted(removed):
        for now in sorted(added):
            share = len(old_subs[gone] & new_subs[now]) / len(old_subs[gone])
            if share >= RENAME_SHARE:
                report["renamed"].append((gone, now, share))
    return report


def describe(report):
    """Markdown lines for the pull request."""
    lines = []
    for kind in KINDS:
        for name, file in report["removed"][kind]:
            lines.append(f"- The {kind} `{name}` is not in the model any more: its row and `{file}` are removed.")
        for name, file in report["added"][kind]:
            lines.append(f"- The {kind} `{name}` is new: it has a row and the blank map `{file}`, which the editor fills.")
    for gone, now, share in report["renamed"]:
        lines.append(f"- The subsystem `{gone}` may have been renamed to `{now}` ({share:.0%} of its reactions are in "
                     "the new one). Its map is removed and a blank map is added; rename by hand if it is the same subsystem.")
    return lines


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    print("\n".join(describe(adapt(*sys.argv[1:]))) or "no changes to the lists")
