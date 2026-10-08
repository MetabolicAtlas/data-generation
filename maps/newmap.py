"""Blank maps in the style of an existing map of the same model, for subsystems or compartments without one.

The template's layers, background band group, title and licence badge are kept; its content is removed.
The editor then fills the blank map (D8 where a KEGG map covers the subsystem, D9 for the rest).
"""
import os
import sys

from lxml import etree

SVG = "{http://www.w3.org/2000/svg}"
WIDTH, HEIGHT = 4000.0, 2500.0
# band colours of the Human-GEM maps, given to new maps in turn
PALETTE = ["#64ff64", "#9ef0ce", "#ffce64", "#ffb199", "#64c4ff", "#f37bff", "#c6bfcf", "#e6b3d7", "#b6bce6",
           "#f8da75", "#08d6b5", "#8864ff", "#fa9005", "#01feee", "#ff64f0", "#64a2ff", "#b8e6c5", "#e0e6a3",
           "#cbb8e6", "#c4dde6", "#e6afb3", "#d5e6ba", "#319131", "#f6b40a"]


def blank(template, title, out, width=WIDTH, height=HEIGHT, colour=None):
    tree = etree.parse(template)
    root = tree.getroot()
    w0, h0 = float(root.get("width")), float(root.get("height"))
    root.set("width", f"{width:g}")
    root.set("height", f"{height:g}")
    root.set("data-new", "1")  # the editor centres the title once the map is filled
    for g in root.iter(SVG + "g"):
        if g.get("id") in ("fes", "fehs", "ees", "nodes"):
            for c in list(g):
                g.remove(c)
    bands = [g for g in root.iter(SVG + "g") if g.get("class") == "subsystem"]
    for g in bands[1:]:
        g.getparent().remove(g)
    for g in [g for g in root.iter(SVG + "g") if g.get("class") == "compartment"]:
        g.getparent().remove(g)
    if bands:
        band = bands[0]
        band.set("id", title)
        band.attrib.pop("data-bounds", None)
        path = band.find(SVG + "path")
        if path is not None:
            path.set("d", "")
            if colour:
                path.set("stroke", colour)
        t = band.find(SVG + "text")
        if t is not None:
            t.text = title
            if colour:
                t.set("fill", colour)
            size = float(t.get("font-size") or 100)
            x = width / 2 if t.get("text-anchor") == "middle" else width / 2 - 0.28 * size * len(title)
            t.set("transform", f"matrix(1,0,0,1,{x:.1f},{size * 2:.1f})")
    for g in root.iter(SVG + "g"):
        if g.get("id") == "license":
            g.set("transform", f"translate({width - w0:g},{height - h0:g})")
    tree.write(out, xml_declaration=False, encoding="utf-8")


if __name__ == "__main__":
    # usage: newmap.py <template.svg> <out dir> <table: filename<TAB>title per line>
    template, outdir, table = sys.argv[1:4]
    os.makedirs(outdir, exist_ok=True)
    for line in open(table):
        if line.strip():
            fname, title = line.rstrip("\n").split("\t")[:2]
            blank(template, title, os.path.join(outdir, fname))
            print("blank map", fname, title)
