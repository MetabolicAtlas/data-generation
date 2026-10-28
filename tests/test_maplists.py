"""Adapting the lists of maps to the subsystems and compartments of a new release."""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "maps"))
import maplists  # noqa: E402

SVG = ('<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100"><g id="fes"/>'
       '<g class="subsystem" id="t"><path d="M0 0"/><text font-size="10">t</text></g></svg>')
TRANSPORT = SVG.replace("<g id=\"fes\"/>", '<g id="fes" data-transport="1"/>')


def model(subsystems, compartments):
    """A model file with one reaction per entry of `subsystems`: {subsystem: [reaction ids]}."""
    lines = ["- metaData:", "    version: \"1.0.0\"", "- compartments:"]
    lines += [f"    {c[0]}: {c}" for c in compartments]
    lines += ["- metabolites:", "  - id: m1", "    name: x", f"    compartment: {compartments[0][0]}", "- reactions:"]
    for sub, ids in subsystems.items():
        for rid in ids:
            lines += [f"  - id: {rid}", f"    subsystem: {sub}", "    metabolites: {m1: -1}"]
    return "\n".join(lines) + "\n"


class Adapt(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.svg = os.path.join(self.dir.name, "svg")
        os.makedirs(self.svg)
        for name, text in (("a.svg", SVG), ("b.svg", SVG), ("cyto.svg", SVG), ("transport_x.svg", TRANSPORT)):
            self.put(os.path.join(self.svg, name), text)
        self.subsystems = self.put(os.path.join(self.dir.name, "s.tsv"),
                                   "# comment\n\n@subsystem name\tmap name\tfilename\nA\tA\ta.svg\nB\tB\tb.svg\n"
                                   "Transport reactions\tTransport: x\ttransport_x.svg\n")
        self.compartments = self.put(os.path.join(self.dir.name, "c.tsv"),
                                     "@compartment name\tMap name\tfilename\nCytosol\tCytosol\tcyto.svg\n")

    def put(self, path, text):
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def adapt(self, old, new):
        old_yml = self.put(os.path.join(self.dir.name, "old.yml"), old)
        new_yml = self.put(os.path.join(self.dir.name, "new.yml"), new)
        return maplists.adapt(old_yml, new_yml, self.subsystems, self.compartments, self.svg)

    def test_the_same_release_changes_nothing(self):
        same = model({"A": ["r1"], "B": ["r2"], "Transport reactions": ["r3"]}, ["Cytosol"])
        before = open(self.subsystems).read()
        report = self.adapt(same, same)
        self.assertEqual(open(self.subsystems).read(), before)
        self.assertEqual(report["removed"], {"subsystem": [], "compartment": []})
        self.assertEqual(report["added"], {"subsystem": [], "compartment": []})

    def test_a_removed_subsystem_loses_its_row_and_map_but_a_transport_row_stays(self):
        old = model({"A": ["r1"], "B": ["r2"]}, ["Cytosol"])
        new = model({"A": ["r1"]}, ["Cytosol"])
        report = self.adapt(old, new)
        self.assertEqual(report["removed"]["subsystem"], [("B", "b.svg")])
        self.assertFalse(os.path.exists(os.path.join(self.svg, "b.svg")))
        self.assertNotIn("b.svg", open(self.subsystems).read())
        self.assertIn("transport_x.svg", open(self.subsystems).read())

    def test_a_new_subsystem_gets_a_row_and_a_blank_map(self):
        old = model({"A": ["r1"], "B": ["r2"]}, ["Cytosol"])
        new = model({"A": ["r1"], "B": ["r2"], "New one": ["r4"], "Transport new": ["r5"]}, ["Cytosol"])
        report = self.adapt(old, new)
        self.assertEqual(report["added"]["subsystem"], [("New one", "new_one.svg")])
        self.assertTrue(os.path.exists(os.path.join(self.svg, "new_one.svg")))
        self.assertIn("New one\tNew one\tnew_one.svg\n", open(self.subsystems).read())

    def test_a_subsystem_that_no_row_lists_and_the_old_model_had_stays_unlisted(self):
        same = model({"A": ["r1"], "B": ["r2"], "Unmapped": ["r4"]}, ["Cytosol"])
        report = self.adapt(same, same)
        self.assertEqual(report["added"]["subsystem"], [])

    def test_a_new_compartment_gets_a_row_and_a_removed_one_loses_it(self):
        old = model({"A": ["r1"]}, ["Cytosol"])
        new = model({"A": ["r1"]}, ["Cytosol", "Golgi Membrane"])
        self.assertEqual(self.adapt(old, new)["added"]["compartment"], [("Golgi Membrane", "golgi_membrane.svg")])
        self.assertIn("Golgi Membrane\tGolgi Membrane\tgolgi_membrane.svg\n", open(self.compartments).read())
        report = self.adapt(new, model({"A": ["r1"]}, ["Golgi Membrane"]))
        self.assertEqual(report["removed"]["compartment"], [("Cytosol", "cyto.svg")])

    def test_a_probable_rename_is_reported_and_not_applied(self):
        old = model({"A": ["r1", "r2", "r3", "r4", "r5"], "B": ["r9"]}, ["Cytosol"])
        new = model({"A2": ["r1", "r2", "r3", "r4", "r5"], "B": ["r9"]}, ["Cytosol"])
        report = self.adapt(old, new)
        self.assertEqual([(a, b) for a, b, _ in report["renamed"]], [("A", "A2")])
        self.assertEqual(report["removed"]["subsystem"], [("A", "a.svg")])

    def test_a_gone_subsystem_without_a_row_is_not_reported_as_renamed(self):
        old = model({"A": ["r1", "r2", "r3", "r4", "r5"], "Unlisted": ["r6", "r7", "r8", "r9", "r10"]}, ["Cytosol"])
        new = model({"A": ["r1", "r2", "r3", "r4", "r5"], "Unlisted 2": ["r6", "r7", "r8", "r9", "r10"]}, ["Cytosol"])
        self.assertEqual(self.adapt(old, new)["renamed"], [])

    def test_a_subsystem_list_on_a_reaction_counts_for_each_subsystem(self):
        old = model({"A": ["r1"], "B": ["r2"]}, ["Cytosol"])
        new = old.replace("subsystem: B", "subsystem: [B, C]")
        report = self.adapt(old, new)
        self.assertEqual(report["added"]["subsystem"], [("C", "c.svg")])
        self.assertEqual(report["removed"]["subsystem"], [])

    def test_description_has_one_line_per_change(self):
        old = model({"A": ["r1"], "B": ["r2"]}, ["Cytosol"])
        new = model({"A": ["r1"], "C": ["r3"]}, ["Cytosol"])
        self.assertEqual(len(maplists.describe(self.adapt(old, new))), 2)


if __name__ == "__main__":
    unittest.main()
