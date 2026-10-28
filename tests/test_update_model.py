"""The parts of the model update that decide whether an update PR is correct."""
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "update"))
import update_model as um  # noqa: E402
import update_pr_body  # noqa: E402

RAVEN3 = """\
- metaData: !!omap
  - id: yeastGEM_v9.1.1
  - name: The Consensus Genome-Scale Metabolic Model of Yeast
  - version: 9.1.1
  - date: 2026-08-28
- metabolites:
  - id: s_0001
"""

OLD = """\
- metaData:
    short_name: "Yeast-GEM"
    full_name: "Consensus model"
    id: "yeastGEM_v9.0.1"
    version: "9.0.1"
    date: "2025-01-01"
    github: "https://github.com/SysBioChalmers/yeast-GEM"
- metabolites:
  - id: s_0001
"""


class MetaData(unittest.TestCase):
    def test_raven3_metadata_is_merged_with_the_fields_only_the_old_copy_had(self):
        head, block, rest = um.split_metadata(RAVEN3)
        _, old_block, _ = um.split_metadata(OLD)
        merged = dict(um.merge_metadata(um.parse_metadata(block), um.parse_metadata(old_block), "Yeast-GEM"))
        self.assertEqual(merged["short_name"], '"Yeast-GEM"')
        self.assertEqual(um.unquote(merged["version"]), "9.1.1")
        self.assertEqual(um.unquote(merged["date"]), "2026-08-28")
        self.assertIn("full_name", merged)
        self.assertIn("github", merged)

    def test_metadata_as_the_last_item_of_the_file(self):
        head, block, rest = um.split_metadata("- metaData:\n    version: 1\n")
        self.assertEqual(rest, [])  # no crash when no item follows
        self.assertEqual(um.parse_metadata(block), [("version", "1")])


class Tables(unittest.TestCase):
    def test_ragged_rows_are_found(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "reactions.tsv")
            with open(path, "w") as fh:
                fh.write("rxns\ta\tb\nr1\t1\t2\nr2\t1\nr3\t1\t2\n")
            self.assertEqual(um.ragged_rows(path), [(3, "r2", 2, 3)])


class Baseline(unittest.TestCase):
    def test_baseline_is_not_the_folder_of_the_new_run(self):
        self.assertEqual(um.baseline_dir("w", "M", "1.0.0", "1.1.0"), os.path.join("w", "M-1.0.0"))
        self.assertNotEqual(um.baseline_dir("w", "M", "1.1.0", "1.1.0"), os.path.join("w", "M-1.1.0"))


class CopyMaps(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        d = self.dir.name
        self.repo = os.path.join(d, "maps-repo")
        self.model_dir = os.path.join(d, "data-files", "integrated-models", "M")
        self.svg_dir = os.path.join(d, "data-files", "svg", "M")
        for kind in ("subsystem", "compartment"):
            os.makedirs(os.path.join(self.repo, kind, "svg"))
        os.makedirs(self.model_dir)
        os.makedirs(self.svg_dir)
        self.write(os.path.join(self.repo, "README.md"), "The maps match **M 2.0.0**.\n")
        self.write(os.path.join(self.repo, "subsystem", "svg", "a.svg"), "<svg>new</svg>")
        self.write(os.path.join(self.repo, "subsystemSVG.tsv"), "@s\tm\tf\nA\tA\ta.svg\n")
        self.write(os.path.join(self.repo, "compartmentSVG.tsv"), "@c\tm\tf\n")
        self.write(os.path.join(self.svg_dir, "a.svg"), "<svg>old</svg>")
        self.write(os.path.join(self.svg_dir, "gone.svg"), "<svg/>")
        self.write(os.path.join(self.svg_dir, "custom.svg"), "<svg/>")
        self.write(os.path.join(self.model_dir, "customSVG.tsv"), "@m\tf\nCustom\tcustom.svg\n")
        self.args = mock.Mock(model="M", maps_repo=self.repo)
        patch = mock.patch.object(um, "DATA_FILES", os.path.join(d, "data-files"))
        patch.start()
        self.addCleanup(patch.stop)
        self.addCleanup(self.dir.cleanup)

    @staticmethod
    def write(path, text):
        with open(path, "w") as fh:
            fh.write(text)

    def test_maps_and_lists_are_copied_and_unlisted_maps_removed_except_custom_ones(self):
        self.assertEqual(um.copy_maps(self.args, self.model_dir, "2.0.0"), [])
        self.assertEqual(sorted(os.listdir(self.svg_dir)), ["a.svg", "custom.svg"])
        self.assertEqual(open(os.path.join(self.svg_dir, "a.svg")).read(), "<svg>new</svg>")
        self.assertIn("A\tA\ta.svg", open(os.path.join(self.model_dir, "subsystemSVG.tsv")).read())

    def test_nothing_is_copied_when_the_maps_repository_holds_another_release(self):
        problems = um.copy_maps(self.args, self.model_dir, "2.1.0")
        self.assertEqual(len(problems), 1)
        self.assertEqual(open(os.path.join(self.svg_dir, "a.svg")).read(), "<svg>old</svg>")
        self.assertTrue(os.path.exists(os.path.join(self.svg_dir, "gone.svg")))

    def test_a_list_row_without_a_map_stops_the_copy(self):
        self.write(os.path.join(self.repo, "subsystemSVG.tsv"), "@s\tm\tf\nA\tA\ta.svg\nB\tB\tb.svg\n")
        self.assertEqual(len(um.copy_maps(self.args, self.model_dir, "2.0.0")), 1)
        self.assertTrue(os.path.exists(os.path.join(self.svg_dir, "gone.svg")))

    def test_no_maps_repository_is_a_problem_unless_the_maps_are_skipped(self):
        self.args.maps_repo = None
        self.assertEqual(len(um.copy_maps(self.args, self.model_dir, "2.0.0")), 1)

    def test_a_file_name_in_both_kinds_is_a_problem(self):
        self.write(os.path.join(self.repo, "compartment", "svg", "a.svg"), "<svg/>")
        self.assertEqual(len(um.copy_maps(self.args, self.model_dir, "2.0.0")), 1)
        self.assertEqual(open(os.path.join(self.svg_dir, "a.svg")).read(), "<svg>old</svg>")

    def test_the_maps_version_is_read_from_the_readme(self):
        self.assertEqual(um.maps_version(self.repo), "2.0.0")


class PullRequestText(unittest.TestCase):
    def test_status_2_covers_maps_that_are_not_up_to_date(self):
        self.assertIn("does not hold the maps of this release", update_pr_body.STATUS[2])


if __name__ == "__main__":
    unittest.main()
