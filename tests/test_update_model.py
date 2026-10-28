"""The parts of the model update that decide whether an update PR is correct."""
import os
import subprocess
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


class MapFailure(unittest.TestCase):
    def update_maps(self, completed):
        with tempfile.TemporaryDirectory() as d:
            model_dir = os.path.join(d, "data-files", "integrated-models", "M")
            svg_dir = os.path.join(d, "data-files", "svg", "M")
            os.makedirs(model_dir)
            os.makedirs(svg_dir)
            with open(os.path.join(svg_dir, "a.svg"), "w") as fh:
                fh.write("<svg/>")
            args = mock.Mock(model="M", gene_label=None, kegg_dir=None)
            summary = os.path.join(d, "work", "maps_summary.md")
            with mock.patch.object(um, "DATA_FILES", os.path.join(d, "data-files")), \
                    mock.patch.object(um.subprocess, "run", return_value=completed):
                problems = um.update_maps(args, model_dir, "new.yml", "old.yml", os.path.join(d, "work"), summary)
            return problems, open(summary).read()

    def test_a_failing_map_editor_is_reported(self):
        failed = subprocess.CompletedProcess([], 1, stdout="", stderr="Traceback: boom")
        problems, summary = self.update_maps(failed)
        self.assertEqual(len(problems), 1)
        self.assertTrue(summary.startswith("**Maps NOT updated"))
        self.assertIn("boom", summary)


class PullRequestText(unittest.TestCase):
    def test_status_2_covers_failed_maps(self):
        self.assertIn("maps could not be updated", update_pr_body.STATUS[2])


if __name__ == "__main__":
    unittest.main()
