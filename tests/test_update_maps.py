"""The checks of update_maps.py that run before a maps repository is changed."""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "maps"))
import update_maps  # noqa: E402

SCRIPT = os.path.join(HERE, "..", "maps", "update_maps.py")


class Readme(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def readme(self, text):
        with open(os.path.join(self.dir.name, "README.md"), "w") as fh:
            fh.write(text)

    def test_the_release_is_read_for_the_right_model(self):
        self.readme("The maps match **Yeast-GEM 9.1.1**.\n")
        self.assertEqual(update_maps.readme_version(self.dir.name, "Yeast-GEM"), "9.1.1")

    def test_a_readme_of_another_model_stops_the_update(self):
        self.readme("The maps match **Human-GEM 2.1.1**.\n")
        with self.assertRaises(SystemExit) as stop:
            update_maps.readme_version(self.dir.name, "Yeast-GEM")
        self.assertEqual(stop.exception.code, 1)

    def test_a_readme_without_the_line_stops_the_update(self):
        self.readme("No version here.\n")
        with self.assertRaises(SystemExit) as stop:
            update_maps.readme_version(self.dir.name, "Yeast-GEM")
        self.assertEqual(stop.exception.code, 1)

    def test_a_work_folder_with_files_is_refused_and_kept(self):
        self.readme("The maps match **Human-GEM 2.1.1**.\n")
        work = os.path.join(self.dir.name, "work")
        os.makedirs(work)
        keep = os.path.join(work, "keep.txt")
        open(keep, "w").close()
        result = subprocess.run([sys.executable, SCRIPT, "--model", "Human-GEM", "--repo", self.dir.name,
                                 "--version", "2.1.1", "--force", "--work-dir", work], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertTrue(os.path.exists(keep))
        self.assertIn("not empty", result.stderr)


if __name__ == "__main__":
    unittest.main()
