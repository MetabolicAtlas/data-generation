"""A reaction is reversible when its lower bound is below 0 and its upper bound above 0."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "maps"))
import mapedit  # noqa: E402


class ReversibleTest(unittest.TestCase):
    def test_bounds(self):
        self.assertTrue(mapedit.reversible({"lower_bound": -1000, "upper_bound": 1000}))
        self.assertFalse(mapedit.reversible({"lower_bound": 0, "upper_bound": 1000}))
        self.assertFalse(mapedit.reversible({"lower_bound": -1000, "upper_bound": 0}))
        self.assertFalse(mapedit.reversible({"lower_bound": -1000}))
        self.assertFalse(mapedit.reversible({}))


if __name__ == "__main__":
    unittest.main()
