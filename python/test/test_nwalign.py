
import unittest

from nwalign import NeedlemanWunsch


class TestNeedlemanWunsch(unittest.TestCase):

    def test_1(self): # Durbin et al. p. 21
        nw = NeedlemanWunsch("HEAGAWGHEE", "PAWHEAE")
        self.assertEqual(nw.alignment(), ('HEAGAWGHE-E', '--P-AW-HEAE'))
        self.assertEqual(nw.identity(), 0.4)
