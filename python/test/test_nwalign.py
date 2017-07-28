
import tempfile
import unittest
import os.path as op

from nwalign import needlemanWunsch, read_fasta


class TestFastaReader(unittest.TestCase):

    def test_1(self):
        f = tempfile.NamedTemporaryFile(suffix=".fasta").name
        fasta_in = ">seq1\nHEAGAWGHEE\n>seq2\nPAWHEAE"
        with open(f, "w") as out:
            out.write(fasta_in)
        records = read_fasta(f)
        self.assertEqual(len(records), 2)
        self.assertEqual("\n".join([str(r) for r in records]), fasta_in)

    def test_2(self):
        fn = op.join(op.dirname(__file__), "../../data/hemoglobinA.fasta")
        records = read_fasta(fn)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "hemoglobin_A")


class TestNeedlemanWunsch(unittest.TestCase):

    def test_1(self): # Durbin et al. p. 21
        nw = needlemanWunsch("HEAGAWGHEE", "PAWHEAE")
        self.assertEqual(nw.alignment(), ('HEAGAWGHE-E', '--P-AW-HEAE'))
        self.assertEqual(nw.length(), 11)
        self.assertAlmostEqual(nw.identity(), 0.4545, 4)

    def test_2(self):
        seq1 = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
        seq2 = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        nw = needlemanWunsch(seq1, seq2)
        self.assertEqual(nw.alignment(), ('V-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR', 'VHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'))
        self.assertAlmostEqual(nw.identity(), 0.4324, 4)
