# Author: Nat Echols
# This source code is in the public domain.

"""
Needleman-Wunsch global pairwise (protein) sequence alignment, as described in
Durbin et al. chapter 2.
"""

from collections import namedtuple
import sys

from .blosum import MATRICES

MatrixCell = namedtuple("MatrixCell", ("score", "pointer"))
UP = (0, -1)
LEFT = (-1, 0)
DIAG = (-1, -1)


class NeedlemanWunsch(object):

    """
    Result of a global pairwise sequence alignment.
    """

    def __init__(self, seq1, seq2, score_matrix="BLOSUM50", gap_penalty=8):
        self._seq1 = seq1
        self._seq2 = seq2
        self._score_matrix = score_matrix
        if isinstance(score_matrix, basestring):
            score_matrix = MATRICES[score_matrix]
        max_x = len(seq1) + 1
        max_y = len(seq2) + 1
        f = [[None for y in range(max_y)] for x in range(max_x)]
        # fill top row
        for i in range(max_x):
            f[i][0] = MatrixCell(-i * gap_penalty, LEFT)
        # fill left column
        for j in range(max_y):
            f[0][j] = MatrixCell(-j * gap_penalty, UP)
        for i in range(1, max_x):
            aa = seq1[i - 1]
            for j in range(1, max_y):
                bb = seq2[j - 1]
                f[i][j] = max(
                    MatrixCell(
                        f[i - 1][j - 1].score + score_matrix(aa, bb), DIAG),
                    MatrixCell(f[i - 1][j].score - gap_penalty, LEFT),
                    MatrixCell(f[i][j - 1].score - gap_penalty, UP),
                    key=lambda c: c.score)
        self._dp_matrix = f
        tmp_seq1 = []
        tmp_seq2 = []
        i = len(seq1)
        j = len(seq2)
        # traceback step
        while i >= 0 and j >= 0:
            cell = f[i][j]
            if cell.pointer == DIAG:
                tmp_seq1.append(seq1[i - 1])
                tmp_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif cell.pointer == LEFT:
                tmp_seq1.append(seq1[i - 1])
                tmp_seq2.append("-")
                i -= 1
            elif cell.pointer == UP:
                tmp_seq1.append("-")
                tmp_seq2.append(seq2[j - 1])
                j -= 1
            else:
                assert False
            if i == j == 0:
                break
        self._seq1_aligned = "".join(tmp_seq1[::-1])
        self._seq2_aligned = "".join(tmp_seq2[::-1])

    def alignment(self):
        return (self._seq1_aligned, self._seq2_aligned)

    def length(self):
        return len(self._seq1_aligned)

    def identity(self):
        n_mm = n_ins = n_del = 0
        for aa, bb in zip(self._seq1_aligned, self._seq2_aligned):
            if aa == "-":
                n_ins += 1
            elif bb == "-":
                n_del += 1
            elif aa != bb:
                n_mm += 1
        return max(0, 1. - float(n_mm + n_ins + n_del) / self.length())

    def show(self, out=sys.stdout, seq_width=50):
        def _get_seq_matches(s1, s2):
            chars = []
            for aa, bb in zip(s1, s2):
                if aa == bb:
                    chars.append('*')
                else:
                    chars.append(" ")
            return "".join(chars)
        k = 0
        rows = []
        id1, id2 = "seq1", "seq2"
        if not isinstance(self._seq1, basestring):
            id1 = self._seq1.id
        if not isinstance(self._seq2, basestring):
            id2 = self._seq2.id
        id_width = max(len(id1), len(id2))
        line_fmt = "%%%ds  %%s" % id_width
        while k < self.length():
            s1 = self._seq1_aligned[k:k+seq_width]
            s2 = self._seq2_aligned[k:k+seq_width]
            seq_matches = _get_seq_matches(s1, s2)
            line1 = line_fmt % (id1, s1)
            line2 = line_fmt % ("", seq_matches)
            line3 = line_fmt % (id2, s2)
            k += seq_width
            rows.append("\n".join([line1,line2,line3]))
        out.write("\n\n\n".join(rows))
        out.write("\n")
