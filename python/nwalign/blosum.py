# Author: Nat Echols
# This source code is in the public domain.

import os.path


def _load_file(file_name):
    path = os.path.join(os.path.dirname(__file__), "resources", file_name)
    return open(path).read()


class BlosumMatrix(object):
    def __init__(self, name, data):
        self._scores = {}
        aa_cols = None
        for line in data.splitlines():
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            if line.startswith(" ") and aa_cols is None:
                aa_cols = fields
                assert len(aa_cols) == 24 # XXX ???
                for aa in aa_cols:
                    self._scores[aa] = {}
            else:
                assert len(fields) == len(aa_cols) + 1
                aa = fields[0]
                for i, score in enumerate(fields[1:]):
                    self._scores[aa][aa_cols[i]] = int(score)

    def __call__(self, a, b):
        return self._scores[a.upper()][b.upper()]


BLOSUM50 = BlosumMatrix("BLOSUM50", _load_file("BLOSUM50.txt"))
BLOSUM62 = BlosumMatrix("BLOSUM62", _load_file("BLOSUM62.txt"))
MATRICES = {
    "BLOSUM50": BLOSUM50,
    "BLOSUM62": BLOSUM62,
}
