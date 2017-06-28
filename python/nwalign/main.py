
"""
Command-line utility for Needleman-Wunsch global alignment of two protein
sequences.

Usage: nwalign seq1.fasta seq2.fasta
"""

from __future__ import print_function
import argparse
import os.path as op
import sys

from .fasta import read_fasta
from .core import NeedlemanWunsch


VERSION = "0.1"


def validate_file(fn):
    if op.isfile(fn):
        return op.realpath(fn)
    else:
        raise IOError("Argument {a} is not a file".format(a=fn))


def _get_parser():
    p = argparse.ArgumentParser(description=__doc__, version=VERSION)
    p.add_argument("seq1", type=validate_file)
    p.add_argument("seq2", type=validate_file)
    p.add_argument("-m", "--matrix", action="store", default="BLOSUM50",
                   help="Name of BLOSUM matrix to use")
    p.add_argument("-g", "--gap-penalty", action="store", type=int,
                   default=8, help="Penalty for opening a gap")
    return p


def main(argv=sys.argv[1:], out=sys.stdout):
    p = _get_parser()
    args = p.parse_args(argv)
    records1 = read_fasta(args.seq1)
    records2 = read_fasta(args.seq2)
    if len(records1) > 1:
        warnings.warn("Multiple records in file 1, will only align the 1st")
    if len(records2) > 1:
        warnings.warn("Multiple records in file 2, will only align the 1st")
    nw = NeedlemanWunsch(
        records1[0], records2[0], args.matrix, args.gap_penalty)
    print("Sequence identity = {i:.2f}%".format(i=nw.identity() * 100),
          file=out)
    out.write("\n")
    nw.show(out=out)
    return 0


if __name__ == "__main__":
    main(sys.argv[1:])
