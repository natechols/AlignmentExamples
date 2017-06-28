Awk implementation of Needleman-Wunsch alignment
------------------------------------------------

This is a somewhat ridiculous example, but it's a completely functional
program with output identical to the other implementations:

  $ nwalign.awk seq1.fasta seq2.fasta

however it is restricted to BLOSUM50 with a gap penalty of 8.  I have only
tested it with GNU Awk 4.1.3; vanilla awk (as found on MacOS for example) may
or may not work.  Running "make" will generate a symlink called "nwalign" for
parity with other builds.

Note that at present the BLOSUM50 matrix is embedded in the source code, so
the entire mess is LGPL.
