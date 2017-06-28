// Author: Nat Echols
// License: public domain

package main

import (
  "libnwalign"
  "fmt"
  "flag"
)

func main() {
  matrixName := flag.String("matrix-name", "BLOSUM50", "Scoring matrix name")
  gapPenalty := flag.Int("gap-penalty", 8, "Penalty for opening a gap")
  flag.Parse()
  args := flag.Args()
  if (len(args) != 2) {
    panic("Exactly two FASTA input files are required")
  }
  fasta1 := libnwalign.FastaReader(args[0])
  fasta2 := libnwalign.FastaReader(args[1])
  if (len(fasta1) > 1) {
    fmt.Printf("WARNING: %d sequences in file 1", len(fasta1))
  }
  if (len(fasta2) > 1) {
    fmt.Printf("WARNING: %d sequences in file 2", len(fasta2))
  }
  nw := libnwalign.NeedlemanWunsch(fasta1[0].Sequence, fasta2[0].Sequence,
                                   *matrixName, *gapPenalty)
  fmt.Printf("Sequence identity = %.2f%%\n\n", nw.Identity()*100)
  nw.Show(fasta1[0].Id, fasta2[0].Id)
}
