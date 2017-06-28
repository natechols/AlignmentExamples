
package libnwalign

import (
  "testing"
  "math"
  "io/ioutil"
)

func TestNeedlemanWunsch (t *testing.T) {
  var a Alignment
  a = NeedlemanWunsch("HEAGAWGHEE", "PAWHEAE", "BLOSUM50", 8)
  if a.seq1Aligned != "HEAGAWGHE-E" {
    t.Error("Expected HEAGAWGHE-E, got ", a.seq1Aligned)
  } else if a.seq2Aligned != "--P-AW-HEAE" {
    t.Error("Expected --P-AW-HEAE, got ", a.seq2Aligned)
  }
  if a.Length() != 11 {
    t.Error("Expected length == 11, got ", a.Length())
  }
  if math.Abs(a.Identity() - 0.4545) > 0.0001 {
    t.Error("Expected identity ~ 0.4, got ", a.Identity())
  }
  //a.show("seq1", "seq2")
}

func TestFastaReader (t *testing.T) {
  f, err := ioutil.TempFile("/tmp", "test-fasta")
  if (err != nil) {
    t.Error("Error creating temp file: " + err.Error())
  }
  f.Close()
  fastaRaw := ">seq1\nHEAGAWGHEE\nQWERTY\n>seq2\nPAWHEAE"
  ioutil.WriteFile(f.Name(), []byte(fastaRaw), 0644)
  records := FastaReader(f.Name())
  if (len(records) != 2) {
    t.Error("Expected two fasta records")
  }
  if (records[0].Id != "seq1") {
    t.Error("Wrong sequence ID: " + records[0].Id)
  }
  if (records[0].Sequence != "HEAGAWGHEEQWERTY") {
    t.Error("Wrong sequence: " + records[0].Sequence)
  }
  if (records[1].Id != "seq2") {
    t.Error("Wrong sequence ID: " + records[1].Id)
  }
  if (records[1].Sequence != "PAWHEAE") {
    t.Error("Wrong sequence: " + records[1].Sequence)
  }
}
