
package nwalign

import "testing"
import "math"

func TestNeedlemanWunsch (t *testing.T) {
  var a Alignment
  a = NeedlemanWunsch("HEAGAWGHEE", "PAWHEAE", "BLOSUM50", 8)
  if a.seq1Aligned != "HEAGAWGHE-E" {
    t.Error("Expected HEAGAWGHE-E, got ", a.seq1Aligned)
  } else if a.seq2Aligned != "--P-AW-HEAE" {
    t.Error("Expected --P-AW-HEAE, got ", a.seq2Aligned)
  }
  if a.length() != 11 {
    t.Error("Expected length == 11, got ", a.length())
  }
  if math.Abs(a.identity() - 0.4545) > 0.0001 {
    t.Error("Expected identity ~ 0.4, got ", a.identity())
  }
  //a.show("seq1", "seq2")
}
