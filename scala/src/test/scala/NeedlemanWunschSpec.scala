
import org.scalatest._
import Matchers._

import com.github.natechols.alignmentexamples.NeedlemanWunsch

class NeedlemanWunschSpec extends FlatSpec {
  "Needleman-Wunsch alignment" should "Align short sequences" in {
    val nw = NeedlemanWunsch("HEAGAWGHEE", "PAWHEAE", "BLOSUM50")
    assert(nw.size == 11)
    assert(nw.format() == "seq1  HEAGAWGHE-E\n          ** ** *\nseq2  --P-AW-HEAE")
    assert(nw.alignment == ("HEAGAWGHE-E", "--P-AW-HEAE"))
    nw.identity shouldBe 0.4 +- 0.00001
  }
}

