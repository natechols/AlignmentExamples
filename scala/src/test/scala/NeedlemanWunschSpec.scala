
import org.scalatest._
import Matchers._

import com.github.natechols.alignmentexamples.NeedlemanWunsch

class NeedlemanWunschSpec extends FlatSpec {
  "Needleman-Wunsch alignment" should "Align short sequences" in {
    val nw = new NeedlemanWunsch("HEAGAWGHEE", "PAWHEAE", "BLOSUM50")
    assert(nw.alignment == ("HEAGAWGHE-E", "--P-AW-HEAE"))
    nw.identity shouldBe 0.4 +- 0.00001
  }
}

