
import org.scalatest.FlatSpec

import com.github.natechols.alignmentexamples.Blosum

class BlosumSpec extends FlatSpec {
  "Blosum matrices" should "Load BLOSUM50" in {
      val m = Blosum.getMatrix("BLOSUM50")
      assert(m(('W','W')) == 15)
      assert(m(('A','X')) == -1)
  }
  it should "Load BLOSUM62" in {
      val m = Blosum.getMatrix("BLOSUM62")
      assert(m(('A','A')) == 4)
      assert(m(('T','*')) == -4)
  }
}
