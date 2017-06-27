
import org.scalatest.FlatSpec

import java.io.{File,FileWriter,BufferedWriter}

import com.github.natechols.alignmentexamples.FastaReader


class FastaSpec extends FlatSpec {
  "FASTA reader" should "load simple FASTA file" in {
    val FASTA = ">seq1 this is a protein\nHEAGAWGHEE\nYQCND\n>seq2\nPAWHEAE"
    val tmpFile = File.createTempFile("test", ".fasta")
    val bw = new BufferedWriter(new FileWriter(tmpFile))
    bw.write(FASTA)
    bw.close()
    val records = FastaReader(tmpFile.toPath)
    assert(records.size == 2)
    assert(records(0).id == "seq1")
    assert(records(0).sequence == "HEAGAWGHEEYQCND")
    assert(records(1).format() == ">seq2\nPAWHEAE")
  }
}
