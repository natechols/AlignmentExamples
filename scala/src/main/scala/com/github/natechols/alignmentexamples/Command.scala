
package com.github.natechols.alignmentexamples

import java.nio.file.{Path,Paths}

import scala.math._

import scopt.OptionParser


trait NWAlignRunner {
  case class NWAlignOptions(seq1: Path,
                            seq2: Path,
                            scoreMatrixName: String = "BLOSUM50",
                            gapPenalty: Int = 8)
  lazy val DEFAULTS = NWAlignOptions(null, null)

  val parser = new OptionParser[NWAlignOptions]("nwalign") {
    head("Needleman-Wunsch alignment (Scala version)")

    arg[String]("seq1") action { (s, c) =>
      c.copy(seq1 = Paths.get(s))
    } text "First FASTA file to align (only the first sequence will be used)"
    arg[String]("seq2") action { (s, c) =>
      c.copy(seq2 = Paths.get(s))
    } text "Second FASTA file to align (only the first sequence will be used)"
    opt[String]('s', "score-matrix-name") action { (n, c) =>
      c.copy(scoreMatrixName = n)
    } text "Scoring matrix to use (default: BLOSUM50)"
    opt[Int]('g', "gap-penalty") action { (g, c) =>
      c.copy(gapPenalty = abs(g))
    } text "Penalty for opening a gap (will be treated as a negative number"
  }

  def runner(args: Array[String]): Unit = {

    val result = parser.parse(args, DEFAULTS) map { c =>
      val rec1 = FastaReader(c.seq1)
      val rec2 = FastaReader(c.seq2)
      val nw = new NeedlemanWunsch(rec1(0), rec2(0), c.scoreMatrixName,
                                   c.gapPenalty)
      println("Sequence identity = %.2f%%".format(nw.identity*100))
      println("")
      println(nw.format())
      0
    }

    val exitCode = result.getOrElse(1)
    System.exit(exitCode)
  }
}

object Program extends App with NWAlignRunner{
  runner(args)
}
