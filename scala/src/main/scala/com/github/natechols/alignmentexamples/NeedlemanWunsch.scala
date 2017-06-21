
package com.github.natechols.alignmentexamples

import scala.collection.mutable.ArrayBuffer
import scala.math._


object Constants {
  val DIAG = (-1, -1)
  val UP = (0, -1)
  val LEFT = (-1, 0)
}

class NeedlemanWunsch(val seq1: String,
                      val seq2: String,
                      val scoreMatrixName: String = "BLOSUM50",
                      val gapPenalty: Int = 8) {

  case class MatrixCell(score: Int, pointer: (Int, Int))

  private val scoreMatrix = Blosum.getMatrix(scoreMatrixName)
  private val maxX = seq1.size
  private val maxY = seq2.size

  private val f: Seq[ArrayBuffer[MatrixCell]] = (0 to maxX).map { i =>
    val row = new ArrayBuffer[MatrixCell]()
    row.append(MatrixCell(-i * gapPenalty, Constants.LEFT))
    row
  }
  (1 to maxY).foreach { j =>
    f(0).append(MatrixCell(-j * gapPenalty, Constants.UP))
  }
  (1 to maxX).foreach { i =>
    val aa: Char = seq1.charAt(i - 1)
    (1 to maxY).foreach { j =>
      val bb: Char = seq2.charAt(j - 1)
      val best = Seq(
        MatrixCell(f(i-1)(j-1).score + scoreMatrix((aa, bb)), Constants.DIAG),
        MatrixCell(f(i-1)(j).score - gapPenalty, Constants.LEFT),
        MatrixCell(f(i)(j-1).score - gapPenalty, Constants.UP)
      ).reduceLeft((x,y) => if (x.score >= y.score) x else y)
      f(i).append(best)
    }
  }
  private val tmpSeq1 = new ArrayBuffer[Char]()
  private val tmpSeq2 = new ArrayBuffer[Char]()
  private var (i, j) = (maxX, maxY)
  // FIXME this could be more elegant...
  while ((i >= 0) && (j >= 0) && !((i == j) && (i == 0 ))) {
    val cell = f(i)(j)
    if (cell.pointer == Constants.DIAG) {
      tmpSeq1.append(seq1.charAt(i-1))
      tmpSeq2.append(seq2.charAt(j-1))
      i -= 1
      j -= 1
    } else if (cell.pointer == Constants.LEFT) {
      tmpSeq1.append(seq1.charAt(i-1))
      tmpSeq2.append('-')
      i -= 1
    } else if (cell.pointer == Constants.UP) {
      tmpSeq1.append('-')
      tmpSeq2.append(seq2.charAt(j-1))
      j -= 1
    } else throw new Exception(s"Can't traceback with ${cell.pointer}")
  }
  val seq1Aligned = tmpSeq1.reverse.mkString
  val seq2Aligned = tmpSeq2.reverse.mkString

  def alignment: (String, String) = (seq1Aligned, seq2Aligned)

  def identity: Double = {
    var (nMm, nIns, nDel) = (0,0,0)
    for ((aa, bb) <- seq1Aligned.toCharArray.zip(seq2Aligned.toCharArray)) {
      if (aa == '-') nIns += 1
      else if (bb == '-') nDel += 1
      else if (aa != bb) nMm += 1
    }
    max(0.0, 1.0 - (nMm+nIns+nDel).toDouble / seq1.size.toDouble)
  }
}
