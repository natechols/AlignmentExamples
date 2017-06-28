
package com.github.natechols.alignmentexamples

import scala.collection.mutable.ArrayBuffer
import scala.math._


class NeedlemanWunsch(val seq1fasta: FastaRecord,
                      val seq2fasta: FastaRecord,
                      val scoreMatrixName: String = "BLOSUM50",
                      val gapPenalty: Int = 8) {

  case class MatrixCell(score: Int, pointer: (Int, Int))
  object Constants {
    val DIAG = (-1, -1)
    val UP = (0, -1)
    val LEFT = (-1, 0)
  }

  private val scoreMatrix = Blosum.getMatrix(scoreMatrixName)
  private val seq1 = seq1fasta.sequence
  private val seq2 = seq2fasta.sequence
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

  def size: Int = seq1Aligned.size

  def identity: Double = {
    var (nMm, nIns, nDel) = (0,0,0)
    seq1Aligned.toCharArray.zip(seq2Aligned.toCharArray).foreach { case (aa,bb) =>
      if (aa == '-') nIns += 1
      else if (bb == '-') nDel += 1
      else if (aa != bb) nMm += 1
    }
    max(0.0, 1.0 - (nMm+nIns+nDel).toDouble / size.toDouble)
  }

  private def matchChars = {
    seq1Aligned.toCharArray.zip(seq2Aligned.toCharArray).map { case (aa, bb) =>
      if (aa == bb) '*' else ' '
    }.mkString
  }

  def format(seqWidth: Int = 50): String = {
    var k = 0
    val idWidth = max(seq1fasta.id.size.toDouble, seq2fasta.id.size.toDouble).toInt
    val fmtString = "%%%ds  %%s".format(idWidth)
    val lines = new ArrayBuffer[String]()
    var matches = matchChars
    while (k < this.size) {
      if (k > 0) lines += "\n"
      lines += fmtString.format(seq1fasta.id, seq1Aligned.slice(k, k+seqWidth))
      lines += fmtString.format("", matches.slice(k, k+seqWidth))
      lines += fmtString.format(seq2fasta.id, seq2Aligned.slice(k, k+seqWidth))
      k += seqWidth
    }
    lines.mkString("\n")
  }

  def show(seqWidth: Int = 50): Unit = {
    println(format(seqWidth))
  }
}

object NeedlemanWunsch {
  def apply(seq1: String,
            seq2: String,
            scoreMatrixName: String = "BLOSUM50",
            gapPenalty: Int = 8) = {
    val seq1fasta = FastaRecord("seq1", seq1)
    val seq2fasta = FastaRecord("seq2", seq2)
    new NeedlemanWunsch(seq1fasta, seq2fasta, scoreMatrixName, gapPenalty)
  }
}
