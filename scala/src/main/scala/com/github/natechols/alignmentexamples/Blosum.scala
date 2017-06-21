
package com.github.natechols.alignmentexamples

import scala.io.Source

object Blosum {

  def getMatrix(scoreMatrixName: String): Map[(Char,Char), Int] = {
    val lines = Source.fromResource(s"${scoreMatrixName}.txt").getLines.toList
    val aaCols = lines.filter(_.startsWith(" ")).toList(0).trim.split(" +")
    lines.filter(l => ! (l.startsWith("#") || l.startsWith(" "))).map { line =>
      val fields = line.split(" +")
      val aa: Char = fields(0).charAt(0)
      fields.drop(1).zipWithIndex.map { case (score, i) =>
        ((aa, aaCols(i).charAt(0)), score.toInt)
      }
    }.flatten.toMap
  }

}
