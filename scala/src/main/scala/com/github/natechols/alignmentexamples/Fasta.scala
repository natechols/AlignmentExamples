// Author: Nat Echols
// License: public domain

package com.github.natechols.alignmentexamples

import java.nio.file.Path

import scala.io.Source
import scala.collection.mutable.ListBuffer


case class FastaRecord(id: String,
                       sequence: String,
                       description: Option[String]) {

  def header = s">$id" + description.map(d => s" $d").getOrElse("")

  def format(width: Int = 60) = {
    val lines = new ListBuffer[String]()
    lines += this.header
    var k = 0
    while (k < sequence.size) {
      lines += sequence.slice(k, k+width)
      k += width
    }
    lines.mkString("\n")
  }

  override def toString = format()
}

object FastaReader {

  private def newRecord(header: String, sequence: Option[Seq[String]]) = {
    val fields = header.stripPrefix(">").split(" ")
    if (fields.size == 0) {
      throw new Exception(s"Cannot parse as FASTA header:\n${header}")
    } else if (sequence.isEmpty || sequence.get.size == 0) {
      throw new Exception(s"Empty sequence for FASTA record:\n${header}")
    } else {
      val desc = if (fields.size > 1) Some(fields.drop(1).mkString(" ")) else None
      FastaRecord(fields(0), sequence.get.mkString(""), desc)
    }
  }

  def apply(path: Path): Seq[FastaRecord] = {
    val s = Source.fromFile(path.toFile)
    val records = new ListBuffer[FastaRecord]()
    var currentHeader: Option[String] = None
    var currentSequence: Option[ListBuffer[String]] = None
    s.getLines.foreach { line =>
      if (line.trim == "") {
        throw new Exception("Blank lines not allowed")
      } else if (line.startsWith(">")) {
        if (currentHeader.isDefined) {
          records += newRecord(currentHeader.get, currentSequence)
        }
        currentHeader = Some(line)
        currentSequence = Some(new ListBuffer[String]())
      } else if (currentHeader.isDefined) {
        ///currentSequence = currentSequence.map(s => s ++ Seq(line))
        currentSequence.get += line
      } else {
        throw new Exception(s"Illegal FASTA line:\n$line")
      }
    }
    if (currentHeader.isEmpty) throw new Exception("No valid FASTA records")
    records += newRecord(currentHeader.get, currentSequence)
    records
  }
}
