// Author: Nat Echols
// License: public domain

var nwalign = (function() {

  var CELL_DIAG = 0;
  var CELL_LEFT = -1;
  var CELL_UP = 1;
  var GAP_CHAR = '-';

  function Alignment (seq1aligned, seq2aligned) {
    this.seq1aligned = seq1aligned;
    this.seq2aligned = seq2aligned;

    this.length = function () {
      return this.seq1aligned.length;
    };

    this.identity = function () {
      var nMismatch = 0, nIns = 0, nDel = 0;
      for (var i = 0; i < this.seq1aligned.length; i++) {
        var aa = this.seq1aligned.charAt(i);
        var bb = this.seq2aligned.charAt(i);
        if (aa == "-") {
          nIns++;
        } else if (bb == "-") {
          nDel++;
        } else if (aa != bb) {
          nMismatch++;
        }
      }
      return 1.0 - (nIns + nDel + nMismatch) / (this.length() * 1.0);
    };

    function getMatches () {
      var matches = "";
      for (var i = 0; i < seq1aligned.length; i++) {
        var aa = seq1aligned.charAt(i);
        var bb = seq2aligned.charAt(i);
        if (aa == bb) {
          matches += "*";
        } else {
          matches += " ";
        }
      }
      return matches;
    }

    // XXX newer JS implementations have String.padEnd()
    function padString (s, targetLen) {
      var pad = "";
      for (var i = s.length; i < targetLen; i++) {
        pad += " ";
      }
      return s + pad;
    };

    this.format = function (seqid1, seqid2) {
      var matches = getMatches();
      var idlen = Math.max(seqid1.length, seqid2.length);
      var headerm = padString("", idlen);
      var header1 = padString(seqid1, idlen);
      var header2 = padString(seqid2, idlen);
      var rows = new Array();
      var k = 0;
      while (k < this.length()) {
        var s1 = this.seq1aligned.substring(k, k + 50);
        var s2 = this.seq2aligned.substring(k, k + 50);
        var m = matches.substring(k, k + 50);
        rows.push([
          header1 + "  " + s1,
          headerm + "  " + m,
          header2 + "  " + s2
        ].join("\n"));
        k += 50;
      }
      return rows.join("\n\n\n");
    };
  };

  function NeedlemanWunsch (seq1, seq2, scoreMatrix, gapPenalty) {
    function MatrixCell (score, cellPointer) {
      this.score = score;
      this.cellPointer = cellPointer;
    }

    var maxX = seq1.length + 1;
    var maxY = seq2.length + 1;
    var f = new Array(maxX);
    for (var i = 0; i < maxX; i++) {
      f[i] = new Array(maxY);
      f[i][0] = new MatrixCell(-i * gapPenalty, CELL_LEFT);
    }
    for (var j = 1; j < maxY; j++) {
      f[0][j] = new MatrixCell(-j * gapPenalty, CELL_UP);
    }
    for (var i = 1; i < maxX; i++) {
      var aa = seq1.charAt(i-1);
      for (var j = 1; j < maxY; j++) {
        var bb = seq2.charAt(j-1);
        f[i][j] = [
          new MatrixCell(f[i-1][j-1].score + scoreMatrix[aa][bb], CELL_DIAG),
          new MatrixCell(f[i-1][j].score - gapPenalty, CELL_LEFT),
          new MatrixCell(f[i][j-1].score - gapPenalty, CELL_UP)
        ].sort(function (a,b) { return b.score - a.score; })[0];
      }
    }
    var i = seq1.length;
    var j = seq2.length;
    var tmpSeq1 = "";
    var tmpSeq2 = "";
    while ((i >= 0) && (j >= 0) && !((i == 0) && (j == 0))) {
      var c = f[i][j];
      if (c.cellPointer == CELL_DIAG) {
        tmpSeq1 += seq1.charAt(i - 1);
        tmpSeq2 += seq2.charAt(j - 1);
        i--;
        j--;
      } else if (c.cellPointer == CELL_LEFT) {
        tmpSeq1 += seq1.charAt(i - 1);
        tmpSeq2 += GAP_CHAR;
        i--;
      } else {
        tmpSeq1 += GAP_CHAR;
        tmpSeq2 += seq2.charAt(j - 1);
        j--;
      }
    }
    var seq1aligned = tmpSeq1.split("").reverse().join("");
    var seq2aligned = tmpSeq2.split("").reverse().join("");
    return new Alignment(seq1aligned, seq2aligned);
  };

  return {
    NeedlemanWunsch: NeedlemanWunsch
  };

}());
