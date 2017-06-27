// Author: Nat Echols
// License: public domain

package nwalign

import "math"
import "fmt"

type Alignment struct {
  seq1Aligned string
  seq2Aligned string
}

type MatrixCell struct {
  score int
  cellPointer int
}

func getMatrix (matrixName string) map[byte]map[byte]int {
  var matrixCols = [24]byte{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'}
  matrices := blosumMatrices()
  var rawData [24][24]int = matrices[matrixName]
  scoreMatrix := make(map[byte]map[byte]int)
  for i, aa := range matrixCols {
    scoreMatrix[aa] = make(map[byte]int)
    for j, bb  := range matrixCols {
      scoreMatrix[aa][bb] = rawData[i][j]
    }
  }
  return scoreMatrix
}

func NeedlemanWunsch (seq1 string,
                      seq2 string,
                      matrixName string,
                      gapPenalty int) Alignment {
  scoreMatrix := getMatrix(matrixName)
  const CELL_LEFT int = -1
  const CELL_UP int = 1
  const CELL_DIAG int = 0
  maxX := len(seq1) + 1
  maxY := len(seq2) + 1
  var f = make([][]MatrixCell, maxX)
  for i:= 0; i < maxX; i++ {
    f[i] = make([]MatrixCell, maxY)
    f[i][0] = MatrixCell{-i * gapPenalty, CELL_LEFT}
  }
  for j:= 1; j < maxY; j++ {
    f[0][j] = MatrixCell{-j * gapPenalty, CELL_UP}
  }
  for i := 1; i < maxX; i++ {
    aa := seq1[i-1]
    for j := 1; j < maxY; j++ {
      bb := seq2[j-1]
      cellDiag := MatrixCell{f[i-1][j-1].score + scoreMatrix[aa][bb], CELL_DIAG}
      cellLeft := MatrixCell{f[i-1][j].score - gapPenalty, CELL_LEFT}
      cellUp := MatrixCell{f[i][j-1].score - gapPenalty, CELL_UP}
      if cellDiag.score >= cellLeft.score && cellDiag.score >= cellUp.score {
        f[i][j] = cellDiag
      } else if cellLeft.score >= cellUp.score {
        f[i][j] = cellLeft
      } else {
        f[i][j] = cellUp
      }
    }
  }
  var tmp1Aligned []byte = make([]byte, 0, maxX+maxY-2)
  var tmp2Aligned []byte = make([]byte, 0, maxX+maxY-2)
  i := len(seq1)
  j := len(seq2)
  const GAP_CHAR byte = '-'
  for {
    if f[i][j].cellPointer == CELL_DIAG {
      tmp1Aligned = append(tmp1Aligned, seq1[i-1])
      tmp2Aligned = append(tmp2Aligned, seq2[j-1])
      i--
      j--
    } else if f[i][j].cellPointer == CELL_LEFT {
      tmp1Aligned = append(tmp1Aligned, seq1[i-1])
      tmp2Aligned = append(tmp2Aligned, GAP_CHAR)
      i--
    } else {
      tmp1Aligned = append(tmp1Aligned, GAP_CHAR)
      tmp2Aligned = append(tmp2Aligned, seq2[j-1])
      j--
    }
    if (i == 0 && j == 0) || i < 0 || j < 0 {
      break
    }
  }
  var n = len(tmp1Aligned)
  var seq1Aligned = make([]byte, n)
  var seq2Aligned = make([]byte, n)
  for i := 0; i < n; i++ {
    seq1Aligned[i] = tmp1Aligned[n-i-1]
  }
  for j := 0; j < n; j++ {
    seq2Aligned[j] = tmp2Aligned[n-j-1]
  }
  return Alignment{string(seq1Aligned[:]), string(seq2Aligned[:])}
}

func (a *Alignment) length() int {
  return len(a.seq1Aligned)
}

func (a *Alignment) identity() float64 {
  nMm, nDel, nIns := 0, 0, 0
  for i := 0; i < a.length(); i++ {
    aa := a.seq1Aligned[i]
    bb := a.seq2Aligned[i]
    if aa == '-' {
      nIns++
    } else if bb == '-' {
      nDel++
    } else if (aa != bb) {
      nMm++
    }
  }
  return 1.0 - float64(nMm + nDel + nIns) / float64(a.length())
}

func getSeqMatches (s1 string, s2 string) string {
  var matchChars = make([]byte, len(s1))
  for i := 0; i < len(s1); i++ {
    if (s1[i] == s2[i]) {
      matchChars = append(matchChars, '*')
    } else {
      matchChars = append(matchChars, ' ')
    }
  }
  return string(matchChars)
}

// Display alignment (with matches highlighted) formatted for a narrow terminal
func (a *Alignment) show (seqid1 string, seqid2 string) {
  const SEQ_WIDTH = 50
  var nRows = int(math.Ceil(float64(a.length()) / float64(SEQ_WIDTH)))
  var idWidth = int(math.Max(float64(len(seqid1)), float64(len(seqid2))))
  var fmtString = fmt.Sprintf("%%%ds %%s\n", idWidth)
  for k := 0; k < nRows; k++ {
    var start = k * SEQ_WIDTH
    var end = int(math.Min(float64((k + 1) * SEQ_WIDTH), float64(a.length())))
    var s1 = a.seq1Aligned[start:end]
    var s2 = a.seq2Aligned[start:end]
    var seqMatches = getSeqMatches(s1, s2)
    if (k > 0) {
      fmt.Println("\n")
    }
    fmt.Printf(fmtString, seqid1, s1)
    fmt.Printf(fmtString, "", seqMatches)
    fmt.Printf(fmtString, seqid2, s2)
  }
}
