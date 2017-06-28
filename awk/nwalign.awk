#!/usr/bin/awk -f
#
# Author: Nat Echols
# License: LGPL due to inclusion of BLOSUM50

function nwalign (seq1, seq2, gapPenalty) {
  # FIXME is there a way to load BLOSUM50 externally?
  CELL_UP = -1
  CELL_DIAG = 0
  CELL_LEFT = 1
  split("A R N D C Q E G H I L K M F P S T W Y V B Z X *", matrixCols, " ");
  for (i = 1; i <= length(matrixCols); i++) {
    aa = matrixCols[i]
    blosum[aa][0] = "";
  }
  split("5 -2 -1 -2 -1 -1 -1 0 -2 -1 -2 -1 -1 -3 -1 1 0 -3 -2 0 -2 -1 -1 -5", blosum["A"], " ");
  split("-2 7 -1 -2 -4 1 0 -3 0 -4 -3 3 -2 -3 -3 -1 -1 -3 -1 -3 -1 0 -1 -5", blosum["R"], " ");
  split("-1 -1 7 2 -2 0 0 0 1 -3 -4 0 -2 -4 -2 1 0 -4 -2 -3 4 0 -1 -5", blosum["N"], " ");
  split("-2 -2 2 8 -4 0 2 -1 -1 -4 -4 -1 -4 -5 -1 0 -1 -5 -3 -4 5 1 -1 -5", blosum["D"], " ");
  split("-1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -3 -2 -5", blosum["C"], " ");
  split("-1 1 0 0 -3 7 2 -2 1 -3 -2 2 0 -4 -1 0 -1 -1 -1 -3 0 4 -1 -5", blosum["Q"], " ");
  split("-1 0 0 2 -3 2 6 -3 0 -4 -3 1 -2 -3 -1 -1 -1 -3 -2 -3 1 5 -1 -5", blosum["E"], " ");
  split("0 -3 0 -1 -3 -2 -3 8 -2 -4 -4 -2 -3 -4 -2 0 -2 -3 -3 -4 -1 -2 -2 -5", blosum["G"], " ");
  split("-2 0 1 -1 -3 1 0 -2 10 -4 -3 0 -1 -1 -2 -1 -2 -3 2 -4 0 0 -1 -5", blosum["H"], " ");
  split("-1 -4 -3 -4 -2 -3 -4 -4 -4 5 2 -3 2 0 -3 -3 -1 -3 -1 4 -4 -3 -1 -5", blosum["I"], " ");
  split("-2 -3 -4 -4 -2 -2 -3 -4 -3 2 5 -3 3 1 -4 -3 -1 -2 -1 1 -4 -3 -1 -5", blosum["L"], " ");
  split("-1 3 0 -1 -3 2 1 -2 0 -3 -3 6 -2 -4 -1 0 -1 -3 -2 -3 0 1 -1 -5", blosum["K"], " ");
  split("-1 -2 -2 -4 -2 0 -2 -3 -1 2 3 -2 7 0 -3 -2 -1 -1 0 1 -3 -1 -1 -5", blosum["M"], " ");
  split("-3 -3 -4 -5 -2 -4 -3 -4 -1 0 1 -4 0 8 -4 -3 -2 1 4 -1 -4 -4 -2 -5", blosum["F"], " ");
  split("-1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -1 -2 -5", blosum["P"], " ");
  split("1 -1 1 0 -1 0 -1 0 -1 -3 -3 0 -2 -3 -1 5 2 -4 -2 -2 0 0 -1 -5", blosum["S"], " ");
  split("0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 2 5 -3 -2 0 0 -1 0 -5", blosum["T"], " ");
  split("-3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1 1 -4 -4 -3 15 2 -3 -5 -2 -3 -5", blosum["W"], " ");
  split("-2 -1 -2 -3 -3 -1 -2 -3 2 -1 -1 -2 0 4 -3 -2 -2 2 8 -1 -3 -2 -1 -5", blosum["Y"], " ");
  split("0 -3 -3 -4 -1 -3 -3 -4 -4 4 1 -3 1 -1 -3 -2 0 -3 -1 5 -4 -3 -1 -5", blosum["V"], " ");
  split("-2 -1 4 5 -3 0 1 -1 0 -4 -4 0 -3 -4 -2 0 0 -5 -3 -4 5 2 -1 -5", blosum["B"], " ");
  split("-1 0 0 1 -3 4 5 -2 0 -3 -3 1 -1 -4 -1 0 -1 -2 -2 -3 2 5 -1 -5", blosum["Z"], " ");
  split("-1 -1 -1 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -1 0 -3 -1 -1 -1 -1 -1 -5", blosum["X"], " ");
  split("-5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 1", blosum["*"], " ");
  for (i = 1; i <= length(matrixCols); i++) {
    aa = matrixCols[i];
    for (j = 1; j <= length(matrixCols); j++) {
      bb = matrixCols[j];
      blosum[aa][bb] = blosum[aa][j];
    }
  }
  maxX = length(seq1) + 1;
  maxY = length(seq2) + 1;
  for (i = 0; i < maxX; i++) {
    scores[i][0] = -i * gapPenalty;
    ptrs[i][0] = CELL_LEFT;
  }
  for (j = 0; j < maxY; j++) {
    scores[0][j] = -j * gapPenalty;
    ptrs[0][j] = CELL_UP;
  }
  split(seq1, seq1a, "");
  split(seq2, seq2a, "");
  for (i = 1; i < maxX; i++) {
    aa = seq1a[i];
    for (j = 1; j < maxY; j++) {
      bb = seq2a[j];
      scoreDiag = scores[i-1][j-1] + blosum[aa][bb];
      scoreLeft = scores[i-1][j] - gapPenalty;
      scoreUp = scores[i][j-1] - gapPenalty;
      if ((scoreDiag >= scoreLeft) && (scoreDiag >= scoreUp)) {
        scores[i][j] = scoreDiag;
        ptrs[i][j] = CELL_DIAG;
      } else if (scoreLeft >= scoreUp) {
        scores[i][j] = scoreLeft;
        ptrs[i][j] = CELL_LEFT;
      } else {
        scores[i][j] = scoreUp;
        ptrs[i][j] = CELL_UP;
      }
    }
  }
  i = length(seq1);
  j = length(seq2);
  k = 0;
  while ((i >= 0) && (j >= 0) && !((i == 0) && (j == 0))) {
    if (ptrs[i][j] == CELL_DIAG) {
      tmp1[k] = seq1a[i];
      tmp2[k] = seq2a[j];
      i--;
      j--;
    } else if (ptrs[i][j] == CELL_LEFT) {
      tmp1[k] = seq1a[i];
      tmp2[k] = "-";
      i--;
    } else {
      tmp1[k] = "-";
      tmp2[k] = seq2a[j];
      j--;
    }
    k++;
  }
  seq1aligned = "";
  seq2aligned = "";
  for (i = 0; i < k; i++) {
    seq1aligned = seq1aligned "" tmp1[k-i-1];
    seq2aligned = seq2aligned "" tmp2[k-i-1];
  }
  return seq1aligned "\n" seq2aligned;
}

function getMatches(s1, s2) {
  split(s1, s1a, "");
  split(s2, s2a, "");
  matches = "";
  for (i = 1; i <= length(s1a); i++) {
    if (s1a[i] == s2a[i]) {
      matches = matches "*";
    } else {
      matches = matches " ";
    }
  }
  return matches;
}

function identity(s1, s2) {
  nMm = nDel = nIns = 0;
  split(s1, s1a, "");
  split(s2, s2a, "");
  matches = "";
  for (i = 1; i <= length(s1a); i++) {
    if (s1a[i] == "-") {
      nIns++;
    } else if (s2a[i] == "-") {
      nDel++;
    } else if (s1a[i] != s2a[i]) {
      nMm++;
    }
  }
  return 1.0 - (nMm + nDel + nIns) / length(s1);
}

BEGIN {
  k = 0;
  ids[0] = "";
  seqs[0] = "";
}

# Main loop: FASTA parsing
{
  if ($1 ~ /^>/) {
    k += 1;
    sub(/^>/, "");
    sub(/ .*/, "");
    ids[k] = $0;
    seqs[k] = "";
  } else if ($0 ~ /[A-Z]{1,}/) {
    if (ids[k] == "") {
      print "Missing header record";
      exit 1;
    }
    seqs[k] = seqs[k] $0;
  } else {
    print "Can't interpret" $0;
    exit 1;
  }
}

END {
  if (k != 2) {
    print "Require two sequences to align";
    exit 1;
  }
  for (i = 1; i <= k; i++) {
    if (length(ids[i]) == 0) {
      print "Missing ID for sequence " i;
      exit 1;
    } else if (length(seqs[i]) == 0) {
      print "Missing sequence " i;
      exit 1;
    }
  } 
  result = nwalign(seqs[1], seqs[2], 8);
  split(result, aln, "\n");
  k = 0;
  idWidth = length(ids[1]);
  if (length(ids[2]) > idWidth) {
    idWidth = length(ids[2]);
  }
  fmt = sprintf("%%-%ds  %%s\n", idWidth);
  # XXX random observation: calling this 'matches' does not work because it
  # will always be empty!
  matchChars = getMatches(aln[1], aln[2]);
  printf("Sequence identity = %.2f%%\n\n", identity(aln[1], aln[2])*100);
  while (k < length(aln[1])) {
    if (k > 0) printf("\n\n");
    printf(fmt, ids[1], substr(aln[1], k+1, 50));
    printf(fmt, "", substr(matchChars, k+1, 50));
    printf(fmt, ids[2], substr(aln[2], k+1, 50));
    k += 50;
  }
  exit 0;
}
