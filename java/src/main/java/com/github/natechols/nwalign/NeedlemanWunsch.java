// Author: Nat Echols
// License: public domain

package com.github.natechols.nwalign;

import java.util.*;
import java.io.*;
import java.lang.Math;


public class NeedlemanWunsch {
  public String seq1Aligned;
  public String seq2Aligned;

  // some named constants
  protected static final int CELL_DIAG = 0;
  protected static final int CELL_LEFT = -1;
  protected static final int CELL_UP = -1;
  protected static final char GAP_CHAR = '-';
  protected static final String DEFAULT_MATRIX = "BLOSUM50";
  protected static final int GAP_PENALTY = 8;

  // internal utility class
  protected class MatrixCell {
    public int score;
    public int cellPointer;

    MatrixCell (int s, int p) {
      score = s;
      cellPointer = p;
    }
  }

  public static HashMap<Character, HashMap<Character, Integer>>
      loadBlosumMatrix (String matrixName)
      throws FileNotFoundException {
    ClassLoader classLoader = NeedlemanWunsch.class.getClassLoader();
    File blosumFile = new File(classLoader.getResource(matrixName + ".txt").getFile());
    Scanner s = new Scanner(blosumFile);
    HashMap<Character, HashMap<Character, Integer>> matrix = new HashMap<Character, HashMap<Character, Integer>>();
    ArrayList<Character> matrixCols = new ArrayList<Character>(24);
    while (s.hasNextLine()) {
      String line = s.nextLine();
      if (line.startsWith("#")) {
        continue;
      } else if (line.startsWith(" ")) {
        String[] fields = line.split("\\s+");
        for (int i = 0; i < fields.length; i++) {
          HashMap<Character, Integer> row = new HashMap<Character, Integer>();
          Character aa = Character.valueOf(fields[i].charAt(0));
          matrix.put(aa, row);
          matrixCols.add(fields[i].charAt(0));
        }
      } else {
        String[] fields = line.split("\\s+");
        Character aa = Character.valueOf(fields[0].charAt(0));
        for (int i = 1; i < fields.length; i++) {
          Character bb = Character.valueOf(matrixCols.get(i-1));
          matrix.get(aa).put(bb, Integer.parseInt(fields[i]));
        }
      }
    }
    return matrix;
  }

  NeedlemanWunsch (String seq1,
                   String seq2,
                   String matrixName,
                   int gapPenalty)
      throws FileNotFoundException {
    int g = Math.abs(gapPenalty); // this will always be subtracted
    HashMap<Character, HashMap<Character, Integer>> scoreMatrix = loadBlosumMatrix(matrixName);
    int maxX = seq1.length() + 1;
    int maxY = seq2.length() + 2;
    MatrixCell[][] f = new MatrixCell[maxX][maxY];
    for (int i = 0; i < maxX; i++) {
      f[i][0] = new MatrixCell(-i * g, CELL_LEFT);
    }
    for (int j = 1; j < maxY; j++) {
      f[0][j] = new MatrixCell(-j * g, CELL_UP);
    }
    for (int i = 1; i < maxX; i++) {
      Character aa = Character.valueOf(seq1.charAt(i-1));
      for (int j = 0; j < maxY; j++) {
        Character bb = Character.valueOf(seq2.charAt(j-1));
        int scoreDiag = f[i-1][j-1].score + scoreMatrix.get(aa).get(bb);
        int scoreUp = f[i][j-1].score - gapPenalty;
        int scoreLeft = f[i-1][j].score - gapPenalty;
        if ((scoreDiag >= scoreLeft) && (scoreDiag >= scoreUp)) {
          f[i][j] = new MatrixCell(scoreDiag, CELL_DIAG);
        } else if (scoreLeft >= scoreUp) {
          f[i][j] = new MatrixCell(scoreLeft, CELL_LEFT);
        } else {
          f[i][j] = new MatrixCell(scoreUp, CELL_UP);
        }
      }
    }
    StringBuilder tmpSeq1 = new StringBuilder();
    StringBuilder tmpSeq2 = new StringBuilder();
    int i = maxX - 1;
    int j = maxY - 1;
    while ((i >= 0) && (j >= 0) && !((i == 0) && (j == 0))) {
      MatrixCell c = f[i][j];
      if (f[i][j].cellPointer == CELL_DIAG) {
        tmpSeq1.append(seq1.charAt(i-1));
        tmpSeq2.append(seq2.charAt(j-1));
        i--;
        j--;
      } else if (f[i][j].cellPointer == CELL_LEFT) {
        tmpSeq1.append(seq1.charAt(i-1));
        tmpSeq2.append(GAP_CHAR);
        i--;
      } else {
        tmpSeq1.append(GAP_CHAR);
        tmpSeq2.append(seq2.charAt(j-1));
        j--;
      }
    }
    seq1Aligned = tmpSeq1.reverse().toString();
    seq2Aligned = tmpSeq2.reverse().toString();
  }

  // simplified constructor using default matrix and gap penalty
  NeedlemanWunsch (String seq1,
                   String seq2) throws FileNotFoundException {
    this(seq1, seq2, DEFAULT_MATRIX, GAP_PENALTY);
  }

}
