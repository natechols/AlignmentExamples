// Author: Nat Echols
// This source code is in the public domain.

#include "nwalign.hpp"
#include "blosum.hpp"

namespace nwalign {

  namespace impl {

    void BlosumMatrix::setData(const int rawData[24][24]) {
      for (size_t i = 0; i < 24; i++) {
        for (size_t j = 0; j < 24; j++) {
          data[i][j] = rawData[i][j];
        }
      }
    }

    BlosumMatrix::BlosumMatrix(std::string name) {
      for (size_t i = 0; i < 24; i++) {
        indexMap[matrix_cols[i]] = i;
      }
      if (name == "BLOSUM50") {
        setData(blosum50_data);
      } else {
        throw new std::runtime_error("Matrix " + name + " not found.");
      }
    }

    int BlosumMatrix::operator()(char aa, char bb) {
      size_t i = indexMap[aa];
      size_t j = indexMap[bb];
      return data[i][j];
    }
   
  }

  NeedlemanWunsch::NeedlemanWunsch(std::string seq1,
                                   std::string seq2,
                                   std::string scoreMatrixName = "BLOSUM50",
                                   int gapPenalty = 8) {
    impl::BlosumMatrix scoreMatrix(scoreMatrixName);
    std::vector<std::vector<MatrixCell>> f;
    for (size_t i = 0; i <= seq1.size(); i++) {
      f.push_back(std::vector<MatrixCell>());
      f[i].push_back(MatrixCell(-i*gapPenalty, CELL_LEFT));
      for (size_t j = 1; j <= seq2.size(); j++) {
        f[i].push_back(MatrixCell(-j*gapPenalty, CELL_UP));
      }
    }
    for (size_t i = 1; i <= seq1.size(); i++) {
      char aa = seq1[i-1];
      for (size_t j = 1; j <= seq2.size(); j++) {
        char bb = seq2[i-1];
        int scoreDiag = f[i-1][j-1].score + scoreMatrix(aa, bb);
        int scoreLeft = f[i-1][j].score - gapPenalty;
        int scoreUp = f[i][j-1].score - gapPenalty;
        if ((scoreDiag >= scoreLeft) && (scoreDiag >= scoreUp)) {
          f[i][j].score = scoreDiag;
          f[i][j].cell_pointer = CELL_DIAG;
        } else if (scoreLeft >= scoreUp) {
          f[i][j].score = scoreLeft;
          f[i][j].cell_pointer = CELL_LEFT;
        } else {
          f[i][j].score = scoreUp;
          f[i][j].cell_pointer = CELL_UP;
        }
      }
    }
    std::vector<char> tmp1Aligned;
    std::vector<char> tmp2Aligned;
    size_t i = seq1.size();
    size_t j = seq2.size();
    while ((i >= 0) && (j >= 0) && !((i == 0) && (j == 0))) {
      if (f[i][j].cell_pointer == CELL_DIAG) {
        tmp1Aligned.push_back(seq1[i-1]);
        tmp2Aligned.push_back(seq2[j-1]);
        i--;
        j--;
      } else if (f[i][j].cell_pointer == CELL_LEFT) {
        tmp1Aligned.push_back(seq1[i-1]);
        tmp2Aligned.push_back('-');
        i--;
      } else {
        tmp1Aligned.push_back('-');
        tmp2Aligned.push_back(seq2[j-1]);
        j--;
      }
    }
    seq1Aligned = std::string(tmp1Aligned.rbegin(), tmp1Aligned.rend());
    seq2Aligned = std::string(tmp2Aligned.rbegin(), tmp2Aligned.rend());
  }

  double NeedlemanWunsch::identity() {
    return 0.0;
  }
}
