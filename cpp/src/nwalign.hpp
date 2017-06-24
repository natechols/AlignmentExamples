// Author: Nat Echols
// This source code is in the public domain.

#ifndef NWALIGN_H
#define NWALIGN_H

#include <string>
#include <vector>
#include <tuple>
#include <map>

#define CELL_DIAG 0
#define CELL_LEFT -1
#define CELL_UP 1

namespace nwalign {

  namespace impl {
    class BlosumMatrix {
      private:
        std::map<char, size_t> indexMap;
        int data[24][24];
        void setData(const int rawData[24][24]);

      public:
        BlosumMatrix(std::string name);
        int operator()(char aa, char bb);
    };
  }

  class NeedlemanWunsch {
    private:
      std::string seq1Aligned;
      std::string seq2Aligned;
      struct MatrixCell {
        int score;
        int cell_pointer;
        MatrixCell(int s, int p): score(s), cell_pointer(p) {}
      };

    public:
      NeedlemanWunsch(std::string seq1,
                      std::string seq2,
                      std::string scoreMatrixName,
                      int gapPenalty);

      std::tuple<std::string, std::string> alignments() {
        return std::tuple<std::string, std::string>(seq1Aligned, seq2Aligned);
      }
      size_t alignmentLength() {
        return seq1Aligned.size();
      }
      double identity();
      void show();
  };
}

#endif
