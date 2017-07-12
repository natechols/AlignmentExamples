
#include "nwalign.h"

int main(void) {
  char seq1[] = "HEAGAWGHEE";
  char seq2[] = "PAWHEAE";
  int **score_matrix = get_blosum_matrix("BLOSUM50");
  sequence_alignment *result = NeedlemanWunsch(seq1, seq2, score_matrix, 8);
  return 0;
}
