
#include "nwalign.h"

int main(void) {
  char seq1[] = "HEAGAWGHEE";
  char seq2[] = "PAWHEAE";
  sequence_alignment *result = NeedlemanWunsch(seq1, seq2, "BLOSUM50", 8);
  return 0;
}
