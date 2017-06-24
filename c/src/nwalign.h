/* Author: Nat Echols
 * This source code is in the public domain.
 */

#ifndef NWALIGN_H
#define NWALIGN_H

#define CELL_DIAG 0
#define CELL_LEFT -1
#define CELL_UP 1

#define BLOSUM_UNDEF -9999999

typedef struct {
  char *seq1aligned;
  char *seq2aligned;
  int aln_size;
} sequence_alignment;

typedef struct {
  char c;
  void *next;
} char_list;

typedef struct {
  int score;
  short cell_pointer;
} matrix_cell;

sequence_alignment *NeedlemanWunsch(char *seq1, char *seq2, char *matrix_name,
                                    int gap_penalty);

int **get_blosum_matrix(char *matrix_name);

#endif
