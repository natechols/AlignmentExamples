
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "nwalign.h"


int **get_blosum_matrix(char *matrix_name) {
  _Bool have_env = true;
  int i, j, k, m, n, score;
  int **full_matrix;
  int raw_data[24][24];
  char matrix_cols[24];
  char *blosum_file, *blosum_dir, line[256], *linep;
  FILE *f;
  blosum_dir = getenv("BLOSUM_DATA");
  if (blosum_dir == NULL) {
    have_env = false;
    blosum_dir = (char *) malloc(sizeof(char) * 2);
    strcpy(blosum_dir, ".");
  }
  blosum_file = (char *) malloc(strlen(matrix_name) + strlen(blosum_dir) + 6);
  sprintf(blosum_file, "%s/%s.txt", blosum_dir, matrix_name);
  f = fopen(blosum_file, "r");
  if (f == NULL) {
    fprintf(stderr, "Can't read file %s\n", blosum_file);
    return NULL;
  }
  free(blosum_file);
  if (!have_env) {
    free(blosum_dir);
  }
  i = 0;
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '#') {
      continue;
    } else if (line[0] == ' ') {
      j = 0;
      for (k = 1; k < strlen(line); k++) {
        if ((line[k] != ' ') && (line[k] != '\n') && (line[k] != '\r')) {
          matrix_cols[j++] = line[k];
        }
      }
      if (j != 24) {
        fprintf(stderr, "Expect 24 columns, got %d\n", j);
        return NULL;
      }
    } else {
      j = 0;
      k = 1;
      n = strlen(line);
      while (k < n) {
        if (line[k] != ' ') {
          break;
        }
        k++;
      }
      while ((j < 24) && (k < n - 1)) {
        linep = &(line[k]);
        m = k + 1;
        while (m < n - 1) {
          if (line[m] == ' ') {
            line[m] = '\0';
            break;
          }
          m++;
        }
        sscanf(linep, "%d", &score);
        raw_data[i][j] = score;
        k = m + 1;
        j++;
      }
      if (j != 24) {
        fprintf(stderr, "Expected 24 integer values, got %d\n", j);
        return NULL;
      }
      i++;
    }
  }
  fclose(f);
  full_matrix = (int **) malloc(sizeof(int *) * 256);
  for (i = 0; i < 256; i++) {
    full_matrix[i] = (int *) malloc(sizeof(int) * 256);
    for (j = 0; j < 256; j++) {
      full_matrix[i][j] = BLOSUM_UNDEF;
    }
  }
  for (i = 0; i < 24; i++) {
    for (j = 0; j < 24; j++) {
      m = (int) matrix_cols[i];
      n = (int) matrix_cols[j];
      full_matrix[m][n] = raw_data[i][j];
    }
  }
  return full_matrix;
}

sequence_alignment *NeedlemanWunsch(char *seq1, char *seq2, int **score_matrix,
                                    int gap_penalty) {
  int i = 0, j = 0;
  int max_x, max_y;
  int aln_size = 0;
  unsigned short aa, bb;
  sequence_alignment *result;
  matrix_cell **f;
  matrix_cell *cell_tmp;
  int score_diag, score_left, score_up;
  char_list *tmp1_aligned, *tmp2_aligned, *last1_aligned, *last2_aligned;
  char_list *seq1_traceback, *seq2_traceback;
  if (score_matrix == NULL) {
    fprintf(stderr, "Scoring matrix cannot be NULL");
    return NULL;
  }
  max_x = strlen(seq1) + 1;
  max_y = strlen(seq2) + 1;
  f = (matrix_cell **) malloc(sizeof(matrix_cell *) * max_x);
  for (i = 0; i < max_x; i++) {
    f[i] = (matrix_cell *) malloc(sizeof(matrix_cell) * max_y);
    f[i][0].score = -i*gap_penalty;
    f[i][0].cell_pointer = CELL_LEFT;
  }
  for (j = 0; j < max_y; j++) {
    f[0][j].score = -j*gap_penalty;
    f[0][j].cell_pointer = CELL_UP;
  }
  for (i = 1; i < max_x; i++) {
    aa = (unsigned short) seq1[i-1];
    for (j = 1; j < max_y; j++) {
      bb = (unsigned short) seq2[j-1];
      score_diag = f[i-1][j-1].score + score_matrix[aa][bb];
      score_left = f[i-1][j].score - gap_penalty;
      score_up = f[i][j-1].score - gap_penalty;
      cell_tmp = &(f[i][j]);
      if ((score_diag >= score_left) && (score_diag >= score_up)) {
        cell_tmp->score = score_diag;
        cell_tmp->cell_pointer = CELL_DIAG;
      } else if (score_left >= score_up) {
        cell_tmp->score = score_left;
        cell_tmp->cell_pointer = CELL_LEFT;
      } else {
        cell_tmp->score = score_up;
        cell_tmp->cell_pointer = CELL_UP;
      }
    }
  }
  seq1_traceback = seq2_traceback = NULL;
  i = max_x - 1;
  j = max_y - 1;
  while ((i >= 0) && (j >= 0) && !((i == 0) && (j == 0))) {
    cell_tmp = &(f[i][j]);
    tmp1_aligned = (char_list *) malloc(sizeof(char_list));
    tmp2_aligned = (char_list *) malloc(sizeof(char_list));
    if (cell_tmp->cell_pointer == CELL_DIAG) {
      tmp1_aligned->c = seq1[i-1];
      tmp2_aligned->c = seq2[j-1];
      i--;
      j--;
    } else if (cell_tmp->cell_pointer == CELL_LEFT) {
      tmp1_aligned->c = seq1[i-1];
      tmp2_aligned->c = '-';
      i--;
    } else {
      tmp1_aligned->c = '-';
      tmp2_aligned->c = seq2[j-1];
      j--;
    }
    if (seq1_traceback == NULL) {
      seq1_traceback = tmp1_aligned;
      seq2_traceback = tmp2_aligned;
    } else {
      last1_aligned->next = (void *) tmp1_aligned;
      last2_aligned->next = (void *) tmp2_aligned;
    }
    last1_aligned = tmp1_aligned;
    last2_aligned = tmp2_aligned;
    aln_size++;
  }
  result = (sequence_alignment *) malloc(sizeof(sequence_alignment));
  result->seq1aligned = (char *) malloc(sizeof(char) * (aln_size+1));
  result->seq2aligned = (char *) malloc(sizeof(char) * (aln_size+1));
  result->aln_size = aln_size;
  tmp1_aligned = seq1_traceback;
  tmp2_aligned = seq2_traceback;
  for (i = 0; i < aln_size; i++) {
    result->seq1aligned[aln_size-i-1] = tmp1_aligned->c;
    result->seq2aligned[aln_size-i-1] = tmp2_aligned->c;
    last1_aligned = tmp1_aligned;
    last2_aligned = tmp2_aligned;
    tmp1_aligned = (char_list *) tmp1_aligned->next;
    tmp2_aligned = (char_list *) tmp2_aligned->next;
    free(last1_aligned);
    free(last2_aligned);
  }
  result->seq1aligned[aln_size] = '\0';
  result->seq2aligned[aln_size] = '\0';
  for (i = 0; i < 256; i++) {
    free(score_matrix[i]);
  }
  free(score_matrix);
  for (i = 0; i < max_x; i++) {
    free(f[i]);
  }
  free(f);
  return result;
}
