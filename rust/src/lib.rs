// Author: Nat Echols
// License: public domain

pub mod fasta;

use std::cmp;

mod blosum {

  use std::collections::HashMap;
  type ScoreMatrix = HashMap<(u8, u8), i32>;

  /// Read a BLOSUM matrix specified by name.  The allowed choices are
  /// hardcoded using files in the top level of this project
  pub fn read_matrix(matrix_name: &str) -> ScoreMatrix {
    let blosum50 = include_str!("../../data/BLOSUM50.txt");
    let blosum62 = include_str!("../../data/BLOSUM62.txt");
    match matrix_name {
      "BLOSUM50" => read_matrix_impl(blosum50),
      "BLOSUM62" => read_matrix_impl(blosum62),
      _ => panic!("Unrecognized matrix name {}", matrix_name)
    }
  }

  fn read_matrix_impl(matrix_str: &str) -> ScoreMatrix {
    let mut scores = ScoreMatrix::new();
    let mut matrix_cols = vec![0 as u8; 24];
    for line in matrix_str.lines() {
      if line.starts_with("#") {
        continue;
      } else if line.starts_with(" ") {
        let fields = line.trim().split_whitespace();
        for (i, aa) in fields.enumerate() {
          matrix_cols[i] = aa.bytes().next().unwrap();
        }
      } else {
        let mut fields = line.trim().split_whitespace();
        let aa = fields.next().unwrap().bytes().next().unwrap();
        for i in 0..24 {
          let bb = matrix_cols[i];
          let score: i32 = fields.next().unwrap().parse().expect(
            &(format!("Can't parse field {} for amino acid '{}'", i+1, aa)));
          scores.insert((aa, bb), score);
        }
      }
    }
    scores
  }
}


/// Data structure for storing pairwise alignment strings
pub struct Alignment {
  pub seq1aligned: String,
  pub seq2aligned: String,
}

impl Alignment {

  /// Length of the alignment including all gaps
  pub fn length(&self) -> usize {
    self.seq1aligned.len()
  }

  /// Calculate sequence identity, including insertions and deletions
  pub fn identity(&self) -> f64 {
    let (mut n_mismatch, mut n_del, mut n_ins) = (0.0, 0.0, 0.0);
    for (aa, bb) in self.seq1aligned.chars().zip(self.seq2aligned.chars()) {
      if aa == '-' {
        n_ins += 1.0;
      } else if bb == '-' {
        n_del += 1.0;
      } else if aa != bb {
        n_mismatch += 1.0;
      }
    }
    1.0 - (n_mismatch + n_del + n_ins) / (self.length() as f64)
  }

  /// Create a string representing matching amino acids as asterisks
  fn get_seq_matches(&self) -> String {
    let mut match_chars = Vec::with_capacity(self.length());
    for (aa, bb) in self.seq1aligned.chars().zip(self.seq2aligned.chars()) {
      if aa == bb {
        match_chars.push(b'*');
      } else {
        match_chars.push(b' ');
      }
    }
    String::from_utf8(match_chars).unwrap()
  }

  pub fn format(&self, seqid1: &str, seqid2: &str) -> String {
    const SEQ_WIDTH: usize = 50;
    let seqid_len = cmp::max(seqid1.len(), seqid2.len());
    let seqid1_padded = format!("{}{}", seqid1, String::from_utf8(vec![b' '; seqid_len-seqid1.len()]).unwrap());
    let seqid2_padded = format!("{}{}", seqid2, String::from_utf8(vec![b' '; seqid_len-seqid2.len()]).unwrap());
    let seqid_matches = String::from_utf8(vec![b' '; seqid_len]).unwrap();
    let mut rows = Vec::with_capacity(3);
    let mut k: usize = 0;
    let matches = self.get_seq_matches();
    while k < self.length() {
      let end = cmp::min(self.length(), k + SEQ_WIDTH);
      let s1 = &self.seq1aligned[k..end];
      let m = &matches[k..end];
      let s2 = &self.seq2aligned[k..end];
      let line1 = format!("{}  {}", seqid1_padded, s1);
      let line2 = format!("{}  {}", seqid_matches, m);
      let line3 = format!("{}  {}", seqid2_padded, s2);
      rows.push([line1,line2,line3].join("\n"));
      k += SEQ_WIDTH;
    }
    rows.join("\n\n\n")
  }
}

struct MatrixCell {
  score: i32,
  cell_ptr: i8,
}

/// Needleman-Wunsch global pairwise protein sequence alignment algorithm,
/// as described in Durbin et al. chapter 2
pub fn needleman_wunsch (seq1: &str,
                         seq2: &str,
                         matrix_name: &str,
                         gap_penalty: i32) -> Alignment {
  const CELL_DIAG: i8 = 0;
  const CELL_LEFT: i8 = -1;
  const CELL_UP: i8 = 1;
  const GAP_CHAR: u8 = b'-';
  let gp = gap_penalty.abs(); // always a positive number, but subtracted
  if seq1.len() == 0 || seq2.len() == 0 {
    panic!("One or both sequences is empty.");
  }
  let score_matrix = blosum::read_matrix(matrix_name);
  let max_x = seq1.len() + 1;
  let max_y = seq2.len() + 1;
  let mut f = Vec::with_capacity(max_x);
  for i in 0..max_x {
    f.push(Vec::with_capacity(max_y));
    f[i].push(MatrixCell{score: -(i as i32) * gp, cell_ptr: CELL_LEFT});
  }
  for j in 1..max_y {
    f[0].push(MatrixCell{score: -(j as i32) * gp, cell_ptr: CELL_UP})
  }
  let seq1b = seq1.as_bytes();
  let seq2b = seq2.as_bytes();
  for i in 1..max_x {
    let aa = seq1b[i-1];
    for j in 1..max_y {
      let bb = seq2b[j-1];
      let score_diag = f[i-1][j-1].score + score_matrix.get(&(aa, bb)).unwrap();
      let score_left = f[i-1][j].score - gp;
      let score_up = f[i][j-1].score - gp;
      let cell = if (score_diag >= score_left) && (score_diag >= score_up) {
        MatrixCell{score: score_diag, cell_ptr: CELL_DIAG}
      } else if score_left >= score_up {
        MatrixCell{score: score_left, cell_ptr: CELL_LEFT}
      } else {
        MatrixCell{score: score_up, cell_ptr: CELL_UP}
      };
      f[i].push(cell);
    }
  }
  let mut tmp1aligned = Vec::with_capacity(max_x+max_y);
  let mut tmp2aligned = Vec::with_capacity(max_x+max_y);
  let mut i = max_x - 1;
  let mut j = max_y - 1;
  loop {
    if f[i][j].cell_ptr == CELL_DIAG {
      tmp1aligned.push(seq1b[i-1]);
      tmp2aligned.push(seq2b[j-1]);
      i -= 1;
      j -= 1;
    } else if f[i][j].cell_ptr == CELL_LEFT {
      tmp1aligned.push(seq1b[i-1]);
      tmp2aligned.push(GAP_CHAR);
      i -= 1;
    } else {
      tmp1aligned.push(GAP_CHAR);
      tmp2aligned.push(seq2b[j-1]);
      j -= 1;
    }
    if (i == 0) && (j == 0) {
      break;
    }
  }
  tmp1aligned.reverse();
  tmp2aligned.reverse();
  Alignment{seq1aligned: String::from_utf8(tmp1aligned).unwrap(),
            seq2aligned: String::from_utf8(tmp2aligned).unwrap()}
}

/*
  macro_rules! needleman_wunsch {
    ($seq1: expr, $seq2: expr, $matrix_name: expr, $gap_penalty: expr) => { needleman_wunsch($seq1, $seq2, $matrix_name, $gap_penalty) };
    ($seq1: expr, $seq2: expr, $matrix_name: expr) => { needleman_wunsch($seq1, $seq2, $matrix_name, 8) };
    ($seq1: expr, $seq2: expr) => { needleman_wunsch($seq1, $seq2, "BLOSUM50", 8) }
  }
*/


#[cfg(test)]
mod tests {
  use super::blosum;
  use super::needleman_wunsch;

  #[test]
  fn test_blosum() {
    let scores = blosum::read_matrix("BLOSUM50");
    assert_eq!(scores.get(&('A' as u8, 'A' as u8)).unwrap(), &5i32);
    assert_eq!(scores.get(&('W' as u8, 'W' as u8)).unwrap(), &15i32);
    assert_eq!(scores.get(&('C' as u8, 'W' as u8)).unwrap(), &-5i32);
    assert_eq!(scores.get(&('W' as u8, 'C' as u8)).unwrap(), &-5i32);
    assert_eq!(scores.get(&('*' as u8, '*' as u8)).unwrap(), &1i32);
    let scores2 = blosum::read_matrix("BLOSUM62");
    assert_eq!(scores2.get(&('A' as u8, 'A' as u8)).unwrap(), &4i32);
    assert_eq!(scores2.get(&('W' as u8, 'W' as u8)).unwrap(), &11i32);
    assert_eq!(scores2.get(&('C' as u8, 'W' as u8)).unwrap(), &-2i32);
    assert_eq!(scores2.get(&('W' as u8, 'C' as u8)).unwrap(), &-2i32);
    assert_eq!(scores2.get(&('*' as u8, '*' as u8)).unwrap(), &1i32);
  }

  #[test]
  #[should_panic]
  fn test_blosum_unknown_matrix() {
    let scores = blosum::read_matrix("foo");
  }

  #[test]
  fn test_nwalign_simple() {
    let a = needleman_wunsch("HEAGAWGHEE", "PAWHEAE", "BLOSUM50", 8);
    assert_eq!(&(a.seq1aligned), "HEAGAWGHE-E");
    assert_eq!(&(a.seq2aligned), "--P-AW-HEAE");
    assert!((a.identity() - 0.4545).abs() < 0.0001);
  }

  #[test]
  fn test_nwalign_negative_penalty() {
    let a = needleman_wunsch("HEAGAWGHEE", "PAWHEAE", "BLOSUM50", -8);
    assert_eq!(&(a.seq1aligned), "HEAGAWGHE-E");
    assert_eq!(&(a.seq2aligned), "--P-AW-HEAE");
    assert!((a.identity() - 0.4545).abs() < 0.0001);
  }

  #[test]
  #[should_panic]
  fn test_nwalign_empty_seq() {
    let a = needleman_wunsch("HEAGAWGHEE", "", "BLOSUM50", 8);
  }

  #[test]
  #[should_panic]
  fn test_nwalign_unknown_matrix() {
    let a = needleman_wunsch("HEAGAWGHEE", "PAWHEAE", "foo", 8);
  }

/*
  #[test]
  fn test_nwalign_macros() {
    let a = needleman_wunsch!("HEAGAWGHEE", "PAWHEAE");
    assert_eq!(&(a.seq1aligned), "HEAGAWGHE-E");
    assert_eq!(&(a.seq2aligned), "--P-AW-HEAE");
  }
*/

  #[test]
  fn test_nwalign_hemoglobin() {
    let seq1 = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    let seq2 = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
    let aln1 = "V-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    let aln2 = "VHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
    let a = needleman_wunsch(seq1, seq2, "BLOSUM50", 8);
    assert_eq!(&(a.seq1aligned), aln1);
    assert_eq!(&(a.seq2aligned), aln2);
    assert!((a.identity() - 0.4324).abs() < 0.0001);
  }
}
