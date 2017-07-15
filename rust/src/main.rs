
extern crate nwalign;

use nwalign::fasta;

use std::env;

fn main () {
  let args: Vec<_> = env::args().collect();
  let mut files = Vec::with_capacity(2);
  let mut matrix_name = "BLOSUM50".to_string();
  let mut gap_penalty = 8;
  let mut k = 1;
  while k < args.len() {
    if args[k] == "-m" || args[k] == "--matrix" {
      matrix_name = (&args[k+1]).to_string();
      k += 1;
    } else if args[k] == "-g" || args[k] == "--gap-penalty" {
      gap_penalty = args[k+1].parse().unwrap();
      k += 1
    } else{
      files.push((&args[k]).to_string());
    }
    k += 1
  }
  if files.len() != 2 {
    panic!("Exactly two files are required as input, got {}:\n{}",
           files.len(), files.join(", "));
  }
  let fasta1 = fasta::read_fasta(&files[0]);
  let fasta2 = fasta::read_fasta(&files[1]);
  let aln = nwalign::needleman_wunsch(&fasta1[0].sequence, &fasta2[0].sequence,
                                      &matrix_name, gap_penalty);
  let pct_id = format!("{:.2}", aln.identity()*100.0);
  println!("Sequence identity = {}%\n", pct_id);
  println!("{}", aln.format(&fasta1[0].id, &fasta2[0].id));
}
