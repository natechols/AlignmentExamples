// Author: Nat Echols
// License: public domain

use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::Path;

/// Representation of a FASTA file record, with header split into ID and
/// description
pub struct FastaRecord {
  pub id: String,
  pub sequence: String,
  pub description: Option<String>,
}

impl FastaRecord {

  /// Return the full record header including description, if any
  pub fn header(&self) -> String {
    match self.description {
      Some(ref desc) => format!(">{} {}", self.id, desc),
      None => format!(">{}", self.id)
    }
  }

  // TODO
  //pub fn format(&self, width: u32) -> String {
  //}
}

fn new_record(header: String, sequence_lines: Vec<String>) -> FastaRecord {
  if sequence_lines.len() == 0 {
    panic!("Empty sequence for {}", header);
  } else if header.len() == 0 {
    panic!("Empty header line not allowed");
  }
  let mut header_fields = header.split_whitespace();
  let id = header_fields.next().unwrap().to_string();
  let mut desc_fields = Vec::with_capacity(1);
  loop {
    let field = header_fields.next();
    match field {
      Some(s) => desc_fields.push(s),
      None => break
    };
  }
  let desc = if desc_fields.len() > 0 {
    Some(desc_fields.join(" "))
  } else {
    None
  };
  FastaRecord{id: id, sequence: sequence_lines.join(""), description: desc}
}

/// Read a FASTA file line by line and return a Vec of FastaRecord objects
pub fn read_fasta(file_name: &str) -> Vec<FastaRecord> {
  let file = match File::open(Path::new(file_name)) {
    Err(_) => panic!("Couldn't open {}", file_name),
    Ok(f) => f
  };
  let mut records = Vec::with_capacity(1);
  let mut current_header: Option<String> = None;
  let mut current_sequence: Vec<String> = Vec::with_capacity(1);
  for line in BufReader::new(file).lines() {
    let l = line.unwrap();
    if l.starts_with(">") {
      if current_header.is_some() {
        records.push(new_record(current_header.unwrap(), current_sequence));
      }
      current_header = Some(l.trim_left_matches('>').trim().to_string());
      current_sequence = Vec::with_capacity(1);
    } else if current_header.is_some() {
      current_sequence.push(l.trim().to_string());
    } else {
      panic!("Illegal FASTA line without header:\n{}", l);
    }
  }
  if current_header.is_none() {
    panic!("No valid FASTA records in file");
  }
  records.push(new_record(current_header.unwrap(), current_sequence));
  records
}


#[cfg(test)]
mod tests {
  //use tempdir::TempDir;
  use super::read_fasta;
  use std::io::Write;
  use std::fs;

  #[test]
  fn test_fasta() {
    // FIXME figure out how to use temp files properly
    //let tmp_dir = TempDir::new("test-fasta");
    //let tmp_file = tmp_dir.path().join("test1.fasta");
    let mut f = fs::File::create("test1.fasta");
    writeln!(f.unwrap(), ">seq1\nHEAGAWGHEE\n>seq2 my sequence\nPAWHEAE");
    let records = read_fasta("test1.fasta");
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "seq1".to_string());
    assert_eq!(records[0].sequence, "HEAGAWGHEE".to_string());
    assert!(records[0].description.is_none());
    assert_eq!(records[0].header(), ">seq1".to_string());
    assert_eq!(records[1].id, "seq2".to_string());
    assert_eq!(records[1].sequence, "PAWHEAE".to_string());
    assert!(records[1].description.is_some());
    assert_eq!(records[1].header(), ">seq2 my sequence".to_string());
    fs::remove_file("test1.fasta");
  }

  // TODO test failure modes
}
