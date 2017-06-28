// Author: Nat Echols
// License: public domain

package libnwalign

import (
  "os"
  "bufio"
  "bytes"
  "strings"
)

type FastaRecord struct {
  Id string
  Sequence string
  Description string
}

func makeRecord (header string, sequence string) FastaRecord {
  headerFields := strings.Split(header, " ")
  return FastaRecord{headerFields[0], sequence, strings.Join(headerFields[1:], " ")}
}

func FastaReader (fileName string) []FastaRecord {
  f, err := os.Open(fileName)
  if (err != nil) {
    panic("Can't read file " + fileName + " (error: " + err.Error() + ")")
  }
  defer f.Close()

  records := make([]FastaRecord, 0)
  scanner := bufio.NewScanner(f)
  var currentHeader string
  var currentSequence bytes.Buffer

  for scanner.Scan() {
    line := scanner.Text()
    if (strings.HasPrefix(line, ">")) {
      if (currentHeader != "") {
        records = append(records, makeRecord(currentHeader, currentSequence.String()))
      }
      currentHeader = strings.Trim(line, ">")
      currentSequence.Reset()
    } else {
      currentSequence.WriteString(line)
    }
  }
  if (currentHeader == "") {
    panic("No valid FASTA records found")
  }
  records = append(records, makeRecord(currentHeader, currentSequence.String()))
  return records
}
