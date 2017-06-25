
import gzip


class FastaRecord(object):

    def __init__(self, id_, sequence, description=None):
        self.id = id_
        self.sequence = sequence
        self.description = description

    @property
    def header(self):
        if not self.description:
            return "".join([">", self.id])
        else:
            return "".join([">", self.id, " ", self.description])

    def format(self, width=60):
        lines = [self.header]
        k = 0
        while k < len(self.sequence):
            lines.append(self.sequence[k:(k + width)])
            k += width
        return "\n".join(lines)

    def __str__(self):
        return self.format()

    def __getitem__(self, idx):
        return self.sequence[idx]

    def __len__(self):
        return len(self.sequence)


def read_fasta(file_name):
    _open_file = lambda f: open(f, "rU")
    if file_name.endswith(".gz"):
        _open_file = lambda f: gzip.open(f, "r")
    records = []

    def new_record(h, s):
        fields = h.lstrip(">").strip().split()
        if len(fields) == 0:
            raise IOError("Cannot parse as FASTA header:\n{h}".format(h=l))
        elif len(s) == 0:
            raise IOError("Empty sequence for FASTA record:\n{h}".format(h=l))
        return FastaRecord(fields[0], "".join(s), " ".join(fields[1:]))
    current_header = None
    current_sequence = None
    with _open_file(file_name) as f:
        for line in f:
            if line.strip() == "":
                raise IOError("Empty lines not allowed")
            elif line.startswith(">"):
                if current_header is not None:
                    records.append(
                        new_record(current_header, current_sequence))
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line.strip())
    if current_header is None or current_sequence is None:
        raise IOError("No valid FASTA records in {f}".format(f=file_name))
    records.append(new_record(current_header, current_sequence))
    return records
