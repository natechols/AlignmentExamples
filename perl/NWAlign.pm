
package NWAlign;

use strict;
use warnings;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 0.1;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(read_matrix needleman_wunsch);

sub read_matrix {
  my ($matrix_name) = @_;
  my $file_name = $ENV{'BLOSUM_DATA'} . "/${matrix_name}.txt";
  open(my $fh, "<", $file_name) or die "Couldn't open $file_name";
  my %matrix_cols;
  my %score_matrix;
  while (my $line = <$fh>) {
    if ($line =~ /^ /) {
      $line =~ s/^\ *//;
      my @fields = split(/\s+/, $line);
      for (my $i = 0; $i < scalar @fields; $i++) {
        $matrix_cols{$i} = $fields[$i];
      }
    } elsif ($line !~ /^#/) {
      my @fields = split(/\s+/, $line);
      my $aa = $fields[0];
      for (my $i = 1; $i < scalar @fields; $i++) {
        my $bb = $matrix_cols{$i - 1};
        $score_matrix{$aa}{$bb} = $fields[$i] + 0;
      }
    }
  }
  return %score_matrix;
}

sub needleman_wunsch {
  my ($seq1, $seq2, $score_matrix, $gap_penalty) = @_;
  my ($CELL_DIAG, $CELL_LEFT, $CELL_UP) = (0, -1, 1);
  my $maxX = length($seq1) + 1;
  my $maxY = length($seq2) + 1;
  my @seq1a = split(//, $seq1);
  my @seq2a = split(//, $seq2);
  my @f;
  for (my $i = 0; $i < $maxX; $i++) {
    $f[$i][0] = {'score' => - $i * $gap_penalty, 'cell_ptr' => $CELL_LEFT};
  }
  for (my $j = 0; $j < $maxY; $j++) {
    $f[0][$j] = {'score' => - $j * $gap_penalty, 'cell_ptr' => $CELL_UP};
  }
  for (my $i = 1; $i < $maxX; $i++) {
    my $aa = $seq1a[$i - 1];
    for (my $j = 1; $j < $maxY; $j++) {
      my $bb = $seq2a[$j - 1];
      my $score_diag = $f[$i-1][$j-1]{'score'} + $$score_matrix{$aa}{$bb};
      my $score_left = $f[$i-1][$j]{'score'} - $gap_penalty;
      my $score_up = $f[$i][$j-1]{'score'} - $gap_penalty;
      if (($score_diag >= $score_left) and ($score_diag >= $score_up)) {
        $f[$i][$j] = {'score' => $score_diag, 'cell_ptr' => $CELL_DIAG};
      } elsif ($score_left >= $score_up) {
        $f[$i][$j] = {'score' => $score_left, 'cell_ptr' => $CELL_LEFT}
      } else {
        $f[$i][$j] = {'score' => $score_up, 'cell_ptr' => $CELL_UP};
      }
    }
  }
  my @tmp1;
  my @tmp2;
  my $i = length($seq1);
  my $j = length($seq2);
  while (($i >= 0) and ($j >= 0) and !(($i == 0) and ($j == 0))) {
    my %cell = %{$f[$i][$j]};
    if ($cell{'cell_ptr'} == $CELL_DIAG) {
      push(@tmp1, $seq1a[$i - 1]);
      push(@tmp2, $seq2a[$j - 1]);
      $i--;
      $j--;
    } elsif ($cell{'cell_ptr'} == $CELL_LEFT) {
      push(@tmp1, $seq1a[$i - 1]);
      push(@tmp2, '-');
      $i--;
    } else {
      push(@tmp1, '-');
      push(@tmp2, $seq2a[$j - 1]);
      $j--;
    }
  }
  my $seq1aligned = join('', reverse(@tmp1));
  my $seq2aligned = join('', reverse(@tmp2));
  return (
    seq1aligned => $seq1aligned,
    seq2aligned => $seq2aligned,
  );
}
