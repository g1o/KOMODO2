use strict;
use warnings;

open (IN, "$ARGV[0]");

#this program takes a multi-proteome-per-line file (since uniprot is 
#gene/protein-centered) and splits it into a single-proteome-per-line
#file. Remember to change $prot_id to point to the right collumn!

my $index = "11";

while (my $line = <IN>) {
  chomp $line;
  my @aux;
  @aux = split (/\t/, $line);
#  print "$#aux\n";
#  my $a = <STDIN>;
  my $prot_id = $aux[$index];
  if (defined $prot_id) {
    my $count = ($prot_id =~ tr/://);
#    print "$prot_id\t$count\n";
#    my $a = <STDIN>;
    if ($count > 1) {
#      print "$prot_id\n";
#      my $a = <STDIN>;
#    print "Here!\n";
      my @aux2 = split (/,\s+UP/, $prot_id);
#      print "\n\n\n@aux2\n\n\n";
      for (my $i = 0; $i <= $#aux2; $i++) {
        if ($i == 0) { #No UP missing at the beginning
          my $element = $aux2[$i];
          splice(@aux, $index, 1, $element);
          my $newline = join ("\t", @aux);
#          $newline =~ s/\t+$//;
          print "$newline\n";
        } else {# UP is missing
          my $element = $aux2[$i];
          $element = "UP".$element;
#          print "\n\n\n$element\n\n\n";
#          my $a = <STDIN>;
          splice(@aux, $index, 1, $element);
          my $newline = join ("\t", @aux);
#          $newline =~ s/\t+$//;
          print "$newline\n";
        }
#        print "$line\n";
#        my $a = <STDIN>;
      }
    }
    else {
      my $newline = join ("\t", @aux);
#    $newline =~ s/\t+$//;
      print $newline."\n";
    }
  } else { #if line does not contain a proteome ID, remove it
#    print "Line $line not defined\n";
  }
}

close IN;
