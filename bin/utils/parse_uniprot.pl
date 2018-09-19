use strict;
use warnings;

#takes as input a sorted-by-proteome file and splits it into individual
#proteome files.

open (IN, "$ARGV[0]");

my $header = <IN>;

my $firstline = <IN>;

my @tmp = split(/\t/, $firstline);

my $index = 10;

my $previous_id = $tmp[$index];

#@tmp = "";

#@tmp = split (/:\s+/, $previous_id);

#$previous_id = $tmp[0];

#$previous_id =~ s/\s+/_/g;

open (OUT, ">$previous_id.txt");

print OUT $header;

print OUT $firstline;

my $i = 0;
while (my $line = <IN>) {
  chomp $line;
  my @aux = split (/\t/, $line);
  my $id = $aux[$index];
#  @tmp = "";
#  @tmp = split (/:\s+/, $id);
#  $id = $tmp[0];
#  $id =~ s/\s+/_/g;
  if ($id eq $previous_id) {
    print OUT "$line\n";
    $previous_id = $id;
  }
  else {
    close OUT;
    open (OUT, ">$id.txt");
    print OUT "$header";
    print OUT "$line\n";
    $previous_id = $id;
    $i++;
    print "$i\n";
  }
}

close IN;

#print "Creating files\n";

#$i = 0;

#foreach my $key (keys %data) {
#  open (OUT, ">$key.txt");
#  my @aux = split (/\*\*/, $data{$key});
#  foreach my $element (@aux) {
#    print OUT "$element\n";
#  }
#  print "$i\n";
#  $i++;
#  close OUT;
#}

