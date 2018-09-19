use strict;
use warnings;

my %flag; #will check for unique IDs

my $indir = $ARGV[0];

opendir(DIR, $indir);

my @files = readdir(DIR);

closedir(DIR);

foreach my $file (@files) {
  next if (($file eq ".")||($file eq ".."));
  print "Analyzing file $file\n";
  my $infile_path = $indir.$file;
  my $outfile_path = $file.".tmp";
  open (IN, $infile_path);
  open (OUT, ">$outfile_path");
  my $firstline = <IN>;
  chomp $firstline;
  my @aux = split(/\t/, $firstline);
#  $firstline = join ("\t", $aux[0], $aux[7]);
  $firstline = join ("\t", $aux[0], $aux[7], $aux[9], $aux[15], $aux[16]);
  print OUT $firstline."\n";
  while (my $line = <IN>) {
    chomp $line;
    my @aux = split(/\t/, $line);
    if (defined $flag{$aux[0]}) {
      print ("Duplicated ID found: $aux[0]\n");
      $flag{$aux[0]}++;
      $aux[0] = join("_",$aux[0],$flag{$aux[0]});
#      $line = join ("\t", $aux[0], $aux[7]);
      $line = join ("\t", $aux[0], $aux[7], $aux[9], $aux[15], $aux[16]);
      print OUT $line."\n";
    } else {
#      $line = join ("\t", $aux[0], $aux[7]);
      $line = join ("\t", $aux[0], $aux[7], $aux[9], $aux[15], $aux[16]);
      print OUT $line."\n";
      $flag{$aux[0]} = 1;
    }
  }
  close IN;
  close OUT;
}
