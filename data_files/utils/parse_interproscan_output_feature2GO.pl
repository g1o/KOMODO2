use strict;
use warnings;

if (!$ARGV[2]) {
  print_help();
}

my $infile = $ARGV[0];

my $feature_name = $ARGV[1];

my $outdir = $ARGV[2];

chomp $outdir;

my @path = split(/\//, $infile);

my $file_name = pop @path;

my @tmp_name = split(/_/, $file_name);

my $final_name = join("_", $tmp_name[0], $tmp_name[1]);

$final_name =~ s/\./_/g;

my $outfile = join("/", $outdir, $final_name);

$outfile = join(".", $outfile, "$feature_name"."2GO.txt");

print ("Generating file $outfile\n");

open(IN, "<$infile") || die($!);

open(OUT, ">$outfile") || die ($!);

#my $header = <IN>;

#my @col_names = split(/\t/, $header);

print OUT ("Feature\t$feature_name\n");

#my %data;

#my $total = 0;

my $i = 0;

while (my $line = <IN>) {
  chomp $line;   
  my @aux = split(/\t/, $line);
  my $actual_feature = $aux[3];
  if ((defined $actual_feature)&&($actual_feature eq $feature_name)) {
    my $GOs = get_GO($line);
    if (defined $GOs) { #new feature has GO terms; print them
      $GOs =~ s/\|/;/g;
      print OUT "$feature_name\_$i\t$GOs\n";
      $i++;
    } else {
      next();
    }
  }
}

close IN;

close OUT;

sub print_help {
  die("Use this program like: perl parse_interproscan_feature2GO.pl <path to interproscan output file> <feature to extract (e.g. \"Pfam\")> <output dir>\n");
}

sub get_GO {
  my $tmp = $_[0];
  my @aux = split(/\t/, $tmp);
  foreach my $element (@aux) {
    if ($element =~ /^GO:/) {
      return $element;
    }
  }
  return undef;
}
