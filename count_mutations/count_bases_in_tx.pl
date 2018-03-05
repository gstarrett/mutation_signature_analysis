#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Fasta;

my $fasta = shift;
my $in = shift;

print "creating fasta db... ";
my $db = Bio::DB::Fasta->new($fasta);
print "complete!\n";


open (IN, "< $in");
while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t",$line);
	my $A=0;
	my $C=0;
	my $G=0;
	my $T=0;
	my $region = uc($db->seq($f[0], $f[1] => $f[2]));
	while ($region =~ /A/g) { $A++ }
	while ($region =~ /C/g) { $C++ }
	while ($region =~ /G/g) { $G++ }
	while ($region =~ /T/g) { $T++ }
	print join("\t",$line,$A,$C,$G,$T), "\n";
}


# open (OUT, "> $name.count.txt");
# print OUT $out;
# close(OUT);


exit;