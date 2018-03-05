#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Fasta;

my $n = 0;
my $fasta = shift;
my $variants = shift;
my $name = shift;
my $Ctrinucs;
my $Gtrinucs;

my %trinucHash;

my %normalization = {
	ACA =>	0.07,
	ACC =>	0.06,
	ACG =>	0.02,
	ACT =>	0.06,
	CCA =>	0.10,
	CCC =>	0.09,
	CCG =>	0.02,
	CCT =>	0.09,
	GCA =>	0.07,
	GCC =>	0.08,
	GCG =>	0.02,
	GCT =>	0.07,
	TCA =>	0.07,
	TCC =>	0.08,
	TCG =>	0.01,
	TCT =>	0.08
};

print "creating fasta db... ";
my $db = Bio::DB::Fasta->new($fasta);
print "complete!\n";

open (VAR, "< $variants");
while (my $line = <VAR>) {
	my @fields = split("\t",$line);
	my $chr = shift @fields;
	my $pos = shift @fields;
	my $ref = shift @fields;
	my $alt = shift @fields;
	#print "$chr $pos $ref $alt\n";
	if ($ref eq "C") {
		my $seq = uc($db->seq($chr, $pos-1+$n => $pos+1+$n));
		my @bases = split('', $seq);
		$seq = $bases[0] . "C" . $bases[2];
		$Ctrinucs .= "$seq\t$ref>$alt\n";
		
# 		if ($seq =~ /[AaTtCcGg][Cc][Gg]/) {
# 			next;
# 		} else {
# 			$Ctrinucs .= "$seq\n";
# 		}
	} elsif ($ref eq "G" ) {
		my $seq = uc($db->seq($chr, $pos-1+$n => $pos+1+$n));
		my @bases = split('', $seq);
		$seq = $bases[2] . "G" . $bases[0];
		$Gtrinucs .= "$seq\t$ref>$alt\n";
		
# 		if ($seq =~ /[Cc][Gg][AaTtCcGg]/) {
# 			next;
# 		} else {
# 			$Gtrinucs .= "$seq\n";
# 		}
	} else { next }
}
open (OUTC, "> $name.C_trinucs_count.txt");
open (OUTG, "> $name.G_trinucs_count.txt");

print OUTC $Ctrinucs;
print OUTG $Gtrinucs;

close(OUTC);
close(OUTG);