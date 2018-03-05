#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Fasta;

my $n = 0;
my $fasta = shift;
my $name = shift;
my $Ctrinucs;
my $Gtrinucs;

# print "creating fasta db... ";
my $db = Bio::DB::Fasta->new($fasta);

my @chrs = $db->ids;
#print Dumper(@chrs);

for my $chr (@chrs) {
	my $seqlen = $db->length($chr);
	for (my $pos = 2; $pos < $seqlen; $pos++) {
		my $ref = uc($db->seq($chr,$pos => $pos));
		#print $ref,"\n";
		if ($ref eq "C") {
			my $seq = uc($db->seq($chr, $pos-2+$n => $pos+1+$n));
			my @bases = split('', $seq);
			$seq = $bases[0] . "C" . $bases[2] . $bases[3];
			$Ctrinucs .= "$seq\n";
			print "$seq\n";
			if ($seq =~ /[AaTtCcGg][Cc][Gg]/) {
				next;
			} else {
				$Ctrinucs .= "$seq\n";
			}
		} elsif ($ref eq "G" ) {
			my $seq = uc($db->seq($chr, $pos-1+$n => $pos+2+$n));
			my @bases = split('', $seq);
			$seq = $bases[3] . $bases[2] . "G" . $bases[0];
			$Gtrinucs .= "$seq\n";
	
			if ($seq =~ /[Cc][Gg][AaTtCcGg]/) {
				next;
			} else {
				$Gtrinucs .= "$seq\n";
			}
		} else { next }
	}
}

open (OUTC, "> $name.C_trinucs_count.txt");
open (OUTG, "> $name.G_trinucs_count.txt");

print OUTC $Ctrinucs;
print OUTG $Gtrinucs;

close(OUTC);
close(OUTG);