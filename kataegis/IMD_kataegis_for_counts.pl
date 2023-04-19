#!/usr/bin/perl -w
use strict;
no warnings 'uninitialized';

my $in = shift;

my %dataHash;
# maybe put in a system sort before this step;
open (IN, "< $in");
#open (OUT, "> $out.imd.txt");

open (OUT, "> $in.imd.txt");
while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);
	my $chr = $f[0];
	my @data = @f[1,5];
	if (exists $dataHash{$chr}) {
		push(@{$dataHash{$chr}}, \@data);
	} else {
		$dataHash{$chr} = [\@data];
	}
}

print OUT "#CHR\tPOS\tMUT\tPROX_IMD\tDIST_IMD\tTYPE\n";

for my $key (keys %dataHash) {
	my $max = scalar @{$dataHash{$key}};
	for my $idx (1 .. $max-2) {
		my $prev = ${$dataHash{$key}}[$idx-1];
		my $dra = ${$dataHash{$key}}[$idx];
		my $next = ${$dataHash{$key}}[$idx+1];
		my $prox = ${$dra}[0] - ${$prev}[0];
		my $dist = ${$next}[0] - ${$dra}[0];

		my $cprev = "";
		my $cprox = "";
		my $cnext = "";
		my $cdist = "";

# 		if (${$dra}[1] eq "C") {
# 			my $n = 1;
# 			until (${$dataHash{$key}}[$idx-$n] eq "C") {
# 				$n++;
# 			}
# 			if (${$dataHash{$key}}[$idx-$n] eq "C") {
# 				$cprev = ${$dataHash{$key}}[$idx-$n];
# 				$cprox = ${$dra}[0] - ${$cprev}[0];
# 			}
#
# 			$n = 1;
# 			until (${$dataHash{$key}}[$idx+$n] eq "C") {
# 				$n++;
# 			}
# 			if (${$dataHash{$key}}[$idx-$n] eq "C") {
# 				$cnext = ${$dataHash{$key}}[$idx+$n];
# 				$cdist = ${$cnext}[0] - ${$dra}[0];
# 			}
# 		} elsif (${$dra}[1] eq "G") {
# 			my $n = 1;
# 			until (${$dataHash{$key}}[$idx-$n] eq "G") {
# 				$n++;
# 			}
# 			if (${$dataHash{$key}}[$idx-$n] eq "G") {
# 				$cprev = ${$dataHash{$key}}[$idx-$n];
# 				$cprox = ${$dra}[0] - ${$cprev}[0];
# 			}
#
# 			$n = 1;
# 			until (${$dataHash{$key}}[$idx+$n] eq "G") {
# 				$n++;
# 			}
# 			if (${$dataHash{$key}}[$idx-$n] eq "G") {
# 				$cnext = ${$dataHash{$key}}[$idx+$n];
# 				$cdist = ${$cnext}[0] - ${$dra}[0];
# 			}
# 		}
		my $type;
		my $mut = ${$dra}[1];
		if ($mut eq "C>T" || $mut eq "T>C" || $mut eq "G>A" || $mut eq "A>G") {
			$type = "transition";
		} else {
			$type = "transversion";
		}
		print OUT "$key\t${$dra}[0]\t$mut\t$prox\t$dist\t$type\n";
	}
}

close(IN);
close(OUT);
