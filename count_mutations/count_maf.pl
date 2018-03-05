#!/usr/bin/perl -w
use strict;

my $infile = shift;

open(FILE, "< $infile");

my %C;
my %C2T;
my %total;

while(my $line = <FILE>) {
	next if $line =~ /^#/ || $line !~ /SNP/;
	chomp($line);
	my @f = split("\t", $line);
	my $sample = $f[15];
	if (exists $total{$sample}) {
		$total{$sample}++;
	} else {
		$total{$sample} = 1;
	}
	if ($f[17] eq "C" && $f[18] eq "C") {
		if (exists $C{$sample}) {
			$C{$sample}++;
		} else {
			$C{$sample} = 1;
		}
		if ($f[12] eq "T") {
			if (exists $C2T{$sample}) {
				$C2T{$sample}++;
			} else {
				$C2T{$sample} = 1;
			}
		}
	}
}

close(FILE);

for my $key (keys %total) {
	print "$key\t$total{$key}\ttotal\n";
}

for my $key (keys %C) {
	print "$key\t$C{$key}\tC\n";
}

for my $key (keys %C2T) {
	print "$key\t$C2T{$key}\tC2T\n";
}