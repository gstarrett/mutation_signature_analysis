#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $mutsin =shift;
my $annin = shift;
my $other = shift;

my @possMuts = ("C>T","C>A","C>G","G>A","G>C","G>T","A>T","A>C","A>G","T>A","T>C","T>G");

my %normCount;
my %annHash;
my %otherHash;

open(ANN, "< $annin");
while (<ANN>) {
	my @f = split("\t",$_);
	$annHash{$f[0]} = $f[1];
}
close(ANN);

open(OTHER, "< $other");
while (<OTHER>) {
	my @f = split("\t",$_);
	$otherHash{$f[1]} = $f[0];
}
#print STDERR Dumper(%otherHash);
close(OTHER);

open(MUTS, "< $mutsin");
while(<MUTS>) {
	my @f = split("\t", $_);
	if (exists $annHash{$f[0]}) {
		my $normMut = 1/$annHash{$f[0]};
		my %mutHash;
		my $mut = "$f[11]>$f[12]";
		my $st = "NA";
		if (exists $otherHash{$f[0]}) {
			$st = $otherHash{$f[0]};
		}
		if (exists $normCount{"$f[0]\t$st"}) {
			if (exists ${$normCount{"$f[0]\t$st"}}{$mut}) {
				${$normCount{"$f[0]\t$st"}}{$mut} += $normMut;
			} else {
				${$normCount{"$f[0]\t$st"}}{$mut} += $normMut;
			}
		} else {
			$mutHash{$mut} = $normMut;
			$normCount{"$f[0]\t$st"} = \%mutHash;
		}
	} else {
		next;
	}
}

#print Dumper(%normCount);
print join("\t", "GENE", "STRAND", (sort @possMuts), "TOTAL"),"\n";
for my $key (sort keys %normCount){
	my @outl;
	my $sum;
	for (sort @possMuts) {
		if (exists ${$normCount{$key}}{$_}) {
			push(@outl, ${$normCount{$key}}{$_});
		} else {
			push(@outl, 0);
		}
	}
	$sum += $_ for @outl;
	print join("\t", $key, @outl, $sum),"\n"; 
}