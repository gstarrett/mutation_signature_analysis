#!/usr/bin/perl -w 
use strict;

my $in = shift;

open(IN, "< $in");
my $n = 0;
my @pclust = ();
my $imdSum = 0;
my $head = <IN>;
my $prev;
while (<IN>) {
	chomp($_);
	my @f = split("\t", $_);
	$imdSum += $f[3];
	if ($imdSum/(scalar(@pclust)+1) <= 5000) {
		if (scalar(@pclust)==0) {
			push(@pclust, $prev);
			$imdSum += 5000;
		}
		push(@pclust, $_);
	} elsif ($imdSum/(scalar(@pclust)+1) > 5000 && scalar(@pclust) > 5) {
		for my $line (@pclust) { 
			print join("\t",$line,$n,scalar(@pclust)), "\n";
		}
		$n++;
		@pclust = ();
		$imdSum = 0;
	} else {
		@pclust = ();
		$imdSum = 0;
	}
		$prev = $_;
}