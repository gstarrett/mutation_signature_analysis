#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $bedtools = "/Applications/bioinf_tools/bedtools2.23/bin/bedtools";

my $type = shift;
my $somatic = shift;
my $cnvDir = shift;

opendir(CNV, $cnvDir);
my @files = readdir(CNV);
closedir(CNV);

open(SOM, "< $somatic");

my $early;
my $late;

open(EARLY, "> early.maf.txt");
open(LATE, "> late.maf.txt");

while (my $line = <SOM>) {
	next if $line =~ /^Hugo/;
	my $found = 0;
	my @f = split("\t", $line);
	#open(TMP, "> tmp.bed");
	my $chr = "chr$f[4]";
	my $snpBed = join("\t", $chr, $f[5]-1, $f[5]);
	#print TMP $snpBed;
	#print $snpBed,"\n";
	#close(TMP);
	my $sample = substr($f[15],0,12);
	foreach (@files) {
		$found = 1 if $_ =~ /$sample/;
	}
	next unless $found == 1;
	my $segment = `echo "$snpBed" | $bedtools intersect -a $cnvDir/$sample.bed -b stdin | head -1`;
	next if length($segment) < 1;
	#print "match\n",$segment;
	my @s = split("\t", $segment);
	if ($s[3] > 2) {
		my $lcutoff = (1/$s[3])+((1/$s[3])*0.1);
		my $ucutoff = (1/$s[3])-((1/$s[3])*0.1);
		my $af = 0;
		if ($type eq "BRCA") {
			$af = $f[63]/100;
		} elsif ($type eq "LUAD") {
			$af = $f[80]/($f[80]+$f[81])
		} elsif ($type eq "SKCM") {
			$af = $f[83]/($f[80]+$f[83])
		} elsif ($type eq "LIHC") {
			$af = $f[43]/($f[42])
		} else {
			die "$!: Invalid cancer type provided.\n\n";
		}

		if ($af < $lcutoff) {
			print LATE $line;
		} elsif ($af > $ucutoff) {
			print EARLY $line;
		}
	} elsif ($s[3] < 2 && $s[3] > 1.8) {
		my $cutoff = (1/$s[3])-((1/$s[3])*0.1);
		my $af = 0;
		if ($type eq "BRCA") {
			$af = $f[63]/100;
		} elsif ($type eq "LUAD") {
			$af = $f[80]/($f[80]+$f[81])
		} elsif ($type eq "SKCM") {
			$af = $f[83]/($f[80]+$f[83])
		} elsif ($type eq "LIHC") {
			$af = $f[43]/($f[42])
		} else {
			die "$!: Invalid cancer type provided.\n\n";
		}

		if ($af < $cutoff) {
			print LATE $line;
		} elsif ($af > $cutoff) {
			print EARLY $line;
		}
	}
}

#print $early;

#print $late;