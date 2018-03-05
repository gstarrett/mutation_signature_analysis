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

open(EARLY, "> $somatic.early");
open(LATE, "> $somatic.late");

while (my $line = <SOM>) {
	next if $line =~ /^Hugo/;
	my $found = 0;
	my @f = split("\t", $line);
	#open(TMP, "> tmp.bed");
	next unless defined $f[4];
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
		my $lcutoff = (1/$s[3])+((1/$s[3])*0.2);
		my $ucutoff = (1/$s[3])-((1/$s[3])*0.2);
		my $af = 0;
		if ($type eq "BRCA") {
			$af = $f[63]/100;
		} elsif ($type eq "LUAD") {
			$af = $f[80]/($f[80]+$f[81])
		} elsif ($type eq "SKCM") {
			$af = $f[83]/($f[80]+$f[83])
		} elsif ($type eq "LIHC") {
			$af = $f[43]/($f[42])
		} elsif ($type eq "LUSC" || $type eq "PRAD" || $type eq "HNSC" || $type eq "THCA") {
			$af = $f[79]/($f[80] + $f[79])
		} elsif ($type eq "LGG") {
			$af = $f[36]/($f[37] + $f[36])
		} elsif ($type eq "BLCA") {
			$af = $f[84]/($f[85] + $f[84])
		} elsif ($type eq "UCEC") {
			$f[42] =~ s/%//;
			$af = $f[42]/100;
		} else {
			die "$!: Invalid cancer type provided.\n\n";
		}

		if ($af < $lcutoff) {
			print LATE $line;
		} elsif ($af > $ucutoff) {
			print EARLY $line;
		}
	} elsif ($s[3] < 2 && $s[3] > 1.8) {
		my $cutoff = (1/$s[3])-((1/$s[3])*0.2);
		my $af = 0;
		if ($type eq "BRCA") {
			$af = $f[63]/100;
		} elsif ($type eq "LUAD") {
			$af = $f[80]/($f[80]+$f[81])
		} elsif ($type eq "SKCM") {
			$af = $f[83]/($f[80]+$f[83])
		} elsif ($type eq "LIHC") {
			$af = $f[43]/($f[42])
		} elsif ($type eq "LUSC" || $type eq "PRAD" || $type eq "HNSC" || $type eq "THCA") {
			$af = $f[79]/($f[80] + $f[79])
		} elsif ($type eq "LGG") {
			$af = $f[36]/($f[37] + $f[36])
		} elsif ($type eq "BLCA") {
			$af = $f[84]/($f[85] + $f[84])
		} elsif ($type eq "UCEC") {
			$f[42] =~ s/%//;
			$af = $f[42]/100;
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