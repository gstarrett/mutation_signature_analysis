#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $bedtools = "/Applications/bioinf_tools/bedtools2.23/bin/bedtools";

my $type = shift;
my $somatic = shift;
my $cnvDir = shift;
my $estimate = shift;
my $pre = shift;

unless (defined $pre) {
	$pre = $somatic;
}

my %purity;

open(EST, "< $estimate");
my $head = <EST>;
while (my $line = <EST>) {
	chomp($line);
	my @f = split("\t", $line);
	my $sample = substr($f[0],0,12);
	my $percent = cos(0.6049872018 + (0.0001467884*$f[3]));
	$purity{$sample} = $percent;
}
close(EST);

opendir(CNV, $cnvDir);
my @files = readdir(CNV);
closedir(CNV);

open(SOM, "< $somatic");

my $early;
my $late;

open(EARLY, "> $pre.early");
open(LATE, "> $pre.late");
open(MID, "> $pre.mid");

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
	
	next unless exists $purity{$sample};
	
	foreach (@files) {
		$found = 1 if $_ =~ /$sample/;
	}
	next unless $found == 1;
	my $segment = `echo "$snpBed" | $bedtools intersect -a $cnvDir/$sample.bed -b stdin | head -1`;
	next if length($segment) < 1;
	#print "match\n",$segment;
	my @s = split("\t", $segment);

	my $cov = 0;

	my $af = 0;
	if ($type eq "BRCA") {
		$af = $f[63]/100;
		$cov = 40;
	} elsif ($type eq "LUAD") {
		$af = $f[80]/($f[80]+$f[81]);
		$cov = ($f[80]+$f[81]);
	} elsif ($type eq "SKCM") {
		$af = $f[83]/($f[80]+$f[83]);
		$cov = ($f[80]+$f[83]);
	} elsif ($type eq "LIHC") {
		$af = $f[43]/($f[42]);
		$cov = ($f[42]+$f[43]);
	} elsif ($type eq "LUSC" || $type eq "PRAD" || $type eq "HNSC" || $type eq "THCA") {
		$af = $f[79]/($f[80] + $f[79]);
		$cov = ($f[80]+$f[79]);
	} elsif ($type eq "LGG") {
		$af = $f[36]/($f[37] + $f[36]);
		$cov = ($f[36]+$f[37]);
	} elsif ($type eq "BLCA") {
		$af = $f[84]/($f[85] + $f[84]);
		$cov = ($f[84]+$f[85]);
	} elsif ($type eq "UCEC") {
		$f[42] =~ s/%//;
		$af = $f[42]/100;
		$cov = 40;
	} elsif ($type eq "CESC") {
		$af = $f[59]/($f[58] + $f[59]);
		$cov = ($f[58]+$f[59]);
	} else {
		die "$!: Invalid cancer type provided.\n\n";
	}
	
	$af = $af/$purity{$sample};
	my $cn = $s[3]/$purity{$sample};
	
	if ($cn > 2.2) {
		my @cutoff = &bootstrap($af,$cov,$cn);

		if ($cutoff[0] < 0.9/$cn) {
			print LATE $line;
		} elsif ($cutoff[1] > 1.8/$cn) {
			print EARLY $line;
		} elsif ($cutoff[1] >= 0.9/$cn) {
			print MID $line;
		}
	} elsif ($cn <= 2.2 && $cn >= 1.8) {
		my @cutoff = &bootstrap($af,$cov,2);

		if ($cutoff[0] < 0.45) {
			print LATE $line;
		} elsif ($cutoff[1] >= 0.45 && $cutoff[1] <= 0.55) {
			print MID $line;
		} elsif ($cutoff[1] > 0.55) {
			print EARLY $line;
		}
	}
}


close(SOM);
close(EARLY);
close(LATE);
#print $early;

#print $late;
sub bootstrap {
	my $af = shift;
	my $cov = shift;
	my $ploidy = shift;

	my $NUM_TOSSES = $cov/2;
	my $RUNS = 10000;

	my @collect = map 0, 0 .. $NUM_TOSSES;
	my $total = 0;
	my $sumDev = 0;

	for ( 1 .. $RUNS )
	{
		my $tailsCt = 0;

		# Toss the coin NUM_TOSSES times.
		#
		for ( 1 .. $NUM_TOSSES )
		{

			$tailsCt++ if rand 2 <= $af;
		}
		$collect[ $tailsCt ] ++;
	}

	foreach my $tailsCt ( 0 .. $NUM_TOSSES )
	{
		$total += $tailsCt * $collect[ $tailsCt ]
	}

	my $mean = $af*($cov/2);

	foreach my $tailsCt ( 0 .. $NUM_TOSSES )
	{
		$sumDev += (($tailsCt- $mean)**2) * $collect[ $tailsCt ];
	}

	my $stdev = sqrt($sumDev/$RUNS);

	my $l95 = ($mean - 1.960 * ($stdev / sqrt($RUNS)))/$NUM_TOSSES;
	my $u95 = ($mean + 1.960 * ($stdev / sqrt($RUNS)))/$NUM_TOSSES;
	return($l95,$u95);
}