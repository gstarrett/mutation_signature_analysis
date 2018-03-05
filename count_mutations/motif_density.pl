#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $in;
my $win=600;
my $int=60;
my $kmer=2;

GetOptions ("window=i" => \$win,
			"interval=i" => \$int,
			"file=s" => \$in,
			"kmer=i" => \$kmer
) or die("Error in command line arguments\n");

print STDERR "Using these settings:\n\tkmer: $kmer\n\tinterval: $int\n\twindow: $win\n";

my $dint = calcMotifs($kmer);
my $trint = calcMotifs($kmer+1);

my %posHash;
#my $lnum = 0;

open(IN, "< $in");
while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t",$line);
	my $seq = $f[4];
	my $pos = 0;
	until ($pos+$win > length($seq)) {
		my %kmercount;
		my $seqw = substr($seq,$pos,$win);
		for (my $i = 0; $i < length($seqw)-$kmer; $i++) {
			my $dn = substr($seqw, $i, $kmer);
			if (exists $kmercount{$dn}) {
				$kmercount{$dn}++;
			} else {
				$kmercount{$dn} = 1;
			}
		}
		for (my $i = 0; $i < length($seqw)-$kmer-1; $i++) {
			my $tn = substr($seqw, $i, $kmer+1);
			if (exists $kmercount{$tn}) {
				$kmercount{$tn}++;
			} else {
				$kmercount{$tn} = 1;
			}
		}
		for my $key (keys %kmercount) {
			unless ($key =~ /-/) {
				print join("\t", @f[0..3,13 ], $pos, $key, $kmercount{$key}),"\n";
			}
		}
		$pos += $int;
	}
	#$lnum++;
}

close(IN);

sub calcMotifs{
	my @out;
	my $kmer = shift;
	my @bases = ("A","T","C","G");
	for my $base1 (@bases) {
		for my $base2 (@bases) {
			push(@out,"$base1$base2");
		}
	}
}