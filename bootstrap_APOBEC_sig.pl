#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $n = shift;
my $list = shift;
my $data = shift;
my $script = "perl /Users/gabestarrett/Google_Drive/Scripts/count_trinuc_muts_v4_STDINOUT.pl maf /Volumes/HD_2/Refs/gatk_b37/human_g1k_v37.fasta";


open(LIST, "< $list");
my @names = <LIST>;
chomp(@names);
close(LIST);

my $bs = 10000;

open(DATA, "< $data");
my @entries = <DATA>;
close(DATA);

for (my $i=0; $i<=$bs; $i++) {
	my %count;
	my $masterCount = 0;
	my @deck = @names;
	my @results;
	my $picks = $n;
	while ($picks > 0) {  # when we have all our picks, stop
		# random number from 0..$num_left-1
		my $rand = int(rand($#deck)-1);
		# pick successful
		if (defined $deck[$rand]) {
			push (@results, $deck[$rand]);
			$picks--;
			delete $deck[$rand];
		}
	}
	#print Dumper(@results);
	my $regex = join('|',@results);
	#print $regex;
	my @matches = grep (/$regex/, @entries);
	my $out = join("",@matches);
	open(TMP, "> tmp.maf");
	print TMP $out;
	close(TMP);
	my $counts = `$script tmp.maf`;
	unlink("tmp.maf");
	my @lines = split("\n", $counts);
	for my $line (@lines) {
		my @f = split("\t", $line);
		next if $#f < 8;
		if (exists $count{$f[6]}) {
			$count{$f[6]}++;
		} else {
			$count{$f[6]} = 1;
		}
		$masterCount++;
	}
	for my $key (sort keys %count) {
		print "$i\t$key\t" . $count{$key}/$masterCount . "\n";
	}
}