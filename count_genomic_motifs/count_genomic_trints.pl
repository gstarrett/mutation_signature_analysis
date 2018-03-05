#!/usr/bin/perl -w
use strict;
use Bio::PrimarySeq;
use base qw(Bio::SeqIO);
use Data::Dumper;
use threads;

my %count;
my $in = shift;
my $n = 10;
my $seqio = Bio::SeqIO->new(-file => $in, -format => "fasta");
while( my $seq = $seqio->next_seq ) {
	my @thr;
	print "Analyzing " . $seq->id . "\n";
	my @nts = split '', $seq->seq();
	my $w = $#nts/$n;
	for (my $start = 0; $start <= $#nts-$w; $start+=$w) {
		my $end = $start + $w;
		my @args = (\@nts, $start, $end);
		push (@thr, threads->create(\&countThread,@args));
	}
	foreach (@thr) {
		my @join = $_->join();
	}
}

print Dumper(%count);

sub countThread {
	my $seqnt = shift;
	my $start = shift;
	my $end = shift;
	print "$start-$end started\n";
	for (my $i = $start; $i <= $end-1; $i++) {
		my $trint = ${$seqnt}[$i-1] . ${$seqnt}[$i] . ${$seqnt}[$i+1];
		next if $trint =~ /N/;
		if (exists $count{$trint}) {
			$count{$trint}++;
		} else {
			$count{$trint} = 1;
		}
	}
	print "$start-$end completed\n";
}