#!/usr/bin/perl -w
use strict;
use Bio::PrimarySeq;
use base qw(Bio::SeqIO);
use Data::Dumper;
use threads;
use Term::ReadKey;

my %count;
my $in = shift;
my $n = 10;
my $seqio = Bio::SeqIO->new(-file => $in, -format => "fasta");
while( my $seq = $seqio->next_seq ) {
	my $m1;
	my $m2;
	print "Analyzing " . $seq->id . "... ";
	#my $s = $seq->seq();

	# while (my $c = chop $s) {
	ReadMode 3;
	while (my $c = <$seq->seq()> ) {
		#print "$c\n";
		my $c = uc($c);
		unless (defined $m2) {
			$m2 = $c;
			next;
		}
		unless (defined $m1) {
			$m1 = $c;
			next;
		}
		my $trint = $c . $m1 . $m2;
		my $dint = $c . $m1;
		#print "$trint\n";
		next if $trint =~ /N/;
		if (exists $count{$trint}) {
			$count{$trint}++;
			$count{$dint}++;
		} else {
			$count{$trint} = 1;
			$count{$dint} = 1;
		}
		$m2 = $m1;
		$m1 = $c;

	}
	print "done\n";
	ReadMode 0;
}

print Dumper(%count);
