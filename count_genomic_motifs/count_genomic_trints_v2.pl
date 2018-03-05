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
	my $m1;
	my $m2;
	my $id = $seq->id;
	next if $id =~ /_|M/;
	print "Analyzing " . $id . "... ";
	my $s = $seq->seq();
	
	while (my $c = chop $s) {
	# while (my $c = $seq->seq() ) {
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
		#print "$trint\n";
		$m2 = $m1;
		$m1 = $c;
		next if $trint =~ /N/;
		if (exists $count{$trint}) {
			$count{$trint}++;
		} else {
			$count{$trint} = 1;
		}

		
	}
	print "done\n"
}

for my $key (keys %count) {
	print $key, "\t", $count{$key},"\n";
}