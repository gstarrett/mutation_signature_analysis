#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;

my $kmer = 2;
# open
my $in = shift;
#my $inseq = Bio::SeqIO->new(-file => $in);
open(IN, "< $in");
my @seqs = <IN>;
close(IN);
#while (my $entry = $inseq->next_seq) {
my $n = 1;
for my $seq (@seqs) {
	chomp($seq);
	#my $sample = $entry->id;
	#my $seq = $entry->seq();
	$seq =~ s/-//g;
	$seq = uc $seq;
	# p obs ratio of the total count of the kmer relative to the total count of all kmers with the same length
	my %kmercount;
	my %kmer1count;
	my %basecount;
	#print length($seq), "\n";
	for (my $i = 0; $i < length($seq)-$kmer; $i++) {
		my $dn = substr($seq, $i, $kmer+1);
		#print $dn, "\n";
		if (exists $kmer1count{$dn}) {
			$kmer1count{$dn}++;
		} else {
			$kmer1count{$dn} = 1;
		}
	}
	for (my $i = 0; $i <= length($seq)-$kmer; $i++) {
		my $dn = substr($seq, $i, $kmer);
		#print $dn, "\n";
		if (exists $kmercount{$dn}) {
			$kmercount{$dn}++;
		} else {
			$kmercount{$dn} = 1;
		}
	}
	for (my $i = 0; $i < length($seq); $i++) {
		my $base = substr($seq, $i, 1);
		#print $dn, "\n";
		if (exists $basecount{$base}) {
			$basecount{$base}++;
		} else {
			$basecount{$base} = 1;
		}
	}
	for my $key (keys %kmercount) {
		my $pobs = $kmercount{$key}/(length($seq)+1-$kmer);
		my $pexp = ($basecount{substr($key,0,1)}/length($seq))*($basecount{substr($key,-1)}/length($seq));
		print $n,"\t",$key,"\t",$pobs/$pexp,"\n";
	}
	for my $key (keys %kmer1count) {
		my $pobs = $kmer1count{$key}/(length($seq)-$kmer);
		#print substr($key,0,$kmer), substr($key,-$kmer), "\n";
		my $pexp = (($kmercount{substr($key,0,$kmer)}/(length($seq)+1-$kmer))*($kmercount{substr($key,-$kmer)}/(length($seq)+1-$kmer)))/($basecount{substr($key,1,1)}/length($seq));
		print $n,"\t",$key,"\t",$pobs/$pexp,"\n";
	}
	$n++;
}