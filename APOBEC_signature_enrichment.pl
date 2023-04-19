#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;

my %results;
my $fmt = shift;
my $kmer = 2;
my $in = shift;

my $n = 1;

if ($fmt eq "genbank" | $fmt eq "fasta") {
	my $inseq = Bio::SeqIO->new(-file => $in, -format => $fmt);
	while (my $entry = $inseq->next_seq) {
		my $id = $entry->id;
		my $seq = $entry->seq();
		my $display = "";
		if ($fmt eq "genbank") {
			$display = $entry->desc();
		}
		my $name = join("\t", $id,$display);
		$results{"full"} .= process($name,$seq,$n);
		$n++;
	}
} elsif ($fmt eq "tab") {
	open(IN, "< $in");
	my @seqs = <IN>;
	close(IN);
	for my $seq (@seqs) {
		chomp($seq);
		$results{"full"} .= process($n,$seq,$n);
		$n++;
	}	
} elsif ($fmt eq "BKdetail") {
	open(IN, "< $in");
	my $head = <IN>;
	while (my $line = <IN>) {
		my @f = split("\t",$line);
		my $id = $f[2];
		my $group = $f[3];
		my $clin = $f[13];	
		my $name = join("\t", $id, $group, $clin);
		
		$results{"full"} .= process($name,$f[4],$n);
		$results{"agno"} .= process($name,$f[5],$n);
		$results{"VP1"} .= process($name,$f[6],$n);
		$results{"VP2"} .= process($name,$f[7],$n);
		$results{"VP3"} .= process($name,$f[8],$n);
		$results{"LT"} .= process($name,$f[9],$n);
		$results{"st"} .= process($name,$f[10],$n);
		$n++;
	}
	close(IN);
} else {
	print "Invalid format!\n";
}

for my $key (keys %results) {
	open(OUT, "> $in.$key.txt");
	print OUT $results{$key};
}

sub process {
	my $sample = shift;
	my $seq = shift;
	my $n = shift;
	my $output;
	
	#next if length($seq) < 1000;
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
		$output .= $n . "\t" . $key . "\t" . $pobs/$pexp . "\t" . "$sample" . "\n";
	}
	for my $key (keys %kmer1count) {
		my $pobs = $kmer1count{$key}/(length($seq)-$kmer);
		#print substr($key,0,$kmer), substr($key,-$kmer), "\n";
		my $pexp = (($kmercount{substr($key,0,$kmer)}/(length($seq)+1-$kmer))*($kmercount{substr($key,-$kmer)}/(length($seq)+1-$kmer)))/($basecount{substr($key,1,1)}/length($seq));
		$output .= $n . "\t" . $key . "\t" . $pobs/$pexp . "\t" . "$sample" . "\n";
	}
	return $output;
}