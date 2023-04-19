#!/usr/bin/perl -w
use strict;
no warnings 'uninitialized';
my $bedtools = "/Applications/bioinf_tools/bedtools-2.17.0/bin/bedtools";
my $vcftools = "/Users/gabestarrett/vcftools_0.1.12b/bin/vcftools ";

my $w10k = "/Volumes/HD_2/Refs/gatk_hg19/hg19.10k.bed";
my $w1k = "/Volumes/HD_2/Refs/gatk_hg19/hg19.1k.bed";

my $in = shift;
my $out = "$in.out";

system("awk  \'BEGIN {OFS=\"\t\"; FS=\"\t\"}; {print \$1,\$2,\$2,\$4\">\"\$5}\' $in > temp.bed");

system("$bedtools intersect -c -a $w1k -b temp.bed > $out.intersect.1k.bed");

system("$bedtools intersect -c -a $w10k -b $in > $out.intersect.10k.bed");

system("awk \'\$4~/^C>/ {print}\' temp.bed > temp.c.bed");
system("$bedtools intersect -c -a $w1k -b temp.c.bed > $out.intersect.1k.c.bed");
system("$bedtools intersect -c -a $w10k -b temp.c.bed > $out.intersect.10k.c.bed");
system("rm temp.c.bed");

system("awk \'\$4~/^G>/ {print}\' temp.bed > temp.g.bed");
system("$bedtools intersect -c -a $w1k -b temp.g.bed > $out.intersect.1k.g.bed");
system("$bedtools intersect -c -a $w10k -b temp.g.bed > $out.intersect.10k.g.bed");
system("rm temp.g.bed");

system("rm temp.bed");

system("awk \'\$4 > 2 {print}\' $out.intersect.10k.c.bed > $out.intersect.10k.cg.tmp.bed");
system("awk \'\$4 > 2 {print}\' $out.intersect.10k.g.bed >> $out.intersect.10k.cg.tmp.bed");

system("$bedtools merge -i $out.intersect.10k.cg.tmp.bed > $out.intersect.10k.cg.bed");
system("rm $out.intersect.10k.cg.tmp.bed");

system("$vcftools --vcf $in --recode --bed $out.intersect.10k.cg.bed --out $in.cg");

exit;

my %dataHash;
# maybe put in a system sort before this step;
#open (IN, "< $in");
#open (OUT, "> $out.imd.txt");

open (OUT, "> $in.imd.txt");
while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);
	my $chr = shift(@f);
	if (exists $dataHash{$chr}) {
		push(@{$dataHash{$chr}}, \@f);
	} else {
		$dataHash{$chr} = [\@f];
	}
}

print OUT "#CHR\tPOS\tMUT\tPROX_IMD\tDIST_IMD\tC_PROX_IMD\tC_DIST_IMD\tTYPE\n";

for my $key (keys %dataHash) {
	my $max = scalar @{$dataHash{$key}};
	for my $idx (1 .. $max-2) {
		my $prev = ${$dataHash{$key}}[$idx-1];
		my $dra = ${$dataHash{$key}}[$idx];
		my $next = ${$dataHash{$key}}[$idx+1];
		my $prox = ${$dra}[0] - ${$prev}[0];
		my $dist = ${$next}[0] - ${$dra}[0];
		
		my $cprev = "";
		my $cprox = "";
		my $cnext = "";
		my $cdist = "";
		
		if (${$dra}[1] eq "C") {
			my $n = 1;
			until (${$dataHash{$key}}[$idx-$n] eq "C") {
				$n++;
			}
			if (${$dataHash{$key}}[$idx-$n] eq "C") {
				$cprev = ${$dataHash{$key}}[$idx-$n];
				$cprox = ${$dra}[0] - ${$cprev}[0];
			}
			
			$n = 1;
			until (${$dataHash{$key}}[$idx+$n] eq "C") {
				$n++;
			}
			if (${$dataHash{$key}}[$idx-$n] eq "C") {
				$cnext = ${$dataHash{$key}}[$idx+$n];
				$cdist = ${$cnext}[0] - ${$dra}[0];
			}
		} elsif (${$dra}[1] eq "G") {
			my $n = 1;
			until (${$dataHash{$key}}[$idx-$n] eq "G") {
				$n++;
			}
			if (${$dataHash{$key}}[$idx-$n] eq "G") {
				$cprev = ${$dataHash{$key}}[$idx-$n];
				$cprox = ${$dra}[0] - ${$cprev}[0];
			}
			
			$n = 1;
			until (${$dataHash{$key}}[$idx+$n] eq "G") {
				$n++;
			}
			if (${$dataHash{$key}}[$idx-$n] eq "G") {
				$cnext = ${$dataHash{$key}}[$idx+$n];
				$cdist = ${$cnext}[0] - ${$dra}[0];
			}
		}
		my $type;
		my $mut = ${$dra}[1] . ">" . ${$dra}[2];
		if ($mut eq "C>T" || $mut eq "T>C" || $mut eq "G>A" || $mut eq "A>G") {
			$type = "transition";
		} else {
			$type = "transversion";
		}
		print OUT "$key\t${$dra}[0]\t$mut\t$prox\t$dist\t$cprox\t$cdist\t$type\n";
	}
}

close(IN);
close(OUT);