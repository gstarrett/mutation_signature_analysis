#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $infile = shift; # Input fastq files for alignment
my $ref = shift; # Reference fasta

####### comment out below to use bwa
#goto START;
#######

### Version 0.4
### Last updated: August 4, 2016
### Developed by starr114@umn.edu
### University of Minnesota-Twin Cities

my $aln;
# Path to bwa mem command
my $bwamem = "/Applications/bioinf_tools/bwa-0.7.5a/bwa mem";

my (@summary, @prefilter, @out, %readLenHash);

# Add ref to first position of input file 
my %phredHash;
my $val = 0;
my @phred_ascii = split('', "!\"#$%&'()*+,-./0123456789:;<=>?\@ABCDEFGHIJ");
foreach (@phred_ascii) {
	$phredHash{$_} = $val;
	$val++;
}

my $format = "fasta";

if ($infile =~ /fasta$|fa$/) {
	print "Input file is in fasta format...\n";
	my $newfile = "tmp.fa";
	system("cat $ref $infile > $newfile");
	$aln = `$bwamem $ref $newfile`;
	system("rm $newfile");
} elsif ($infile =~ /fastq$|fq$/) {
	print "Input file is in fastq format...\n";
	open(REF, "< $ref");
	my @reflines = <REF>;
	chomp($reflines[1]);
	my $qual = "+\n" . "I" x length($reflines[1]);
	open(QUAL, "> ref_tmp.fq");
	print QUAL join("\n", "\@reference_sequence",$reflines[1],$qual,"");
	system("cat ref_tmp.fq $infile > tmp.fq");
	close(QUAL);
	close(REF);
	$aln = `$bwamem $ref tmp.fq`;
	system("rm ref_tmp.fq");
	system("rm tmp.fq");
	$format ="fastq";
	open(SAM, "> $infile.sam");
	print SAM $aln;
	close(SAM);
} else {die "Error: wrong file format extension.\n\n"}

# Comment out below for sam input

my @alns = grep(!/^@/,split("\n",$aln));

my @seqs;
my @quals;
my $refline = shift @alns;
my @reff = split("\t", $refline);
my $refseq = $reff[9];
my $reflen = length($refseq);
my $i=0;
my $bp=0;
my %baseHash;
for my $line (@alns) {
	my %cigarstr;
	my @f = split("\t", $line);
	my @cigar = ("M","I","D","N","S","H","P","=","X");
	foreach(@cigar) {
		if ($f[5] =~ /(\d+)$_/) {
			$cigarstr{$_} = $1; 
		} else {
			$cigarstr{$_} = 0;
		}
	}
	my $seq = "-" x ($f[3] - 1) . substr($f[9], $cigarstr{"S"}, $cigarstr{"M"}+1);
	my $qual = "-" x ($f[3] - 1) . substr($f[10], $cigarstr{"S"}, $cigarstr{"M"}+1);
	unless ($cigarstr{"M"} > (0.50*$reflen)) {next} else {
		if ($format eq "fastq") {
			push(@quals,$qual);
		}
		push(@seqs, $seq);
		$i++;
		$bp=$bp+$cigarstr{"M"};
	}
}
open(SUMMARY, "> $infile.summary.txt");
print SUMMARY "Number of reads used:\t$i\n";


############# For reading in tab file
#START:

#open(IN, "< $infile");
#@seqs = <IN>;
#$refseq = shift(@seqs);

#############

my %posHash;
my %readHash;

my @ref = split('', $refseq);
my @qscore;
for (my $idx=0; $idx <= $#seqs; $idx++) {
	chomp($seqs[$idx]);
	my @bases = split('',$seqs[$idx]);
	if ($format eq "fastq") {
		@qscore = split('', $quals[$idx]);
	}
	my $prev = "";
	my %counthash;
	my %readBaseHash;
	my $prevContext = "";
	for (my $sidx=0; $sidx < $#bases; $sidx++) {
		next unless $ref[$sidx] =~ /[ATCG]/ && $bases[$sidx] ne "-";
		
		if (exists $readBaseHash{$ref[$sidx]}) {
			$readBaseHash{$ref[$sidx]}++;
		} else {
			$readBaseHash{$ref[$sidx]} = 1;
		}
		
		if (exists $baseHash{$ref[$sidx]}) {
			$baseHash{$ref[$sidx]}++;
		} else {
			$baseHash{$ref[$sidx]} = 1;
		}
		if ($bases[$sidx] ne $ref[$sidx]) {
			next unless $bases[$sidx] =~ /[ATCG]/;
			if ($format eq "fastq") {
				next unless $phredHash{$qscore[$sidx]} > 30;
			}
			
			my $pos = $sidx + 1;
			
			## Remove DNPs (likely misalignments)
			
			if (exists $posHash{$sidx}) {
				$posHash{$sidx}--;
				delete $counthash{$prevContext} if exists $counthash{$prevContext};
				next;
			}
			
			if (exists $posHash{$pos}) {
				$posHash{$pos}++;
			} else {
				$posHash{$pos} = 1;
			}
			
			if (exists $readHash{$idx}) {
				$readHash{$idx}++;
			} else {
				$readHash{$idx} = 1;
			}
			
			my $context = "$pos\t$ref[$sidx-1]$ref[$sidx]\t$ref[$sidx]$ref[$sidx+1]\t$ref[$sidx-1]$ref[$sidx]$ref[$sidx+1]\t$ref[$sidx]>$bases[$sidx]";
			if (exists $counthash{$context}) {
				$counthash{$context}++;
			} else {
				$counthash{$context} = 1;
			}
			$prevContext = $context;
		}
	$readLenHash{$idx} = \%readBaseHash;
	}
	for my $key (keys %counthash) {
		push(@prefilter, "$idx\t$key\t$counthash{$key}\t$#bases");
	}
}


my $posOutliers = calculate(\%posHash);
my $readOutliers = calculate(\%readHash);

my @posFlt;

print SUMMARY "===Unfiltered base counts===\n";
for my $baseKey (keys %baseHash) {
	print SUMMARY "$baseKey:\t$baseHash{$baseKey}\n";
}
print SUMMARY "===Filtered base counts===\n";
foreach (@prefilter) {
	if (scalar @$posOutliers == 0) {
		push(@posFlt, "$_\tPASS");
	} else {
		my $bool = 0;
		for my $posOutlier (@$posOutliers) {
			if ($bool == 0) {
				if ($_ =~ /\t$posOutlier\t/) {
					$bool = 1;
				} else {
					next;
				}
			} else {
				next;
			}
		}
		if ($bool == 0 ) {
			push(@posFlt, "$_\tPASS");
		}
		else {
			push(@posFlt, "$_\tOUTLIER");
		}
	}
}

foreach (@posFlt) {
	if (scalar @$readOutliers == 0) {
		push(@out, "$_\tPASS");
	} else {
		my $bool = 0;
		for my $readOutlier (@$readOutliers) {
			if ($bool == 0) {
				if ($_ =~ /^$readOutlier\t/) {
					$bool = 1;
				} else {
					next;
				}
			} else {
				next;
			}
		}
		if ($bool == 0) {
			push(@out, "$_\tPASS");
		} else {
			push(@out, "$_\tOUTLIER");
		}		
	}
}

if (scalar @$readOutliers == 0) {
	for my $baseKey (keys %baseHash) {
		print SUMMARY "$baseKey:\t$baseHash{$baseKey}\n";
	}
} else {
	for my $readOutlier (@$readOutliers) {
		my $int = $readLenHash{$readOutlier};
		for my $baseKey (keys %{$int}) {
			$baseHash{$baseKey} = $baseHash{$baseKey} - ${$int}{$baseKey};
		}
	}
	for my $baseKey (keys %baseHash) {
		print SUMMARY "$baseKey:\t$baseHash{$baseKey}\n";
	}
}

open(OUT, "> $infile.counts.txt");

foreach (@out) {
	print OUT $_, "\n";
}


close(OUT);
close(SUMMARY);

# Calculate if any read has unusually high (u+3sd) frequency of mutations and remove
# Calculate if any position has unusually high (u+3sd) frequency of mutations

sub calculate {
	my $hash = shift;
	#print Dumper(%{$hash});
	my $sum = 0;
	my $sumDev = 0;
	my @outliers;
	for my $key (keys %$hash) {
		$sum += $$hash{$key};
	}
	my $mean = $sum/(scalar keys %$hash);
	for my $key (keys %$hash) {
		$sumDev += ($$hash{$key} - $mean) ** 2;
	}
	my $stdev = sqrt($sum/(scalar keys %$hash));
	for my $key (keys %$hash) {
		if ($$hash{$key} > ($mean + 3*$stdev)) {
			push(@outliers, $key);
		}
	}
	return(\@outliers);
}
