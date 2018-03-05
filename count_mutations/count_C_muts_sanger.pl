#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $infile = shift;
my $ref = shift;

####### comment out below to use bwa
#goto START;
#######

my $aln;
#my $clustalo = "/Users/gabestarrett/Desktop/clustal-omega-1.2.0-macosx";
#my $trim = "T";
my $bwamem = "/Applications/bioinf_tools/bwa-0.7.5a/bwa mem";

# if ($trim eq "T") {
# 	open(FILE, "< $infile");
# 	my $trimmedFile = <FILE>;
# 	$trimmedFile .= <FILE>;
# 	while (my $line = <FILE>) {
# 		if ($line =~ /^>/) {
# 			$trimmedFile .= $line;
# 		} else {
# 			my $subseq = substr $line, 300, 500;
# 			$trimmedFile .= $subseq . "\n";
# 		}
# 	}
# 	close(FILE);
# 	#print $trimmedFile;
# 	$aln = `echo "$trimmedFile" | $clustalo -i - --output-order=input-order --wrap=2000`;
# 	
# } else {
# 
# 	$aln = `$clustalo -i $infile --output-order=input-order --wrap=2000`;
# 
# }

# add ref to input file 
my %phredHash;
my $val = 0;
my @phred_ascii = split('', "!\"#$%&'()*+,-./0123456789:;<=>?\@ABCDEFGHIJ");
foreach (@phred_ascii) {
	$phredHash{$_} = $val;
	$val++;
}

my $format = "fasta";

if ($infile =~ /fasta/) {
	
	my $newfile = "tmp.fa";
	system("cat $ref $infile > $newfile");
	$aln = `$bwamem $ref $newfile`;
	system("rm $newfile");
} elsif ($infile =~ /fastq/) {
	print "Input file is in fastq format...\n";
	open(REF, "< $ref");
	my @reflines = <REF>;
	chomp($reflines[1]);
	my $qual = "+\n" . "I" x length($reflines[1]);
	open(QUAL, "> ref_tmp.fq");
	print QUAL join('\n', "\@reference_sequence",$reflines[1],$qual);
	system("cat ref_tmp.fq $infile > tmp.fq");
	close(QUAL);
	close(REF);
	$aln = `$bwamem $ref ref_tmp.fq`;
	system("rm ref_tmp.fq");
	$format ="fastq";
}

# Comment out below for sam input

my @alns = grep(!/^@/,split("\n",$aln));


my @seqs;
my @quals;
my $refline = shift @alns;
print $refline,"\n";	
my @reff = split("\t", $refline);
my $refseq = $reff[9];
my $reflen = length($refseq);
#print $refseq,"\n";
my $i=0;
my $bp=0;
my %baseHash;
for my $line (@alns) {
	my %cigarstr;
	#print $line,"\n";
	my @f = split("\t", $line);
	#print "$f[5]\n";
	#push(@seqs,$f[9]);
	if ($format eq "fastq") {
		push(@quals,$f[10]);
	}
	my @cigar = ("M","I","D","N","S","H","P","=","X");
	foreach(@cigar) {
		if ($f[5] =~ /(\d+)$_/) {
			$cigarstr{$_} = $1; 
		} else {
			$cigarstr{$_} = 0;
		}
		#print "$_\t" . "$1\n";
	}
	my $seq = "-" x ($f[3] - 1) . substr($f[9], $cigarstr{"S"}, $cigarstr{"M"}+1);
	unless ($cigarstr{"M"} > (0.50*$reflen)) {next} else {
		push(@seqs, $seq);
		$i++;
		$bp=$bp+$cigarstr{"M"};
	}
}
print STDERR "==============================\n$infile\n";
print STDERR "Number of reads used:\t$i\n";
print STDERR "Basepairs analyzed:\t$bp\n";


#print $aln,"\n";

#my @seqs = split(/>.+\n/, $aln);

#shift(@seqs);
#print Dumper(@seqs);
#print $#seqs,"\n";
#chomp($seqs[0]);

############# For reading in tab file
#START:

#open(IN, "< $infile");
#@seqs = <IN>;
#$refseq = shift(@seqs);

#############

my @ref = split('', $refseq);
my @qscore;

for (my $idx=0; $idx <= $#seqs; $idx++) {
	chomp($seqs[$idx]);
	#print "$seqs[$idx]";
	my @bases = split('',$seqs[$idx]);
	if ($format eq "fastq") {
		@qscore = split('', $quals[$idx]);
	}
	my $prev = "";
	my %counthash;
	for (my $sidx=0; $sidx < $#bases; $sidx++) {
		next unless $ref[$sidx] =~ /[ATCG]/ && $bases[$sidx] ne "-";
		if (exists $baseHash{$ref[$sidx]}) {
			$baseHash{$ref[$sidx]}++;
		} else {
			$baseHash{$ref[$sidx]} = 1;
		}
		if ($bases[$sidx] ne $ref[$sidx]) {
			next unless $bases[$sidx] =~ /[ATCG]/;
			if ($format eq "fastq") {
				next unless $phredHash{$qscore[$sidx]} > 20;
			}
			my $pos = $sidx + 1;
			my $context = "$pos\t$ref[$sidx-1]$ref[$sidx]\t$ref[$sidx]$ref[$sidx+1]\t$ref[$sidx-1]$ref[$sidx]$ref[$sidx+1]\t$ref[$sidx-1]$ref[$sidx]$ref[$sidx+1]$ref[$sidx+2]\t$ref[$sidx]>$bases[$sidx]";
			if (exists $counthash{$context}) {
				$counthash{$context}++;
			} else {
				$counthash{$context} = 1;
			}
		}
	}
	for my $key (keys %counthash) {
		print "$idx\t$key\t$counthash{$key}\t$#bases\n";
	}
}
for my $baseKey (keys %baseHash) {
	print STDERR "$baseKey:\t$baseHash{$baseKey}\n";
}
print STDERR "==============================\n";
