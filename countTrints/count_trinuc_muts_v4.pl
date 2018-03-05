#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Fasta;
use Statistics::R;

my $n = -1;
my $type = shift;
my $fasta = shift;
my $variants = shift;
my $stdout = 0;

unless (defined $type && defined $fasta && defined $variants) {
	my $usage="Error: Not enough inputs\n\n\tUsage: perl count_trinuc_muts <maf, vcf, or pvcf> <reference> <variants>\\nn";
	die $usage;
}

my $name = "$variants";
my $trinucs = "#chr\tpos\t5'tetranuc\t3'tetranuc\ttrinuc\tmut\ttrinuc_mut\tstrand\tTCWcount\tCGcount\n";


my %trinucHash;

my %compDict = (
	A => "T",
	T => "A",
	C => "G",
	G => "C"
);

my %normalization = (
	ACA =>	0.07,
	ACC =>	0.06,
	ACG =>	0.02,
	ACT =>	0.06,
	CCA =>	0.10,
	CCC =>	0.09,
	CCG =>	0.02,
	CCT =>	0.09,
	GCA =>	0.07,
	GCC =>	0.08,
	GCG =>	0.02,
	GCT =>	0.07,
	TCA =>	0.07,
	TCC =>	0.08,
	TCG =>	0.01,
	TCT =>	0.08
);

my @norm = (
	0.07,
	0.06,
	0.02,
	0.06,
	0.10,
	0.09,
	0.02,
	0.09,
	0.07,
	0.08,
	0.02,
	0.07,
	0.07,
	0.08,
	0.01,
	0.08
);

#system("awk 'BEGIN {FS="\t"}; {print "chr"$5"\t"$6"\t"$18"\t"$12"\t"$16}' $filein");

print "creating fasta db... ";
my $db = Bio::DB::Fasta->new($fasta);
print "complete!\n";

my %repeatHash;

open (VAR, "< $variants");
while (my $line = <VAR>) {
	chomp($line);
	next if $line =~ /^$/;
	my ($chr,$pos,$ref,$alt);
	my $sample = "";
	my @fields = split("\t",$line);
	unless ($type eq "maf") {
		$chr = shift @fields;
		$pos = shift @fields;
		$pos++;
		my $id = shift @fields if $type eq "vcf"; # for full vcf not pvcf
		$ref = shift @fields;
		$alt = shift @fields;
		$sample = shift @fields if $type eq "pvcf"; # if sample name is available
	} else {
		next unless $fields[9] eq "SNP";
		$chr = $fields[4];
		$pos = $fields[5];
		$pos++;
		$ref = $fields[10];
		$alt = $fields[12];
		$sample = substr($fields[15],0,12);
	}
	next unless length($alt) == 1;
	#print "$chr $pos $ref $alt\n";
	if ($ref =~ /^[TC]$/) {
		my $CGcxt=0;
		my $TCWcxt=0;
		my $WGAcxt=0;
		my $flank = uc($db->seq($chr, $pos-10+$n => $pos+10+$n));
		my $flank40 = uc($db->seq($chr, $pos-20+$n => $pos+20+$n));
		while ($flank40 =~ /TC[AT]/g) { $TCWcxt++ }
		while ($flank40 =~ /[AT]GA/g) { $WGAcxt++ }
		my $A3cxt = $WGAcxt + $TCWcxt;
		while ($flank40 =~ /[CG]/g) { $CGcxt++ }
		if (exists $repeatHash{$flank}) {
			next;
		} else {
			$repeatHash{$flank} = 1;
		}
		my $seq = uc($db->seq($chr, $pos-2+$n => $pos+2+$n));
		my @bases = split('', $seq);
		$seq = $bases[1] . $bases[2] . $bases[3];
		my $tetraseq5 = $bases[0] . $bases[1] . $bases[2] . $bases[3];
		my $tetraseq3 = $bases[1] . $bases[2] . $bases[3] . $bases[4];
		my $fullstring = "$bases[1]\[$ref>$alt\]$bases[3]";
		$trinucs .= "$chr\t$pos\t$tetraseq5\t$tetraseq3\t$seq\t$ref>$alt\t$fullstring\t1\t$A3cxt\t$CGcxt\t$sample\n";
		
	} elsif ($ref =~ /^[AG]$/ ) {
		my $CGcxt=0;
		my $TCWcxt=0;
		my $WGAcxt=0;
		my $flank = uc($db->seq($chr, $pos-10+$n => $pos+10+$n));
		my $flank40 = uc($db->seq($chr, $pos-20+$n => $pos+20+$n));
		while ($flank40 =~ /TC[AT]/g) { $TCWcxt++ }
		while ($flank40 =~ /[AT]GA/g) { $WGAcxt++ }
		my $A3cxt = $WGAcxt + $TCWcxt;
		while ($flank40 =~ /[CG]/g) { $CGcxt++ }
		my @fbases = split('', $flank);
		my @revcomp;
		foreach(@fbases) {
			unshift(@revcomp, $compDict{$_});
		}
		my $rcflank = join('', @revcomp);
		if (exists $repeatHash{$rcflank}) {
			next;
		} else {
			$repeatHash{$rcflank} = 1;
		}
		my $seq = uc($db->seq($chr, $pos-2+$n => $pos+3+$n));
		my @bases = split('', $seq);
		$seq = $compDict{$bases[3]} . $compDict{$bases[2]} . $compDict{$bases[1]};
		my $tetraseq3 = $compDict{$bases[3]} . $compDict{$bases[2]} . $compDict{$bases[1]} . $compDict{$bases[0]};
		my $tetraseq5 = $compDict{$bases[4]} . $compDict{$bases[3]} . $compDict{$bases[2]} . $compDict{$bases[1]};
		my $fullstring = "$compDict{$bases[3]}\[$compDict{$ref}>$compDict{$alt}\]$compDict{$bases[1]}";
		$trinucs .= "$chr\t$pos\t$tetraseq5\t$tetraseq3\t$seq\t$compDict{$ref}>$compDict{$alt}\t$fullstring\t-1\t$A3cxt\t$CGcxt\t$sample\n";
	
	} else { next }
}

if ($stdout == 0) {
	open (OUT, "> $name.count.txt");
	print OUT $trinucs;
	close(OUT);
} else {
	print $trinucs;
}
# count all TCW:WGA and sum flank C:Gs

# count all C:Gs and sum flank TCW:WGA

### Process maf file
# data.w <- mapply("*", data.c, data.radj, SIMPLIFY=F)

exit;

my $R = Statistics::R->new();
$R -> start();
#print $R->is_started(),"\n";
#$R -> set('filename', "$name.C_trinucs_count.txt");
#$R -> set('norm', @norm);
#$R -> set('fileout', "$name.C_trinucs.pdf");
#print $R->get('filename'),"\n";
$R->run( qq`setwd("/Users/gabestarrett")` );
$R->run(qq`data <- read.table("$name.trinucs_count.txt", header=T, sep="\t")`);
$R->run(qq`fileout <- "$name.C_trinucs.pdf"`);
my $cmds = <<'EOF';
library("reshape")

data.m <- melt(data, id=c(5:6), measure=c(abs(8)))
data.c <- cast(data.m, V1~V2, sum)
rownames(data.c) <- data.c\$V1
data.c\$V1 <- NULL
data.rfreq <- rowsum(data.c)/sum(data.c)
data.radj <- mapply("/",data.rfreq,norm,SIMPLIFY=F)

par(pin=c(12,5))
pdf(fileout)
barplot(t(as.matrix(data.c/sum(data.c))), col=c("blue","black","red"), ylim=c(0,0.5), border=NA)
dev.off()
print("success!")
EOF

$R->stop();
my $run = $R->run($cmds);
print $run, "\n";



