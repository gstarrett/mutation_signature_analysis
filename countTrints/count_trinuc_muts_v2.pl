#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Fasta;
use Statistics::R;

my $n = 0;
my $fasta = shift;
my $variants = shift;
my $name = shift;
my $Ctrinucs;
my $Gtrinucs;

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
	my @fields = split("\t",$line);
	my $chr = shift @fields;
	my $pos = shift @fields;
	my $id = shift @fields; # for full vcf not pvcf
	my $ref = shift @fields;
	my $alt = shift @fields;
	#print "$chr $pos $ref $alt\n";
	if ($ref eq "C") {
		my $flank = uc($db->seq($chr, $pos-10+$n => $pos+10+$n));
		if (exists $repeatHash{$flank}) {
			next;
		} else {
			$repeatHash{$flank} = 1;
		}
		my $seq = uc($db->seq($chr, $pos-1+$n => $pos+1+$n));
		my @bases = split('', $seq);
		$seq = $bases[0] . $bases[1] . $bases[2];
		$Ctrinucs .= "$seq\t$ref>$alt\t1\n";
		
# 		if ($seq =~ /[AaTtCcGg][Cc][Gg]/) {
# 			next;
# 		} else {
# 			$Ctrinucs .= "$seq\n";
# 		}
	} elsif ($ref eq "G" ) {
		my $flank = uc($db->seq($chr, $pos-10+$n => $pos+10+$n));
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
		my $seq = uc($db->seq($chr, $pos-1+$n => $pos+1+$n));
		my @bases = split('', $seq);
		$seq = $compDict{$bases[2]} . $compDict{$bases[1]} . $compDict{$bases[0]};
		$Gtrinucs .= "$seq\t$compDict{$ref}>$compDict{$alt}\t1\n";
		
# 		if ($seq =~ /[Cc][Gg][AaTtCcGg]/) {
# 			next;
# 		} else {
# 			$Gtrinucs .= "$seq\n";
# 		}
	} else { next }
}
open (OUTC, "> $name.C_trinucs_count.txt");
open (OUTG, "> $name.G_trinucs_count.txt");

print OUTC $Ctrinucs;
print OUTG $Gtrinucs;

close(OUTC);
close(OUTG);

### Process maf file
# data.w <- mapply("*", data.c, data.radj, SIMPLIFY=F)

my $R = Statistics::R->new();
$R -> start();
print $R->is_started(),"\n";
#$R -> set('filename', "$name.C_trinucs_count.txt");
#$R -> set('norm', @norm);
#$R -> set('fileout', "$name.C_trinucs.pdf");
#print $R->get('filename'),"\n";
$R->run(qq`data <- read.table("$name.C_trinucs_count.txt")`);
$R->run(qq`fileout <- "$name.C_trinucs.pdf"`);
my $cmds = <<'EOF';
library("reshape")

data.m <- melt(data, id=c(1:2), measure=c(3))
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