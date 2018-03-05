#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $bedtools = "/Applications/bioinf_tools/bedtools2.23/bin/bedtools";

my $type = shift;
my $somatic = shift;
my $cnvDir = shift;
my $estimate = shift;
my $pre = shift;

unless (defined $pre) {
	$pre = $somatic;
}

my %purity;

open(EST, "< $estimate");
my $head = <EST>;
while (my $line = <EST>) {
	chomp($line);
	my @f = split("\t", $line);
	my $sample = substr($f[0],0,12);
	my $percent = cos(0.6049872018 + (0.0001467884*$f[3]));
	$purity{$sample} = $percent;
}
close(EST);

opendir(CNV, $cnvDir);
my @files = readdir(CNV);
closedir(CNV);

open(SOM, "< $somatic");

my $early;
my $late;

open(EARLY, "> $pre.early");
open(LATE, "> $pre.late");
open(MID, "> $pre.mid");

while (my $line = <SOM>) {
	next if $line =~ /^Hugo/;
	my $found = 0;
	my @f = split("\t", $line);
	#open(TMP, "> tmp.bed");
	next unless defined $f[4];
	my $chr = "chr$f[4]";
	my $snpBed = join("\t", $chr, $f[5]-1, $f[5]);
	#print TMP $snpBed;
	#print $snpBed,"\n";
	#close(TMP);
	my $sample = substr($f[15],0,12);
	
	next unless exists $purity{$sample};
	
	foreach (@files) {
		$found = 1 if $_ =~ /$sample/;
	}
	next unless $found == 1;
	my $segment = `echo "$snpBed" | $bedtools intersect -a $cnvDir/$sample.bed -b stdin | head -1`;
	next if length($segment) < 1;
	#print "match\n",$segment;
	my @s = split("\t", $segment);

	my $cov = 0;

	my $af = 0;
	if ($type eq "BRCA") {
		$af = $f[63]/100;
		$cov = 40;
	} elsif ($type eq "LUAD") {
		$af = $f[80]/($f[80]+$f[81]);
		$cov = ($f[80]+$f[81]);
	} elsif ($type eq "SKCM") {
		$af = $f[83]/($f[80]+$f[83]);
		$cov = ($f[80]+$f[83]);
	} elsif ($type eq "LIHC") {
		$af = $f[43]/($f[42]);
		$cov = ($f[42]+$f[43]);
	} elsif ($type eq "LUSC" || $type eq "PRAD" || $type eq "HNSC" || $type eq "THCA") {
		$af = $f[79]/($f[80] + $f[79]);
		$cov = ($f[80]+$f[79]);
	} elsif ($type eq "LGG") {
		$af = $f[36]/($f[37] + $f[36]);
		$cov = ($f[36]+$f[37]);
	} elsif ($type eq "BLCA") {
		$af = $f[84]/($f[85] + $f[84]);
		$cov = ($f[84]+$f[85]);
	} elsif ($type eq "UCEC") {
		$f[42] =~ s/%//;
		$af = $f[42]/100;
		$cov = 40;
	} elsif ($type eq "CESC") {
		$af = $f[59]/($f[58] + $f[59]);
		$cov = ($f[58]+$f[59]);
	} else {
		die "$!: Invalid cancer type provided.\n\n";
	}
	
	$af = $af/$purity{$sample};
	my $cn = $s[3]/$purity{$sample};
	
	if ($cn > 2.2) {
		my @cutoff = &bootstrap($af,$cov,$cn);

		if ($cutoff[0] < 0.9/$cn) {
			print LATE $line;
		} elsif ($cutoff[1] > 1.8/$cn) {
			print EARLY $line;
		} elsif ($cutoff[1] >= 0.9/$cn) {
			print MID $line;
		}
	} elsif ($cn <= 2.2 && $cn >= 1.8) {
		my @cutoff = &bootstrap($af,$cov,2);

		if ($cutoff[0] < 0.45) {
			print LATE $line;
		} elsif ($cutoff[1] >= 0.45 && $cutoff[1] <= 0.55) {
			print MID $line;
		} elsif ($cutoff[1] > 0.55) {
			print EARLY $line;
		}
	}
}


close(SOM);
close(EARLY);
close(LATE);
#print $early;

#print $late;
sub bootstrap {
	my $af = shift;
	my $cov = shift;
	my $ploidy = shift;

	my $NUM_TOSSES = $cov/2;
	my $RUNS = 10000;

	my @collect = map 0, 0 .. $NUM_TOSSES;
	my $total = 0;
	my $sumDev = 0;

	for ( 1 .. $RUNS )
	{
		my $tailsCt = 0;

		# Toss the coin NUM_TOSSES times.
		#
		for ( 1 .. $NUM_TOSSES )
		{

			$tailsCt++ if rand 2 <= $af;
		}
		$collect[ $tailsCt ] ++;
	}

	foreach my $tailsCt ( 0 .. $NUM_TOSSES )
	{
		$total += $tailsCt * $collect[ $tailsCt ]
	}

	my $mean = $af*($cov/2);

	foreach my $tailsCt ( 0 .. $NUM_TOSSES )
	{
		$sumDev += (($tailsCt- $mean)**2) * $collect[ $tailsCt ];
	}

	my $stdev = sqrt($sumDev/$RUNS);

	my $l95 = ($mean - 1.960 * ($stdev / sqrt($RUNS)))/$NUM_TOSSES;
	my $u95 = ($mean + 1.960 * ($stdev / sqrt($RUNS)))/$NUM_TOSSES;
	return($l95,$u95);
}

exit;

# plot the signatures in R

my $script = "

getAging <- function (x) {
    y <- x$raw[which((x$raw$trinuc == \"AxG\" |  x$raw$trinuc == \"CxG\" |  x$raw$trinuc == \"GxG\") & (x$raw$V6 == \"C>T\")),]
    result <- dcast(y, V11 ~ ., value.var=".", fun.aggregate=sum)
    return(result)
}

getAPOBEC <- function (x) {
    y <- x$raw[which((x$raw$trinuc == \"TxA\" |  x$raw$trinuc == \"TxT\") & (x$raw$V6 == \"C>T\")),]
    result <- dcast(y, V11 ~ ., value.var=".", fun.aggregate=sum)
    return(result)
}

getAPOBECall <- function (x) {
    y <- x$raw[which((x$raw$trinuc == \"TxA\" |  x$raw$trinuc == \"TxT\") & (x$raw$V6 == \"C>T\" | x$raw$V6 == \"C>G\")),]
    result <- dcast(y, V11 ~ ., value.var=".", fun.aggregate=sum)
    return(result)
}

getAPOBECrev1 <- function (x) {
    y <- x$raw[which((x$raw$trinuc == \"TxA\" |  x$raw$trinuc == \"TxT\") & (x$raw$V6 == \"C>G\")),]
    result <- dcast(y, V11 ~ ., value.var=".", fun.aggregate=sum)
    return(result)
}

getSmoking <- function (x) {
    y <- x$raw[which(x$raw$V6 == \"C>A\" & (x$raw$trinuc != \"TxA\" | x$raw$trinuc != \"TxC\" | x$raw$trinuc != \"TxG\" | x$raw$trinuc != \"TxT\" )),]
    result <- dcast(y, V11 ~ ., value.var=".", fun.aggregate=sum)
    return(result)
}

mutSigByIndiv <- function (x) {
    MUTS.total <- dcast(x, V11 ~ ., value.var=\"V11\", fun.aggregate=length)
    keep <- MUTS.total[which(MUTS.total$.>3),\"V11\"]
    MUTS.keep <- x[x$V11 %in% keep,]
    MUTS.d <- dcast(MUTS.keep, V11 + V5 + V6 ~ ., value.var=\"V11\", fun.aggregate=length)
    MUTS.d$prob <- apply(MUTS.d, 1, function(x) as.numeric(x["."])/sum(MUTS.d[which(MUTS.d$V11 == x[\"V11\"]),\".\"]))
    MUTS.d$trinuc <- paste(substr(MUTS.d\$V5,1,1), \"x\",substr(MUTS.d\$V5,3,3), sep=\"\")
    MUTS.empty <- expand.grid(levels(MUTS.d\$V11),\"x\",0,trint,base.subs,0)
    colnames(MUTS.empty) <- c(\"V11\",\"V5\",\".\",\"trinuc\",\"V6\",\"prob\")
    MUTS.d.2 <- dcast(rbind(MUTS.d,MUTS.empty), V11 + trinuc + V6 ~ ., value.var=\"prob\", fun.aggregate=sum)
    MUTS.d.d <- dcast(MUTS.d, trinuc + V6 ~ ., value.var = \"prob\", fun.aggregate = function(x) sum(x)/length(levels(MUTS.d\$V11)))
    MUTS.p <- ggplot(MUTS.d.d, aes(x=trinuc,y=.,fill=V6)) + facet_grid(. ~ V6) + geom_bar(stat=\"identity\") + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    result <- list(data=MUTS.d.d,plot=MUTS.p,raw=MUTS.d.2)
    return(result)
}

getSigs <- function(x,y,z) {
    aging <- getAging(x)
    apobec <- getAPOBEC(x)
    apobecrev1 <- getAPOBECrev1(x)
    apobecall <- getAPOBECall(x)
    smoking <- getSmoking(x)
    aging$hap <- y
    aging$shap <- z
    apobec$hap <- y
    apobec$shap <- z
    apobecrev1$hap <- y
    apobecrev1$shap <- z
    apobecall$hap <- y
    apobecall$shap <- z
    smoking$hap <- y
    smoking$shap <- z
    result <- list(aging=aging,apobec=apobec,apobecrev1=apobecrev1,apobecall=apobecall,smoking=smoking)
    return(result)
}

cancer.a3hI.late <- read.table(\"$late1in\",sep=\"\t\")
cancer.noa3hI.late <- read.table(\"$late2in\",sep=\"\t\")
cancer.a3hI.early <- read.table(\"$early1in\",sep=\"\t\")
cancer.noa3hI.early <- read.table(\"$early2in\",sep=\"\t\")
cancer.a3hI.mid <- read.table(\"$mid1in\",sep=\"\t\")
cancer.noa3hI.mid <- read.table(\"$mid2in\",sep=\"\t\")

cancer.noa3hI.early.results <- mutSigByIndiv(cancer.noa3hI.early)
cancer.a3hI.early.results <- mutSigByIndiv(cancer.a3hI.early)
cancer.noa3hI.mid.results <- mutSigByIndiv(cancer.noa3hI.mid)
cancer.a3hI.mid.results <- mutSigByIndiv(cancer.a3hI.mid)
cancer.noa3hI.late.results <- mutSigByIndiv(cancer.noa3hI.late)
cancer.a3hI.late.results <- mutSigByIndiv(cancer.a3hI.late)
cancer.a3hI.early.results <- mutSigByIndiv(cancer.a3hI.early)
cancer.noa3hI.mid.results <- mutSigByIndiv(cancer.noa3hI.mid)

cancer.noa3hI.late.results.sigs <- getSigs(cancer.noa3hI.late.results,\"non-A3H-I\",\"non-A3H-I\")
cancer.noa3hI.early.results.sigs <- getSigs(cancer.noa3hI.early.results,\"non-A3H-I\",\"non-A3H-I\")
cancer.a3hI.early.results.sigs <- getSigs(cancer.a3hI.early.results,\"A3H-I\",\"A3H-I\")
cancer.a3hI.late.results.sigs <- getSigs(cancer.a3hI.late.results,\"A3H-I\",\"A3H-I\")
cancer.noa3hI.mid.results.sigs <- getSigs(cancer.noa3hI.mid.results,\"non-A3H-I\",\"non-A3H-I\")
cancer.a3hI.mid.results.sigs <- getSigs(cancer.a3hI.mid.results,\"A3H-I\",\"A3H-I\")

cancer.early.smoking.p1 <- ggplot(rbind(cancer.noa3hI.early.results.sigs\$smoking,cancer.a3hI.early.results.sigs\$smoking), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()
cancer.early.apobec.p1 <- ggplot(rbind(cancer.noa3hI.early.results.sigs\$apobec,cancer.a3hI.early.results.sigs\$apobec), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()
cancer.early.aging.p1 <- ggplot(rbind(cancer.noa3hI.early.results.sigs\$aging,cancer.a3hI.early.results.sigs\$aging), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()

cancer.mid.smoking.p1 <- ggplot(rbind(cancer.noa3hI.mid.results.sigs\$smoking,cancer.a3hI.mid.results.sigs\$smoking), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()
cancer.mid.apobec.p1 <- ggplot(rbind(cancer.noa3hI.mid.results.sigs\$apobec,cancer.a3hI.mid.results.sigs\$apobec), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()
cancer.mid.aging.p1 <- ggplot(rbind(cancer.noa3hI.mid.results.sigs\$aging,cancer.a3hI.mid.results.sigs\$aging), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()

cancer.late.smoking.p1 <- ggplot(rbind(cancer.noa3hI.late.results.sigs\$smoking,cancer.a3hI.late.results.sigs\$smoking), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()
cancer.late.apobec.p1 <- ggplot(rbind(cancer.noa3hI.late.results.sigs\$apobec,cancer.a3hI.late.results.sigs\$apobec), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()
cancer.late.aging.p1 <- ggplot(rbind(cancer.noa3hI.late.results.sigs\$aging,cancer.a3hI.late.results.sigs\$aging), aes(shap, .)) + stat_summary(fun.y=mean, fun.ymin = min, fun.ymax = max, aes(fill = shap), geom=\"bar\") + stat_sum_df(\"mean_se\", geom = \"errorbar\") + stat_summary(fun.data = give.n, geom = \"text\", fun.y = mean) + scale_y_continuous(expand = c(0, 0)) + guides(fill=F) + theme_classic()

pdf(\"$out.pdf\")
multiplot(cancer.early.apobec.p1, cancer.early.smoking.p1, cancer.early.aging.p1, cancer.mid.apobec.p1, cancer.mid.smoking.p1, cancer.mid.aging.p1, cancer.late.apobec.p1, cancer.late.smoking.p1, cancer.late.aging.p1, cols=9)
dev.off()

cancer.pvalues = data.frame(
	cancer = \"cancer\",
	late.smoking = t.test(x = cancer.a3hI.late.results.sigs\$smoking\$., y = cancer.noa3hI.late.results.sigs\$smoking\$., paired = F, var.equal = F)\$p.value,
	early.smoking = t.test(x = cancer.a3hI.early.results.sigs\$smoking\$., y = cancer.noa3hI.early.results.sigs\$smoking\$., paired = F, var.equal = F)\$p.value,
	late.aging = t.test(x = cancer.a3hI.late.results.sigs\$aging\$., y = cancer.noa3hI.late.results.sigs\$aging\$., paired = F, var.equal = F)\$p.value,
	early.aging = t.test(x = cancer.a3hI.early.results.sigs\$aging\$., y = cancer.noa3hI.early.results.sigs\$aging\$., paired = F, var.equal = F)\$p.value,
	late.apobec = t.test(x = cancer.a3hI.late.results.sigs\$apobec\$., y = cancer.noa3hI.late.results.sigs\$apobec\$., paired = F, var.equal = F)\$p.value,
	early.apobec = t.test(x = cancer.a3hI.early.results.sigs\$apobec\$., y = cancer.noa3hI.early.results.sigs\$apobec\$., paired = F, var.equal = F)\$p.value
)
write.table(\"$out.txt\",sep=\"\t\")
";

system("Rscript script$n.R");