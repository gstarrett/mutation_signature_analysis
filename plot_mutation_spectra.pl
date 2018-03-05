#!/usr/bin/perl -w
use strict;

my $file = shift;

my $n = time;

my $Rscript = "
library(\"reshape2\")
library(\"ggplot2\")
sample <- read.table(\"$file\", sep =\"\t\")
sample\$trint <- paste(substr(sample\$V5,1,1),\"x\",substr(sample\$V5,3,3),sep=\"\")
sample.cast <- dcast(sample,trint + V6 ~ .)
colnames(sample.cast) <- c(\"trint\",\"mut\",\"count\")
sample.cast\$per <- sample.cast\$count/sum(sample.cast\$count)
pdf(\"$file.all.pdf\", width=8, height=4)
ggplot(sample.cast,aes(trint, per)) + geom_bar(aes(fill=mut),stat = \"identity\") + facet_grid(. ~ mut) + theme_classic() + ylim(c(0,0.2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=6))
dev.off()
pdf(\"$file.C.pdf\", width=8, height=4)
ggplot(sample.cast[grepl(sample.cast\$mut,pattern = \"^C>[ATCG]\"),],aes(trint, per)) + geom_bar(aes(fill=mut),stat = \"identity\") + theme_classic() + scale_fill_manual(values=c(\"blue\",\"black\",\"red\")) + ylim(c(0,0.2))
dev.off()
";
open(OUT, "> script$n.R");
print OUT $Rscript;
close(OUT);
system("Rscript script$n.R");
