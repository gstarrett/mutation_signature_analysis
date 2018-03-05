#!/usr/bin/perl -w
use strict;

my $infile = shift;

open (IN, "< $infile");

my @lines = <IN>;

for (my $i=1; $i<$#lines; $i++) {
	if ($lines[$i] =~ /^#/ || $lines[$i] !~ /LOH/) {
		print $lines[$i];
		next;
	}
	my @prev = split("\t", $lines[$i-1]);
	my @curr = split("\t",$lines[$i]);
	my @next = split("\t",$lines[$i+1]);
	
	if ($curr[0] ne $prev[0] || $curr[0] ne $next[0]) {
		print join("\t",@curr);
		next;
	}
	
	my $prevd = $curr[1]-$prev[1];
	my $nextd = $next[1]-$curr[1];
	
	if ($prev[12] eq "Somatic" && $next[12] eq "Somatic") {
		$curr[12] = "Somatic_candidate";
	} elsif ($prevd > 50000 && $nextd > 50000) {
		$curr[12] = "Somatic_candidate";
	} elsif (($prev[12] eq "Somatic" && 3*$prevd < $nextd) || ($next[12] eq "Somatic" && 3*$nextd < $prevd)) {
		$curr[12] = "Somatic_candidate";
	} 
		
	print join("\t", @curr);
		
	#print "$curr[0]\t$curr[1]\t$prev[12]\t$prevd\t$next[12]\t$nextd\n";
	
	# if prev or next is somatic within 50000bp and LOH is more than 50000bp
	
}

