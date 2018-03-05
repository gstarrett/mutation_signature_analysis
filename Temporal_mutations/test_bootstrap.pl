#!/usr/bin/perl -w
use strict;

my $af = .60;
my $cov = 80;

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
print "$l95 ... $u95\n";