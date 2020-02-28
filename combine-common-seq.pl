#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

my %seqs;
my $id;
while (<>) {
    chomp;
    if (/^>(.+)/) { $id = $1 }
    else {
        my $s = $_;
        my @ids;
        @ids = @{$seqs{$s}} if $seqs{$s};
        push @ids, $id;
        $seqs{$s} = \@ids
    }
}
#print Dumper(\%seqs); exit;

open(OUT1,'>', 'group.txt');
open(OUT2,'>', 'common.fas');

my $n = 1;
foreach (keys %seqs) {
    print OUT2 '>g_', $n, "\n", $_, "\n";
    my @ids = @{$seqs{$_}};
    print OUT1 $_, "\t", 'g_', $n, "\n" foreach @ids;
    $n++
}

close OUT1;
close OUT2;
