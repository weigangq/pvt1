#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

my @sampleList;
my @seqs;
while (<>) {
    chomp;
    my @data = split;
    my $ref = shift @data;
    my $alt = shift @data;
    if (/^REF/) { @sampleList = @data; next }

    for (my $i=0; $i<=$#sampleList; $i++){
        my @arr = split /\|/, $data[$i];
        my @ar = map {$_? $alt : $ref} @arr;
        if (!$seqs[$i]){ push @seqs, \@ar; next }
        my @seq = @{$seqs[$i]};
        $seqs[$i] = [$seq[0].$ar[0], $seq[1].$ar[1]]
    }
}
#print Dumper(\@seqs); exit;

for (my $i=0; $i<=$#sampleList; $i++){
    my $sample = $sampleList[$i];
    my @seq = @{$seqs[$i]};
    for (my $j=0; $j<2; $j++){
        print '>', $sample, '|', $j+1, "\n", $seq[$j], "\n"
    }
}
