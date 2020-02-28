#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my %collect;
while (<>) {
	chomp;
	my @data = split;
	my $type = $data[5] eq 'archaic'? 'A' : 'a';
	$collect{$data[4]}{$data[0]}{$type}++
}
#print Dumper(\%collect); exit;

my %counts;
foreach my $pop (keys %collect){
	my %population = %{$collect{$pop}};
	foreach (keys %population){
		my %obj = %{$population{$_}};
		my @arr = keys %obj;
		my $type = @arr==1? $arr[0].$arr[0] : 'Aa';
		$counts{$pop}{$type}++
	}
}
#print Dumper(\%counts); exit;

print "pop\tAA\tAa\taa\tAA(exp)\tAa(exp)\taa(exp)\n";
foreach my $pop (sort keys %counts){
	my %c = %{$counts{$pop}};
#	print join "\t", ($pop, $c{'AA'}? $c{'AA'} : 0, $c{'Aa'}? $c{'Aa'} : 0, $c{'aa'}? $c{'aa'} : 0, &hardy);
	my @obs = ($c{'AA'}? $c{'AA'} : 0, $c{'Aa'}? $c{'Aa'} : 0, $c{'aa'}? $c{'aa'} : 0);
	print join ($pop, @obs, &hardy(@obs));
	print "\n";
}
