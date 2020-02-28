#!/usr/bin/perl -w
use strict;
#use JSON;
use Data::Dumper;

my ($totalFile, $popFile) = @ARGV;

open(totalF, $totalFile);
my %allele;
while (<totalF>) {
	chomp;
	next if /CHROM/;
	my @data=split;
	my @f1 = split /:/, $data[4];
	my @f2 = split /:/, $data[5];
	$allele{$data[1]} = $f1[1]<$f2[1]? 4 : 5
}
close totalF;
#print Dumper(\%allele);exit;

my @pop = split /\./, $popFile;
print $pop[1];

open(popF, $popFile);
while (<popF>) {
	chomp;
	next if /CHROM/;
	my @data=split;
	my @f = split /:/, $data[$allele{$data[1]}];
	print "\t", $f[1]
}
close totalF;
print "\n";
