use strict;
use warnings;

my @file;
while (<>) {
	chomp;
	push @file, $_
}

for (my $i=0; $i<$#file; $i++) {
	for (my $j=$i+1; $j<=$#file; $j++) {
		my $f1 = $1 if $file[$i]=~/(.+)\.+/;
		my $f2 = $1 if $file[$j]=~/(.+)\.+/;
		print 'vcftools --weir-fst-pop ', $file[$i], ' --weir-fst-pop ', $file[$j], ' --vcf pvt1.recode.vcf --out ', $f1, '-', $f2, '  --remove-indels', "\n"
	}
}
