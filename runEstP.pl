#!/usr/bin/perl
#
# run estpEM for each population
#

foreach $in (@ARGV){
	$in =~ m/^([A-Z]+)/ or die "failed here: $in\n";
	$out = "p_$1.txt";
	system "estpEM -i $in -o $out -e 0.001 -m 50 -h 2\n";
}
close(IN);
