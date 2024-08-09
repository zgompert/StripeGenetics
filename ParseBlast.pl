#!/usr/bin/perl
#
# parse local blast
# for each query report gene, length of gene, best hist

## get file from command line
$in = shift(@ARGV);
open(IN, $in) or die "failed to open $in for reading";
open(OUT, "> summary_$in") or die "failed to write outfile\n";

## read file
$reading = 0;
while(<IN>){
	chomp;
	if(m/^Query= *(\S+)/){
		$gene = $1;
		$reading = 1;
	} elsif ($reading==1 and m/Length=(\d+)/){
		$length = $1;
	} elsif ($reading==1 and m/^(Scaffold\S+)\s+\S+\s+(\S+)/){
		print OUT "$gene $length $1 $2\n";
		$reading = 0;
	} elsif ($reading==1 and m/No hits found/){
		print OUT "$gene $length NA NA\n";
		$reading = 0;
	}
}
close(IN);
close(OUT);
