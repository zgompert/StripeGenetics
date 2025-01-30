#!/usr/bin/perl
#
# generates a table of GO terms per gene ID

# usage: out_GO_clean_genes* uni_GO_genes_*
#

$in1 = shift(@ARGV); ## go terms

open(IN, $in1) or die "failed to read $in1\n";

while(<IN>){
	chomp;
	m/(GO:\d+)\s+(\S+)/ or die "failed here: $_\n";
	$go = $1;
	$type = $2;
	m/D\d+\s(.+)$/ or die "failed term: $_\n";
	$term = $1;
	$gos{$go}{$type} = $term;
}

close(IN);


$in2 = shift(@ARGV); ## genes with GO numbers
open(IN, $in2) or die "failed to read $in2\n";
open(OUT, "> table_$in2") or die "failed to write table for $in2\n";
while(<IN>){
	chomp;
	m/(\d+)\s+(\S+)/ or die "what is this? $_\n";
	$gene = $1;
	$go = $2;
	foreach $type (sort keys %{$gos{$go}}){
		print OUT "$gene $go $type $gos{$go}{$type}\n";
	}
}
close(IN);
close(OUT);
