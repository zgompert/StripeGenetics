#!/usr/bin/perl
#

$min = 887; ## includes bounds, SV is 888 to 935
$max = 936;

open(IN, "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h1/braker/interpro_aa.gff3") or die "failed to read\n";
open(OUT, "> GO_genes_h154_gs_h1.txt") or die "failed to write\n";

while(<IN>){
	chomp;
	if(m/^g(\d+)/){
		$gid = $1;
		if(($gid >= $min) and ($gid <= $max)){
	                while(s/(GO:\d+)//){
				print OUT "$gid $1\n";
			}
		}
	}
}
close(IN);
close(OUT);
