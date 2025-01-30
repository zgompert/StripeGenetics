#!/usr/bin/perl
#

$min = 3262; ## includes bounds, SV is 3263-3346 = 84 genes
$max = 3345;

open(IN, "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h2/braker/interpro_aa.gff3") or die "failed to read\n";
open(OUT, "> GO_genes_h154_gus_h2.txt") or die "failed to write\n";

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
