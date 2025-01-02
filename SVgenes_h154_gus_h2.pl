#!/usr/bin/perl
#

$min = 3260; ## includes bounds, SV is 3261-3318 = 58 genes
$max = 3319;

open(IN, "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h2/braker/interpro_aa.gff3") or die "failed to read\n";
open(OUT, "> SV_genes_h154_gus_h2.txt") or die "failed to write\n";

while(<IN>){
	chomp;
	if(m/^g(\d+)/){
		$gid = $1;
		if(($gid >= $min) and ($gid <= $max)){
			if(m/signature_desc=([a-zA-Z0-9 \-]+);/){
				print OUT "$gid $1\n";
			}
		}
	}
}
close(IN);
close(OUT);
