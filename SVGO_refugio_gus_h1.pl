#!/usr/bin/perl
#

$min = 4059; ## includes bounds, SV is 4060 to 4302
$max = 4303;

open(IN, "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h1/braker/interpro_aa.gff3") or die "failed to read\n";
open(OUT, "> GO_genes_refug_gus_h1.txt") or die "failed to write\n";

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
