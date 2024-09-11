#!/usr/bin/perl
#
# conver sam to bam, then sort and index 
#


use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $sam (@ARGV){
	$pm->start and next FILES; ## fork
	$sam =~ m/^([A-Za-z0-9_\-]+\.)sam/ or die "failed to match $fq\n";
	$base = $1;
	system "samtools view -b -O BAM -o $base"."bam $sam\n";
        system "samtools sort -O BAM -o $base"."sorted.bam $base"."bam\n";
        system "samtools index -b $base"."sorted.bam\n";
        $pm->finish;
}

$pm->wait_all_children;



