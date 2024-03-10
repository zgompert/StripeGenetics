# StripeGenetics
Genomic analyses of color-pattern in *Timema cristinae* with a focus on structural variation and Refugio

# Genomes

We have two new pairs of phased genomes for *T. cristinae*. Both are from Refugio, one is for a striped morph and one is for a green (unstriped) morph. This were both generated from PacBio data combined with Illumina/Omni-C.

The data are here:

Refugio striped (CEN4122) [report](https://github.com/zgompert/StripeGenetics/blob/main/CEN4122_report.html)

Haplotype 1: `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1`

Haplotype 2: `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap2`

Refugio green (CEN4120) [report](https://github.com/zgompert/StripeGenetics/blob/main/CEN4120_report.html)

Haplotype 1: `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap1`

Haplotype 2: `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap2`


# Comparative alignments

We then made synteny and alignment plots in R, see [SynPlotsRefugio.R](https://github.com/zgompert/StripeGenetics/blob/main/SynPlotsRefugioPlus.R).

# Demographic inference

# Tests of admixture and introgression

We wanted to know whether the stripe translocation introgressed from another species. We used previously published whole genome sequence data for this. This includes whole genomes from 20 *T. californicum* (SM), 20 *T. chumash* (BS), 20 *T. curi* (CR), 16 *T. knulli* (BCTUR), 19 *T. landelsensis* (BCBOG), 19 *T. poppensis* (SM), 20 *T. cristinae* from Hwy154 (HVA), and 21 *T. cristinae* from Refugio (R12A) (155 total). 

First, the DNA fastq data were aligned to the Refugio striped genome haplotype 1.

```bash
#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bwa-tcr
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load bwa 

cd /scratch/general/nfs1/u6000989/timema_wgs_abba_baba
## for cristinae
perl BwaAlignWgsFork.pl tcr_R12A*_1.fq
## for others
perl BwaAlignWgsFork2.pl 
```
The alignment script for *T. cristinae* (working form fq files directly) `BwaAlignWgsFork.pl`:

```perl
#!/usr/bin/perl
#
# alignment with bwa
#


use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta";

foreach $file1 (@ARGV){
	$pm->start and next; ## fork
	$file2 = $file1;
	$file2 =~ s/_1\.fq/_2.fq/ or die "failed file sub\n";
 	$file1 =~ m/^tcr_([A-Z0-9]+_\d+)/ or die "failed id match\n";
	$id = $1;
   
        
	system "bwa mem -t 1 -M -r 1.3 -k 19 -R \'\@RG\\tID:tcr-"."$id\\tPL:ILLUMINA\\tLB:tcr-"."$id\\tSM:tcr-"."$id"."\' -o $id".".sam $genome $file1 $file2\n";

	$pm->finish;
}

$pm->wait_all_children;
```

For the others, working from zipped files `BwaAlignWgsFork2.pl`: 

```perl
#!/usr/bin/perl
#
# bwa mem
#


use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);


my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta";

open(IN,"fqA_1.txt");
while(<IN>){
    $pm->start and next; ## fork
    chomp;
    $file1 = $_;
    $file2 = $file1;
    $file2 =~ s/_1\.fa/_2.fa/ or die "failed file sub\n";
    $file1 =~ m/fq\/([a-z]+_[A-Z]+_[A-Z]+_moe[a-zA-Z0-9]+_WTCHG_\d+)/ or die "failed id match\n";
    $id = $1;
   
    $fq1 = $file1;
    $fq2 = $file2;
    $fq1 = "t$id"."_1.fq";
    $fq2 = "t$id"."_2.fq";
    system "gunzip --stdout $file1 > $fq1\n";
    system "gunzip --stdout $file2 > $fq2\n";
    	
    system "bwa mem -t 1 -M -r 1.3 -k 19 -R \'\@RG\\tID:t"."$id\\tPL:ILLUMINA\\tLB:t"."$id\\tSM:t"."$id"."\' -o t$id".".sam $genome $fq1 $fq2\n";
    $pm->finish;
}
 

$pm->wait_all_children;
```

All alignments used `bwa mem` version (0.7.17-r1198-dirty).

The alignments were then compressed, sorted and indexed with `samtools` (version 1.16):

```bash
#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=sam2bam
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
## samtools 1.16

cd /scratch/general/nfs1/u6000989/timema_wgs_abba_baba
perl Sam2BamFork.pl *sam
```

```perl
#!/usr/bin/perl
#
# conver sam to bam, then sort and index 
#


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $sam (@ARGV){
	$pm->start and next FILES; ## fork
	$sam =~ m/^([A-Za-z0-9_\-]+\.)sam/ or die "failed to match $sam\n";
	$base = $1;
	system "samtools view -b -O BAM -o $base"."bam $sam\n";
        system "samtools sort -O BAM -o $base"."sorted.bam $base"."bam\n";
        system "samtools index -b $base"."sorted.bam\n";
        $pm->finish;
}

$pm->wait_all_children;
```

Next, PCR duplicates were masked and removed with `samtools` (version 1.16):

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=dedup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /scratch/general/nfs1/u6000989/timema_wgs_abba_baba
perl RemoveDupsFork.pl *bam
```

```perl
#!/usr/bin/perl
#
# PCR duplicate removal with samtools
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $bam (@ARGV){
	$pm->start and next FILES; ## fork
	$bam =~ m/^([A-Za-z0-9_]+)/ or die "failed to match $bam\n";
	$base = $1;
	system "samtools collate -o co_$base.bam $bam /scratch/general/nfs1/u6000989/timema_wgs_abba_baba/Dedup/t$bam\n";
	system "samtools fixmate -m co_$base.bam fix_$base.bam\n";
	system "samtools sort -o sort_$base.bam fix_$base.bam\n";
	## using default definition of dups
	## measure positions based on template start/end (default). = -m t
	system "samtools markdup -T /scratch/general/nfs1/u6000989/timema_wgs_abba_baba/Dedup/ -r sort_$base.bam dedup_$base.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```

For the non *T. cristinae*, alignments for each individual were spread over multiple (4) bam files. These were merged with `samtools` (version 1.16):

```perl
#!/usr/bin/perl
#
# merge bam alignments for each individual 
#


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $file (@ARGV){
        $pm->start and next FILES; ## fork
        system "samtools merge -r -c -p -b $file -O BAM $file.bam\n";        
        system "samtools sort -O BAM -o $file.sorted.bam $file.bam\n";
        system "samtools index $file.sorted.bam\n";
        $pm->finish;
}

$pm->wait_all_children;
```

Variants were then called with `bcftools` (version 1.16) call. This was done in parallel for each of the 13 chromosomes:

```bash
#!/bin/sh
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bcf_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /scratch/general/nfs1/u6000989/timema_wgs_abba_baba

perl BcfForkLg.pl chrom*list 
```

```perl
#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#


use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);

my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta";


foreach $chrom (@ARGV){
	$pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
	$out = "o_timema_chrom$1";
	system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o $out"."vcf\n";

	$pm->finish;
}

$pm->wait_all_children;

```

Variants were filtered with `GATK` (version 4.1.4.1) based on depth (1.5X per individual, or 233X total) and  mapping quality (30):

```perl
#!/usr/bin/perl
#
# filter vcf files 
#


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);


foreach $vcf (@ARGV){
	$pm->start and next; ## fork
	$o = $vcf;
	$o =~ s/o_// or die "failed sub $vcf\n";
	system "bgzip $vcf\n";
	system "tabix $vcf.gz\n";
	system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar IndexFeatureFile -I $vcf.gz\n";
	system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration -R /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta -V $vcf.gz -O filt_$o.gz --filter-name \"depth\" --filter-expression \"DP < 233\" --filter-name \"mapping\" --filter-expression \"MQ < 30\"\n";
	system "bgzip -d filt_$o.gz\n";

	$pm->finish;

}

$pm->wait_all_children;
```
```bash
## keeps biallelic SNPs that passed the above depth and mapping filters
grep ^Sc filt_o_timema1.vcf | grep PASS | grep -v [ATCG],[ATCG] > clean_o_timema1.vcf"
```


# GWA of stripe

# Cline analyses

A collection of notes and links that need cleaning:

* In general, cline width is proportional to sigma/sqrt(s).
* Slope from linear regression of logit(p) on location is four times the allele frequency gradient (1/width) ([Barton & Hewitt 1989](), [Kruuk et al. 1999]()).
* [Mallet et al 1990](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203983/pdf/ge1244921.pdf) provides a nice discussion of estimating dispersal (sigma) for a hybrid zone based on LD. This includes an equation for cases where the cline widths of a pair of loci are not the same, D = sigma2 / r w1 w2, note R (correlation) = D/ sqrt(p1 q1 p2 q2), so this can be re-written in terms of R.
* [Haldane 1948](https://www.ias.ac.in/article/fulltext/jgen/048/03/0277-0284) provides a model for an ecotonal cline that estimates selection (k and K) based on the interquartile range. The equations all assume some form of dominance (I think), but still give selection on par with what I get using alternatives.
* [Szymura & Barton 1986](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1558-5646.1986.tb05740.x) gives a justification for using genetic estimates of dispersal from LD, and provides an equation for doing this, it is the same as [Mallet et al 1990](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203983/pdf/ge1244921.pdf) but assumes all allele frequencies at the center are 0.5. This paper also gives an equation for wdith as width (l) = sqrt(8 sigma2)/s, but the 8 is specific to underdominance.
* [Szymura & Barton 1991](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1558-5646.1991.tb04400.x) makes similar arguments and includes the sigmoidal cline function, P = [1 + tanh(2(x-c)/w]/2]. Good one to cite for this I think.
* Need to carefully read this one in general, [Nagylaki 1974](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1213362/pdf/595.pdf).
* Cline width depends on s and sigma in slightly different ways for different types of selection. This is summarized in [Barton & Gale](https://books.google.com/books?hl=en&lr=&id=aFJFkVKskYIC&oi=fnd&pg=PA13&ots=MFk-dnK9NI&sig=r7KfAnHvJyLgyUJsMeLs7jpp33c#v=onepage&q&f=false). Let s* be the difference in mean fitness between populations at the center and edge of the zone, then w = 1.732 sigma/srt(s*) for selection favoring alternative alleles on sides of an ecotone with no dominance, 1.782 with dominance.
