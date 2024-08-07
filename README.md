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

We used RepeatMasker (version 4.0.7) and the *T. cristinae* repeat library to mask repeates for each genome.

```{bash}
#!/bin/sh 
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=wolf-kp
#SBATCH --partition=wolf-kp
#SBATCH --job-name=repeatm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load repeatmasker

#version 4.0.7
cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/repeat_mask

## run repeat masker on each genome sequence, uses library from the 2020 Science paper 
## developed by Victor

RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta
RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap2/ojincantatabio-cen4122-hap2-mb-hirise-14fv0__08-10-2023__final_assembly.fasta
RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap1/ojincantatabio-cen4120-hap1-mb-hirise-wlbll__08-15-2023__final_assembly.fasta
RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap2/ojincantatabio-cen4120-hap2-mb-hirise-bn0ko__08-15-2023__final_assembly.fasta
```

# Comparative alignments
We used `cactus` for align the new genomes to each other and to our existing striped and green genomes (Hi-C/Omni-C, but not phased, from Hwy154).

```{bash}
#!/bin/sh 
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=cactus-master
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

cd /scratch/general/nfs1/cactusNp

module load cactus

cactus jobStoreGUSR1_GUSR2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGUSR1_TcrGUS2.txt cactusTcrGUSR1_TcrGUSR2.hal --maxCores 80 
cactus jobStoreGSR1_GSR2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGSR1_TcrGS2.txt cactusTcrGSR1_TcrGSR2.hal --maxCores 80 

cactus jobStoreGSR1_GUSR1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGSR1_TcrGUS1.txt cactusTcrGSR1_TcrGUSR1.hal --maxCores 80 
cactus jobStoreGSR1_GUSR2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGSR1_TcrGUS2.txt cactusTcrGSR1_TcrGUSR2.hal --maxCores 80 
cactus jobStoreGSR2_GUSR2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGSR2_TcrGUS2.txt cactusTcrGSR2_TcrGUSR2.hal --maxCores 80 
cactus jobStoreGSR2_GUSR1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGSR2_TcrGUS1.txt cactusTcrGSR2_TcrGUSR1.hal --maxCores 80 

cactus jobStoreGSR1_GSM /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGSR1_TcrMainGS.txt cactusTcrGSR1_TcrGS.hal --maxCores 80 
cactus jobStoreGUSR1_GUSM /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusStripe_TcrGUSR1_TcrMainGUS.txt cactusTcrGUSR1_TcrGUS.hal --maxCores 80 
```
Next, we used `halSyntency` to extract synteny blocks from the genome alignments.

```{bash}
#!/bin/sh 
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=cactus-syn
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

cd /scratch/general/nfs1/cactusNp

perl forkHal.pl *hal
```

```{perl}
#!/usr/bin/perl
#
# hal syntency batch run 
#


use Parallel::ForkManager;
my $max = 10;
my $pm = Parallel::ForkManager->new($max);


@target = ("t_cris_r_gs1","t_cris_r_gs1","t_cris_r_gs1","t_cris_r_gs1","t_cris_r_gs2","t_cris_r_gs2","t_cris_r_gus1","t_cris_r_gus1");
@query = ("t_cris_gs","t_cris_r_gs2","t_cris_r_gus1","t_cris_r_gus2","t_cris_r_gus1","t_cris_r_gus2","t_cris_gus","t_cris_r_gus2");

foreach $hal (@ARGV){
	$g1 = shift(@target);
	$g2 = shift(@query);
	$pm->start and next;
	$out = "out_$hal";
	$out =~ s/hal/psl/ or die;
	system "~/source/hal/bin/halSynteny --queryGenome $g2 --targetGenome $g1 $hal $out\n"; 
	$pm->finish;
}


$pm->wait_all_children;
```

We then made synteny and alignment plots in R, see [SynPlotsRefugio.R](https://github.com/zgompert/StripeGenetics/blob/main/SynPlotsRefugioPlus.R).

# Genome annotations

We annotated the eight *Timema* genomes using [braker3](https://github.com/Gaius-Augustus/BRAKER) (structural annotation) and [interproscan5](https://academic.oup.com/bioinformatics/article/30/9/1236/237988?login=false) (functional annotation). 

Annotations are in /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/ with one subdirectory for each genome. 

Structural annotation with `braker` (version 3.0.8) used protein and RNA sequence data. The protain data comprised 2,601,995 arthropod proteins from OrthDB (arthropod data set version 10 (odb10_arthropoda_fasta.tar.gz) (the input file is named proteins.fasta; this is the same input used for the annotation of *T. knulli* in [Nosil et al. 2023](https://www.pnas.org/doi/abs/10.1073/pnas.2300673120)). The RNA data are *T. cristinae* (N = 18) from 2017, which were originally published as part of a study of DNA methylation, see [de Carvalho et al. 2023](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.17165) ([BioProject PRJNA1010130](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1010130)). It took some work to get this to run on the cluster (see comments in script for things that had to be done initially to make things work). Here is an example (for the main mountain green haplotype 1) of the script I used for the run. This includes a BUSO assessment (based on insecta_odb10) of the inferred genes and of the original genome. All of the results are in the braker subdirectory for each genome.

```{bash}
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=braker
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

source ~/.bashrc

ml braker/3.0.8
ml busco
#mkdir ~/augustus
#containerShell
#cp /usr/share/augustus/config ~/augustus
#exit

#cd
#mv ~/augustus/config ~/augustus/config-old
#ml braker/3.0.8
#containerShell
#cp -r /opt/Augustus/config ~/augustus/
#exit

cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h1

## run braker
braker.pl --genome=../../t_crist_gus_hap_cen4280/HiRise/Hap1/ojincantatabio-cen4280-hap1-mb-hirise-ig5ps__01-30-2024__hic_output.fasta.masked --prot_seq=../proteins.fasta --rnaseq_sets_ids=tcr135.17_0003_R,tcr137.17_0006_R,tcr139.17_0012_R,tcr140.17_0015_R,tcr141.17_0019_R,tcr142.17_0043_R,tcr143.17_0045_R,tcr144.17_0049_R,tcr145.17_0051_R,tcr146.17_0057_R,tcr148.17_0062_R,tcr149.17_0065_R,tcr150.17_0067_R,tcr151.17_0070_R,tcr152.17_0074_R,tcr173.17_0075_R,tcr174.17_0081_R,tcr175.17_0082_R --rnaseq_sets_dirs=/scratch/general/nfs1/u6000989/timema_rna_2020/ --AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts --AUGUSTUS_CONFIG_PATH=/uufs/chpc.utah.edu/common/home/u6000989/augustus/config --threads=48 --gff3

## run busco, genome and aa
cd braker
## genome
busco -i ../../../t_crist_gus_hap_cen4280/HiRise/Hap1/ojincantatabio-cen4280-hap1-mb-hirise-ig5ps__01-30-2024__hic_output.fasta.masked -m geno -o busco_genome_out -l insecta_odb10

## amino acids
busco -i braker.aa -m prot -o busco_aa_out -l insecta_odb10
```
We have very high BUSCO scores for the genomes, and more modest scores for the found genes, which is to be expected (braker is not designed to maximize BUSCOs, but rather to identify genes with high confidence) (more details are in busco_genome_out and busco_aa_out).

| Genome | Genome % complete | AA % complete |
|--------|-------------------|---------------|
| t_crist_hyw154_green_h1 | 99.1 | 67.0 |
| t_crist_hyw154_green_h2 | 99.5 | 70.0 |
| t_crist_refug_green_h1 | 99.4 | 70.0 |
| t_crist_refug_green_h2 | 99.2 | 70.9 |
| t_crist_refug_stripe_h1 | 99.2 | 69.0 |
| t_crist_refug_stripe_h2 | 90.3 | 65.9 |

I am using `interproscan` (version 5.63) for functional annotation. This is being run on the amino acid sequence predictions from `braker` (the output files require a small bit of reformatting to get rid of the stop codons). Here is the version of the script for one of the genomes.

```{bash}
#!/bin/bash 
#SBATCH --time=140:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=braker
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

source ~/.bashrc

module load interproscan

cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h1/braker

cat braker.aa | perl -p -i -e 's/\*//g' > braker_aa.fasta

/uufs/chpc.utah.edu/sys/installdir/r8/interproscan/5.63/interproscan.sh -cpu 48 -i braker_aa.fasta -goterms -b interpro_aa
````

# Demographic inference

We used `moments` to fit and compare four historical demographic models for *T. cristinae* Refugio versus Hwy 154.

```{bash}
#!/bin/sh 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=wolf-kp
#SBATCH --partition=wolf-kp
#SBATCH --job-name=momments
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load miniconda3/latest

cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/timema_SV_balance/demog

max=16
count=0

for mfile in *m_sub_mom*py;
do
    python3 ${mfile} &
    ((count++))

    if ((count >= max)); then
        wait -n
        ((count--))
    fi
done

wait
```

[Strict isolation = SI](https://github.com/zgompert/StripeGenetics/blob/main/si_sub_moments.py)

[Isolation-with-migration = IM](https://github.com/zgompert/StripeGenetics/blob/main/im_sub_moments.py)

[Ancient migration = AM](https://github.com/zgompert/StripeGenetics/blob/main/am_sub_moments.py)

[Secondary contact = SC](https://github.com/zgompert/StripeGenetics/blob/main/sc_sub_moments.py)


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

We conducted a GWA mapping analysis based on 238 *T. cristinae* from Refugio (all collected in 2016). The fastq files are in `/uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/reads_refugio`. 

We aligned these data to haplotype 1 of both the green and striped genome from Refugio using `bwa` (version 0.7.17-r1198-dirty).

```perl
#!/usr/bin/perl
#
# alignment with bwa
#

use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

## green refugio hap 1
#my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap1/ojincantatabio-cen4120-hap1-mb-hirise-wlbll__08-15-2023__final_assembly.fasta";

## stripe refugio hap 1
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta";

FILES:
foreach $fq (@ARGV){
        $pm->start and next FILES; ## fork
        if ($fq =~ m/(16_[a-zA-Z0-9_]+)/){
                $ind = $1;
        }
        else {
                die "Failed to match $file\n";
        }
        system "bwa aln -n 4 -l 20 -k 2 -t 1 -q 10 -f aln"."$ind".".sai $genome $fq\n";
        system "bwa samse -n 1 -r \'\@RG\\tID:tcr-"."$ind\\tPL:ILLUMINA\\tLB:tcr-"."$ind\\tSM:tcr-"."$ind"."\' -f aln"."$ind".".sam $genome aln"."$ind".".sai $fq\n";
        $pm->finish;
}

$pm->wait_all_children;

```
We then  compressed, sorted and indexed the alignment files with `samtools` (version 1.16). Here is the version for the striped genome (green is bascially the same thing).

```bash
#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=gompert
#SBATCH --partition=notchpeak
#SBATCH --job-name=sam2bam
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
## samtools 1.16

cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/align_refugio_gs
perl Sam2BamFork.pl *sam

```

```perl
#!/usr/bin/perl
#
# conver sam to bam, then sort and index 
#


use Parallel::ForkManager;
my $max = 16;
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

````

We then used `samtools` and `bcftools` (version 1.16 for both) for variant calling (again with each genome).

```bash
#!/bin/sh 
#SBATCH --time=196:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=wolf-kp
#SBATCH --partition=wolf-kp
#SBATCH --job-name=bcfCall
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load bcftools
module load samtools
## bcftools 1.16
## samtools 1.16


cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/align_refugio_gs

## index genome just needs to be done once
#samtools faidx /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta

bcftools mpileup -b bams -C 50 -d 500 -f /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta -q 20 -Q 30 -I -Ou -a DP,AD,ADF,ADR | bcftools call -v -c -p 0.01 -P 0.001 -O v -o tcr_refugio_variants_gs.vcf 
```
I used our standard perl scripts (vcfFilter.pl) for filtering the vcf files. The one for the striped genome is. 

```perl
#!/usr/bin/perl

use warnings;
use strict;

# this program filters a vcf file based on overall sequence coverage, number of non-reference reads, number of alleles, and reverse orientation reads

# usage vcfFilter.pl infile.vcf
#
# change the marked variables below to adjust settings
#

#### stringency variables, edits as desired
## 238 inds, 2x
my $minCoverage = 476; # minimum number of sequences; DP
my $minAltRds = 10; # minimum number of sequences with the alternative allele; AC
my $notFixed = 1.0; # removes loci fixed for alt; AF
my $bqrs = 0.005; # p-value base quality rank sum test; BaseQRankSum
my $mqrs = 0.005; # p-value mapping quality rank sum test; MQRankSum
my $rprs = 0.005; # p-value read position rank sum test; ReadPosRankSum
my $qd = 2; # minimum ratio of variant confidenct to non reference read depth; QD
my $mq = 30; # minimum mapping quality; MQ
my $miss = 47; # maximum number of individuals with no data
##### this set is for GBS
my $d;

my @line;

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z_0-9\-]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> filtered2x_$1") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
        chomp;
        $flag = 1;
        if (m/^\#/){ ## header row, always write
                $flag = 1;
        }
        elsif (m/^Sc/){ ## this is a sequence line, you migh need to edit this reg. expr.
                $flag = 1;
                $d = () = (m/\d\/\d:0,0,0:0/g); ## for bcftools call
                if ($d >= $miss){
                        $flag = 0;
                        ##print "fail missing : ";
                }
                if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
                        $flag = 0;
                        #print "fail allele : ";
                }
                @line = split(/\s+/,$_);
                if(length($line[3]) > 1 or length($line[4]) > 1){
                        $flag = 0;
                        #print "fail INDEL : ";
                }
                m/DP=(\d+)/ or die "Syntax error, DP not found\n";
                if ($1 < $minCoverage){
                        $flag = 0;
                        #print "fail DP : ";
               }
## bcftools call version

                m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error DP4 not found\n";
                if(($1 + $2) < $minAltRds){
                        $flag = 0;
                }
                m/AF1*=([0-9\.e\-]+)/ or die "Syntax error, AF not found\n";
                if ($1 == $notFixed){
                        $flag = 0;
                #       print "fail AF : ";
                }

## bcftools call verions, these are p-values, use 0.01
                if(m/BQB=([0-9e\-\.]*)/){
                        if ($1 < 0.005){
                                $flag = 0;
#                               print "fail BQRS : ";
                        }
                }
                if(m/MQB=([0-9e\-\.]*)/){
                        if ($1 < 0.005){
                                $flag = 0;
#                               print "fail MQRS : ";
                        }
                }
                if(m/RPB=([0-9e\-\.]*)/){
                        if ($1 < 0.005){
                                $flag = 0;
#                               print "fail RPRS : ";
                        }
                }
                if(m/MQ=([0-9\.]+)/){
                        if ($1 < $mq){
                                $flag = 0;
#                               print "fail MQ : ";
                        }
                }
                else{
                        $flag = 0;
                        print "faile no MQ : ";
                }
                if ($flag == 1){
                        $cnt++; ## this is a good SNV
                }
        }
        else{
                print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
                $flag = 0;
        }
        if ($flag == 1){
                print OUT "$_\n";
        }
}
close (IN);
close (OUT);                                             
``
I then applied the depth filter from `filterSomeMore.pl`.

```perl
#!/usr/bin/perl

use warnings;
use strict;

# filter vcf files


### stringency variables, edits as desired
my $maxCoverage =  1814; # maximum depth to avoid repeats, mean + 2sd

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z0-9_]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> morefilter_$1") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
        chomp;
        $flag = 1;
        print "\n";
        if (m/^\#/){ ## header row, always write
                $flag = 1;
        }
        elsif (m/^Sc/){ ## this is a sequence line, you migh need to edit this reg. expr.
                $flag = 1;
                m/DP=(\d+)/ or die "Syntax error, DP not found\n";
                if ($1 > $maxCoverage){
                        $flag = 0;
                        print "fail DP\n";
                }
                if ($flag == 1){
                        $cnt++; ## this is a good SNV
                }
        }
        else{
                print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
                $flag = 0;
        }
        if ($flag == 1){
                print OUT "$_\n";
        }
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
```

Lastly for variant generation, I converted the vcf to gl format and applied a minor allele frequency (0.01) filter:

```bash
perl vcf2glSamt.pl 0.01 morefilter_filtered2x_tcr_refugio_variants.vcf
```

```perl
#!/usr/bin/perl

my @line = ();
my $word;
my $nind = 0;
my $nloc = 0;
my $first = 1; ## first vcf file, get ids from here
my $out = "filtered_tcr_refugio_variants_gs.gl";

open (OUT, "> $out") or die "Could not write the outfile\n";

if ($out =~ s/gl/txt/){
	open (OUT2, "> af_$out") or die "Count not write the alt. af file\n";
}

my $maf = shift (@ARGV);

foreach my $in (@ARGV){
	open (IN, $in) or die "Could not read the vcf file\n";
	while (<IN>){
		chomp;
		## get individual ids
		if (m/^#CHROM/ & ($first == 1)){
			@line = split(m/\s+/, $_);	
			foreach $word (@line){
				if ($word =~ m/tcr/){
					push (@inds, $word);
					$nind++;
				}
			}
			print OUT "$nind $nloc\n";
			$word = join (" ", @inds);
			print OUT "$word\n";
		}
		## read genetic data lines, write gl
		elsif (m/^Sclu3Hs_(\d+);HRSCAF_\d+\s+(\d+)/){
			$word = "$1".":"."$2";
			if (m/AF1=([0-9\.\-e]+)/){
				$palt = $1;
				if ($palt > 0.5){
					$p = 1 - $palt;
				}
				else {
					$p = $palt;
				}
				#print "$word = $p\n";
				if ($p >= $maf){ ## keep this locus
					$nloc++;
					print OUT "$word";
					@line = split(m/\s+/, $_);
					$i = 0;
					foreach $word (@line){
						if ($word =~ m/^\d\/\d\:/){
							$word =~ m/(\d+),(\d+),(\d+)/ or die "failed match $word \n";
							print OUT " $1 $2 $3";
						}
				
					}
					print OUT "\n";
					print OUT2 "$palt\n"; ## print p before converting to maf
				}
			}
		}	
		
	}
}
close (OUT);
close (OUT2);
```
This results in 85,558 SNPs and 238 individuals for the striped genome and 87,202 SNPs and 238 individuals for the green genome. 
 
Next, I used `entropy` (version 1.2) to estimate genotypes. This was done with each genome and with 2 or 3 source populations (5 chains each). Initial values for MCMC were generated with [initq.R](initq.R), based on simple genotype point estimates from [gl2genest.pl](gl2genest.pl) and [runEstP.pl](runEstP.pl).

```bash
#!/bin/sh 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=entropy
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module purge
module load gcc/8.5.0 hdf5/1.10.7

## entropy 1.2
#~/bin/entropy

cd /uufs/chpc.utah.edu/common/home/gompert-group4/projects/timema_SV_balance/gen_refugio
perl forkEntropy.pl
```
```perl
#!/usr/bin/perl
# fork script for entropy 

use Parallel::ForkManager;
my $max = 30;
my $pm = Parallel::ForkManager->new($max);


foreach $k (2..3){
	foreach $ch (0..4){ 
		foreach $g ("gs","gus"){
			sleep 2;
			$pm->start and next;
			$out = "out_$g"."_k$k"."_ch$ch".".hdf5";
			system "entropy -i filtered_tcr_refugio_variants_$g".".gl -w 0 -m 1 -l 2000 -b 1000 -t 5 -k $k -o $out -q $g"."_ldak$k".".txt -s 20\n";
			$pm->finish;
		}
	}
}
$pm->wait_all_children;
```

I prepared the genotypic and phenotypic data for mapping with [FormatGeno.pl](FormatGeno.pl) and [FormatPheno.R](FormatPheno.R), respectively.

I used the LMM in gemma (version 0.95a) for genome-wide association mapping of stripe (and color) using haplotype 1 from both the green and striped genomes and with and without including chromosome 8 for the kinship matrix. Here are the full set of commands I used:

```bash
#!/bin/sh 

module load gemma
# 0.95a

## pattern gs, kinship includes ch8
gemma -g  pattern_g_tcr_refugio_gs.geno -p ph_pattern.txt -gk 1 -o o_ref_pattern_gs -maf 0
gemma -g  pattern_g_tcr_refugio_gs.geno -p ph_pattern.txt -k output/o_ref_pattern_gs.cXX.txt -lmm 4 -n 1 -o o_ref_pattern_gs -maf 0

## color gs, kinship includes ch8
gemma -g  color_g_tcr_refugio_gs.geno -p ph_color.txt -gk 1 -o o_ref_color_gs -maf 0
gemma -g  color_g_tcr_refugio_gs.geno -p ph_color.txt -k output/o_ref_color_gs.cXX.txt -lmm 4 -n 1 -o o_ref_color_gs -maf 0

## pattern gs, kinship excludes ch8
gemma -g  pattern_no8_g_tcr_refugio_gs.geno -p ph_pattern.txt -gk 1 -o o_ref_pattern_no8_gs -maf 0
gemma -g  pattern_g_tcr_refugio_gs.geno -p ph_pattern.txt -k output/o_ref_pattern_no8_gs.cXX.txt -lmm 4 -n 1 -o o_ref_pattern_no8_gs -maf 0

## color gs, kinship excludes ch8
gemma -g  color_no8_g_tcr_refugio_gs.geno -p ph_color.txt -gk 1 -o o_ref_color_no8_gs -maf 0
gemma -g  color_g_tcr_refugio_gs.geno -p ph_color.txt -k output/o_ref_color_no8_gs.cXX.txt -lmm 4 -n 1 -o o_ref_color_no8_gs -maf 0

## pattern gus, kinship includes ch8
gemma -g  pattern_g_tcr_refugio_gus.geno -p ph_pattern.txt -gk 1 -o o_ref_pattern_gus -maf 0
gemma -g  pattern_g_tcr_refugio_gus.geno -p ph_pattern.txt -k output/o_ref_pattern_gus.cXX.txt -lmm 4 -n 1 -o o_ref_pattern_gus -maf 0

## color gus, kinship includes ch8
gemma -g  color_g_tcr_refugio_gus.geno -p ph_color.txt -gk 1 -o o_ref_color_gus -maf 0
gemma -g  color_g_tcr_refugio_gus.geno -p ph_color.txt -k output/o_ref_color_gus.cXX.txt -lmm 4 -n 1 -o o_ref_color_gus -maf 0

## pattern gus, kinship excludes ch8
gemma -g  pattern_no8_g_tcr_refugio_gus.geno -p ph_pattern.txt -gk 1 -o o_ref_pattern_no8_gus -maf 0
gemma -g  pattern_g_tcr_refugio_gus.geno -p ph_pattern.txt -k output/o_ref_pattern_no8_gus.cXX.txt -lmm 4 -n 1 -o o_ref_pattern_no8_gus -maf 0

## color gus, kinship excludes ch8
gemma -g  color_no8_g_tcr_refugio_gus.geno -p ph_color.txt -gk 1 -o o_ref_color_no8_gus -maf 0
gemma -g  color_g_tcr_refugio_gus.geno -p ph_color.txt -k output/o_ref_color_no8_gus.cXX.txt -lmm 4 -n 1 -o o_ref_color_no8_gus -maf 0
```
I then summarized the results (for pattern) in R with [summarize_gemma_gs.R](summarize_gemma_gs.R) and [summarize_gemma_gus.R](summarize_gemma_gus.R).

# Cline analyses

We fit geographic clines for stripe (phenotype and SV genotype) and background SNPs. This all uses the same genotype estiamtes (from `entropy`) used for the GWA mapping.

The core analyses are all described in [ClineAnalyses.R](ClineAnalyses.R). Stan model specifications are in [cline.stan](cline.stan), [cline_con.stan](cline_con.stan) and [llcline.stan](llcline.stan).

General notes that were useful for developing the ideas in the code are below.

* In general, cline width is proportional to sigma/sqrt(s).
* Slope from linear regression of logit(p) on location is four times the allele frequency gradient (1/width) ([Barton & Hewitt 1989](), [Kruuk et al. 1999]()).
* [Mallet et al 1990](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203983/pdf/ge1244921.pdf) provides a nice discussion of estimating dispersal (sigma) for a hybrid zone based on LD. This includes an equation for cases where the cline widths of a pair of loci are not the same, D = sigma2 / r w1 w2, note R (correlation) = D/ sqrt(p1 q1 p2 q2), so this can be re-written in terms of R.
* [Haldane 1948](https://www.ias.ac.in/article/fulltext/jgen/048/03/0277-0284) provides a model for an ecotonal cline that estimates selection (k and K) based on the interquartile range. The equations all assume some form of dominance, but still give selection on par with what I get using alternatives.
* [Szymura & Barton 1986](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1558-5646.1986.tb05740.x) gives a justification for using genetic estimates of dispersal from LD, and provides an equation for doing this, it is the same as [Mallet et al 1990](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203983/pdf/ge1244921.pdf) but assumes all allele frequencies at the center are 0.5. This paper also gives an equation for wdith as width (l) = sqrt(8 sigma2)/s, but the 8 is specific to underdominance.
* [Szymura & Barton 1991](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1558-5646.1991.tb04400.x) makes similar arguments and includes the sigmoidal cline function, P = [1 + tanh(2(x-c)/w]/2]. 
* Cline width depends on s and sigma in slightly different ways for different types of selection. This is summarized in [Barton & Gale](https://books.google.com/books?hl=en&lr=&id=aFJFkVKskYIC&oi=fnd&pg=PA13&ots=MFk-dnK9NI&sig=r7KfAnHvJyLgyUJsMeLs7jpp33c#v=onepage&q&f=false). Let s* be the difference in mean fitness between populations at the center and edge of the zone, then w = 1.732 sigma/srt(s*) for selection favoring alternative alleles on sides of an ecotone with no dominance, 1.782 with dominance.

# Comparative functional genomics

I am using a number of analyses for functional characterization of the SVs based on the comparative alignments and annotations. 

I am using `GENESPACE` (version 1.3.1) (see [Lovell et al. 2022](https://elifesciences.org/articles/78526) and [this GitHub page](https://github.com/jtlovell/GENESPACE)) to look at syteny in terms of protein/gene sequences. Here is code for an analysis in progress (I need to update this with new versions of files).

```{bash}
module load orthofinder
```

```{R}
library(GENESPACE)
###############################################
genomeRepo<-"/scratch/general/nfs1/u6000989/rawGenomes"
wd<-"/scratch/general/nfs1/u6000989/"
path2mcscanx<-"~/bin/"

# -- download raw data from NCBI for human and chicken genomes UPDATE
dir.create(genomeRepo)
rawFiles <- download_exampleData(filepath = genomeRepo)

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("human", "chicken"),
  genomeIDs = c("human", "chicken"),
  presets = "ncbi",
  genespaceWd = wd)

# -- initalize the run and QC the inputs


gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = path2mcscanx)

## need to set this
gpar$shellCalls$orthofinder<-"orthofinder"

# -- accomplish the run
out <- run_genespace(gpar)
```
