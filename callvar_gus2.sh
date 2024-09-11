#!/bin/sh 
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert
#SBATCH --partition=notchpeak
#SBATCH --job-name=bcfCall
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load bcftools
module load samtools
## bcftools 1.16
## samtools 1.16


cd /scratch/general/nfs1/u6000989/t_cris_fha
## index genome just needs to be done once
#samtools faidx /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gs_hap_cen4119/HiRise/Hap1/final_assembly.fasta
#samtools faidx /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap2/ojincantatabio-cen4280-hap2-mb-hirise-i2xb7__01-30-2024__hic_output.fasta

bcftools mpileup -b bams_gus2 -C 50 -d 500 -f /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap2/ojincantatabio-cen4280-hap2-mb-hirise-i2xb7__01-30-2024__hic_output.fasta -q 20 -Q 30 -I -Ou -a DP,AD,ADF,ADR | bcftools call -v -c -p 0.01 -P 0.001 -O v -o tcr_h154_variants_gus.vcf 
