#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 14-00:00:00
#SBATCH --job-name=genotypeGVCF_cvir
#SBATCH --mem=248G
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/07/21                  Last Modified: 09/07/21 ###
###########################################################################
###########################################################################
###                     genotypeGVCF.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load bioinfo
module load GATK/4.1.3.0 #gatk/4.1.3.0
module load vcftools

## Step 1: CombineGVCFs 


cd /scratch/bell/mathur20/osu/rnaseq/popgen/gatk/gvcfs/

gatk --java-options "-Xmx248g -XX:+UseParallelGC -XX:ParallelGCThreads=128" CombineGVCFs  \
-R /scratch/bell/mathur20/osu/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
--variant SRS5767847.cvir.raw.variants.g.vcf \
--variant SRS5767848.cvir.raw.variants.g.vcf \
--variant SRS5767849.cvir.raw.variants.g.vcf \
--variant SRS5767850.cvir.raw.variants.g.vcf \
--variant SRS5767851.cvir.raw.variants.g.vcf \
--variant SRS5767852.cvir.raw.variants.g.vcf \
--variant SRS5767853.cvir.raw.variants.g.vcf \
--variant SRS5767854.cvir.raw.variants.g.vcf \
--variant SRS5767855.cvir.raw.variants.g.vcf \
--variant SRS5767856.cvir.raw.variants.g.vcf \
--variant SRS5767857.cvir.raw.variants.g.vcf \
--variant SRS5767859.cvir.raw.variants.g.vcf \
--variant SRS5767870.cvir.raw.variants.g.vcf \
--variant SRS5767880.cvir.raw.variants.g.vcf \
--variant SRS5767881.cvir.raw.variants.g.vcf \
--variant SRS5767882.cvir.raw.variants.g.vcf \
--variant SRS5767883.cvir.raw.variants.g.vcf \
--variant SRS5767884.cvir.raw.variants.g.vcf \
--variant SRS5767885.cvir.raw.variants.g.vcf \
-O all.cvir.g.vcf.gz

# Step2: GenotypeGVCF

cd /scratch/bell/mathur20/osu/rnaseq/popgen/gatk/variants/

gatk --java-options "-Xmx248g -XX:+UseParallelGC -XX:ParallelGCThreads=128" GenotypeGVCFs \
-R /scratch/bell/mathur20/osu/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
-V /scratch/bell/mathur20/osu/rnaseq/popgen/gatk/gvcfs/all.cvir.g.vcf.gz \
-O all.cvir.raw.variants.vcf.gz

# Step3: Select SNPs

gatk --java-options "-Xmx248g -XX:+UseParallelGC -XX:ParallelGCThreads=128" SelectVariants \
-R /scratch/bell/mathur20/osu/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
-V all.cvir.raw.variants.vcf.gz \
--select-type-to-include SNP \
-O all.cvir.raw.SNPs.vcf.gz 

#Filter SNPs
# From grossen et al 2020 QD < 2.0, FS > 40.0, SOR > 5.0, MQ< 20.0, −3.0 > MQRandkSum >3.0, 
# −3.0 > ReadPosRankSum >3.0 and AN < 62 (80% of all Alpine ibex individuals)

# using AN < 78 (80% of all samples (2X49 alleles))

gatk --java-options "-Xmx248g" VariantFiltration \
-R /scratch/bell/mathur20/osu/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
-V all.cvir.raw.SNPs.vcf.gz  \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "MQ < 20.0" --filter-name "MQ20" \
-filter "MQRankSum < -3.0" --filter-name "MQRankSum-3.0" \
-filter "MQRankSum > 3.0" --filter-name "MQRankSum3.0" \
-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3.0" \
-filter "ReadPosRankSum > 3.0" --filter-name "ReadPosRankSum3.0" \
-filter "SOR > 5.0" --filter-name "SOR5" \
-filter "FS > 40.0" --filter-name "FS40" \
-filter "AN < 78.0" --filter-name "AN78" \
-filter "AF < 0.05" --filter-name "MAF0.05" \
-O all.cvir.filterflag.SNPs.vcf.gz

gatk SelectVariants \
-R /scratch/bell/mathur20/osu/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
-V all.cvir.filterflag.SNPs.vcf.gz \
-select 'vc.isNotFiltered()' \
-O all.cvir.filtered.SNPs.vcf.gz


# keep only biallelic sites

vcftools --gzvcf all.cvir.filtered.SNPs.vcf.gz \
--recode --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 \
--out all.cvir.final.SNPs


