#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --job-name=getSites_DEGbyPop
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/13/21                  Last Modified: 09/13/21 ###
###########################################################################
###########################################################################
###                     variant_paritioning.sh              			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools

#No. of sites:
# Whole genome = 7,454,238
# Genic = 2,089,052 (No. of genes = 18,538)
# Intergenic = 5,365,186
# Exonic = 196,617 (No. of exons = 157,177)
# Synonymous = 67,323
# Nonsynonymous = 44,517
# Missense = 43,929
# Nonsense = 588
# DEGs = 103,463 (No. of DEGs = 842)
# DEGs_exons = 10,374
# DEGs_Nonsyn = 2,361

cd /scratch/bell/mathur20/osu/rnaseq/popgen/vcf/

# Genes

#vcftools --vcf all49.cvir.wholegenome.SNPs.vcf \
#--recode --recode-INFO-all \
#--bed /scratch/bell/mathur20/osu/rnaseq/popgen/bedfiles/cvir_justgenes.bed \
#--out all.cvir.genic.SNPs

# Intergenic

#vcftools --vcf all49.cvir.wholegenome.SNPs.vcf \
#--recode --recode-INFO-all \
#--exclude-bed /scratch/bell/mathur20/osu/rnaseq/popgen/bedfiles/cvir_justgenes.bed \
#--out all.cvir.intergenic.SNPs

# Exons

#vcftools --vcf all49.cvir.wholegenome.SNPs.vcf \
#--recode --recode-INFO-all \
#--bed /scratch/bell/mathur20/osu/rnaseq/popgen/bedfiles/cvir_justexons.bed \
#--out all.cvir.exonic.SNPs

# Exonic mutations

#for i in synonymous nonsynonymous missense nonsense
#do
#	vcftools --vcf all.cvir.wholegenome.SNPs.vcf \
#	--recode --recode-INFO-all \
#	--positions /scratch/bell/mathur20/osu/rnaseq/popgen/bedfiles/snps/all.cvir.$i.sites \
#	--out all.cvir.$i.SNPs
#done

#### By differentially Expressed Genes (DEGs) ####

# Exons
cd /scratch/bell/mathur20/osu/rnaseq/popgen/vcf/

#vcftools --vcf all.cvir.exonic.SNPs.vcf \
#--recode --recode-INFO-all \
#--bed /scratch/bell/mathur20/osu/rnaseq/popgen/degs/genelists/unique_DEGs.bed \
#--out all.cvir.DEGs_exons.SNPs

# Nonsynonymous

vcftools --vcf all.cvir.nonsynonymous.SNPs.vcf \
--recode --recode-INFO-all \
--bed /scratch/bell/mathur20/osu/rnaseq/popgen/degs/genelists/unique_DEGs.bed \
--out all.cvir.DEGs_nonsyn.SNPs
