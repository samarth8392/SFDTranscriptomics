#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --job-name=snpeff
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/13/21                  Last Modified: 09/13/21 ###
###########################################################################
###########################################################################
###                     variant_annotation.sh                        	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load snpEff/4.3

# SNPeff

# Steps:
# 1. Copy snpeff.config file from /group/bioinfo/apps/apps/snpEff-4.3/ to local directory
# 2. Make a directory in your local directory called "data" and another directory within data with your species name. I call it cvir (C. viridis)
# 3. copy ref fasta and annotation gff as sequences.fa and genes.gff into "cvir" folder
# 4. Change the following in config file:
# 	(a) Add in the first two lines:
#			# C.viridis genome
#			cvir.genome : C. viridis
#	(b) Change data.dir to /scratch/bell/mathur20/osu/rnaseq/popgen/snpeff/data/
# 5. Build the database

#cd /scratch/bell/mathur20/osu/rnaseq/popgen/snpeff/
#snpEff build -c snpEff.config -gff3 -v cvir &> build.logfile.txt

# If the database builds without error, you should see snpEffectPredictor.bin inside your cvir folder

# Step6: Annotate your variants
cd /scratch/bell/mathur20/osu/rnaseq/popgen/snpeff/

snpEff ann -stats -c snpEff.config \
-no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v cvir \
/scratch/bell/mathur20/osu/rnaseq/popgen/vcf/all.cvir.final.wholegenome.SNPs.vcf \
> SNPeff.all.final.vcf

