#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 300:00:00
#SBATCH --job-name=depth_fungus
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/20/21                  Last Modified: 09/23/21 ###
###########################################################################
###########################################################################
###                     aligment.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR

#cd /scratch/bell/mathur20/osu/rnaseq/fungus/align/hisat/
#hisat2-build -p 64 \
#/scratch/bell/mathur20/osu/rnaseq/ref/ophidiomyces_ophidiicola/GCA_002167195.1_ASM216719v1_genomic.fna \
#oo_fungus

module load bioinfo
#module load bwa
module load samtools

#cd /scratch/bell/mathur20/osu/rnaseq/fungus/ref/
#bwa index GCA_002167195.1_ASM216719v1_genomic.fna
#samtools faidx GCA_002167195.1_ASM216719v1_genomic.fna




## map reads to fungal reference genome using Hisat
while read -a line
do
	echo "#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 300:00:00
#SBATCH --job-name=${line}_hisat
#SBATCH -e /scratch/bell/mathur20/osu/rnaseq/fungus/align/hisat/log/%x
#SBATCH -o /scratch/bell/mathur20/osu/rnaseq/fungus/align/hisat/log/%x

cd $SLURM_SUBMIT_DIR
module load bioinfo
module load hisat2/2.1.0

cd /scratch/bell/mathur20/osu/rnaseq/fungus/align/hisat/
hisat2 --dta --time --threads 8 \
-s ${line}.fungus.wgr.sam \
-x oo_fungus \
-1 /scratch/bell/mathur20/osu/rnaseq/data/after/final/${line}.paired.R1.fastq \
-2 /scratch/bell/mathur20/osu/rnaseq/data/after/final/${line}.paired.R2.fastq \
> ${line}.fungus.hisat.sam" \
> /scratch/bell/mathur20/osu/rnaseq/jobcodes/rnaseq/per_treatment/${line}/${line}_align_hisat.sh
done < /scratch/bell/mathur20/osu/rnaseq/data/after/treatment.txt

## map reads to fungal reference genome using BWA
while read -a line
do
	echo "#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 300:00:00
#SBATCH --job-name=${line}_bwa
#SBATCH -e ${line[0]}_bwa
#SBATCH -o ${line[0]}_bwa

cd $SLURM_SUBMIT_DIR
module load bioinfo
module load bwa
module load samtools

cd /scratch/bell/mathur20/osu/rnaseq/fungus/align/bwa/

bwa mem -t 50 -M -R \"@RG\tID:group1\tSM:${line}\tPL:illumina\tLB:lib1\tPU:unit1\" \
/scratch/bell/mathur20/osu/rnaseq/fungus/ref/GCA_002167195.1_ASM216719v1_genomic.fna \
/scratch/bell/mathur20/osu/rnaseq/data/after/final/${line}.paired.R1.fastq \
/scratch/bell/mathur20/osu/rnaseq/data/after/final/${line}.paired.R2.fastq \
> ${line}.fungus.bwa.sam

samtools flagstat ${line}.fungus.bwa.sam > ${line}.fungus.mapstats.txt" \
> /scratch/bell/mathur20/osu/rnaseq/jobcodes/rnaseq/per_treatment/${line}/${line}_align_bwa.sh
done < /scratch/bell/mathur20/osu/rnaseq/data/after/treatment.txt

#Submit jobs #

#while read -a line
#do
#	cd /scratch/bell/mathur20/osu/rnaseq/errors/per_treatment/${line}
#	sbatch /scratch/bell/mathur20/osu/rnaseq/jobcodes/rnaseq/per_treatment/${line}/${line}_align_hisat.sh
#	sbatch /scratch/bell/mathur20/osu/rnaseq/jobcodes/rnaseq/per_treatment/${line}/${line}_align_bwa.sh
#done < /scratch/bell/mathur20/osu/rnaseq/data/after/treatment.txt


# Get depth statistics for the Infected skin at low temp

cd /scratch/bell/mathur20/osu/rnaseq/fungus/align/stats/

for i in CV01_F_20C_Inf CV03_F_20C_Inf CV16_M_20C_Inf CV21_M_20C_Inf
do
	samtools view -@ 2 -Su /scratch/bell/mathur20/osu/rnaseq/fungus/align/hisat/${i}_Skin.fungus.hisat.sam | samtools sort -@ 2 -o ${i}_Skin.fungus.sorted.bam
	samtools depth -a ${i}_Skin.fungus.sorted.bam > ${i}_Inf_Skin.fungus.depth.txt
done