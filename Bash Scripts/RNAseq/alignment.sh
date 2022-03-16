#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=aligment.noSE
#SBATCH -e %x_%j
#SBATCH -o %x_%j

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 01/27/21                  Last Modified: 01/27/21 ###
###########################################################################
###########################################################################
###                     aligment.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR

#module load hisat2/2.1.0
#module load samtools

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/

#hisat2-build -p 64 \
#/fs/ess/scratch/PAS1533/smathur/rnaseq/ref/cvir/CroVir_genome_L77pg_16Aug2017.final_rename.fasta \
#cvir_wholegenome

## map reads to reference genome
while read -a line
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line[0]}.align
#SBATCH -e %x_%j
#SBATCH -o %x_%j

cd $SLURM_SUBMIT_DIR
module load hisat2/2.1.0
module load samtools

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/

hisat2 --dta --time --threads 64 \
-s ${line[0]}.cvir.wgr.sam \
-x ../cvir_wholegenome \
-1 /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/trimmed/${line[0]}.R1.fastq.paired \
-2 /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/trimmed/${line[0]}.R2.fastq.paired \
> ${line[0]}.cvir.wgr.sam


samtools view -@ 2 -Su ${line[0]}.cvir.wgr.sam | samtools sort -@ 2 -o ${line[0]}.sorted.bam" > /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/per_treatment/${line[0]}.align.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/treatment.list


#Submit jobs

while read -a line
do
	cd /fs/ess/scratch/PAS1533/smathur/rnaseq/errors/per_treatment/
	sbatch /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/per_treatment/${line[0]}.align.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/treatment.list

# Sort SAM


#while read -a line
#do
#	samtools view -@ 2 -Su ${line[0]}.cvir.wgr.noSE.sam | samtools sort -@ 2 -o ${line[0]}.sorted.noSE.bam
#done < /fs/ess/scratch/PAS1533/smathur/rnaseq/data/treatment.list