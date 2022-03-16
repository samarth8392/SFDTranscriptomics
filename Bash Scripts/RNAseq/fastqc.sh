#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=fastqc
#SBATCH -e %x_%j
#SBATCH -o %x_%j

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 01/26/21                  Last Modified: 04/08/21 ###
###########################################################################
###########################################################################
###                     fastqc.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load fastqc/0.11.8


cd /fs/ess/scratch/PAS1533/smathur/rnaseq/fastqc/trimmed/
while read -a line
do
	mkdir ${line[0]}
	fastqc -o ${line[0]}/ \
	/fs/ess/scratch/PAS1533/smathur/rnaseq/data/final/paired/${line[0]}.R1.fastq.paired \
	/fs/ess/scratch/PAS1533/smathur/rnaseq/data/final/paired/${line[0]}.R2.fastq.paired 
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/data/final/treatment.all.list