#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=adapter_removal
#SBATCH -e adapter_removal
#SBATCH -o adapter_removal

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 01/26/21                  Last Modified: 04/08/21 ###
###########################################################################
###########################################################################
###                     adapter_removal.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR

#Step0: Create folders for jobcodes and errors

#while read -a line
#do
#	mkdir /scratch/bell/mathur20/osu/rnaseq/jobcodes/per_ind/${line}
#	mkdir /scratch/bell/mathur20/osu/rnaseq/errors/per_ind/${line}
#done < /scratch/bell/mathur20/osu/rnaseq/data/samples.list

#run for each sample
while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line[0]}_adptrem
#SBATCH -e ${line[0]}_adptrem
#SBATCH -o ${line[0]}_adptrem

cd $SLURM_SUBMIT_DIR

#module load trimmomatic

### Step1: Adapter removal 

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/trimmed/

java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar \
PE /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/raw/${line[0]}.R1.fastq /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/raw/${line[0]}.R2.fastq \
${line[0]}.R1.fastq.paired ${line[0]}.R1.fastq.unpaired \
${line[0]}.R2.fastq.paired ${line[0]}.R2.fastq.unpaired \
LEADING:20 TRAILING:20 MINLEN:30 -threads 20 \
ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/TruSeq2-PE.fa:2:40:10" \
>/fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/per_treatment/${line[0]}_adptrem.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/treatment.list



#Submit jobs #

while read -a line
do
	cd /fs/ess/scratch/PAS1533/smathur/rnaseq/errors/per_treatment/
	sbatch  /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/per_treatment/${line[0]}_adptrem.sh
done </fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/treatment.list


