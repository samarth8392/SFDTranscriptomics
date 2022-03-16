#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=adapter_removal
#SBATCH -e adapter_removal
#SBATCH -o adapter_removal

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 08/12/21                  Last Modified: 08/14/21 ###
###########################################################################
###########################################################################
###                     fastqc.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR


# FASTQC and adapter removal on raw Cvir reads

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=${line}.fastqc
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
#module load fastqc/0.11.8
module load trimmomatic

# Step1: Fastqc on raw reads

#cd /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/fastqc/
#mkdir ${line}

#fastqc -o ${line}/ \
#/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/${line}_1.fastq \
#/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/${line}_2.fastq

#cd /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/fastqc/${line}
#unzip *zip

# Step2: adapter removal

java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar \
PE /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/raw/${line}_1.fastq /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/raw/${line}_2.fastq \
/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/trimmed/${line}.R1.paired /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/trimmed/${line}.R1.unpaired \
/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/trimmed/${line}.R2.paired /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/trimmed/${line}.R2.unpaired \
LEADING:20 TRAILING:20 MINLEN:30 -threads 20 \
ILLUMINACLIP:/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/fastqc/NexteraPE-PE.fa:2:40:10" \
> /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/popgen/per_ind/${line}_fastqc.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir_SRA.list

# Submit jobs

while read -a line
do
	cd /fs/ess/scratch/PAS1533/smathur/rnaseq/errors/per_ind/${line}
	sbatch /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/popgen/per_ind/${line}_fastqc.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir_SRA.list
