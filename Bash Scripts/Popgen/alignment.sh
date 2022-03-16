#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 50:00:00
#SBATCH --job-name=alignment
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 08/07/21                  Last Modified: 08/18/21 ###
###########################################################################
###########################################################################
###                     aligment.sh                        				###
###########################################################################


cd $SLURM_SUBMIT_DIR

# step0a: Index reference and dictionary

#module load bwa/0.7.17-r1198 #bwa/0.7.17
#module load picard/2.3.0
#module load gatk/4.1.2.0
#module load samtools/1.9 

#bwa index /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta
#samtools faidx /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta

#java -jar /usr/local/picard/picard-tools-2.3.0/picard.jar CreateSequenceDictionary \
#reference=/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
#output=/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.dict

#step0b: Create directories per sample to store jobfiles and error files

#while read -a line
#do
#	mkdir /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/popgen/per_ind/${line}
#	mkdir /fs/ess/scratch/PAS1533/smathur/rnaseq/errors/per_ind/${line}
#done < /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir.list

# Step1: Align cvir samples to reference genome

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line}.align
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge
module load bwa/0.7.17-r1198 #bwa/0.7.17
module load picard/2.3.0
module load gatk/4.1.2.0
module load samtools/1.9 

# STEP_1: Align to Cvir genome 

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/alignment/step1/

bwa mem -t 50 -M -R \"@RG\tID:group1\tSM:${line}\tPL:illumina\tLB:lib1\tPU:unit1\" \
/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/trimmed/${line}.R1.paired \
/fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir/trimmed/${line}.R1.paired \
> ${line}.cvir.sam

java -jar /usr/local/picard/picard-tools-2.3.0/picard.jar ValidateSamFile \
I=${line}.cvir.sam MODE=SUMMARY O=${line}.cvir.sam.txt

java -jar /usr/local/picard/picard-tools-2.3.0/picard.jar SortSam \
INPUT=${line}.cvir.sam OUTPUT=sorted_${line}.cvir.bam SORT_ORDER=coordinate

java -jar /usr/local/picard/picard-tools-2.3.0/picard.jar MarkDuplicates \
INPUT=sorted_${line}.cvir.bam OUTPUT=dedup_${line}.cvir.bam METRICS_FILE=metrics_${line}.cvir.bam.txt

java -jar /usr/local/picard/picard-tools-2.3.0/picard.jar BuildBamIndex \
INPUT=dedup_${line}.cvir.bam

# STEP_2: Fix mate pair info in BAM

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/alignment/step2/

java -jar /usr/local/picard/picard-tools-2.3.0/picard.jar FixMateInformation \
INPUT=../step1/dedup_${line}.cvir.bam \
OUTPUT=${line}.cvir.fixmate.bam \
SO=coordinate \
CREATE_INDEX=true

# STEP_3: Get mapping stats

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/alignment/stats/

samtools depth -a ../step2/${line}.cvir.fixmate.bam \
| awk '{c++;s+=\$3}END{print s/c}' \
> ${line}.cvir.meandepth.txt

samtools depth -a ../step2/${line}.cvir.fixmate.bam \
| awk '{c++; if(\$3>=10) total+=1}END{print (total/c)*100}' \
> ${line}.cvir.10xbreadth.txt

samtools flagstat ../step2/${line}.cvir.fixmate.bam \
> ${line}.cvir.mapstats.txt

samtools stats ../step2/${line}.cvir.fixmate.bam \
> ${line}.cvir.stats.txt" \
> /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/popgen/per_ind/${line}/${line}.align.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir_SRA.list


# Submit alignment jobs

while read -a line
do
	cd /fs/ess/scratch/PAS1533/smathur/rnaseq/errors/per_ind/${line}
	sbatch /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/popgen/per_ind/${line}/${line}.align.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/cvir_SRA.list
