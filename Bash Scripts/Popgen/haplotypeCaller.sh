#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=HaplotypeCaller
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 08/30/21                  Last Modified: 08/30/21 ###
###########################################################################
###########################################################################
###                     HaplotypeCaller.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR

samples=( SRS5767847 SRS5767848 SRS5767850 SRS5767851 SRS5767852 SRS5767853 SRS5767854 SRS5767855 SRS5767856 SRS5767857 SRS5767859 SRS5767870 SRS5767880 SRS5767882 SRS5767883 SRS5767884 SRS5767885 )

#while read -a line

for line in "${samples[@]}"
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line}_hapCall
#SBATCH -e ${line}_hapCall
#SBATCH -o ${line}_hapCall

cd $SLURM_SUBMIT_DIR
module load gatk/4.1.2.0

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/gatk/gvcf/

gatk --java-options \"-Xmx100g -XX:+UseParallelGC -XX:ParallelGCThreads=20\" HaplotypeCaller  \
-R /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/ref/CroVir_genome_L77pg_16Aug2017.final_rename2.fasta \
-I /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/bamfiles/${line}.*.bam \
-ERC GVCF --native-pair-hmm-threads 20 \
-O ${line}.cvir.raw.variants.g.vcf" \
> /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/popgen/per_ind/${line}/${line}_HapCall.sh
done
#done < /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/all_samples.list

#while read -a line
for line in "${samples[@]}"
do
	cd /fs/ess/scratch/PAS1533/smathur/rnaseq/errors/per_ind/${line}
	sbatch /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/popgen/per_ind/${line}/${line}_HapCall.sh
done
#done < /fs/ess/scratch/PAS1533/smathur/rnaseq/popgen/reads/all_samples.list