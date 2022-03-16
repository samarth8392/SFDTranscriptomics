#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=stringtie_merge
#SBATCH -e %x_%j
#SBATCH -o %x_%j

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 03/29/21                  Last Modified: 04/08/21 ###
###########################################################################
###########################################################################
###                     stringtie.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR


# create assembly per sample
cd /fs/ess/scratch/PAS1533/smathur/rnaseq/stringtie/


while read -a line
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line[0]}.stringtie
#SBATCH -e %x_%j
#SBATCH -o %x_%j

cd $SLURM_SUBMIT_DIR

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/stringtie/r2/

/fs/ess/scratch/PAS1533/smathur/softwares/stringtie/stringtie \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/${line[0]}.sorted.bam \
-l ${line[0]} -p 20 \
-G /fs/ess/scratch/PAS1533/smathur/rnaseq/ref/cvir/CroVir_rnd1.all.maker.final.homologIDs.gff \
-o ${line[0]}.sorted.bam.gtf" > /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/per_treatment/${line[0]}.stringtie.sh
done < /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/treatment.list

#while read -a line
#do
#	cd /fs/ess/scratch/PAS1533/smathur/rnaseq/errors/per_treatment/
#	sbatch /fs/ess/scratch/PAS1533/smathur/rnaseq/jobcodes/per_treatment/${line[0]}.stringtie.sh
#done < /fs/ess/scratch/PAS1533/smathur/rnaseq/data/r2/treatment.list

## once all sample assemblies have been generated:
# merge all transcripts from all samples + reference annotation

/fs/ess/scratch/PAS1533/smathur/softwares/stringtie/stringtie --merge -p 20 \
-G /fs/ess/scratch/PAS1533/smathur/rnaseq/ref/cvir/CroVir_rnd1.all.maker.final.homologIDs.gff \
-o stringtie_all_merged.gtf merge_all

# compare assembled transcripts to known transcripts

/fs/ess/scratch/PAS1533/smathur/softwares/gffcompare/gffcompare -r \
/fs/ess/scratch/PAS1533/smathur/rnaseq/ref/cvir/CroVir_rnd1.all.maker.final.homologIDs.gff \
-G -o stringtie_all_merged stringtie_all_merged.gtf
