#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=featureCounts
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

# -p = is paired-end; -B = require that fragments have both ends successfully aligned; -C
# = chimeric fragments (ends mapped to 2 different chroms) not counted

cd /fs/ess/scratch/PAS1533/smathur/rnaseq/featurecounts/

/fs/ess/scratch/PAS1533/smathur/softwares/subread-2.0.1-source/bin/featureCounts \
-p -B -C -T 20 -t exon -g transcript_id \
-a /fs/ess/scratch/PAS1533/smathur/rnaseq/stringtie/stringtie_all_merged.gtf \
-o featureCount_counts_all.txt \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV01_F_20C_Inf_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV01_F_20C_Inf_Kidney.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV01_F_20C_Inf_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV03_F_20C_Inf_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV03_F_20C_Inf_Kidney.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV03_F_20C_Inf_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV05_F_20C_Con_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV05_F_20C_Con_Kidney.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV05_F_20C_Con_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV06_F_26C_Inf_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV06_F_26C_Inf_Kidney.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV06_F_26C_Inf_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV13_M_20C_Con_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV13_M_20C_Con_Kidney.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV13_M_20C_Con_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV16_M_20C_Inf_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV16_M_20C_Inf_Kidney.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV16_M_20C_Inf_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV17_M_26C_Con_Kidney.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV17_M_26C_Con_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV17_M_26C_Con_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV18_M_26C_Inf_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV18_M_26C_Inf_Kidney.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV19_M_26C_Inf_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV19_M_26C_Inf_Kidney.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV19_M_26C_Inf_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV21_M_20C_Inf_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV21_M_20C_Inf_Kidney.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV21_M_20C_Inf_Skin.sorted.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV24_F_26C_Con_Kidney.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r1/CV24_F_26C_Con_Liver.sorted.noSE.bam \
/fs/ess/scratch/PAS1533/smathur/rnaseq/align/r2/CV24_F_26C_Con_Skin.sorted.bam 