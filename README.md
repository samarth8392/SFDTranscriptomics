# SFD Transcriptomics 🐍
>This Github repository contains scripts developed for the transcriptomic study of host response to snake fungal disease (ophidiomycosis) in different temperature conditions

**Publication:**

Mathur, S., Haynes, E., Allender, M. C., & Gibbs, H. L. (2024). Genetic mechanisms and biological processes underlying host response to ophidiomycosis (snake fungal disease) inferred from tissue‐specific transcriptome analyses. *Molecular Ecology*, 33(2), e17210. https://doi.org/10.1111/mec.17210

## Folders ##

### The directory structure of the repository is:
> - ***_Bash Scripts_***: Contains scripts for data analysis on OSC clusters.
>> _Popgen_ : Scripts for population genomics of _C.viridis_ genomes download from NCBI SRA database
>> _RNASeq_ : Scripts for RNASeq data processing and mapping
> - ***_RScripts_***: Contains scripts for statistical analysis and output data visualization in R.
> - ***_DataFies_***: Contains miscellenious datafiles used/generated for the study
***_Bash Scripts_***

- RNASeq
1. adapter_removal.sh
> To remove adapter and low base quality sequences from raw RNASeq reads using triommomatic.

2. fastqc.sh
> To quality check filtered RNASeq reads

3. alignment.sh
> To build reference genome database and map filtered RNAseq reads using hisat2

4. stringtie.sh
> To assemble transcripts for each sample using _C. viridis_ reference annotation file to guide assembly and merged sample transcripts 

5. featureCounts.sh
> To generate a transcript count matrix 

- Popgen
1. adapter_removal.sh
> To remove adapter and low base quality sequences from raw whole genome sequence reads using triommomatic.

2. alignment.sh
> To align filtered genomic reads to reference genome. Script also sort, mark duplications, fix mate pair information, and get mapping stats.

3. base_recalibration.sh
> To identify known variants from our dataset and recalibrate base quality scores for genomic reads for each sample

- RScripts

1. Differential gene expression (DGE) using DESeq2
> DGE_liver.R : For liver samples
> DGE_kidney.R : For kidney samples
> DGE_skin.R : For skin samples

2. DGE_Compare.R
> To compare DGE results for each tissue and extract differentially expressed genes (DEGs) due to fixed effects only

3. DGE_Contrast_Compare.R
> To compare DGE results for each tissue and extract differentially expressed genes (DEGs) due to interaction between infection and temperature

4. Weighted gene co-expression network analysis (WGCNA)
> WGCNA_liver.R : For liver samples
> WGCNA_kidney.R : For kidney samples
> WGCNA_skin.R : For skin samples

5. Compare_WGCNA.R
> To compare WGCNA results for each tissue and extract genes within each module significantly associated with infection or the interaction between infection and temperature.

6. gene_list.R
> To get the final gene lists from DGE and WGCNA analysis.
