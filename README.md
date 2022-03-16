# SFD Transcriptomics
Scripts developed for manuscript **_"Genetic mechanisms and biological processes underlying host response to ophiodiomycosis (Snake Fungal Disease) inferred from tissue-specific transcriptomes"._**

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

3. alignment.sh
> To align filtered genomic reads to reference genome. Script also sort, mark duplications, fix mate pair information, and get mapping stats.

