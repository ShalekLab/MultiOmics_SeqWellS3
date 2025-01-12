# Processing Multi-omic data generated utilizing SeqWell S3
Workflow scripts for data processing of multi-omics data generated using Seq-Well S3 and analysis scripts for downstream data analysis as presented in _A scalable, low-cost, sample hashing workflow for multi-omic single-cell analysis using the Seq-Well S3 platform_
Data pre-processing wrappers should be run prior to analysis scripts to ensure all necessary files are available. 

## data_preprocessing_wrappers
Wrapper scripts to aggregate and parallelize multi-step processes required to achieve analysis ready files of genetic hashing, antibody hashing, and cell surface protein quantifications. 

gatk_preprocessing.sh: wrapper script for the systematic running of main steps required for GATK pipelines. Input should be an unmapped bam files from each sample specific bulk-RNA sequencing. If starting with fastq files, use Picard's FastqToSam to convert to unmapped bam file format first.
haplotype_caller.sh: wrapper script for the systematic running of HaplotypeCaller pipeline steps. Input should be the unmapped bam files that have already undergone the GATK pre-processing steps. 
run_CITESeq_Counts.sh: wrapper for running multiple samples through CITE-seq-Count. Input is a list of sample names and a list of ID numbers which give fastq files their unique naming property. 

## analysis_scripts 

Script for the processing of scRNA-seq data and integration of genetic hashing, antibody hashing, and cell surface protein data. Scripts should be run based on the chronological order of the script titles. 
