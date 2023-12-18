#!/bin/bash

##### ANTIBODY DEMULTIPLEXING ######
##

# Overview: 
# This workflow will allow you to perform demultiplexing on samples that have been hashed using antibody tags. 
# It assumes you are using Seq-Well + Biolegend or Seq-Well + BD for generating multiplexed, single-cell data. 
# This workflow can also be used to demultiplex CITE-Seq data, just swap the HTO.csv file with an ADT.csv file which would contain barcodes for surface protein markers. 

#

# Understanding this file: 
# Distinct sections are contained within sets of two hashes (##) 
# Separation within sections, usually between explanations and scripts occur between single (#)

##

# transfer STARSolo/CellRanger output to VM alongside HTO data and sample sheet in .csv with antibody indices 
# create the HTO.csv file using your favourite bash editor (I like nano, but some prefer vim), if you try to do this with excel, 
# it probably won't work due to "extra characters", but do whatever works for you
# You can create a single HTO.csv file for all arrays, or separate files for each array. If you have a unique file for each array, name these files HTO_"sample_name".csv


run_CITESeq_Counts() {
    sample=$1
    idx=$2

    celln=$(cat "./starsolo/${sample}/Solo.out/Gene/filtered/barcodes.tsv" | wc -l)

    # Determine the HTO file to use
    if [ -f "HTO_${sample}.csv" ]; then
        hto_file="HTO_${sample}.csv"
    elif [ -f "HTO.csv" ]; then
        hto_file="HTO.csv"
    else
        echo "Error: HTO file not found for sample ${sample}."
        return 1
    fi

    nohup CITE-seq-Count -R1 "fastq/HTO_${sample}_S${idx}_R1_001.fastq.gz" \
    -R2 "fastq/HTO_${sample}_S${idx}_R2_001.fastq.gz" \
    -t "$hto_file" \
    -cbf 1 -cbl 12 \
    -umif 13 -umil 20 \
    -cells $celln \
    -wl "./starsolo/${sample}/Solo.out/Gene/filtered/barcodes.tsv" -o "HTO_${sample}" > "HTO_${sample}.out" &
}

# Convert string arguments to arrays
IFS=' ' read -r -a samples <<< "$1"
IFS=' ' read -r -a idxs <<< "$2"

# Check if the number of samples and indices match
if [ "${#samples[@]}" -ne "${#idxs[@]}" ]; then
    echo "Error: The number of samples and indices must match."
    exit 1
fi

for i in "${!samples[@]}"; do
    run_CITESeq_Counts "${samples[i]}" "${idxs[i]}"
done


# Output of this code provides files needed for Seurat analysis
