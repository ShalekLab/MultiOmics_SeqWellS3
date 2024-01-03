#!/bin/bash

##### ANTIBODY DEMULTIPLEXING ######
##

# Overview: 
# This workflow will allow you to perform demultiplexing on samples that have been hashed using antibody tags. 
# It assumes you are using Seq-Well + Biolegend or Seq-Well + BD for generating multiplexed, single-cell data. 
# This workflow can also be used to demultiplex CITE-Seq data, just swap the HTO.csv file with an ADT.csv file which would contain barcodes for surface protein markers. 

##

# Understanding this file: 
# Distinct sections are contained within sets of two hashes (##) 
# Separation within sections, usually between explanations and scripts occur between single (#)

##

# transfer STARSolo or CellRanger output to VM alongside HTO data and sample sheet in .csv with antibody indices 
# create the HTO.csv file using your favourite bash editor (I like nano, but some prefer vim), if you try to do this with excel, 
# it probably won't work due to "extra characters", but do whatever works for you
# You can create a single HTO.csv file for all arrays, or separate files for each array. If you have a unique file for each array, name these files HTO_"sample_name".csv

##

# Running the script 
# For HTO only: ./run_CITESeq_Counts.sh HTO "sample1 sample2" "idx1 idx2"
# For ADT only: ./run_CITESeq_Counts.sh ADT "sample1 sample2" "idx1 idx2"
# For both HTO and ADT: ./run_CITESeq_Counts.sh Both "sample1_HTO sample2_HTO" "idx1_HTO idx2_HTO" "sample1_ADT sample2_ADT" "idx1_ADT idx2_ADT"

##

# Function to run CITESeq Counts for HTO
run_CITESeq_Counts_HTO() {
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

# Function to run CITESeq Counts for ADT
run_CITESeq_Counts_ADT() {
    sample=$1
    idx=$2

    celln=$(cat "./starsolo/${sample}/Solo.out/Gene/filtered/barcodes.tsv" | wc -l)

    # Determine the ADT file to use
    if [ -f "ADT_${sample}.csv" ]; then
        ADT_file="ADT_${sample}.csv"
    elif [ -f "ADT.csv" ]; then
        ADT_file="ADT.csv"
    else
        echo "Error: ADT file not found for sample ${sample}."
        return 1
    fi

    nohup CITE-seq-Count -R1 "fastq/ADT_${sample}_S${idx}_R1_001.fastq.gz" \
    -R2 "fastq/ADT_${sample}_S${idx}_R2_001.fastq.gz" \
    -t "$ADT_file" \
    -cbf 1 -cbl 12 \
    -umif 13 -umil 20 \
    -cells $celln \
    -wl "./starsolo/${sample}/Solo.out/Gene/filtered/barcodes.tsv" -o "ADT_${sample}" > "ADT_${sample}.out" &
}

# Main function to process the input type and call the respective function
process_input() {
    type=$1
    shift # Shift the arguments to exclude the first one (type)

    if [ "$type" == "HTO" ]; then
        IFS=' ' read -r -a samples <<< "$1"
        IFS=' ' read -r -a idxs <<< "$2"
        for i in "${!samples[@]}"; do
            run_CITESeq_Counts_HTO "${samples[i]}" "${idxs[i]}"
        done
    elif [ "$type" == "ADT" ]; then
        IFS=' ' read -r -a samples <<< "$1"
        IFS=' ' read -r -a idxs <<< "$2"
        for i in "${!samples[@]}"; do
            run_CITESeq_Counts_ADT "${samples[i]}" "${idxs[i]}"
        done
    elif [ "$type" == "Both" ]; then
        IFS=' ' read -r -a samples_HTO <<< "$1"
        IFS=' ' read -r -a idxs_HTO <<< "$2"
        IFS=' ' read -r -a samples_ADT <<< "$3"
        IFS=' ' read -r -a idxs_ADT <<< "$4"
        for i in "${!samples_HTO[@]}"; do
            run_CITESeq_Counts_HTO "${samples_HTO[i]}" "${idxs_HTO[i]}"
        done
        for i in "${!samples_ADT[@]}"; do
            run_CITESeq_Counts_ADT "${samples_ADT[i]}" "${idxs_ADT[i]}"
        done
    else
        echo "Error: Invalid type. Please choose HTO, ADT, or Both."
        exit 1
    fi
}

# Call the main function with all arguments
process_input "$@"
