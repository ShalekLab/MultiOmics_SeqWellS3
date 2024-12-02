#!/bin/bash
#
# Author: Daniela D. Russo <danieladrusso@gmail.com>
#
# Allows for the passing of one or more sample BAM files to be passed. Samples will undergo each step of GATK preprocessing for variant discovery as outlined here: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
# Samples are processed individually, in parallel
#
# Enable error tracing and exit on error
set -e
set -o pipefail

# Function for logging
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] $1" | tee -a "${LOG_FILE}"
}

# Function to show usage
usage() {
    echo "Usage: $0 [options]"
    echo
    echo "Required arguments:"
    echo "  -s sample1 [sample2 sample3 ...]   Specify one or more sample names"
    echo "  -r REF_NAME                        Reference name"
    echo "  -f REF_FASTA                       Reference FASTA file"
    echo "  -i REF_FASTA_INDEX                 Reference FASTA index file"
    echo "  -d REF_DICT                        Reference dictionary file"
    echo "  -v DBSNP_VCF                       dbSNP VCF file"
    echo "  -x DBSNP_VCF_INDEX                 dbSNP VCF index file"
    echo "  -b UNMAPPED_BAM_LIST               List of unmapped BAM files"
    echo "  -u UNMAPPED_BAM_SUFFIX             Suffix for unmapped BAM files"
    echo "  -k KNOWN_INDELS_SITES_VCF          VCF file with known indels sites"
    echo
    echo "Optional arguments:"
    echo "  -h                                 Display this help message"
    echo
    exit 1
}

# Parse command line arguments
declare -a SAMPLE_NAMES
while getopts ":s:r:f:i:d:v:x:b:u:k:h" opt; do
    case $opt in
        s) 
            # Handle multiple sample names
            for sample in ${OPTARG}; do
                SAMPLE_NAMES+=("$sample")
            done
            ;;
        r) REF_NAME="$OPTARG" ;;
        f) REF_FASTA="$OPTARG" ;;
        i) REF_FASTA_INDEX="$OPTARG" ;;
        d) REF_DICT="$OPTARG" ;;
        v) DBSNP_VCF="$OPTARG" ;;
        x) DBSNP_VCF_INDEX="$OPTARG" ;;
        b) UNMAPPED_BAM_LIST="$OPTARG" ;;
        u) UNMAPPED_BAM_SUFFIX="$OPTARG" ;;
        k) KNOWN_INDELS_SITES_VCF="$OPTARG" ;;
        h) usage ;;  # Display usage if -h is passed
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check if at least one sample name is provided
if [ ${#SAMPLE_NAMES[@]} -eq 0 ]; then
    echo "Error: At least one sample name must be provided with -s flag"
    usage
fi

# Check other required arguments
required_args=("REF_NAME" "REF_FASTA" "REF_FASTA_INDEX" "REF_DICT" 
              "DBSNP_VCF" "DBSNP_VCF_INDEX" "UNMAPPED_BAM_LIST" 
              "UNMAPPED_BAM_SUFFIX" "KNOWN_INDELS_SITES_VCF")

for arg in "${required_args[@]}"; do
    if [[ -z "${!arg}" ]]; then
        echo "Error: Required argument $arg is not set"
        usage
    fi
done

# Configuration variables
GATK="/gatk/gatk"
BWA_PATH="/usr/gitc/"
COMPRESSION_LEVEL=5
LOG_FILE="pipeline_$(date '+%Y%m%d_%H%M%S').log"

# Function to process a single sample
process_sample() {
    local sample_name=$1
    local sample_log_file="${sample_name}_pipeline.log"
    local checkpoint_file=".${sample_name}_checkpoint"
    local base_filename="${sample_name}.${REF_NAME}"
    
    log "Starting processing for sample: ${sample_name}"
    
    # Get current checkpoint or start new
    local current_checkpoint="start"
    if [[ -f "${checkpoint_file}" ]]; then
        current_checkpoint=$(cat "${checkpoint_file}")
        log "Resuming from checkpoint: ${current_checkpoint}"
    fi
    
    # Process unmapped BAMs
    if [[ "$current_checkpoint" == "start" ]]; then
        log "Starting alignment for ${sample_name}"
        while read -r unmapped_bam; do
            if [[ "$unmapped_bam" == *"${sample_name}"* ]]; then
                bam_basename=$(basename ${unmapped_bam} ${UNMAPPED_BAM_SUFFIX})
                align_reads ${unmapped_bam} ${bam_basename}
                echo "aligned_${bam_basename}" > "${checkpoint_file}"
            fi
        done < "${UNMAPPED_BAM_LIST}"
    fi
    
    # Mark duplicates
    if [[ "$current_checkpoint" == "aligned_"* || "$current_checkpoint" == "start" ]]; then
        log "Marking duplicates for ${sample_name}"
        mark_duplicates "${base_filename}.aligned.unsorted.duplicates_marked" \
                       "${base_filename}.duplicate_metrics"
        echo "marked_duplicates" > "${checkpoint_file}"
    fi
    
    # Sort and fix tags
    if [[ "$current_checkpoint" == "marked_duplicates" || "$current_checkpoint" == "aligned_"* || "$current_checkpoint" == "start" ]]; then
        log "Sorting BAM and fixing tags for ${sample_name}"
        ${GATK} --java-options "-Dsamjdk.compression_level=${COMPRESSION_LEVEL} -Xms9G" \
            SortSam \
            --INPUT "${base_filename}.aligned.unsorted.duplicates_marked.bam" \
            --OUTPUT /dev/stdout \
            --SORT_ORDER "coordinate" \
            --CREATE_INDEX false \
            --CREATE_MD5_FILE false \
            | \
            ${GATK} --java-options "-Dsamjdk.compression_level=${COMPRESSION_LEVEL} -Xms1G" \
            SetNmMdAndUqTags \
            --INPUT /dev/stdin \
            --OUTPUT "${base_filename}.aligned.duplicate_marked.sorted.bam" \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --REFERENCE_SEQUENCE ${REF_FASTA}
        echo "sorted_fixed_tags" > "${checkpoint_file}"
    fi
    
    # Base recalibration
    if [[ "$current_checkpoint" == "sorted_fixed_tags" || "$current_checkpoint" == "marked_duplicates" || "$current_checkpoint" == "aligned_"* || "$current_checkpoint" == "start" ]]; then
        log "Performing base recalibration for ${sample_name}"
        ${GATK} --java-options "-Xms4G" \
            BaseRecalibrator \
            -R ${REF_FASTA} \
            -I "${base_filename}.aligned.duplicate_marked.sorted.bam" \
            --use-original-qualities \
            -O "${base_filename}.recal_data.csv" \
            --known-sites ${DBSNP_VCF} \
            --known-sites ${KNOWN_INDELS_SITES_VCF}
        echo "recalibrated" > "${checkpoint_file}"
    fi
    
    # Apply BQSR
    if [[ "$current_checkpoint" == "recalibrated" || "$current_checkpoint" == "sorted_fixed_tags" || "$current_checkpoint" == "marked_duplicates" || "$current_checkpoint" == "aligned_"* || "$current_checkpoint" == "start" ]]; then
        log "Applying BQSR for ${sample_name}"
        ${GATK} --java-options "-Xms3G" \
            ApplyBQSR \
            -R ${REF_FASTA} \
            -I "${base_filename}.aligned.duplicate_marked.sorted.bam" \
            -O "${base_filename}.aligned.duplicate_marked.recalibrated.bam" \
            -bqsr "${base_filename}.recal_data.csv" \
            --static-quantized-quals 10 \
            --static-quantized-quals 20 \
            --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --create-output-bam-md5 \
            --use-original-qualities
        echo "completed" > "${checkpoint_file}"
    fi
    
    # Clean up intermediate files if processing completed successfully
    if [[ "$current_checkpoint" == "completed" ]]; then
        log "Cleaning up intermediate files for ${sample_name}"
        rm -f "${base_filename}.aligned.unsorted.duplicates_marked.bam"
        rm -f "${base_filename}.aligned.duplicate_marked.sorted.bam"
        rm -f "${base_filename}.aligned.duplicate_marked.sorted.bai"
        rm -f "${base_filename}.aligned.duplicate_marked.sorted.bam.md5"
        rm -f "${checkpoint_file}"
    fi
    
    log "Completed processing for sample: ${sample_name}"
}

export -f process_sample
export -f log
export -f align_reads
export -f mark_duplicates
export -f base_recalibration
export -f apply_bqsr

# Main execution
main() {
    log "Starting pipeline execution for ${#SAMPLE_NAMES[@]} samples"
    
    # Process samples in parallel using GNU parallel
    printf '%s\n' "${SAMPLE_NAMES[@]}" | parallel --jobs 80% process_sample {}
    
    log "Pipeline completed for all samples"
}

# Set up trap for cleanup
trap cleanup EXIT
trap 'log "Error occurred. Check logs for details."; exit 1' ERR

# Execute main workflow
main "$@"
