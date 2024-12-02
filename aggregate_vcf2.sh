#!/bin/bash
# Author: Daniela D. Russo <danieladrusso@gmail.com>
#
# Aggregate vcf files from samples run through haplotype caller pipeline
# Output is a single vcf file containing information on samples from the same experiment 
#
set -e
set -o pipefail

# Initialize logging
LOG_FILE="vcf_aggregation_$(date '+%Y%m%d_%H%M%S').log"

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
    echo "  -i vcf_list.txt                    Text file containing paths to VCF files (one per line)"
    echo "  -s SAMPLE_SET_ID                   Sample set ID for output naming"
    echo "  -o MERGE_OPTION                    BCFtools merge option (e.g., 'none', 'snps', 'indels', 'both', 'all')"
    echo "  -m MEMORY_GB                       Memory in GB (default: 4)"
    echo "  -t THREADS                         Number of threads (default: 1)"
    echo "  -d DISK_SPACE_GB                   Disk space in GB (default: 100)"
    echo
    exit 1
}

# Default values
MEMORY_GB=4
THREADS=1
DISK_SPACE_GB=100

# Parse command line arguments
while getopts ":i:s:o:m:t:d:" opt; do
    case $opt in
        i) VCF_LIST="$OPTARG" ;;
        s) SAMPLE_SET_ID="$OPTARG" ;;
        o) MERGE_OPTION="$OPTARG" ;;
        m) MEMORY_GB="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) DISK_SPACE_GB="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check required arguments
if [[ -z "${VCF_LIST}" || -z "${SAMPLE_SET_ID}" || -z "${MERGE_OPTION}" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate input file
if [[ ! -f "${VCF_LIST}" ]]; then
    echo "Error: VCF list file ${VCF_LIST} not found"
    exit 1
fi

# Function to check bcftools installation
check_bcftools() {
    if ! command -v bcftools &> /dev/null; then
        echo "Error: bcftools is not installed or not in PATH"
        exit 1
    fi
    log "Using bcftools version: $(bcftools --version | head -n1)"
}

# Function to validate VCF files
validate_vcfs() {
    local vcf_files=()
    while IFS= read -r vcf; do
        if [[ ! -f "${vcf}" ]]; then
            echo "Error: VCF file not found: ${vcf}"
            exit 1
        fi
        # Check if file is compressed and has index
        if [[ "${vcf}" == *.gz ]]; then
            if [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
                echo "Error: Index file not found for compressed VCF: ${vcf}"
                exit 1
            fi
        fi
        vcf_files+=("${vcf}")
    done < "${VCF_LIST}"
    
    if [[ ${#vcf_files[@]} -eq 0 ]]; then
        echo "Error: No VCF files found in ${VCF_LIST}"
        exit 1
    fi
    
    log "Found ${#vcf_files[@]} VCF files to merge"
}

# Function to merge VCF files
merge_vcfs() {
    local output_vcf="${SAMPLE_SET_ID}.called.vcf"
    
    log "Starting VCF merge"
    log "Output file: ${output_vcf}"
    log "Using merge option: ${MERGE_OPTION}"
    log "Using ${THREADS} threads"
    
    bcftools merge \
        -m "${MERGE_OPTION}" \
        --force-samples \
        --threads "${THREADS}" \
        -o "${output_vcf}" \
        --file-list "${VCF_LIST}"
    
    local exit_code=$?
    if [[ ${exit_code} -eq 0 ]]; then
        log "VCF merge completed successfully"
        log "Output written to: ${output_vcf}"
    else
        log "Error: VCF merge failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
}

# Function to check disk space
check_disk_space() {
    local available_space=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//')
    if [[ ${available_space} -lt ${DISK_SPACE_GB} ]]; then
        echo "Error: Insufficient disk space. Required: ${DISK_SPACE_GB}GB, Available: ${available_space}GB"
        exit 1
    fi
    log "Disk space check passed. Available: ${available_space}GB"
}

# Main execution
main() {
    log "Starting VCF aggregation pipeline"
    
    # Check requirements
    check_bcftools
    check_disk_space
    
    # Validate input files
    validate_vcfs
    
    # Merge VCF files
    merge_vcfs
    
    log "Pipeline completed successfully"
}

# Set up trap for cleanup
trap 'echo "Error occurred. Check logs for details."; exit 1' ERR

# Execute main workflow
main "$@"