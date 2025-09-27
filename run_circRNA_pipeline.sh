#!/bin/bash

# --- Stricter Shell Safety ---
set -e
set -o pipefail
IFS=$'\n\t'

# ==============================================================================
#
#           Single-File circRNA Pipeline with CIRI2 and CircTest
#
# Description:
#   This self-contained script executes a complete circRNA-seq workflow.
#   It uses CIRI2 for detection and the CircTest R package for differential
#   analysis. The necessary Python and R helper scripts are automatically

#   generated at runtime.
#
# ==============================================================================


# --- 1. Default Parameters & Help Function ---
THREADS=8
BSJ_FILTER=2 # Default filter: require at least 2 BSJ reads for initial extraction

HELP_MSG="Usage: $0 -s <sample_sheet.csv> -o <out_dir> -r <ref_dir> -c <circRNA.sif> [OPTIONS]

Required:
  -s  Sample sheet (CSV). Format: sample_name,condition,fastq1_path,fastq2_path
  -o  Main output directory.
  -r  Reference data directory (must contain BWA indices and annotation).
  -c  Path to the Singularity container (circRNA.sif).

Optional:
  -b  Minimum BSJ count for filtering circRNAs before aggregation (Default: 2).
  -t  Number of threads to use (Default: 8).
  -h  Display this help message.
"

# --- 2. Parse Command-line Arguments ---
while getopts "s:o:r:c:b:t:h" opt; do
  case ${opt} in
    s ) SAMPLE_SHEET=$(realpath "${OPTARG}") ;;
    o ) OUT_DIR=$(realpath "${OPTARG}") ;;
    r ) REF_DIR=$(realpath "${OPTARG}") ;;
    c ) SIF_PATH=$(realpath "${OPTARG}") ;;
    b ) BSJ_FILTER=${OPTARG} ;;
    t ) THREADS=${OPTARG} ;;
    h ) echo "${HELP_MSG}"; exit 0 ;;
    \? ) echo "Invalid option: -${OPTARG}" >&2; echo "${HELP_MSG}"; exit 1 ;;
  esac
done

# Check for mandatory arguments
if [ -z "${SAMPLE_SHEET}" ] || [ -z "${OUT_DIR}" ] || [ -z "${REF_DIR}" ] || [ -z "${SIF_PATH}" ]; then
    echo "Error: Missing mandatory arguments." >&2; echo "${HELP_MSG}"; exit 1
fi
if [ ! -f "${SIF_PATH}" ]; then
    echo "Error: Singularity container not found at: ${SIF_PATH}" >&2; exit 1
fi

# --- 3. Setup Environment and Generate Helper Scripts ---
mkdir -p "${OUT_DIR}"
SINGULARITY_BASE_CMD="singularity exec --cleanenv -B ${OUT_DIR}:/output -B ${REF_DIR}:/reference"

echo "==================================================="
echo "====== circRNA Pipeline with CircTest Started ======"
echo "==================================================="
echo "Sample Sheet: ${SAMPLE_SHEET}"
echo "Output Directory: ${OUT_DIR}"
echo "Reference Directory: ${REF_DIR}"
echo "BSJ Filter Threshold: ${BSJ_FILTER}"
echo "Threads: ${THREADS}"
echo "Singularity Container: ${SIF_PATH}"
echo "==================================================="

# --- 3a. Generate Python Aggregation Script ---
AGGREGATION_SCRIPT_HOST_PATH="${OUT_DIR}/aggregate_counts.py"
AGGREGATION_SCRIPT_CONT_PATH="/output/aggregate_counts.py"
echo "Generating Python helper script at ${AGGREGATION_SCRIPT_HOST_PATH}"
cat <<'EOF' > "${AGGREGATION_SCRIPT_HOST_PATH}"
import pandas as pd
import argparse
import os

def aggregate_counts(sample_sheet_path, input_dir, suffix, output_path):
    """
    Aggregates individual count files into a single matrix.
    
    Args:
        sample_sheet_path (str): Path to the sample sheet (CSV: sample_name,condition).
        input_dir (str): Directory containing sample subdirectories with count files.
        suffix (str): Suffix of the count files to aggregate (e.g., '.circ.counts').
        output_path (str): Path to save the final aggregated matrix.
    """
    try:
        sample_sheet = pd.read_csv(sample_sheet_path, header=0)
        sample_names = sample_sheet.iloc[:, 0].tolist()
    except Exception as e:
        print(f"Error reading sample sheet {sample_sheet_path}: {e}")
        return

    # Initialize the master dataframe with the first sample
    first_sample_name = sample_names[0]
    first_file_path = os.path.join(input_dir, first_sample_name, f"{first_sample_name}{suffix}")
    
    if not os.path.exists(first_file_path):
        print(f"Error: Count file for the first sample not found at {first_file_path}")
        return
        
    master_df = pd.read_csv(first_file_path, sep='\t', header=None, names=['CircID', first_sample_name])
    
    # Loop through the rest of the samples and merge
    for sample_name in sample_names[1:]:
        file_path = os.path.join(input_dir, sample_name, f"{sample_name}{suffix}")
        if os.path.exists(file_path):
            sample_df = pd.read_csv(file_path, sep='\t', header=None, names=['CircID', sample_name])
            master_df = pd.merge(master_df, sample_df, on='CircID', how='outer')
        else:
            print(f"Warning: Count file not found for sample {sample_name}. Skipping.")
            master_df[sample_name] = 0 # Add column with zeros if file is missing
            
    # Fill NaN values with 0 (for circRNAs not found in some samples) and set index
    master_df = master_df.fillna(0)
    
    # Convert count columns to integer
    for col in master_df.columns[1:]:
        master_df[col] = master_df[col].astype(int)
        
    master_df.to_csv(output_path, sep='\t', index=False)
    print(f"Successfully created count matrix at: {output_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Aggregate circRNA count files into a matrix.")
    parser.add_argument('--sample_sheet', required=True, help='Path to the sample sheet (CSV format: sample_name,condition).')
    parser.add_argument('--input_dir', required=True, help='Directory containing sample subdirectories.')
    parser.add_argument('--suffix', required=True, help='Suffix of the count files (e.g., .circ.counts).')
    parser.add_argument('--output', required=True, help='Path for the output matrix file.')
    
    args = parser.parse_args()
    
    aggregate_counts(args.sample_sheet, args.input_dir, args.suffix, args.output)
EOF

# --- 3b. Generate R CircTest Script ---
CIRCTEST_SCRIPT_HOST_PATH="${OUT_DIR}/run_CircTest.R"
CIRCTEST_SCRIPT_CONT_PATH="/output/run_CircTest.R"
echo "Generating R helper script at ${CIRCTEST_SCRIPT_HOST_PATH}"
cat <<'EOF' > "${CIRCTEST_SCRIPT_HOST_PATH}"
#!/usr/bin/env Rscript

# Load necessary libraries. Ensure these are installed in your container.
# install.packages(c("optparse", "devtools", "ggplot2"))
# devtools::install_github('dieterich-lab/CircTest')
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(CircTest))
suppressPackageStartupMessages(library(ggplot2))

# --- Argument Parsing ---
option_list <- list(
  make_option(c("-c", "--circ_matrix"), type="character", help="Path to the circular counts matrix"),
  make_option(c("-l", "--linear_matrix"), type="character", help="Path to the linear counts matrix"),
  make_option(c("-s", "--sample_sheet"), type="character", help="Path to the sample sheet (CSV: sample,condition)"),
  make_option(c("-o", "--out_dir"), type="character", help="Path to the output directory")
)
args <- parse_args(OptionParser(option_list=option_list))

# --- Main Analysis ---

# 1. Create output directory
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Read input files
message("Reading input matrices and sample sheet...")
circ_counts <- read.delim(args$circ_matrix, header = TRUE, as.is = TRUE, sep = "\t")
linear_counts <- read.delim(args$linear_matrix, header = TRUE, as.is = TRUE, sep = "\t")
sample_info <- read.csv(args$sample_sheet, header = TRUE)

# 3. Prepare data for CircTest
# Ensure row names are CircIDs for easy filtering
rownames(circ_counts) <- circ_counts[, 1]
rownames(linear_counts) <- linear_counts[, 1]

# Ensure sample order matches between matrices and sample sheet
circ_counts <- circ_counts[, c("CircID", sample_info[,1])]
linear_counts <- linear_counts[, c("CircID", sample_info[,1])]

# Get group vector and number of replicates per condition
conditions <- unique(sample_info$condition)
group_vector <- as.numeric(factor(sample_info$condition, levels = conditions))
n_replicates <- as.numeric(table(group_vector)[1]) # Assumes equal replicates

message(paste("Detected conditions:", paste(conditions, collapse=", ")))
message(paste("Number of replicates per condition:", n_replicates))

# 4. Filter data
message("Filtering count tables...")
circ_filtered <- Circ.filter(
  circ = circ_counts,
  linear = linear_counts,
  Nreplicates = n_replicates,
  filter.sample = n_replicates,
  filter.count = 0,
  percentage = 0.01,
  circle_description = 1
)

if (nrow(circ_filtered) == 0) {
  stop("No circRNAs passed the filtering criteria. Analysis cannot continue.")
}

# Filter the linear table to match the filtered circular table
linear_filtered <- linear_counts[rownames(circ_filtered), ]
message(paste(nrow(circ_filtered), "circRNAs passed the filter."))


# 5. Run Circ.test
message("Running Circ.test for differential expression...")
test_results <- Circ.test(
  circ_filtered,
  linear_filtered,
  group = group_vector,
  alpha=1, plotsig=F,
  circle_description = 1
)
final <- data.frame(test_results$summary_table[,1],test_results$p.val,test_results$p.adj,test_results$ratios)
colnames(final) <- c("CircID","p.val","p.adj","group_1_ratio_mean","group_2_ratio_mean")
# 6. Save results
message("Saving results...")
if (!is.null(test_results$summary_table) && nrow(test_results$summary_table) > 0) {
  write.csv(final,file.path(args$out_dir, "CircTest_results.csv"), row.names = FALSE)
} else {
  message("No significant differentially expressed circRNAs found.")
  write("No significant differentially expressed circRNAs found.", file.path(args$out_dir, "CircTest_no_significant_results.txt"))
}

message("CircTest analysis complete!")
EOF

# --- 4. Loop Through and Process Each Sample ---
# Skip header in sample sheet
tail -n +2 "${SAMPLE_SHEET}" | while IFS=, read -r SAMPLE_NAME CONDITION FQ1 FQ2; do
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    CONDITION=$(echo "${CONDITION}" | tr -d '[:space:]')
    FQ1=$(echo "${FQ1}" | tr -d '[:space:]')
    FQ2=$(echo "${FQ2}" | tr -d '[:space:]')

    echo -e "\n\n======== Starting to process sample: [${SAMPLE_NAME}] ========"
    
    FQ_DIR=$(dirname "${FQ1}")
    SAMPLE_OUT_DIR_HOST="${OUT_DIR}/${SAMPLE_NAME}"
    SAMPLE_OUT_DIR_CONT="/output/${SAMPLE_NAME}" # Path inside the container
    mkdir -p "${SAMPLE_OUT_DIR_HOST}"
    SINGULARITY_CMD="singularity exec --cleanenv -B ${OUT_DIR}:/output -B ${REF_DIR}:/reference -B ${FQ_DIR}:/reads ${SIF_PATH}"

    # --- 4.1. QC and Trimming ---
    echo "======== [${SAMPLE_NAME}] Step 1: Trimming and Raw QC ========"
    eval "$SINGULARITY_CMD fastqc -t ${THREADS} /reads/$(basename ${FQ1}) /reads/$(basename ${FQ2}) -o ${SAMPLE_OUT_DIR_CONT}"
    eval "$SINGULARITY_CMD trim_galore --paired --cores ${THREADS} -o ${SAMPLE_OUT_DIR_CONT} /reads/$(basename ${FQ1}) /reads/$(basename ${FQ2})"

    FQ1_BASENAME=$(basename "${FQ1}")
    FQ2_BASENAME=$(basename "${FQ2}")
    FQ1_TRIMMED_BASENAME=$(echo "${FQ1_BASENAME}" | sed -E 's/(\.fq|\.fastq)(\.gz)?$/_val_1.fq.gz/')
    FQ2_TRIMMED_BASENAME=$(echo "${FQ2_BASENAME}" | sed -E 's/(\.fq|\.fastq)(\.gz)?$/_val_2.fq.gz/')
    TRIMMED_FQ1="${SAMPLE_OUT_DIR_CONT}/${FQ1_TRIMMED_BASENAME}"
    TRIMMED_FQ2="${SAMPLE_OUT_DIR_CONT}/${FQ2_TRIMMED_BASENAME}"
    
    echo "======== [${SAMPLE_NAME}] Step 2: Post-trimming FastQ QC ========"
    eval "$SINGULARITY_CMD fastqc -t ${THREADS} ${TRIMMED_FQ1} ${TRIMMED_FQ2} -o ${SAMPLE_OUT_DIR_CONT}"

    # --- 4.2. Alignment for CIRI2 ---
    echo "======== [${SAMPLE_NAME}] Step 3a: Alignment for CIRI2 (BWA-MEM) ========"
    SAM_PREFIX="${SAMPLE_OUT_DIR_CONT}/${SAMPLE_NAME}"
    eval "$SINGULARITY_CMD bash -c 'bwa mem -T ${THREADS} /reference/GRCh38.primary_assembly.genome.fa ${TRIMMED_FQ1} ${TRIMMED_FQ2} > ${SAM_PREFIX}.sam'"
    
    # --- 4.3. Run circRNA Detector ---
    echo "======== [${SAMPLE_NAME}] Step 3b: circRNA Detection with CIRI2 ========"
    CIRI2_OUT_CONT="${SAMPLE_OUT_DIR_CONT}/ciri2_output.tsv"
    eval "$SINGULARITY_CMD perl /mnt/software/CIRI2.pl -I ${SAM_PREFIX}.sam -O ${CIRI2_OUT_CONT} -F /reference/GRCh38.primary_assembly.genome.fa -A /reference/annotation.gtf -T ${THREADS}"
    
    # --- 4.4. Extract Counts for CircTest ---
    echo "======== [${SAMPLE_NAME}] Step 4: Extracting counts for CircTest ========"
    CIRI2_OUT_HOST="${SAMPLE_OUT_DIR_HOST}/ciri2_output.tsv"
    CIRI2_CIRC_COUNTS="${SAMPLE_OUT_DIR_HOST}/${SAMPLE_NAME}.circ.counts"
    CIRI2_LINEAR_COUNTS="${SAMPLE_OUT_DIR_HOST}/${SAMPLE_NAME}.linear.counts"

    # Extract circRNA_ID (col 1), back-spliced reads (col 5), and linear reads (col 7)
    awk -F'\t' -v OFS='\t' 'NR > 1 && $5 >= '${BSJ_FILTER}' {print $1, $5}' "${CIRI2_OUT_HOST}" > "${CIRI2_CIRC_COUNTS}"
    awk -F'\t' -v OFS='\t' 'NR > 1 && $5 >= '${BSJ_FILTER}' {print $1, $7}' "${CIRI2_OUT_HOST}" > "${CIRI2_LINEAR_COUNTS}"

done

echo -e "\n\n======= All samples processed. Aggregating results and running CircTest... ======="

# --- 5. Aggregate Counts and Run CircTest ---
R_SAMPLE_SHEET_HOST="${OUT_DIR}/samples_for_r.csv"
R_SAMPLE_SHEET_CONT="/output/samples_for_r.csv"

# Create a clean sample sheet for R (name, condition)
cut -d, -f1,2 "${SAMPLE_SHEET}" > "${R_SAMPLE_SHEET_HOST}"

# Aggregate circular (back-spliced) counts
echo "======== Step 5a: Aggregating circular read counts ========"
eval "$SINGULARITY_BASE_CMD ${SIF_PATH} python3 ${AGGREGATION_SCRIPT_CONT_PATH} \
    --sample_sheet ${R_SAMPLE_SHEET_CONT} \
    --input_dir /output \
    --suffix .circ.counts \
    --output /output/circ_counts_matrix.tsv"

# Aggregate linear (host gene) counts
echo "======== Step 5b: Aggregating linear read counts ========"
eval "$SINGULARITY_BASE_CMD ${SIF_PATH} python3 ${AGGREGATION_SCRIPT_CONT_PATH} \
    --sample_sheet ${R_SAMPLE_SHEET_CONT} \
    --input_dir /output \
    --suffix .linear.counts \
    --output /output/linear_counts_matrix.tsv"

# Run R script for CircTest
echo "======== Step 5c: Running CircTest for differential analysis ========"
eval "$SINGULARITY_BASE_CMD ${SIF_PATH} Rscript ${CIRCTEST_SCRIPT_CONT_PATH} \
    --circ_matrix /output/circ_counts_matrix.tsv \
    --linear_matrix /output/linear_counts_matrix.tsv \
    --sample_sheet ${R_SAMPLE_SHEET_CONT} \
    --out_dir /output/CircTest_results"

# --- 6. Generate MultiQC Report ---
echo "======== Step 6: Generating Final MultiQC Report ========"
MULTIQC_DIR_HOST="${OUT_DIR}/multiqc_report"
mkdir -p "${MULTIQC_DIR_HOST}"
eval "$SINGULARITY_BASE_CMD ${SIF_PATH} multiqc /output -o /output/multiqc_report --force"
while IFS=, read -r SAMPLE_NAME CONDITION FQ1 FQ2; do
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    find "${OUT_DIR}/${SAMPLE_NAME}" -type f \
        ! -name "ciri2*" \
        -delete
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')
rm -f "${OUT_DIR}/aggregate_counts.py"
rm -f "${OUT_DIR}/run_CircTest.R"
rm -f "${OUT_DIR}/samples_for_r.csv"
echo -e "\n\n==================================================="
echo "====== circRNA Pipeline Completed Successfully! ======"
echo "==================================================="
echo "Final results are in: ${OUT_DIR}"
echo "CircTest results are in: ${OUT_DIR}/CircTest_results"
echo "View QC report at: ${MULTIQC_DIR_HOST}/multiqc_report.html"
