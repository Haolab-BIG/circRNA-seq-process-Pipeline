# circRNA-seq Pipeline

This unified circRNA-seq pipeline processes raw paired-end FASTQ files through to differential expression results, using Singularity for reproducibility. It leverages the CIRI2 detector for accurate circRNA identification and the CircTest for statistical analysis.

## Workflow

<img width="2184" height="666" alt="CleanShot 2025-09-21 at 21 55 17@2x" src="https://github.com/user-attachments/assets/943e5529-1099-4ae3-b639-0f0522501c8a" />


## Features

  * **Single Command Execution**: Executes the entire workflow—from FASTQ input and QC, through circRNA quantification, to differential expression analysis—with a single command.
  * **Reproducible**: All software (FastQC, Trim Galore, BWA, Samtools, CIRI2, Python, R/CircTest) is encapsulated within a Singularity container (`circRNA.sif`), ensuring analysis is fully reproducible.
  * **Automated Reporting**: Generates a final, interactive MultiQC report summarizing quality control metrics across all samples and steps for easy assessment.

## Requirements

1.  **Recommended System Configuration**:

      * 8-core CPU
      * 64 GB RAM

2.  **Singularity**: Must be installed on your system. Below are detailed steps for installing on an Ubuntu 22.04 system. For other operating systems, please refer to the [official installation guide](https://www.google.com/search?q=https://docs.sylabs.io/guides/latest/user-guide/installation.html).

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
        	libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go (check for the latest version)
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Navigate to a suitable directory for downloading source code
        cd /tmp

        # Download the Singularity CE source code (check for the latest version)
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure, build, and install Singularity
        ./mconfig
        cd builddir
        make
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version
        ```


3.  **Pipeline Files**:

      * `run_circRNA_pipeline.sh`
      * `circRNA.sif` (The Singularity container)

4.  **Reference Data**: A directory containing all necessary reference files.

## Setup

### 1\. Prepare the Sample Sheet

This is the most critical input file. Create a CSV file named `samplesheet.csv`. The pipeline is designed for **paired-end** sequencing data.

  * `sample_name`: A unique identifier for the sample (e.g., `Control_Rep1`).
  * `condition`: The experimental group for the sample (e.g., `Control`, `Treated`).
  * `fastq1_path`: The **absolute path** to the Read 1 FASTQ file.
  * `fastq2_path`: The **absolute path** to the Read 2 FASTQ file.

**Example `samplesheet.csv`:**

```csv
sample_name,condition,fastq1_path,fastq2_path
Control_Rep1,Control,/path/to/data/C1_R1.fastq.gz,/path/to/data/C1_R2.fastq.gz
Control_Rep2,Control,/path/to/data/C2_R1.fastq.gz,/path/to/data/C2_R2.fastq.gz
Treated_Rep1,Treated,/path/to/data/T1_R1.fastq.gz,/path/to/data/T1_R2.fastq.gz
Treated_Rep2,Treated,/path/to/data/T2_R1.fastq.gz,/path/to/data/T2_R2.fastq.gz
```

### 2\. Prepare the Reference Data

The pipeline requires a reference genome, annotation, and a BWA index.

#### Create Reference Directory

Create a dedicated directory for all reference data:

```bash
mkdir -p reference_data
cd reference_data
```

#### Download Required Files (Human hg38/GRCh38 example)

```bash
# 1. Download Genome FASTA (from GENCODE)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# 2. Download GTF Annotation (from GENCODE)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz
gunzip gencode.v46.primary_assembly.annotation.gtf.gz
mv gencode.v46.primary_assembly.annotation.gtf annotation.gtf
```

#### Build BWA Index

CIRI2 requires reads to be mapped with BWA-MEM. You must build the index using the `bwa` version inside the container.

```bash
# The output prefix should match the FASTA filename for simplicity
singularity exec circRNA.sif bwa index GRCh38.primary_assembly.genome.fa
```

#### Final Reference Structure

Your `reference_data` directory should now contain the following files:

```
reference_data/
├── GRCh38.primary_assembly.genome.fa
├── GRCh38.primary_assembly.genome.fa.amb
├── GRCh38.primary_assembly.genome.fa.ann
├── GRCh38.primary_assembly.genome.fa.bwt
├── GRCh38.primary_assembly.genome.fa.pac
├── GRCh38.primary_assembly.genome.fa.sa
└── annotation.gtf
```

## Running

Execute the pipeline using a single command.

### Command Parameters

  * `-s`: Path to the sample sheet CSV file (required).
  * `-o`: Output directory path where results will be saved (required).
  * `-r`: Reference data directory (required).
  * `-c`: Path to the `circRNA.sif` Singularity container file (required).
  * `-b`: Minimum BSJ count for initial filtering (optional, default: 2).
  * `-t`: Number of threads to use for processing (optional, default: 8).
  * `-h`: Display the help message.

### Example Command

```bash
bash run_circRNA_pipeline.sh \
  -s ./samplesheet.csv \
  -o ./circrna_project_results \
  -r ./reference_data \
  -c ./circRNA.sif \
  -t 16
```

## Output Structure and Interpretation

After the pipeline completes successfully, the output directory (`circrna_project_results/`) will be organized as follows.

```
./circrna_project_results/
├── Control_Rep1/
│   ├── ciri2_output.tsv
│   └── ciri2_output.tsv.log
├── Control_Rep2/
│   └── ... (same structure as Control_Rep1)
├── ...
├── CircTest_results/
│   └── CircTest_significant_results.tsv
├── multiqc_report/
│   └── multiqc_report.html
├── circ_counts_matrix.tsv
└── linear_counts_matrix.tsv
```

-----

### Aggregate Result Files

These files represent the final, combined analysis results from all samples.

  * **`CircTest_results/CircTest_significant_results.tsv`**

      * **Content**: This is the primary result file, containing the list of statistically significant differentially expressed circRNAs from CircTest. Key columns include:
          * `CircID`: The unique identifier of the circRNA.
          * `pval`: The raw p-value from the statistical test.
          * `padj`: The p-value adjusted for multiple testing.
          * `group_1_ratio_mean`: The mean circular-to-linear ratio of group 1.
          * `group_2_ratio_mean`: The mean circular-to-linear ratio of group 2.
      * **Application**: Use this file to identify your top candidate circRNAs. An orthogonal validation method must be used to validate a predicted circRNA; qPCR validation on its own is not sufficient, at least qPCR + RNase R treatment or preferably qPCR + amplicon sequencing should be used.

  * **`circ_counts_matrix.tsv`** and **`linear_counts_matrix.tsv`**

      * **Content**: Two matrices of raw read counts. Rows are circRNAs, and columns are samples. `circ_counts_matrix.tsv` contains the back-spliced junction (BSJ) read counts, while `linear_counts_matrix.tsv` contains reads supporting the linear host gene at the same locus.
      * **Application**: These matrices are the direct inputs for the CircTest analysis. They are useful for custom quality control checks or for creating visualizations like PCA plots or heatmaps.

* **`multiqc_report`**: Open `multiqc_report.html` in a web browser to explore all sections interactively.

  	- **Application**: This is the first file you should check to assess the overall quality of your sequencing data and the alignment process. It helps identify problematic samples (e.g., low alignment rate, high duplication) early on.

    	- **General Statistics**: A combined table summarizing important metrics for each sample:

      <img width="1914" height="1184" alt="CleanShot 2025-09-21 at 21 56 18@2x" src="https://github.com/user-attachments/assets/e8b46539-a5cb-4ef5-9517-f2e6f6eac34f" />


    	- **FastQC**: Quality-control metrics on raw and trimmed reads, including  
      'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores',  
      'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content',  
      'Sequence Length Distribution', 'Sequence Duplication Levels',  
      'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content'.

          - **Sequence Quality Histograms**: The mean quality value across each base position in the read.  

          <img width="1912" height="1322" alt="CleanShot 2025-09-21 at 21 56 33@2x" src="https://github.com/user-attachments/assets/500a5cd6-20c0-4126-a4b0-9fdc3a1d4378" />


          - **Adapter Content**: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.  

          <img width="1926" height="1222" alt="CleanShot 2025-09-21 at 21 56 49@2x" src="https://github.com/user-attachments/assets/85f72bd6-4a52-43f0-99fd-4c3576220595" />

-----

### Per-Sample Files (`Control_Rep1/`, etc.)

These directories contain the direct output from the CIRI2 detection step for each sample.

  * **`ciri2_output.tsv`**

      * **Content**: A detailed, tab-separated file listing every circRNA detected in that sample. The columns are:
          * `circRNA_ID`: Unique identifier (chr:start|end).
          * `chr`, `circRNA_start`, `circRNA_end`: Genomic coordinates.
          * `#junction_reads`: The number of back-spliced junction (BSJ) reads. **This is the circular count.**
          * `SM_MS_SMS`: Information on CIGAR strings of junction reads.
          * `#non_junction_reads`: The number of reads spanning the junction site consistent with linear splicing. **This is the linear count.**
          * `junction_reads_ratio`: A metric of circular vs. linear expression.
          * `circRNA_type`: Annotation of the circRNA (e.g., exon, intron).
          * `gene_id`, `strand`: Host gene information from the GTF file.
          * `junction_reads_ID`: A comma-separated list of all BSJ read IDs.
      * **Application**: This is the raw data file for a single sample. It's useful for deep dives into specific samples or for manual inspection of all detected circRNAs before filtering. The `#junction_reads` and `#non_junction_reads` columns are extracted by the pipeline to build the final count matrices.
   
      <img width="2272" height="102" alt="CleanShot 2025-09-21 at 21 57 54@2x" src="https://github.com/user-attachments/assets/6cb4ca67-952a-4091-8424-fb13471a894a" />


  * **`ciri2_output.tsv.log`**

      * **Content**: A log file from the CIRI2 run, containing runtime information, parameters used, and summary statistics.
      * **Application**: Useful for debugging a failed run or for recording the exact parameters and version of the software used for a specific sample.
      <img width="1484" height="444" alt="CleanShot 2025-09-21 at 21 58 32@2x" src="https://github.com/user-attachments/assets/ca09f4e7-fd48-40c9-9aeb-0f21cd240219" />
