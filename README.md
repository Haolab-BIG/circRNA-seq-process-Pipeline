# circRNA-seq-Processing-Pipeline
This pipeline provides a fully containerized Singularity environment that bundles all required tools and dependencies. With a single command, the entire circRNA-seq workflow—from raw FASTQ input through trimming, quality control, genome alignment (BWA-MEM), linear gene counting (featureCounts), circRNA detection (CIRI3), and differential expression analysis (edgeR & GLM)—can be executed reproducibly on any compatible system.

# Part I Workflow
<img width="798" height="248" alt="CleanShot 2025-11-06 at 18 47 03@2x" src="https://github.com/user-attachments/assets/615f63b9-3408-48ba-b4bd-6d6026580566" />


# Part II Requirements
1.  **Recommended Specs**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

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
        # Download and install Go
        wget [https://go.dev/dl/go1.21.3.linux-amd64.tar.gz](https://go.dev/dl/go1.21.3.linux-amd64.tar.gz)
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget [https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz](https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz)

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

      ```bash
      pip install snakemake
      ```

4.  **Reference Data**: A directory containing the BWA index (Detailed steps for the human hg38 genome are below. For other reference genomes, please download the corresponding files and replace as needed).
      ```bash
      mkdir reference_data
      cd reference_data
      
      # Download Genome FASTA
      wget [https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz)
      
      # Download Genome GTF
      wget [https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz)
      
      # Unzip the files
      gunzip GRCh38.primary_assembly.genome.fa.gz
      gunzip gencode.v46.primary_assembly.annotation.gtf.gz
      
      # (Optional) Rename for clarity
      mv GRCh38.primary_assembly.genome.fa GRCh38.primary_assembly.genome.fa
      mv gencode.v46.primary_assembly.annotation.gtf annotation.gtf
      
      # Build BWA index
      # This command will use the Singularity container to ensure 'bwa' is available
      # Replace 'circRNA_modified.sif' with the path to your container
      singularity exec --cleanenv /path/to/circRNA_modified.sif bwa index GRCh38.primary_assembly.genome.fa
      ```

5.  **Required Files**:

      ```bash
      project_directory/
      ├── snakemake/
            ├── config.yaml
            ├── circRNA-seq.smk
            └── de_ciri3_all.R
      ├── Containers/
            └── circRNA_modified.sif
      ├── reference_data/
            ├── illumina_adapter.fa
            ├── GRCh38.primary_assembly.genome.fa
            ├── GRCh38.primary_assembly.genome.fa.amb
            ├── GRCh38.primary_assembly.genome.fa.ann
            ├── GRCh38.primary_assembly.genome.fa.bwt
            ├── GRCh38.primary_assembly.genome.fa.pac
            ├── GRCh38.primary_assembly.genome.fa.sa
            └── annotation.gtf
      └── fq/
            ├── Control_Rep1_R1.fastq.gz
            ├── Control_Rep1_R2.fastq.gz
            ├── Treated_Rep1_R1.fastq.gz
            └── ... (etc)
      ```
      
      - **circRNA-seq.smk** — The main Snakemake workflow script.
      - **config.yaml** — Configuration file containing paths, parameters, and sample information.
      - **de_ciri3_all.R** — R script for all three types of circRNA differential expression analysis.
      - **circRNA_modified.sif** — Singularity container image with all required software (FastQC, TrimGalore, BWA, samtools, featureCounts, CIRI3, R-libraries) and dependencies pre-installed.
      - **illumina_adapter.fa** — FASTA file containing Illumina adapter sequences; replace with your own if needed.
      - **GRCh38.primary_assembly.genome.fa** — Reference genome FASTA file.
      - **annotation.gtf** — Genome annotation GTF file.
      - **BWA Index Files** (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) — Genome index generated by `bwa index`.

# Part III Running

   * **Example code**

      * **Step 1: Edit `config.yaml`**

        ```yaml
        # --- 1. Container ---
        # Absolute path to Singularity image (.sif file)
        # Must contain: fastqc, trim_galore, bwa, samtools, subread (featureCounts), java, multiqc
        singularity_image: "/mnt/guanli/zh.4.cirRNAseq.singularity/snakemake/circRNA_modified.sif"

        # --- 2. Output directory ---
        output_dir: "/mnt/guanli/zh.4.cirRNAseq.singularity/snakemake/circleRNA_results"

        # --- 3. Reference ---
        ref:
          fasta: "/mnt/guanli/zh.1.rnaseq.singularity/reference_data/GRCh38.primary_assembly.genome.fa"
          gtf: "/mnt/guanli/zh.1.rnaseq.singularity/reference_data/annotation.gtf"
          bwa_index: "/mnt/guanli/zh.1.rnaseq.singularity/reference_data/GRCh38.primary_assembly.genome.fa"

        # --- 4. Adapter and trimming ---
        trim:
          adapter_fa: "/mnt/guanli/zh.1.rnaseq.singularity/reference_data/illumina_adapter.fa"
          min_length: 20 # Minimum read length to keep

        # --- 5. BWA-MEM parameters ---
        bwa:
          min_score: 19 # BWA-MEM -T option

        # --- 6. featureCounts parameters ---
        featurecounts:
          threads: 8
          opts: "-p -s 0 -g gene_id" # -p (paired-end), -s 0 (unstranded), -g (group by gene_id)

        # --- 7. CIRI3 parameters ---
        ciri3:
          threads: 8 # Threads for multi-sample run
          threads_single: 4 # Threads for single-sample runs
          max_span: 200000
          min_span: 140
          mapq_uni: 10
          strigency: 2 # 2: >2 PCC signals; 1: >2 junction reads; 0: all

        # --- 8. Samples ---
        # rnase: 0 = untreated, 1 = RNase R treated, 2 = unknown
        samples:
          Control_Rep1:
            R1: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Control_Rep1_R1.fastq.gz"
            R2: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Control_Rep1_R2.fastq.gz"
            condition: "Control"
            rnase: 0
          Control_Rep2:
            R1: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Control_Rep2_R1.fastq.gz"
            R2: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Control_Rep2_R2.fastq.gz"
            condition: "Control"
            rnase: 0
          Treated_Rep1:
            R1: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Treated_Rep1_R1.fastq.gz"
            R2: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Treated_Rep1_R2.fastq.gz"
            condition: "Treated"
            rnase: 0
          Treated_Rep2:
            R1: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Treated_Rep2_R1.fastq.gz"
            R2: "/mnt/guanli/zh.4.cirRNAseq.singularity/fq/Treated_Rep2_R2.fastq.gz"
            condition: "Treated"
            rnase: 0

        # --- 9. Path to unified DE script ---
        de_script: "/mnt/guanli/zh.4.cirRNAseq.singularity/snakemake/de_ciri3_all.R"
        ```

      * **Step 2: run snakemake**
        (Ensure you run this from the directory containing `circRNA-seq.smk` and `config.yaml`)

        ```bash
        snakemake -s circRNA-seq.smk --use-singularity --cores 8 \
          --singularity-args "--bind /mnt/guanli"
        ```
        Then delete the intermediate files and folders:

        ```bash
        rm -r circleRNA_results/qc/ circleRNA_results/trimmed/ circleRNA_results/alignment/
        ```

        *Note*: The `--bind` parameter is crucial. It must include all parent directories referenced in your `config.yaml` (e.g., paths for `output_dir`, `ref`, `samples` R1/R2, `de_script`, and the `singularity_image` itself).

   * **Command Parameters**

      **edit `config.yaml`**
      - `singularity_image`: Absolute path to the Singularity image (`.sif` file) (required).
      - `output_dir`: Path to the directory where all output will be stored (required).
      - `ref`:
          - `fasta`: Path to the reference genome FASTA file (required).
          - `gtf`: Path to the genome annotation GTF file (required).
          - `bwa_index`: Prefix for the BWA index (usually same as FASTA path) (required).
      - `trim`:
          - `adapter_fa`: Path to the adapter sequences FASTA file (required).
          - `min_length`: Minimum read length to keep after trimming (required).
      - `bwa`:
          - `min_score`: BWA-MEM `-T` option, alignments scoring lower will be discarded (required).
      - `featurecounts`:
          - `threads`: Number of threads for featureCounts (required).
          - `opts`: Extra arguments for featureCounts (e.g., `-p` for paired-end, `-s 0` for unstranded) (required).
      - `ciri3`:
          - `threads`: Threads for CIRI3 multi-sample mode (required).
          - `threads_single`: Threads for CIRI3 single-sample mode (required).
          - `max_span`, `min_span`, `mapq_uni`, `strigency`: CIRI3 detection parameters (required).
      - `samples`: Sample definition block (required).
          - `Sample_Name`: A unique identifier for each sample.
          - `R1`: Path to the R1 FASTQ file (required).
          - `R2`: Path to the R2 FASTQ file (leave empty if not paired-end).
          - `condition`: Experimental group for the sample (e.g., "Control", "Treated") (required).
          - `rnase`: RNase R treatment status (0=untreated, 1=treated, 2=unknown) (required).
      - `de_script`: Absolute path to the `de_ciri3_all.R` script (required).

      **run snakemake**
      - `--use-singularity`: Enables execution of rules within a Singularity container to ensure a fully reproducible environment.
      - `--singularity-args`: Allows passing additional arguments to the Singularity runtime (e.g., `--bind`, `--nv`, or custom options).
      - `--cores`: Specifies the maximum number of CPU cores (threads) that Snakemake can use in parallel when executing workflow rules.
      - `--bind`: Specifies the directories to be mounted within the Singularity container. Include all required paths such as raw data, scripts, container images, and references.

# Part IV Output

   * **Output Structure**
      ```bash
      circleRNA_results/
      ├── ciri3_de/
            ├── circ_Gene.txt
            ├── DE_BSJ_results.txt
            ├── DE_Ratio_results.txt
            ├── DE_Relative_results.txt
            ├── infor_bsj_relative.tsv
            └── infor_ratio.tsv
      ├── ciri3_multi/
            ├── results.txt
            ├── results.txt.BSJ_Matrix
            ├── results.txt.FSJ_Matrix
            └── samples.tsv
      ├── ciri3_single/
            ├── Control_Rep1.ciri3.txt
            ├── ... (files for each sample)
      ├── featureCounts/
            ├── counts.txt
            ├── counts.txt.summary
            └── Gene_Expression.txt
      ├── multiqc/
            ├── multiqc_data/
            └── multiqc_report.html
      ```
      
   * **Output Interpretation**
      - **`featureCounts/Gene_Expression.txt`**

        - **Content**: A tab-delimited matrix where rows are gene IDs and columns are samples. Values are the raw read counts assigned to each linear gene.
        - **Application**: Used as input for the `DE_Relative` analysis (to calculate total host gene counts) and can also be used for standard linear mRNA differential expression analysis (e.g., with DESeq2 or edgeR).

      - **`ciri3_multi/results.txt.BSJ_Matrix`**

        - **Content**: A tab-delimited matrix where rows are circRNA IDs (e.g., `chr1:100|200`) and columns are samples. Values are the raw counts of back-spliced junction (BSJ) reads supporting each circRNA.
        - **Application**: The primary input file for `DE_BSJ` (differential expression of circRNAs) analysis.

      - **`ciri3_multi/results.txt.FSJ_Matrix`**

        - **Content**: A tab-delimited matrix in the same format as the BSJ_Matrix, but values are the counts of forward-spliced junction (FSJ) reads corresponding to each circRNA's region.
        - **Application**: Used in conjunction with the BSJ_Matrix for `DE_Ratio` (differential splicing ratio) analysis to model circRNA splicing efficiency.

      - **`ciri3_de/DE_BSJ_results.txt`**

        - **Content**: Differential expression results from the `edgeR` QLF test. Columns include `circRNA_ID`, `logFC`, `PValue`, and `FDR`.
        - **Application**: Identifies circRNAs that are differentially expressed based on their BSJ counts, normalized for total library size.

      - **`ciri3de/DE_Ratio_results.txt`**

        - **Content**: Differential ratio results from the quasibinomial GLM model. Columns include `circRNA_ID`, `Beta` (log-odds ratio), `Mean_Ratio_C1` (Group 1 mean ratio), `Mean_Ratio_C2` (Group 2 mean ratio), `Delta` (difference in ratios), `PValue`, and `FDR`.
        - **Application**: Identifies *differential usage* of a circRNA. It models the ratio of BSJ / (BSJ + FSJ), which is largely independent of the host gene's expression level.

      - **`ciri3_de/DE_Relative_results.txt`**

        - **Content**: Differential relative expression results from the quasibinomial GLM model. Columns include `circRNA_ID`, `gene_id` (host gene), `Beta`, `Mean_Rel_C1`, `Mean_Rel_C2`, `Delta`, `PValue`, and `FDR`.
        - **Application**: Identifies circRNAs that change in abundance *relative to their host gene*. It models the ratio of a specific circRNA's BSJ count to the total BSJ counts of all circRNAs from that same host gene (circ_i / (sum(all_circs_from_gene))).

      - **`multiqc_report.html`** : Open `multiqc_report.html` in a web browser to explore all sections interactively.

        - **General Statistics**: A combined table summarizing important metrics for each sample (e.g., FastQC quality, % GC, Total Sequences, TrimGalore trimmed reads, featureCounts assignment rate).
          <img width="1948" height="1068" alt="CleanShot 2025-11-06 at 18 45 10@2x" src="https://github.com/user-attachments/assets/6dcf5576-dfcb-40e8-be66-b725db4136ce" />

	  
        - **FastQC**: Quality-control metrics on raw and trimmed reads, including 'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores', 'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content', 'Sequence Length Distribution', 'Sequence Duplication Levels', 'Adapter Content'.
          <img width="1922" height="1436" alt="CleanShot 2025-11-06 at 18 44 40@2x" src="https://github.com/user-attachments/assets/44b203c1-1b5f-4b17-a12a-6365d6e03182" />

	      
        - **Trim Galore!**: Reports the percentage and number of adapter sequences trimmed from each sample.
          <img width="1956" height="1288" alt="CleanShot 2025-11-06 at 18 44 57@2x" src="https://github.com/user-attachments/assets/38d104ce-831b-4a49-a6d6-f1b153830ecb" />

	  
        - **featureCounts**: Visualizes the percentage of reads assigned to genomic features (e.g., Exon, Intron, Intergenic). A high "Assigned" percentage is desirable.
          <img width="1912" height="1312" alt="CleanShot 2025-11-06 at 18 44 24@2x" src="https://github.com/user-attachments/assets/b1289371-d1bd-45a6-84d7-06840aff1599" />


# Reference
Zheng, X., Zhang, J., Song, L. et al. Detecting and quantifying circular RNAs in terabyte-scale RNA-seq datasets with CIRI3. Nat Biotechnol (2025). https://doi.org/10.1038/s41587-025-02835-1
