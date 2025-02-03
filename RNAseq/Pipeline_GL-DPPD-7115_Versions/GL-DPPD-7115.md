# GeneLab bioinformatics processing pipeline for Illumina RNA-sequencing data

> **This page holds an overview and instructions for how GeneLab processes RNAseq datasets. Exact processing commands, GL-DPPD-7101 version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** January 28, 2025  
**Document Number:** GL-DPPD-7115   

**Submitted by:**  
Alexis Torres (GeneLab Data Processing Team)  
Crystal Han (GeneLab Data Processing Team)

**Approved by:**  
Barbara Novak (GeneLab Data Processing Lead)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Samrawit Gebre (GeneLab Project Manager)  
Danielle Lopez (GeneLab Deputy Project Manager)  
Lauren Sanders (GeneLab Project Scientist)

---

## Updates from previous version  

This initial release of the RNAseq pipeline documentation details the workflow steps executed when the `--microbes` parameter is used. In this mode, Bowtie 2 is used for sequence alignment and featureCounts is used for gene-level quantification, replacing the default STAR and RSEM which are typically used for eukaryotic organisms.

Differences with default workflow:
- STAR is replaced with Bowtie 2 for alignment
- RSEM is replaced with featureCounts for gene quantification
- kentUtils gtfToGenePred and genePredToBed are replaced with gtfToBed.py
- rRNA genes are removed from featureCounts results on a dataset-wide basis. rRNA removal logs are all reported in the same file.
- Rather than importing RSEM Genes.Results files for each sample using tximport, the entire FeatureCounts table is imported into R for normalization and DGE analysis.


---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [**2. Trim/Filter Raw Data and Trimmed Data QC**](#2-trimfilter-raw-data-and-trimmed-data-qc)
    - [2a. Trim/Filter Raw Data](#2a-trimfilter-raw-data)
    - [2b. Trimmed Data QC](#2b-trimmed-data-qc)
    - [2c. Compile Trimmed Data QC](#2c-compile-trimmed-data-qc)
  - [**3. Build Bowtie 2 Reference**](#3-build-bowtie-2-reference)
  - [**4. Align Reads to Reference Genome then Sort and Index**](#4-align-reads-to-reference-genome-then-sort-and-index)
    - [4a. Align Reads to Reference Genome with Bowtie 2](#4a-align-reads-to-reference-genome-with-bowtie-2)
    - [4b. Compile Alignment Logs](#4b-compile-alignment-logs)
    - [4c. Sort Aligned Reads](#4c-sort-aligned-reads)
    - [4d. Index Sorted Aligned Reads](#4d-index-sorted-aligned-reads)
  - [**5. Create Reference BED File**](#5-create-reference-bed-file)
  - [**6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC**](#6-assess-strandedness-genebody-coverage-inner-distance-and-read-distribution-with-rseqc)
    - [6a. Determine Read Strandedness](#6a-determine-read-strandedness)
    - [6b. Compile Strandedness Reports](#6b-compile-strandedness-reports)
    - [6c. Evaluate GeneBody Coverage](#6c-evaluate-genebody-coverage)
    - [6d. Compile GeneBody Coverage Reports](#6d-compile-genebody-coverage-reports)
    - [6e. Determine Inner Distance (For Paired End Datasets)](#6e-determine-inner-distance-for-paired-end-datasets-only)
    - [6f. Compile Inner Distance Reports](#6f-compile-inner-distance-reports)
    - [6g. Assess Read Distribution](#6g-assess-read-distribution)
    - [6h. Compile Read Distribution Reports](#6h-compile-read-distribution-reports)
  - [**7. Quantitate Aligned Reads**](#7-quantitate-aligned-reads)
    - [7a. Count Aligned Reads with FeatureCounts](#7a-count-aligned-reads-with-featurecounts)
    - [7b. Compile FeatureCounts Logs](#7b-compile-featurecounts-logs)
    - [7c. Calculate Total Number of Genes Expressed Per Sample in R](#7c-calculate-total-number-of-genes-expressed-per-sample-in-r)
    - [7d. Remove rRNA Genes from FeatureCounts](#7d-remove-rrna-genes-from-featurecounts)
      - [7d.1 Extract rRNA Gene IDs from GTF](#7d1-extract-rrna-gene-ids-from-gtf)
      - [7d.2 Filter rRNA Genes from FeatureCounts](#7d2-filter-rrna-genes-from-featurecounts)
  - [**8. Normalize Read Counts and Perform Differential Gene Expression Analysis**](#8-normalize-read-counts-and-perform-differential-gene-expression-analysis)
    - [8a. Create Sample RunSheet](#8a-create-sample-runsheet)
    - [8b. Environment Set Up](#8b-environment-set-up)
    - [8c. Configure Metadata, Sample Grouping, and Group Comparisons](#8c-configure-metadata-sample-grouping-and-group-comparisons)
    - [8d. Import FeatureCounts Data](#8d-import-featurecounts-data)
    - [8e. Perform DGE Analysis](#8e-perform-dge-analysis)
    - [8f. Add Statistics and Gene Annotations to DGE Results](#8f-add-statistics-and-gene-annotations-to-dge-results)
    - [8g. Export DGE Tables](#8g-export-dge-tables)

  - [**9. Evaluate ERCC Spike-In Data**](#9-evaluate-ercc-spike-in-data)
    - [9a. Evaluate ERCC Count Data in Python](#9a-evaluate-ercc-count-data-in-python)
    - [9b. Perform DESeq2 Analysis of ERCC Counts in R](#9b-perform-deseq2-analysis-of-ercc-counts-in-r)
    - [9c. Analyze ERCC DESeq2 Results in Python](#9c-analyze-ercc-deseq2-results-in-python)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.12.1|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.26|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|4.2|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|0.6.10|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|Bowtie 2|2.5.4|[https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)|
|subread|2.0.8|[https://subread.sourceforge.net/](https://subread.sourceforge.net/)|
|Samtools|1.21|[http://www.htslib.org/](http://www.htslib.org/)|
|infer_experiment|5.0.4|[http://rseqc.sourceforge.net/#infer-experiment-py](http://rseqc.sourceforge.net/#infer-experiment-py)|
|geneBody_coverage|5.0.4|[http://rseqc.sourceforge.net/#genebody-coverage-py](http://rseqc.sourceforge.net/#genebody-coverage-py)|
|inner_distance|5.0.4|[http://rseqc.sourceforge.net/#inner-distance-py](http://rseqc.sourceforge.net/#inner-distance-py)|
|read_distribution|5.0.4|[http://rseqc.sourceforge.net/#read-distribution-py](http://rseqc.sourceforge.net/#read-distribution-py)|
|R|4.4.2|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.20|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|1.46.0|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tidyverse|2.0.0|[https://www.tidyverse.org](https://www.tidyverse.org)|
|stringr|1.5.1|[https://github.com/tidyverse/stringr](https://github.com/tidyverse/stringr)|
|dp_tools|1.3.5|[https://github.com/J-81/dp_tools](https://github.com/J-81/dp_tools)|
|pandas|2.2.3|[https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)|
|seaborn|0.13.2|[https://seaborn.pydata.org/](https://seaborn.pydata.org/)|
|matplotlib|3.10.0|[https://matplotlib.org/stable](https://matplotlib.org/stable)|
|jupyter notebook|7.3.2|[https://jupyter-notebook.readthedocs.io/](https://jupyter-notebook.readthedocs.io/)|
|numpy|2.2.1|[https://numpy.org/](https://numpy.org/)|
|scipy|1.15.1|[https://scipy.org/](https://scipy.org/)|
|singularity|3.9|[https://sylabs.io/](https://sylabs.io/)|

---

# General processing overview with example commands  

> Exact processing commands and output files listed in **bold** below are included with each RNAseq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). 

---

## 1. Raw Data QC

<br>

### 1a. Raw Data QC  

```bash
fastqc -o /path/to/raw_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- *fastqc.html (FastQC output html summary)
- *fastqc.zip (FastQC output data)

<br>

### 1b. Compile Raw Data QC  

```bash
multiqc --interactive -n raw_multiqc_GLbulkRNAseq -o /path/to/raw_multiqc/output/raw_multiqc_GLbulkRNAseq_report /path/to/directory/containing/raw_fastqc/files
zip -r raw_multiqc_GLbulkRNAseq_report.zip raw_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

* **raw_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **raw_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **raw_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

---

## 2. Trim/Filter Raw Data and Trimmed Data QC

<br>

### 2a. Trim/Filter Raw Data  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores NumberOfThreads \
  --phred33 \
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \ # only for PE studies, remove this parameter if raw data are SE
  sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz
# if SE, replace the last line with only the forward reads (R1) of each sample

```

**Parameter Definitions:**

- `--gzip` – compress the output files with `gzip`
- `--path_to_cutadapt` – specify path to cutadapt software if it is not in your `$PATH`
- `--cores` – specify the number of threads available on the server node to perform trimming
- `--phred33` – instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
- `--output_dir` – the output directory to store results
- `--paired` – indicates paired-end reads - both reads, forward (R1) and reverse (R2) must pass length threshold or else both reads are removed
- `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- **\*fastq.gz** (trimmed reads)
- **\*trimming_report.txt** (trimming report)

<br>

### 2b. Trimmed Data QC  

```bash
fastqc -o /path/to/trimmed_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**

- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *fastqc.html (FastQC output html summary)
- *fastqc.zip (FastQC output data)

<br>

### 2c. Compile Trimmed Data QC  

```bash
multiqc --interactive -n trimmed_multiqc_GLbulkRNAseq -o /path/to/trimmed_multiqc/output/trimmed_multiqc_GLbulkRNAseq_report /path/to/directory/containing/trimmed_fastqc/files
zip -r trimmed_multiqc_GLbulkRNAseq_report.zip /path/to/trimmed_multiqc/output/trimmed_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/trimmed_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 2b](#2b-trimmed-data-qc))

**Output Data:**

* **trimmed_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **trimmed_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **trimmed_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

---

## 3. Build Bowtie 2 Reference  

```bash
bowtie2-build --threads NumberOfThreads \
    -f /path/to/genome/fasta/file \
    bt2_base
```

**Parameter Definitions:**

- `--threads` – number of threads available on server node to create Bowtie 2 reference
- `-f` – specifies one or more fasta file(s) containing the genome reference sequences
- `bt2_base` – specifies the basename of of the Bowtie 2 index files to write

**Input Data:**

- *.fasta (genome sequence, this scRCP version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)


**Output Data:**

Bowtie 2 genome reference, which consists of the following files:

- bt2_base.1.bt2 (contains half of the forward strand BWT)
- bt2_base.2.bt2 (contains half of the forward strand BWT) 
- bt2_base.3.bt2 (contains BWT ranks)
- bt2_base.4.bt2 (contains BWT ranks)
- bt2_base.rev.1.bt2 (contains half of the reverse strand BWT)
- bt2_base.rev.2.bt2 (contains half of the reverse strand BWT)

<br>

---

## 4. Align Reads to Reference Genome then Sort and Index

<br>

### 4a. Align Reads to Reference Genome with Bowtie 2

```bash
bowtie2 -x /path/to/bowtie2/index \
 --threads NumberOfThreads \
 --minins 0 \
 --maxins 500 \
 -k 1 \
 -1 /path/to/trimmed_forward_reads \
 -2 /path/to/trimmed_reverse_reads \
 --un /path/to/bowtie2/output/directory/<sample_id>.unmapped.fastq \
 -S /path/to/bowtie2/output/directory/<sample_id>.sam \
 2> /path/to/bowtie2/output/directory/<sample_id>.bowtie2.log
```

**Parameter Definitions:**

- `-x` - specifies the path to the Bowtie2 index prefix
- `--threads` - number of threads to use for alignment
- `--minins` - minimum fragment length for valid paired-end alignments (0 for no minimum)
- `--maxins` - maximum fragment length for valid paired-end alignments
- `-k` - report up to k alignments per read (1 for unique mapping only)
- `-1` - path to input forward reads (R1)
- `-2` - path to input reverse reads (R2) (omit -1/-2 and use -U for single-end reads)
- `--un` - write unmapped reads to specified file
- `-S` - write alignments to SAM format file
- `2>` - redirect stderr containing alignment statistics to log file

**Input Data:**

- Bowtie2 genome reference (output from [Step 3](#3-build-bowtie-2-reference))
- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- **\*.sam** (alignments in SAM format)
- **\*.bowtie2.log** (log file containing alignment statistics)
- \*.unmapped.fastq (unmapped reads in FASTQ format)

<br>

### 4b. Compile Alignment Logs

```bash
multiqc --interactive -n align_multiqc_GLbulkRNAseq -o /path/to/align_multiqc/output/align_multiqc_GLbulkRNAseq_report /path/to/*.bowtie2.log/files
zip -r align_multiqc_GLbulkRNAseq_report.zip /path/to/align_multiqc/output/align_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*.bowtie2.log/files` – the directory holding the *.bowtie2.log output files from the [Bowtie2 alignment step](#4a-align-reads-to-reference-genome-with-bowtie-2), provided as a positional argument

**Input Data:**

- *.bowtie2.log (log file containing alignment info/stats such as reads mapped, etc., output from [Step 4a](#4a-align-reads-to-reference-genome-with-bowtie-2))

**Output Data:**

* **align_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **align_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **align_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

### 4c. Sort Aligned Reads

```bash
samtools sort -m 3G \
  --threads NumberOfThreads \
  -o /path/to/*_sorted.bam \
  /path/to/*.sam
```

**Parameter Definitions:**

- `-m` – memory available per thread, `3G` indicates 3 gigabytes, this can be changed based on user resources
- `--threads` – number of threads available on server node to sort genome alignment files
- `-o` - output file name (sorted BAM)
- `/path/to/*.sam` – path to the SAM files output from the [Bowtie2 alignment step](#4a-align-reads-to-reference-genome-with-bowtie-2), provided as a positional argument

**Input Data:**

- *.sam (alignments in SAM format, output from [Step 4a](#4a-align-reads-to-reference-genome-with-bowtie-2))

**Output Data:**

- **\*_sorted.bam** (samtools sorted genome aligned bam file)

<br>

### 4d. Index Sorted Aligned Reads

```bash
samtools index -@ NumberOfThreads /path/to/*_sorted.bam
```

**Parameter Definitions:**

- `-@` – number of threads available on server node to index the sorted alignment files
- `/path/to/*_sorted.bam` – the path to the sorted BAM files from the [sorting step](#4c-sort-aligned-reads), provided as a positional argument

**Input Data:**

- *_sorted.bam (sorted mapping to genome file, output from [Step 4c](#4c-sort-aligned-reads))

**Output Data:**

- **\*_sorted.bam.bai** (index of sorted mapping to genome file)

<br>

---

## 5. Create Reference BED File

```bash
gtf_to_bed.py \
  /path/to/annotation/gtf/file \
  /path/to/output/bed/file

```

**Parameter Definitions:**

- `/path/to/annotation/gtf/file` – specifies the file(s) containing annotated reference transcripts in the standard gtf format, provided as a positional argument
- `/path/to/output/bed/file` – specifies the location and name of the output BED file(s), provided as a positional argument

**Input Data:**

- *.gtf (genome annotation, this scRCP version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

**Output Data:**

- *.bed (genome annotation in BED format)

<br>

---

## 6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC

<br>

### 6a. Determine Read Strandedness

```bash
infer_experiment.py -r /path/to/annotation/BED/file \
 -i /path/to/*_sorted.bam \
 -s 15000000 > /path/to/*infer_expt.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `-s` – specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `>` – redirects standard output to specified file
- `/path/to/*infer_expt.out` – specifies the location and name of the file containing the infer_experiment standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5](#5-create-reference-bed-file))
- *_sorted.bam (sorted mapping to genome file, output from [Step 4c](#4c-sort-aligned-reads))
- *_sorted.bam.bai (index of sorted mapping to genome file, output from [Step 4d](#4d-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*_sorted.bam file)

**Output Data:**

- *infer_expt.out (file containing the infer_experiment standard output)

<br>

### 6b. Compile Strandedness Reports

```bash
multiqc --interactive -n infer_exp_multiqc_GLbulkRNAseq -o /path/to/infer_exp_multiqc/output/infer_exp_multiqc_GLbulkRNAseq_report /path/to/*infer_expt.out/files
zip -r infer_exp_multiqc_GLbulkRNAseq_report.zip /path/to/infer_exp_multiqc/output/infer_exp_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*infer_expt.out/files` – the directory holding the *infer_expt.out output files from the [read strandedness step](#6a-determine-read-strandedness), provided as a positional argument

**Input Data:**

- *infer_expt.out (file containing the infer_experiment standard output, output from [Step 6a](#6a-determine-read-strandedness))

**Output Data:**

* **infer_exp_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **infer_exp_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **infer_exp_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

### 6c. Evaluate GeneBody Coverage

```bash
geneBody_coverage.py -r /path/to/annotation/BED/file \
 -i /path/to/*_sorted.bam \
 -o /path/to/geneBody_coverage/output/directory/<sample_id>
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `-o` – specifies the path to the output directory
- `/path/to/geneBody_coverage/output/directory/<sample_id>` – specifies the location and name of the directory containing the geneBody_coverage output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5](#5-create-reference-bed-file))
- *_sorted.bam (sorted mapping to genome file, output from [Step 4c](#4c-sort-aligned-reads))
- *_sorted.bam.bai (index of sorted mapping to genome file, output from [Step 4d](#4d-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*_sorted.bam file)

**Output Data:**

- *.geneBodyCoverage.curves.pdf (genebody coverage line plot)
- *.geneBodyCoverage.r (R script that generates the genebody coverage line plot)
- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values used to generate the line plot)

<br>

### 6d. Compile GeneBody Coverage Reports

```bash
multiqc --interactive -n genebody_cov_multiqc_GLbulkRNAseq -o /path/to/geneBody_cov_multiqc/output/geneBody_cov_multiqc_GLbulkRNAseq_report /path/to/geneBody_coverage/output/files
zip -r genebody_cov_multiqc_GLbulkRNAseq_report.zip /path/to/genebody_cov_multiqc/output/genebody_cov_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/geneBody_coverage/output/files` – the directory holding the geneBody_coverage output files from [step 6c](#6c-evaluate-genebody-coverage), provided as a positional argument

**Input Data:**

- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values, output from [Step 6c](#6c-evaluate-genebody-coverage))

**Output Data:**

* **genebody_cov_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **genebody_cov_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **genebody_cov_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

### 6e. Determine Inner Distance (For Paired End Datasets ONLY)

```bash
inner_distance.py -r /path/to/annotation/BED/file \
 -i /path/to/*_sorted.bam \
 -k 15000000 \
 -l -150 \
 -u 350 \
 -o  /path/to/inner_distance/output/directory
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `-k` – specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `-l` – specifies the lower bound of inner distance (bp), set to -150 or negative of maximum read length if read length is greater than 150
- `-u` – specifies the upper bound of inner distance (bp)
- `/path/to/inner_distance/output/directory` – specifies the location and name of the directory containing the inner_distance output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5](#5-create-reference-bed-file))
- *_sorted.bam (sorted mapping to genome file, output from [Step 4c](#4c-sort-aligned-reads))
- *_sorted.bam.bai (index of sorted mapping to genome file, output from [Step 4d](#4d-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*_sorted.bam file)

**Output Data:**

- *.inner_distance.txt (log of read-wise inner distance results)
- *.inner_distance_freq.txt (tab delimited table of inner distances mapped to number of reads with that distance)
- *.inner_distance_plot.pdf (histogram plot of inner distance distribution)
- *.inner_distance_plot.r (R script that generates the histogram plot)

<br>

### 6f. Compile Inner Distance Reports

```bash
multiqc --interactive -n inner_dist_multiqc_GLbulkRNAseq /path/to/align_multiqc/output/inner_dist_multiqc_GLbulkRNAseq_report /path/to/inner_dist/output/files
zip -r inner_dist_multiqc_GLbulkRNAseq_report.zip /path/to/align_multiqc/output/inner_dist_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/inner_dist/output/files` – the directory holding the inner_distance output files from [Step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only), provided as a positional argument

**Input Data:**

- *.inner_distance_freq.txt (tab delimited table of inner distances from [step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only))

**Output Data:**

* **inner_dist_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **inner_dist_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **inner_dist_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

### 6g. Assess Read Distribution

```bash
read_distribution.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam > /path/to/*read_dist.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `>` – redirects standard output to specified file
- `/path/to/*read_dist.out` – specifies the location and name of the file containing the read_distribution standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5](#5-create-reference-bed-file))
- *_sorted.bam (sorted mapping to genome file, output from [Step 4c](#4c-sort-aligned-reads))
- *_sorted.bam.bai (index of sorted mapping to genome file, output from [Step 4d](#4d-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*_sorted.bam file)

**Output Data:**

- *read_dist.out (file containing the read_distribution standard output)

<br>

### 6h. Compile Read Distribution Reports

```bash
multiqc --interactive -n read_dist_multiqc_GLbulkRNAseq -o /path/to/read_dist_multiqc/output/read_dist_multiqc_GLbulkRNAseq_report /path/to/*read_dist.out/files
zip -r read_dist_multiqc_GLbulkRNAseq_report.zip /path/to/read_dist_multiqc/output/read_dist_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*read_dist.out/files` – the directory holding the *read_dist.out output files from [Step 6g](#6g-assess-read-distribution) provided as a positional argument

**Input Data:**

- *read_dist.out (files containing the read_distribution standard output, output from [Step 6g](#6g-assess-read-distribution))

**Output Data:**

* **read_dist_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **read_dist_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **read_dist_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

---

## 7. Quantitate Aligned Reads

<br>

### 7a. Count Aligned Reads with FeatureCounts

```bash
featureCounts -p \
  -T NumberOfThreads \
  -a /path/to/annotation/gtf/file \
  -s 0|1|2 \
  -t exon \
  -g gene_id \
  -o /path/to/featurecounts/output/directory/FeatureCounts_GLbulkRNAseq.csv \
  /path/to/*_sorted.bam
```

**Parameter Definitions:**

- `-p` - indicates paired-end reads (omit for single-end data)
- `-T` - number of threads to use
- `-a` - path to genome annotation GTF file
- `-s` - specifies strandedness: 0=unstranded, 1=stranded (forward), 2=stranded (reverse)
- `-t` - specify feature type (default: exon)
- `-g` - specify attribute type used to group features (default: gene_id)
- `-o` - output file path and name
- `/path/to/*_sorted.bam` - input sorted BAM files, can specify multiple files separated by spaces

**Input Data:**

- *.gtf (genome annotation, this version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)
- *_sorted.bam (sorted mapping to genome file, output from [Step 4c](#4c-sort-aligned-reads))
- *_sorted.bam.bai (index of sorted mapping to genome file, output from [Step 4d](#4d-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*_sorted.bam file)

**Output Data:**

- **FeatureCounts_GLbulkRNAseq.csv** (table containing raw read counts per gene for each sample)
- **FeatureCounts_GLbulkRNAseq.csv.summary** (table containing counting statistics for each sample)

<br>

### 7b. Compile FeatureCounts Logs

```bash
multiqc --interactive -n FeatureCounts_GLbulkRNAseq -o /path/to/FeatureCounts_multiqc/output/FeatureCounts_multiqc_GLbulkRNAseq_report /path/to/*stat/files
zip -r FeatureCounts_multiqc_GLbulkRNAseq_report.zip /path/to/FeatureCounts_multiqc/output/FeatureCounts_multiqc_GLbulkRNAseq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*stat/files` – the directories holding the *stat output files from the [RSEM Counts step](#8a-count-aligned-reads-with-rsem), provided as a positional argument

**Input Data:**

- *stat (directory containing the following stats files, output from [Step 7a](#7a-count-aligned-reads-with-featurecounts))
  - *cnt
  - *model
  - *theta

**Output Data:**

* **RSEM_count_multiqc_GLbulkRNAseq_report.zip** (zip containing the following)
  * **FeatureCounts_multiqc_GLbulkRNAseq.html** (multiqc output html summary)
  * **FeatureCounts_multiqc_GLbulkRNAseq_data** (directory containing multiqc output data)

<br>

### 7c. Calculate Total Number of Genes Expressed Per Sample in R

```R
### Load required packages
library(tidyverse)

### Set paths
counts_dir="/path/to/directory/containing/featurecounts/output"

### Read FeatureCounts data, skipping the first line which contains command info
data <- read.table(file.path(counts_dir, "FeatureCounts_GLbulkRNAseq.csv"), 
                  header=TRUE, 
                  skip=1,
                  sep="\t",
                  row.names=1)

### Remove metadata columns not needed for counting
data <- data %>% 
  select(-Chr, -Start, -End, -Strand, -Length)

### Remove '.bam' from column names if present
colnames(data) <- gsub("\\.bam$", "", colnames(data))

### Count number of genes with non-zero counts for each sample
NumNonZeroGenes <- colSums(data > 0)
NumNonZeroGenes <- as.data.frame(NumNonZeroGenes)
colnames(NumNonZeroGenes) <- "Number of genes with non-zero counts"

### Export results
write.csv(NumNonZeroGenes, 
          file=file.path(counts_dir, 'NumNonZeroGenes_GLbulkRNAseq.csv'),
          quote=TRUE)

### Print session info
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- FeatureCounts_GLbulkRNAseq.csv (table containing raw read counts per gene for each sample, output from [Step 7a](#7a-count-aligned-reads-with-featurecounts))

**Output Data:**

- **NumNonZeroGenes_GLbulkRNAseq.csv** (A samplewise table of the number of genes expressed)

<br>

### 7d. Remove rRNA Genes from FeatureCounts

<br>

#### 7d.1 Extract rRNA Gene IDs from GTF

```bash
### Extract unique rRNA ENSEMBL gene IDs from GTF file ###
grep "rRNA" /path/to/annotation/gtf/file \
    | grep -o 'gene_id "[^"]*"' \
    | sed 's/gene_id "\(.*\)"/\1/' \
    | sort -u > *_rrna_ensembl_ids.txt
```

**Input Data:**
- *.gtf (genome annotation)

**Output Data:**
- *rrna_ensembl_ids.txt (list of unique rRNA ENSEMBL gene IDs in a GTF file)

<br>

#### 7d.2 Filter rRNA Genes from FeatureCounts

```bash
### Filter out rRNA entries ###
awk 'NR==1; NR==2; NR==FNR {next} FNR==NR {ids[$1]=1; next} !($1 in ids)' \
    FeatureCounts_GLbulkRNAseq.csv \
    *rrna_ensembl_ids.txt \
    FeatureCounts_GLbulkRNAseq.csv > FeatureCounts_rRNA_removed_GLbulkRNAseq.csv

### Count removed rRNA entries ###
rRNA_count=$(awk 'NR==1; NR==2 {next} NR==FNR {next} FNR==NR {ids[$1]=1; next} $1 in ids' \
    FeatureCounts_GLbulkRNAseq.csv \
    *rrna_ensembl_ids.txt \
    FeatureCounts_GLbulkRNAseq.csv | wc -l)
echo "FeatureCounts: ${rRNA_count} rRNA entries removed." > FeatureCounts_rRNA_counts.txt
```

**Input Data:**
- FeatureCounts_GLbulkRNAseq.csv (table containing raw read counts per gene for each sample, output from [Step 7a](#7a-count-aligned-reads-with-featurecounts))
- *rrna_ensembl_ids.txt (file containing list of gene IDs with rRNA features, output from [Step 7d.1](#7d1-extract-rrna-gene-ids-from-gtf))

**Output Data:**
- **FeatureCounts_rRNA_removed_GLbulkRNAseq.csv** (FeatureCounts output for all samples with rRNA entries removed)
- *rRNA_counts.txt (Summary of number of rRNA entries removed)

<br>

---

## 8. Normalize Read Counts and Perform Differential Gene Expression Analysis

> **Note:** DGE Analysis is performed twice with different sets of input files:
> 1. Using raw FeatureCounts output
> 2. Using rRNA-removed FeatureCounts output

<br>

### 8a. Create Sample RunSheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/NF_RCP-F/examples/runsheet/README.md).

```bash
### Download the *ISA.zip file from the GeneLab Repository ###

dpt-get-isa-archive \
 --accession GLDS-###

### Parse the metadata from the *ISA.zip file to create a sample runsheet ###

dpt-isa-to-runsheet --accession GLDS-### \
 --config-type bulkRNASeq \
 --config-version Latest \
 --isa-archive *ISA.zip
```

**Parameter Definitions:**

- `--accession GLDS-###` – GLDS accession ID (replace ### with the GLDS number being processed), used to retrieve the urls for the ISA archive and raw reads hosted on the GeneLab Repository
- `--config-type` – Instructs the script to extract the metadata required for `bulkRNAseq` processing from the ISA archive
- `--config-version` – Specifies the `dp-tools` configuration version to use, a value of `Latest` will specify the most recent version
- `--isa-archive` – Specifies the *ISA.zip file for the respective GLDS dataset, downloaded in the `dpt-get-isa-archive` command


**Input Data:**

- No input data required but the GLDS accession ID needs to be indicated, which is used to download the respective ISA archive 

**Output Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, used to define sample groups - the *ISA.zip file is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'Study Files' -> 'metadata')

- **{GLDS-Accession-ID}_bulkRNASeq_v{version}_runsheet.csv** (table containing metadata required for processing, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

### 8b. Environment Set Up

```R
### Install and load required packages ###

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# List of required packages
cran_packages <- c("stringr", "knitr", "yaml", "dplyr")
bioc_packages <- c("tximport", "DESeq2", "BiocParallel")

# Install missing packages
for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
    }
}
for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg)
    }
}

# Load all packages
library(stringr)
library(knitr)
library(yaml)
library(dplyr)
library(tximport)
library(DESeq2)
library(BiocParallel)

### Define which organism is used in the study - this should be consistent with the name in the "name" column of the GL-DPPD-7110_annotations.csv file, which matches the abbreviations used in the Panther database for each organism ###

organism <- "organism_that_samples_were_derived_from"


### Define the location of the input data and where the output data will be printed to ###

runsheet_path="/path/to/directory/containing/runsheet.csv/file" ## This is the runsheet created in Step 9a above
work_dir="/path/to/working/directory/where/script/is/executed/from" 
counts_dir="/path/to/directory/containing/RSEM/counts/files"
norm_output="/path/to/normalized/counts/output/directory"
DGE_output="/path/to/DGE/output/directory"


### Pull in the GeneLab annotation table (GL-DPPD-7110_annotations.csv) file ###

org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv"

org_table <- read.table(org_table_link, sep = ",", header = TRUE)


### Define the link to the GeneLab annotation table for the organism of interest ###

annotations_link <- org_table[org_table$name == organism, "genelab_annots_link"]


### Set your working directory to the directory where you will execute your DESeq2 script from ###

setwd(file.path(work_dir))

```

<br>

### 8c. Configure Metadata, Sample Grouping, and Group Comparisons

```R
### Pull all factors for each sample in the study from the runsheet created in Step 9a ###

compare_csv_from_runsheet <- function(runsheet_path) {
    df = read.csv(runsheet_path)
    # get only Factor Value columns
    factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
    colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")
    result = data.frame(sample_id = df[,c("Sample.Name")], factors)	
    return(result)
}


### Load metadata from runsheet csv file ###

compare_csv <- compare_csv_from_runsheet(runsheet_path)


### Create data frame containing all samples and respective factors ###

study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
rownames(study) <- compare_csv[,1]


### Format groups and indicate the group that each sample belongs to ###

if (dim(study) >= 2){
    group<-apply(study,1,paste,collapse = " & ") ## concatenate multiple factors into one condition per sample
} else{
    group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") ## human readable group names
group <- sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", group))) # group naming compatible with R models, this maintains the default behaviour of make.names with the exception that 'X' is never prepended to group names
names(group) <- group_names
rm(group_names)


### Format contrasts table, defining pairwise comparisons for all groups ###

contrast.names <- combn(levels(factor(names(group))),2) ## generate matrix of pairwise group combinations for comparison
contrasts <- apply(contrast.names, MARGIN=2, function(col) sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))) # limited make.names call for each group (also removes leading parentheses)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) ## format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 

```

<br>

### 8d. Import FeatureCounts Data

```R
### Import FeatureCounts data ###
counts_file <- "/path/to/FeatureCounts_GLbulkRNAseq.csv"

# Load featureCounts data
featurecounts_data <- read.csv(file = counts_file, 
                              header = TRUE, 
                              sep = "\t", 
                              skip = 1, 
                              stringsAsFactors = FALSE, 
                              check.names = FALSE)

# Identify metadata columns and sample columns
metadata_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_cols <- setdiff(colnames(featurecounts_data), metadata_cols)

# Remove the ".bam" suffix from sample columns
sample_cols <- sub("\\.bam$", "", sample_cols)

# Reorder sample columns to match the sample order in the study
samples <- rownames(study)
sample_col_indices <- match(samples, sample_cols)

# Create counts matrix
counts <- featurecounts_data[, sample_col_indices, drop = FALSE]
counts <- as.data.frame(lapply(counts, as.numeric))
colnames(counts) <- samples
rownames(counts) <- featurecounts_data$Geneid

# Create DESeq2 input object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sampleTable,
    design = ~condition
)
```

<br>

### 8e. Perform DGE Analysis

```R
### Create sample table ###
sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)

### Handle technical replicates in sample table ###
# Replicates are identified by the presence of a "_techrepX" suffix in the sample name.
# Replicates are removed from the sampleTable such that the minimal number of equal technical replicates are kept across all samples.
handle_technical_replicates <- function(sampleTable) {
    # Extract base names and tech rep numbers
    sample_info <- data.frame(
        full_name = rownames(sampleTable),
        condition = sampleTable$condition,
        stringsAsFactors = FALSE
    )
    
    # Identify technical replicates
    tech_rep_pattern <- "_techrep(\\d+)$"
    has_tech_rep <- grepl(tech_rep_pattern, sample_info$full_name)
    
    if (!any(has_tech_rep)) {
        return(sampleTable)  # No technical replicates found
    }
    
    # Extract base names and tech rep numbers
    sample_info$base_name <- sub(tech_rep_pattern, "", sample_info$full_name)
    sample_info$tech_rep <- as.numeric(sub(".*_techrep", "", sample_info$full_name))
    sample_info$tech_rep[!has_tech_rep] <- 1
    
    # Group samples by base name and keep track of which ones to keep
    samples_to_keep <- c()
    unique_base_names <- unique(sample_info$base_name)
    
    for (base_name in unique_base_names) {
        base_samples <- which(sample_info$base_name == base_name)
        if (length(base_samples) == 1) {
            # Single sample, keep as is
            samples_to_keep <- c(samples_to_keep, base_samples)
        } else {
            # Multiple samples (tech reps), keep only the minimum number
            tech_rep_samples <- base_samples[order(sample_info$tech_rep[base_samples])]
            samples_to_keep <- c(samples_to_keep, tech_rep_samples[1])  # Keep only first tech rep
        }
    }
    
    # Create new sample table with only the kept samples
    new_sampleTable <- data.frame(
        condition = sample_info$condition[samples_to_keep],
        row.names = sample_info$base_name[samples_to_keep],
        stringsAsFactors = FALSE
    )
    
    return(new_sampleTable)
}

# Apply the technical replicate handling
sampleTable <- handle_technical_replicates(sampleTable)

# Update the counts matrix to match the new sample table
txi.rsem$counts <- txi.rsem$counts[, rownames(sampleTable)]

### Build dds object ###
dds <- DESeqDataSetFromTximport(
    txi = txi.rsem,
    colData = sampleTable,
    design = ~condition
)

### Collapse technical replicates if present ###
if (any(grepl("_techrep\\d+$", rownames(sampleTable)))) {
    tech_rep_groups <- sub("_techrep\\d+$", "", rownames(sampleTable))
    dds <- collapseReplicates(dds, tech_rep_groups)
}

### Filter low count genes ###
keep <- rowSums(counts(dds)) > 10
print(sprintf("Removed %d genes with dataset-wide count sum less than 10", sum(!keep)))
dds <- dds[keep,]

### Remove ERCC spike-in genes if present ###
if (length(grep("ERCC-", rownames(dds))) != 0) {
    dds <- dds[-c(grep("ERCC-", rownames(dds))), ]
}

### Perform DESeq analysis ###
dds <- DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)

### Generate normalized counts ###
normCounts <- as.data.frame(counts(dds, normalized = TRUE))
VSTCounts <- as.data.frame(assay(vst(dds)))

### Add 1 to normalized counts for log calculations ###
normCounts <- normCounts + 1

### Calculate LRT statistics ###
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res_lrt <- results(dds_lrt)

```

<br>

### 8f. Add Statistics and Gene Annotations to DGE Results

```R
### Initialize output table with normalized counts ###
output_table <- tibble::rownames_to_column(normCounts, var = "ENSEMBL")

### Add LRT p-values ###
output_table$LRT.p.value <- res_lrt@listData$padj

### Iterate through Wald Tests to generate pairwise comparisons of all groups ###
compute_contrast <- function(i) {
    res_1 <- results(
        dds_1,
        contrast = c("condition", contrasts[1, i], contrasts[2, i]),
        parallel = FALSE  # Disable internal parallelization
    )
    res_1_df <- as.data.frame(res_1@listData)[, c(2, 4, 5, 6)]
    colnames(res_1_df) <- c(
        paste0("Log2fc_", colnames(contrasts)[i]),
        paste0("Stat_", colnames(contrasts)[i]),
        paste0("P.value_", colnames(contrasts)[i]),
        paste0("Adj.p.value_", colnames(contrasts)[i])
    )
    return(res_1_df)
}

### Use bplapply to compute results in parallel ###
res_list <- bplapply(1:dim(contrasts)[2], compute_contrast, BPPARAM = BPPARAM)

### Combine the list of data frames into a single data frame ###
res_df <- do.call(cbind, res_list)

### Combine with the existing output_table ###
output_table <- cbind(output_table, res_df)

### Add summary statistics ###
output_table$All.mean <- rowMeans(normCounts, na.rm = TRUE)
output_table$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE)
output_table$LRT.p.value <- res_1_lrt@listData$padj

### Add group-wise statistics ###
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- names(group)

# Calculate group means and standard deviations
group_means <- as.data.frame(t(aggregate(. ~ group, data = tcounts, mean)))
group_stdev <- as.data.frame(t(aggregate(. ~ group, data = tcounts, sd)))

# Remove group name rows
group_means <- group_means[-1,]
group_stdev <- group_stdev[-1,]

# For each group, add mean and stdev columns
for (group_name in names(group)) {
    mean_col <- paste0("Group.Mean_(", group_name, ")")
    stdev_col <- paste0("Group.Stdev_(", group_name, ")")
    output_table[[mean_col]] <- group_means[, paste0("Group.Mean_", group_means['group',])]
    output_table[[stdev_col]] <- group_stdev[, paste0("Group.Stdev_", group_stdev['group',])]
}

### Read in GeneLab annotation table for the organism of interest ###
annot <- read.table(annotations_link, 
    sep = "\t", 
    header = TRUE, 
    quote = "", 
    comment.char = "", 
    row.names = 1
)

### Combine annotations table and the DGE table ###
output_table <- merge(annot, output_table, by='row.names', all.y=TRUE)
output_table <- output_table %>% 
  rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )

```

### 8g. Export DGE Tables

```R
### Export unnormalized and normalized counts tables ###
write.csv(txi.rsem$counts, 
    file.path(norm_output, "RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"))

write.csv(normCounts,
    file.path(norm_output, "Normalized_Counts_GLbulkRNAseq.csv"))

### Export sample grouping and contrasts tables ###
write.csv(sampleTable,
    file.path(DGE_output, "SampleTable_GLbulkRNAseq.csv"))

write.csv(contrasts,
    file.path(DGE_output, "contrasts_GLbulkRNAseq.csv"))

### Export DGE Results table ###
write.csv(output_table,
    file.path(DGE_output, "differential_expression_GLbulkRNAseq.csv"), 
    row.names = FALSE)

### print session info ###

print("Session Info below: ")
sessionInfo()

```


**Input Data:**

- *runsheet.csv file (table containing metadata required for analysis, output from [step 9a](#9a-create-sample-runsheet))
- [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) (csv file containing link to GeneLab annotations) 
- FeatureCounts_GLbulkRNAseq.csv (table containing raw read counts per gene for each sample, output from [Step 7a](#7a-count-aligned-reads-with-featurecounts))

**Output Data:**

- **RSEM_Unnormalized_Counts_GLbulkRNAseq.csv** (table containing raw RSEM gene counts for each sample)
- **Normalized_Counts_GLbulkRNAseq.csv** (table containing normalized gene counts for each sample)
- **VST_Counts_GLbulkRNAseq.csv** (table containing VST normalized gene counts for each sample)
- **SampleTable_GLbulkRNAseq.csv** (table containing samples and their respective groups)
- **differential_expression_GLbulkRNAseq.csv** (table containing normalized counts for each sample, group statistics, DESeq2 DGE results for each pairwise comparison, and gene annotations) 
- **contrasts_GLbulkRNAseq.csv** (table containing all pairwise comparisons)

> Note: Datasets with technical replicates are handled by collapsing them such that the minimum number of equal technical replicates is retained across all samples. Before normalization, the counts of technical replicates are summed to combine them into a single sample representing the biological replicate.

> Note: RNAseq processed data interactive tables and plots are found in the [GLDS visualization portal](https://visualization.genelab.nasa.gov/data/studies).

<br>

---

## 9. Evaluate ERCC Spike-In Data 

> Note: This is only applicable for datasets with ERCC spike-in

### 9a. Evaluate ERCC Count Data in Python

```python
### Setting up the notebook

# import python packages
import pandas as pd
pd.set_option('mode.chained_assignment', None) # suppress chained indexing warnings
import numpy as np
from json import loads
from re import search
import zipfile
import seaborn as sns
from scipy.stats import linregress
import matplotlib.pyplot as plt


### Get and parse data and metadata

# Get and unzip ISA.zip to extract metadata.

accession = 'GLDS-NNN' # Replace Ns with GLDS number
isaPath = '/path/to/GLDS-NNN_metadata_GLDS-NNN-ISA.zip' # Replace with path to ISA archive file
zip_file_object =  zipfile.ZipFile(isaPath, "r")
list_of_ISA_files = zip_file_object.namelist() # Print contents of zip file. Pick relevant one from list
UnnormalizedCountsPath = '/path/to/GLDS-NNN_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv'

GENE_ID_PREFIX = "ENSMU" # change according to a common prefix for all Gene IDs associated with the subject organism

# Print contents of ISA zip file to view file order
list_of_ISA_files

# There are datasets that have multiple assays (including microarray), so the RNAseq ISA files from the above output must be selected. 
# Txt files outputted above are indexed as 0, 1, 2, etc. Fill in the indexed number corresponding to the sample (s_*txt) and assay files for RNAseq (a_*_(RNA-Seq).txt) in the code block below.

# Extract metadata from the sample file (s_*txt)
sample_file = list_of_ISA_files[2] # replace [2] with index corresponding to the (s_*txt) file
file = zip_file_object.open(sample_file)
sample_table = pd.read_csv(zip_file_object.open(sample_file), sep='\t')

# Extract metadata from the assay (a_*_(RNA-Seq).txt) file
assay_file = list_of_ISA_files[0] # replace [0] with index corresponding to the (a_*_(RNA-Seq).txt) file
file = zip_file_object.open(assay_file)
assay_table = pd.read_csv(zip_file_object.open(assay_file), sep='\t')

# Check the sample table
pd.set_option('display.max_columns', None)
print(sample_table.head(n=3))

# Check the assay table
pd.set_option('display.max_columns', None)
assay_table.head(n=3)

### Get raw counts table

raw_counts_table = pd.read_csv(UnnormalizedCountsPath, index_col=0) 
raw_counts_table.index.rename('Gene_ID', inplace=True)
print(raw_counts_table.head(n=3))

raw_counts_transcripts = raw_counts_table[raw_counts_table.index.str.contains(f"^{GENE_ID_PREFIX}")]
assert len(raw_counts_transcripts) != 0, f"Looks like {GENE_ID_PREFIX} matched no genes, probably the wrong prefix"
raw_counts_transcripts = raw_counts_transcripts.sort_values(by=list(raw_counts_transcripts), ascending=False)
print(raw_counts_transcripts)

### Get ERCC counts

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts.reset_index(inplace=True)
ercc_counts = ercc_counts.rename(columns={'Gene_ID':'ERCC ID'})
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())

### Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
ercc_table = pd.read_csv(ercc_url, sep = '\t')
print(ercc_table.head(n=3))

## Calculate the number of ERCC genes detected in each of the 4 (A, B, C and D) groups for each sample
### Extract ERCC counts and calculate the log(2)

meltERCC = ercc_counts.melt(id_vars=['ERCC ID'])
meltERCC['log2 Count'] = meltERCC['value']+1
meltERCC['log2 Count'] = np.log2(meltERCC['log2 Count'])
meltERCC = meltERCC.rename(columns={'variable':'Sample Name', 'value':'Count'})
print(meltERCC.head(n=3))

### Build Mix dictionary to link sample name to mix added and read depth using the assay table

mix_dict = assay_table.filter(['Sample Name','Parameter Value[Spike-in Mix Number]', 
                       'Parameter Value[Read Depth]'])
mix_dict = mix_dict.rename(columns={'Parameter Value[Spike-in Mix Number]':'Mix',
                                    'Parameter Value[Read Depth]':
                                    'Total Reads'})
print(mix_dict.head(n=3))

### Make combined ercc counts and assay table

merged_ercc = meltERCC.merge(mix_dict, on='Sample Name')
print(merged_ercc)

### Read ERCC info including concentrations from merged_ercc table

groupA = ercc_table.loc[ercc_table['subgroup'] == 'A']['ERCC ID']
groupB = ercc_table.loc[ercc_table['subgroup'] == 'B']['ERCC ID']
groupC = ercc_table.loc[ercc_table['subgroup'] == 'C']['ERCC ID']
groupD = ercc_table.loc[ercc_table['subgroup'] == 'D']['ERCC ID']

### Make a dictionary for ERCC groups

group_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['subgroup']))

### Calculate ERCC counts per million and log(2) counts per million

merged_ercc['Count per million'] = merged_ercc['Count'] / (merged_ercc['Total Reads'] / 1000000.0)
merged_ercc['log2 Count per million'] = np.log2(merged_ercc['Count per million']+1)

### Add ERCC group column

merged_ercc['ERCC group'] = merged_ercc['ERCC ID'].map(group_dict)
merged_ercc = merged_ercc.sort_values(by=['Mix'], ascending=True)
print(merged_ercc)

## Filter and calculate mean counts per million of Mix1 and Mix2 spiked samples in each of the 4 groups
### Filter Mix1 CPM and Mix2 CPM in group A 

Adf = merged_ercc.loc[merged_ercc['ERCC group'] == 'A']
Amix1df = Adf.loc[Adf['Mix']=='Mix 1']
Amix1df['Mix1 CPM'] = Amix1df[Amix1df['Count per million'] > 0]['Count per million'].dropna()
Amix1df = Amix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Amix1df = Amix1df.to_frame()
Amix2df = Adf.loc[Adf['Mix']=='Mix 2']
Amix2df['Mix2 CPM'] = Amix2df[Amix2df['Count per million'] > 0]['Count per million'].dropna()
Amix2df = Amix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Amix2df = Amix2df.to_frame()

adf = Amix1df.merge(Amix2df, on='ERCC ID', suffixes=('', '_2'))
adf = adf.reset_index()
adf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (adf['Avg Mix1 CPM'] / adf['Avg Mix2 CPM'])

### Filter Mix1 CPM and Mix2 CPM in group B

Bdf = merged_ercc.loc[merged_ercc['ERCC group'] == 'B']
Bmix1df = Bdf.loc[Bdf['Mix']=='Mix 1']
Bmix1df['Mix1 CPM'] = Bmix1df[Bmix1df['Count per million'] > 0]['Count per million'].dropna()
Bmix1df = Bmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Bmix1df = Bmix1df.to_frame()
Bmix2df = Bdf.loc[Bdf['Mix']=='Mix 2']
Bmix2df['Mix2 CPM'] = Bmix2df[Bmix2df['Count per million'] > 0]['Count per million'].dropna()
Bmix2df = Bmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Bmix2df = Bmix2df.to_frame()

bdf = Bmix1df.merge(Bmix2df, on='ERCC ID')
bdf = bdf.reset_index()
bdf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (bdf['Avg Mix1 CPM'] / bdf['Avg Mix2 CPM'])

### Filter Mix1 CPM and Mix2 CPM in group C

Cdf = merged_ercc.loc[merged_ercc['ERCC group'] == 'C']
Cmix1df = Cdf.loc[Cdf['Mix']=='Mix 1']
Cmix1df['Mix1 CPM'] = Cmix1df[Cmix1df['Count per million'] > 0]['Count per million'].dropna()
Cmix1df = Cmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Cmix1df = Cmix1df.to_frame()
Cmix2df = Cdf.loc[Cdf['Mix']=='Mix 2']
Cmix2df['Mix2 CPM'] = Cmix2df[Cmix2df['Count per million'] > 0]['Count per million'].dropna()
Cmix2df = Cmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Cmix2df = Cmix2df.to_frame()

cdf = Cmix1df.merge(Cmix2df, on='ERCC ID')
cdf = cdf.reset_index()
cdf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (cdf['Avg Mix1 CPM'] / cdf['Avg Mix2 CPM'])

### Filter Mix1 CPM and Mix2 CPM in group D

Ddf = merged_ercc.loc[merged_ercc['ERCC group'] == 'D']
Dmix1df = Ddf.loc[Ddf['Mix']=='Mix 1']
Dmix1df['Mix1 CPM'] = Dmix1df[Dmix1df['Count per million'] > 0]['Count per million'].dropna()
Dmix1df = Dmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Dmix1df = Dmix1df.to_frame()
Dmix2df = Ddf.loc[Ddf['Mix']=='Mix 2']
Dmix2df['Mix2 CPM'] = Dmix2df[Dmix2df['Count per million'] > 0]['Count per million'].dropna()
Dmix2df = Dmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Dmix2df = Dmix2df.to_frame()

ddf = Dmix1df.merge(Dmix2df, on='ERCC ID')
ddf = ddf.reset_index()
ddf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (ddf['Avg Mix1 CPM'] / ddf['Avg Mix2 CPM'])

## Multi-sample ERCC analyses
### Create box and whisker plots of the log(2) CPM for each ERCC detected in group A in Mix 1 and Mix 2 spiked samples

a = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupA[::-1], hue="Mix",data=merged_ercc[merged_ercc['ERCC ID'].isin(groupA)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
a.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 4")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group A ERCC genes (for group A we expect Mix 1 CPM / Mix 2 CPM = 4)

a1 = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupA[::-1], palette="rocket_r", data=adf, kind="bar", height=5, aspect=1, linewidth=0.5)
a1.set_xticklabels(rotation=90)
plt.title("ERCC Group A")
a1.set(ylim=(0, 6))
a1.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group A (out of 23) =', adf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group B in Mix 1 and Mix 2 spiked samples

b = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupB[::-1], hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupB)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
b.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 1")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group B ERCC genes (for group B we expect Mix 1 CPM / Mix 2 CPM = 1)

b = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupB[::-1], palette="rocket_r", data=bdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
b.set_xticklabels(rotation=90)
plt.title("ERCC Group B")
b.set(ylim=(0, 2))
b.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group B (out of 23) =', bdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group C in Mix 1 and Mix 2 spiked samples

c = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupC[::-1], hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupC)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
c.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.67")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group C ERCC genes (for group C we expect Mix 1 CPM / Mix 2 CPM = 0.67)

c = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupC[::-1], palette="rocket_r", data=cdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
c.set_xticklabels(rotation=90)
plt.title("ERCC Group C")
c.set(ylim=(0, 2))
c.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group C (out of 23) =', cdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group D in Mix 1 and Mix 2 spiked samples

d = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupD[::-1], hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupD)], col="ERCC group", kind="box", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
d.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.5")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group D ERCC genes (for group D we expect Mix 1 CPM / Mix 2 CPM = 0.5)

d = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupD[::-1], palette="rocket_r", data=ddf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
d.set_xticklabels(rotation=90)
plt.title("ERCC Group D")
d.set(ylim=(0, 1))
d.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group D (out of 23) =', ddf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

## Individual sample ERCC analyses
# Calculate and plot ERCC metrics from individual samples, including limit of detection, dynamic range, and R^2 of counts vs. concentration.

#Calculate and plot ERCC metrics from individual samples, including limit of detection, dynamic range, and R^2 of counts vs. concentration.
print(ercc_table.head(n=3))

# Make a dictionary for ERCC concentrations for each mix

mix1_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 1 (attomoles/ul)']))
mix2_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 2 (attomoles/ul)']))

# Check assay_table header to identify the 'Sample Name' column and the column title indicating the 'Spike-in Mix Nmber' if it's indicated in the metadata.

pd.set_option('display.max_columns', None)
print(assay_table.head(n=3))

# Get samples that use mix 1 and mix 2

mix1_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 1']['Sample Name']
mix2_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 2']['Sample Name']

# Get ERCC counts for all samples

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())

# Get ERCC counts for Mix 1 spiked samples

ercc_counts_mix_1 = ercc_counts[mix1_samples]
ercc_counts_mix_1['ERCC conc (attomoles/ul)'] = ercc_counts_mix_1.index.map(mix1_conc_dict)
print(ercc_counts_mix_1.head(n=3))

# Get ERCC counts for Mix 2 spiked samples

ercc_counts_mix_2 = ercc_counts[mix2_samples]
ercc_counts_mix_2['ERCC conc (attomoles/ul)'] = ercc_counts_mix_2.index.map(mix2_conc_dict)
print(ercc_counts_mix_2.head(n=3))

# Create a scatter plot of log(2) ERCC counts versus log(2) ERCC concentration for each sample

columns_mix_1 = ercc_counts_mix_1.columns.drop(['ERCC conc (attomoles/ul)'])
columns_mix_2 = ercc_counts_mix_2.columns.drop(['ERCC conc (attomoles/ul)'])
all_columns = columns_mix_1.to_list() + columns_mix_2.to_list()
total_columns = len(columns_mix_1) + len(columns_mix_2) 
side_size = np.int32(np.ceil(np.sqrt(total_columns)))# calculate grid side size. take sqrt of total plots and round up.
fig, axs = plt.subplots(side_size, side_size, figsize=(22,26), sharex='all', sharey='all'); #change figsize x,y labels if needed.
fig.tight_layout(pad=1, w_pad=2.5, h_pad=3.5)

# Iterate over subplot positions, if subplot does not have data, hide the subplot with ax.set_visible(False)

for i in range(side_size):
    for j in range(side_size):
        ax = axs[i, j]
        index = i * side_size + j

        if index < total_columns:
            if index < len(columns_mix_1):
                ax.scatter(x=np.log2(ercc_counts_mix_1['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_1[all_columns[index]]+1), s=7)
                ax.set_title(all_columns[index][-45:], fontsize=9)
            else:
                ax.scatter(x=np.log2(ercc_counts_mix_2['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_2[all_columns[index]]+1), s=7)
                ax.set_title(all_columns[index][-45:], fontsize=9)

            ax.set_xlabel('log2 ERCC conc (attomoles/ ul)', fontsize=9)
            ax.set_ylabel('log2 Counts per million', fontsize=9)
            ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True)
        else:
            ax.set_visible(False)  # Hide the subplot if it's not needed

plt.show()

## Calculate and plot linear regression of log(2) ERCC counts versus log(2) ERCC concentration for each sample

# Filter counts > 0

nonzero_counts_list_1 = []
for i in range(0, len(ercc_counts_mix_1.columns)-1):
  counts = ercc_counts_mix_1[columns_mix_1[i]]
  counts.index.rename('Gene_ID', inplace=True)
  countsdf = pd.DataFrame(counts)
  nonzero_counts = countsdf[ercc_counts_mix_1[columns_mix_1[i]] > 0.0]
  nonzero_counts['Conc'] = nonzero_counts.index.map(mix1_conc_dict)
  nonzero_counts.columns = ['Counts','Conc']
  nonzero_counts_sorted = nonzero_counts.sort_values('Conc')
  nonzero_counts_list_1.append(nonzero_counts_sorted)

nonzero_counts_list_2 = []
for i in range(0, len(ercc_counts_mix_2.columns)-1):
  counts = ercc_counts_mix_2[columns_mix_2[i]]
  counts.index.rename('Gene_ID', inplace=True)
  countsdf = pd.DataFrame(counts)
  nonzero_counts = countsdf[ercc_counts_mix_2[columns_mix_2[i]] > 0.0]
  nonzero_counts['Conc'] = nonzero_counts.index.map(mix2_conc_dict)
  nonzero_counts.columns = ['Counts','Conc']
  nonzero_counts_sorted = nonzero_counts.sort_values('Conc')
  nonzero_counts_list_2.append(nonzero_counts_sorted)


# Plot each sample using linear regression of scatter plot with x = log2 Conc and y = log2 Counts.
# Return min, max, R^2 and dynamic range (max / min) values.

samples = []
mins = []
maxs = []
dyranges = []
rs = []

fig, axs = plt.subplots(side_size, side_size, figsize=(22,26), sharex='all', sharey='all');
fig.tight_layout(pad=1, w_pad=2.5, h_pad=3.5)

counter = 0
list2counter = 0
for i in range(side_size):
    for j in range(side_size):
        ax = axs[i, j]
        index = i * side_size + j

        if index < len(columns_mix_1):
            nonzero_counts = nonzero_counts_list_1[index]
            xvalues = nonzero_counts['Conc']
            yvalues = nonzero_counts['Counts']

            if len(xvalues) > 0:
                sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax)
                ax.set_title(all_columns[index][-47:], fontsize=9)
                ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9)
                ax.set_ylabel('log2 Counts per million', fontsize=9)
                ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True)
                samples.append(all_columns[index])

                min_val = xvalues.min()
                max_val = xvalues.max()
                dynamic_range = max_val / min_val
                slope, intercept, r_value, p_value, std_err = linregress(np.log2(xvalues), np.log2(yvalues))

                ax.text(0.02, 0.98, f'Min:{min_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.88, f'Max:{max_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.78, f'Dyn:{dynamic_range:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.68, f'R:{r_value:.2f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)

                mins.append(min_val)
                maxs.append(max_val)
                dyranges.append(dynamic_range)
                rs.append(r_value)
            else:
                ax.set_visible(False)

        elif index < total_columns:
            nonzero_counts = nonzero_counts_list_2[list2counter]
            xvalues = nonzero_counts['Conc']
            yvalues = nonzero_counts['Counts']

            if len(xvalues) > 0:
                sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax)
                ax.set_title(all_columns[index][-47:], fontsize=9)
                ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9)
                ax.set_ylabel('log2 Counts per million', fontsize=9)
                ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True)
                samples.append(all_columns[index])

                min_val = xvalues.min()
                max_val = xvalues.max()
                dynamic_range = max_val / min_val
                slope, intercept, r_value, p_value, std_err = linregress(np.log2(xvalues), np.log2(yvalues))

                ax.text(0.02, 0.98, f'Min:{min_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.88, f'Max:{max_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.78, f'Dyn:{dynamic_range:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.68, f'R:{r_value:.2f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)

                mins.append(min_val)
                maxs.append(max_val)
                dyranges.append(dynamic_range)
                rs.append(r_value)
            else:
                ax.set_visible(False)

            list2counter += 1
        else:
            ax.set_visible(False)

plt.show()

# Create directory for saved files

import os
os.makedirs(name="ERCC_analysis", exist_ok=True)

# Print tables containing the dynamic range and R^2 values for each sample.
# Remember to change file names to GLDS# analyzing

stats = pd.DataFrame(list(zip(samples, mins, maxs, dyranges, rs)))
stats.columns = ['Samples', 'Min', 'Max', 'Dynamic range', 'R']
stats.to_csv(f'ERCC_analysis/ERCC_stats_{accession}_GLbulkRNAseq.csv', index = False)
stats.filter(items = ['Samples', 'Dynamic range']).to_csv(f'ERCC_analysis/ERCC_dynrange_{accession}_mqc_GLbulkRNAseq.csv', index = False)
stats.filter(items = ['Samples', 'R']).to_csv(f'ERCC_analysis/ERCC_rsq_{accession}_mqc_GLbulkRNAseq.csv', index = False)

## Generate data and metadata files needed for ERCC DESeq2 analysis

# ERCC Mix 1 and Mix 2 are distributed so that half the samples receive Mix 1 spike-in and half receive Mix 2 spike-in. Transcripts in Mix 1 and Mix 2 are present at a known ratio, so we can determine how well these patterns are revealed in the dataset.

# Get sample table

combined = sample_table.merge(assay_table, on='Sample Name')
combined = combined.set_index(combined['Sample Name'])
pd.set_option('display.max_columns', None)
print(combined)

# Create metadata table containing samples and their respective ERCC spike-in Mix number
# Sometimes Number in [Spike-in Mix Number] is spelled 'number' and this could cause error in mismatch search 

ERCCmetadata = combined[['Parameter Value[Spike-in Mix Number]']]
ERCCmetadata.index = ERCCmetadata.index.str.replace('-','_')
ERCCmetadata.columns = ['Mix']
#ERCCmetadata = ERCCmetadata.rename(columns={'Parameter Value[Spike-in Mix Number]':'Mix'})
print(ERCCmetadata)

# Export ERCC sample metadata

ERCCmetadata.to_csv('ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv') 

# Export ERCC count data

ercc_counts.columns = ercc_counts.columns.str.replace('-','_')
ERCCcounts = ercc_counts.loc[:,ERCCmetadata.index]
ERCCcounts.head()

ERCCcounts.to_csv('ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv') 
```


**Input Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, output from [Step 9a](#9a-create-sample-runsheet))
- RSEM_Unnormalized_Counts.csv (RSEM raw counts table, output from [Step 9](#ERCCspikeOut))

**Output Data:**

- ERCC_analysis/ERCC_stats_GLDS-*_GLbulkRNAseq.csv (Samplewise counts statistics table containing 'Min', 'Max', 'Dynamic range', 'R')
- ERCC_analysis/ERCC_dynrange_GLDS-*_GLbulkRNAseq.csv (Samplewise counts statistics subset table containing 'Dynamic range')
- ERCC_analysis/ERCC_rsq_GLDS-*_GLbulkRNAseq.csv (Samplewise counts statistics subset table containing 'R')
- ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv (Samplewise metadata table inlcuding ERCC mix number)
- ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv (Samplewise ERCC counts table)

<br>

### 9b. Perform DESeq2 Analysis of ERCC Counts in R

```R

## Install R packages if not already installed

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

## Import DESeq2 library

library("DESeq2")

## Import and format ERCC count data and metadata

cts <- as.matrix(read.csv('ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv',sep=",",row.names="Gene_ID")) #INPUT
coldata <- read.csv('ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv', row.names=1) #INPUT

coldata$Mix <- factor(coldata$Mix)
all(rownames(coldata) == colnames(cts))

## Make DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Mix)
dds

## Filter out ERCC genes with counts of less than 10 in all samples #####

keepGenes <- rowSums(counts(dds)) > 10
dds <- dds[keepGenes,]

dds

## Run DESeq2 analysis and calculate results

## Try first to use the default type="median", but if there is an error (usually due to zeros in genes), use type="poscounts"
## From DESeq2 manual: "The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts."
dds <- tryCatch(
      expr = { estimateSizeFactors(dds) },
      error = function(e) { estimateSizeFactors(dds, type="poscounts")}
)


dds <- DESeq(dds)
res <- results(dds, contrast=c("Mix","Mix 1","Mix 2")) # remove space before mix number if needed
res

## Export DESeq2 results table and normalized ERCC counts table

write.csv(res, 'ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv') #OUTPUT
normcounts = counts(dds, normalized=TRUE)
write.csv(normcounts, 'ERCC_analysis/ERCC_normcounts_GLbulkRNAseq.csv') #OUTPUT
```

**Input Data:**

- ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv (samplewise metadata table inlcuding ERCC mix number, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))
- ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv (samplewise ERCC counts table, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))

**Output Data:**

- ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv (DESeq2 results table)
- ERCC_analysis/ERCC_normcounts_GLbulkRNAseq.csv (normalized ERCC Counts table)

<br>

### 9c. Analyze ERCC DESeq2 Results in Python

```python

# Import python packages

## Import python packages

import pandas as pd
from urllib.request import urlopen, quote, urlretrieve
import seaborn as sns
import matplotlib.pyplot as plt

accession = 'GLDS-NNN' # Replace Ns with GLDS number

## Import ERCC DESeq2 results

deseq2out = pd.read_csv('ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv', index_col=0) # INPUT
#deseq2out.index = deseq2out.index.str.replace('_','-')
deseq2out.rename(columns ={'baseMean' : 'meanNormCounts'}, inplace = True)
print(deseq2out.head())

## Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
ercc_table = pd.read_csv(ercc_url, sep='\t', index_col='ERCC ID')
print(ercc_table.head(n=3))

## Combine ERCC DESeq2 results and ercc_table

combined = deseq2out.merge(ercc_table, left_index=True, right_index=True)
print(combined.head())

## Filter p-value and adj. p-value cutoff at 10^-3

combined['cleaned_padj'] = combined['padj']
combined.loc[(combined.cleaned_padj < 0.001),'cleaned_padj']=0.001

combined['cleaned_pvalue'] = combined['pvalue']
combined.loc[(combined.cleaned_pvalue < 0.001),'cleaned_pvalue']=0.001

print(combined.head())

## Export the filtered combined ERCC DESeq2 results and ercc_table
### Remember to change file name to GLDS# analyzing

combined.filter(items = ['ERCC ID', 'meanNormCounts', 'cleaned_pvalue','cleaned_padj']).to_csv(f'ERCC_analysis/ERCC_lodr_{accession}_mqc_GLbulkRNAseq.csv') 

## Plot p-value vs. mean normalized ERCC counts

fig, ax = plt.subplots(figsize=(10, 7))

sns.scatterplot(data=combined, x="meanNormCounts", y="cleaned_pvalue",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

sns.lineplot(data=combined, x="meanNormCounts", y="cleaned_pvalue",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

#g.set_xscale("log", base=2)
ax.set_xscale("linear");
ax.set_yscale("log");

## Plot Adjp-value vs. mean normalized ERCC counts

fig, ax = plt.subplots(figsize=(10, 7))

sns.scatterplot(data=combined, x="meanNormCounts", y="cleaned_padj",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

sns.lineplot(data=combined, x="meanNormCounts", y="cleaned_padj",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

#g.set_xscale("log", base=2)
ax.set_xscale("linear");
ax.set_yscale("log");
```

**Input Data:**

- ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv (ERCC DESeq2 results table, output from [Step 10b](#10b-perform-deseq2-analysis-of-ercc-counts-in-r))

**Output Data:**

- ERCC_analysis/ERCC_lodr_*_GLbulkRNAseq.csv (ERCC Gene Table including mean counts, adjusted p-value and p-value, and filtered to genes with both adj. p-value and p-value < 0.001)

> All steps of the ERCC Spike-In Data Analysis are performed in a Jupyter Notebook (JN) and the completed JN is exported as an html file (**ERCC_analysis_GLbulkRNAseq.html**) and published in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/) for the respective dataset.