# **ATAC-seq Snakemake Pipeline**

## **Overview**
This repository contains an **automated ATAC-seq analysis pipeline** built using **Snakemake**. The pipeline is designed for efficient and reproducible processing of **Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq)** data. It follows a structured workflow to **process raw FASTQ files, align reads, perform quality control, call peaks, and compare different peak calling methods**.

The pipeline allows for the **benchmarking of peak calling tools** by systematically comparing their outputs, including **peak overlap analysis and Jaccard similarity calculations**. This ensures that the most suitable peak caller can be selected for specific datasets.

### **Key Features**
- Fully automated workflow using **Snakemake**.
- Supports multiple peak calling methods (**MACS2, HMMRATAC, Genrich, SEACR**).
- Implements **quality control** at multiple stages, including **read filtering and duplicate removal**.
- Provides **comparison and benchmarking** of different peak calling tools to assess performance.
- Generates **bigWig tracks** for visualization.
- Configurable through a **config.yaml** file.

---

## **Pipeline Workflow**
The pipeline consists of the following major steps:

1. **Pre-processing and Quality Control**
   - Adapter trimming and quality filtering using **fastp**.
   - Initial sequence quality assessment using **FastQC**.
   - Reference genome indexing using **BWA**.

2. **Read Alignment and Post-processing**
   - Aligning reads to the reference genome using **BWA**.
   - Sorting and indexing BAM files using **samtools**.
   - Removing low-quality reads and duplicate reads using **Picard MarkDuplicates**.
   - Removing mitochondrial and blacklisted reads using **bedtools**.

3. **Peak Calling and Benchmarking**
   - Running multiple peak calling tools: **MACS2, HMMRATAC, Genrich, and SEACR**.
   - Computing **Fraction of Reads in Peaks (FRiP) Score**.
   - Computing **Jaccard similarity between peak sets**.
   - **Comparing the performance of peak callers through overlap analysis**.

4. **Signal Track Generation**
   - Generating **bigWig tracks** for visualization using **deepTools bamCoverage**.

5. **Comparison of Peak Callers**
   - Evaluating **peak overlaps** between different peak calling tools.
   - Calculating **Jaccard similarity scores** to assess the similarity of peak sets.
   - Counting the **total number of peaks called by each tool**.

---

## **Installation and Dependencies**
### Install Snakemake
Ensure Snakemake is installed:
```bash
conda install -c bioconda snakemake
```

### Clone the Repository
```bash
git clone https://github.com/your-repo/your-project.git
cd your-project
```

### Set Up a Conda Environment
```bash
conda env create -f environment.yml
conda activate atac-seq-pipeline
```

---

## **Pipeline Execution**
### Running the Full Workflow
To execute the entire pipeline using 8 CPU cores:
```bash
snakemake --cores 8 --use-conda
```

### Running a Specific Rule
To execute only a specific step:
```bash
snakemake --cores 4 output/sample_1.bam
```

### Generating a Workflow Visualization
```bash
snakemake --dag | dot -Tpng > workflow.png
```

---

## **Detailed Explanation of Each Step**

### **Pre-processing and Quality Control**
- Adapter trimming using `fastp` to remove sequencing adapters and low-quality bases.
- Quality control assessment using `FastQC` to inspect sequence quality.
- Indexing the reference genome using `BWA index`.

```bash
bwa index reference.fa
```

### **Read Alignment and Post-processing**
- Aligning reads using `BWA`.
- Sorting and indexing BAM files using `samtools`.
- Marking and removing duplicate reads using `Picard MarkDuplicates`.
- Filtering out low-quality reads and blacklisted regions using `bedtools`.

```bash
bwa mem -t 8 reference.fa trimmed_R1.fastq.gz trimmed_R2.fastq.gz | samtools sort -o sample.sorted.bam
picard MarkDuplicates I=sample.sorted.bam O=sample.dedup.bam M=metrics.txt REMOVE_DUPLICATES=true
```

### **Peak Calling and Benchmarking**
The pipeline runs multiple peak calling tools and compares their performance.

#### MACS2 (Standard ATAC-seq Peak Caller)
```bash
macs2 callpeak -t sample.filtered.bam -f BAMPE -g hs -n sample --outdir peaks/
```

#### HMMRATAC (Specialized ATAC-seq Peak Caller)
```bash
HMMRATAC -b sample.filtered.bam -i sample.filtered.bam.bai -g genome.info -o sample
```

#### Genrich (Optimized for Removing Duplicates)
```bash
Genrich -t sample.filtered.bam -o sample_peaks.narrowPeak
```

#### SEACR (Sparse Enrichment Peak Caller)
```bash
bash SEACR.sh sample.filtered.bam peaks/
```

---

## **Benchmarking and Comparison of Peak Callers**
To evaluate and benchmark peak calling methods, the pipeline performs the following analyses:

### **Peak Overlap Analysis**
Identifying peaks that are called by multiple tools:
```bash
bedtools intersect -a macs2_peaks.narrowPeak -b hmmratac_peaks.narrowPeak > macs2_vs_hmmratac.bed
```

### **Jaccard Similarity Calculation**
Quantifying the similarity between peak sets:
```bash
bedtools jaccard -a macs2_peaks.narrowPeak -b genrich_peaks.narrowPeak
```

### **Counting the Number of Peaks Detected**
Checking how many peaks are detected by each tool:
```bash
wc -l peaks/*.narrowPeak
```

---

## **Signal Track Generation**
Generating signal tracks for visualization:
```bash
bamCoverage --binSize 10 --normalizeUsing BPM --bam sample.filtered.bam -o sample.bw
```

---

## **Citations and References**
If you use this pipeline, please cite:
- **Snakemake**: Mölder et al., *Bioinformatics* (2021)
- **MACS2**: Zhang et al., *Genome Biology* (2008)
- **HMMRATAC**: Tarbell et al., *Nucleic Acids Research* (2019)
- **Genrich**: Developed by John Stamatoyannopoulos Lab
- **SEACR**: Meers et al., *Genome Biology* (2019)
- **BWA**: Li & Durbin, *Bioinformatics* (2009)
- **FeatureCounts**: Liao et al., *Bioinformatics* (2014)
- **Bedtools**: Quinlan & Hall, *Bioinformatics* (2010)

---

## **Contact Information**
For questions or bug reports, open an issue on GitHub or contact:
📧 **happystone@**

---

## **Next Steps**
- Verify the `config.yaml` file
- Run the pipeline with test data
- Analyze and interpret results

This pipeline is continuously updated. If you have suggestions or feature requests, feel free to contribute!
