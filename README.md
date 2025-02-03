# ATAC-seq-pipline
Benchmarking methods of ATAC-seq pipline

# ATAC-seq Snakemake Pipeline

This repository contains an **automated ATAC-seq analysis pipeline** built with **Snakemake**, designed for efficient and reproducible processing of **Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq)** data. The pipeline performs **FASTQ processing, alignment, quality control, peak calling, and benchmarking of different peak calling methods** while also providing **batch correction strategies** for high-quality downstream analysis.

---

## **📌 Key Features**
### **Preprocessing & Alignment**
✅ **Quality Control & Adapter Trimming** (`fastp`)  
✅ **Read Alignment with BWA-MEM** (`bwa mem`)  

### **Post-Alignment QC & Filtering**
✅ **Remove mitochondrial reads** (`samtools idxstats`)  
✅ **Mark and remove PCR duplicates** (`picard MarkDuplicates`)  
✅ **Filter low-quality reads** (`samtools view -q 30`)  
✅ **Remove ENCODE Blacklist regions** (`bedtools intersect`)  

### **Peak Calling & Signal Track Generation**
✅ **Multiple peak calling tools** (`MACS2`, `Genrich`, `HMMRATAC`, `SEACR`)  
✅ **Generate signal tracks for visualization** (`bamCoverage`)  
✅ **Benchmark peak calling performance using peak overlap analysis** (`bedtools jaccard`)  

### **Quality Control Metrics**
✅ **Fraction of Reads in Peaks (FRiP) Score Calculation** (`featureCounts`)  
✅ **Fragment size distribution analysis** (`ATACseqQC`)  

### **Batch Effect Correction**
✅ **Compare batch correction methods** (`Signac LSI` vs. `ArchR IterativeLSI`)  

---

## **📥 Installation**
### **1️⃣ Clone the Repository**
```bash
git clone https://github.com/YOUR_USERNAME/ATAC-seq-pipeline.git
cd ATAC-seq-pipeline
