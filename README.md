# **ATAC-seq Snakemake Pipeline**
Benchmarking methods of ATAC-seq pipeline

## **Overview**
This repository contains an **automated ATAC-seq analysis pipeline** built with **Snakemake**, designed for efficient and reproducible processing of **Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq)** data. The pipeline performs **FASTQ processing, alignment, quality control, peak calling, and benchmarking of different peak calling methods**.

### **Pipeline Workflow**
The pipeline consists of the following key steps:

1. **Pre-alignment processing**
   - **Quality Control (QC):** Adapter trimming and quality filtering
   - **Reference genome indexing:** Preparing genome index for alignment

2. **Alignment, QC, and track visualization**
   - **Read Alignment:** Mapping trimmed reads to the reference genome
   - **Post-alignment QC:** Filtering low-quality reads and removing duplicates
   - **Signal track generation:** Converting BAM files to bigWig for visualization

3. **Peak calling and quality control**
   - **Peak calling:** MACS2, HMMRATAC, Genrich, and SEACR
   - **Quality assessment:** FRiP score, Jaccard similarity, and peak overlap analysis

4. **Differential accessibility and motif analysis**
   - **Differential analysis:** Identifying differential accessibility regions
   - **Motif analysis:** Discovering transcription factor binding motifs

---

## **Installation and Dependencies**
### **1️⃣ Install Snakemake**
Ensure you have Snakemake installed:
```bash
conda install -c bioconda snakemake
```

### **2️⃣ Clone this repository**
```bash
git clone https://github.com/your-repo/your-project.git
cd your-project
```

### **3️⃣ Set up a Conda environment**
Create a Conda environment for this pipeline:
```bash
conda env create -f environment.yml
conda activate atac-seq-pipeline
```

---

## **Pipeline Execution**
### **Run the full workflow**
To execute the full pipeline with 8 CPU cores:
```bash
snakemake --cores 8 --use-conda
```

### **Run a specific rule**
If you want to execute only a specific step, specify the output file:
```bash
snakemake --cores 4 output/sample_1.bam
```

### **Visualize the workflow**
To generate a **workflow DAG visualization**:
```bash
snakemake --dag | dot -Tpng > workflow.png
```

---

## **Pipeline Steps**

### **1️⃣ Quality Control (Pre-alignment QC)**
#### FastQC Report
The raw sequencing reads are first assessed for quality using `fastqc`. Run the following command manually to inspect the quality:
```bash
fastqc <sample>_R1.fastq.gz <sample>_R2.fastq.gz -o qc_reports/
```

#### Adapter Trimming
Adapters and low-quality bases are removed using `fastp`:
```bash
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -j <sample>.fastp.json -h <sample>.fastp.html
```

---

### **2️⃣ Read Alignment**
Reads are aligned to the reference genome using `BWA`:
```bash
bwa mem -t 8 reference/hg38.fa <sample>_R1.trimmed.fastq.gz <sample>_R2.trimmed.fastq.gz | samtools sort -o <sample>.sorted.bam
```

---

### **3️⃣ Post-alignment Quality Control**
#### Remove mitochondrial reads
```bash
samtools view -h <sample>.sorted.bam | grep -v "chrM" | samtools sort -o <sample>.rmChrM.bam
```

#### Mark and Remove Duplicates
```bash
picard MarkDuplicates I=<sample>.rmChrM.bam O=<sample>.dedup.bam M=<sample>.dup.metrics REMOVE_DUPLICATES=true
```

#### Remove Low-Quality Alignments
```bash
samtools view -h -b -q 30 <sample>.dedup.bam > <sample>.filtered.bam
samtools index <sample>.filtered.bam
```

---

### **4️⃣ Peak Calling**
#### MACS2 Peak Calling
```bash
macs2 callpeak -t <sample>.filtered.bam -f BAMPE -g hs -n <sample> --outdir peaks/
```

#### HMMRATAC Peak Calling
```bash
HMMRATAC -b <sample>.filtered.bam -i <sample>.filtered.bam.bai -g genome.info -o <sample>
```

---

### **5️⃣ Peak Quality Control**
#### Calculate FRiP Score
```bash
featureCounts -p -a <sample>_peaks.narrowPeak -F SAF -o <sample>-readCountInPeaks.txt <sample>.filtered.bam
```

---

### **6️⃣ Visualization: Generate Signal Tracks**
```bash
bamCoverage --binSize 10 --normalizeUsing BPM --bam <sample>.filtered.bam -o <sample>_coverage.bw
```

---

## **Citation & References**
If you use this pipeline, please cite:
- **Snakemake**: Mölder et al., *Bioinformatics* (2021)
- **MACS2**: Zhang et al., *Genome Biology* (2008)
- **BWA**: Li & Durbin, *Bioinformatics* (2009)

---

## **Contact**
For questions or bug reports, open an issue on GitHub or contact:  
📧 **your_email@example.com**

---

### 🔥 **Next Steps**
✅ Verify the `config.yaml` file  
✅ Run the workflow with test data  
✅ Analyze the results  

Let me know if you need further modifications! 🚀
