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
   - Initial sequence quality assessment using **FastQC**.
   - Adapter trimming and quality filtering using **fastp**.
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


### **1. Pre-alignment QC**
- Adapter trimming using `fastp` to remove sequencing adapters and low-quality bases.
- Quality control assessment using `FastQC` to inspect sequence quality.
- Indexing the reference genome using `BWA index`.

Independent replicates should be processed separately.

#### **Initial QC Report**

The raw sequence data should first be assessed for quality. [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content, and adapter contamination. In ATAC-seq data, it is likely that Nextera sequencing adapters will be over-represented. As described by [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3), base quality should be high although it may drop slightly at the 3' end, while GC content and read length should be consistent with expected values. For paired-end reads, run FastQC on both files, with the results output to the current directory:

```bash
fastqc <sample>_R1.fastq.gz -d . -o .
fastqc <sample>_R2.fastq.gz -d . -o .
```

#### **Adapter Trimming**

Adapters and low-quality reads/bases should be trimmed using one of several programs, such as [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), or [fastp](https://github.com/OpenGene/fastp). Adapter contamination can be seen in the FastQC report.

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/adapters.png" width="800">
For this pipeline, fastp is used to remove adapter sequences.

```bash
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -j <sample>.fastp.json -h <sample>.fastp.html
```

---

### **Reference Genome Indexing**
Before aligning reads, the reference genome must be indexed:
```bash
bwa index reference.fa
```

### **2. Read Alignment**

- Aligning reads using `BWA`.
- Sorting and indexing BAM files using `samtools`.
- Marking and removing duplicate reads using `Picard MarkDuplicates`.
- Filtering out low-quality reads and blacklisted regions using `bedtools`.

```bash
bwa mem -t 8 reference.fa trimmed_R1.fastq.gz trimmed_R2.fastq.gz | samtools sort -o sample.sorted.bam
picard MarkDuplicates I=sample.sorted.bam O=sample.dedup.bam M=metrics.txt REMOVE_DUPLICATES=true
```

### 3. Post-alignment QC

The post-alignment QC involves several steps:

- [Remove mitochondrial reads](#remove-mitochondrial-reads)
- [Remove duplicates & low-quality alignments](#remove-duplicates-&-low-quality-alignments) (including non-uniquely mapped reads)
- [Calculate library complexity and QC](#calculate-library-complexity-and-QC)
- [Remove ENCODE blacklist regions](#remove-encode-blacklist-regions)
- [Shift read coordinates](#shift-read-coordinates)

For an ATAC-seq experiment, the number of uniquely mapped reads ***after these steps*** is recommended to be 25 million for single-end or 50 million paired-end reads [(Buenrostro et al. 2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4374986/). Specific to ATAC-seq, an additional QC step is to check the fragment size distribution, which is expected to correspond to the length of nucleosomes:

- [Assess fragment size distribution](#assess-fragment-size-distribution)

#### Remove mitochondrial reads

ATAC-seq experiments commonly include a high proportion of mitochondrial reads. These should be removed. To assess the total % of mitochondrial reads, `samtools idxstats` can be run to report the total number of reads mapped to each chromosome. `samtools flagstat` provides a short report including the total number of DNA reads as the first line (halve this number for the total number of fragments). 

```bash
#Generate the idxstats report
samtools idxstats <sample>_sorted.bam > <sample>_sorted.idxstats

#Check the number of reads mapped to the mitochondria (chrM)
grep "chrM" <sample>_sorted.idxstats
```

The second column is the length of the chromosome and the third column is the total number of reads aligned to the chromosome (chrM). To see the total number of DNA fragments, run:

```bash
#Generate the flagstat report
samtools flagstat <sample>_sorted.bam > <sample>_sorted.flagstat

#check the total number of aligned fragments
head <sample>_sorted.flagstat
```

The % of DNA fragments aligned to chrM can be calculated as a % of the total DNA fragments. To remove any mitocondrial DNA, run the following:

```bash
#Remove reads aligned to the mitochondria
samtools view -h <sample>_sorted.bam | grep -v chrM | samtools sort -O bam -o <sample>.rmChrM.bam -T .
```

#### Mark duplicates 

To mark duplicate reads and view the % of duplicates:

```bash
picard MarkDuplicates QUIET=true INPUT=<sample>.rmChrM.bam OUTPUT=<sample>.marked.bam METRICS_FILE=<sample>.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

#View the % of duplicates
head -n 8 <sample>.dup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2
```

#### Remove duplicates & low-quality alignments 

The output `sam/bam` files contain several measures of quality. First, the alignment quality score. Reads which are uniquely mapped are assigned a high alignment quality score and one genomic position. If reads can map to more than one location, Bowtie2 reports one position and assigns a low quality score. The proportion of uniquely mapped reads can be assessed. In general, >70% uniquely mapped reads is expected, while <50% may be a cause for concern [(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf). Secondly, the 'flag' reports information such as whether the read is mapped as a pair or is a PCR duplicate. The individual flags are reported [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html) and are combined in a `sam/bam` file to one score, which can be deconstructed back to the original flags using [online interpretation tools](https://broadinstitute.github.io/picard/explain-flags.html). In this pipeline, the bowtie2 parameters `--no-mixed` and `--no-discordant` prevent the mapping of only one read in a pair, so these flags will not be present. All flags reported in a `sam` file can optionally be viewed using  `grep -v ^@ <sample>.sam | cut -f 2 | sort | uniq`.

**Multi-mapping:** the user should decide whether or not to retain reads which have multi-mapped, i.e. aligned to more than one position in the reference genome. When using paired-end data, it may be the case that one read aligns to a repetitive region (and therefore can map elsewhere), while the mate aligns to a unique sequence with a high quality. The bowtie2 parameters used above required reads to align within 50-700bp, so there should be no reads incorrectly aligned outside this distance. As such, the user may decide to keep multi-mapping reads on the assumption that they are likely to be mapped to the correct sequence, within the length of the DNA fragment. This may, however, cause incorrect alignments in extended repetitive regions where a read could map to multiple positions within the length of the DNA fragment. This should be minimised by the downstream removal of the [ENCODE Blacklisted regions](https://www.nature.com/articles/s41598-019-45839-z).

If a read is multi-mapped, it is assigned a low quality score by bowtie2. To view how many DNA reads align with a quality score >30, run the following (divide this number by 2 to calculate the # of DNA fragments):

```bash
samtools view -q 30 -c <sample>.marked.bam
```

A low % of uniquely mapped reads map result from short reads, excessive PCR amplification or problems with the PCR (Bailey et al. 2013). The following code uses the `sam/bam` flags to retain properly mapped pairs (`-f 2`) and to remove reads which fail the platform/vendor QC checks (`-F 512`), duplicate reads (`-F 1024`) and those which are unmapped (`-F 12`). The three flags to be removed can be combined into `-F 1548`, which will remove reads which meet any of the three individual flags

Here we will ***remove*** multi-mapped reads:

```bash
samtools view -h -b -f 2 -F 1548 -q 30 <sample>.marked.bam | samtools sort -o <sample>.filtered.bam

samtools index <sample>.filtered.bam

#To retain multi-mapped reads:
#samtools view -h -b -f 2 -F 1548 <sample>.rmChrM.bam | samtools sort -n -o <sample>.filtered.bam 
```

#### Remove ENCODE blacklist regions

The [ENCODE blacklist regions](https://github.com/Boyle-Lab/Blacklist/), most recently reported by [Amemiya et al. (2019)](https://www.nature.com/articles/s41598-019-45839-z) are defined as 'a comprehensive set of regions in the human, mouse, worm, and fly genomes that have anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment.' These problematic regions should be removed before further analysis. Download the blacklist files for your chosen reference genome from the [Boyle Lab github repository](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). Details regarding the identification of blacklist regions are reported [here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist-README.pdf).

```bash
#Remove reads within the blacklist regions
bedtools intersect -nonamecheck -v -abam <sample>.filtered.bam -b hg19-blacklist.v2.bed > <sample>.tmp.bam

#Sort and index the bam file
samtools sort -O bam -o <sample>.blacklist-filtered.bam <sample>.tmp.bam
samtools index <sample>.blacklist-filtered.bam

rm <sample>.tmp.bam
```

#### Shift read coordinates

An optional step in analysing data generated using the Tn5 transposase (such as ATAC-seq, ChIPmentation etc.) is to account for a small DNA insertion, introducted as repair of the transposase-induced nick introduces a 9bp insertion. Reads aligning to the + strand should be offset by +4bp and reads aligned to the -ve strand should be offset by -5bp. For references, see the first ATAC-seq paper by [Buenrostro et al., (2013)](https://www.nature.com/articles/nmeth.2688) and the analysis by [Adey et al., (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r119) which showed this insertion bias. Shifting coordinates is only really important if single-base resolution is required, for example in the analysis of transcription factor motifs in ATAC-seq peak footprints. Be aware that some tools do this shifting themselves (so double check manuals!).

We can use the deeptools package. 

```bash
alignmentSieve --numberOfProcessors max --ATACshift --blackListFileName hg38-blacklist.v2.bed --bam blacklist-filtered.bam -o ${prefix}.shifted.bam

# samtools index ${prefix}.shifted.bam
```







#### **4. Peak Calling and Benchmarking**
- Running multiple peak calling tools: **MACS2, HMMRATAC, Genrich, and SEACR**.
- Computing **Fraction of Reads in Peaks (FRiP) Score**.
- Computing **Jaccard similarity between peak sets**.
- **Comparing the performance of peak callers through overlap analysis**.

##### MACS2 (Standard ATAC-seq Peak Caller)
```bash
macs2 callpeak -t sample.filtered.bam -f BAMPE -g hs -n sample --outdir peaks/
```

##### HMMRATAC (Specialized ATAC-seq Peak Caller)
```bash
HMMRATAC -b sample.filtered.bam -i sample.filtered.bam.bai -g genome.info -o sample
```

##### Genrich (Optimized for Removing Duplicates)
```bash
Genrich -t sample.filtered.bam -o sample_peaks.narrowPeak
```

##### SEACR (Sparse Enrichment Peak Caller)
```bash
bash SEACR.sh sample.filtered.bam peaks/
```

---

#### **5. Benchmarking and Comparison of Peak Callers**
- Evaluating **peak overlaps** between different peak calling tools.
- Calculating **Jaccard similarity scores** to assess the similarity of peak sets.
- Counting the **total number of peaks called by each tool**.

##### **Peak Overlap Analysis**
Identifying peaks that are called by multiple tools:
```bash
bedtools intersect -a macs2_peaks.narrowPeak -b hmmratac_peaks.narrowPeak > macs2_vs_hmmratac.bed
```

##### **Jaccard Similarity Calculation**
Quantifying the similarity between peak sets:
```bash
bedtools jaccard -a macs2_peaks.narrowPeak -b genrich_peaks.narrowPeak
```

##### **Counting the Number of Peaks Detected**
Checking how many peaks are detected by each tool:
```bash
wc -l peaks/*.narrowPeak
```

---

#### **Signal Track Generation**
Generating signal tracks for visualization:
```bash
bamCoverage --binSize 10 --normalizeUsing BPM --bam sample.filtered.bam -o sample.bw
```

---

#### **Citations and References**
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

#### **Contact Information**
For questions or bug reports, open an issue on GitHub or contact:
📧 **suki5976@yuhs.ac**

---

#### **Next Steps**
- Verify the `config.yaml` file
- Run the pipeline with test data
- Analyze and interpret results

This pipeline is continuously updated. If you have suggestions or feature requests, feel free to contribute!
