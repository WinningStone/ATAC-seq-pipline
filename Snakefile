# Load configuration file
configfile: "config.yaml"

# Load variables from config
GENOME = config["genome"]
OUTPUT_DIR = config["output_dir"]
PEAK_CALLERS = config["peak_callers"]
GROUPS = list(config["samples"].keys())
SAMPLES = [(group, sample) for group in GROUPS for sample in config["samples"][group]]
BLACKLIST = config["blacklist"]

# Rule to generate all expected output files
rule all:
    input:
        [f"{OUTPUT_DIR}/{group}/{sample}_{caller}_peaks.narrowPeak"
         for group in GROUPS
         for sample in config["samples"][group]
         for caller in PEAK_CALLERS],

        [f"{OUTPUT_DIR}/{group}/{sample}.bigWig"
         for group in GROUPS
         for sample in config["samples"][group]],

        [f"benchmarking/{group}/{sample}_peak_jaccard.txt"
         for group in GROUPS
         for sample in config["samples"][group]],

        [f"{OUTPUT_DIR}/{group}/{sample}_{caller}_FRiP.txt"
         for group in GROUPS
         for sample in config["samples"][group]
         for caller in PEAK_CALLERS]

# Reference genome indexing
rule Index_reference:
    output:
        f"genome_index/{GENOME}.bwt"
    shell:
        "bwa index {GENOME}"


# Fastq preprocess

rule Pre_alignment_QC:
    input:
        fastq1 = lambda wildcards: f"data/{wildcards.group}/{wildcards.sample}_1.fastq.gz",
        fastq2 = lambda wildcards: f"data/{wildcards.group}/{wildcards.sample}_2.fastq.gz"
    output:
        trimmed_fastq1 = f"{OUTPUT_DIR}/{{group}}/{{sample}}_1.trimmed.fastq.gz",
        trimmed_fastq2 = f"{OUTPUT_DIR}/{{group}}/{{sample}}_2.trimmed.fastq.gz",
        json_report = f"{OUTPUT_DIR}/{{group}}/{{sample}}.fastp.json",
        html_report = f"{OUTPUT_DIR}/{{group}}/{{sample}}.fastp.html"
    shell:
        """
        fastp -i {input.fastq1} -I {input.fastq2} \
              -o {output.trimmed_fastq1} -O {output.trimmed_fastq2} \
              --detect_adapter_for_pe \
              -j {output.json_report} -h {output.html_report}
        """

# Read Alignment
rule Align_reads:
    input:
        fastq1 = f"{OUTPUT_DIR}/{{group}}/{{sample}}_1.trimmed.fastq.gz",
        fastq2 = f"{OUTPUT_DIR}/{{group}}/{{sample}}_2.trimmed.fastq.gz",
        index = f"genome_index/{GENOME}.bwt"
    output:
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input.index} {input.fastq1} {input.fastq2} | samtools view -Sb - > {output.bam}"


#remove mito reads
rule remove_mito_reads:
    input:
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.bam"
    output:
        bam_rmChrM = f"{OUTPUT_DIR}/{{group}}/{{sample}}.rmChrM.bam"
    shell:
        "samtools view -h {input.bam} | grep -v chrM | samtools sort -O bam -o {output.bam_rmChrM}"


# Mark Duplicates
rule mark_duplicates:
    input:
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.rmChrM.bam"
    output:
        marked_bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.marked.bam",
        metrics = f"{OUTPUT_DIR}/{{group}}/{{sample}}.dup.metrics"
    shell:
        """
        picard MarkDuplicates QUIET=true \
            INPUT={input.bam} \
            OUTPUT={output.marked_bam} \
            METRICS_FILE={output.metrics} \
            REMOVE_DUPLICATES=false \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT
        """


# Filter Low-Quality Reads
rule filter_low_quality:
    input:
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.marked.bam"
    output:
        filtered_bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.filtered.bam"
    shell:
        """
        samtools view -h -b -f 2 -F 1548 -q 30 {input.bam} | samtools sort -o {output.filtered_bam}
        samtools index {output.filtered_bam}
        """

# Remove Blacklist Regions
rule remove_blacklist:
    input:
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.filtered.bam"
    output:
        blacklist_filtered_bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.blacklist-filtered.bam"
    params:
        blacklist = BLACKLIST
    shell:
        """
        bedtools intersect -nonamecheck -v -abam {input.bam} -b {params.blacklist} > {output.blacklist_filtered_bam}
        samtools sort -O bam -o {output.blacklist_filtered_bam} {output.blacklist_filtered_bam}
        samtools index {output.blacklist_filtered_bam}
        """

# Shift Reads for Tn5 Bias Correction
rule shift_reads:
    input:
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.blacklist-filtered.bam"
    output:
        shifted_bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.shifted.bam"
    params:
        blacklist = BLACKLIST
    shell:
        """
        alignmentSieve --numberOfProcessors max --ATACshift --blackListFileName {params.blacklist} \
            --bam {input.bam} -o {output.shifted_bam}
        samtools index {output.shifted_bam}
        """

# Fragment Size QC & Library Complexity
rule fragment_size_qc:
    input:
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.shifted.bam"
    output:
        frag_size_plot = f"{OUTPUT_DIR}/{{group}}/{{sample}}_fragment_size.png"
    shell:
        """
        plotFingerprint -b {input.bam} --plotFile {output.frag_size_plot} --labels {wildcards.sample} \
            --numberOfSamples 50000
        """


# BAM Processing
#rule Filter_dedup:
#    input:
#        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.marked.bam"
#    output:
#        filtered_bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.filtered.bam"
#    shell:
#        """
#        samtools view -h -q 30 -F 4 -F 256 -F 1024 {input.bam} | samtools view -Sb - > {output.filtered_bam}
#        samtools rmdup {output.filtered_bam} {output.filtered_bam}
#        """


# Peak Calling
rule Call_peaks:
    input:
        filtered_bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.shifted.bam"
    output:
        peaks = f"{OUTPUT_DIR}/{{group}}/{{sample}}_{{caller}}_peaks.narrowPeak"
    params:
        caller = lambda wildcards: wildcards.caller
    shell:
        """
        if [[ "{params.caller}" == "MACS2" ]]; then
            macs2 callpeak -t {input.filtered_bam} -f BAM -g hs --outdir {OUTPUT_DIR}/{wildcards.group} --name {wildcards.sample}_{params.caller};

        elif [[ "{params.caller}" == "Genrich" ]]; then
            Genrich -t {input.filtered_bam} -o {output.peaks} -j -y -v;

        elif [[ "{params.caller}" == "HMMRATAC" ]]; then
            mkdir -p {OUTPUT_DIR}/{wildcards.group}/HMMRATAC
            HMMRATAC -b {input.filtered_bam} -i {GENOME} -o {OUTPUT_DIR}/{wildcards.group}/HMMRATAC/{wildcards.sample}_{params.caller}_peaks.narrowPeak
            cp {OUTPUT_DIR}/{wildcards.group}/HMMRATAC/{wildcards.sample}_{params.caller}_peaks.narrowPeak {output.peaks};

        elif [[ "{params.caller}" == "SEACR" ]]; then
            mkdir -p {OUTPUT_DIR}/{wildcards.group}/SEACR
            bash SEACR.sh {input.filtered_bam} 0.01 {OUTPUT_DIR}/{wildcards.group}/SEACR/{wildcards.sample}_{params.caller}_peaks.narrowPeak
            cp {OUTPUT_DIR}/{wildcards.group}/SEACR/{wildcards.sample}_{params.caller}_peaks.narrowPeak {output.peaks};
        fi
        """

# BigWig Generation (Signal Track)
rule macs2_signal_track:
    input:
        filtered_bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.shifted.bam"
    output:
        bigwig = f"{OUTPUT_DIR}/{{group}}/{{sample}}.bigWig"
    shell:
        "bamCoverage --bam {input.filtered_bam} --outFileName {output.bigwig} --binSize 10 --normalizeUsing RPKM"

# Convert BED to SAF format
rule convert_bed_to_saf:
    input:
        peaks = f"{OUTPUT_DIR}/{{group}}/{{sample}}_{{caller}}_peaks.narrowPeak"
    output:
        saf = f"{OUTPUT_DIR}/{{group}}/{{sample}}_{{caller}}_peaks.saf"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}} \
        {{print $4, $1, $2+1, $3, "."}}' {input.peaks} > {output.saf}
        """

# Count reads in peaks using featureCounts
rule count_reads_in_peaks:
    input:
        saf = f"{OUTPUT_DIR}/{{group}}/{{sample}}_{{caller}}_peaks.saf",
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.shifted.bam"
    output:
        counts = f"{OUTPUT_DIR}/{{group}}/{{sample}}_{{caller}}_readCountInPeaks.txt"
    shell:
        """
        featureCounts -p -a {input.saf} -F SAF -o {output.counts} {input.bam}
        """

# Calculate FRiP score
rule calculate_frip_score:
    input:
        counts = f"{OUTPUT_DIR}/{{group}}/{{sample}}_{{caller}}_readCountInPeaks.txt",
        bam = f"{OUTPUT_DIR}/{{group}}/{{sample}}.shifted.bam"
    output:
        frip_score = f"{OUTPUT_DIR}/{{group}}/{{sample}}_{{caller}}_FRiP.txt"
    shell:
        """
        total_reads=$(samtools view -c {input.bam})
        peak_reads=$(awk 'NR>2 {{sum+=$7}} END {{print sum}}' {input.counts})
        frip=$(echo "scale=4; $peak_reads / $total_reads" | bc)
        echo -e "Total Reads: $total_reads\\nReads in Peaks: $peak_reads\\nFRiP Score: $frip" > {output.frip_score}
        """

# Peak Overlap Analysis
rule compare_peaks:
    input:
        macs2 = f"{OUTPUT_DIR}/{{group}}/{{sample}}_MACS2_peaks.narrowPeak",
        genrich = f"{OUTPUT_DIR}/{{group}}/{{sample}}_Genrich_peaks.narrowPeak",
        hmmratac = f"{OUTPUT_DIR}/{{group}}/{{sample}}_HMMRATAC_peaks.narrowPeak",
        seacr = f"{OUTPUT_DIR}/{{group}}/{{sample}}_SEACR_peaks.narrowPeak"
    output:
        jaccard_results = f"benchmarking/{{group}}/{{sample}}_peak_jaccard.txt"
    shell:
        """
        echo "MACS2 vs Genrich" > {output.jaccard_results}
        bedtools jaccard -a {input.macs2} -b {input.genrich} >> {output.jaccard_results}

        echo "MACS2 vs HMMRATAC" >> {output.jaccard_results}
        bedtools jaccard -a {input.macs2} -b {input.hmmratac} >> {output.jaccard_results}

        echo "MACS2 vs SEACR" >> {output.jaccard_results}
        bedtools jaccard -a {input.macs2} -b {input.seacr} >> {output.jaccard_results}
        """

