#!/bin/bash
set -e
set -u
set -o pipefail

while read -r sample <&3
do

echo "processing sample: $sample"

# Detect sequencing layout
layout=$(esearch -db sra -query "$sample" | efetch -format runinfo | csvcut -c LibraryLayout | tail -n 1)

if [ -z "$layout" ]; then
    echo "ERROR: Could not determine sequencing layout"
    exit 1
elif [ "$layout" = "PAIRED" ]; then
    echo "This sample is paired-end"
elif [ "$layout" = "SINGLE" ]; then
    echo "This sample is single-end"
else
    echo "ERROR: Unexpected sequencing layout: $layout"
    exit 1
fi

# Download SRA data
prefetch "$sample"

if [ "$layout" = "PAIRED" ]; then

    if [ ! -f "data/raw/${sample}_1.fastq" ] || [ ! -f "data/raw/${sample}_2.fastq" ]; then
        echo "Creating paired-end FASTQ files..."
        fasterq-dump "$sample" --split-files --outdir data/raw --threads 8
    else
        echo "Paired-end FASTQ files already exist. Skipping conversion."
    fi

elif [ "$layout" = "SINGLE" ]; then

    if [ ! -f "data/raw/${sample}.fastq" ]; then
        echo "Creating single-end FASTQ file..."
        fasterq-dump "$sample" --outdir data/raw --threads 8
    else
        echo "Single-end FASTQ file already exists. Skipping conversion."
    fi

fi

# Raw-read quality control with Falco
if [ "$layout" = "PAIRED" ]; then
    echo "Running Falco on paired-end raw reads..."

    falco "data/raw/${sample}_1.fastq" --no-tiles -o qc/raw
    falco "data/raw/${sample}_2.fastq" --no-tiles -o qc/raw

elif [ "$layout" = "SINGLE" ]; then
    echo "Running Falco on single-end raw reads..."

    falco "data/raw/${sample}.fastq" -o qc/raw

fi

# Read trimming with Trimmomatic
if [ "$layout" = "PAIRED" ]; then

    echo "Running Trimmomatic on paired-end reads..."

    trimmomatic PE -threads 8 \
        "data/raw/${sample}_1.fastq" \
        "data/raw/${sample}_2.fastq" \
        "data/trimmed/${sample}_1.paired.fastq" \
        "data/trimmed/${sample}_1.unpaired.fastq" \
        "data/trimmed/${sample}_2.paired.fastq" \
        "data/trimmed/${sample}_2.unpaired.fastq" \
        SLIDINGWINDOW:4:30

elif [ "$layout" = "SINGLE" ]; then

    echo "Running Trimmomatic on single-end reads..."

    trimmomatic SE -threads 8 \
        "data/raw/${sample}.fastq" \
        "data/trimmed/${sample}.trimmed.fastq" \
        SLIDINGWINDOW:4:30

fi

# Quality control after trimming
if [ "$layout" = "PAIRED" ]; then
    echo "Running Falco on trimmed paired-end reads..."

    falco "data/trimmed/${sample}_1.paired.fastq" --no-tiles -o qc/trimmed
    falco "data/trimmed/${sample}_2.paired.fastq" --no-tiles -o qc/trimmed

elif [ "$layout" = "SINGLE" ]; then
    echo "Running Falco on trimmed single-end reads..."

    falco "data/trimmed/${sample}.trimmed.fastq" -o qc/trimmed

fi

# Host alignment with HISAT2
if [ -f "results/alignment/${sample}.hisat2.log" ] && \
   [ -f "results/alignment/${sample}.sam" ]; then

    echo "HISAT2 outputs already exist. Skipping alignment."

else
    if [ "$layout" = "PAIRED" ]; then
        echo "Running HISAT2 on paired-end reads..."

        hisat2 -p 8 \
            -x index/grch38/genome \
            -1 "data/trimmed/${sample}_1.paired.fastq" \
            -2 "data/trimmed/${sample}_2.paired.fastq" \
            --un-conc "results/unmapped/${sample}.unmapped.%.fastq" \
            -S "results/alignment/${sample}.sam" \
            2> "results/alignment/${sample}.hisat2.log"

    elif [ "$layout" = "SINGLE" ]; then
        echo "Running HISAT2 on single-end reads..."

        hisat2 -p 8 \
            -x index/grch38/genome \
            -U "data/trimmed/${sample}.trimmed.fastq" \
            --un "results/unmapped/${sample}.unmapped.fastq" \
            -S "results/alignment/${sample}.sam" \
            2> "results/alignment/${sample}.hisat2.log"
    fi
fi

samtools view -bS "results/alignment/${sample}.sam" > "results/alignment/${sample}.bam"

samtools sort "results/alignment/${sample}.bam" -o "results/alignment/${sample}.sorted.bam"

samtools index "results/alignment/${sample}.sorted.bam"

samtools flagstat "results/alignment/${sample}.sorted.bam" > "results/alignment/${sample}.flagstat.txt"

echo "completed alignment successfully"

# Quality control of host-unmapped reads
if [ "$layout" = "PAIRED" ]; then

    echo "Running Falco on paired-end host-unmapped reads..."

    falco "results/unmapped/${sample}.unmapped.1.fastq" --no-tiles -o qc/unmapped
    falco "results/unmapped/${sample}.unmapped.2.fastq" --no-tiles -o qc/unmapped

elif [ "$layout" = "SINGLE" ]; then

    echo "Running Falco on single-end host-unmapped reads..."

    falco "results/unmapped/${sample}.unmapped.fastq" -o qc/unmapped

fi

# Kraken2 classification
echo "Running Kraken2 on host-unmapped reads..."

if [ "$layout" = "PAIRED" ]; then

    echo "Running Kraken2 separately on forward and reverse host-unmapped reads..."

    # Forward read
    kraken2 \
        --db database/k2_standard_20240904 \
        --threads 8 \
        --memory-mapping \
        --confidence 0.00 \
        --minimum-hit-groups 2 \
        --report "results/kraken2/${sample}_1.report" \
        --output /dev/null \
        "results/unmapped/${sample}.unmapped.1.fastq"

    # Reverse read
    kraken2 \
        --db database/k2_standard_20240904 \
        --threads 8 \
        --memory-mapping \
        --confidence 0.00 \
        --minimum-hit-groups 2 \
        --report "results/kraken2/${sample}_2.report" \
        --output /dev/null \
        "results/unmapped/${sample}.unmapped.2.fastq"

elif [ "$layout" = "SINGLE" ]; then

    echo "Running Kraken2 on single-end host-unmapped reads..."

    kraken2 \
        --db database/k2_standard_20240904 \
        --threads 8 \
        --memory-mapping \
        --confidence 0.00 \
        --minimum-hit-groups 2 \
        --report "results/kraken2/${sample}.report" \
        --output /dev/null \
        "results/unmapped/${sample}.unmapped.fastq"

fi

echo "Kraken2 completed."

# Bracken species-level abundance estimation
echo "Running Bracken species-level abundance estimation..."

if [ "$layout" = "PAIRED" ]; then

    echo "Running Bracken separately on forward and reverse Kraken2 reports..."

    # Forward read
    bracken \
        -d database/k2_standard_20240904 \
        -i "results/kraken2/${sample}_1.report" \
        -o "results/bracken/${sample}_1.bracken" \
        -w "results/bracken/${sample}_1.bracken.kraken" \
        -r 100 \
        -l S


    # Reverse read
    bracken \
        -d database/k2_standard_20240904 \
        -i "results/kraken2/${sample}_2.report" \
        -o "results/bracken/${sample}_2.bracken" \
        -w "results/bracken/${sample}_2.bracken.kraken" \
        -r 100 \
        -l S

elif [ "$layout" = "SINGLE" ]; then

    echo "Running Bracken on single-end Kraken2 report..."

    bracken \
        -d database/k2_standard_20240904 \
        -i "results/kraken2/${sample}.report" \
        -o "results/bracken/${sample}.bracken" \
        -w "results/bracken/${sample}.bracken.kraken" \
        -r 100 \
        -l S

fi

echo "Bracken completed."

# Clean up large intermediate files
echo "Cleaning up large intermediate files..."

if [ "$layout" = "PAIRED" ]; then

    rm -f "data/raw/${sample}_1.fastq"
    rm -f "data/raw/${sample}_2.fastq"

    rm -f "data/trimmed/${sample}_1.paired.fastq"
    rm -f "data/trimmed/${sample}_1.unpaired.fastq"
    rm -f "data/trimmed/${sample}_2.paired.fastq"
    rm -f "data/trimmed/${sample}_2.unpaired.fastq"

elif [ "$layout" = "SINGLE" ]; then

    rm -f "data/raw/${sample}.fastq"
    rm -f "data/trimmed/${sample}.trimmed.fastq"

fi

rm -f "results/alignment/${sample}.sam"
rm -f "results/alignment/${sample}.bam"
rm -f "$HOME/ncbi/sra/${sample}.sra"

echo "Cleanup completed."

echo "completed processing: $sample safely and clean"

done 3< samples.txt

echo "======================================"
echo "ALL SAMPLES FINISHED"
echo "======================================"

echo "Running final MultiQC report..."

multiqc qc/raw qc/trimmed qc/unmapped results/alignment -o qc/multiqc

echo "MultiQC completed."

echo "Removing sorted BAM files to save disk space..."

rm -f results/alignment/*.sorted.bam
rm -f results/alignment/*.sorted.bam.bai

echo "======================================"
echo "ALL SAMPLES COMPLETED SUCCESSFULLY"
echo "======================================"
