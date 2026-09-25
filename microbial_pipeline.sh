#!/bin/bash
set -e
set -u
set -o pipefail

sample=$1

echo "processing sample: $sample"

if [ ! -f "data/raw/${sample}.fastq" ]; then

    echo "Raw FASTQ not found. Downloading and converting $sample..."

    prefetch "$sample"

    fasterq-dump "$sample" --outdir data/raw --threads 8

    echo "Raw FASTQ created."

else

    echo "Raw FASTQ already exists. Skipping download."

fi

fastqc "data/raw/${sample}.fastq" -o qc

if [ ! -f "data/trimmed/${sample}.trimmed.fastq" ]; then
    echo "Running fastp..."

    fastp \
    -i "data/raw/${sample}.fastq" \
    -o "data/trimmed/${sample}.trimmed.fastq" \
    -h qc/fastp.html \
    -j qc/fastp.json

    echo "fastp completed."
else
    echo "Trimmed FASTQ already exists. Skipping fastp."
fi

fastqc "data/trimmed/${sample}.trimmed.fastq" -o  qc

if [ ! -f "results/unmapped/${sample}.unmapped.fastq" ]; then
    echo "Running HISAT2 host alignment..."

    hisat2 -p 8 \
        -x index/grch38/genome \
        -U "data/trimmed/${sample}.trimmed.fastq" \
        --un "results/unmapped/${sample}.unmapped.fastq" \
        -S "results/alignment/${sample}.sam" \
        2> "results/alignment/${sample}.hisat2.log"

    echo "HISAT2 completed."
else
    echo "Unmapped FASTQ already exists. Skipping HISAT2."
fi

samtools view -bS "results/alignment/${sample}.sam" > "results/alignment/${sample}.bam"

samtools sort "results/alignment/${sample}.bam" -o "results/alignment/${sample}.sorted.bam"

samtools index "results/alignment/${sample}.sorted.bam"

samtools flagstat "results/alignment/${sample}.sorted.bam" > "results/alignment/${sample}.flagstat.txt"

echo "completed alignment successfully"

echo "Running Kraken2 on unmapped reads..."

kraken2 --db database/k2_standard_20240904 --threads 2 --memory-mapping --confidence 0.00  --minimum-hit-groups 2 --report "results/kraken2/${sample}.report" --output "results/kraken2/${sample}.kraken" "results/unmapped/${sample}.unmapped.fastq"

echo "Kraken2 completed."

echo "Running Bracken species-level abundance estimation..."

bracken -d database/k2_standard_20240904 -i "results/kraken2/${sample}.report" -o "results/bracken/${sample}.bracken" -w "results/bracken/${sample}.bracken.kraken" -r 100 -l S -t 10

echo "Bracken completed."

echo "Running MultiQC..."

multiqc qc results/alignment -o qc/multiqc

echo "MultiQC report completed."

echo "Cleaning up large intermediate files..."

rm -f "data/raw/${sample}.fastq"
rm -f "data/trimmed/${sample}.trimmed.fastq"
rm -f "results/alignment/${sample}.sam"
rm -f "results/alignment/${sample}.bam"
rm -f "results/alignment/${sample}.sorted.bam"
rm -f "results/alignment/${sample}.sorted.bam.bai"
rm -f "$HOME/ncbi/sra/${sample}.sra"

echo "Cleanup completed."

echo "completed processing: $sample safely and clean"
