#!/bin/bash
set -euo pipefail

# ===== CONFIG =====
GENOME_INDEX="/Data/Austin/workdir/genome/hg38/hg38_bt2/GRCh38"
GENOME_SIZE=2913022398
THREADS=10
ORIGINAL_FASTQ_PATH="/Data/Ana/Jackie/ATACAug2025/fastq/"
SAMPLE_ID="D5_ATAC"

# ===== PATHS =====
FASTQ_DIR="fastq"
TRIM_DIR="trimmedFastq"
BAM_DIR="BAM"
BW_DIR="BW"
PEAKS_DIR="Peaks"
STATS_DIR="stats"

# ===== PREP =====
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$BW_DIR" "$PEAKS_DIR" \
         "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed"

# ===== SYM LINK =====
for file in $(ls $ORIGINAL_FASTQ_PATH/${SAMPLE_ID}*.fastq.gz);
	do ln -s $file $FASTQ_DIR/
done

# Make list of pairs
ls --color=never "$FASTQ_DIR"/*fastq.gz | sort | sed "s|$FASTQ_DIR/||g" | paste -sd';\n' - > "$FASTQ_DIR/filePairs.txt"

# ===== FASTQC (raw data) =====
echo "[FastQC] Executing in raw data..."
fastqc -t "$THREADS" "$FASTQ_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_raw"

# ===== ADAPTER TRIMMING =====
while IFS=";" read -r F1 F2; do
    echo "[Cutadapt] $F1 / $F2"
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --minimum-length=25 -j $THREADS \
        -o "$TRIM_DIR/trimmed_${F1}" -p "$TRIM_DIR/trimmed_${F2}" \
        "$FASTQ_DIR/$F1" "$FASTQ_DIR/$F2" > "$STATS_DIR/${F1}_cutadapt_report.txt"
done < "$FASTQ_DIR/filePairs.txt"

# ===== FASTQC after trimming =====
echo "[FastQC] Executing in trimmed fastq..."
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq* -o "$STATS_DIR/fastqc_trimmed"

# ===== ALIGNMENT =====
while IFS=";" read -r F1 F2; do
    OUT_BAM="$BAM_DIR/${F1/.fastq.gz/.bam}"
    echo "[Bowtie2] $F1 / $F2"
    bowtie2 --local --very-sensitive-local \
        --no-unal --no-mixed --no-discordant \
        -x "$GENOME_INDEX" -I 10 -X 1000 \
        -1 "$TRIM_DIR/trimmed_${F1}" -2 "$TRIM_DIR/trimmed_${F2}" \
        --threads "$THREADS" \
    | samtools view -@ "$THREADS" -F 780 -f 2 -bh - \
    | samtools sort -@ "$THREADS" -T "${OUT_BAM%.bam}" \
    | samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o "$OUT_BAM" -
done < "$FASTQ_DIR/filePairs.txt"

# ===== MARK DUPLICATES =====
for BAMFILE in "$BAM_DIR"/*.bam; do
    echo "[Picard] $BAMFILE"
    picard MarkDuplicates \
        I="$BAMFILE" \
        O="${BAMFILE/.bam/_dupsMarked.bam}" \
        M="${BAMFILE/.bam/_dupsMarkedStats.txt}" \
        VALIDATION_STRINGENCY=SILENT
done

# ===== REMOVE DUPLICATES =====
for BAMFILE in "$BAM_DIR"/*_dupsMarked.bam; do
    echo "[Samtools rmdup] $BAMFILE"
    samtools view -@ "$THREADS" -F 1804 -f 2 -b "$BAMFILE" \
	    | samtools sort -@ "$THREADS" -o "${BAMFILE/_dupsMarked.bam/_nodups.bam}"
    samtools index "${BAMFILE/_dupsMarked.bam/_nodups.bam}"
done

# ===== BIGWIG =====
for BAMFILE in "$BAM_DIR"/*_nodups.bam; do
    bamCoverage --bam "$BAMFILE" \
        --outFileName "$BW_DIR/$(basename "$BAMFILE").bw" \
        --outFileFormat bigwig \
        --binSize 5 \
        --numberOfProcessors "$THREADS" \
        --normalizeUsing RPKM \
        --extendReads
done

# ===== MACS2 =====
for BAMFILE in "$BAM_DIR"/*_nodups.bam; do
    SAMPLE=$(basename "$BAMFILE" .bam)
    macs2 callpeak -t "$BAMFILE" \
        -n "$SAMPLE" \
        -f BAMPE -g "$GENOME_SIZE" \
        -q 0.05 \
        --call-summits \
        --outdir "$PEAKS_DIR"
done

# ===== CLEAN INTERMEDIATES =====
echo "[Cleanup] Removing intermediate files..."
rm -f "$TRIM_DIR"/*.fastq.gz
for file in $(ls --color=never "$BAM_DIR"/*.bam "$BAM_DIR"/*.bai | grep -E -v "nodups|dupsMarked" ); do
      rm -rf $file
done

echo "Pipeline complete!"
