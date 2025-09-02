#!/bin/bash
set -euo pipefail
source /etc/profile.d/apps.sh

### Effective genome size for Human: 2.7e9
### Effective genome size for Mouse: 1.87e9
### Effective genome size for Chicken: 1.2e9
# ===== CONFIG =====
GENOME_INDEX="/data/Megan/genome_data/galGal6/galGal6_bt2/galGal6_bt2"
GENOME_SIZE="1.2e9"
THREADS=10
# the path ANA will send you
ORIGINAL_FASTQ_PATH="/Data/Fjodor/MRE_PROJECT/osteoblast_atac_seq/functional_exp/fastq_og"
# your samples id as in SampleSheet
SAMPLE_ID=("day3_osteo_control_rep1" "day3_osteo_control_rep2" "day3_osteo_rotenone_rep1" "day3_osteo_rotenone_rep2" 
		"day6_osteo_control_rep1" "day6_osteo_control_rep2" "day6_osteo_rotenone_rep1" "day6_osteo_rotenone_rep2" 
		"pgo91_control_rep1" "pgo91_control_rep2")

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
for sample in "${SAMPLE_ID[@]}"; do
	ln -s $ORIGINAL_FASTQ_PATH/${sample}*.fastq.gz $FASTQ_DIR
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
        --minimum-length=25 -j "$THREADS" \
        -o "$TRIM_DIR/trimmed_${F1}" -p "$TRIM_DIR/trimmed_${F2}" \
        "$FASTQ_DIR/$F1" "$FASTQ_DIR/$F2" > "$STATS_DIR/${F1}_cutadapt_report.txt"
done < "$FASTQ_DIR/filePairs.txt"

# ===== FASTQC after trimming =====
echo "[FastQC] Executing in trimmed fastq..."
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_trimmed"

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
    PICARD MarkDuplicates \
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
        --normalizeUsing RPGC \
	--effectiveGenomeSize "$GENOME_SIZE" \
        --extendReads
done

# ===== MACS2 =====
for BAMFILE in "$BAM_DIR"/*_nodups.bam; do
    SAMPLE=$(basename "$BAMFILE" .bam)
    macs2 callpeak -t "$BAMFILE" \
        -n "$SAMPLE" \
        -f BAM -g "$GENOME_SIZE" \
        -q 0.05 \
        --call-summits \
        --outdir "$PEAKS_DIR"
done

# ===== CLEAN INTERMEDIATES =====
echo "[Cleanup] Removing intermediate files..."
rm -f "$TRIM_DIR"/*.fastq.gz
for file in $(ls --color=never "$BAM_DIR"/*.bam "$BAM_DIR"/*.bai | grep -E -v "nodups|dupsMarked"); do
      rm -rf $file
done

echo "Pipeline complete!"

