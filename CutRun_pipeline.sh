#!/bin/bash
set -euo pipefail
# =============================================== PUT YOUR SETTINGS HERE, THE SCRIPT WILL DO THE REST :) ===============================================
# ===== CONFIG =====
### possible GENOME choices = galgal6, hg38 or mm39
GENOME=galgal6
THREADS=10

# the path ANA will send you
ORIGINAL_FASTQ_PATH="/Data/Ana/NextSeq1000/SCL_08112025_100bp_AAHUSLTK/"

# your samples id as in SampleSheet
SAMPLE_ID=("HH9_AbcamSnai2_rep2"
"HH9_CSTSnai2_rep1"
"HH9_CSTSnai2_rep2"
"HH9_H3K23Ac_rep2"
"HH9_H3K9Ac_rep2"
"HH9_PanklaMono_rep2"
"HH9_UCSCZeb2_rep1"
"HH9_UCSCZeb2_rep2")
# ======================================================================================================================================================



# ==== GENOME SPECS ====
## DO NOT CHANGE THIS!
case "$GENOME" in
  galgal6)
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal6/bowtie2_index/galgal6"
    GENOME_SIZE=1004243013
    MACS2_GENOME="1.0e9"
    ;;
  hg38)
    GENOME_INDEX="/Data/GENOMES/HomoSapiens/hg38.p14/bowtie2_index/hg38.p14"
    GENOME_SIZE=2509528065
    MACS2_GENOME="2.5e9"
    ;;
  mm39)
    GENOME_INDEX="/Data/GENOMES/MusMusculus/mm39/bowtie2_index/mm39"
    GENOME_SIZE=2208238602
    MACS2_GENOME="2.2e9"
    ;;
  *)
    echo "Error: Unknown genome '$GENOME'. Choose from {galgal6, hg38, mm39}."
    exit 1
    ;;
esac

# ===== PATHS =====
FASTQ_DIR="fastq"
TRIM_DIR="trimmedFastq"
BAM_DIR="BAM"
BW_DIR="BW"
PEAKS_DIR="Peaks"
STATS_DIR="stats"

# ===== TOOL CHECK (fail fast if missing) =====
for cmd in fastqc cutadapt bowtie2 samtools picard bamCoverage macs2; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# ===== PREP =====
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$BW_DIR" "$PEAKS_DIR" \
         "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed"

# ensure shell globbing of missing files produces empty list instead of literal
shopt -s nullglob

# ===== SYM LINK =====
for sample in "${SAMPLE_ID[@]}"; do
    matches=("$ORIGINAL_FASTQ_PATH"/${sample}*.fastq.gz)
    if (( ${#matches[@]} == 0 )); then
        echo "Warning: no fastq files found for pattern ${sample}*" >&2
        continue
    fi
    for fq in "${matches[@]}"; do
        ln -sfn "$fq" "$FASTQ_DIR/$(basename "$fq")"
    done
done

# ===== Make list of pairs (R1;R2) =====
filepairs="$FASTQ_DIR/filePairs.txt"
: > "$filepairs"
for R1path in "$FASTQ_DIR"/*_R1_*.fastq.gz; do
    [[ -e "$R1path" ]] || continue
    R2path="${R1path/_R1_/_R2_}"
    if [[ -f "$R2path" ]]; then
        echo "$(basename "$R1path");$(basename "$R2path")" >> "$filepairs"
    else
        echo "Warning: no R2 found for $(basename "$R1path")" >&2
    fi
done

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
        "$FASTQ_DIR/$F1" "$FASTQ_DIR/$F2" \
    > "$STATS_DIR/${F1%.fastq.gz}_cutadapt_report.txt" 2>&1
done < "$filepairs"

# ===== FASTQC after trimming =====
echo "[FastQC] Executing in trimmed fastq..."
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_trimmed"

# ===== ALIGNMENT =====
while IFS=";" read -r F1 F2; do
    BASE=$(basename "$F1" .fastq.gz)
    if echo "$BASE" | grep -qP '_S[0-9]+'; then
        SNUM=$(echo "$BASE" | grep -oP '_S[0-9]+' | sed 's/^_S//')
    else
        SNUM="NA"
    fi
    SAMPLE=$(echo "$BASE" | sed -E "s/_S[0-9]+.*//")
    OUT_BAM="$BAM_DIR/${SAMPLE}.${GENOME}.bam"

    echo "[Bowtie2] $SAMPLE ($F1 / $F2)"
    set -o pipefail
    bowtie2 --local --very-sensitive-local \
        --no-unal --no-mixed --no-discordant \
        -x "$GENOME_INDEX" -I 10 -X 1000 \
        -1 "$TRIM_DIR/trimmed_${F1}" -2 "$TRIM_DIR/trimmed_${F2}" \
        --threads "$THREADS" \
    | samtools view -@ "$THREADS" -F 780 -f 2 -bh - \
    | samtools sort -@ "$THREADS" -T "${OUT_BAM%.bam}" \
    | samtools addreplacerg \
        -r "@RG\tID:${SAMPLE}.${GENOME}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPU:S${SNUM}" \
        -o "$OUT_BAM" -
    set +o pipefail
done < "$filepairs"

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
    echo "[Samtools remove duplicates] $BAMFILE"
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
    [[ -e "$BAMFILE" ]] || continue
    SAMPLE=$(basename "$BAMFILE" .bam)
    echo "[MACS2] $SAMPLE"
    macs2 callpeak -t "$BAMFILE" \
        -n "$SAMPLE" \
        -f BAMPE -g "$MACS2_GENOME" \
        -q 0.05 \
        --call-summits \
        --outdir "$PEAKS_DIR"
done

# ===== CLEAN INTERMEDIATES =====
echo "[Cleanup] Removing intermediate files..."
rm -rf "$TRIM_DIR"

# remove intermediate bams and bai that are not nodups or dupsMarked
for f in "$BAM_DIR"/*; do
    case "$(basename "$f")" in
        *_nodups.bam|*_nodups.bam.bai|*_dupsMarked.bam|*_dupsMarked.bai) continue ;;
        *.bam|*.bai) rm -f "$f" ;;
    esac
done

echo "Pipeline complete!"
