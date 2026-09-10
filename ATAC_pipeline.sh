#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# USER CONFIGURATION — edit only this section
# ==============================================================================
GENOME="mm39"                    # galgal7, galgal6, hg38, or mm39
THREADS=10                       # CPU threads
MAPQ_THRESHOLD=30                # Minimum MAPQ retained for downstream analysis
ORIGINAL_FASTQ_PATH="/Data/Fjodor/MRE_PROJECT/osteoblast_atac_seq/functional_exp/fastq_og"
SAMPLE_ID=("day3_osteo_control_rep1" "day3_osteo_control_rep2" "day3_osteo_rotenone_rep1" "day3_osteo_rotenone_rep2" 
		"day6_osteo_control_rep1" "day6_osteo_control_rep2" "day6_osteo_rotenone_rep1" "day6_osteo_rotenone_rep2" 
		"pgo91_control_rep1" "pgo91_control_rep2")
R1_PATTERN="_R1_001.fastq.gz"
R2_PATTERN="_R2_001.fastq.gz"
TRACK_BIN_SIZE=5                # BigWig resolution in base pairs
MACS2_QVALUE=0.05
KEEP_TRIMMED_FASTQ=false         # true keeps trimmed reads; false removes them
START_AT="all"                   # all, trim, align, duplicates, tracks, or peaks

# ==============================================================================
# END USER CONFIGURATION
# ==============================================================================

case "$START_AT" in
  all)        START_STAGE=1 ;;
  trim)       START_STAGE=2 ;;
  align)      START_STAGE=3 ;;
  duplicates) START_STAGE=4 ;;
  tracks)     START_STAGE=5 ;;
  peaks)      START_STAGE=6 ;;
  *) echo "ERROR: START_AT must be all, trim, align, duplicates, tracks, or peaks." >&2; exit 1 ;;
esac

# SAM flags excluded from aligned BAMs:
#   0x4 unmapped, 0x8 mate unmapped, 0x100 secondary,
#   0x200 QC failed, 0x800 supplementary; final filtering also excludes 0x400 duplicates.
PRIMARY_PAIR_EXCLUDE_FLAGS=2828
FINAL_EXCLUDE_FLAGS=3852


# ==============================================================================
# GENOME SETTINGS — do not edit unless reference paths change
# ==============================================================================
case "$GENOME" in
  galgal7)
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal7/bowtie2_index_main/galgal7_mainchr"
    GENOME_SIZE=1041139641
    MACS2_GENOME="1.0e9"
    ;;
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
    echo "ERROR: GENOME='$GENOME'; choose galgal7, galgal6, hg38, or mm39." >&2
    exit 1
    ;;
esac

# ==============================================================================
# OUTPUT PATHS
# ==============================================================================
FASTQ_DIR="fastq"
TRIM_DIR="trimmedFastq"
BAM_DIR="BAM"
BW_DIR="BW"
PEAKS_DIR="Peaks"
STATS_DIR="stats"

# ==============================================================================
# REQUIRED TOOLS
# ==============================================================================
for cmd in fastqc cutadapt bowtie2 samtools picard bamCoverage macs2; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

[[ -d "$ORIGINAL_FASTQ_PATH" ]] || {
    echo "ERROR: FASTQ directory does not exist: $ORIGINAL_FASTQ_PATH" >&2
    exit 1
}

# ==============================================================================
# OUTPUT DIRECTORY SETUP
# ==============================================================================
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$BW_DIR" "$PEAKS_DIR" \
         "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed"

# ensure shell globbing of missing files produces empty list instead of literal
shopt -s nullglob

# ==============================================================================
# FASTQ DISCOVERY AND PAIRING
# ==============================================================================
filepairs="$FASTQ_DIR/filePairs.txt"
: > "$filepairs"

for sample in "${SAMPLE_ID[@]}"; do
    # Find all R1 files belonging to this sample
    R1_matches=("$ORIGINAL_FASTQ_PATH"/${sample}*"${R1_PATTERN}")

    if (( ${#R1_matches[@]} == 0 )); then
        echo "WARNING: No R1 for ${sample}*${R1_PATTERN}" >&2
        continue
    fi

    for R1path in "${R1_matches[@]}"; do

        R1name=$(basename "$R1path")

        # Remove R1_PATTERN and replace it with R2_PATTERN
        prefix="${R1name%${R1_PATTERN}}"
        R2name="${prefix}${R2_PATTERN}"
        R2path="${ORIGINAL_FASTQ_PATH}/${R2name}"

        # Check that matching R2 exists
        if [[ ! -f "$R2path" ]]; then
            echo "WARNING: No R2 for $R1name" >&2
            continue
        fi

        # Create symlinks
        ln -sfn "$R1path" "$FASTQ_DIR/$R1name"
        ln -sfn "$R2path" "$FASTQ_DIR/$R2name"

        # Add pair to filePairs.txt
        echo "${R1name};${R2name}" >> "$filepairs"

    done
done

if [[ ! -s "$filepairs" ]]; then
    echo "ERROR: No complete FASTQ pairs were found. Check SAMPLE_ID and read patterns." >&2
    exit 1
fi

# ==============================================================================
# FASTQ VALIDATION
# ==============================================================================
N_R1=$(find "$FASTQ_DIR" -maxdepth 1 -name "*${R1_PATTERN}" | wc -l)
N_R2=$(find "$FASTQ_DIR" -maxdepth 1 -name "*${R2_PATTERN}" | wc -l)

echo "[Input] R1=${N_R1}, R2=${N_R2}"

if [[ "$N_R1" -eq 0 ]]; then
    echo "ERROR: No R1 FASTQ files matched the specified pattern."
    exit 1
fi

if [[ "$N_R1" -ne "$N_R2" ]]; then
    echo "WARNING: Number of R1 and R2 files differs."
fi

# ==============================================================================
# RAW READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 1 )); then
echo "[QC:raw] ${N_R1} pair(s)"
fastqc -t "$THREADS" "$FASTQ_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_raw"
fi

# ==============================================================================
# ADAPTER TRIMMING
# ==============================================================================
if (( START_STAGE <= 2 )); then
while IFS=";" read -r F1 F2; do
    echo "[Trim] $F1 + $F2"
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --minimum-length=25 -j "$THREADS" \
        -o "$TRIM_DIR/trimmed_${F1}" -p "$TRIM_DIR/trimmed_${F2}" \
        "$FASTQ_DIR/$F1" "$FASTQ_DIR/$F2" \
    > "$STATS_DIR/${F1%.fastq.gz}_cutadapt_report.txt" 2>&1
done < "$filepairs"
fi

# ==============================================================================
# TRIMMED READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 2 )); then
echo "[QC:trimmed]"
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_trimmed"
fi

# ==============================================================================
# ALIGNMENT
# ==============================================================================
if (( START_STAGE <= 3 )); then
trimmed_fastqs=("$TRIM_DIR"/*.fastq.gz)
(( ${#trimmed_fastqs[@]} > 0 )) || { echo "ERROR: No trimmed FASTQ files found for alignment." >&2; exit 1; }
while IFS=";" read -r F1 F2; do
    SAMPLE="${F1%${R1_PATTERN}}"
    OUT_BAM="$BAM_DIR/${SAMPLE}.${GENOME}.bam"

    echo "[Align:Bowtie2] $SAMPLE (primary proper pairs, MAPQ >= $MAPQ_THRESHOLD)"
    set -o pipefail
    bowtie2 --local --very-sensitive-local \
        --no-unal --no-mixed --no-discordant --phred33 \
        -x "$GENOME_INDEX" -I 10 -X 1000 \
        -1 "$TRIM_DIR/trimmed_${F1}" -2 "$TRIM_DIR/trimmed_${F2}" \
        --threads "$THREADS" \
    | samtools view -@ "$THREADS" -q "$MAPQ_THRESHOLD" \
        -F "$PRIMARY_PAIR_EXCLUDE_FLAGS" -f 2 -bh - \
    | samtools sort -@ "$THREADS" -T "${OUT_BAM%.bam}" \
    | samtools addreplacerg \
        -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
        -o "$OUT_BAM" -
    set +o pipefail
done < "$filepairs"
fi

# ==============================================================================
# DUPLICATE METRICS AND FINAL ALIGNMENT FILTERING
# ==============================================================================
if (( START_STAGE <= 4 )); then
    raw_bams=("$BAM_DIR"/*."$GENOME".bam)
    (( ${#raw_bams[@]} > 0 )) || {
        echo "ERROR: No aligned BAM files found." >&2
        exit 1
    }

    for BAMFILE in "${raw_bams[@]}"; do
        SAMPLE=$(basename "$BAMFILE" ".${GENOME}.bam")
        MARKED="$BAM_DIR/${SAMPLE}.${GENOME}_dupsMarked.bam"
        FILTERED="$BAM_DIR/${SAMPLE}.${GENOME}_nodups.bam"

        echo "[Duplicates:mark] $SAMPLE"
        picard MarkDuplicates \
            I="$BAMFILE" \
            O="$MARKED" \
            M="$BAM_DIR/${SAMPLE}.${GENOME}_dupsMarkedStats.txt" \
            VALIDATION_STRINGENCY=SILENT

        echo "[Filter] $SAMPLE (marked duplicates removed)"

        samtools view -@ "$THREADS" -q "$MAPQ_THRESHOLD" \
          -F "$FINAL_EXCLUDE_FLAGS" -f 2 -b "$MARKED" \
        | samtools sort -@ "$THREADS" -o "$FILTERED"
        samtools index -@ "$THREADS" "$FILTERED"
        samtools flagstat -@ "$THREADS" "$FILTERED" > "$STATS_DIR/${SAMPLE}.filtered.flagstat.txt"
    done
fi

# ==============================================================================
# BIGWIG TRACK GENERATION
# ==============================================================================
if (( START_STAGE <= 5 )); then
nodup_bams=("$BAM_DIR"/*_nodups.bam)
(( ${#nodup_bams[@]} > 0 )) || { echo "ERROR: No duplicate-filtered BAM files found for track generation." >&2; exit 1; }
for BAMFILE in "${nodup_bams[@]}"; do
    TRACK_NAME=$(basename "$BAMFILE" .bam)
    echo "[Track:CPM] $TRACK_NAME"
    bamCoverage --bam "$BAMFILE" \
        --outFileName "$BW_DIR/${TRACK_NAME}.CPM.bw" \
        --outFileFormat bigwig \
        --binSize "$TRACK_BIN_SIZE" \
        --numberOfProcessors "$THREADS" \
        --normalizeUsing CPM \
        --extendReads
    echo "[Track:RPGC] $TRACK_NAME"
    bamCoverage --bam "$BAMFILE" \
        --outFileName "$BW_DIR/${TRACK_NAME}.RPGC.bw" \
        --outFileFormat bigwig \
        --binSize "$TRACK_BIN_SIZE" \
        --numberOfProcessors "$THREADS" \
        --normalizeUsing RPGC \
        --effectiveGenomeSize "$GENOME_SIZE" \
        --extendReads
done
fi

# ==============================================================================
# PEAK CALLING
# ==============================================================================
if (( START_STAGE <= 6 )); then
nodup_bams=("$BAM_DIR"/*_nodups.bam)
(( ${#nodup_bams[@]} > 0 )) || { echo "ERROR: No duplicate-filtered BAM files found for peak calling." >&2; exit 1; }
for BAMFILE in "${nodup_bams[@]}"; do
    [[ -e "$BAMFILE" ]] || continue
    SAMPLE=$(basename "$BAMFILE" .bam)
    echo "[Peaks:MACS2] $SAMPLE"
    macs2 callpeak -t "$BAMFILE" \
        -n "$SAMPLE" \
        -f BAMPE -g "$MACS2_GENOME" \
        -q "$MACS2_QVALUE" \
        --call-summits \
        --keep-dup all \
        --outdir "$PEAKS_DIR"
done
fi

# ==============================================================================
# CLEANUP
# ==============================================================================
echo "[Cleanup]"
if [[ "$KEEP_TRIMMED_FASTQ" == false ]] && (( START_STAGE <= 3 )); then
    rm -rf "$TRIM_DIR"
fi

# remove intermediate bams and bai that are not nodups or dupsMarked
for f in "$BAM_DIR"/*; do
    case "$(basename "$f")" in
        *_nodups.bam|*_nodups.bam.bai|*_dupsMarked.bam|*_dupsMarked.bai) continue ;;
        *.bam|*.bai) rm -f "$f" ;;
    esac
done

echo "[Done] Pipeline complete"
