#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# USER CONFIGURATION — edit only this section
# ==============================================================================
GENOME="hg38"                 # galgal7, galgal6, hg38, or mm39
THREADS=10
MAPQ_THRESHOLD=0                # minimum MAPQ retained for downstream analysis
                                # I recommend 20 for histone marks, but 0 for TFs
ORIGINAL_FASTQ_PATH="/Data/Ana/NextSeq1000/SCL_11252024_100bp_AASLVB/fastq/"
SAMPLE_ID=("cJun_EpiCypher" "cJun_CST")
R1_PATTERN="_R1_001.fastq.gz"
R2_PATTERN="_R2_001.fastq.gz"

REMOVE_DUPLICATES=true            # option to remove or keep duplicated reads. Default remove, use carefully!!
MACS2_QVALUE=0.01
TRACK_BIN_SIZE=5                  # BigWig resolution in base pairs
KEEP_TRIMMED_FASTQ=true
START_AT="all"                    # all, trim, align, duplicates, tracks, or peaks

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
#   0x200 QC failed, and 0x800 supplementary.
# Duplicate removal additionally excludes 0x400.
PRIMARY_PAIR_EXCLUDE_FLAGS=2828
FINAL_EXCLUDE_FLAGS=2828
if [[ "$REMOVE_DUPLICATES" == true ]]; then
    FINAL_EXCLUDE_FLAGS=3852
fi

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
  *) echo "ERROR: GENOME='$GENOME'; choose galgal7, galgal6, hg38, or mm39." >&2; exit 1 ;;
esac

FASTQ_DIR="fastq"
TRIM_DIR="trimmedFastq"
BAM_DIR="BAM"
BW_DIR="BW"
PEAKS_DIR="Peaks"
STATS_DIR="stats"
filepairs="$FASTQ_DIR/filePairs.txt"

# ==============================================================================
# REQUIRED TOOLS
# ==============================================================================
required_tools=()
(( START_STAGE <= 1 )) && required_tools+=(fastqc)
(( START_STAGE <= 2 )) && required_tools+=(fastqc cutadapt)
(( START_STAGE <= 3 )) && required_tools+=(bowtie2 samtools)
(( START_STAGE <= 4 )) && required_tools+=(picard samtools)
(( START_STAGE <= 5 )) && required_tools+=(bamCoverage)
(( START_STAGE <= 6 )) && required_tools+=(samtools bedtools macs2)
for cmd in "${required_tools[@]}"; do
    command -v "$cmd" &>/dev/null || {
        echo "ERROR: required command '$cmd' not found in PATH" >&2
        exit 1
    }
done

# ==============================================================================
# OUTPUT DIRECTORY SETUP
# ==============================================================================
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$BW_DIR" "$PEAKS_DIR" \
         "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed"
shopt -s nullglob

echo "[Configuration]; MAPQ>=$MAPQ_THRESHOLD; remove_duplicates=$REMOVE_DUPLICATES; q=$MACS2_QVALUE"

# ==============================================================================
# FASTQ DISCOVERY AND PAIRING
# ==============================================================================
if (( START_STAGE <= 3 )); then
    [[ -d "$ORIGINAL_FASTQ_PATH" ]] || {
        echo "ERROR: FASTQ directory does not exist: $ORIGINAL_FASTQ_PATH" >&2
        exit 1
    }

    : > "$filepairs"
    for sample in "${SAMPLE_ID[@]}"; do
        R1_matches=("$ORIGINAL_FASTQ_PATH"/${sample}*"${R1_PATTERN}")

        if (( ${#R1_matches[@]} == 0 )); then
            echo "WARNING: No R1 for ${sample}*${R1_PATTERN}" >&2
            continue
        fi

        for R1path in "${R1_matches[@]}"; do
            R1name=$(basename "$R1path")
            prefix="${R1name%${R1_PATTERN}}"
            R2name="${prefix}${R2_PATTERN}"
            R2path="${ORIGINAL_FASTQ_PATH}/${R2name}"

            if [[ ! -f "$R2path" ]]; then
                echo "WARNING: No R2 for $R1name" >&2
                continue
            fi

            ln -sfn "$R1path" "$FASTQ_DIR/$R1name"
            ln -sfn "$R2path" "$FASTQ_DIR/$R2name"
            echo "${R1name};${R2name}" >> "$filepairs"
        done
    done

    if [[ ! -s "$filepairs" ]]; then
        echo "ERROR: No complete FASTQ pairs were found. Check SAMPLE_ID and read patterns." >&2
        exit 1
    fi

    N_PAIRS=$(wc -l < "$filepairs" | tr -d ' ')
    echo "[Input] ${N_PAIRS} FASTQ pair(s)"
fi

# ==============================================================================
# RAW READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 1 )); then
    echo "[QC:raw]"
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

    echo "[QC:trimmed]"
    fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_trimmed"
fi

# ==============================================================================
# END-TO-END, DOVETAIL-AWARE ALIGNMENT
# ==============================================================================
if (( START_STAGE <= 3 )); then
    trimmed_fastqs=("$TRIM_DIR"/*.fastq.gz)
    (( ${#trimmed_fastqs[@]} > 0 )) || {
        echo "ERROR: No trimmed FASTQ files found for alignment." >&2
        exit 1
    }

    while IFS=';' read -r F1 F2; do
        SAMPLE="${F1%${R1_PATTERN}}"
        OUT_BAM="$BAM_DIR/${SAMPLE}.${GENOME}.bam"

        echo "[Align: CUT&RUN] $SAMPLE (primary proper pairs, MAPQ >= $MAPQ_THRESHOLD)"
        bowtie2 --local --very-sensitive-local --dovetail \
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

        samtools index -@ "$THREADS" "$OUT_BAM"
        samtools flagstat -@ "$THREADS" "$OUT_BAM" > "$STATS_DIR/${SAMPLE}.aligned.flagstat.txt"
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

        if [[ "$REMOVE_DUPLICATES" == true ]]; then
            echo "[Filter] $SAMPLE (marked duplicates removed)"
        else
            echo "[Filter] $SAMPLE (marked duplicates retained)"
        fi

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

echo "[Done]  CUT&RUN pipeline complete"
