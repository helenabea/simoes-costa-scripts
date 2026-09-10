#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# USER CONFIGURATION — edit only this section
# ==============================================================================
GENOME="galgal7"                 # galgal7, galgal6, hg38, or mm39
THREADS=10                       # CPU threads
ORIGINAL_FASTQ_PATH="/Data/Ana/NextSeq1000/SCL_10102025_100bp_AA/fastq"
SAMPLE_ID=("CRISPR_3RNA")
R1_PATTERN="_R1_001.fastq.gz"    # Exact filename suffix for single-end reads
MIN_LEN=20                       # Minimum read length after trimming
FEATURECOUNTS_STRANDED=1         # 0=unstranded, 1=forward, 2=reverse
FEATURECOUNTS_MODE="matrix"      # matrix or per_sample
TRACK_BIN_SIZE=10               # BigWig resolution in base pairs
KEEP_TRIMMED_FASTQ=false         # true keeps trimmed reads; false removes them
START_AT="all"                   # all, trim, align, or count

# ==============================================================================
# END USER CONFIGURATION
# ==============================================================================

case "$START_AT" in
  all)       START_STAGE=1 ;;
  trim)      START_STAGE=2 ;;
  align)     START_STAGE=3 ;;
  count)     START_STAGE=4 ;;
  *) echo "ERROR: START_AT must be all, trim, align, or count." >&2; exit 1 ;;
esac


# ==============================================================================
# GENOME SETTINGS — do not edit unless reference paths change
# ==============================================================================
case "$GENOME" in
  galgal7)
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal7/star_index_main"
    ANNOTATION="/Data/GENOMES/GallusGallus/galgal7/galgal7.mainchr.gtf"
    ;;
  galgal6)
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal6/star_index_100bp/"
    ANNOTATION="/Data/GENOMES/GallusGallus/galgal6/galGal6.ncbiRefSeq.gtf"
    ;;
  hg38)
    GENOME_INDEX="/Data/GENOMES/HomoSapiens/hg38.p14/star_index_100bp/"
    ANNOTATION="/Data/GENOMES/HomoSapiens/hg38.p14/gencode.v49.annotation.gtf"
    ;;
  mm39)
    GENOME_INDEX="/Data/GENOMES/MusMusculus/mm39/star_index_100bp/"
    ANNOTATION="/Data/GENOMES/MusMusculus/mm39/gencode.vM38.annotation.gtf"
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
STATS_DIR="stats"
COUNTS_DIR="counts"
filelist="$FASTQ_DIR/fileR1.txt"

# ==============================================================================
# REQUIRED TOOLS
# ==============================================================================
required_tools=(featureCounts)
if (( START_STAGE <= 3 )); then
    required_tools+=(STAR samtools bamCoverage)
fi
if (( START_STAGE <= 2 )); then
    required_tools+=(fastqc cutadapt)
fi
for cmd in "${required_tools[@]}"; do
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
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$BW_DIR" "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed" "$COUNTS_DIR"
shopt -s nullglob

# ==============================================================================
# FASTQ DISCOVERY
# ==============================================================================
for sample in "${SAMPLE_ID[@]}"; do
    matches=("$ORIGINAL_FASTQ_PATH"/${sample}*"$R1_PATTERN")
    if (( ${#matches[@]} == 0 )); then
        echo "WARNING: No FASTQ for ${sample}*${R1_PATTERN}" >&2
        continue
    fi
    for fq in "${matches[@]}"; do
        ln -sfn "$fq" "$FASTQ_DIR/$(basename "$fq")"
    done
done

# ==============================================================================
# FASTQ MANIFEST
# ==============================================================================
: > "$filelist"
for R1path in "$FASTQ_DIR/"*"$R1_PATTERN"; do
    [[ -e "$R1path" ]] || continue
    echo "$(basename "$R1path")" >> "$filelist"
done

# Exit if no fastqs
if [[ ! -s "$filelist" ]]; then
    echo "ERROR: No R1 FASTQ files in $FASTQ_DIR." >&2
    exit 1
fi

# ==============================================================================
# RAW READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 1 )); then
echo "[QC:raw]"
fastqc -t "$THREADS" "$FASTQ_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_raw" || true
fi

# ==============================================================================
# ADAPTER AND POLY(A) TRIMMING
# ==============================================================================
if (( START_STAGE <= 2 )); then
echo "[Trim] minimum length: $MIN_LEN"
while IFS= read -r F1; do
    echo "[Trim] $F1"
    sample_name="${F1%${R1_PATTERN}}"
    REPORT="$STATS_DIR/${sample_name}_cutadapt_report.txt"

    # Illumina 3' adapter + trim long polyA tails of length >=10
    cutadapt -j "$THREADS" \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -a A{10} \
      --trim-n \
      --minimum-length "$MIN_LEN" \
      -o "$TRIM_DIR/${sample_name}_R1.trim.fastq.gz" \
      "$FASTQ_DIR/$F1" > "$REPORT"
done < "$filelist"
fi

# ==============================================================================
# TRIMMED READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 2 )); then
echo "[QC:trimmed]"
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_trimmed" || true
fi

# ==============================================================================
# STAR ALIGNMENT
# ==============================================================================
if (( START_STAGE <= 3 )); then
    trimmed_fastqs=("$TRIM_DIR"/*.fastq.gz)
    (( ${#trimmed_fastqs[@]} > 0 )) || { echo "ERROR: No trimmed FASTQ files found for alignment." >&2; exit 1; }
while IFS= read -r F1; do
    sample_name="${F1%${R1_PATTERN}}"
    echo "[Align:STAR] $sample_name"

    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$GENOME_INDEX" \
      --readFilesIn "$TRIM_DIR/${sample_name}_R1.trim.fastq.gz" \
      --readFilesCommand zcat \
      --sjdbGTFfile "$ANNOTATION" \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "$BAM_DIR/${sample_name}.${GENOME}_" \
      --quantMode TranscriptomeSAM

    # index BAM
    samtools index "$BAM_DIR/${sample_name}.${GENOME}_Aligned.sortedByCoord.out.bam"
    echo "[Track:CPM] $sample_name"
    bamCoverage \
        --bam "$BAM_DIR/${sample_name}.${GENOME}_Aligned.sortedByCoord.out.bam" \
        --outFileName "$BW_DIR/${sample_name}.${GENOME}.CPM.bw" \
        --outFileFormat bigwig \
        --binSize "$TRACK_BIN_SIZE" \
        --normalizeUsing CPM \
        --numberOfProcessors "$THREADS"
done < "$filelist"
fi
# ==============================================================================
# FEATURECOUNTS
# ==============================================================================
# Using annotation GTF. featureCounts will count exons grouped by gene_id.
bam_list=("$BAM_DIR"/*_Aligned.sortedByCoord.out.bam)
if (( ${#bam_list[@]} == 0 )); then
    echo "ERROR: STAR produced no coordinate-sorted BAM files." >&2
    exit 1
fi
case "$FEATURECOUNTS_MODE" in
  matrix)
    echo "[Count:matrix] ${#bam_list[@]} sample(s)"
    featureCounts \
        -T "$THREADS" -a "$ANNOTATION" \
        -o "$COUNTS_DIR/all_samples.${GENOME}.featureCounts.txt" \
        -t exon -g gene_id -s "$FEATURECOUNTS_STRANDED" \
        "${bam_list[@]}"
    ;;
  per_sample)
    for bam in "${bam_list[@]}"; do
        sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
        echo "[Count:sample] $sample_name"
        featureCounts \
            -T "$THREADS" -a "$ANNOTATION" \
            -o "$COUNTS_DIR/${sample_name}.featureCounts.txt" \
            -t exon -g gene_id -s "$FEATURECOUNTS_STRANDED" \
            "$bam"
    done
    ;;
  *)
    echo "ERROR: FEATURECOUNTS_MODE must be 'matrix' or 'per_sample'." >&2
    exit 1
    ;;
esac

# ==============================================================================
# CLEANUP
# ==============================================================================
if [[ "$KEEP_TRIMMED_FASTQ" == false ]] && (( START_STAGE <= 3 )); then
    echo "[Cleanup]"
    rm -rf "$TRIM_DIR"
fi
find "$BAM_DIR" -type f -name "*Log.out" -delete
find "$BAM_DIR" -type f -name "*Log.progress.out" -delete
find "$BAM_DIR" -type f -name "*SJ.out.tab" -delete
find "$BAM_DIR" -type d -name "*_STARgenome" -exec rm -rf {} +

echo "[Done] Pipeline complete"
