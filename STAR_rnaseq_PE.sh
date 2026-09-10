#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# USER CONFIGURATION — edit only this section
# ==============================================================================
GENOME="galgal7"                 # galgal7, galgal6, hg38, or mm39
THREADS=10                       # CPU threads
ORIGINAL_FASTQ_PATH="/mnt/rcfs/Public/FASTQS/raw_data/NextSeq1000/SCL_09102024_NAAASL/"
SAMPLE_ID=("HH4_anterior_ectoderm_RNA")
R1_PATTERN="_R1_001.fastq.gz"
R2_PATTERN="_R2_001.fastq.gz"
FEATURECOUNTS_STRANDED=2         # 0=unstranded, 1=forward, 2=reverse
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
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal6/star_index_76bp/"
    ANNOTATION="/Data/GENOMES/GallusGallus/galgal6/galGal6.ncbiRefSeq.gtf"
    ;;
  hg38)
    GENOME_INDEX="/Data/GENOMES/HomoSapiens/hg38.p14/star_index/"
    ANNOTATION="/Data/GENOMES/HomoSapiens/hg38.p14/gencode.v49.annotation.gtf"
    ;;
  mm39)
    GENOME_INDEX="/Data/GENOMES/MusMusculus/mm39/star_index/"
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
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$BW_DIR" \
         "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed" "$COUNTS_DIR"

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
    REPORT="$STATS_DIR/$(basename "$F1" .fastq.gz)_cutadapt_report.txt"
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --minimum-length=25 -j "$THREADS" \
        -o "$TRIM_DIR/trimmed_${F1}" -p "$TRIM_DIR/trimmed_${F2}" \
        "$FASTQ_DIR/$F1" "$FASTQ_DIR/$F2" > "$REPORT"
done < "$filepairs"
fi

# ==============================================================================
# TRIMMED READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 2 )); then
echo "[QC:trimmed]"
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq* -o "$STATS_DIR/fastqc_trimmed"
fi

# ==============================================================================
# STAR ALIGNMENT
# ==============================================================================
if (( START_STAGE <= 3 )); then
    trimmed_fastqs=("$TRIM_DIR"/*.fastq.gz)
    (( ${#trimmed_fastqs[@]} > 0 )) || { echo "ERROR: No trimmed FASTQ files found for alignment." >&2; exit 1; }
while IFS=";" read -r F1 F2; do
    sample_name="${F1%${R1_PATTERN}}"
    echo "[Align:STAR] $sample_name"
    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$GENOME_INDEX" \
      --readFilesIn "$TRIM_DIR/trimmed_${F1}" "$TRIM_DIR/trimmed_${F2}" \
      --readFilesCommand zcat \
      --sjdbGTFfile "$ANNOTATION" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix "$BAM_DIR/${sample_name}.${GENOME}_"

    samtools index "$BAM_DIR/${sample_name}.${GENOME}_Aligned.sortedByCoord.out.bam"
    echo "[Track:CPM] $sample_name"
    bamCoverage \
        --bam "$BAM_DIR/${sample_name}.${GENOME}_Aligned.sortedByCoord.out.bam" \
        --outFileName "$BW_DIR/${sample_name}.${GENOME}.CPM.bw" \
        --outFileFormat bigwig \
        --binSize "$TRACK_BIN_SIZE" \
        --normalizeUsing CPM \
        --numberOfProcessors "$THREADS"
done < "$filepairs"
fi
# ==============================================================================
# FEATURECOUNTS
# ==============================================================================
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
        -p --countReadPairs "${bam_list[@]}"
    ;;
  per_sample)
    for bam in "${bam_list[@]}"; do
        sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
        echo "[Count:sample] $sample_name"
        featureCounts \
            -T "$THREADS" -a "$ANNOTATION" \
            -o "$COUNTS_DIR/${sample_name}.featureCounts.txt" \
            -t exon -g gene_id -s "$FEATURECOUNTS_STRANDED" \
            -p --countReadPairs "$bam"
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
    rm -rf "$TRIM_DIR"
fi
find "$BAM_DIR" -type f -name "*Log.out" -delete
find "$BAM_DIR" -type f -name "*Log.progress.out" -delete
find "$BAM_DIR" -type f -name "*SJ.out.tab" -delete
find "$BAM_DIR" -type d -name "*_STARgenome" -exec rm -rf {} +

echo "[Done] Pipeline complete"
