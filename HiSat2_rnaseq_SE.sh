#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# USER CONFIGURATION — edit only this section
# ==============================================================================
GENOME="galgal7"                 # galgal7, galgal6, hg38, or mm39
THREADS=10                       # CPU threads
ORIGINAL_FASTQ_PATH="/Data/Ana/NextSeq1000/SCL_10102025_100bp_AA/fastq"
SAMPLE_ID=("CRISPR_3RNA")
R1_PATTERN=".fastq.gz"           # Exact filename suffix for single-end reads
MIN_LEN=20                       # Minimum read length after trimming
HISAT_STRANDNESS="F"             # F, R, or "" for unstranded
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
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal7/hisat2_index_main/galgal7_mainchr"
    ANNOTATION="/Data/GENOMES/GallusGallus/galgal7/galgal7.mainchr.gtf"
    ;;
  galgal6)
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal6/hisat2_index/galgal6"
    ANNOTATION="/Data/GENOMES/GallusGallus/galgal6/galGal6.ncbiRefSeq.gtf"
    ;;
  hg38)
    GENOME_INDEX="/Data/GENOMES/HomoSapiens/hg38.p14/hisat2_index/hg38.p14"
    ANNOTATION="/Data/GENOMES/HomoSapiens/hg38.p14/gencode.v49.annotation.gtf"
    ;;
  mm39)
    GENOME_INDEX="/Data/GENOMES/MusMusculus/mm39/hisat2_index/mm39"
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
TRIM_SUFFIX="_R1.trim.fastq.gz"

# ==============================================================================
# OUTPUT DIRECTORY SETUP
# ==============================================================================
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$BW_DIR" "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed" "$COUNTS_DIR"
shopt -s nullglob

hisat_strand_args=()
if [[ -n "$HISAT_STRANDNESS" ]]; then
    hisat_strand_args=(--rna-strandness "$HISAT_STRANDNESS")
fi

# ==============================================================================
# REQUIRED TOOLS FOR THE SELECTED START STAGE
# ==============================================================================
required_tools=(featureCounts)
if (( START_STAGE <= 3 )); then
    required_tools+=(hisat2 samtools bamCoverage)
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

# ==============================================================================
# FASTQ DISCOVERY AND MANIFEST
# ==============================================================================
fastq_list=()
trimmed_fastqs=()
bam_list=()
if (( START_STAGE <= 2 )); then
    [[ -d "$ORIGINAL_FASTQ_PATH" ]] || {
        echo "ERROR: FASTQ directory does not exist: $ORIGINAL_FASTQ_PATH" >&2
        exit 1
    }

    : > "$filelist"
    for sample in "${SAMPLE_ID[@]}"; do
        matches=("$ORIGINAL_FASTQ_PATH"/${sample}*"$R1_PATTERN")
        if (( ${#matches[@]} == 0 )); then
            echo "WARNING: No FASTQ for ${sample}*${R1_PATTERN}" >&2
            continue
        fi
        for fq in "${matches[@]}"; do
            fq_name=$(basename "$fq")
            if ! grep -Fqx -- "$fq_name" "$filelist"; then
                ln -sfn "$fq" "$FASTQ_DIR/$fq_name"
                echo "$fq_name" >> "$filelist"
                fastq_list+=("$FASTQ_DIR/$fq_name")
            fi
        done
    done

    if [[ ! -s "$filelist" ]]; then
        echo "ERROR: No FASTQ files were found. Check SAMPLE_ID and R1_PATTERN." >&2
        exit 1
    fi
fi

# ==============================================================================
# RAW READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 1 )); then
echo "[QC:raw]"
fastqc -t "$THREADS" "${fastq_list[@]}" -o "$STATS_DIR/fastqc_raw"
fi

# ==============================================================================
# READ TRIMMING
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
      -o "$TRIM_DIR/${sample_name}${TRIM_SUFFIX}" \
      "$FASTQ_DIR/$F1" > "$REPORT"
    trimmed_fastqs+=("$TRIM_DIR/${sample_name}${TRIM_SUFFIX}")
done < "$filelist"
fi

# ==============================================================================
# TRIMMED READ QUALITY CONTROL
# ==============================================================================
if (( START_STAGE <= 2 )); then
echo "[QC:trimmed]"
fastqc -t "$THREADS" "${trimmed_fastqs[@]}" -o "$STATS_DIR/fastqc_trimmed"
fi

# ==============================================================================
# HISAT2 ALIGNMENT
# ==============================================================================
if (( START_STAGE <= 3 )); then
if (( ${#trimmed_fastqs[@]} == 0 )); then
    trimmed_fastqs=("$TRIM_DIR"/*"$TRIM_SUFFIX")
fi
(( ${#trimmed_fastqs[@]} > 0 )) || { echo "ERROR: No trimmed FASTQ files found for alignment." >&2; exit 1; }
for fq in "${trimmed_fastqs[@]}"; do
    fq_name=$(basename "$fq")
    if [[ "$fq_name" != *"$TRIM_SUFFIX" ]]; then
        echo "ERROR: Unexpected trimmed FASTQ name: $fq_name" >&2
        exit 1
    fi
    sample="${fq_name%${TRIM_SUFFIX}}"
    echo "[Align:HISAT2] $sample"
    hisat2 \
        -p "$THREADS" \
        --phred33 \
        "${hisat_strand_args[@]}" \
        --dta \
        --no-unal \
        -x "$GENOME_INDEX" \
        -U "$fq" \
        --summary-file "$STATS_DIR/${sample}.${GENOME}_hisatSummary.txt" \
        | samtools view -b -@ "$THREADS" - \
        | samtools sort -@ "$THREADS" -o "$BAM_DIR/${sample}.${GENOME}.bam"
    samtools index -@ "$THREADS" "$BAM_DIR/${sample}.${GENOME}.bam"
    echo "[Track:CPM] $sample"
    bamCoverage \
        --bam "$BAM_DIR/${sample}.${GENOME}.bam" \
        --outFileName "$BW_DIR/${sample}.${GENOME}.CPM.bw" \
        --outFileFormat bigwig \
        --binSize "$TRACK_BIN_SIZE" \
        --normalizeUsing CPM \
        --numberOfProcessors "$THREADS"
    bam_list+=("$BAM_DIR/${sample}.${GENOME}.bam")
done
fi

# ==============================================================================
# FEATURECOUNTS
# ==============================================================================
if (( ${#bam_list[@]} == 0 )); then
    bam_list=("$BAM_DIR"/*.bam)
fi
if (( ${#bam_list[@]} == 0 )); then
    echo "ERROR: HISAT2 produced no BAM files." >&2
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
        sample_name=$(basename "$bam" .bam)
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

if [[ "$KEEP_TRIMMED_FASTQ" == false ]] && (( START_STAGE <= 3 )); then
    echo "[Cleanup]"
    rm -rf "$TRIM_DIR"
fi

echo "[Done] Pipeline complete"
