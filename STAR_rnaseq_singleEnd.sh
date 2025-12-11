#!/bin/bash
set -euo pipefail

# ===============================================
# STAR pipeline + Lexogen QuantSeq 3' FWD
# Single-end sequencing (demultiplex: 100x12x12x0 -> R1 = 100bp)
# ===============================================

# =============================================== PUT YOUR SETTINGS HERE, THE SCRIPT WILL DO THE REST :) ===============================================
# ===== CONFIG =====
### possible GENOME choices = galgal6, hg38 or mm39
GENOME=galgal6
THREADS=16

# path where demultiplexed fastq are (R1 files only)
ORIGINAL_FASTQ_PATH="/Data/Ana/NextSeq1000/SCL_10102025_100bp_AA/fastq"
# sample ids as in your SampleSheet (patterns that match R1 fastq basename)
SAMPLE_ID=("CRISPR_3RNA_CTCF")
# ======================================================================================================================================================


# ===== GENOME SETTINGS (do not change names without checking paths) =====
case "$GENOME" in
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
    echo "Error: Unknown genome '$GENOME'. Choose from {galgal6, hg38, mm39}."
    exit 1
    ;;
esac

# ===== PATHS =====
FASTQ_DIR="fastq"
TRIM_DIR="trimmedFastq"
BAM_DIR="BAM"
STATS_DIR="stats"
COUNTS_DIR="counts"
filelist="$FASTQ_DIR/fileR1.txt"

# minimum trimmed read length to keep
MIN_LEN=20

# featureCounts strandedness default for Lexogen QuantSeq FWD
# 1 -> stranded (reads are in the same orientation as transcripts)
# 2 -> reverse stranded (if your data behaves like that, change to 2)
# 0 -> unstranded
FEATURECOUNTS_STRANDED=1

# ===== TOOL CHECK (fail fast if missing) =====
for cmd in fastqc cutadapt STAR samtools featureCounts multiqc; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# ===== PREP =====
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed" "$COUNTS_DIR"
shopt -s nullglob

# ===== SYMLINK fastq files =====
for sample in "${SAMPLE_ID[@]}"; do
    matches=("$ORIGINAL_FASTQ_PATH"/${sample}*R1*.fastq.gz)
    if (( ${#matches[@]} == 0 )); then
        echo "Warning: no fastq files found for pattern ${sample}*R1* in $ORIGINAL_FASTQ_PATH" >&2
        continue
    fi
    for fq in "${matches[@]}"; do
        ln -sfn "$fq" "$FASTQ_DIR/$(basename "$fq")"
    done
done

# ===== Make list of single-end R1 files =====
: > "$filelist"
for R1path in "$FASTQ_DIR/"*R1*.fastq.gz; do
    [[ -e "$R1path" ]] || continue
    echo "$(basename "$R1path")" >> "$filelist"
done

# Exit if no fastqs
if [[ ! -s "$filelist" ]]; then
    echo "No R1 fastq files found in $FASTQ_DIR. Exiting." >&2
    exit 1
fi

# ===== FASTQC (raw data) =====
echo "[FastQC] raw fastq..."
fastqc -t "$THREADS" "$FASTQ_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_raw" || true

# ===== ADAPTER + polyA TRIMMING (single-end) =====
echo "[Cutadapt] trimming adapters + polyA tails (min length = $MIN_LEN)"
while IFS= read -r F1; do
    echo "[Cutadapt] $F1"
    sample_name=$(basename "$F1" | sed -E 's/_R1.*//')
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

# ===== FASTQC after trimming =====
echo "[FastQC] trimmed fastq..."
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_trimmed" || true

# ===== STAR ALIGNMENT (single-end) =====
while IFS= read -r F1; do
    sample_name=$(basename "$F1" | sed -E 's/_R1.*//')
    echo "[STAR] Aligning $sample_name..."

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
done < "$filelist"

# ===== featureCounts (gene counts) =====
# Using annotation GTF. featureCounts will count exons grouped by gene_id.
echo "[featureCounts] computing gene counts..."
# build list of BAMs
for sample in "${SAMPLE_ID[@]}"; do
	bamlist=()
	for bam in "$BAM_DIR"/${sample}*Aligned.sortedByCoord.out.bam; do
    		[[ -e "$bam" ]] || continue
    		bamlist+=("$bam")
	done

	# Run featureCounts
	featureCounts -T "$THREADS" -a "$ANNOTATION" -o "$COUNTS_DIR"/${sample}.${GENOME}.featureCounts_counts.txt \
	  -t exon -g gene_id -s "$FEATURECOUNTS_STRANDED" "${bamlist[@]}"
done

# ===== CLEANUP (highly recommended) =====
rm -rf "$TRIM_DIR"
find "$BAM_DIR" -type f -name "*Aligned.toTranscriptome.out.bam" -delete
find "$BAM_DIR" -type f -name "*Log.out" -delete
find "$BAM_DIR" -type f -name "*Log.progress.out" -delete
find "$BAM_DIR" -type f -name "*SJ.out.tab" -delete
find "$BAM_DIR" -type d -name "*_STARgenome" -exec rm -rf {} +

echo "Pipeline complete!"
