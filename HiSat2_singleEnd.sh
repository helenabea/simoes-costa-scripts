#!/bin/bash
set -euo pipefail

# ===============================================
# CONFIG
# ===============================================
GENOME=galgal6
THREADS=10
ORIGINAL_FASTQ_PATH="/Data/FASTQS/raw_data/NextSeq1000/SCL_09102024_NAAASL/"
SAMPLE_ID=("HH4_anterior_ectoderm_RNA_1")
MIN_LEN=20
FEATURECOUNTS_STRANDED=1
### For single end use F or R For paired-end use FR OR RF - Lexogen is FOWARD (F)
HISAT_STRANDNESS=F

# ===============================================
# GENOME PATHS
# ===============================================
case "$GENOME" in
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
    echo "Error: Unknown genome '$GENOME'." >&2
    exit 1
    ;;
esac

# ===============================================
# PATHS
# ===============================================
FASTQ_DIR="fastq"
TRIM_DIR="trimmedFastq"
BAM_DIR="BAM"
STATS_DIR="stats"
COUNTS_DIR="counts"

mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" "$STATS_DIR/fastqc_raw" "$STATS_DIR/fastqc_trimmed" "$COUNTS_DIR"

# ===============================================
# LINK FASTQ FILES
# ===============================================
for sample in "${SAMPLE_ID[@]}"; do
    matches=("$ORIGINAL_FASTQ_PATH"/${sample}*.fastq.gz)
    if (( ${#matches[@]} == 0 )); then
        echo "Warning: no fastq files found for ${sample}" >&2
        continue
    fi
    for fq in "${matches[@]}"; do
        ln -sfn "$fq" "$FASTQ_DIR/$(basename "$fq")"
    done
done

# ===============================================
# FASTQC (raw)
# ===============================================
echo "[FastQC] Raw reads..."
fastqc -t "$THREADS" "$FASTQ_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_raw"

# ===============================================
# TRIMMING (fastp)
# ===============================================
for fq in "$FASTQ_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "[fastp] Trimming $sample..."
    fastp \
        -i "$fq" \
        -o "$TRIM_DIR/${sample}_trimmed.fastq.gz" \
        --qualified_quality_phred 20 \
        --length_required "$MIN_LEN" \
        --thread "$THREADS" \
        --html "$STATS_DIR/${sample}_fastp.html" \
        --json "$STATS_DIR/${sample}_fastp.json"
done

# ===============================================
# FASTQC (trimmed)
# ===============================================
echo "[FastQC] Trimmed reads..."
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq.gz -o "$STATS_DIR/fastqc_trimmed"

# ===============================================
# ALIGNMENT (HISAT2)
# ===============================================
for fq in "$TRIM_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" _trimmed.fastq.gz)
    echo "[HISAT2] Aligning $sample..."
    hisat2 \
        -p "$THREADS" \
        --phred33 \
        --rna-strandness $HISAT_STRANDNESS \
        --dta \
        -x "$GENOME_INDEX" \
        -U "$fq" \
        --summary-file "$STATS_DIR/${sample}.${GENOME}_hisatSummary.txt" \
        | samtools view -bS -@ "$THREADS" - \
        | samtools sort -@ "$THREADS" -o "$BAM_DIR/${sample}.${GENOME}.bam"
    samtools index -@ "$THREADS" "$BAM_DIR/${sample}.${GENOME}.bam"
done

# ===============================================
# FEATURE COUNTS
# ===============================================
for bam in "$BAM_DIR"/*.bam; do
    base=$(basename "$bam" .bam)
    echo "[featureCounts] Counting for $base..."
    featureCounts \
        -T "$THREADS" \
        -t exon \
        -s "$FEATURECOUNTS_STRANDED" \
        -g gene_id \
        -a "$ANNOTATION" \
        -o "$COUNTS_DIR/${base}_featureCounts.txt" \
        "$bam"
done

echo "Pipeline complete!"
