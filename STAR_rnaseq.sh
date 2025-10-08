#!/bin/bash
set -euo pipefail
# =============================================== PUT YOUR SETTINGS HERE, THE SCRIPT WILL DO THE REST :) ===============================================
# ===== CONFIG =====
### possible GENOME choices = galgal6, hg38 or mm39
GENOME=galgal6
THREADS=10

# the path ANA will send you
ORIGINAL_FASTQ_PATH="/Data/FASTQS/raw_data/NextSeq1000/SCL_09102024_NAAASL/"
# your samples id as in SampleSheet
SAMPLE_ID=("HH4_anterior_ectoderm_RNA_1")

# ======================================================================================================================================================



# ==== GENOME SPECS ====
## DO NOT CHANGE THIS!
case "$GENOME" in
  galgal6)
    GENOME_INDEX="/Data/GENOMES/GallusGallus/galgal6/star_index_76bp/"
    ANNOTATION="/Data/GENOMES/GallusGallus/galgal6/galGal6.refGene.gtf"
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
    echo "Error: Unknown genome '$GENOME'. Choose from {galgal6, hg38, mm39}."
    exit 1
    ;;
esac

# ===== PATHS =====
FASTQ_DIR="fastq"
TRIM_DIR="trimmedFastq"
BAM_DIR="BAM"
STATS_DIR="stats"

# ===== TOOL CHECK (fail fast if missing) =====
for cmd in fastqc cutadapt STAR samtools; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# ===== PREP =====
mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$BAM_DIR" \
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
    REPORT="$STATS_DIR/$(basename "$F1" .fastq.gz)_cutadapt_report.txt"
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        --minimum-length=25 -j "$THREADS" \
        -o "$TRIM_DIR/trimmed_${F1}" -p "$TRIM_DIR/trimmed_${F2}" \
        "$FASTQ_DIR/$F1" "$FASTQ_DIR/$F2" > "$REPORT"
done < "$filepairs"

# ===== FASTQC after trimming =====
echo "[FastQC] Executing in trimmed fastq..."
fastqc -t "$THREADS" "$TRIM_DIR"/*.fastq* -o "$STATS_DIR/fastqc_trimmed"

# ===== STAR ALIGNMENT =====
while IFS=";" read -r F1 F2; do
    sample_name=$(basename "$F1" | sed 's/_R1_.*//')
    echo "[STAR] Aligning $sample_name..."
    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$GENOME_INDEX" \
      --readFilesIn "$TRIM_DIR/trimmed_${F1}" "$TRIM_DIR/trimmed_${F2}" \
      --readFilesCommand zcat \
      --sjdbGTFfile "$ANNOTATION" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix "$BAM_DIR/${sample_name}_"
    
    samtools index "$BAM_DIR/${sample_name}_Aligned.sortedByCoord.out.bam"
done < "$filepairs"


# ===== CLEAN INTERMEDIATES =====
echo "[Cleanup] Removing intermediate files..."
rm -rf "$TRIM_DIR"

echo "Pipeline complete!"
