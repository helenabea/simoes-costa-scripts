echo "Program Started: $(date)" > timelog.txt

mkdir stats
mkdir trimmedFastq
mkdir BAM
mkdir counts

fastqc ./fastq/*.fastq -o ./stats -t 8

for FILE in ./fastq/*.fastq
    do
        NUM=$(wc -l $FILE | cut -f1 -d" ")
        NUMR=$(echo "$NUM"/4 | bc)
        echo -e "${FILE}\t${NUMR}" >> ./stats/readCounts.txt
    done
wait

for FILEPAIR in $(cat filePairs.txt)
    do
        F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
        F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
        echo "Paired End mode ${FILEPAIR}"
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --minimum-length=25 -j 8 \
        -o "trimmed_${F1}"  -p "trimmed_${F2}" \
        ${F1} ${F2} >"${F1}_cutadapt_report.txt"
    done
wait

echo "Adapters Trimmed: $(date)" >> timelog.txt

for FILE in ./fastq/trimmed*
	do
	   {
		mv $FILE ./trimmedFastq
	   } &
	done
wait

for FILE in ./fastq/*_cutadapt_report.txt
	do
	{
		mv $FILE ./stats
	} &
	done
wait

cat ./stats/*_cutadapt_report.txt > ./stats/CombinedTrimmingStats.txt

for FILE in ./trimmedFastq/*.fastq
    do
        NUM=$(wc -l $FILE | cut -f1 -d" ")
        NUMR=$(echo "$NUM"/4 | bc)
        echo -e "${FILE}\t${NUMR}" >> ./stats/readCountsTrimmed.txt
    done
wait



for FILEPAIR in $(cat filePairs.txt)
    do
    	F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
        F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
    	hisat2 \
    	-p 8 \
    	--phred33 \
    	--rna-strandness R \
    	--dta \
    	--no-unal \
        -1 "trimmed_${F1}" \
        -2 "trimmed_${F2}" \
    	-x /data/fm449/MRE_PROJECT/mouse_genome/ens_mus_grc39_hisat2_index/mus_grc39_hisat2_index \
        --summary-file ../stats/${FN/hisatSummary} |
    	samtools view -bS -@ 8 - |
    	samtools sort -@ 8 - > "../BAM/${F1/fastq/bam}" &&
    	samtools index -@ 8 "../BAM/${F1/fastq/bam}"
    done
wait

echo "Reads Aligned: $(date)" >> timelog.txt


     featureCounts \
	./BAM/*.bam \
	-T 8 \
	-t exon \
	-s 2 \
	-g gene_id \
	-a /data/fm449/MRE_PROJECT/mouse_genome/Mus_musculus.GRCm39.111.gtf \
	-o ./counts/featureCounts.txt \
	-p \

echo "Reads Counted at $(date)" >> timelog.txt

multiqc ./stats --cl_config "extra_fn_clean_exts: { '_trimmed.BAM'}"

echo "MultiQC Report Generated at $(date)...Complete!" >> timelog.txt
