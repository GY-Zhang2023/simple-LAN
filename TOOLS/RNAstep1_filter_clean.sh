#!/bin/bash
THREADS_DEFAULT=8
CORES_DEFAULT=2
MEM_DEFAULT=50G
WORK_DIR=$(pwd)
MT_INDEX_HUMAN="$HOME/ribo_references/rrna_trna/human_mtDNA"
MT_INDEX_MOUSE="$HOME/reference/mouse_transcriptome/mouse_mtDNA"
RRNA_TRNA_INDEX_HUMAN="$HOME/ribo_references/rrna_trna/human_rrna_tRNA_index"
RRNA_TRNA_INDEX_MOUSE="$HOME/reference/mouse_transcriptome/mouse_rrna_tRNA_index"
THREADS_BOWTIE2=${THREADS_DEFAULT}
THREADS_FASTQC=${THREADS_DEFAULT}
THREADS_CUTADAPT=${THREADS_DEFAULT}
THREADS_TRIMGALORE=${THREADS_DEFAULT}
MEM_LIMIT_BYTES=$(( ${MEM_DEFAULT%G} * 1024 * 1024 * 1024 ))
ulimit -v ${MEM_LIMIT_BYTES}
declare -A SPECIES_MT_INDEX_MAP
SPECIES_MT_INDEX_MAP["human"]="${MT_INDEX_HUMAN}"
SPECIES_MT_INDEX_MAP["mouse"]="${MT_INDEX_MOUSE}"
declare -A SPECIES_RRRNA_TRNA_INDEX_MAP
SPECIES_RRRNA_TRNA_INDEX_MAP["human"]="${RRNA_TRNA_INDEX_HUMAN}"
SPECIES_RRRNA_TRNA_INDEX_MAP["mouse"]="${RRNA_TRNA_INDEX_MOUSE}"
if [ $# -ne 1 ]; then
    echo "Usage: $0 <CSV filename>"
    echo "CSV file should be placed under 1.raw_fastq/ in current directory"
    echo "Example: $0 RNA-seq_sample_list.csv"
    exit 1
fi
CSV_ABS_PATH="$1"
CSV="$1"
if [ ! -f "${CSV_ABS_PATH}" ]; then
    echo "ERROR: CSV file not found: ${CSV_ABS_PATH}"
    exit 1
fi
if [ ! -s "${CSV_ABS_PATH}" ]; then
    echo "ERROR: CSV file is empty: ${CSV_ABS_PATH}"
    exit 1
fi
echo -e "\n==================== Species selection ===================="
echo "Select species (input number and press enter):"
echo "1 human"
echo "2 mouse"
read -p "Input selection (1/2): " SPECIES_CHOICE
case ${SPECIES_CHOICE} in
    1)
        SELECTED_SPECIES="human"
        ;;
    2)
        SELECTED_SPECIES="mouse"
        ;;
    *)
        echo "ERROR: invalid selection, only 1 (human) or 2 (mouse) allowed"
        exit 1
        ;;
esac
mtCONTAMINANT_INDEX=${SPECIES_MT_INDEX_MAP[${SELECTED_SPECIES}]}
rRNA_tRNA_CONTAMINANT_INDEX=${SPECIES_RRRNA_TRNA_INDEX_MAP[${SELECTED_SPECIES}]}
echo -e "\nConfirmed species: ${SELECTED_SPECIES}"
echo "mtRNA contaminant index: ${mtCONTAMINANT_INDEX}"
echo "rRNA+tRNA contaminant index: ${rRNA_tRNA_CONTAMINANT_INDEX}"
read -p "Continue? (y/n): " CONFIRM
if [ "${CONFIRM}" != "y" ] && [ "${CONFIRM}" != "Y" ]; then
    echo "User cancelled, exit"
    exit 0
fi
RAW_DATA_DIR="${WORK_DIR}/1.raw_fastq"
CLEAN_DATA_DIR="${WORK_DIR}/2.clean_fastq"
STATS_FILE="${CLEAN_DATA_DIR}/rna_seq_filter_stats.csv"
FASTQC_CLEAN_DIR="${CLEAN_DATA_DIR}/fastqc_clean"
CONTAM_STATS_FILE="${CLEAN_DATA_DIR}/contaminant_removal_stats.csv"
echo -e "\n==================== Run confirmation ===================="
echo "CSV absolute path: ${CSV_ABS_PATH}"
echo "Raw data directory: ${RAW_DATA_DIR}"
echo "Selected species: ${SELECTED_SPECIES}"
echo "mtRNA index: ${mtCONTAMINANT_INDEX}"
echo "rRNA-tRNA index: ${rRNA_tRNA_CONTAMINANT_INDEX}"
echo "Memory limit: ${MEM_DEFAULT}"
mkdir -p "$CLEAN_DATA_DIR" "$FASTQC_CLEAN_DIR" || { echo "Failed to create directories!"; exit 1; }
echo "sample_name,sequencing_type,raw_reads,after_trim_qc_reads,after_length_filter_reads,after_contaminant_filter_reads,final_reads,min_length,max_length,avg_length,rRNA_removal_rate,mtRNA_removal_rate" > "$STATS_FILE"
echo "sample_name,total_contaminant_reads,rRNA_reads,mtRNA_reads,tRNA_reads,contaminant_removal_rate" > "$CONTAM_STATS_FILE"
declare -A sample_info
declare -A is_paired_end
sample_list=()
dos2unix "$CSV" 2>/dev/null
while IFS= read -r line; do
    line=$(echo "$line" | sed -e 's/^[ \t\r]*//' -e 's/[ \t\r]*$//')
    [ -z "$line" ] && continue
    if [[ "$line" =~ ^[Ss][Aa][Mm][Pp][Ll][Ee]_?[Nn][Aa][Mm][Ee],.* ]]; then
        continue
    fi
    IFS=',' read -r sample adapter3 adapter5 <<< "$line"
    sample=$(echo "$sample" | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//' -e "s/'//g")
    adapter3=$(echo "$adapter3" | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//' -e "s/'//g" -e 's/^NA$//i')
    adapter5=$(echo "$adapter5" | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//' -e "s/'//g" -e 's/^NA$//i')
    [ -z "$sample" ] && continue
    if [[ "$sample" =~ _1$ ]] || [[ "$sample" =~ _2$ ]]; then
        is_paired_end["$sample"]="true"
    else
        is_paired_end["$sample"]="false"
    fi
    sample_info["$sample"]="$adapter3,$adapter5"
    sample_list+=("$sample")
done < "$CSV"
total=${#sample_list[@]}
paired_samples=()
single_samples=()
for sample in "${sample_list[@]}"; do
    if [[ "${is_paired_end[$sample]}" == "true" ]]; then
        paired_samples+=("$sample")
    else
        single_samples+=("$sample")
    fi
done
echo "===== RNA-seq filter pipeline start ====="
echo "Start time: $(date +'%Y-%m-%d %H:%M:%S')"
echo "Total samples: $total | single-end: ${#single_samples[@]} | paired-end: ${#paired_samples[@]}"
echo "Output directory: $CLEAN_DATA_DIR"
echo "Core rules:"
echo "  1. trim_galore -> cutadapt: adapter & QC, then length filter 36-150nt;"
echo "  2. Contaminant removal target: rRNA + tRNA + mtRNA;"
echo "  3. Switch to single-end decontamination if R1/R2 read count ratio > 5;"
echo "  4. QC threshold: q20;"
echo "=========================="
current=0
success=0
failed=0
process_single_end() {
    local sample=$1
    local adapter3=$2
    local adapter5=$3
    local raw_reads=0 after_trim_qc_reads=0 after_length_filter_reads=0 after_contaminant_filter_reads=0 final_reads=0
    local min_length=0 max_length=0 avg_length=0
    local contaminant_reads=0 rRNA_reads=0 mtRNA_reads=0 tRNA_reads=0 contaminant_removal_rate=0 rRNA_removal_rate=0 mtRNA_removal_rate=0
    echo -e "\n[$current/$total] Processing single-end sample: $sample"
    raw_file="$RAW_DATA_DIR/$sample.fastq.gz"
    if [ ! -f "$raw_file" ]; then
        echo "  ERROR: raw file not found: $raw_file"
        echo "$sample,single,0,0,0,0,0,,,,0,0" >> "$STATS_FILE"
        echo "$sample,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
        return 1
    fi
    raw_reads=$(zcat "$raw_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    echo "  Raw reads: $raw_reads"
    qc_trimmed_file="$CLEAN_DATA_DIR/${sample}_qc_trimmed.fastq.gz"
    length_filtered_file="$CLEAN_DATA_DIR/${sample}_length_filtered.fastq.gz"
    trim_ok=0
    length_ok=0
    if [ -n "$adapter3" ] || [ -n "$adapter5" ]; then
        echo "  Adapter detected, first step: cutadapt adapter trimming + QC (3'=$adapter3, 5'=$adapter5)..."
        cutadapt_args=()
        [ -n "$adapter3" ] && cutadapt_args+=("-a" "$adapter3")
        [ -n "$adapter5" ] && cutadapt_args+=("-g" "$adapter5")
        cutadapt "${cutadapt_args[@]}" -j $THREADS_CUTADAPT --trim-n -q 20 --trimmed-only -o "$qc_trimmed_file" "$raw_file" 2>/dev/null
        if [ -f "$qc_trimmed_file" ] && [ $(stat -c %s "$qc_trimmed_file") -gt 0 ]; then
            after_trim_qc_reads=$(zcat "$qc_trimmed_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            echo "  After trim & QC reads: $after_trim_qc_reads"
            echo "  Second step: cutadapt length filter 36-150nt..."
            cutadapt --minimum-length 36 --maximum-length 150 -j $THREADS_CUTADAPT -o "$length_filtered_file" "$qc_trimmed_file" 2>/dev/null
            if [ -f "$length_filtered_file" ] && [ $(stat -c %s "$length_filtered_file") -gt 0 ]; then
                after_length_filter_reads=$(zcat "$length_filtered_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                echo "  After QC + length filter reads: $after_length_filter_reads (36-150nt)"
                length_ok=1
                rm -f $qc_trimmed_file
            fi
        fi
    else
        echo "  No adapter, first step: trim_galore QC only..."
        trim_galore --quality 20 --trim-n --cores $THREADS_TRIMGALORE -o "$CLEAN_DATA_DIR" "$raw_file" 2>/dev/null
        qc_trimmed_tmp="${CLEAN_DATA_DIR}/${sample}_trimmed.fq.gz"
        [ ! -f "$qc_trimmed_tmp" ] && qc_trimmed_tmp="${CLEAN_DATA_DIR}/${sample}_val.fq.gz"
        if [ -f "$qc_trimmed_tmp" ]; then
            mv "$qc_trimmed_tmp" "$qc_trimmed_file"
            after_trim_qc_reads=$(zcat "$qc_trimmed_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            echo "  After QC reads: $after_trim_qc_reads"
            echo "  Second step: cutadapt length filter 36-150nt..."
            cutadapt --minimum-length 36 --maximum-length 150 -j $THREADS_CUTADAPT -o "$length_filtered_file" "$qc_trimmed_file" 2>/dev/null
            if [ -f "$length_filtered_file" ] && [ $(stat -c %s "$length_filtered_file") -gt 0 ]; then
                after_length_filter_reads=$(zcat "$length_filtered_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                echo "  After QC + length filter reads: $after_length_filter_reads (36-150nt)"
                length_ok=1
                rm -f $qc_trimmed_file
            fi
        fi
    fi
    if [ $length_ok -eq 0 ]; then
        echo "  ERROR: empty after length filter, mark sample failed"
        echo "$sample,single,$raw_reads,$after_trim_qc_reads,0,-1,0,,,0,0" >> "$STATS_FILE"
        echo "$sample,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
        rm -f "$qc_trimmed_file" "$length_filtered_file"
        return 1
    fi
    nomt_contam_file="$CLEAN_DATA_DIR/${sample}_nomt_contam.fastq"
    mtcontam_sam="$CLEAN_DATA_DIR/${sample}_mtcontaminant_mapped.sam"
    echo "  Running single-end mtRNA filtering..."
    bowtie2 -x "$mtCONTAMINANT_INDEX" -U "$length_filtered_file" --un "$nomt_contam_file" --no-unal -L 20 -p $THREADS_BOWTIE2 -S "$mtcontam_sam" 2>/dev/null
    rm -f $length_filtered_file
    if [ ! -f "$nomt_contam_file" ] || [ $(stat -c %s "$nomt_contam_file") -eq 0 ]; then
        echo "  ERROR: empty after mtRNA filter, mark sample failed"
        echo "$sample,single,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,0,0,,,0,0" >> "$STATS_FILE"
        echo "$sample,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
        rm -f "$mtcontam_sam" "$nomt_contam_file"
        return 1
    fi
    after_mt_filter_reads=$(cat "$nomt_contam_file" | wc -l | awk '{print int($1/4)}')
    mt_removed_reads=$((after_length_filter_reads - after_mt_filter_reads))
    echo "  After mtRNA filter reads: $after_mt_filter_reads (mtRNA removed: $mt_removed_reads)"
    no_contam_file="$CLEAN_DATA_DIR/${sample}_no_contam.fastq"
    rrna_trna_contam_sam="$CLEAN_DATA_DIR/${sample}_rrna_trna_contaminant_mapped.sam"
    echo "  Running single-end rRNA+tRNA filtering..."
    bowtie2 -x "$rRNA_tRNA_CONTAMINANT_INDEX" -U "$nomt_contam_file" --un "$no_contam_file" --no-unal -L 20 -p $THREADS_BOWTIE2 -S "$rrna_trna_contam_sam" 2>/dev/null
    rm -f $nomt_contam_file
    if [ ! -f "$no_contam_file" ] || [ $(stat -c %s "$no_contam_file") -eq 0 ]; then
        echo "  ERROR: empty after rRNA+tRNA filter, mark sample failed"
        echo "$sample,single,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,0,0,,,0,0" >> "$STATS_FILE"
        echo "$sample,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
        rm -f "$rrna_trna_contam_sam" "$no_contam_file"
        return 1
    fi
    after_contaminant_filter_reads=$(cat "$no_contam_file" | wc -l | awk '{print int($1/4)}')
    rrna_trna_removed_reads=$((after_mt_filter_reads - after_contaminant_filter_reads))
    echo "  After rRNA+tRNA filter reads: $after_contaminant_filter_reads (rRNA+tRNA removed: $rrna_trna_removed_reads)"
    echo "  Total after decontamination: $after_contaminant_filter_reads"
    contaminant_reads=$((mt_removed_reads + rrna_trna_removed_reads))
    rRNA_reads=$rrna_trna_removed_reads
    mtRNA_reads=$mt_removed_reads
    tRNA_reads=0
    if [ $after_length_filter_reads -gt 0 ]; then
        contaminant_removal_rate=$(echo "scale=4; $contaminant_reads / $after_length_filter_reads * 100" | bc)
        rRNA_removal_rate=$(echo "scale=4; $rRNA_reads / $after_length_filter_reads * 100" | bc)
        mtRNA_removal_rate=$(echo "scale=4; $mtRNA_reads / $after_length_filter_reads * 100" | bc)
    else
        contaminant_removal_rate=0
        rRNA_removal_rate=0
        mtRNA_removal_rate=0
    fi
    echo "  Contaminant summary: total=$contaminant_reads (mtRNA=$mtRNA_reads, rRNA+tRNA=$rRNA_reads)"
    echo "  Removal rate: $contaminant_removal_rate% (mtRNA: $mtRNA_removal_rate%, rRNA+tRNA: $rRNA_removal_rate%)"
    final_file="$CLEAN_DATA_DIR/${sample}_final.fastq.gz"
    temp_length_file="$CLEAN_DATA_DIR/${sample}_lengths.txt"
    gzip -c "$no_contam_file" > "$final_file"
    rm -f $no_contam_file
    if [ ! -f "$final_file" ] || [ $(stat -c %s "$final_file") -eq 0 ]; then
        echo "  WARNING: final file empty, proceed with stats"
        final_reads=0
        min_length=0
        max_length=0
        avg_length=0
        touch "$final_file"
    else
        final_reads=$(zcat "$final_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
        zcat "$final_file" 2>/dev/null | awk 'NR%4==2 {print length($0)}' > "$temp_length_file"
        min_length=$(sort -n "$temp_length_file" | head -n1 2>/dev/null || echo 0)
        max_length=$(sort -n "$temp_length_file" | tail -n1 2>/dev/null || echo 0)
        avg_length=$(awk '{sum+=$1} END {if(NR>0) print sum/NR; else print 0}' "$temp_length_file")
        echo "  Final clean reads: $final_reads (${min_length}-${max_length}nt, mean ${avg_length}nt)"
    fi
    echo "  Generating FastQC report..."
    fastqc -o "$FASTQC_CLEAN_DIR" -t $THREADS_FASTQC "$final_file" > /dev/null 2>&1
    if [ -f "$FASTQC_CLEAN_DIR/${sample}_final_fastqc.html" ]; then
        echo "  FastQC report generated"
    else
        echo "  WARNING: FastQC report failed (empty file)"
    fi
    echo "$sample,single,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,$after_contaminant_filter_reads,$final_reads,$min_length,$max_length,$avg_length,$rRNA_removal_rate,$mtRNA_removal_rate" >> "$STATS_FILE"
    echo "$sample,$contaminant_reads,$rRNA_reads,$mtRNA_reads,$tRNA_reads,$contaminant_removal_rate" >> "$CONTAM_STATS_FILE"
    rm -f "$temp_length_file" "$mtcontam_sam" "$rrna_trna_contam_sam"
    echo "  Single-end sample finished"
    return 0
}
process_paired_end() {
    local sample1=$1
    local sample2=$2
    local adapter3_1=$3
    local adapter5_1=$4
    local adapter3_2=$5
    local adapter5_2=$6
    local base_name=${sample1%_1}
    local raw_reads=0 after_trim_qc_reads=0 after_length_filter_reads=0 after_contaminant_filter_reads=0 final_reads=0
    local min_length=0 max_length=0 avg_length=0
    local raw_reads1=0 raw_reads2=0
    local after_trim_qc_reads1=0 after_trim_qc_reads2=0
    local after_length_filter_reads1=0 after_length_filter_reads2=0
    local after_contaminant_filter_reads1=0 after_contaminant_filter_reads2=0 final_reads1=0 final_reads2=0
    local contaminant_reads=0 rRNA_reads=0 mtRNA_reads=0 tRNA_reads=0 contaminant_removal_rate=0 rRNA_removal_rate=0 mtRNA_removal_rate=0
    echo -e "\n[$current/$total] Processing paired-end sample: $base_name (${sample1}/${sample2})"
    raw_file1="$RAW_DATA_DIR/$sample1.fastq.gz"
    raw_file2="$RAW_DATA_DIR/$sample2.fastq.gz"
    if [ ! -f "$raw_file1" ] || [ ! -f "$raw_file2" ]; then
        echo "  ERROR: raw file missing: $raw_file1 or $raw_file2"
        echo "$base_name,paired,0,0,0,0,0,,,0,0" >> "$STATS_FILE"
        echo "$base_name,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
        return 1
    fi
    raw_reads1=$(zcat "$raw_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    raw_reads2=$(zcat "$raw_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    raw_reads=$((raw_reads1 + raw_reads2))
    echo "  Raw reads: $raw_reads (R1: $raw_reads1, R2: $raw_reads2)"
    qc_trimmed_file1="$CLEAN_DATA_DIR/${sample1}_qc_trimmed.fastq.gz"
    qc_trimmed_file2="$CLEAN_DATA_DIR/${sample2}_qc_trimmed.fastq.gz"
    length_filtered_file1="$CLEAN_DATA_DIR/${sample1}_length_filtered.fastq.gz"
    length_filtered_file2="$CLEAN_DATA_DIR/${sample2}_length_filtered.fastq.gz"
    trim_ok=0
    length_ok=0
    if [ -n "$adapter3_1" ] || [ -n "$adapter5_1" ] || [ -n "$adapter3_2" ] || [ -n "$adapter5_2" ]; then
        echo "  Adapter detected, first step: cutadapt paired adapter trimming + QC..."
        cutadapt_args=()
        [ -n "$adapter3_1" ] && cutadapt_args+=("-a" "$adapter3_1")
        [ -n "$adapter5_1" ] && cutadapt_args+=("-g" "$adapter5_1")
        [ -n "$adapter3_2" ] && cutadapt_args+=("-A" "$adapter3_2")
        [ -n "$adapter5_2" ] && cutadapt_args+=("-G" "$adapter5_2")
        cutadapt "${cutadapt_args[@]}" -j $THREADS_CUTADAPT --trim-n -q 20 --trimmed-only -o "$qc_trimmed_file1" -p "$qc_trimmed_file2" "$raw_file1" "$raw_file2" 2>/dev/null
        if [ -f "$qc_trimmed_file1" ] && [ -f "$qc_trimmed_file2" ] && [ $(stat -c %s "$qc_trimmed_file1") -gt 0 ] && [ $(stat -c %s "$qc_trimmed_file2") -gt 0 ]; then
            after_trim_qc_reads1=$(zcat "$qc_trimmed_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_qc_reads2=$(zcat "$qc_trimmed_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_qc_reads=$((after_trim_qc_reads1 + after_trim_qc_reads2))
            echo "  After paired trim & QC reads: $after_trim_qc_reads (R1: $after_trim_qc_reads1, R2: $after_trim_qc_reads2)"
            echo "  Second step: cutadapt paired length filter 36-150nt..."
            cutadapt --minimum-length 36 --maximum-length 150 --paired -j $THREADS_CUTADAPT -o "$length_filtered_file1" -p "$length_filtered_file2" "$qc_trimmed_file1" "$qc_trimmed_file2" 2>/dev/null
            if [ -f "$length_filtered_file1" ] && [ -f "$length_filtered_file2" ] && [ $(stat -c %s "$length_filtered_file1") -gt 0 ] && [ $(stat -c %s "$length_filtered_file2") -gt 0 ]; then
                after_length_filter_reads1=$(zcat "$length_filtered_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_filter_reads2=$(zcat "$length_filtered_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_filter_reads=$((after_length_filter_reads1 + after_length_filter_reads2))
                echo "  After paired QC + length filter reads: $after_length_filter_reads (R1: $after_length_filter_reads1, R2: $after_length_filter_reads2) (36-150nt)"
                length_ok=1
                rm -f $qc_trimmed_file1 $qc_trimmed_file2
            fi
        fi
    else
        echo "  No adapter, first step: trim_galore paired QC..."
        trim_galore --quality 20 --trim-n --paired --cores $THREADS_TRIMGALORE -o "$CLEAN_DATA_DIR" "$raw_file1" "$raw_file2" 2>/dev/null
        qc_trimmed_tmp1="${CLEAN_DATA_DIR}/${sample1}_val_1.fq.gz"
        qc_trimmed_tmp2="${CLEAN_DATA_DIR}/${sample2}_val_2.fq.gz"
        [ ! -f "$qc_trimmed_tmp1" ] && qc_trimmed_tmp1="${CLEAN_DATA_DIR}/${base_name}_1_trimmed.fq.gz"
        [ ! -f "$qc_trimmed_tmp2" ] && qc_trimmed_tmp2="${CLEAN_DATA_DIR}/${base_name}_2_trimmed.fq.gz"
        if [ -f "$qc_trimmed_tmp1" ] && [ -f "$qc_trimmed_tmp2" ]; then
            mv "$qc_trimmed_tmp1" "$qc_trimmed_file1"
            mv "$qc_trimmed_tmp2" "$qc_trimmed_file2"
            after_trim_qc_reads1=$(zcat "$qc_trimmed_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_qc_reads2=$(zcat "$qc_trimmed_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_qc_reads=$((after_trim_qc_reads1 + after_trim_qc_reads2))
            echo "  After paired QC reads: $after_trim_qc_reads (R1: $after_trim_qc_reads1, R2: $after_trim_qc_reads2)"
            echo "  Second step: cutadapt paired length filter 36-150nt..."
            cutadapt --minimum-length 36 --maximum-length 150 -j $THREADS_CUTADAPT -o "$length_filtered_file1" "$qc_trimmed_file1" 2>/dev/null
            cutadapt --minimum-length 36 --maximum-length 150 -j $THREADS_CUTADAPT -o "$length_filtered_file2" "$qc_trimmed_file2" 2>/dev/null
            if [ -f "$length_filtered_file1" ] && [ -f "$length_filtered_file2" ] && [ $(stat -c %s "$length_filtered_file1") -gt 0 ] && [ $(stat -c %s "$length_filtered_file2") -gt 0 ]; then
                after_length_filter_reads1=$(zcat "$length_filtered_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_filter_reads2=$(zcat "$length_filtered_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_filter_reads=$((after_length_filter_reads1 + after_length_filter_reads2))
                echo "  After paired QC + length filter reads: $after_length_filter_reads (R1: $after_length_filter_reads1, R2: $after_length_filter_reads2) (36-150nt)"
                length_ok=1
                rm -f $qc_trimmed_file1 $qc_trimmed_file2
            fi
        fi
    fi
    if [ $length_ok -eq 0 ]; then
        echo "  ERROR: empty after length filter, mark sample failed"
        echo "$base_name,paired,$raw_reads,$after_trim_qc_reads,0,-1,0,,,0,0" >> "$STATS_FILE"
        echo "$base_name,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
        rm -f "$qc_trimmed_file1" "$qc_trimmed_file2" "$length_filtered_file1" "$length_filtered_file2"
        return 1
    fi
    local diff_ratio=0
    local keep_paired="true"
    local single_sample=""
    local single_length_file=""
    if [ $after_length_filter_reads1 -eq 0 ] || [ $after_length_filter_reads2 -eq 0 ]; then
        echo "  WARNING: one of R1/R2 empty, switch to single-end decontamination"
        keep_paired="false"
        if [ $after_length_filter_reads1 -gt $after_length_filter_reads2 ]; then
            single_sample="$sample1"
            single_length_file="$length_filtered_file1"
        else
            single_sample="$sample2"
            single_length_file="$length_filtered_file2"
        fi
    else
        if [ $after_length_filter_reads1 -gt $after_length_filter_reads2 ]; then
            diff_ratio=$(echo "scale=2; $after_length_filter_reads1 / $after_length_filter_reads2" | bc)
        else
            diff_ratio=$(echo "scale=2; $after_length_filter_reads2 / $after_length_filter_reads1" | bc)
        fi
        if (( $(echo "$diff_ratio > 5" | bc -l) )); then
            echo "  WARNING: R1/R2 ratio=$diff_ratio (>5), switch to single-end decontamination"
            keep_paired="false"
            if [ $after_length_filter_reads1 -gt $after_length_filter_reads2 ]; then
                single_sample="$sample1"
                single_length_file="$length_filtered_file1"
            else
                single_sample="$sample2"
                single_length_file="$length_filtered_file2"
            fi
        else
            echo "  R1/R2 ratio=$diff_ratio (<=5), keep paired-end decontamination"
            keep_paired="true"
        fi
    fi
    nomt_contam_file1="$CLEAN_DATA_DIR/${sample1}_nomt_contam.fastq"
    nomt_contam_file2="$CLEAN_DATA_DIR/${sample2}_nomt_contam.fastq"
    mtcontam_sam="$CLEAN_DATA_DIR/${base_name}_mtcontaminant_mapped.sam"
    no_contam_file1="$CLEAN_DATA_DIR/${sample1}_no_contam.fastq"
    no_contam_file2="$CLEAN_DATA_DIR/${sample2}_no_contam.fastq"
    rrna_trna_contam_sam="$CLEAN_DATA_DIR/${base_name}_rrna_trna_contaminant_mapped.sam"
    if [ "$keep_paired" == "true" ]; then
        echo "  Running paired-end mtRNA filtering..."
        bowtie2 -x "$mtCONTAMINANT_INDEX" -1 "$length_filtered_file1" -2 "$length_filtered_file2" \
            --un-conc "$CLEAN_DATA_DIR/${base_name}_nomt_contam_.fastq" --no-unal -L 20 -p $THREADS_BOWTIE2 -S "$mtcontam_sam" 2>/dev/null
        rm -f $length_filtered_file1 $length_filtered_file2
        if [ -f "$CLEAN_DATA_DIR/${base_name}_nomt_contam_.1.fastq" ]; then
            mv "$CLEAN_DATA_DIR/${base_name}_nomt_contam_.1.fastq" "$nomt_contam_file1"
        else
            echo "  ERROR: R1 missing after paired mtRNA filter"
            echo "$base_name,paired,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,0,0,,,0,0" >> "$STATS_FILE"
            echo "$base_name,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
            rm -f "$mtcontam_sam" "$nomt_contam_file1" "$nomt_contam_file2"
            return 1
        fi
        if [ -f "$CLEAN_DATA_DIR/${base_name}_nomt_contam_.2.fastq" ]; then
            mv "$CLEAN_DATA_DIR/${base_name}_nomt_contam_.2.fastq" "$nomt_contam_file2"
        else
            echo "  ERROR: R2 missing after paired mtRNA filter"
            echo "$base_name,paired,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,0,0,,,0,0" >> "$STATS_FILE"
            echo "$base_name,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
            rm -f "$mtcontam_sam" "$nomt_contam_file1" "$nomt_contam_file2"
            return 1
        fi
        after_mt_filter_reads1=$(cat "$nomt_contam_file1" | wc -l | awk '{print int($1/4)}')
        after_mt_filter_reads2=$(cat "$nomt_contam_file2" | wc -l | awk '{print int($1/4)}')
        after_mt_filter_reads=$((after_mt_filter_reads1 + after_mt_filter_reads2))
        mt_removed_reads=$((after_length_filter_reads - after_mt_filter_reads))
        echo "  After mtRNA filter reads: $after_mt_filter_reads (mtRNA removed: $mt_removed_reads)"
        echo "  Running paired-end rRNA+tRNA filtering..."
        bowtie2 -x "$rRNA_tRNA_CONTAMINANT_INDEX" -1 "$nomt_contam_file1" -2 "$nomt_contam_file2" \
            --un-conc "$CLEAN_DATA_DIR/${base_name}_no_contam_.fastq" --no-unal -L 20 -p $THREADS_BOWTIE2 -S "$rrna_trna_contam_sam" 2>/dev/null
        rm -f $nomt_contam_file1 $nomt_contam_file2
        if [ -f "$CLEAN_DATA_DIR/${base_name}_no_contam_.1.fastq" ]; then
            mv "$CLEAN_DATA_DIR/${base_name}_no_contam_.1.fastq" "$no_contam_file1"
        else
            echo "  ERROR: R1 missing after paired rRNA+tRNA filter"
            echo "$base_name,paired,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,0,0,,,0,0" >> "$STATS_FILE"
            echo "$base_name,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
            rm -f "$rrna_trna_contam_sam" "$no_contam_file1" "$no_contam_file2"
            return 1
        fi
        if [ -f "$CLEAN_DATA_DIR/${base_name}_no_contam_.2.fastq" ]; then
            mv "$CLEAN_DATA_DIR/${base_name}_no_contam_.2.fastq" "$no_contam_file2"
        else
            echo "  ERROR: R2 missing after paired rRNA+tRNA filter"
            echo "$base_name,paired,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,0,0,,,0,0" >> "$STATS_FILE"
            echo "$base_name,0,0,0,0,0" >> "$CONTAM_STATS_FILE"
            rm -f "$rrna_trna_contam_sam" "$no_contam_file1" "$no_contam_file2"
            return 1
        fi
        after_contaminant_filter_reads1=$(cat "$no_contam_file1" | wc -l | awk '{print int($1/4)}')
        after_contaminant_filter_reads2=$(cat "$no_contam_file2" | wc -l | awk '{print int($1/4)}')
        after_contaminant_filter_reads=$((after_contaminant_filter_reads1 + after_contaminant_filter_reads2))
        rrna_trna_removed_reads=$((after_mt_filter_reads - after_contaminant_filter_reads))
        echo "  After rRNA+tRNA filter reads: $after_contaminant_filter_reads (rRNA+tRNA removed: $rrna_trna_removed_reads)"
        echo "  Total after decontamination: $after_contaminant_filter_reads (R1: $after_contaminant_filter_reads1, R2: $after_contaminant_filter_reads2)"
    else
        echo "  Running single-end mtRNA filtering (keep $single_sample)..."
        single_nomt_contam_file="$CLEAN_DATA_DIR/${single_sample}_nomt_contam.fastq"
        bowtie2 -x "$mtCONTAMINANT_INDEX" -U "$single_length_file" --un "$single_nomt_contam_file" --no-unal -L 20 -p $THREADS_BOWTIE2 -S "$mtcontam_sam" 2>/dev/null
        rm -f $length_filtered_file1 $length_filtered_file2
        after_mt_filter_reads_single=$(cat "$single_nomt_contam_file" | wc -l | awk '{print int($1/4)}')
        mt_removed_reads=$((after_length_filter_reads - after_mt_filter_reads_single))
        echo "  After mtRNA filter reads: $after_mt_filter_reads_single (mtRNA removed: $mt_removed_reads)"
        echo "  Running single-end rRNA+tRNA filtering (keep $single_sample)..."
        single_no_contam_file="$CLEAN_DATA_DIR/${single_sample}_no_contam.fastq"
        bowtie2 -x "$rRNA_tRNA_CONTAMINANT_INDEX" -U "$single_nomt_contam_file" --un "$single_no_contam_file" --no-unal -L 20 -p $THREADS_BOWTIE2 -S "$rrna_trna_contam_sam" 2>/dev/null
        rm -f $single_nomt_contam_file
        if [ "$single_sample" == "$sample1" ]; then
            after_contaminant_filter_reads1=$(cat "$single_no_contam_file" | wc -l | awk '{print int($1/4)}')
            after_contaminant_filter_reads2=0
            touch "$no_contam_file2"
            mv "$single_no_contam_file" "$no_contam_file1"
        else
            after_contaminant_filter_reads1=0
            after_contaminant_filter_reads2=$(cat "$single_no_contam_file" | wc -l | awk '{print int($1/4)}')
            touch "$no_contam_file1"
            mv "$single_no_contam_file" "$no_contam_file2"
        fi
        after_contaminant_filter_reads=$((after_contaminant_filter_reads1 + after_contaminant_filter_reads2))
        rrna_trna_removed_reads=$((after_mt_filter_reads_single - after_contaminant_filter_reads))
        echo "  After rRNA+tRNA filter reads: $after_contaminant_filter_reads (rRNA+tRNA removed: $rrna_trna_removed_reads)"
        echo "  Total after decontamination: $after_contaminant_filter_reads (only $single_sample retained: $after_contaminant_filter_reads)"
    fi
    contaminant_reads=$((mt_removed_reads + rrna_trna_removed_reads))
    rRNA_reads=$rrna_trna_removed_reads
    mtRNA_reads=$mt_removed_reads
    tRNA_reads=0
    if [ $after_length_filter_reads -gt 0 ]; then
        contaminant_removal_rate=$(echo "scale=4; $contaminant_reads / $after_length_filter_reads * 100" | bc)
        rRNA_removal_rate=$(echo "scale=4; $rRNA_reads / $after_length_filter_reads * 100" | bc)
        mtRNA_removal_rate=$(echo "scale=4; $mtRNA_reads / $after_length_filter_reads * 100" | bc)
    else
        contaminant_removal_rate=0
        rRNA_removal_rate=0
        mtRNA_removal_rate=0
    fi
    echo "  Contaminant summary: total=$contaminant_reads (mtRNA=$mtRNA_reads, rRNA+tRNA=$rRNA_reads)"
    echo "  Removal rate: $contaminant_removal_rate% (mtRNA: $mtRNA_removal_rate%, rRNA+tRNA: $rRNA_removal_rate%)"
    if [ $after_contaminant_filter_reads -eq 0 ]; then
        echo "  ERROR: empty after decontamination, mark sample failed"
        echo "$base_name,paired,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,0,0,,,0,0" >> "$STATS_FILE"
        echo "$base_name,$contaminant_reads,$rRNA_reads,$mtRNA_reads,$tRNA_reads,$contaminant_removal_rate" >> "$CONTAM_STATS_FILE"
        rm -f "$mtcontam_sam" "$rrna_trna_contam_sam" "$no_contam_file1" "$no_contam_file2"
        return 1
    fi
    final_file1="$CLEAN_DATA_DIR/${sample1}_final.fastq.gz"
    final_file2="$CLEAN_DATA_DIR/${sample2}_final.fastq.gz"
    temp_length_file="$CLEAN_DATA_DIR/${base_name}_lengths.txt"
    if [ -f "$no_contam_file1" ] && [ $(stat -c %s "$no_contam_file1") -gt 0 ]; then
        gzip -c "$no_contam_file1" > "$final_file1"
        final_reads1=$(zcat "$final_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    else
        touch "$final_file1"
        final_reads1=0
        echo "  WARNING: R1 final file empty"
    fi
    if [ -f "$no_contam_file2" ] && [ $(stat -c %s "$no_contam_file2") -gt 0 ]; then
        gzip -c "$no_contam_file2" > "$final_file2"
        final_reads2=$(zcat "$final_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    else
        touch "$final_file2"
        final_reads2=0
        echo "  WARNING: R2 final file empty"
    fi
    final_reads=$((final_reads1 + final_reads2))
    rm -f $no_contam_file1 $no_contam_file2
    if [ $final_reads -gt 0 ]; then
        cat <(zcat "$final_file1" 2>/dev/null) <(zcat "$final_file2" 2>/dev/null) | awk 'NR%4==2 {print length($0)}' > "$temp_length_file"
        min_length=$(sort -n "$temp_length_file" | head -n1 2>/dev/null || echo 0)
        max_length=$(sort -n "$temp_length_file" | tail -n1 2>/dev/null || echo 0)
        avg_length=$(awk '{sum+=$1} END {if(NR>0) print sum/NR; else print 0}' "$temp_length_file")
    fi
    echo "  Final clean reads: $final_reads (R1: $final_reads1, R2: $final_reads2) (${min_length}-${max_length}nt, mean ${avg_length}nt)"
    echo "  Generating paired FastQC reports..."
    fastqc -o "$FASTQC_CLEAN_DIR" -t $THREADS_FASTQC "$final_file1" "$final_file2" > /dev/null 2>&1
    if [ -f "$FASTQC_CLEAN_DIR/${sample1}_final_fastqc.html" ] && [ -f "$FASTQC_CLEAN_DIR/${sample2}_final_fastqc.html" ]; then
        echo "  Paired FastQC reports generated"
    else
        echo "  WARNING: FastQC reports incomplete (some empty files)"
    fi
    echo "$base_name,paired,$raw_reads,$after_trim_qc_reads,$after_length_filter_reads,$after_contaminant_filter_reads,$final_reads,$min_length,$max_length,$avg_length,$rRNA_removal_rate,$mtRNA_removal_rate" >> "$STATS_FILE"
    echo "$base_name,$contaminant_reads,$rRNA_reads,$mtRNA_reads,$tRNA_reads,$contaminant_removal_rate" >> "$CONTAM_STATS_FILE"
    rm -f "$temp_length_file" "$mtcontam_sam" "$rrna_trna_contam_sam"
    echo "  Paired-end sample finished"
    return 0
}
declare -A paired_groups
declare -A processed_paired
for sample in "${sample_list[@]}"; do
    if [[ "${is_paired_end[$sample]}" == "true" ]] && [[ -z "${processed_paired[$sample]}" ]]; then
        if [[ "$sample" =~ _1$ ]]; then
            base_name=${sample%_1}
            sample2="${base_name}_2"
            if [[ " ${sample_list[@]} " =~ " $sample2 " ]] && [[ "${is_paired_end[$sample2]}" == "true" ]] && [[ -z "${processed_paired[$sample2]}" ]]; then
                adapter3_1=$(echo "${sample_info[$sample]}" | cut -d',' -f1)
                adapter5_1=$(echo "${sample_info[$sample]}" | cut -d',' -f2)
                adapter3_2=$(echo "${sample_info[$sample2]}" | cut -d',' -f1)
                adapter5_2=$(echo "${sample_info[$sample2]}" | cut -d',' -f2)
                paired_groups["$base_name"]="$sample,$sample2,$adapter3_1,$adapter5_1,$adapter3_2,$adapter5_2"
                processed_paired[$sample]="true"
                processed_paired[$sample2]="true"
                echo "  Paired detected: $sample <-> $sample2"
            else
                echo "  WARNING: $sample missing _2 pair, downgrade to single-end"
                is_paired_end["$sample"]="false"
            fi
        elif [[ "$sample" =~ _2$ ]]; then
            echo "  WARNING: $sample missing _1 pair, downgrade to single-end"
            is_paired_end["$sample"]="false"
        fi
    fi
done
single_samples=()
paired_samples=()
for sample in "${sample_list[@]}"; do
    if [[ "${is_paired_end[$sample]}" == "true" ]]; then
        paired_samples+=("$sample")
    else
        single_samples+=("$sample")
    fi
done
echo "  After pair validation - single-end: ${#single_samples[@]} | paired-end: ${#paired_samples[@]}"
for sample in "${single_samples[@]}"; do
    if [[ -z "${processed_paired[$sample]}" ]]; then
        current=$((current + 1))
        adapter3=$(echo "${sample_info[$sample]}" | cut -d',' -f1)
        adapter5=$(echo "${sample_info[$sample]}" | cut -d',' -f2)
        if process_single_end "$sample" "$adapter3" "$adapter5"; then
            success=$((success + 1))
        else
            failed=$((failed + 1))
        fi
    fi
done
for base_name in "${!paired_groups[@]}"; do
    current=$((current + 1))
    IFS=',' read -r sample1 sample2 adapter3_1 adapter5_1 adapter3_2 adapter5_2 <<< "${paired_groups[$base_name]}"
    if process_paired_end "$sample1" "$sample2" "$adapter3_1" "$adapter5_1" "$adapter3_2" "$adapter5_2"; then
        success=$((success + 1))
    else
        failed=$((failed + 1))
    fi
done
echo -e "\n===== RNA-seq filter pipeline finished ====="
echo "End time: $(date +'%Y-%m-%d %H:%M:%S')"
echo "Total processed samples: $current | Success: $success | Failed: $failed"
echo "Main statistics: $STATS_FILE"
echo "Contaminant removal statistics: $CONTAM_STATS_FILE"
echo "Output file description:"
echo "  - Final clean reads: ${CLEAN_DATA_DIR}/*_final.fastq.gz"
echo "  - FastQC reports: $FASTQC_CLEAN_DIR"
echo "Stat definition:"
echo "  - after_trim_qc_reads: reads after adapter trimming and QC"
echo "  - after_length_filter_reads: reads after QC + 36-150nt length filter"
echo "  - rRNA_removal_rate/mtRNA_removal_rate: rRNA/mtRNA removal rate (%)"
echo "=========================="
