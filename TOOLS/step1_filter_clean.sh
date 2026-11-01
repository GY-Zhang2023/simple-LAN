#!/bin/bash
# step1 preprocessing script
# default resource: 2 cores + 8 threads, 50G memory
# only human and mouse supported
# all comments removed except this header
# unified English output
# adjustable thread params at top

THREADS_BOWTIE2=8
THREADS_FASTQC=8
THREADS_CUTADAPT=8
THREADS_TRIMGALORE=8

FILTER_MONOSOME_LENGTH=true

declare -A SPECIES_REF_MAP
SPECIES_REF_MAP["human"]="$HOME/ribo_references/rrna_trna/human_rrna_tRNA_index"
SPECIES_REF_MAP["mouse"]="$HOME/reference/mouse_transcriptome/mouse_rrna_tRNA_index"

if [ $# -ne 3 ]; then
    echo "Usage: $0 <CSV_ABS_PATH> <MIN_LENGTH> <MAX_LENGTH>"
    echo "Example: $0 /home/data/WTC-11.csv 21 40"
    exit 1
fi

CSV_ABS_PATH="$1"
MIN_LENGTH="$2"
MAX_LENGTH="$3"

if ! [[ "$MIN_LENGTH" =~ ^[0-9]+$ && "$MAX_LENGTH" =~ ^[0-9]+$ && "$MAX_LENGTH" -ge "$MIN_LENGTH" ]]; then
    echo "ERROR: invalid length parameters, must be positive integer and MAX >= MIN"
    exit 1
fi

if [ "$MIN_LENGTH" -eq 21 ] && [ "$MAX_LENGTH" -eq 40 ]; then
    SELECTED_ANALYSIS="Monosome"
elif [ "$MIN_LENGTH" -eq 21 ] && [ "$MAX_LENGTH" -eq 80 ]; then
    SELECTED_ANALYSIS="Disome"
else
    SELECTED_ANALYSIS="Custom"
fi

WORK_DIR=$(pwd)
CSV="$CSV_ABS_PATH"

if [ ! -f "${CSV_ABS_PATH}" ]; then
    echo "ERROR: CSV file not found: ${CSV_ABS_PATH}"
    exit 1
fi
if [ ! -s "${CSV_ABS_PATH}" ]; then
    echo "ERROR: CSV file is empty: ${CSV_ABS_PATH}"
    exit 1
fi

echo -e "\n==================== Species Selection ===================="
echo "Select species:"
echo "1 human"
echo "2 mouse"
read -p "Input choice (1/2): " SPECIES_CHOICE

case ${SPECIES_CHOICE} in
    1)
        SELECTED_SPECIES="human"
        ;;
    2)
        SELECTED_SPECIES="mouse"
        ;;
    *)
        echo "ERROR: invalid choice, only 1(human) or 2(mouse)"
        exit 1
        ;;
esac

RRNA_TRNA_INDEX=${SPECIES_REF_MAP[${SELECTED_SPECIES}]}

echo -e "\nConfirmed species: ${SELECTED_SPECIES}"
echo "rRNA-tRNA index path: ${RRNA_TRNA_INDEX}"
read -p "Continue? (y/n): " CONFIRM
if [ "${CONFIRM}" != "y" ] && [ "${CONFIRM}" != "Y" ]; then
    echo "User aborted"
    exit 0
fi

RAW_DATA_DIR="${WORK_DIR}/1.raw_fastq"
CLEAN_DATA_DIR="${WORK_DIR}/2.clean_fastq"
STATS_FILE="${CLEAN_DATA_DIR}/filter_stats.csv"
FASTQC_CLEAN_DIR="${CLEAN_DATA_DIR}/fastqc_clean"

echo -e "\n==================== Run Confirmation ===================="
echo "CSV path: ${CSV_ABS_PATH}"
echo "Raw fastq dir: ${RAW_DATA_DIR}"
echo "Species: ${SELECTED_SPECIES}"
echo "rRNA-tRNA index: ${RRNA_TRNA_INDEX}"
echo "Analysis type: ${SELECTED_ANALYSIS}"
echo "Length filter: ${MIN_LENGTH}-${MAX_LENGTH} nt"

mkdir -p "$CLEAN_DATA_DIR" "$FASTQC_CLEAN_DIR" || { echo "mkdir failed"; exit 1; }

echo "sample_name,sequencing_type,analysis_type,length_filter_rule,raw_reads,after_trim_adapter_reads,after_length_constraint_reads,after_rtrna_filter_reads,final_reads,min_length,max_length,avg_length" > "$STATS_FILE"

check_dependency() {
    if [ ! -f "$1.1.bt2" ]; then echo "ERROR: bowtie2 index missing: $1"; exit 1; fi
}
check_dependency "$RRNA_TRNA_INDEX"
for tool in bowtie2 cutadapt fastqc trim_galore; do
    if ! command -v $tool &> /dev/null; then echo "ERROR: tool not found: $tool"; exit 1; fi
done

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

echo "===== Preprocessing Start ====="
echo "Start time: $(date +'%Y-%m-%d %H:%M:%S')"
echo "Total samples: $total | single-end: ${#single_samples[@]} | paired-end: ${#paired_samples[@]}"
echo "Output dir: $CLEAN_DATA_DIR"
echo "Rule list:"
echo "  1. trim_galore / cutadapt two-step workflow: adapter removal first, then length filter"
echo "  2. separate stats: after adapter trimming / after length constraint"
echo "  3. if R1/R2 read count ratio >3 or one empty: single-end decontamination"
echo "  4. if R1/R2 read count ratio <=3: paired-end decontamination"
echo "  5. empty file after length filter: warning only, pipeline continue"
echo "  6. analysis type: ${SELECTED_ANALYSIS}, length parameter passed from master script"
echo "=========================="

current=0
success=0
failed=0

process_single_end() {
    local sample=$1
    local adapter3=$2
    local adapter5=$3
    local raw_reads=0 after_trim_adapter_reads=0 after_length_constraint_reads=0 after_rtrna_filter_reads=0 final_reads=0
    local min_length=0 max_length=0 avg_length=0
    echo -e "\n[$current/$total] Process single-end sample: $sample"
    raw_file="$RAW_DATA_DIR/$sample.fastq.gz"
    if [ ! -f "$raw_file" ]; then
        echo "  ERROR: raw file missing: $raw_file"
        echo "$sample,single,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},0,0,0,0,0,,,," >> "$STATS_FILE"
        return 1
    fi
    raw_reads=$(zcat "$raw_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    echo "  Raw reads: $raw_reads"
    adapter_trimmed_file="$CLEAN_DATA_DIR/${sample}_adapter_trimmed.fastq.gz"
    length_filtered_file="$CLEAN_DATA_DIR/${sample}_length_filtered.fastq.gz"
    length_ok=0
    if [ -n "$adapter3" ] || [ -n "$adapter5" ]; then
        echo "  Adapters detected, step1: cutadapt adapter trimming only"
        cutadapt_args=()
        [ -n "$adapter3" ] && cutadapt_args+=("-a" "$adapter3")
        [ -n "$adapter5" ] && cutadapt_args+=("-g" "$adapter5")
        cutadapt "${cutadapt_args[@]}" -j $THREADS_CUTADAPT --trim-n -q 28 --trimmed-only -o "$adapter_trimmed_file" "$raw_file" 2>/dev/null
        if [ -f "$adapter_trimmed_file" ] && [ $(stat -c %s "$adapter_trimmed_file") -gt 0 ]; then
            after_trim_adapter_reads=$(zcat "$adapter_trimmed_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            echo "  Reads after adapter trim: $after_trim_adapter_reads"
            echo "  Step2: cutadapt length filter ${MIN_LENGTH}-${MAX_LENGTH}nt"
            if [[ "${SELECTED_ANALYSIS}" != "Monosome" || "${FILTER_MONOSOME_LENGTH}" == "true" ]]; then
                cutadapt --minimum-length ${MIN_LENGTH} --maximum-length ${MAX_LENGTH} -j $THREADS_CUTADAPT -o "$length_filtered_file" "$adapter_trimmed_file" 2>/dev/null
            else
                cp "$adapter_trimmed_file" "$length_filtered_file"
            fi
            if [ -f "$length_filtered_file" ] && [ $(stat -c %s "$length_filtered_file") -gt 0 ]; then
                after_length_constraint_reads=$(zcat "$length_filtered_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                echo "  Reads after length filter: $after_length_constraint_reads (${MIN_LENGTH}-${MAX_LENGTH}nt)"
                length_ok=1
            fi
        fi
    else
        echo "  No adapter, step1: trim_galore quality trim only"
        trim_galore --quality 28 --trim-n --cores $THREADS_TRIMGALORE -o "$CLEAN_DATA_DIR" "$raw_file" 2>/dev/null
        adapter_trimmed_tmp="${CLEAN_DATA_DIR}/${sample}_trimmed.fq.gz"
        [ ! -f "$adapter_trimmed_tmp" ] && adapter_trimmed_tmp="${CLEAN_DATA_DIR}/${sample}_val.fq.gz"
        if [ -f "$adapter_trimmed_tmp" ]; then
            mv "$adapter_trimmed_tmp" "$adapter_trimmed_file"
            after_trim_adapter_reads=$(zcat "$adapter_trimmed_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            echo "  Reads after quality trim: $after_trim_adapter_reads"
            echo "  Step2: cutadapt length filter ${MIN_LENGTH}-${MAX_LENGTH}nt"
            if [[ "${SELECTED_ANALYSIS}" != "Monosome" || "${FILTER_MONOSOME_LENGTH}" == "true" ]]; then
                cutadapt --minimum-length ${MIN_LENGTH} --maximum-length ${MAX_LENGTH} -j $THREADS_CUTADAPT -o "$length_filtered_file" "$adapter_trimmed_file" 2>/dev/null
            else
                cp "$adapter_trimmed_file" "$length_filtered_file"
            fi
            if [ -f "$length_filtered_file" ] && [ $(stat -c %s "$length_filtered_file") -gt 0 ]; then
                after_length_constraint_reads=$(zcat "$length_filtered_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                echo "  Reads after length filter: $after_length_constraint_reads (${MIN_LENGTH}-${MAX_LENGTH}nt)"
                length_ok=1
            fi
        fi
    fi
    if [ $length_ok -eq 0 ]; then
        echo "  WARNING: empty file after length filter"
        echo "$sample,single,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,0,-1,0,,," >> "$STATS_FILE"
        rm -f "$adapter_trimmed_file" "$length_filtered_file"
        return 1
    fi
    no_contam_file="$CLEAN_DATA_DIR/${sample}_no_rrna_tRNA.fastq"
    echo "  rRNA/tRNA decontamination"
    bowtie2 -x "$RRNA_TRNA_INDEX" -U "$length_filtered_file" --un "$no_contam_file" --no-unal -L 15 -p $THREADS_BOWTIE2 > /dev/null 2>&1
    if [ ! -f "$no_contam_file" ] || [ $(stat -c %s "$no_contam_file") -eq 0 ]; then
        echo "  WARNING: empty file after decontamination"
        echo "$sample,single,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,$after_length_constraint_reads,0,0,,," >> "$STATS_FILE"
        rm -f "$adapter_trimmed_file" "$length_filtered_file"
        return 1
    fi
    after_rtrna_filter_reads=$(cat "$no_contam_file" | wc -l | awk '{print int($1/4)}')
    echo "  Reads after decontamination: $after_rtrna_filter_reads"
    final_file="$CLEAN_DATA_DIR/${sample}_final.fastq.gz"
    temp_length_file="$CLEAN_DATA_DIR/${sample}_lengths.txt"
    gzip -c "$no_contam_file" > "$final_file"
    if [ ! -f "$final_file" ] || [ $(stat -c %s "$final_file") -eq 0 ]; then
        echo "  WARNING: final file empty"
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
        echo "  Final clean reads: $final_reads (${min_length}-${max_length}nt, avg ${avg_length}nt)"
    fi
    echo "  Run FastQC"
    fastqc -o "$FASTQC_CLEAN_DIR" -t $THREADS_FASTQC "$final_file" > /dev/null 2>&1
    if [ -f "$FASTQC_CLEAN_DIR/${sample}_final_fastqc.html" ]; then
        echo "  FastQC report generated"
    else
        echo "  WARNING: FastQC report missing"
    fi
    echo "$sample,single,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,$after_length_constraint_reads,$after_rtrna_filter_reads,$final_reads,$min_length,$max_length,$avg_length" >> "$STATS_FILE"
    rm -f "$temp_length_file" "$adapter_trimmed_file" "$length_filtered_file"
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
    local raw_reads=0 after_trim_adapter_reads=0 after_length_constraint_reads=0 after_rtrna_filter_reads=0 final_reads=0
    local min_length=0 max_length=0 avg_length=0
    local raw_reads1=0 raw_reads2=0
    local after_trim_adapter_reads1=0 after_trim_adapter_reads2=0
    local after_length_constraint_reads1=0 after_length_constraint_reads2=0
    local after_rtrna_filter_reads1=0 after_rtrna_filter_reads2=0 final_reads1=0 final_reads2=0
    echo -e "\n[$current/$total] Process paired-end sample: $base_name (${sample1}/${sample2})"
    raw_file1="$RAW_DATA_DIR/$sample1.fastq.gz"
    raw_file2="$RAW_DATA_DIR/$sample2.fastq.gz"
    if [ ! -f "$raw_file1" ] || [ ! -f "$raw_file2" ]; then
        echo "  ERROR: raw file missing: $raw_file1 or $raw_file2"
        echo "$base_name,paired,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},0,0,0,0,0,,,," >> "$STATS_FILE"
        return 1
    fi
    raw_reads1=$(zcat "$raw_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    raw_reads2=$(zcat "$raw_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    raw_reads=$((raw_reads1 + raw_reads2))
    echo "  Raw reads: $raw_reads (R1: $raw_reads1, R2: $raw_reads2)"
    adapter_trimmed_file1="$CLEAN_DATA_DIR/${sample1}_adapter_trimmed.fastq.gz"
    adapter_trimmed_file2="$CLEAN_DATA_DIR/${sample2}_adapter_trimmed.fastq.gz"
    length_filtered_file1="$CLEAN_DATA_DIR/${sample1}_length_filtered.fastq.gz"
    length_filtered_file2="$CLEAN_DATA_DIR/${sample2}_length_filtered.fastq.gz"
    length_ok=0
    if [ -n "$adapter3_1" ] || [ -n "$adapter5_1" ] || [ -n "$adapter3_2" ] || [ -n "$adapter5_2" ]; then
        echo "  Adapters detected, step1: cutadapt paired adapter trimming only"
        cutadapt_args=()
        [ -n "$adapter3_1" ] && cutadapt_args+=("-a" "$adapter3_1")
        [ -n "$adapter5_1" ] && cutadapt_args+=("-g" "$adapter5_1")
        [ -n "$adapter3_2" ] && cutadapt_args+=("-A" "$adapter3_2")
        [ -n "$adapter5_2" ] && cutadapt_args+=("-G" "$adapter5_2")
        cutadapt "${cutadapt_args[@]}" -j $THREADS_CUTADAPT --trim-n -q 28 --trimmed-only -o "$adapter_trimmed_file1" -p "$adapter_trimmed_file2" "$raw_file1" "$raw_file2" 2>/dev/null
        if [ -f "$adapter_trimmed_file1" ] && [ -f "$adapter_trimmed_file2" ] && [ $(stat -c %s "$adapter_trimmed_file1") -gt 0 ] && [ $(stat -c %s "$adapter_trimmed_file2") -gt 0 ]; then
            after_trim_adapter_reads1=$(zcat "$adapter_trimmed_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_adapter_reads2=$(zcat "$adapter_trimmed_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_adapter_reads=$((after_trim_adapter_reads1 + after_trim_adapter_reads2))
            echo "  Reads after adapter trim: $after_trim_adapter_reads (R1: $after_trim_adapter_reads1, R2: $after_trim_adapter_reads2)"
            echo "  Step2: cutadapt paired length filter ${MIN_LENGTH}-${MAX_LENGTH}nt"
            if [[ "${SELECTED_ANALYSIS}" != "Monosome" || "${FILTER_MONOSOME_LENGTH}" == "true" ]]; then
                cutadapt --minimum-length ${MIN_LENGTH} --maximum-length ${MAX_LENGTH} -j $THREADS_CUTADAPT -o "$length_filtered_file1" "$adapter_trimmed_file1" 2>/dev/null
                cutadapt --minimum-length ${MIN_LENGTH} --maximum-length ${MAX_LENGTH} -j $THREADS_CUTADAPT -o "$length_filtered_file2" "$adapter_trimmed_file2" 2>/dev/null
            else
                cp "$adapter_trimmed_file1" "$length_filtered_file1"
                cp "$adapter_trimmed_file2" "$length_filtered_file2"
            fi
            if [ -f "$length_filtered_file1" ] && [ -f "$length_filtered_file2" ] && [ $(stat -c %s "$length_filtered_file1") -gt 0 ] && [ $(stat -c %s "$length_filtered_file2") -gt 0 ]; then
                after_length_constraint_reads1=$(zcat "$length_filtered_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_constraint_reads2=$(zcat "$length_filtered_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_constraint_reads=$((after_length_constraint_reads1 + after_length_constraint_reads2))
                echo "  Reads after length filter: $after_length_constraint_reads (R1: $after_length_constraint_reads1, R2: $after_length_constraint_reads2) (${MIN_LENGTH}-${MAX_LENGTH}nt)"
                length_ok=1
            fi
        fi
    else
        echo "  No adapter, step1: trim_galore paired quality trim only"
        trim_galore --quality 28 --trim-n --paired --cores $THREADS_TRIMGALORE -o "$CLEAN_DATA_DIR" "$raw_file1" "$raw_file2" 2>/dev/null
        adapter_trimmed_tmp1="${CLEAN_DATA_DIR}/${sample1}_val_1.fq.gz"
        adapter_trimmed_tmp2="${CLEAN_DATA_DIR}/${sample2}_val_2.fq.gz"
        [ ! -f "$adapter_trimmed_tmp1" ] && adapter_trimmed_tmp1="${CLEAN_DATA_DIR}/${base_name}_1_trimmed.fq.gz"
        [ ! -f "$adapter_trimmed_tmp2" ] && adapter_trimmed_tmp2="${CLEAN_DATA_DIR}/${base_name}_2_trimmed.fq.gz"
        if [ -f "$adapter_trimmed_tmp1" ] && [ -f "$adapter_trimmed_tmp2" ]; then
            mv "$adapter_trimmed_tmp1" "$adapter_trimmed_file1"
            mv "$adapter_trimmed_tmp2" "$adapter_trimmed_file2"
            after_trim_adapter_reads1=$(zcat "$adapter_trimmed_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_adapter_reads2=$(zcat "$adapter_trimmed_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
            after_trim_adapter_reads=$((after_trim_adapter_reads1 + after_trim_adapter_reads2))
            echo "  Reads after quality trim: $after_trim_adapter_reads (R1: $after_trim_adapter_reads1, R2: $after_trim_adapter_reads2)"
            echo "  Step2: cutadapt paired length filter ${MIN_LENGTH}-${MAX_LENGTH}nt"
            cutadapt --minimum-length ${MIN_LENGTH} --maximum-length ${MAX_LENGTH} -j $THREADS_CUTADAPT -o "$length_filtered_file1" "$adapter_trimmed_file1" 2>/dev/null
            cutadapt --minimum-length ${MIN_LENGTH} --maximum-length ${MAX_LENGTH} -j $THREADS_CUTADAPT -o "$length_filtered_file2" "$adapter_trimmed_file2" 2>/dev/null
            if [ -f "$length_filtered_file1" ] && [ -f "$length_filtered_file2" ] && [ $(stat -c %s "$length_filtered_file1") -gt 0 ] && [ $(stat -c %s "$length_filtered_file2") -gt 0 ]; then
                after_length_constraint_reads1=$(zcat "$length_filtered_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_constraint_reads2=$(zcat "$length_filtered_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
                after_length_constraint_reads=$((after_length_constraint_reads1 + after_length_constraint_reads2))
                echo "  Reads after length filter: $after_length_constraint_reads (R1: $after_length_constraint_reads1, R2: $after_length_constraint_reads2) (${MIN_LENGTH}-${MAX_LENGTH}nt)"
                length_ok=1
            fi
        fi
    fi
    if [ $length_ok -eq 0 ]; then
        echo "  WARNING: empty file after length filter"
        echo "$base_name,paired,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,0,-1,0,,," >> "$STATS_FILE"
        rm -f "$adapter_trimmed_file1" "$adapter_trimmed_file2" "$length_filtered_file1" "$length_filtered_file2"
        return 1
    fi
    local diff_ratio=0
    local keep_paired="true"
    local single_sample=""
    local single_length_file=""
    if [ $after_length_constraint_reads1 -eq 0 ] || [ $after_length_constraint_reads2 -eq 0 ]; then
        echo "  WARNING: one end empty, switch to single-end decontamination"
        keep_paired="false"
        if [ $after_length_constraint_reads1 -gt $after_length_constraint_reads2 ]; then
            single_sample="$sample1"
            single_length_file="$length_filtered_file1"
        else
            single_sample="$sample2"
            single_length_file="$length_filtered_file2"
        fi
    else
        if [ $after_length_constraint_reads1 -gt $after_length_constraint_reads2 ]; then
            diff_ratio=$(echo "scale=2; $after_length_constraint_reads1 / $after_length_constraint_reads2" | bc)
        else
            diff_ratio=$(echo "scale=2; $after_length_constraint_reads2 / $after_length_constraint_reads1" | bc)
        fi
        if (( $(echo "$diff_ratio > 3" | bc -l) )); then
            echo "  WARNING: ratio=$diff_ratio (>3), switch to single-end decontamination, keep larger end"
            keep_paired="false"
            if [ $after_length_constraint_reads1 -gt $after_length_constraint_reads2 ]; then
                single_sample="$sample1"
                single_length_file="$length_filtered_file1"
            else
                single_sample="$sample2"
                single_length_file="$length_filtered_file2"
            fi
        else
            echo "  Ratio=$diff_ratio (<=3), keep paired-end decontamination"
            keep_paired="true"
        fi
    fi
    no_contam_file1="$CLEAN_DATA_DIR/${sample1}_no_rrna_tRNA.fastq"
    no_contam_file2="$CLEAN_DATA_DIR/${sample2}_no_rrna_tRNA.fastq"
    if [ "$keep_paired" == "true" ]; then
        echo "  Paired rRNA/tRNA decontamination"
        bowtie2 -x "$RRNA_TRNA_INDEX" -1 "$length_filtered_file1" -2 "$length_filtered_file2" --un-conc "$CLEAN_DATA_DIR/${base_name}_no_rrna_tRNA_.fastq" --no-unal -L 15 -p $THREADS_BOWTIE2 > /dev/null 2>&1
        if [ -f "$CLEAN_DATA_DIR/${base_name}_no_rrna_tRNA_.1.fastq" ]; then
            mv "$CLEAN_DATA_DIR/${base_name}_no_rrna_tRNA_.1.fastq" "$no_contam_file1"
        else
            echo "  ERROR: paired decontamination R1 missing"
            echo "$base_name,paired,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,$after_length_constraint_reads,0,0,,," >> "$STATS_FILE"
            rm -f "$adapter_trimmed_file1" "$adapter_trimmed_file2" "$length_filtered_file1" "$length_filtered_file2"
            return 1
        fi
        if [ -f "$CLEAN_DATA_DIR/${base_name}_no_rrna_tRNA_.2.fastq" ]; then
            mv "$CLEAN_DATA_DIR/${base_name}_no_rrna_tRNA_.2.fastq" "$no_contam_file2"
        else
            echo "  ERROR: paired decontamination R2 missing"
            echo "$base_name,paired,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,$after_length_constraint_reads,0,0,,," >> "$STATS_FILE"
            rm -f "$adapter_trimmed_file1" "$adapter_trimmed_file2" "$length_filtered_file1" "$length_filtered_file2"
            return 1
        fi
        after_rtrna_filter_reads1=$(cat "$no_contam_file1" | wc -l | awk '{print int($1/4)}')
        after_rtrna_filter_reads2=$(cat "$no_contam_file2" | wc -l | awk '{print int($1/4)}')
        after_rtrna_filter_reads=$((after_rtrna_filter_reads1 + after_rtrna_filter_reads2))
        echo "  Reads after decontamination: $after_rtrna_filter_reads (R1: $after_rtrna_filter_reads1, R2: $after_rtrna_filter_reads2)"
    else
        echo "  Single-end rRNA/tRNA decontamination, keep $single_sample"
        single_no_contam_file="$CLEAN_DATA_DIR/${single_sample}_no_rrna_tRNA.fastq"
        bowtie2 -x "$RRNA_TRNA_INDEX" -U "$single_length_file" --un "$single_no_contam_file" --no-unal -L 15 -p $THREADS_BOWTIE2 > /dev/null 2>&1
        if [ "$single_sample" == "$sample1" ]; then
            after_rtrna_filter_reads1=$(cat "$single_no_contam_file" | wc -l | awk '{print int($1/4)}')
            after_rtrna_filter_reads2=0
            touch "$no_contam_file2"
        else
            after_rtrna_filter_reads1=0
            after_rtrna_filter_reads2=$(cat "$single_no_contam_file" | wc -l | awk '{print int($1/4)}')
            touch "$no_contam_file1"
        fi
        after_rtrna_filter_reads=$((after_rtrna_filter_reads1 + after_rtrna_filter_reads2))
        echo "  Reads after decontamination: $after_rtrna_filter_reads (only $single_sample kept: $after_rtrna_filter_reads)"
    fi
    if [ $after_rtrna_filter_reads -eq 0 ]; then
        echo "  WARNING: empty after decontamination"
        echo "$base_name,paired,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,$after_length_constraint_reads,0,0,,," >> "$STATS_FILE"
        rm -f "$adapter_trimmed_file1" "$adapter_trimmed_file2" "$length_filtered_file1" "$length_filtered_file2"
        return 1
    fi
    final_file1="$CLEAN_DATA_DIR/${sample1}_final.fastq.gz"
    final_file2="$CLEAN_DATA_DIR/${sample2}_final.fastq.gz"
    final_merged_file="$CLEAN_DATA_DIR/${base_name}_final_merged.fastq.gz"
    temp_length_file="$CLEAN_DATA_DIR/${base_name}_lengths.txt"
    if [ -f "$no_contam_file1" ] && [ $(stat -c %s "$no_contam_file1") -gt 0 ]; then
        gzip -c "$no_contam_file1" > "$final_file1"
        final_reads1=$(zcat "$final_file1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    else
        touch "$final_file1"
        final_reads1=0
        echo "  WARNING: R1 final empty"
    fi
    if [ -f "$no_contam_file2" ] && [ $(stat -c %s "$no_contam_file2") -gt 0 ]; then
        gzip -c "$no_contam_file2" > "$final_file2"
        final_reads2=$(zcat "$final_file2" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    else
        touch "$final_file2"
        final_reads2=0
        echo "  WARNING: R2 final empty"
    fi
    final_reads=$((final_reads1 + final_reads2))
    if [ $final_reads -gt 0 ]; then
        cat <(zcat "$final_file1" 2>/dev/null) <(zcat "$final_file2" 2>/dev/null) | awk 'NR%4==2 {print length($0)}' > "$temp_length_file"
        min_length=$(sort -n "$temp_length_file" | head -n1 2>/dev/null || echo 0)
        max_length=$(sort -n "$temp_length_file" | tail -n1 2>/dev/null || echo 0)
        avg_length=$(awk '{sum+=$1} END {if(NR>0) print sum/NR; else print 0}' "$temp_length_file")
    fi
    echo "  Final clean reads: $final_reads (R1: $final_reads1, R2: $final_reads2) (${min_length}-${max_length}nt, avg ${avg_length}nt)"
    if [ $final_reads -eq 0 ]; then
        echo "  WARNING: both ends empty, create empty merged file"
        touch "$final_merged_file"
    else
        zcat "$final_file1" "$final_file2" 2>/dev/null | gzip > "$final_merged_file"
    fi
    echo "  Merged file: $final_merged_file"
    echo "  Run FastQC"
    fastqc -o "$FASTQC_CLEAN_DIR" -t $THREADS_FASTQC "$final_file1" "$final_file2" > /dev/null 2>&1
    if [ -f "$FASTQC_CLEAN_DIR/${sample1}_final_fastqc.html" ] && [ -f "$FASTQC_CLEAN_DIR/${sample2}_final_fastqc.html" ]; then
        echo "  FastQC reports generated"
    else
        echo "  WARNING: FastQC report incomplete"
    fi
    echo "$base_name,paired,${SELECTED_ANALYSIS},${MIN_LENGTH}-${MAX_LENGTH},$raw_reads,$after_trim_adapter_reads,$after_length_constraint_reads,$after_rtrna_filter_reads,$final_reads,$min_length,$max_length,$avg_length" >> "$STATS_FILE"
    rm -f "$temp_length_file" "$adapter_trimmed_file1" "$adapter_trimmed_file2" "$length_filtered_file1" "$length_filtered_file2"
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
                echo "  Paired detected: $sample ↔ $sample2"
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

echo -e "\n===== Preprocessing Complete ====="
echo "End time: $(date +'%Y-%m-%d %H:%M:%S')"
echo "Total processed: $current | Success: $success | Failed: $failed"
echo "Analysis type: ${SELECTED_ANALYSIS}"
echo "Length filter: ${MIN_LENGTH}-${MAX_LENGTH} nt"
echo "Stats file: $STATS_FILE"
echo "Output file list:"
echo "  * _adapter_trimmed.fastq.gz: only adapter trimmed"
echo "  * _length_filtered.fastq.gz: adapter trimmed + length filtered"
echo "  * _final.fastq.gz: clean decontaminated reads, paired includes merged file"
echo "  fastqc_clean/: FastQC html reports"
echo "Stats column description:"
echo "  after_trim_adapter_reads: reads after adapter / quality trim"
echo "  after_length_constraint_reads: reads after adapter trim + length filter"
echo "=========================="
