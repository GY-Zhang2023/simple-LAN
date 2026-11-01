#!/bin/bash
set -euo pipefail

export DISPLAY=
export JAVA_OPTS="-Djava.awt.headless=true"
export QUALIMAP_JAVA_OPTS="-Djava.awt.headless=true"
export LC_ALL=C

THREADS=8
CORES=2
MEM=50G
WORK_DIR=$(pwd)
HOME_REF="$HOME/reference"
HOME_TOOLS="$HOME/rna_simple_LAN/TOOLS"
HOME_BIN="$HOME/tools"

declare -A SPECIES_CONFIG
SPECIES_CONFIG["human_ref_genome"]="${HOME_REF}/human_transcriptome/GRCh38.p14.canonical_genome.fa"
SPECIES_CONFIG["human_gtf"]="${HOME_REF}/human_transcriptome/gencode.v49.annotation.gtf"
SPECIES_CONFIG["human_star_index"]="${HOME_TOOLS}/STAR_index_GRCh38_v49"
SPECIES_CONFIG["human_hisat2_index"]="${HOME_TOOLS}/HISAT2_index_GRCh38_v49/HISAT2_index_GRCh38_v49"

SPECIES_CONFIG["mouse_ref_genome"]="${HOME_REF}/mouse_transcriptome/GRCm39.p6.canonical_genome.fa"
SPECIES_CONFIG["mouse_gtf"]="${HOME_REF}/mouse_transcriptome/gencode.vM35.annotation.gtf"
SPECIES_CONFIG["mouse_star_index"]="${HOME_TOOLS}/STAR_index_GRCm39_vM35"
SPECIES_CONFIG["mouse_hisat2_index"]="${HOME_TOOLS}/HISAT2_index_GRCm39_vM35/HISAT2_index_GRCm39_vM35"

TOOL_PATHS=(
    "HISAT2_BIN:${HOME_BIN}/hisat2/hisat2"
    "HISAT2_BUILD_BIN:${HOME_BIN}/hisat2/hisat2-build"
    "SAMTOOLS_BIN:samtools"
    "STRINGTIE_BIN:stringtie"
    "QUALIMAP_BIN:qualimap"
    "BEDTOOLS_BIN:bedtools"
    "STAR_BIN:STAR"
    "R_SCRIPT_BIN:Rscript"
)
for tool in "${TOOL_PATHS[@]}"; do
    key=$(echo "$tool" | cut -d':' -f1)
    value=$(echo "$tool" | cut -d':' -f2)
    declare "$key=$value"
done
OUTPUT_BASE="${WORK_DIR}/3.rna_seq_analysis"
MAPQ_THRESHOLD=20
MULTIMAP_MAX=20
SORT_RAM=50000000000
GENOME_RAM_LIMIT=300000000000
SIZE_RATIO_THRESHOLD=3
SJDB_OVERHANG=149
READ_LENGTH=150
STATS_REPORT="${OUTPUT_BASE}/RNAseq_Analysis_Report.txt"

ALIGN_INTRON_MAX=""
ALIGN_MATES_GAP_MAX=""
HISAT2_MP="6,2"
HISAT2_SCORE_MIN="L,0,-0.6"
HISAT2_INTRON_MAX=""
HISAT2_MATES_MAX=""
PE_OVERLAP_MM=0.1
MATCH_MIN_OVER_LREAD=0.3
SCORE_MIN_OVER_LREAD=0.3

SELECTED_TOOL=""
SELECTED_SPECIES=""
STRAND_MODE=""
STRINGTIE_STRAND=""
QUALIMAP_STRAND=""

INPUT_DIR="${WORK_DIR}/2.clean_fastq"
DIR_STRUCTURE=(
    "QUANT_DIR:${OUTPUT_BASE}/01_quantification"
    "ALIGN_DIR:${OUTPUT_BASE}/02_alignments"
    "ASSEMBLY_DIR:${OUTPUT_BASE}/03_assembly"
    "STATS_DIR:${OUTPUT_BASE}/04_statistics"
    "VIS_DIR:${OUTPUT_BASE}/05_visualization"
    "TMP_DIR:${OUTPUT_BASE}/00_tmp"
)
for dir in "${DIR_STRUCTURE[@]}"; do
    key=$(echo "$dir" | cut -d':' -f1)
    path=$(echo "$dir" | cut -d':' -f2)
    declare "$key=$path"
    mkdir -p "$path"
done

COUNT_MATRIX="${OUTPUT_BASE}/gene_count_matrix.tsv"
FPKM_MATRIX="${OUTPUT_BASE}/gene_fpkm_matrix.tsv"
TPM_MATRIX="${OUTPUT_BASE}/gene_tpm_matrix.tsv"
TABLE_REPORT="${STATS_DIR}/RNAseq_Alignment_Stats_Table.tsv"

echo -e "\n====================================="
echo "          Select analysis tool              "
echo "====================================="
echo "1 HISAT2+StringTie"
echo "2 STAR"
read -p "Enter choice (1/2): " TOOL_CHOICE

case ${TOOL_CHOICE} in
    1) SELECTED_TOOL="HISAT2";;
    2) SELECTED_TOOL="STAR";;
    *) echo -e "\nInvalid choice, only 1 or 2"; exit 1 ;;
esac

echo -e "\n====================================="
echo "          Select species              "
echo "====================================="
echo "1 human"
echo "2 mouse"
read -p "Enter choice (1/2): " SPECIES_CHOICE

case ${SPECIES_CHOICE} in
    1)
        SELECTED_SPECIES="human"
        ALIGN_INTRON_MAX=100000
        ALIGN_MATES_GAP_MAX=200000
        HISAT2_INTRON_MAX=100000
        HISAT2_MATES_MAX=200000
        CHROM_LEVELS="c(paste0('chr', 1:22), 'chrX', 'chrY')"
        ;;
    2)
        SELECTED_SPECIES="mouse"
        ALIGN_INTRON_MAX=50000
        ALIGN_MATES_GAP_MAX=100000
        HISAT2_INTRON_MAX=50000
        HISAT2_MATES_MAX=100000
        CHROM_LEVELS="c(paste0('chr', 1:19), 'chrX', 'chrY')"
        ;;
    *) echo -e "\nInvalid choice, only 1 or 2"; exit 1 ;;
esac

echo -e "\n====================================="
echo "          Select library type              "
echo "====================================="
echo "1 Unstranded"
echo "2 Stranded"
read -p "Enter choice (1/2): " STRAND_CHOICE

case ${STRAND_CHOICE} in
    1)
        STRAND_MODE="None"
        STRINGTIE_STRAND=""
        QUALIMAP_STRAND="non-strand-specific"
        ;;
    2)
        STRAND_MODE="intronMotif"
        echo -e "\nSelect strand direction:"
        echo "1 --rf"
        echo "2 --fr"
        read -p "Enter choice (1/2): " STRAND_DIR
        if [ "$STRAND_DIR" = "1" ]; then
            STRINGTIE_STRAND="--rf"
            QUALIMAP_STRAND="strand-specific-reverse"
        else
            STRINGTIE_STRAND="--fr"
            QUALIMAP_STRAND="strand-specific-forward"
        fi
        ;;
    *) echo -e "\nInvalid choice, only 1 or 2"; exit 1 ;;
esac

echo -e "\nFinal config:"
echo "   Aligner: $(if [ $SELECTED_TOOL = "STAR" ]; then echo "STAR"; else echo "HISAT2+StringTie"; fi)"
echo "   Species: ${SELECTED_SPECIES}"
echo "   Library: $(if [ $STRAND_MODE = "None" ]; then echo "Unstranded"; else echo "Stranded"; fi)"
read -p "Continue? (y/n): " CONFIRM
if [ "${CONFIRM}" != "y" ] && [ "${CONFIRM}" != "Y" ]; then
    echo -e "\nUser cancelled"; exit 0
fi

REF_GENOME=${SPECIES_CONFIG[${SELECTED_SPECIES}_ref_genome]}
GTF_FILE=${SPECIES_CONFIG[${SELECTED_SPECIES}_gtf]}
if [ "$SELECTED_TOOL" = "STAR" ]; then
    INDEX=${SPECIES_CONFIG[${SELECTED_SPECIES}_star_index]}
else
    INDEX=${SPECIES_CONFIG[${SELECTED_SPECIES}_hisat2_index]}
fi

echo -e "\nLoaded reference:"
echo "   Index: ${INDEX}"
echo "   Genome: ${REF_GENOME}"
echo "   Annotation: ${GTF_FILE}"

judge_sequencing_type() {
    local sample_prefix=$1
    local single_file="${INPUT_DIR}/${sample_prefix}_final.fastq.gz"
    local file1="${INPUT_DIR}/${sample_prefix}_1_final.fastq.gz"
    local file2="${INPUT_DIR}/${sample_prefix}_2_final.fastq.gz"

    if [ -f "$single_file" ]; then
        echo "single"
    elif [ -f "$file1" ] && [ -f "$file2" ]; then
        echo "paired"
    elif [ -f "$file1" ] || [ -f "$file2" ]; then
        echo "single"
    else
        echo "unknown"
    fi
}

check_paired_size_ratio() {
    local sample_prefix=$1
    local file1="${INPUT_DIR}/${sample_prefix}_1_final.fastq.gz"
    local file2="${INPUT_DIR}/${sample_prefix}_2_final.fastq.gz"

    local size1=$(stat -c%s "$file1" 2>/dev/null || echo 0)
    local size2=$(stat -c%s "$file2" 2>/dev/null || echo 0)

    if [ $size1 -eq 0 ] || [ $size2 -eq 0 ]; then
        echo "single"; return
    fi

    local ratio=$(if [ $size1 -ge $size2 ]; then echo "scale=2; $size1/$size2" | bc; else echo "scale=2; $size2/$size1" | bc; fi)
    if (( $(echo "$ratio > $SIZE_RATIO_THRESHOLD" | bc -l) )); then
        echo "single"
    else
        echo "paired"
    fi
}

get_sample_process_files() {
    declare -A sample_process_map
    declare -A sample_type_map

    for fastq in "$INPUT_DIR"/*_final.fastq.gz; do
        sample_base=$(basename "$fastq" "_final.fastq.gz")
        sample_prefix=${sample_base%_*}
        sample_process_map["$sample_prefix"]="processed"
    done

    for prefix in "${!sample_process_map[@]}"; do
        seq_type=$(judge_sequencing_type "$prefix")
        if [[ "$seq_type" == "unknown" ]]; then
            echo "WARNING: sample $prefix no valid file, skip" >&2
            unset sample_process_map["$prefix"]
            continue
        elif [[ "$seq_type" == "single" ]]; then
            single_file="${INPUT_DIR}/${prefix}_final.fastq.gz"
            file1="${INPUT_DIR}/${prefix}_1_final.fastq.gz"
            file2="${INPUT_DIR}/${prefix}_2_final.fastq.gz"

            if [ -f "$single_file" ]; then
                sample_process_map["$prefix"]="$single_file"
            elif [ -f "$file1" ]; then
                sample_process_map["$prefix"]="$file1"
            elif [ -f "$file2" ]; then
                sample_process_map["$prefix"]="$file2"
            fi
            sample_type_map["$prefix"]="single"
        else
            process_type=$(check_paired_size_ratio "$prefix")
            file1="${INPUT_DIR}/${prefix}_1_final.fastq.gz"
            file2="${INPUT_DIR}/${prefix}_2_final.fastq.gz"
            if [[ "$process_type" == "single" ]]; then
                size1=$(stat -c%s "$file1")
                size2=$(stat -c%s "$file2")
                sample_process_map["$prefix"]=$(if [ $size1 -ge $size2 ]; then echo "$file1"; else echo "$file2"; fi)
                sample_type_map["$prefix"]="paired_as_single"
            else
                sample_process_map["$prefix"]="$file1,$file2"
                sample_type_map["$prefix"]="paired"
            fi
        fi
    done

    for prefix in "${!sample_process_map[@]}"; do
        echo "$prefix:${sample_process_map[$prefix]}:${sample_type_map[$prefix]}"
    done
}

count_fastq_reads() {
    local fastq_file="$1"
    local read_count=0

    if [ ! -f "$fastq_file" ]; then
        echo "0"
        return
    fi

    read_count=$(
        if [[ "$fastq_file" == *.gz ]]; then
            zcat "$fastq_file" 2>/dev/null
        else
            cat "$fastq_file" 2>/dev/null
        fi | wc -l | awk '{print int($1/4)}' | tr -d '[:space:]' | tr -d '[:punct:]'
    )

    if [[ "$read_count" =~ ^[0-9]+$ ]]; then
        echo "$read_count"
    else
        echo "0"
    fi
}

run_hisat2_alignment() {
    local prefix=$1
    local files=$2
    local process_type=$3
    local out_prefix=$4

    echo "Start HISAT2 alignment..."
    if [[ "$process_type" == "single" || "$process_type" == "paired_as_single" ]]; then
        "$HISAT2_BIN" -p "$THREADS" \
        --mp "$HISAT2_MP" \
        --score-min "$HISAT2_SCORE_MIN" \
        --max-intronlen "$HISAT2_INTRON_MAX" \
        -x "$INDEX" \
        -U "$files" \
        | "$SAMTOOLS_BIN" sort -@ "$THREADS" -m "$((SORT_RAM/THREADS))" -o "${out_prefix}Aligned.sortedByCoord.out.bam" -
    else
        IFS=',' read -r file1 file2 <<< "$files"
        "$HISAT2_BIN" -p "$THREADS" \
        --mp "$HISAT2_MP" \
        --score-min "$HISAT2_SCORE_MIN" \
        --max-intronlen "$HISAT2_INTRON_MAX" \
        -X "$HISAT2_MATES_MAX" \
        -x "$INDEX" \
        -1 "$file1" -2 "$file2" \
        | "$SAMTOOLS_BIN" sort -@ "$THREADS" -m "$((SORT_RAM/THREADS))" -o "${out_prefix}Aligned.sortedByCoord.out.bam" -
    fi
    "$SAMTOOLS_BIN" index -@ "$THREADS" "${out_prefix}Aligned.sortedByCoord.out.bam"
    echo "HISAT2 alignment done!"
}

run_stringtie_quant() {
    local prefix=$1
    local bam_file=$2
    local out_gtf=$3
    local out_quant=$4

    echo "Start StringTie assembly and quantification..."
    if [ -n "$STRINGTIE_STRAND" ]; then
        "$STRINGTIE_BIN" -p "$THREADS" \
        $STRINGTIE_STRAND \
        -G "$GTF_FILE" \
        -o "$out_gtf" \
        -A "$out_quant" \
        "$bam_file"
    else
        "$STRINGTIE_BIN" -p "$THREADS" \
        -G "$GTF_FILE" \
        -o "$out_gtf" \
        -A "$out_quant" \
        "$bam_file"
    fi
    echo "StringTie assembly done!"
}

TEMP_GTF="${TMP_DIR}/gencode.annotation.chr.gtf"
mkdir -p "${TMP_DIR}"
cp -p "$GTF_FILE" "$TEMP_GTF"
sed -i.bak -E 's/^([0-9XYM]+)\t/chr\1\t/' "$TEMP_GTF"
rm -f "${TEMP_GTF}.bak"

run_qualimap_rnaseq() {
    local prefix=$1
    local bam_file=$2
    local qualimap_out="${ALIGN_DIR}/${prefix}_qualimap"

    echo "Start Qualimap RNA-seq statistics..."
    local pe_param=""
    if [[ "${sample_type[$prefix]}" == "paired" ]]; then
        pe_param="-pe"
    fi

    if [ ! -f "$bam_file" ] || [ ! -f "$TEMP_GTF" ]; then
        echo "BAM or GTF missing, skip Qualimap"
        return 1
    fi

    "$QUALIMAP_BIN" rnaseq \
    --java-mem-size=${MEM} \
    -bam "$bam_file" \
    -gtf "$TEMP_GTF" \
    -outdir "$qualimap_out" \
    -outformat HTML \
    -p "$QUALIMAP_STRAND" \
    $pe_param

    qualimap_stats="${qualimap_out}/genome_results.txt"
    if [ -f "$qualimap_stats" ]; then
        exon_rate=$(grep "Exonic" "$qualimap_stats" | awk '{print $6}' | sed 's/%//g')
        exon_reads=$(grep "Exonic" "$qualimap_stats" | awk '{print $5}')
        intron_rate=$(grep "Intronic" "$qualimap_stats" | awk '{print $6}' | sed 's/%//g')
        intron_reads=$(grep "Intronic" "$qualimap_stats" | awk '{print $5}')
        intergenic_rate=$(grep "Intergenic" "$qualimap_stats" | awk '{print $6}' | sed 's/%//g')
        intergenic_reads=$(grep "Intergenic" "$qualimap_stats" | awk '{print $5}')

        echo "Exon: $exon_reads reads | ${exon_rate}%"
        echo "Intron: $intron_reads reads | ${intron_rate}%"
        echo "Intergenic: $intergenic_reads reads | ${intergenic_rate}%"

        sample_stats="${STATS_DIR}/${prefix}_stats.txt"
        echo "Exonic reads: ${exon_reads}(${exon_rate}%)" >> "$sample_stats"
        echo "Intronic reads: ${intron_reads}(${intron_rate}%)" >> "$sample_stats"
        echo "Intergenic reads: ${intergenic_reads}(${intergenic_rate}%)" >> "$sample_stats"
    fi
}

run_bedtools_chrom_dist() {
    local prefix=$1
    local bam_file=$2
    local chrom_stats="${ALIGN_DIR}/${prefix}_chrom_dist.txt"
    local chrom_bed="${TMP_DIR}/${prefix}_chrom.bed"

    echo "Start BEDTools chromosome distribution..."
    "$BEDTOOLS_BIN" bamtobed -i "$bam_file" > "$chrom_bed"
    awk '{print $1}' "$chrom_bed" | sort | uniq -c | awk '{print $2 "\t" $1}' > "$chrom_stats"

    echo "Top10 chromosome reads:"
    head -10 "$chrom_stats" | while read chrom count; do
        echo "   $chrom: $count reads"
    done

    draw_chrom_density "$prefix" "$chrom_stats" "$VIS_DIR"
}

draw_chrom_density() {
    local prefix=$1
    local chrom_stats=$2
    local out_dir=$3
    local r_script="${TMP_DIR}/${prefix}_chrom_density.R"
    local out_png="${out_dir}/${prefix}_chrom_density.png"

    cat > "$r_script" << EOF
library(ggplot2)
library(dplyr)
chrom_data <- read.table("$chrom_stats", header = FALSE, col.names = c("chrom", "count"))
chrom_data <- chrom_data %>% filter(grepl("^chr[0-9XY]+$", chrom)) %>% arrange(chrom)
chrom_data\$chrom <- factor(chrom_data\$chrom, levels = $CHROM_LEVELS)
ggplot(chrom_data, aes(x = chrom, y = count/1000000)) +
  geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
  geom_smooth(aes(group=1), method = "loess", se = FALSE, color = "#A23B72", linewidth = 1) +
  labs(title = "Reads Density Distribution on Chromosomes", subtitle = "Sample: $prefix", x = "Chromosome", y = "Reads Count (Million)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("$out_png", width = 12, height = 6, dpi = 300)
EOF

    echo "Plot chromosome density..."
    "$R_SCRIPT_BIN" "$r_script"
    echo "Chromosome density plot saved: $out_png"
}

echo -e "\n====================================="
echo "          Pipeline init                  "
echo "====================================="
echo "========================================" > "$STATS_REPORT"
echo "poly-A RNA-seq analysis report" >> "$STATS_REPORT"
echo "Time: $(date +"%Y-%m-%d %H:%M:%S")" >> "$STATS_REPORT"
echo "Aligner: $SELECTED_TOOL" >> "$STATS_REPORT"
echo "Species: ${SELECTED_SPECIES}" >> "$STATS_REPORT"
echo "Library: $(if [ $STRAND_MODE = "None" ]; then echo "Unstranded"; else echo "Stranded"; fi)" >> "$STATS_REPORT"
echo "========================================" >> "$STATS_REPORT"

echo -e "\n====================================="
echo "          Build index                    "
echo "====================================="
if [ "$SELECTED_TOOL" = "STAR" ]; then
    if [ ! -f "$INDEX/Genome" ]; then
        echo "Building STAR index..."
        "$STAR_BIN" --runMode genomeGenerate \
        --genomeDir "$INDEX" \
        --genomeFastaFiles "$REF_GENOME" \
        --sjdbGTFfile "$GTF_FILE" \
        --sjdbOverhang "$SJDB_OVERHANG" \
        --runThreadN "$THREADS" \
        --genomeSAindexNbases 14 \
        --limitGenomeGenerateRAM "$GENOME_RAM_LIMIT"
        echo "STAR index built!"
    else
        echo "STAR index exists, skip"
    fi
else
    if [ ! -f "$INDEX.1.ht2" ]; then
        echo "Building HISAT2 index..."
        "$HISAT2_BUILD_BIN" -p "$THREADS" "$REF_GENOME" "$INDEX"
        echo "HISAT2 index built!"
    else
        echo "HISAT2 index exists, skip"
    fi
fi

echo -e "\n====================================="
echo "          Sample processing                    "
echo "====================================="
sample_list=$(get_sample_process_files)
declare -A sample_files=()
declare -A sample_type=()

while IFS=: read -r prefix files type; do
    sample_files["$prefix"]="$files"
    sample_type["$prefix"]="$type"
done <<< "$sample_list"

if [ ${#sample_files[@]} -eq 0 ]; then
    echo -e "\nNo valid samples found, exit"; exit 0
fi
echo -e "\nTotal ${#sample_files[@]} samples:"
for prefix in "${!sample_files[@]}"; do
    echo "   - $prefix (type: ${sample_type[$prefix]})"
done

for prefix in "${!sample_files[@]}"; do
    echo -e "\nProcess sample: $prefix"
    files="${sample_files[$prefix]}"
    process_type="${sample_type[$prefix]}"
    out_prefix="${ALIGN_DIR}/${prefix}_"
    quant_file="${QUANT_DIR}/${prefix}_FPKM_TPM.tsv"
    assembly_gtf="${ASSEMBLY_DIR}/${prefix}_assembly.gtf"
    bam_file="${out_prefix}Aligned.sortedByCoord.out.bam"
    filtered_bam="${ALIGN_DIR}/${prefix}_filtered.bam"

    echo -e "\n【Sample $prefix】" >> "$STATS_REPORT"
    echo "Process type: $process_type" >> "$STATS_REPORT"
    echo "Input: $files" >> "$STATS_REPORT"

    echo "Count raw reads..."
    if [[ "$process_type" == "single" || "$process_type" == "paired_as_single" ]]; then
        raw_reads=$(count_fastq_reads "$files")
        echo "Raw reads: $raw_reads"
        echo "Raw reads: $raw_reads" >> "$STATS_REPORT"
    else
        IFS=',' read -r file1 file2 <<< "$files"
        raw_reads1=$(count_fastq_reads "$file1")
        raw_reads2=$(count_fastq_reads "$file2")
        raw_reads1=$((raw_reads1 + 0))
        raw_reads2=$((raw_reads2 + 0))
        raw_reads=$((raw_reads1 + raw_reads2))
        echo "Raw reads R1: $raw_reads1 | R2: $raw_reads2 | total: $raw_reads"
        echo "Raw reads R1: $raw_reads1" >> "$STATS_REPORT"
        echo "Raw reads R2: $raw_reads2" >> "$STATS_REPORT"
        echo "Raw reads total: $raw_reads" >> "$STATS_REPORT"
    fi

    if [ "$SELECTED_TOOL" = "STAR" ]; then
        if [ ! -f "${out_prefix}Log.final.out" ]; then
            echo "Start STAR alignment..."
            if [[ "$process_type" == "single" || "$process_type" == "paired_as_single" ]]; then
                "$STAR_BIN" --genomeDir "$INDEX" \
                --readFilesIn "$files" \
                --readFilesCommand zcat \
                --outFileNamePrefix "$out_prefix" \
                --runThreadN "$THREADS" \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped None \
                --outSAMprimaryFlag AllBestScore \
                --outFilterMultimapNmax "$MULTIMAP_MAX" \
                --outSAMmapqUnique "$MAPQ_THRESHOLD" \
                --quantMode GeneCounts \
                --twopassMode Basic \
                --outFilterMatchNminOverLread "$MATCH_MIN_OVER_LREAD" \
                --outFilterScoreMinOverLread "$SCORE_MIN_OVER_LREAD" \
                --alignIntronMax "$ALIGN_INTRON_MAX" \
                --outSAMstrandField "$STRAND_MODE" \
                --limitBAMsortRAM "$SORT_RAM"
            else
                "$STAR_BIN" --genomeDir "$INDEX" \
                --readFilesIn "$file1" "$file2" \
                --readFilesCommand zcat \
                --outFileNamePrefix "$out_prefix" \
                --runThreadN "$THREADS" \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped None \
                --outSAMprimaryFlag AllBestScore \
                --outFilterMultimapNmax "$MULTIMAP_MAX" \
                --outSAMmapqUnique "$MAPQ_THRESHOLD" \
                --quantMode GeneCounts \
                --peOverlapNbasesMin 10 \
                --peOverlapMMp "$PE_OVERLAP_MM" \
                --alignMatesGapMax "$ALIGN_MATES_GAP_MAX" \
                --alignIntronMax "$ALIGN_INTRON_MAX" \
                --twopassMode Basic \
                --outSAMstrandField "$STRAND_MODE" \
                --limitBAMsortRAM "$SORT_RAM" \
                --outFilterMatchNminOverLread "$MATCH_MIN_OVER_LREAD" \
                --outFilterScoreMinOverLread "$SCORE_MIN_OVER_LREAD"
            fi
            echo "STAR alignment done!"
        else
            echo "STAR finished, skip"
        fi

        log_file="${out_prefix}Log.final.out"
        cp "${out_prefix}ReadsPerGene.out.tab" "${QUANT_DIR}/${prefix}_ReadsPerGene.out.tab"
        echo "STAR quantification file saved"
    else
        if [ ! -f "$bam_file" ]; then
            run_hisat2_alignment "$prefix" "$files" "$process_type" "$out_prefix"
        else
            echo "HISAT2 finished, skip"
        fi

        if [ ! -f "$assembly_gtf" ]; then
            run_stringtie_quant "$prefix" "$bam_file" "$assembly_gtf" "$quant_file"
        else
            echo "StringTie finished, skip"
        fi

        if [ ! -f "$filtered_bam" ]; then
            "$SAMTOOLS_BIN" view -b -q "$MAPQ_THRESHOLD" -@ "$THREADS" "$bam_file" -o "$filtered_bam"
            "$SAMTOOLS_BIN" index -@ "$THREADS" "$filtered_bam"
            filtered_reads=$("$SAMTOOLS_BIN" view -c -F 4 "$filtered_bam")
            echo "Filtered high-confidence reads: $filtered_reads"
            echo "Filtered high-confidence reads (MAPQ≥20): $filtered_reads" >> "$STATS_REPORT"
        else
            echo "Filtered BAM exists, skip"
        fi
    fi

    run_qualimap_rnaseq "$prefix" "$bam_file"
    run_bedtools_chrom_dist "$prefix" "$bam_file"
done

echo -e "\n====================================="
echo "          Generate quantification matrices              "
echo "====================================="
sample_names=$(echo "${!sample_files[@]}" | tr ' ' '\t')
if [ "$SELECTED_TOOL" = "STAR" ]; then
    first_sample=$(echo "${!sample_files[@]}" | awk '{print $1}')
    first_quant="${QUANT_DIR}/${first_sample}_ReadsPerGene.out.tab"
    awk '$1 !~ /^N_/ {print $1}' "$first_quant" > "${TMP_DIR}/gene_ids.tmp"

cat > "${TMP_DIR}/generate_matrices.R" << EOF
library(dplyr)
library(rtracklayer)
sample_list <- c($(printf '"%s",' "${!sample_files[@]}" | sed 's/,$//'))
quant_dir <- "${QUANT_DIR}"
gene_ids <- read.table("${TMP_DIR}/gene_ids.tmp", header=F, stringsAsFactors=F, col.names = "gene_id")[,1]
strand_mode <- "${STRAND_MODE}"

count_mat <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
fpkm_mat <- count_mat
tpm_mat <- count_mat

for (sample in sample_list) {
  quant_file <- paste0(quant_dir, "/", sample, "_ReadsPerGene.out.tab")
  dat <- read.table(quant_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                    fill = TRUE, quote = "", col.names = c("gene_id", "pos", "neg", "unstrand"))
  dat <- dat %>% filter(!grepl("^N_", gene_id)) %>% select(gene_id, pos, neg, unstrand)
  dat <- dat[match(gene_ids, dat\$gene_id), ]

  if (strand_mode == "None") {
    dat\$raw_count <- dat\$unstrand
  } else {
    dat\$raw_count <- dat\$neg
  }

  dat\$raw_count <- as.numeric(dat\$raw_count)
  dat\$raw_count[is.na(dat\$raw_count)] <- 0
  dat\$raw_count <- pmax(dat\$raw_count, 0)
  dat\$raw_count <- as.integer(dat\$raw_count)

  count_mat[[sample]] <- dat\$raw_count
}

gtf_file <- "${GTF_FILE}"
gtf <- import(gtf_file)
gtf_gene <- gtf[gtf\$type == "gene"]
gene_len <- data.frame(
  gene_id = gtf_gene\$gene_id,
  length = width(gtf_gene),
  stringsAsFactors = FALSE
) %>% distinct(gene_id, .keep_all = TRUE)
gene_len <- gene_len[match(gene_ids, gene_len\$gene_id), ]
gene_length <- gene_len\$length
gene_length[is.na(gene_length)] <- 0
gene_length <- pmax(gene_length, 1)

for (sample in sample_list) {
  total_mapped <- sum(count_mat[[sample]]) / 1e6
  fpkm_mat[[sample]] <- count_mat[[sample]] / (gene_length / 1000) / total_mapped
  fpkm_mat[[sample]] <- pmax(fpkm_mat[[sample]], 0)
}

for (sample in sample_list) {
  fpkm_sum <- sum(fpkm_mat[[sample]], na.rm = TRUE)
  tpm_mat[[sample]] <- (fpkm_mat[[sample]] / fpkm_sum) * 1e6
  tpm_mat[[sample]] <- pmax(tpm_mat[[sample]], 0)
}

count_mat <- count_mat %>% mutate(across(-gene_id, ~pmax(as.integer(.), 0)))

write.table(count_mat, "${COUNT_MATRIX}", sep = "\t", row.names = FALSE, quote = FALSE, na = "0")
write.table(fpkm_mat, "${FPKM_MATRIX}", sep = "\t", row.names = FALSE, quote = FALSE, na = "0")
write.table(tpm_mat, "${TPM_MATRIX}", sep = "\t", row.names = FALSE, quote = FALSE, na = "0")
EOF
    "$R_SCRIPT_BIN" "${TMP_DIR}/generate_matrices.R"
else
    first_sample=$(echo "${!sample_files[@]}" | awk '{print $1}')
    first_quant="${QUANT_DIR}/${first_sample}_FPKM_TPM.tsv"
    grep -v "^Gene ID" "$first_quant" | awk '{print $1}' > "${TMP_DIR}/gene_ids.tmp"

cat > "${TMP_DIR}/generate_matrices.R" << EOF
library(dplyr)
sample_list <- c($(printf '"%s",' "${!sample_files[@]}" | sed 's/,$//'))
quant_dir <- "${QUANT_DIR}"
gene_ids <- read.table("${TMP_DIR}/gene_ids.tmp", header=F, stringsAsFactors=F)[,1]
read_length <- ${READ_LENGTH}

count_mat <- data.frame(gene_id = gene_ids)
fpkm_mat <- data.frame(gene_id = gene_ids)
tpm_mat <- data.frame(gene_id = gene_ids)

for (sample in sample_list) {
  quant_file <- paste0(quant_dir, "/", sample, "_FPKM_TPM.tsv")
  dat <- read.table(quant_file, header=T, sep="\t", stringsAsFactors=F)

  dat <- dat[match(gene_ids, dat\$Gene.ID),]

  dat\$gene_length <- pmax(dat\$End - dat\$Start + 1, 1)
  dat\$coverage_clean <- pmax(dat\$Coverage, 0)
  dat\$raw_count <- round(dat\$coverage_clean * dat\$gene_length / read_length)

  dat\$raw_count <- as.integer(pmax(dat\$raw_count, 0))

  dat\$raw_count[is.na(dat\$raw_count)] <- 0
  dat\$FPKM[is.na(dat\$FPKM)] <- 0
  dat\$TPM[is.na(dat\$TPM)] <- 0

  count_mat[[sample]] <- dat\$raw_count
  fpkm_mat[[sample]] <- dat\$FPKM
  tpm_mat[[sample]] <- dat\$TPM
}

count_mat <- count_mat %>% mutate(across(-gene_id, ~pmax(., 0)))
fpkm_mat <- fpkm_mat %>% mutate(across(-gene_id, ~pmax(., 0)))
tpm_mat <- tpm_mat %>% mutate(across(-gene_id, ~pmax(., 0)))

write.table(count_mat, "${COUNT_MATRIX}", sep="\t", row.names=F, quote=F)
write.table(fpkm_mat, "${FPKM_MATRIX}", sep="\t", row.names=F, quote=F)
write.table(tpm_mat, "${TPM_MATRIX}", sep="\t", row.names=F, quote=F)
EOF
    "$R_SCRIPT_BIN" "${TMP_DIR}/generate_matrices.R"
fi

rm -f "${TMP_DIR}/gene_ids.tmp"
echo "Raw count matrix: $COUNT_MATRIX"
echo "FPKM matrix: $FPKM_MATRIX"
echo "TPM matrix: $TPM_MATRIX"
echo "Three matrices aligned by gene ID, ready for downstream analysis!"

echo -e "\n====================================="
echo "          Pipeline finished              "
echo "====================================="
echo -e "\nOutput directory: $OUTPUT_BASE"
echo -e "\nMain quantification matrices:"
echo "   1. $COUNT_MATRIX → integer raw count for DESeq2/edgeR"
echo "   2. $FPKM_MATRIX → FPKM for intra-sample comparison"
echo "   3. $TPM_MATRIX → TPM for inter-sample clustering, heatmap, PCA"
echo -e "\nOther files:"
if [ "$SELECTED_TOOL" = "STAR" ]; then
    echo "   - BAM: $ALIGN_DIR/*_Aligned.sortedByCoord.out.bam"
    echo "   - Quant files: $QUANT_DIR/*_ReadsPerGene.out.tab"
else
    echo "   - Filtered BAM: $ALIGN_DIR/*_filtered.bam"
    echo "   - Assembled GTF: $ASSEMBLY_DIR/*_assembly.gtf"
    echo "   - Gene quant: $QUANT_DIR/*_FPKM_TPM.tsv"
fi
echo "   - Report: $STATS_REPORT"
echo "   - Chromosome density plot: $VIS_DIR/*_chrom_density.png"
echo -e "\nAnalysis suggestions:"
echo "   ✔️ DEG: $COUNT_MATRIX"
echo "   ✔️ Heatmap/PCA: $TPM_MATRIX"
echo "   ✔️ Intra-sample rank: $FPKM_MATRIX"
echo "   ✔️ Alternative splicing: GTF in $ASSEMBLY_DIR"

rm -rf "${TMP_DIR}"
