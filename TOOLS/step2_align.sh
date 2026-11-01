#!/bin/bash
set -euo pipefail

##############################################################################
# Adjustable parameters
# Default: 8 threads, 50G RAM
##############################################################################
THREADS=8
SORT_RAM=50000000000
GENOME_RAM_LIMIT=300000000000
MAPQ_THRESH=20
SJDB_OVERHANG=100
SIZE_RATIO_THRESHOLD=3

declare -A SPECIES_CONFIG
# human
SPECIES_CONFIG["human_ref_genome"]="$HOME/reference/human_transcriptome/GRCh38.p14.canonical_genome.fa"
SPECIES_CONFIG["human_gtf"]="$HOME/reference/human_transcriptome/gencode.v49.pc_annotation.gtf"
SPECIES_CONFIG["human_transcript_fa"]="$HOME/reference/human_transcriptome/gencode.v49.pc_transcripts.fa"
SPECIES_CONFIG["human_star_index"]="$HOME/ribo_simple_LAN/TOOLS/star_genome_index"

# mouse
SPECIES_CONFIG["mouse_ref_genome"]="$HOME/reference/mouse_transcriptome/GRCm39.p6.canonical_genome.fa"
SPECIES_CONFIG["mouse_gtf"]="$HOME/reference/mouse_transcriptome/gencode.vM35.pc_annotation.gtf"
SPECIES_CONFIG["mouse_transcript_fa"]="$HOME/reference/mouse_transcriptome/gencode.vM35.pc_transcripts.fa"
SPECIES_CONFIG["mouse_star_index"]="$HOME/ribo_simple_LAN/TOOLS/star_mus_genome_index"

##############################################################################
# Species selection
##############################################################################
echo -e "\n====================================="
echo "          Select species              "
echo "====================================="
echo "1 human (gencode.v49)"
echo "2 mouse (gencode.vM35)"
read -p "Input choice (1/2): " SPECIES_CHOICE

case ${SPECIES_CHOICE} in
    1) SELECTED_SPECIES="human" ;;
    2) SELECTED_SPECIES="mouse" ;;
    *) echo -e "\nERROR: Invalid choice! Only 1 (human) or 2 (mouse) supported"; exit 1 ;;
esac

echo -e "\nWARNING: Species ${SELECTED_SPECIES} will be used for analysis"
read -p "Confirm continue? (y/n): " CONFIRM
if [ "${CONFIRM}" != "y" ] && [ "${CONFIRM}" != "Y" ]; then
    echo -e "\nUser cancelled, exit"; exit 0
fi

##############################################################################
# Path definition
##############################################################################
WORK_DIR=$(pwd)
INPUT_DIR="${WORK_DIR}/2.clean_fastq"
OUTPUT_DIR="${WORK_DIR}/3.alignments"
STATS_REPORT="$OUTPUT_DIR/Report.txt"
PLOT_DIR="${OUTPUT_DIR}/RPF_length_plots"
GENOME_BAM_SUFFIX="_genome_Aligned.toTranscriptome.out.bam"

STAR_INDEX=${SPECIES_CONFIG[${SELECTED_SPECIES}_star_index]}
REF_GENOME=${SPECIES_CONFIG[${SELECTED_SPECIES}_ref_genome]}
GTF_FILE=${SPECIES_CONFIG[${SELECTED_SPECIES}_gtf]}
TRANSCRIPT_FASTA=${SPECIES_CONFIG[${SELECTED_SPECIES}_transcript_fa]}

echo -e "\nLoaded config for ${SELECTED_SPECIES}:"
echo "   STAR index: ${STAR_INDEX}"
echo "   Reference genome: ${REF_GENOME}"
echo "   GTF: ${GTF_FILE}"
echo "   Transcript fasta: ${TRANSCRIPT_FASTA}"
echo -e "\nInput dir: ${INPUT_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "RPF plot dir: ${PLOT_DIR}"
echo -e "=====================================================\n"

##############################################################################
# Functions
##############################################################################
judge_sequencing_type() {
    local sample_prefix=$1
    local single_file="${INPUT_DIR}/${sample_prefix}_final.fastq.gz"
    local file1="${INPUT_DIR}/${sample_prefix}_1_final.fastq.gz"
    local file2="${INPUT_DIR}/${sample_prefix}_2_final.fastq.gz"
    if [ -f "$single_file" ]; then echo "single";
    elif [ -f "$file1" ] && [ -f "$file2" ]; then echo "paired";
    elif [ -f "$file1" ] || [ -f "$file2" ]; then echo "single";
    else echo "unknown"; fi
}

check_paired_size_ratio() {
    local sample_prefix=$1
    local file1="${INPUT_DIR}/${sample_prefix}_1_final.fastq.gz"
    local file2="${INPUT_DIR}/${sample_prefix}_2_final.fastq.gz"
    local size1=$(stat -c%s "$file1" 2>/dev/null || echo 0)
    local size2=$(stat -c%s "$file2" 2>/dev/null || echo 0)
    if [ $size1 -eq 0 ] || [ $size2 -eq 0 ]; then echo "single"; return; fi
    local ratio=$(echo "scale=2; if($size1>=$size2) $size1/$size2 else $size2/$size1" | bc)
    echo "  Pair file size: $(basename $file1)=$size1, $(basename $file2)=$size2, ratio=$ratio" >> "$STATS_REPORT"
    if (( $(echo "$ratio > $SIZE_RATIO_THRESHOLD" | bc -l) )); then echo "single"; else echo "paired"; fi
}

get_sample_process_files() {
    declare -A sample_process_map sample_type_map
    for fastq in "$INPUT_DIR"/*_final.fastq.gz; do
        sample_base=$(basename "$fastq" "_final.fastq.gz")
        sample_prefix=${sample_base%_*}
        [[ "$sample_base" =~ _1$ || "$sample_base" =~ _2$ ]] || sample_prefix=$sample_base
        sample_process_map["$sample_prefix"]="processed"
    done
    for prefix in "${!sample_process_map[@]}"; do
        seq_type=$(judge_sequencing_type "$prefix")
        if [[ "$seq_type" == "unknown" ]]; then unset sample_process_map["$prefix"]; continue; fi
        if [[ "$seq_type" == "single" ]]; then
            single_file="${INPUT_DIR}/${prefix}_final.fastq.gz"
            file1="${INPUT_DIR}/${prefix}_1_final.fastq.gz"
            file2="${INPUT_DIR}/${prefix}_2_final.fastq.gz"
            [ -f "$single_file" ] && sample_process_map["$prefix"]="$single_file" || sample_process_map["$prefix"]="$file1"
            sample_type_map["$prefix"]="single"
        else
            process_type=$(check_paired_size_ratio "$prefix")
            if [[ "$process_type" == "single" ]]; then
                file1="${INPUT_DIR}/${prefix}_1_final.fastq.gz"
                file2="${INPUT_DIR}/${prefix}_2_final.fastq.gz"
                sample_process_map["$prefix"]=$( [ $(stat -c%s "$file1") -ge $(stat -c%s "$file2") ] && echo "$file1" || echo "$file2" )
                sample_type_map["$prefix"]="paired_as_single"
            else
                sample_process_map["$prefix"]="$file1,$file2"
                sample_type_map["$prefix"]="paired"
            fi
        fi
    done
    for prefix in "${!sample_process_map[@]}"; do echo "$prefix:${sample_process_map[$prefix]}:${sample_type_map[$prefix]}"; done
}

count_fastq_reads() { zcat "$1" 2>/dev/null | wc -l | awk '{print int($1/4)}'; }

star_paired_alignment() {
    local prefix=$1 file1=$2 file2=$3 out_prefix=$4
    echo "Start paired alignment: $prefix"
    [ -f "${out_prefix}Log.final.out" ] && { echo "$prefix paired alignment finished, skip"; return; }
    STAR --genomeDir "$STAR_INDEX" --readFilesIn "$file1" "$file2" --readFilesCommand zcat --outFileNamePrefix "$out_prefix" --runThreadN "$THREADS" --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMprimaryFlag AllBestScore --outFilterMultimapNmax 1 --outSAMmapqUnique 20 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMstrandField None --limitBAMsortRAM "$SORT_RAM"
    echo "$prefix paired alignment done"
}

##############################################################################
# Step 1: Environment setup
##############################################################################
echo "========================================"
echo "Step 1: Environment setup"
echo "========================================"
mkdir -p "$OUTPUT_DIR" "$STAR_INDEX" "$PLOT_DIR"
chmod -R 777 "$(dirname "$OUTPUT_DIR")" "$(dirname "$STAR_INDEX")" "$PLOT_DIR"
ulimit -n 4096
echo "Output dir: $OUTPUT_DIR"
echo "Index dir: $STAR_INDEX"
echo "Plot dir: $PLOT_DIR"
echo "Current ulimit -n: $(ulimit -n)"

##############################################################################
# Step 0: Initialize report
##############################################################################
echo "========================================"
echo "Step 0: Initialize report"
echo "========================================"
echo "========================================" > "$STATS_REPORT"
echo "Ribo-seq alignment report" >> "$STATS_REPORT"
echo "Time: $(date +"%Y-%m-%d %H:%M:%S")" >> "$STATS_REPORT"
echo "Reference genome: $(basename "$REF_GENOME")" >> "$STATS_REPORT"
echo "GTF: $(basename "$GTF_FILE")" >> "$STATS_REPORT"
echo "MAPQ threshold: $MAPQ_THRESH" >> "$STATS_REPORT"
echo "Region: 5UTR,5UTR-CDS border,CDS,CDS-3UTR border,3UTR" >> "$STATS_REPORT"
echo "========================================" >> "$STATS_REPORT"
echo -e "\n" >> "$STATS_REPORT"
echo "Report initialized"

##############################################################################
# Step 2: Build STAR index
##############################################################################
echo -e "\n========================================"
echo "Step 2: Build STAR index (run once)"
echo "========================================"
if [ ! -f "$STAR_INDEX/Genome" ]; then
  echo "Start index building..."
  STAR --runMode genomeGenerate --genomeDir "$STAR_INDEX" --genomeFastaFiles "$REF_GENOME" --sjdbGTFfile "$GTF_FILE" --sjdbOverhang "$SJDB_OVERHANG" --runThreadN "$THREADS" --genomeSAindexNbases 14 --genomeSAsparseD 2 --limitGenomeGenerateRAM "$GENOME_RAM_LIMIT"
  echo "Index built: $STAR_INDEX"
else
  echo "Index exists, skip"
fi

##############################################################################
# Step 3: Sample detection & alignment
##############################################################################
echo -e "\n========================================"
echo "Step 3: Detect samples and alignment"
echo "========================================"
echo "Detect samples..."
sample_list=$(get_sample_process_files)
declare -A sample_files=() sample_process_type=()
while IFS=: read -r prefix files type; do
    [ -z "$prefix" ] || [ -z "$files" ] || [ -z "$type" ] && continue
    sample_files["$prefix"]="$files"
    sample_process_type["$prefix"]="$type"
done <<< "$sample_list"

[ ${#sample_files[@]} -eq 0 ] && { echo -e "\nWARNING: No valid samples found"; exit 0; }
echo -e "\nSample list (${#sample_files[@]} samples):"
for prefix in "${!sample_files[@]}"; do echo "  - $prefix: type=${sample_process_type[$prefix]}, file=${sample_files[$prefix]}"; done

for prefix in "${!sample_files[@]}"; do
    files="${sample_files[$prefix]}"
    process_type="${sample_process_type[$prefix]}"
    echo -e "\nProcess sample: $prefix (type: $process_type)"
    echo -e "\n【Sample $prefix start】" >> "$STATS_REPORT"
    echo "Sample name: $prefix" >> "$STATS_REPORT"
    echo "Process type: $process_type" >> "$STATS_REPORT"

    if [[ "$process_type" == "single" || "$process_type" == "paired_as_single" ]]; then
        fastq="$files"
        echo "Input file: $fastq" >> "$STATS_REPORT"
        raw_reads=$(count_fastq_reads "$fastq")
        echo "Clean reads before alignment: $raw_reads"
        echo "Clean reads before alignment: $raw_reads" >> "$STATS_REPORT"
    else
        IFS=',' read -r file1 file2 <<< "$files"
        echo "R1: $file1" >> "$STATS_REPORT"
        echo "R2: $file2" >> "$STATS_REPORT"
        raw_reads1=$(count_fastq_reads "$file1")
        raw_reads2=$(count_fastq_reads "$file2")
        raw_reads=$((raw_reads1 + raw_reads2))
        echo "Clean reads R1: $raw_reads1"
        echo "Clean reads R2: $raw_reads2"
        echo "Clean reads total: $raw_reads"
        echo "Clean reads R1: $raw_reads1" >> "$STATS_REPORT"
        echo "Clean reads R2: $raw_reads2" >> "$STATS_REPORT"
        echo "Clean reads total: $raw_reads" >> "$STATS_REPORT"
    fi

    out_prefix="$OUTPUT_DIR/${prefix}_genome_"
    genome_bam="${out_prefix}Aligned.sortedByCoord.out.bam"
    if [ ! -f "${out_prefix}Log.final.out" ]; then
        if [[ "$process_type" == "single" || "$process_type" == "paired_as_single" ]]; then
            echo "Run single-end alignment..."
            STAR --genomeDir "$STAR_INDEX" --readFilesIn "$fastq" --readFilesCommand zcat --outFileNamePrefix "$out_prefix" --runThreadN "$THREADS" --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMprimaryFlag AllBestScore --outFilterMultimapNmax 1 --outSAMmapqUnique 20 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMstrandField None --limitBAMsortRAM "$SORT_RAM"
        else
            star_paired_alignment "$prefix" "$file1" "$file2" "$out_prefix"
        fi
    else
        echo "$prefix alignment finished, skip"
    fi

    log_file="${out_prefix}Log.final.out"
    [ ! -f "$log_file" ] && { echo "WARNING: Log missing for $prefix, skip stats"; continue; }
    aligned_reads=$(grep "Uniquely mapped reads number" "$log_file" | awk -F'\t' '{print $2}')
    aligned_rate=$(grep "Uniquely mapped reads %" "$log_file" | awk -F'\t' '{print $2}')
    multimap_reads=$(grep "Number of reads mapped to multiple loci" "$log_file" | awk -F'\t' '{print $2}')
    echo "Unique mapped reads: $aligned_reads"
    echo "Unique mapping rate: $aligned_rate"
    echo "Multi-mapped reads: $multimap_reads"
    echo "Unique mapped reads: $aligned_reads" >> "$STATS_REPORT"
    echo "Unique mapping rate: $aligned_rate" >> "$STATS_REPORT"
    echo "Multi-mapped reads: $multimap_reads" >> "$STATS_REPORT"
    total_aligned_rate=$([ $raw_reads -gt 0 ] && echo $((aligned_reads * 100 / raw_reads)) || echo 0)
    echo "Unique mapping rate (raw): $total_aligned_rate%" >> "$STATS_REPORT"
    echo "$prefix alignment stats done"
    echo -e "【Sample $prefix alignment end】\n" >> "$STATS_REPORT"
    echo "$aligned_reads" > "$OUTPUT_DIR/${prefix}.genome_uniq.tmp"
done

##############################################################################
# Step4: Annotation generation & best-hit BAM filter
##############################################################################
echo -e "\n========================================"
echo "Step4: Generate annotation and filter BAM (best hit by FLAG)"
echo "========================================"

ANNOT_DIR="$OUTPUT_DIR"
ANNOT_PATH="$ANNOT_DIR/ribowaltz_annotation.tsv"
CLASSIFY_ANNOT="$ANNOT_DIR/tx_annotation_for_classify.tsv"
mkdir -p "$ANNOT_DIR"

if [ ! -f "$CLASSIFY_ANNOT" ]; then
  echo "Run R script to generate annotation..."
  Rscript -e "
assignInNamespace(
  'makeTxDbFromGFF',
  function(...) {
    txdbmaker::makeTxDbFromGFF(...)
  },
  ns = 'GenomicFeatures'
)

library(S4Vectors)
library(txdbmaker)
library(riboWaltz)
library(data.table)

gtf_path <- '$GTF_FILE'
anno_path <- '$ANNOT_PATH'

if (!file.exists(anno_path)) {
  cat('Create riboWaltz annotation...\n')
  stopifnot(file.exists(gtf_path))
  anno_dt <- create_annotation(gtf = gtf_path)
  fwrite(anno_dt, anno_path, sep='\t')
} else {
  anno_dt <- fread(anno_path)
}

annot_out <- '$CLASSIFY_ANNOT'
fwrite(anno_dt[, .(transcript, l_tr, l_utr5, l_cds, l_utr3)],
       annot_out, sep='\t')
cat('Annotation saved: ', annot_out, '\n')
  "
else
  echo "Annotation exists, skip"
fi

echo -e "\nFilter BAM: keep record with minimal FLAG per QNAME"
for prefix in "${!sample_files[@]}"; do
    [ -z "$prefix" ] && continue
    BAM="$OUTPUT_DIR/${prefix}_genome_Aligned.toTranscriptome.out.bam"
    [ ! -f "$BAM" ] && continue

    echo "Filter sample: $prefix"
    samtools view -h "$BAM" | awk '
    /^@/ { print; next }
    {
        if (!best[$1] || $2 < best_flag[$1]) {
            best[$1] = $0
            best_flag[$1] = $2
        }
    }
    END { for (q in best) print best[q] }' | \
    samtools view -Sb -@ $THREADS - > "${BAM}.tmp" && mv -f "${BAM}.tmp" "$BAM"
    echo "Overwritten best-hit BAM: $BAM"
done

##############################################################################
# Step4 cont: 5-region stats
##############################################################################
echo -e "\n========================================"
echo "5-region stats (best hit)"
echo "================================================\n"

export LC_ALL=C
ANNOT_FILE="$CLASSIFY_ANNOT"

for prefix in "${!sample_files[@]}"; do
    [ -z "$prefix" ] && continue
    BAM="$OUTPUT_DIR/${prefix}_genome_Aligned.toTranscriptome.out.bam"
    [ ! -f "$BAM" ] && continue
    [ ! -f "$ANNOT_FILE" ] && { echo "ERROR: annotation missing"; exit 1; }

    echo -e "\nSample: $prefix"

stats=$(samtools view -F4 "$BAM" 2>/dev/null | awk '
BEGIN { pure5=0; border5c=0; purecds=0; borderc3=0; pure3=0; total=0; }
NR == FNR {
    sub(/_[0-9]+$/, "", $1);
    utr5[$1]=$3;
    cds_len[$1]=$4;
    next;
}
{
    tx=$3;
    sub(/_[0-9]+$/, "", tx);
    start=$4;
    rlen=length($10);
    end=start+rlen-1;
    if (!utr5[tx]) next;
    u5=utr5[tx];
    cds_end=u5+cds_len[tx];
    total++;
    if (end <= u5) pure5++;
    else if (start>u5 && end<=cds_end) purecds++;
    else if (start>cds_end) pure3++;
    else if (start<=u5 && end>u5) border5c++;
    else if (start<=cds_end && end>cds_end) borderc3++;
}
END { print pure5, border5c, purecds, borderc3, pure3, total }
' "$ANNOT_FILE" -)

    read pure5 border5c purecds borderc3 pure3 total <<< "$stats"

    p_pure5=$(echo "scale=2; $pure5*100/$total" | bc)
    p_border5c=$(echo "scale=2; $border5c*100/$total" | bc)
    p_purecds=$(echo "scale=2; $purecds*100/$total" | bc)
    p_borderc3=$(echo "scale=2; $borderc3*100/$total" | bc)
    p_pure3=$(echo "scale=2; $pure3*100/$total" | bc)
    genome_uniq=$(cat "$OUTPUT_DIR/${prefix}.genome_uniq.tmp")
    convert_rate=$(echo "scale=2; $total * 100 / $genome_uniq" | bc)

    echo "1. pure 5UTR      : $pure5 (${p_pure5}%)"
    echo "2. 5UTR-CDS border: $border5c (${p_border5c}%)"
    echo "3. pure CDS       : $purecds (${p_purecds}%)"
    echo "4. CDS-3UTR border: $borderc3 (${p_borderc3}%)"
    echo "5. pure 3UTR      : $pure3 (${p_pure3}%)"
    echo "Total valid reads : $total"
    echo "Genome unique mapped : $genome_uniq"
    echo "Genome -> transcript conversion rate : ${convert_rate}%"

    echo -e "\n【Sample $prefix 5-region stats】" >> "$STATS_REPORT"
    echo "pure 5UTR      : $pure5 (${p_pure5}%)" >> "$STATS_REPORT"
    echo "5UTR-CDS border: $border5c (${p_border5c}%)" >> "$STATS_REPORT"
    echo "pure CDS       : $purecds (${p_purecds}%)" >> "$STATS_REPORT"
    echo "CDS-3UTR border: $borderc3 (${p_borderc3}%)" >> "$STATS_REPORT"
    echo "pure 3UTR      : $pure3 (${p_pure3}%)" >> "$STATS_REPORT"
    echo "Total valid reads : $total" >> "$STATS_REPORT"
    echo "Genome unique mapped reads : $genome_uniq" >> "$STATS_REPORT"
    echo "Genome -> transcript conversion rate : ${convert_rate}%" >> "$STATS_REPORT"
    echo -e "----------------------------------------" >> "$STATS_REPORT"
done

##############################################################################
# Step5: Coverage & RPF length distribution
##############################################################################
echo -e "\n========================================"
echo "Step5: Transcript coverage and RPF length stats"
echo "========================================"
echo -e "\n【Coverage & RPF stats】" >> "$STATS_REPORT"

for prefix in "${!sample_files[@]}"; do
    if [ -z "$prefix" ]; then continue; fi
    process_type="${sample_process_type[$prefix]}"
    echo -e "\nSample: $prefix (type: $process_type)"

    transcript_bam="$OUTPUT_DIR/${prefix}_genome_Aligned.toTranscriptome.out.bam"
    if [ ! -f "$transcript_bam" ]; then echo "Skip: BAM missing"; continue; fi

    total_mapped_nucleotides=$(samtools view -F4 -@ $THREADS "$transcript_bam" | awk '{sum+=length($10)}END{print sum+0}')
    bam_header_dict="$OUTPUT_DIR/${prefix}_transcript_header.dict"
    samtools view -H "$transcript_bam" | grep ^@SQ | awk -F'\t' '{split($2,sn,":");split($3,ln,":");print sn[2]"\t"ln[2]}' > "$bam_header_dict"

    expressed_transcripts_total_length=$(samtools view -F4 -@ $THREADS "$transcript_bam" | cut -f3 | sort -u | awk -F'\t' -v dict="$bam_header_dict" '
    BEGIN{while((getline<dict)>0)len_dict[$1]=$2;close(dict)}{total+=len_dict[$1]}END{print total+0}')
    expressed_transcripts_count=$(samtools view -F4 -@ $THREADS "$transcript_bam" | cut -f3 | sort -u | wc -l)

    if [ "$expressed_transcripts_total_length" -eq 0 ] || [ "$total_mapped_nucleotides" -eq 0 ]; then
        global_coverage="0.00"
    else
        global_coverage=$(echo "scale=2; $total_mapped_nucleotides / $expressed_transcripts_total_length" | bc)
    fi

    rpf_length_stats=$(samtools view -F4 -@ $THREADS "$transcript_bam" | awk '{print length($10)}' | sort -n | awk '
    BEGIN{min_len=1000;max_len=0;}
    {len=$1;count[len]++;total_reads++;if(len<min_len)min_len=len;if(len>max_len)max_len=len;}
    END{printf "%d,%d,%d,",min_len,max_len,total_reads;for(l in count)printf "%d:%.2f%%,",l,count[l]/total_reads*100;}')
    IFS=',' read -r min_rpf_len max_rpf_len total_rpf_reads rpf_distribution <<< "$rpf_length_stats"
    rpf_distribution=${rpf_distribution%,}

    rpf_temp_tsv="$OUTPUT_DIR/${prefix}_rpf_temp_stats.tsv"
    samtools view -F4 -@ $THREADS "$transcript_bam" | awk '{print length($10)}' | sort -n | uniq -c | awk -v total=$total_rpf_reads '{print $2"\t"$1"\t"sprintf("%.2f",$1/total*100)}' > $rpf_temp_tsv
    Rscript -e "
    library(ggplot2)
    df <- read.delim('$rpf_temp_tsv',header=F,col.names=c('Length','Count','Percent'))
    png('$OUTPUT_DIR/${prefix}_RPF_length_distribution.png',width=1000,height=600,res=150)
    p <- ggplot(df,aes(x=as.factor(Length),y=Count,fill=as.factor(Length)))+
         geom_col(alpha=0.8,color='black')+geom_text(aes(label=paste0(Percent,'%')),vjust=-0.2,size=2)+
         labs(title=paste0('$prefix RPF Length Distribution (Best Hit)'),x='RPF Length',y='Count')+theme_bw()+theme(legend.position='none')
    print(p);dev.off()"
    rm -f $rpf_temp_tsv

    echo "Generate CDS-only RPF length..."
    cds_rpf_tsv="$OUTPUT_DIR/${prefix}_CDS_RPF.tsv"

samtools view -F4 "$transcript_bam" 2>/dev/null | awk '
BEGIN {
    while ((getline < "'"$CLASSIFY_ANNOT"'") > 0) {
        sub(/_[0-9]+$/, "", $1);
        utr5[$1]=$3;
        cds_len[$1]=$4;
    }
}
{
    tx=$3;
    sub(/_[0-9]+$/, "", tx);
    start=$4;
    rlen=length($10);
    end=start+rlen-1;
    if (!utr5[tx]) next;
    u5=utr5[tx];
    cds_end=u5+cds_len[tx];
    if (start>u5 && end<=cds_end) print length($10);
}' | sort -n | uniq -c | awk '{print $2"\t"$1}' > $cds_rpf_tsv

    Rscript --no-environ -e "
    library(ggplot2)
    df <- read.delim('$cds_rpf_tsv',header=F,col.names=c('Length','Count'))
    if(nrow(df)>0){
        df\$Percent = sprintf('%.2f', df\$Count/sum(df\$Count)*100)
        p <- ggplot(df,aes(x=as.factor(Length),y=Count,fill=as.factor(Length)))+
             geom_col(alpha=0.8,color='black')+
             geom_text(aes(label=Percent),vjust=-0.2,size=2)+
             labs(title=paste0('$prefix CDS-only RPF Length (Best Hit)'),x='RPF Length',y='Count')+theme_bw()+theme(legend.position='none')
    }else{
        p <- ggplot() +
             annotate('text',x=1,y=1,label='No CDS RPF Reads Detected',size=6)+
             labs(title=paste0('$prefix CDS-only RPF Length (Best Hit)'),x='RPF Length',y='Count')+theme_bw()
    }
    png('$OUTPUT_DIR/${prefix}_CDS_RPF_distribution.png',width=1000,height=600,res=150)
    print(p);dev.off()"
    rm -f $cds_rpf_tsv

    echo -e "\n【Coverage】"
    echo "Expressed transcripts: $expressed_transcripts_count"
    echo "Total transcript length: $expressed_transcripts_total_length bp"
    echo "Total mapped bases: $total_mapped_nucleotides bp"
    echo "Global coverage: $global_coverage ×"
    echo -e "\n【RPF length】"
    echo "Total RPF: $total_rpf_reads"
    echo "Range: $min_rpf_len ~ $max_rpf_len bp"
    echo "Distribution: $rpf_distribution"
    echo -e "\nOutput plots:"
    echo "All transcript RPF: ${prefix}_RPF_length_distribution.png"
    echo "CDS-only RPF: ${prefix}_CDS_RPF_distribution.png"

    echo -e "\n【Sample $prefix coverage & RPF stats】" >> "$STATS_REPORT"
    echo "Process type: $process_type" >> "$STATS_REPORT"
    echo "Expressed transcripts count: $expressed_transcripts_count" >> "$STATS_REPORT"
    echo "Total transcript length: $expressed_transcripts_total_length bp" >> "$STATS_REPORT"
    echo "Total mapped bases: $total_mapped_nucleotides bp" >> "$STATS_REPORT"
    echo "Global coverage: $global_coverage ×" >> "$STATS_REPORT"
    echo "Total RPF reads: $total_rpf_reads" >> "$STATS_REPORT"
    echo "RPF length range: $min_rpf_len - $max_rpf_len bp" >> "$STATS_REPORT"
    echo "RPF length distribution: $rpf_distribution" >> "$STATS_REPORT"
    echo "All transcript RPF plot: ${prefix}_RPF_length_distribution.png" >> "$STATS_REPORT"
    echo "CDS-only RPF plot: ${prefix}_CDS_RPF_distribution.png" >> "$STATS_REPORT"
    echo -e "----------------------------------------" >> "$STATS_REPORT"

    rm -f "$bam_header_dict"
done

##############################################################################
# Final summary
##############################################################################
echo -e "\n========================================"
echo "Pipeline finished"
echo "========================================"
echo -e "\nResults: $OUTPUT_DIR"
echo "Report: $STATS_REPORT"
echo -e "\nReport preview:"
cat "$STATS_REPORT" | grep -E "Sample|5UTR|CDS|Total valid reads|Global coverage|RPF" | head -70
echo "...full report in $STATS_REPORT"
