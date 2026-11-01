#!/bin/bash
# ribo_LAN: Global executable ribo-seq analysis command-line tool
# Usage: ribo_LAN --step1 <CSV name/path> | --step2 | --step3 | --step4 | --step5 | --step6 <group CSV name> <gene list txt name>
# Note: For --step1/--step6, CSV/txt basenames are automatically looked up under ./1.raw_fastq/ first
# Final version: step3 supports Monosome/Disome dual workflow; Disome trim length + RPF length must be manually entered
# Only modified step1: add Monosome/Disome selection, filter length (21-40/21-80) passed to step1 script
# New step6: codon occupancy analysis (reuses species/workflow/length/trim parameters from step3; group CSV in 1.raw_fastq)
# Revision log: Disome workflow changed from 3'-end trimming to 5'-end trimming
# New: step6 accepts gene list txt for gene-specific footprint density plot (no binning)

# ==================== Fixed path configuration (modify according to your environment!) ====================
CONDA_ENV_NAME="ribo_simple"                # conda environment name
PROJECT_TOOLS_DIR="$HOME/ribo_simple_LAN/TOOLS"  # Directory storing analysis scripts
# Absolute paths for each analysis script
STEP1_SCRIPT="$PROJECT_TOOLS_DIR/step1_filter_clean.sh"
STEP2_SCRIPT="$PROJECT_TOOLS_DIR/step2_align.sh"
STEP3_SCRIPT="$PROJECT_TOOLS_DIR/step3_riboWaltz.R"
STEP4_SCRIPT="$PROJECT_TOOLS_DIR/step4_A.P.E_around_motif.R"
STEP5_SCRIPT="$PROJECT_TOOLS_DIR/step5_personalize.R"
# ========== New: step6 script path ==========
STEP6_SCRIPT="$PROJECT_TOOLS_DIR/step6_ribo_group_diff.R"

# ==================== New for step1: default fragment filter length global variables ====================
MONO_MIN=21  # Monosome default minimum filter length
MONO_MAX=40  # Monosome default maximum filter length
DISO_MIN=21  # Disome default minimum filter length
DISO_MAX=80  # Disome default maximum filter length

# ==================== Environment check ====================
# Check whether target conda environment is activated
if [ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV_NAME" ]; then
    echo "❌ ERROR: Target conda environment '$CONDA_ENV_NAME' is NOT activated!"
    echo "Please run first: source $HOME/ribo_simple_LAN/ribo_LAN_env.sh"
    exit 1
fi

# Function to check existence of script file
check_script_exists() {
    if [ ! -f "$1" ]; then
        echo "❌ ERROR: Script not found! Path: $1"
        exit 1
    fi
}

# Function to check R environment
check_r_env() {
    if ! command -v Rscript &> /dev/null; then
        echo "❌ ERROR: Rscript not found! Please confirm R is installed in conda env '$CONDA_ENV_NAME'"
        exit 1
    fi
}

# Function to resolve file path (CSV/TXT, ★ priority lookup under 1.raw_fastq ★)
complete_file_path() {
    local input="$1"
    local work_dir=$(pwd)
    local raw_fastq_file="$work_dir/1.raw_fastq/${input}"
    local abs_path=""

    # 1. First try 1.raw_fastq regardless of whether input carries path
    if [ -f "$raw_fastq_file" ]; then
        abs_path=$(realpath "$raw_fastq_file")
    # 2. If not found in 1.raw_fastq, use input path directly
    elif [ -f "$input" ]; then
        abs_path=$(realpath "$input")
    fi

    # 3. Validate existence
    if [ -n "$abs_path" ]; then
        echo "$abs_path"
        return 0
    else
        echo ""
        return 1
    fi
}

# ==================== New for step1: workflow selection function (only called by step1) ====================
select_ribo_type_for_step1() {
    echo -e "\n====================================="
    echo "          Select workflow (step1)     "
    echo "====================================="
    echo "1) Monosome standard workflow (filter 21-40nt)"
    echo "2) Disome disome workflow (filter 21-80nt)"
    
    while true; do
        read -p "Enter choice (1/2): " RIBO_CHOICE
        case "$RIBO_CHOICE" in
            1)
                RIBO_TYPE="monosome"
                FILT_MIN=$MONO_MIN
                FILT_MAX=$MONO_MAX
                break
                ;;
            2)
                RIBO_TYPE="disome"
                FILT_MIN=$DISO_MIN
                FILT_MAX=$DISO_MAX
                break
                ;;
            *)
                echo "❌ ERROR: Invalid choice! Only 1 (Monosome) or 2 (Disome)"
                ;;
        esac
    done

    # Confirmation
    echo -e "\n⚠️  Selected workflow: ${RIBO_TYPE^^}, fragment filter range ${FILT_MIN}-${FILT_MAX}nt"
    read -p "Confirm to proceed? (y/n): " RIBO_CONFIRM
    if [[ "$RIBO_CONFIRM" != "y" && "$RIBO_CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
        exit 0
    fi

    echo -e "✅ Workflow selected for step1: ${RIBO_TYPE^^}, fragment filter: ${FILT_MIN}-${FILT_MAX}nt\n"
}

# ==================== Species selection (Removed CHmito) ====================
select_species() {
    echo -e "\n====================================="
    echo "          Select species              "
    echo "====================================="
    echo "1) human"
    echo "2) mouse"
    
    while true; do
        read -p "Enter choice (1/2): " SPECIES_CHOICE
        case "$SPECIES_CHOICE" in
            1)
                SPECIES="human"
                break
                ;;
            2)
                SPECIES="mouse"
                break
                ;;
            *)
                echo "❌ ERROR: Invalid choice! Only 1 (human) or 2 (mouse)"
                ;;
        esac
    done

    # Confirmation
    echo -e "\n⚠️  Reference files for species [${SPECIES}] will be used"
    read -p "Confirm to proceed? (y/n): " CONFIRM
    if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
        exit 0
    fi

    echo -e "✅ Species selected: ${SPECIES}\n"
}

# Select Monosome / Disome workflow (shared by step3/step6)
select_ribo_type() {
    echo -e "\n====================================="
    echo "          Select ribosome workflow         "
    echo "====================================="
    echo "1) Monosome standard workflow (no trimming)"
    echo "2) Disome disome workflow (trim fixed length from 5' end of all RPFs)"
    
    while true; do
        read -p "Enter choice (1/2): " RIBO_CHOICE
        case "$RIBO_CHOICE" in
            1)
                RIBO_TYPE="monosome"
                TRIM_LENGTH=0  # Monosome trim length fixed to 0
                break
                ;;
            2)
                RIBO_TYPE="disome"
                TRIM_LENGTH=0  # Initialize to avoid empty variable
                break
                ;;
            *)
                echo "❌ ERROR: Invalid choice! Only 1 (Monosome) or 2 (Disome)"
                ;;
        esac
    done

    echo -e "\n⚠️  Selected workflow: ${RIBO_TYPE^^}"
    if [ "$RIBO_TYPE" == "disome" ];then
        echo "⚠️  Disome note: 5' trim length and RPF length parameters must be manually provided, no default"
    fi
    read -p "Confirm this workflow? (y/n): " RIBO_CONFIRM
    if [[ "$RIBO_CONFIRM" != "y" && "$RIBO_CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
    fi

    echo -e "✅ Workflow selected: ${RIBO_TYPE^^}\n"
}

# Disome only - interactive 5' trim length input (nt)
select_trim_length() {
    if [ "$RIBO_TYPE" != "disome" ]; then
        return 0
    fi

    echo -e "\n====================================="
    echo "        Disome - trim length setting        "
    echo "====================================="
    echo "💡 Recommendation: typical disome uses 33nt trim, adjust based on your experiment"
    echo "💡 Note: trimming from 5' end, retained 3' sequence for downstream analysis"
    
    while true; do
        read -p "Input 5' trim length in nt (positive integer required): " TRIM_LENGTH
        if [ -z "$TRIM_LENGTH" ]; then
            echo "❌ ERROR: Trim length cannot be empty!"
        elif [[ "$TRIM_LENGTH" =~ ^[0-9]+$ && "$TRIM_LENGTH" -gt 0 ]]; then
            break
        else
            echo "❌ ERROR: Must enter positive integer!"
        fi
    done

    echo -e "\n⚠️  Confirm: Disome will trim ${TRIM_LENGTH}nt from 5' end of all RPFs"
    read -p "Confirm trim length? (y/n): " TRIM_CONFIRM
    if [[ "$TRIM_CONFIRM" != "y" && "$TRIM_CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
        exit 0
    fi

    echo -e "✅ Confirmed Disome trim length: ${TRIM_LENGTH}nt (5'-end)\n"
}

# Shared by step3/step6 - interactive RPF length parameters
# $1: workflow type (monosome/disome)
select_rpf_length() {
    local ribo_type=$1
    local prompt_prefix=""
    
    if [ "$ribo_type" == "disome" ]; then
        prompt_prefix="[Disome - no default]"
        echo -e "\n====================================="
        echo "  ${prompt_prefix} Input RPF length after 5' ${TRIM_LENGTH}nt trim  "
        echo "  Check RPF distribution in 3. alignments/Report.txt after trimming!  "
        echo "====================================="
    else
        prompt_prefix="[Monosome - default supported]"
        echo -e "\n====================================="
        echo "  ${prompt_prefix} Check RPF distribution in 3. alignments/Report.txt  "
        echo "  1. Preliminary (default 26-34nt) 2. Strict high-depth (default 28-32nt)  "
        echo "====================================="
    fi
    
    # 1. Min length, preliminary
    while true; do
        if [ "$ribo_type" == "monosome" ]; then
            read -p "${prompt_prefix} Preliminary min RPF length (default 26): " MIN_LEN_DEFAULT
            if [ -z "$MIN_LEN_DEFAULT" ]; then
                MIN_LEN_DEFAULT=26
                break
            fi
        else
            read -p "${prompt_prefix} Preliminary min RPF length after trim (REQUIRED): " MIN_LEN_DEFAULT
            if [ -z "$MIN_LEN_DEFAULT" ]; then
                echo "❌ ERROR: Disome requires manual length input, cannot be empty!"
                continue
            fi
        fi
        if [[ "$MIN_LEN_DEFAULT" =~ ^[0-9]+$ ]]; then
            break
        else
            echo "❌ ERROR: Enter positive integer!"
        fi
    done

    # 2. Max length, preliminary
    while true; do
        if [ "$ribo_type" == "monosome" ]; then
            read -p "${prompt_prefix} Preliminary max RPF length (default 34): " MAX_LEN_DEFAULT
            if [ -z "$MAX_LEN_DEFAULT" ]; then
                MAX_LEN_DEFAULT=34
                break
            fi
        else
            read -p "${prompt_prefix} Preliminary max RPF length after trim (REQUIRED): " MAX_LEN_DEFAULT
            if [ -z "$MAX_LEN_DEFAULT" ]; then
                echo "❌ ERROR: Disome requires manual length input, cannot be empty!"
                continue
            fi
        fi
        if [[ "$MAX_LEN_DEFAULT" =~ ^[0-9]+$ ]]; then
            if [ "$MAX_LEN_DEFAULT" -ge "$MIN_LEN_DEFAULT" ]; then
                break
            else
                echo "❌ ERROR: Max length cannot be smaller than min length ${MIN_LEN_DEFAULT}!"
            fi
        else
            echo "❌ ERROR: Enter positive integer!"
        fi
    done

    # 3. Min length, strict
    while true; do
        if [ "$ribo_type" == "monosome" ]; then
            read -p "${prompt_prefix} Strict min RPF length for high depth (default 28): " MIN_LEN_STRICT
            if [ -z "$MIN_LEN_STRICT" ]; then
                MIN_LEN_STRICT=28
                break
            fi
        else
            read -p "${prompt_prefix} Strict min RPF length after trim for high depth (REQUIRED): " MIN_LEN_STRICT
            if [ -z "$MIN_LEN_STRICT" ]; then
                echo "❌ ERROR: Disome requires manual length input, cannot be empty!"
                continue
            fi
        fi
        if [[ "$MIN_LEN_STRICT" =~ ^[0-9]+$ ]]; then
            if [ "$MIN_LEN_STRICT" -ge "$MIN_LEN_DEFAULT" ] && [ "$MIN_LEN_STRICT" -le "$MAX_LEN_DEFAULT" ]; then
                break
            else
                echo "❌ ERROR: Strict min length must be within [${MIN_LEN_DEFAULT}, ${MAX_LEN_DEFAULT}]!"
            fi
        else
            echo "❌ ERROR: Enter positive integer!"
        fi
    done

    # 4. Max length, strict
    while true; do
        if [ "$ribo_type" == "monosome" ]; then
            read -p "${prompt_prefix} Strict max RPF length for high depth (default 32): " MAX_LEN_STRICT
            if [ -z "$MAX_LEN_STRICT" ]; then
                MAX_LEN_STRICT=32
                break
            fi
        else
            read -p "${prompt_prefix} Strict max RPF length after trim for high depth (REQUIRED): " MAX_LEN_STRICT
            if [ -z "$MAX_LEN_STRICT" ]; then
                echo "❌ ERROR: Disome requires manual length input, cannot be empty!"
                continue
            fi
        fi
        if [[ "$MAX_LEN_STRICT" =~ ^[0-9]+$ ]]; then
            if [ "$MAX_LEN_STRICT" -ge "$MIN_LEN_STRICT" ] && [ "$MAX_LEN_STRICT" -le "$MAX_LEN_DEFAULT" ]; then
                break
            else
                echo "❌ ERROR: Strict max length must be within [${MIN_LEN_STRICT}, ${MAX_LEN_DEFAULT}]!"
            fi
        else
            echo "❌ ERROR: Enter positive integer!"
        fi
    done

    # 5. Confirm all length parameters
    echo -e "\n⚠️  Length parameter summary:"
    echo "   Preliminary range: ${MIN_LEN_DEFAULT}-${MAX_LEN_DEFAULT}nt"
    echo "   Strict high-depth range: ${MIN_LEN_STRICT}-${MAX_LEN_STRICT}nt"
    if [ "$ribo_type" == "disome" ]; then
        echo "   ⚠️  These lengths are AFTER 5' trimming ${TRIM_LENGTH}nt!"
    fi
    read -p "Confirm these parameters? (y/n): " LEN_CONFIRM
    if [[ "$LEN_CONFIRM" != "y" && "$LEN_CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
        exit 0
    fi

    echo -e "✅ RPF length parameters confirmed\n"
}

# ==================== Argument parsing & execution ====================
if [ $# -eq 0 ]; then
    echo "📚 Usage: ribo_LAN --step1 <CSV name> | --step2 | --step3 | --step4 | --step5 | --step6 <group CSV name> <gene list txt name>"
    echo "💡 Note: For --step1/--step6, only basename needed; automatically search ./1.raw_fastq/"
    echo "💡 --step1: Monosome(21-40nt)/Disome(21-80nt) selection, filter lengths passed to step1 script"
    echo "💡 step3/step6: Monosome supports default lengths; Disome requires manual 5'-trim length + RPF lengths"
    echo "💡 --step6: codon occupancy analysis, reuses step3 parameters; group CSV columns: sample_name/group_name, comma separated; gene list txt one gene per line"
    exit 1
fi

case "$1" in
    --step1)
        if [ $# -ne 2 ]; then
            echo "❌ ERROR: --step1 requires CSV basename!"
            echo "📚 Usage: ribo_LAN --step1 <CSV name>"
            echo "💡 Example: ribo_LAN --step1 108T.csv (auto lookup under 1.raw_fastq/)"
            exit 1
        fi
        
        CSV_INPUT="$2"
        CSV_PATH=$(complete_file_path "$CSV_INPUT")
        
        if [ -z "$CSV_PATH" ]; then
            echo "❌ ERROR: ${CSV_INPUT} not found under 1.raw_fastq!"
            echo "🔍 Lookup path: $(pwd)/1.raw_fastq/${CSV_INPUT}"
            exit 1
        fi
        
        echo "✅ CSV found: ${CSV_PATH}"

        select_ribo_type_for_step1
        
        check_script_exists "$STEP1_SCRIPT"
        echo -e "\n🚀 Running step1: Preprocessing (adapter trimming - QC - tRNA/rRNA removal)【${RIBO_TYPE^^}】($(date))\n"
        "$STEP1_SCRIPT" "$CSV_PATH" "$FILT_MIN" "$FILT_MAX"
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step1 finished ($(date))\n"
        else
            echo -e "\n❌ step1 failed ($(date))\n"
            exit 1
        fi
        ;;
    --step2)
        check_script_exists "$STEP2_SCRIPT"
        echo -e "\n🚀 Running step2: Read alignment (STAR) ($(date))\n"
        "$STEP2_SCRIPT"
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step2 finished ($(date))\n"
        else
            echo -e "\n❌ step2 failed ($(date))\n"
            exit 1
        fi
        ;;
    --step3)
        check_r_env
        check_script_exists "$STEP3_SCRIPT"
        select_species
        select_ribo_type
        select_trim_length
        select_rpf_length "$RIBO_TYPE"
        echo -e "\n🚀 Running step3: riboWaltz P-site calling (A/E-site included)【${RIBO_TYPE^^}】($(date))\n"
        Rscript "$STEP3_SCRIPT" "$SPECIES" "$MIN_LEN_DEFAULT" "$MAX_LEN_DEFAULT" "$MIN_LEN_STRICT" "$MAX_LEN_STRICT" "$RIBO_TYPE" "$TRIM_LENGTH"
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step3 finished ($(date))\n"
        else
            echo -e "\n❌ step3 failed ($(date))\n"
            exit 1
        fi
        ;;
     --step4)
        check_r_env
        check_script_exists "$STEP4_SCRIPT"
        select_species
        echo -e "\n🚀 Running step4: Motif analysis around A/P/E-sites ($(date))\n"
        Rscript "$STEP4_SCRIPT" "$SPECIES"
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step4 finished ($(date))\n"
        else
            echo -e "\n❌ step4 failed ($(date))\n"
            exit 1
        fi
        ;;
    --step5)
        check_r_env
        check_script_exists "$STEP5_SCRIPT"
        select_species
        echo -e "\n🚀 Running step5: Ribo-seq translation efficiency analysis (P-site based) ($(date))\n"
        Rscript "$STEP5_SCRIPT" "$SPECIES"
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step5 finished ($(date))\n"
        else
            echo -e "\n❌ step5 failed ($(date))\n"
            exit 1
        fi
        ;;
    --step6)
        if [ $# -ne 3 ]; then
            echo "❌ ERROR: --step6 requires group CSV basename AND gene list txt basename!"
            echo "📚 Usage: ribo_LAN --step6 <group CSV name> <gene list txt name>"
            echo "💡 Example: ribo_LAN --step6 group.csv list.txt (auto lookup under 1.raw_fastq/)"
            echo "💡 Requirement: 1) CSV columns sample_name/group_name, comma separated; 2) TXT one gene per line, no header"
            exit 1
        fi
        
        GROUP_CSV_INPUT="$2"
        GROUP_CSV_PATH=$(complete_file_path "$GROUP_CSV_INPUT")
        if [ -z "$GROUP_CSV_PATH" ]; then
            echo "❌ ERROR: Group CSV ${GROUP_CSV_INPUT} not found under 1.raw_fastq!"
            echo "🔍 Lookup path: $(pwd)/1.raw_fastq/${GROUP_CSV_INPUT}"
            exit 1
        fi
        echo "✅ Group CSV found: ${GROUP_CSV_PATH}"

        GENE_LIST_INPUT="$3"
        GENE_LIST_PATH=$(complete_file_path "$GENE_LIST_INPUT")
        if [ -z "$GENE_LIST_PATH" ]; then
            echo "❌ ERROR: Gene list ${GENE_LIST_INPUT} not found under 1.raw_fastq!"
            echo "🔍 Lookup path: $(pwd)/1.raw_fastq/${GENE_LIST_INPUT}"
            exit 1
        fi
        echo "✅ Gene list TXT found: ${GENE_LIST_PATH}"

        check_r_env
        check_script_exists "$STEP6_SCRIPT"
        select_species
        select_ribo_type
        select_trim_length
        select_rpf_length "$RIBO_TYPE"

        echo -e "\n🚀 Running step6: Codon occupancy analysis (PCA + KS test + footprint density + gene-specific plot)【${RIBO_TYPE^^}】($(date))\n"
        Rscript "$STEP6_SCRIPT" "$SPECIES" "$MIN_LEN_DEFAULT" "$MAX_LEN_DEFAULT" "$MIN_LEN_STRICT" "$MAX_LEN_STRICT" "$RIBO_TYPE" "$TRIM_LENGTH" "$GROUP_CSV_PATH" "$GENE_LIST_PATH"
        
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step6 finished ($(date))\n"
        else
            echo -e "\n❌ step6 failed ($(date))\n"
            exit 1
        fi
        ;;
    *)
        echo "❌ ERROR: Invalid argument!"
        echo "📚 Usage: ribo_LAN --step1 <CSV name> | --step2 | --step3 | --step4 | --step5 | --step6 <group CSV name> <gene list txt name>"
        echo "💡 Note: For --step1/--step6, only basename needed; automatically search ./1.raw_fastq/"
        exit 1
        ;;
esac
