#!/bin/bash
# rna_LAN: Global executable RNA-seq analysis command-line tool
# Usage: rna_LAN --step1 <CSV name/path> | --step2 | --step3 <group CSV name/path> | --step4 | --step5 <group CSV name/path>
# Note: When only a basename is supplied for --step1/--step3/--step5, the script first looks under ./1.raw_fastq/
# Feature: Interactive species selection in bash layer for step4; choose TE normalization method + interactive Ribo-seq root directory
# Consistency: Identical argument passing logic for step3/step5 and step1; basenames auto-resolved, full paths accepted directly
# Compatibility: Fully compatible with the latest RNAstep4_TE.R + RNAstep5_TE_diff_analysis.R
# Fix: Resolve tilde (~) path parsing bug, supports all path formats

# ==================== Fixed path configuration (modify according to your environment!) ====================
CONDA_ENV_NAME="ribo_simple"                # Reuse ribo-seq conda environment (can be changed)
PROJECT_TOOLS_DIR="$HOME/rna_simple_LAN/TOOLS"  # Directory storing RNA-seq scripts
# Absolute paths for analysis scripts
STEP1_SCRIPT="$PROJECT_TOOLS_DIR/RNAstep1_filter_clean.sh"
STEP2_SCRIPT="$PROJECT_TOOLS_DIR/RNAstep2_align.sh"
STEP3_SCRIPT="$PROJECT_TOOLS_DIR/RNAstep3_countPCA.R"
STEP4_SCRIPT="$PROJECT_TOOLS_DIR/RNAstep4_TE.R"  # TE calculation script
STEP5_SCRIPT="$PROJECT_TOOLS_DIR/RNAstep5_TE_diff_analysis.R" # TE differential analysis script

# ==================== Environment check ====================
# Check whether target conda environment is activated
if [ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV_NAME" ]; then
    echo "❌ ERROR: Target conda environment '$CONDA_ENV_NAME' is NOT activated!"
    echo "Please run first: source $HOME/rna_simple_LAN/rna_LAN_env.sh"
    exit 1
fi

# Function to check existence of script file
check_script_exists() {
    if [ ! -f "$1" ]; then
        echo "❌ ERROR: Script not found! Path: $1"
        exit 1
    fi
}

# Check Rscript availability for R scripts
check_r_env() {
    if ! command -v Rscript &> /dev/null; then
        echo "❌ ERROR: Rscript not found! Please install R inside the conda environment."
        exit 1
    fi
}

# Core helper: expand tilde paths
expand_tilde_path() {
    local input_path="$1"
    echo "${input_path/#\~/$HOME}"
}

# Resolve CSV path (shared by step1/step3/step5), tilde-aware
complete_csv_path() {
    local input="$1"
    local work_dir=$(pwd)
    local abs_path=""
    input=$(expand_tilde_path "$input")

    # Pure basename → look under 1.raw_fastq
    if [[ "$input" != */* ]]; then
        abs_path=$(realpath "${work_dir}/1.raw_fastq/${input}" 2>/dev/null)
    # Full path supplied
    else
        abs_path=$(realpath "${input}" 2>/dev/null)
    fi

    if [ -f "$abs_path" ]; then
        echo "$abs_path"
        return 0
    else
        echo ""
        return 1
    fi
}

# ==================== Interactive selection functions ====================
# Select species (human / mouse only, consistent with ribo_LAN)
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

    echo -e "\n⚠️  Reference files for species [${SPECIES}] will be used"
    read -p "Confirm to proceed? (y/n): " CONFIRM
    if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
        exit 0
    fi

    echo -e "✅ Species selected: ${SPECIES}\n"
}

# Step4 only: select TE normalization method
select_te_normalization() {
    echo -e "\n====================================="
    echo "        Select TE normalization method      "
    echo "====================================="
    echo "1) RPM (Reads Per Million, recommended)"
    echo "2) TPM (Transcripts Per Million)"
    echo "3) VST (DESeq2 variance stabilizing transformation)"
    
    while true; do
        read -p "Enter choice (1/2/3): " NORM_CHOICE
        case "$NORM_CHOICE" in
            1)
                NORM_METHOD="RPM"
                break
                ;;
            2)
                NORM_METHOD="TPM"
                break
                ;;
            3)
                NORM_METHOD="VST"
                break
                ;;
            *)
                echo "❌ ERROR: Invalid choice! Only 1/2/3"
                ;;
        esac
    done

    echo -e "\n⚠️  TE will be calculated using [${NORM_METHOD}]"
    read -p "Confirm to proceed? (y/n): " NORM_CONFIRM
    if [[ "$NORM_CONFIRM" != "y" && "$NORM_CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
        exit 0
    fi

    echo -e "✅ TE normalization method selected: ${NORM_METHOD}\n"
}

# Step4 only: interactive Ribo-seq root directory input & validation, tilde-aware
select_ribo_seq_dir() {
    echo -e "\n====================================="
    echo "      Enter Ribo-seq analysis root directory     "
    echo "====================================="
    echo "💡 Requirement: folder must contain subdirectory 4.riboWaltz_output"
    echo "💡 Supports: absolute path / relative path / ~ tilde path / current directory (./)"
    echo "💡 Example: ~/ribo_simple_LAN/YourData/FeiLan_25.12.12  or  /home/user/xxx  or  ./xxx\n"
    
    local ribo_dir=""
    while true; do
        read -p "Enter Ribo-seq root path: " ribo_dir
        ribo_dir=$(echo "$ribo_dir" | sed 's/^[ \t]*//;s/[ \t]*$//')
        ribo_dir=$(expand_tilde_path "$ribo_dir")
        local abs_ribo_dir=$(realpath "${ribo_dir}" 2>/dev/null)
        
        if [ ! -d "$abs_ribo_dir" ]; then
            echo "❌ ERROR: Path does not exist! Please re-enter."
            continue
        fi
        
        local target_subdir="${abs_ribo_dir}/4.riboWaltz_output/2.translational expression"
        if [ -d "$target_subdir" ]; then
            RIBO_SEQ_DIR="$abs_ribo_dir"
            break
        else
            echo "❌ ERROR: Required subdirectory not found!"
            echo "🔍 Missing path: $target_subdir"
            echo "💡 Make sure you supply the Ribo-seq root folder (e.g. FeiLan_25.12.12), NOT 4.riboWaltz_output!"
        fi
    done

    echo -e "\n⚠️  Ribo-seq directory to use: ${RIBO_SEQ_DIR}"
    read -p "Confirm this path? (y/n): " RIBO_CONFIRM
    if [[ "$RIBO_CONFIRM" != "y" && "$RIBO_CONFIRM" != "Y" ]]; then
        echo -e "\n🔴 User cancelled, exit\n"
        exit 0
    fi

    echo -e "✅ Ribo-seq directory selected: ${RIBO_SEQ_DIR}\n"
}

# Step5 only: auto-locate TE result file produced by step4
find_te_result_file() {
    echo -e "\n====================================="
    echo "        Auto-locate TE result file        "
    echo "====================================="
    local work_dir=$(pwd)
    local te_result_root="${work_dir}/3.rna_seq_analysis/TE_analysis"
    local te_file_path=""

    for norm_type in RPM TPM VST; do
        local curr_te_file="${te_result_root}/TE_${norm_type}/TE_results_ALL_samples_${norm_type}.tsv"
        if [ -f "$curr_te_file" ]; then
            te_file_path="$curr_te_file"
            NORM_METHOD_USED="$norm_type"
            break
        fi
    done

    if [ -z "$te_file_path" ]; then
        echo "❌ ERROR: TE result file from step4 not found!"
        echo "🔍 Search root: ${te_result_root}"
        echo "💡 Run rna_LAN --step4 first before step5!"
        exit 1
    else
        te_file_path=$(realpath "$te_file_path")
        echo "✅ TE file located: ${te_file_path}"
        echo "✅ Detected normalization method: ${NORM_METHOD_USED}"
        TE_RESULT_FILE="$te_file_path"
    fi
}

# ==================== Argument parsing & execution ====================
if [ $# -eq 0 ]; then
    echo "📚 Usage: rna_LAN --step1 <CSV name/path> | --step2 | --step3 <group CSV name/path> | --step4 | --step5 <group CSV name/path>"
    echo "💡 Note: For --step1/--step3/--step5, basenames without / are automatically searched under ./1.raw_fastq/"
    exit 1
fi

case "$1" in
    --step1)
        if [ $# -ne 2 ]; then
            echo "❌ ERROR: --step1 requires CSV name/path!"
            echo "📚 Usage: rna_LAN --step1 <CSV name/path>"
            echo "💡 Example1 (basename): rna_LAN --step1 108T.csv"
            echo "💡 Example2 (path): rna_LAN --step1 /home/data/108T.csv or ~/data/108T.csv"
            exit 1
        fi
        
        CSV_INPUT="$2"
        CSV_PATH=$(complete_csv_path "$CSV_INPUT")
        
        if [ -z "$CSV_PATH" ]; then
            echo "❌ ERROR: CSV file not found!"
            if [[ "$CSV_INPUT" != */* ]]; then
                echo "🔍 Searched: $(pwd)/1.raw_fastq/${CSV_INPUT}"
            else
                echo "🔍 Searched: ${CSV_INPUT}"
            fi
            exit 1
        fi
        
        echo "✅ CSV found: ${CSV_PATH}"
        
        check_script_exists "$STEP1_SCRIPT"
        echo -e "\n🚀 Running step1: RNA-seq preprocessing (adapter trimming - QC) ($(date))\n"
        "$STEP1_SCRIPT" "$CSV_PATH"
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step1 finished ($(date))\n"
        else
            echo -e "\n❌ step1 failed ($(date))\n"
            exit 1
        fi
        ;;
    --step2)
        check_script_exists "$STEP2_SCRIPT"
        echo -e "\n🚀 Running step2: RNA-seq alignment (HISAT2 + StringTie) ($(date))\n"
        "$STEP2_SCRIPT"
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step2 finished ($(date))\n"
        else
            echo -e "\n❌ step2 failed ($(date))\n"
            exit 1
        fi
        ;;
    --step3)
        if [ $# -ne 2 ]; then
            echo "❌ ERROR: --step3 requires group CSV name/path!"
            echo "📚 Usage: rna_LAN --step3 <group CSV name/path>"
            echo "💡 Example1 (basename): rna_LAN --step3 group.csv"
            echo "💡 Example2 (path): rna_LAN --step3 /home/data/group.csv or ~/data/group.csv"
            exit 1
        fi
        
        GROUP_CSV_INPUT="$2"
        GROUP_CSV_PATH=$(complete_csv_path "$GROUP_CSV_INPUT")
        
        if [ -z "$GROUP_CSV_PATH" ]; then
            echo "❌ ERROR: Group CSV file not found!"
            if [[ "$GROUP_CSV_INPUT" != */* ]]; then
                echo "🔍 Searched: $(pwd)/1.raw_fastq/${GROUP_CSV_INPUT}"
            else
                echo "🔍 Searched: ${GROUP_CSV_INPUT}"
            fi
            exit 1
        fi
        
        check_r_env
        check_script_exists "$STEP3_SCRIPT"
        select_species
        
        echo "✅ Group CSV found: ${GROUP_CSV_PATH}"
        echo -e "\n🚀 Running step3: RNA-seq TPM quantification + PCA + DEG + KEGG enrichment ($(date))\n"
        
        Rscript "$STEP3_SCRIPT" "$SPECIES" "$GROUP_CSV_PATH"
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
        select_te_normalization
        select_ribo_seq_dir
        
        echo -e "\n🚀 Running step4: RNA-seq + Ribo-seq integration for TE calculation ($(date))\n"
        Rscript "$STEP4_SCRIPT" "$SPECIES" "$NORM_METHOD" "$RIBO_SEQ_DIR"
        
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step4 finished ($(date))\n"
        else
            echo -e "\n❌ step4 failed ($(date))\n"
            exit 1
        fi
        ;;
    --step5)
        if [ $# -ne 2 ]; then
            echo "❌ ERROR: --step5 requires group CSV name/path!"
            echo "📚 Usage: rna_LAN --step5 <group CSV name/path>"
            echo "💡 Example1 (basename): rna_LAN --step5 group.csv"
            echo "💡 Example2 (path): rna_LAN --step5 /home/data/group.csv or ~/data/group.csv"
            exit 1
        fi
        
        GROUP_CSV_INPUT="$2"
        GROUP_CSV_PATH=$(complete_csv_path "$GROUP_CSV_INPUT")
        
        if [ -z "$GROUP_CSV_PATH" ]; then
            echo "❌ ERROR: Group CSV file not found!"
            if [[ "$GROUP_CSV_INPUT" != */* ]]; then
                echo "🔍 Searched: $(pwd)/1.raw_fastq/${GROUP_CSV_INPUT}"
            else
                echo "🔍 Searched: ${GROUP_CSV_INPUT}"
            fi
            exit 1
        fi

        check_r_env
        check_script_exists "$STEP5_SCRIPT"
        find_te_result_file

        echo "✅ Group CSV found: ${GROUP_CSV_PATH}"
        echo -e "\n🚀 Running step5: TE differential analysis (volcano + boxplot + heatmap) ($(date))\n"
        
        Rscript "$STEP5_SCRIPT" "$TE_RESULT_FILE" "$GROUP_CSV_PATH"
        
        if [ $? -eq 0 ]; then
            echo -e "\n✅ step5 finished ($(date))\n"
            echo "📁 Differential TE output directory: $(dirname ${TE_RESULT_FILE})/TE_differential_analysis"
        else
            echo -e "\n❌ step5 failed ($(date))\n"
            exit 1
        fi
        ;;
    *)
        echo "❌ ERROR: Invalid argument!"
        echo "📚 Usage: rna_LAN --step1 <CSV name/path> | --step2 | --step3 <group CSV name/path> | --step4 | --step5 <group CSV name/path>"
        exit 1
        ;;
esac

