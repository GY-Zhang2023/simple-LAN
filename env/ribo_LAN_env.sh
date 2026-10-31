#!/bin/bash
# ribo_LAN_env.sh: environment activator for ribo‑seq & RNA‑seq pipeline
# Usage: source ribo_LAN_env.sh

CONDA_ENV_NAME="ribo_simple"

# Auto‑detect conda init script
CONDA_INIT_PATH=""
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    CONDA_INIT_PATH="$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    CONDA_INIT_PATH="$HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "[ERROR] Cannot find conda init script, check miniconda/anaconda install path."
    return 1
fi

source "$CONDA_INIT_PATH"

if conda info --envs | grep -q "^$CONDA_ENV_NAME\s"; then
    conda activate "$CONDA_ENV_NAME"
    echo "[OK] Conda environment '$CONDA_ENV_NAME' activated."
else
    echo "[ERROR] Conda environment '$CONDA_ENV_NAME' not found."
    echo "Please run: conda env create -f environment.yml at repository root."
    return 1
fi

# Check global command symlink
if command -v ribo_LAN >/dev/null 2>&1 && command -v rna_LAN >/dev/null 2>&1; then
    echo "[OK] Global commands 'rna_LAN' 'ribo_LAN' available."
else
    echo "[WARNING] global commands missing, please re‑run bash install.sh under repository root."
fi

echo -e "\n🎉 simple‑LAN environment ready."
echo "Usage examples:"
echo "  rna_LAN --step1 sample.csv"
echo "  ribo_LAN --step1 ribo_sample.csv"
