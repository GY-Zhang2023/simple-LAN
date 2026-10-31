#!/bin/bash
# rna_LAN_env.sh: environment activator for bulk RNA‑seq only
# Usage: source rna_LAN_env.sh

CONDA_ENV_NAME="ribo_simple"

CONDA_INIT_PATH=""
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    CONDA_INIT_PATH="$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    CONDA_INIT_PATH="$HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "[ERROR] Cannot find conda init script."
    return 1
fi

source "$CONDA_INIT_PATH"

if conda info --envs | grep -q "^$CONDA_ENV_NAME\s"; then
    conda activate "$CONDA_ENV_NAME"
    echo "[OK] Conda environment '$CONDA_ENV_NAME' activated."
else
    echo "[ERROR] Conda environment '$CONDA_ENV_NAME' not found."
    echo "Please run: conda env create -f environment.yml"
    return 1
fi

if command -v rna_LAN >/dev/null 2>&1; then
    echo "[OK] Global command 'rna_LAN' available."
else
    echo "[WARNING] global command missing, re‑run bash install.sh"
fi

echo -e "\n🎉 simple‑LAN RNA‑seq environment ready."
echo "Usage example:"
echo "  rna_LAN --step1 sample.csv"
