#!/bin/bash
set -euo pipefail

echo "=================================================="
echo "          simple‑LAN One‑Click Installer          "
echo "     Bulk RNA‑seq & Ribo‑seq analysis pipeline    "
echo "=================================================="

# only support Linux
if [[ "$(uname -s)" != "Linux" ]]; then
    echo "[ERROR] This pipeline only supports Linux system."
    exit 1
fi

# get the root path
SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
REPO_ROOT=$(dirname "${SCRIPT_PATH}")
TOOLS_DIR="${REPO_ROOT}/TOOLS"
LOCAL_BIN="${HOME}/.local/bin"

echo "[INFO] Repository root: ${REPO_ROOT}"
echo "[INFO] Tools directory: ${TOOLS_DIR}"

# check TOOLS 
if [[ ! -d "${TOOLS_DIR}" ]]; then
    echo "[ERROR] TOOLS folder missing, please check repository completeness."
    exit 1
fi

# build ~/.local/bin
mkdir -p "${LOCAL_BIN}"

# set all scripts and authenity
echo "[INFO] Adding execute permission for all scripts..."
chmod +x "${TOOLS_DIR}"/*.sh "${TOOLS_DIR}"/*.R
chmod +x "${REPO_ROOT}/env/"*.sh

# soft-link
echo "[INFO] Create global command symlinks under ${LOCAL_BIN}"
ln -sf "${TOOLS_DIR}/rna_LANSummary.sh" "${LOCAL_BIN}/rna_LAN"
ln -sf "${TOOLS_DIR}/ribo_LAN.sh" "${LOCAL_BIN}/ribo_LAN"

#  ~/.local/bin in PATH
if ! echo "${PATH}" | grep -q "${LOCAL_BIN}"; then
    echo -e "\n[WARNING] ${LOCAL_BIN} is NOT in your PATH variable!"
    echo "Please append this line into your ~/.bashrc"
    echo "export PATH=\"\$HOME/.local/bin:\$PATH\""
    echo "Then run: source ~/.bashrc"
fi

echo -e "\n✅ Install finished!"
echo "-----------------------------------------"
echo "Next steps:"
echo "1. Create conda environment:"
echo "   conda env create -f ${REPO_ROOT}/environment.yml"
echo ""
echo "2. Activate environment:"
echo "   source ${REPO_ROOT}/env/ribo_LAN_env.sh"
echo ""
echo "3. Test command availability:"
echo "   rna_LAN -h"
echo "   ribo_LAN -h"
echo ""
echo "📌 Example quick start:"
echo "cd your_rna_project"
echo "rna_LAN --step1 sample.csv"
echo ""
echo "📌 Reference prepare document: ${REPO_ROOT}/docs/reference_prepare.md"
echo "=================================================="
