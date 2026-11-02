#!/bin/bash
set -euo pipefail

echo "=================================================="
echo "          simple‑ribo One‑Click Installer          "
echo "     Bulk RNA‑seq & Ribo‑seq analysis pipeline    "
echo "=================================================="

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "[ERROR] This pipeline only supports Linux system."
    exit 1
fi

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
REPO_ROOT=$(dirname "${SCRIPT_PATH}")
TOOLS_DIR="${REPO_ROOT}/TOOLS"
LOCAL_BIN="${HOME}/.local/bin"

echo "[INFO] Repository root: ${REPO_ROOT}"
echo "[INFO] Tools directory: ${TOOLS_DIR}"

if [[ ! -d "${TOOLS_DIR}" ]]; then
    echo "[ERROR] TOOLS folder missing, check repository completeness."
    exit 1
fi

mkdir -p "${LOCAL_BIN}"

echo "[INFO] Set execute permission for scripts..."
chmod +x "${TOOLS_DIR}"/*.sh "${TOOLS_DIR}"/*.R
chmod +x "${REPO_ROOT}/env/"*.sh

echo "[INFO] Create global command symlinks under ${LOCAL_BIN}"
ln -sf "${TOOLS_DIR}/rna_LAN.sh" "${LOCAL_BIN}/rna_LAN"
ln -sf "${TOOLS_DIR}/ribo_LAN.sh" "${LOCAL_BIN}/ribo_LAN"

if ! echo "${PATH}" | grep -q "${LOCAL_BIN}"; then
    echo -e "\n[WARNING] ${LOCAL_BIN} is NOT in your PATH!"
    echo "Append this line to ~/.bashrc"
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
echo "Reference doc: ${REPO_ROOT}/docs/reference_prepare.md"
echo "=================================================="
