# simple-LAN
# simple‑LAN
> Automated Ribo‑seq(mainly) & Bulk RNA‑seq end‑to‑end analysis pipeline for human/mouse/CHmito.
> Open source pipeline for ribosome profiling and polyA RNA‑seq, support TE(Translation Efficiency) calculation and differential analysis.

## ⚡ Quick install (One‑click)
```bash
# clone repository
git clone https://github.com/YOUR_USERNAME/simple-LAN.git
cd simple-LAN

# run one‑click install
bash install.sh

# create conda environment
conda env create -f environment.yml
```

### Activate environment
```bash
source $HOME/.local/bin/../env/ribo_LAN_env.sh
# test commands
rna_LAN -h
ribo_LAN -h
```

## Workflow
- Bulk RNA‑seq: QC → trim → HISAT2/STAR → quantification → PCA/DESeq2/edgeR → GO/KEGG
- Ribo‑seq: Trim‑galore → bowtie2 rRNA/tRNA depletion → STAR alignment → riboWaltz P‑site calibration
- Combined: RNA+Ribo compute TE and differential TE analysis

## Input example
See `example/` folder for csv template.

## Reference genome preparation
Please read `docs/reference_prepare.md`.
> Reference fasta/gtf/index are NOT included in this repository due to large file size, user need download and build index manually.

## Output structure
See `docs/output_structure.md`

## Supported species
- human GRCh38.p14 gencode v49
- mouse GRCm39.p6 gencode vM35
- CHmito Chaetomium thermophilum mitochondrial

> Test environment: Linux, conda, R‑4.4
