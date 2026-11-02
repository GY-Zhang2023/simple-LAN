# simple‑ribo
Automated Ribo‑seq & Bulk RNA‑seq end‑to‑end analysis pipeline.

Workflow
- Ribo‑seq: QC, adapter‑trim, rRNA‑tRNA depletion, STAR alignment, riboWaltz P‑site calibration, periodicity, codon occupancy
- Bulk RNA‑seq: QC → adapter trimming → alignment → gene quantification → PCA / DEG / GO‑KEGG enrichment
- Combined mode: compute Translation Efficiency(TE) and differential TE analysis

> Tested environment: Linux, conda.
> Reference genome / index files are NOT included in this repository due to large size.

## ⚡ Quick install
```bash
git clone https://github.com/YOUR_USERNAME/simple-ribo.git
cd simple-ribo
bash install.sh
conda env create -f environment.yml
```
## Activate environmentbash
```bash
source $HOME/.local/bin/../env/ribo_LAN_env.sh
ribo_LAN -h
rna_LAN -h
```
## Custom resource limit (environment variable override)
```bash
export N_CORES=8
export MEMORY_LIMIT_GB=50
# run pipeline commands ...
```
## Supported species
- human (GRCh38.p14, GENCODE v49)
- mouse (GRCm39, GENCODE vM35)
- other species in future development

## Reference preparation
Read document: [`docs/reference_prepare.md`](docs/reference_prepare.md)

> 
> You must download genome fasta, GTF annotation and build STAR / bowtie2 index manually, fill absolute paths in pipeline script `SPECIES_CONFIG` dictionary.

## Output directory structure
Read document: docs/output_structure.md
## Input template
Sample csv templates are under folder `example/`

- `rna_sample.csv` for bulk RNA‑seq
- `ribo_sample.csv` for ribo‑seq
- `group.csv` for differential analysis
## Usage exampleBulk RNA‑seq
```bash
cd your_rna_project
rna_LAN --step1 sample.csv
rna_LAN --step2
rna_LAN --step3 group.csv
rna_LAN --step4
rna_LAN --step5 group.csv
```
Ribo‑seq
```bash
cd your_ribo_project
ribo_LAN --step1 ribo_sample.csv
ribo_LAN --step2
ribo_LAN --step3
ribo_LAN --step4
ribo_LAN --step5
ribo_LAN --step6 group.csv gene_list.txt
```
## Software version list
| Software  | Version   |
|-----------|-----------|
| R         | 4.4.x     |
| cutadapt  | >=4.0     |
| STAR      | >=2.7.11  |
| HISAT2    | >=2.2.1   |
| bowtie2   | >=2.5     |
| samtools  | >=1.18    |
| riboWaltz | >=1.2.0   |

## Notes

- Do NOT run analysis inside source‑code git repository folder, use independent project workspace.
- All log / plot outputs are in project sub‑folders.
