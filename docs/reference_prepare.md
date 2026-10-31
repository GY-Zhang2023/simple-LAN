# Reference Genome & Index Preparation
This document describes how to download genome sequences, annotations and build indexes required by simple‑LAN.

> Reference files are NOT distributed in this repository due to large file size.
> Users need to prepare them manually and fill absolute paths in the `SPECIES_CONFIG` dictionary inside pipeline shell scripts.

## Supported species
- human (GRCh38.p14, GENCODE v49)
- mouse (GRCm39, GENCODE vM35)
- CHmito (Chaetomium thermophilum mitochondrial, custom mitochondrial reference)

---

## 1. Human GRCh38.p14 | GENCODE v49
### Download genome & annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz

gunzip GRCh38.p14.genome.fa.gz
gunzip gencode.v49.annotation.gtf.gz

### Build STAR genome index
STAR --runMode genomeGenerate \
--genomeDir ./star_index_human_v49 \
--genomeFastaFiles GRCh38.p14.genome.fa \
--sjdbGTFfile gencode.v49.annotation.gtf \
--sjdbOverhang 149 \
--runThreadN ${N_CORES}

### Build bowtie2 index for rRNA‑tRNA contamination filtering
1. Extract rRNA and tRNA entries from GTF file, export sequences to `rrna_trna.fasta`
2. Build bowtie2 index:
bowtie2‑build rrna_trna.fasta human_rrna_tRNA_index

---

## 2. Mouse GRCm39 | GENCODE vM35
### Download genome & annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz

gunzip *.gz

STAR index and bowtie2 rRNA‑tRNA index construction commands are identical to human workflow.

---

## 3. CHmito (Chaetomium thermophilum mitochondrial)
This is a custom mitochondrial reference.
Prepare fasta and gtf manually, then build STAR index and bowtie2 contamination index following exactly the same workflow above.

---

## 4. Fill paths in pipeline config block
After all indexes are ready, open pipeline shell script and modify `SPECIES_CONFIG` associative array:

declare -A SPECIES_CONFIG
SPECIES_CONFIG["human_ref_genome"]="/data/ref/human/GRCh38.p14.genome.fa"
SPECIES_CONFIG["human_gtf"]="/data/ref/human/gencode.v49.annotation.gtf"
SPECIES_CONFIG["human_star_index"]="/data/ref/human/star_index_human_v49"
SPECIES_CONFIG["human_rrna_trna_index"]="/data/ref/human/bowtie2_rrna_tRNA/human_rrna_tRNA_index"

> Pipeline will abort immediately with error if reference file missing.

---

## Recommended reference folder layout
your_reference_root/
├── human/
│   ├── GRCh38.p14.genome.fa
│   ├── gencode.v49.annotation.gtf
│   ├── star_index_human_v49/
│   └── bowtie2_rrna_tRNA/
├── mouse/
└── CHmito/
