# Reference Genome & Index Preparation
This document describes how to download genome sequences, annotations and build indexes required by simple‑LAN.

> Reference files are NOT distributed in this repository due to large file size.
> Users need to prepare them manually and fill absolute paths in the `SPECIES_CONFIG` and `SPECIES_REF_MAP` dictionaries inside pipeline shell scripts.
> This pipeline uses **protein-coding only GTF/transcript FASTA** and canonical genome sequence instead of the full annotation set.

## Supported species
- human (GRCh38.p14 canonical, GENCODE v49 protein-coding)
- mouse (GRCm39 p6 canonical, GENCODE vM35 protein-coding)
---

## 1. Human GRCh38.p14 | GENCODE v49
### Download genome & annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz

gunzip GRCh38.p14.genome.fa.gz
gunzip gencode.v49.annotation.gtf.gz

### Filter protein-coding annotations
Filter the full `gencode.v49.annotation.gtf` to retain only protein-coding entries to generate `gencode.v49.pc_annotation.gtf`.
Extract corresponding protein-coding transcript sequences to `gencode.v49.pc_transcripts.fa`.

### Build STAR genome index
STAR --runMode genomeGenerate \
--genomeDir ./star_genome_index \
--genomeFastaFiles GRCh38.p14.canonical_genome.fa \
--sjdbGTFfile gencode.v49.pc_annotation.gtf \
--sjdbOverhang 149 \
--runThreadN ${N_CORES}

### Build bowtie2 index for rRNA‑tRNA contamination filtering
1. Extract rRNA and tRNA entries from GTF file, export sequences to `rrna_trna.fasta`
2. Build bowtie2 index:
bowtie2-build rrna_trna.fasta human_rrna_tRNA_index

---

## 2. Mouse GRCm39 p6 | GENCODE vM35
### Download genome & annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz

gunzip *.gz

### Filter protein-coding annotations
Filter the full `gencode.vM35.annotation.gtf` to retain only protein-coding entries to generate `gencode.vM35.pc_annotation.gtf`.
Extract corresponding protein-coding transcript sequences to `gencode.vM35.pc_transcripts.fa`.

STAR index and bowtie2 rRNA‑tRNA index construction commands are identical to human workflow, using the pc-filtered GTF.
STAR output index folder name: `star_mus_genome_index`

---

## 3. Fill paths in pipeline config block
After all indexes are ready, open pipeline shell script and modify `SPECIES_CONFIG` and `SPECIES_REF_MAP` associative arrays:

```bash
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

declare -A SPECIES_REF_MAP
SPECIES_REF_MAP["human"]="$HOME/ribo_references/rrna_trna/human_rrna_tRNA_index"
SPECIES_REF_MAP["mouse"]="$HOME/reference/mouse_transcriptome/mouse_rrna_tRNA_index"

> Pipeline will abort immediately with error if reference file missing.

---

## Recommended reference folder layout
$HOME/
├── reference/
│   ├── human_transcriptome/
│   │   ├── GRCh38.p14.canonical_genome.fa
│   │   ├── gencode.v49.pc_annotation.gtf
│   │   └── gencode.v49.pc_transcripts.fa
│   └── mouse_transcriptome/
│       ├── GRCm39.p6.canonical_genome.fa
│       ├── gencode.vM35.pc_annotation.gtf
│       └── gencode.vM35.pc_transcripts.fa
├── ribo_simple_LAN/
│   └── TOOLS/
│       ├── star_genome_index/
│       └── star_mus_genome_index/
└── ribo_references/
    └── rrna_trna/
        ├── mouse_rrna_tRNA_index.*
        └── human_rrna_tRNA_index.*
