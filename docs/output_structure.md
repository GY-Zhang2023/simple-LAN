# Output Directory Structure
simple‑LAN separates analysis workspace from source‑code repository.
> ⚠️ Each analysis project should be an independent folder. DO NOT run analysis inside git source‑code repository directory.

## Bulk RNA‑seq project
project_root/
├── 1.raw_fastq/          # Raw fastq + sample sheet csv
├── 2.clean_fastq/        # Trimmed, filtered fastq, FastQC reports
├── 3.rna_seq_analysis/   # STAR/HISAT2 BAM, gene count / FPKM / TPM matrix
└── 4.rna_seq_qc_diff/    # PCA, sample correlation heatmap, DESeq2/edgeR DEG, GO‑KEGG enrichment tables and figures

## Ribo‑seq project
project_root/
├── 1.raw_fastq/               # Raw RPF fastq + sample sheet csv
├── 2.clean_fastq/             # Adapter trimmed, length filtered, rRNA‑tRNA depleted fastq
├── 3.alignments/              # STAR genome/transcriptome BAM, RPF‑length stats, alignment summary Report.txt
└── 4.riboWaltz_output/        # P‑site calibration, periodicity, footprint density, codon occupancy tables and plots

## Combined RNA‑seq + Ribo‑seq: Translation Efficiency (TE) analysis
TE results locate under `3.rna_seq_analysis/TE_analysis/`

TE_analysis/
├── TE_results_ALL_samples_*.tsv          # Full translation efficiency table
├── TE_RPM/TPM/VST/                       # Per‑method sub‑folder
│   └── per‑sample subfolders: TE distribution, TE‑vs‑RNA scatter PDF plots
└── TE_differential_analysis/
    ├── diff_TE_raw_FC_full_anno.tsv
    ├── TE volcano plot
    ├── TE boxplot
    └── Z‑score normalized TE heatmap

## File notes
- Temporary intermediate files are ignored by `.gitignore`.
- Log files are generated alongside each step for debugging.
- Publication‑ready figures are exported in PDF format.
- All matrix output are tab‑separated tsv for downstream R/Python analysis.
