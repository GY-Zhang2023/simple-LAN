#!/usr/bin/env Rscript
#SBATCH --job-name=ribo_diff
#SBATCH --output=ribo_diff_%j.out
#SBATCH --error=ribo_diff_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=2
#SBATCH --mem=50G

N_CORES <- 8
N_THREADS <- 2
MEMORY_LIMIT <- "50G"

Sys.setenv(LANGUAGE = "en")
Sys.setenv(R_MAX_NUM_DLLS = 999)
options(stringsAsFactors = FALSE, scipen = 999)
rm(list = ls())

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
selected_species    <- args[1]
min_len_default     <- as.integer(args[2])
max_len_default     <- as.integer(args[3])
min_len             <- as.integer(args[4])
max_len             <- as.integer(args[5])
selected_analysis   <- args[6]
trim_nt             <- as.integer(args[7])
csv_path            <- args[8]
gene_list_path      <- if (length(args) >= 9) args[9] else "list.txt"

codon2aa <- data.table(
  codon = c("TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
            "ATT", "ATC", "ATA", "ATG",
            "GTT", "GTC", "GTA", "GTG",
            "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
            "CCT", "CCC", "CCA", "CCG",
            "ACT", "ACC", "ACA", "ACG",
            "GCT", "GCC", "GCA", "GCG",
            "TAT", "TAC", "TAA", "TAG",
            "CAT", "CAC", "CAA", "CAG",
            "AAT", "AAC", "AAA", "AAG",
            "GAT", "GAC", "GAA", "GAG",
            "TGT", "TGC", "TGA", "TGG",
            "CGT", "CGC", "CGA", "CGG",
            "AGA", "AGG",
            "GGT", "GGC", "GGA", "GGG"),
  aa    = c("F", "F", "L", "L", "L", "L", "L", "L",
            "I", "I", "I", "M",
            "V", "V", "V", "V",
            "S", "S", "S", "S", "S", "S",
            "P", "P", "P", "P",
            "T", "T", "T", "T",
            "A", "A", "A", "A",
            "Y", "Y", "*", "*",
            "H", "H", "Q", "Q",
            "N", "N", "K", "K",
            "D", "D", "E", "E",
            "C", "C", "*", "W",
            "R", "R", "R", "R",
            "R", "R",
            "G", "G", "G", "G")
)

PATH_HUMAN_FASTA <- "~/reference/human_transcriptome/gencode.v49.pc_transcripts.fa"
PATH_MOUSE_FASTA  <- "~/reference/mouse_transcriptome/gencode.vM35.pc_transcripts.fa"

species_config <- list(
  human = list(
    gtf   = "~/reference/human_transcriptome/gencode.v49.pc_annotation.gtf",
    fasta = PATH_HUMAN_FASTA
  ),
  mouse = list(
    gtf   = "~/reference/mouse_transcriptome/gencode.vM35.pc_annotation.gtf",
    fasta = PATH_MOUSE_FASTA
  )
)

gtf_path           <- species_config[[selected_species]]$gtf
transcriptome_fasta <- species_config[[selected_species]]$fasta

cat("Loading pre-generated transcript-to-gene mapping...\n")

gtf_path <- species_config[[selected_species]]$gtf
gtf_dir  <- dirname(gtf_path)

tx2gene_filename <- ifelse(selected_species == "human",
                            "tx2gene_human.txt",
                            "tx2gene_mouse.txt")
tx2gene_path <- file.path(gtf_dir, tx2gene_filename)

if (!file.exists(tx2gene_path)) {
  stop(paste0(
    "Mapping file not found!\n",
    "   Species: ", selected_species, "\n",
    "   Expected path: ", tx2gene_path, "\n",
    "   Please confirm tx2gene file is in the same folder as GTF"
  ))
}

tx2gene <- data.table::fread(
  tx2gene_path,
  header      = FALSE,
  stringsAsFactors = FALSE,
  col.names   = c("transcript", "gene_name")
)
tx2gene <- unique(tx2gene, by = "transcript")

cat(paste0("Successfully loaded ", nrow(tx2gene), " transcript-gene mappings\n"))
cat(paste0("   Example: ", tx2gene[1, transcript], " -> ", tx2gene[1, gene_name], "\n\n"))

if (selected_analysis == "disome") {
  flow_text <- paste0("Disome pipeline - 5' trim ", trim_nt, "nt + 5' Rear Codon inference + 3nt periodicity correction")
} else {
  flow_text <- "Standard Monosome pipeline (Load TXT + RDS for footprint)"
}

cat(paste0("\nLoaded reference config for species: ", selected_species, "\n"))
cat(paste0("   Pipeline: ", flow_text, "\n"))
cat(paste0("   Codon base: U -> T (DNA standard)\n"))
cat(paste0("   PCA logic: sample=row | codon=feature | avg_occupancy=expression\n"))
cat(paste0("   Differential test: KS test per codon\n"))
cat(paste0("   Visualization: sample curve plots, no heatmap\n\n"))

WORK_DIR      <- getwd()
output_root   <- paste0(WORK_DIR, "/5.ribo_group_diff_output")
rds_root_dir  <- paste0(WORK_DIR, "/4.riboWaltz_output")
report_path   <- file.path(output_root, "batch_run_report.txt")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_root)) {
  dir.create(output_root, recursive = TRUE)
  cat("Output directory created: ", output_root, "\n")
}

sink(report_path, split = TRUE)

cat("============================================================\n")
cat(paste0("Codon Occupancy Analysis Report - ", flow_text))
cat("\nAnalysis time: ", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

library(dplyr)
library(data.table)
library(ggplot2)
library(Biostrings)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(tidyr)
library(FactoMineR)
library(factoextra)

threads <- N_CORES
rpf_count_threshold <- 1e6
bin_size <- 0.5
length_range <- min_len:max_len
p_threshold <- 0.1

get_group_similar_color <- function(sample_names, group_dt) {
  group_list <- unique(group_dt$group_name)
  color_palette <- c("#E63946", "#457B9D", "#2A9D8F", "#F4A261",
                     "#264653", "#8D99AE", "#EF476F", "#06D6A0")
  names(color_palette) <- group_list
  sample_color <- c()
  for (s in sample_names) {
    g <- group_dt[group_dt$sample_name == s, ]$group_name
    idx <- which(sample_names[sapply(sample_names,
                                     function(x) group_dt[group_dt$sample_name == x, ]$group_name) == g] == s)
    shade <- 0.2 * (idx - 1)
    rgb_val <- col2rgb(color_palette[g]) / 255
    new_rgb <- pmin(rgb_val + shade, 1)
    sample_color[s] <- rgb(new_rgb[1], new_rgb[2], new_rgb[3])
  }
  return(sample_color)
}

cat("1. Analysis Configuration\n")
cat("------------------------------------------------------------\n")
cat(paste0("- Cores: ", N_CORES, "\n"))
cat(paste0("- Threads per core: ", N_THREADS, "\n"))
cat(paste0("- Memory limit: ", MEMORY_LIMIT, "\n"))
cat("- System language: English\n")
cat("- Codon base correction: U -> T (DNA standard)\n")
cat("- Differential threshold: P <", p_threshold, "\n\n")

cat("2. Load Sample Group Info & Match File Paths\n")
cat("------------------------------------------------------------\n")

group_dt <- fread(csv_path, sep = ",")
group_dt <- group_dt[!is.na(group_dt$sample_name) & !is.na(group_dt$group_name), ]
sample_names <- unique(group_dt$sample_name)

txt_paths <- sapply(sample_names, function(x) {
  file.path(rds_root_dir, x, "all_site_codon_occupancy_final.txt")
})
rds_paths <- sapply(sample_names, function(x) {
  file.path(rds_root_dir, x, "reads_psite_raw_version.rds")
})
names(txt_paths) <- sample_names
names(rds_paths) <- sample_names

missing_txt <- names(txt_paths)[!file.exists(txt_paths)]
missing_rds <- names(rds_paths)[!file.exists(rds_paths)]
missing_samples <- unique(c(missing_txt, missing_rds))

if (length(missing_samples) > 0) {
  cat(paste0("Missing files for samples: ", paste(missing_samples, collapse = ", "), "\n"))
  sample_names <- sample_names[!sample_names %in% missing_samples]
  txt_paths    <- txt_paths[!names(txt_paths) %in% missing_samples]
  rds_paths    <- rds_paths[!names(rds_paths) %in% missing_samples]
}

if (length(sample_names) == 0) {
  sink()
  stop("No valid samples with complete files! Exit analysis.")
}

cat(paste0("Matched ", length(sample_names), " valid samples\n"))
group_summary <- group_dt[, list(N = .N), by = group_name]
group_detail <- paste0(group_summary$group_name, "(", group_summary$N, ")") %>% paste(collapse = " | ")
cat(paste0("Groups: ", group_detail, "\n\n"))

group_dt <- group_dt[group_dt$sample_name %in% sample_names, ]
sample_color_vec <- get_group_similar_color(sample_names, group_dt)

cat("3. Load Codon Occupancy Data & Footprint Data\n")
cat("------------------------------------------------------------\n")

all_codon_occupancy_raw <- data.table()
all_psite_raw_dt        <- data.table()

for (sample in sample_names) {
  cat(paste0("   Processing sample: ", sample, "\n"))

  codon_occupancy_raw <- fread(txt_paths[sample], sep = "\t", stringsAsFactors = FALSE)
  codon_occupancy_raw$codon_type <- paste0(codon_occupancy_raw$site, "_codon")
  codon_occupancy_raw$sample_name <- sample
  codon_occupancy_raw$group_name <- group_dt[group_dt$sample_name == sample, ]$group_name
  all_codon_occupancy_raw <- rbindlist(list(all_codon_occupancy_raw, codon_occupancy_raw), fill = TRUE)

  reads_psite_dt <- readRDS(rds_paths[sample])
  setDT(reads_psite_dt)
  reads_psite_dt$sample_name <- sample
  reads_psite_dt$group_name  <- group_dt[group_dt$sample_name == sample, ]$group_name
  reads_psite_dt <- merge(reads_psite_dt, tx2gene, by = "transcript", all.x = TRUE)
  all_psite_raw_dt <- rbindlist(list(all_psite_raw_dt, reads_psite_dt), fill = TRUE)

  cat(paste0("     Done: codon rows=", nrow(codon_occupancy_raw),
             " | footprint rows=", nrow(reads_psite_dt), "\n"))
}

all_codon_occupancy_group <- all_codon_occupancy_raw[, list(
  mean_occupancy   = mean(as.numeric(avg_occupancy)),
  Normalized_mean  = mean(as.numeric(Normalized_avg)),
  sd_occupancy     = sd(as.numeric(avg_occupancy)),
  se_occupancy     = sd(as.numeric(avg_occupancy)) / sqrt(.N),
  total_reads      = sum(as.numeric(total_reads)),
  transcript_count = mean(as.numeric(transcript_count))
), by = list(group_name, codon_type, codon, site)]

all_codon_occupancy_group$norm_occupancy <-
  all_codon_occupancy_group$mean_occupancy /
  ave(all_codon_occupancy_group$mean_occupancy,
      all_codon_occupancy_group$codon_type,
      FUN = mean)

fwrite(all_codon_occupancy_raw,
       file.path(output_root, "all_sample_APEcodon_occupancy_raw.txt"), sep = "\t")
fwrite(all_codon_occupancy_group,
       file.path(output_root, "all_sample_APEcodon_occupancy_group.txt"), sep = "\t")

cat("All data loaded successfully\n\n")

plot_Psite_PCA <- function(all_codon_occupancy_raw, group_dt, output_dir) {
  cat("4. Generate PCA Plot (P-site Codon Occupancy)\n")

  p_raw <- all_codon_occupancy_raw[codon_type == "P_codon",
                                   c("sample_name", "codon", "avg_occupancy")]
  pca_matrix <- dcast(p_raw, sample_name ~ codon,
                      value.var = "avg_occupancy", fill = 0)
  sample_grp <- group_dt[match(pca_matrix$sample_name,
                               group_dt$sample_name), ]$group_name

  pca_res <- PCA(pca_matrix[, -1], scale.unit = TRUE, ncp = 5, graph = FALSE)

  p <- fviz_pca_ind(pca_res,
                    geom.ind = c("point", "text"),
                    pointshape = 21, pointsize = 4,
                    fill.ind = sample_grp,
                    col.ind = "black",
                    palette = "Set1",
                    label = "all", repel = TRUE,
                    text.size = 3,
                    title = "PCA Analysis - P-site Codon Occupancy",
                    xlab = paste0("PC1 (", sprintf("%.1f", pca_res$eig[1, 2]), "%)"),
                    ylab = paste0("PC2 (", sprintf("%.1f", pca_res$eig[2, 2]), "%)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 9),
          legend.position = "top",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))

  ggsave(file.path(output_dir, "0.Psite_all_sample_PCA.pdf"),
         p, width = 10, height = 8, dpi = 300)
  cat("PCA plot generated\n")

  mat <- as.matrix(pca_matrix[, -1])
  rownames(mat) <- pca_matrix$sample_name
  sample_cor <- cor(t(mat), method = "pearson")

  annotation_df <- data.frame(Group = sample_grp)
  rownames(annotation_df) <- pca_matrix$sample_name

  pdf(file.path(output_dir, "0.Psite_sample_correlation_heatmap.pdf"),
      width = 9, height = 7)
  pheatmap(sample_cor,
           scale = "none",
           color = colorRampPalette(c("#457B9D", "white", "#E63946"))(100),
           main = "Sample Correlation (P-site Codon Occupancy, Pearson)",
           display_numbers = TRUE,
           number_color = "black",
           fontsize_number = 7,
           treeheight_row = 15,
           treeheight_col = 15,
           annotation_col = annotation_df,
           annotation_row = annotation_df,
           annotation_names_col = FALSE,
           annotation_names_row = FALSE)
  dev.off()

  cat("Sample correlation heatmap generated\n\n")
}

plot_group_footprint <- function(all_psite_raw_dt, output_dir) {
  cat("5. Generate Ribosome Footprint Density Plots\n")
  cat("------------------------------------------------------------\n")

  psite_cds_dt <- all_psite_raw_dt[!is.na(cds_start) &
                                     !is.na(cds_stop) &
                                     !is.na(psite)]
  psite_cds_dt$l_cds        <- psite_cds_dt$cds_stop - psite_cds_dt$cds_start + 1
  psite_cds_dt$cds_relative_pos <-
    (psite_cds_dt$psite - psite_cds_dt$cds_start + 1) /
    psite_cds_dt$l_cds * 100

  cds_density_dt <- psite_cds_dt[!is.na(cds_relative_pos)]
  cds_density_dt$cds_bin <- round(cds_density_dt$cds_relative_pos / bin_size) * bin_size
  cds_density_dt <- cds_density_dt[cds_bin >= -25 & cds_bin <= 125]
  cds_density_dt <- cds_density_dt[, list(psite_count = .N),
                                   by = list(group_name, cds_bin)]
  cds_density_dt <- cds_density_dt[order(group_name, cds_bin)]
  cds_density_dt$global_avg <- ave(cds_density_dt$psite_count,
                                   cds_density_dt$group_name, FUN = mean)
  cds_density_dt$relative_density <-
    cds_density_dt$psite_count / cds_density_dt$global_avg

  final_plot <- ggplot(cds_density_dt,
                       aes(x = cds_bin, y = relative_density,
                           color = group_name, group = group_name)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "#333f50",
               linetype = "dashed", linewidth = 1, alpha = 0.7) +
    geom_vline(xintercept = 100, color = "#333f50",
               linetype = "dashed", linewidth = 1, alpha = 0.7) +
    annotate("text", x = 2,
             y = max(cds_density_dt$relative_density) * 0.9,
             label = "Start Codon (ATG)",
             color = "#333f50", size = 3.5, hjust = 0) +
    annotate("text", x = 98,
             y = max(cds_density_dt$relative_density) * 0.9,
             label = "Stop Codon",
             color = "#333f50", size = 3.5, hjust = 1) +
    labs(x = "Normalized CDS Length (%)",
         y = "Relative Footprint Density (a.u.)",
         title = "Ribosome Footprint Density - Group Mean",
         subtitle = paste0("Binsize = ", bin_size,
                           "% | Length filter: ",
                           min_len, "-", max_len, "nt")) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

  ggsave(file.path(output_dir, "1.group_footprint_global_CDS_percent.pdf"),
         final_plot, width = 12, height = 6, dpi = 300)

  psite_cds_dt$start_nt_pos <- psite_cds_dt$psite - psite_cds_dt$cds_start
  start_codon_dt <- psite_cds_dt[!is.na(start_nt_pos) &
                                   start_nt_pos >= -25 &
                                   start_nt_pos <= 125]
  start_codon_dt <- start_codon_dt[, list(psite_count = .N),
                                   by = list(group_name, start_nt_pos)]
  start_codon_dt <- start_codon_dt[order(group_name, start_nt_pos)]
  start_codon_dt$global_avg <- ave(start_codon_dt$psite_count,
                                   start_codon_dt$group_name, FUN = mean)
  start_codon_dt$relative_density <-
    start_codon_dt$psite_count / start_codon_dt$global_avg

  start_plot <- ggplot(start_codon_dt,
                       aes(x = start_nt_pos, y = relative_density,
                           color = group_name, group = group_name)) +
    geom_line(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "#e84a5f",
               linetype = "solid", linewidth = 1.2, alpha = 0.9) +
    annotate("text", x = 5,
             y = max(start_codon_dt$relative_density) * 0.9,
             label = "ATG (Start Codon)",
             color = "#e84a5f", size = 4, hjust = 0, fontface = "bold") +
    labs(x = "Position relative to Start Codon (nt)",
         y = "Relative Footprint Density (a.u.)",
         title = "Footprint Density - Start Codon Region") +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

  ggsave(file.path(output_dir, "2.group_footprint_StartCodon.pdf"),
         start_plot, width = 10, height = 5, dpi = 300)

  psite_cds_dt$stop_nt_pos <- psite_cds_dt$psite - psite_cds_dt$cds_stop
  stop_codon_dt <- psite_cds_dt[!is.na(stop_nt_pos) &
                                  stop_nt_pos >= -125 &
                                  stop_nt_pos <= 25]
  stop_codon_dt <- stop_codon_dt[, list(psite_count = .N),
                                 by = list(group_name, stop_nt_pos)]
  stop_codon_dt <- stop_codon_dt[order(group_name, stop_nt_pos)]
  stop_codon_dt$global_avg <- ave(stop_codon_dt$psite_count,
                                   stop_codon_dt$group_name, FUN = mean)
  stop_codon_dt$relative_density <-
    stop_codon_dt$psite_count / stop_codon_dt$global_avg

  stop_plot <- ggplot(stop_codon_dt,
                      aes(x = stop_nt_pos, y = relative_density,
                          color = group_name, group = group_name)) +
    geom_line(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "#e84a5f",
               linetype = "solid", linewidth = 1.2, alpha = 0.9) +
    annotate("text", x = -5,
             y = max(stop_codon_dt$relative_density) * 0.9,
             label = "Stop Codon",
             color = "#e84a5f", size = 4, hjust = 1, fontface = "bold") +
    labs(x = "Position relative to Stop Codon (nt)",
         y = "Relative Footprint Density (a.u.)",
         title = "Footprint Density - Stop Codon Region") +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

  ggsave(file.path(output_dir, "3.group_footprint_StopCodon.pdf"),
         stop_plot, width = 10, height = 5, dpi = 300)

  cat("3 footprint density plots generated\n\n")
}

plot_gene_specific_footprint <- function(all_psite_raw_dt, gene_list_path, tx2gene, output_dir) {
  cat("5.1 Generate Gene-specific Ribosome Footprint Plots\n")
  cat("------------------------------------------------------------\n")

  required_cols <- c("psite", "cds_start", "cds_stop",
                     "transcript", "group_name", "sample_name")
  missing_cols <- setdiff(required_cols, colnames(all_psite_raw_dt))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing required columns: ",
                paste(missing_cols, collapse = ", ")))
  }

  global_group_total <- all_psite_raw_dt[, .(global_total = .N), by = group_name]

  target_genes <- fread(gene_list_path, header = FALSE)$V1
  target_genes <- unique(target_genes[!is.na(target_genes) & target_genes != ""])
  if (length(target_genes) == 0) return()

  tx2gene_filtered <- tx2gene[gene_name %in% target_genes]
  gene2tx <- split(tx2gene_filtered$transcript, tx2gene_filtered$gene_name)

  psite_dt <- all_psite_raw_dt[, ..required_cols]
  psite_dt[, atg_pos      := cds_start]
  psite_dt[, cds_length   := cds_stop - cds_start + 1]
  psite_dt[, frame        := (psite - cds_start) %% 3]
  psite_dt[, in_cds       := psite >= cds_start & psite <= cds_stop]

  gene_plot_root <- file.path(output_dir, "gene_specific_footprint")
  dir.create(gene_plot_root, recursive = TRUE, showWarnings = FALSE)

  total_tx_plots  <- 0
  total_sum_plots <- 0
  total_mean_plots <- 0

  for (gene in target_genes) {
    if (!gene %in% names(gene2tx)) next

    target_transcripts <- gene2tx[[gene]]
    gene_dir <- file.path(gene_plot_root, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    gene_all_dt <- psite_dt[transcript %in% target_transcripts]
    if (nrow(gene_all_dt) == 0) next

    tx_info_all <- unique(gene_all_dt[, .(transcript, atg_pos, cds_length)])
    setorder(tx_info_all, -cds_length)

    longest_tx       <- tx_info_all$transcript[1]
    longest_tx_atg   <- tx_info_all$atg_pos[1]
    longest_tx_cds_len <- tx_info_all$cds_length[1]
    longest_tx_cds_end <- longest_tx_atg + longest_tx_cds_len - 1

    for (tx in target_transcripts) {
      tx_dt <- psite_dt[transcript == tx]
      if (nrow(tx_dt) == 0) next

      tx_atg     <- unique(tx_dt$atg_pos)
      tx_cds_len <- unique(tx_dt$cds_length)
      tx_cds_end <- tx_atg + tx_cds_len - 1

      tx_density <- tx_dt[, .(raw = .N), by = .(group_name, psite, frame)]
      tx_density <- merge(tx_density, global_group_total, by = "group_name")
      tx_density[, rpm := (raw / global_total) * 1e6]

      fwrite(tx_density[, .(group_name, psite, frame, raw, global_total, rpm)],
             file.path(gene_dir, paste0("gene_", gene, "_tx_", tx, "_rpm_table.txt")),
             sep = "\t")

      tx_sum <- tx_density[, .(gene_reads = sum(raw)), by = group_name]
      tx_sum <- merge(tx_sum, global_group_total, by = "group_name")
      tx_sum[, rpm_total := (gene_reads / global_total) * 1e6]
      tx_sum[, label := paste0(group_name, ": ",
                               format(gene_reads, big.mark = ","),
                               " (", format(global_total, big.mark = ","),
                               ") | RPM: ", round(rpm_total, 2))]
      tx_label <- paste(tx_sum$label, collapse = "\n")

      p_tx <- ggplot(tx_density, aes(x = psite, y = rpm, fill = group_name)) +
        geom_col(width = 0.8, alpha = 0.8, position = "identity") +
        geom_vline(xintercept = tx_atg, color = "#e84a5f", linewidth = 1.2) +
        {if (tx == longest_tx)
          annotate("text",
                   x = max(tx_density$psite),
                   y = max(tx_density$rpm) * 1.1,
                   label = "* Longest TX",
                   color = "#2E8B57", size = 4, hjust = 1)} +
        geom_vline(xintercept = tx_cds_end,
                   color = "#2A9D8F", linewidth = 1, linetype = "dashed") +
        annotate("text", x = Inf, y = Inf, label = tx_label,
                 hjust = 1.1, vjust = 1.1, size = 3.3, fontface = "bold") +
        labs(x = "PSite",
             y = "RPM (per group total reads)",
             title = paste0(gene, " | ", tx),
             fill = "Group") +
        scale_fill_brewer(palette = "Set1") +
        theme_bw()

      ggsave(file.path(gene_dir, paste0("gene_", gene, "_tx_", tx, "_rpm.pdf")),
             p_tx, width = 14, height = 6)
      total_tx_plots <- total_tx_plots + 1

      tx_cds_frame_dt <- tx_dt[in_cds == TRUE]
      if (nrow(tx_cds_frame_dt) > 0) {
        tx_cds_frame_dt[, codon_idx := floor((psite - cds_start) / 3)]

        group_frame_raw <- tx_cds_frame_dt[, .(cnt = .N),
                                           by = .(group_name, codon_idx, frame)]
        group_frame_cast <- dcast(group_frame_raw,
                                  group_name + codon_idx ~ frame,
                                  value.var = "cnt", fill = 0)

        frame_cols <- intersect(as.character(c(0, 1, 2)),
                                colnames(group_frame_cast))
        f12_cols   <- intersect(c("1", "2"), frame_cols)

        group_frame_cast[, total := rowSums(.SD, na.rm = TRUE),
                         .SDcols = frame_cols]
        group_frame_cast[, sum12 := rowSums(.SD, na.rm = TRUE),
                         .SDcols = f12_cols]
        group_frame_cast[, non0_frac := ifelse(total == 0, NA,
                                               sum12 / total * 100)]

        fwrite(group_frame_cast,
               file.path(gene_dir,
                         paste0("gene_", gene, "_tx_", tx,
                                "_non0_frame_pct_table.txt")),
               sep = "\t")

        p_tx_frame <- ggplot(group_frame_cast,
                             aes(x = codon_idx, y = non0_frac,
                                 fill = group_name)) +
          geom_col(position = position_dodge(width = 2), width = 1.8) +
          geom_hline(yintercept = 66.67,
                     linetype = "dotted", color = "gray50") +
          geom_vline(xintercept = 0,
                     color = "#e84a5f", linetype = "solid", linewidth = 1.2) +
          labs(x = "Codon index (starting from ATG)",
               y = "Out-of-frame fraction (%)",
               title = paste0(gene, " | ", tx,
                              " | Out-of-frame Percentage per codon"),
               fill = "Group") +
          theme_bw()

        ggsave(file.path(gene_dir,
                         paste0("gene_", gene, "_tx_", tx,
                                "_non0_frame_pct.pdf")),
               p_tx_frame, width = 16, height = 6)
      }
    }

    gene_sum <- gene_all_dt[, .(raw_sum = .N),
                            by = .(group_name, psite, frame)]
    gene_sum <- merge(gene_sum, global_group_total, by = "group_name")
    gene_sum[, rpm_sum := (raw_sum / global_total) * 1e6]

    gene_sum_agg <- gene_sum[, .(raw_sum = sum(raw_sum),
                                 rpm_sum = sum(rpm_sum)),
                            by = .(group_name, psite)]

    gene_sum_total <- gene_sum[, .(gene_reads = sum(raw_sum)), by = group_name]
    gene_sum_total <- merge(gene_sum_total, global_group_total, by = "group_name")
    gene_sum_total[, rpm_total := (gene_reads / global_total) * 1e6]
    gene_sum_total[, label := paste0(group_name, ": ",
                                     format(gene_reads, big.mark = ","),
                                     " (", format(global_total, big.mark = ","),
                                     ") | RPM: ", round(rpm_total, 2))]
    sum_label <- paste(gene_sum_total$label, collapse = "\n")

    p_sum <- ggplot(gene_sum_agg, aes(x = psite, y = rpm_sum, fill = group_name)) +
      geom_col(width = 0.8, alpha = 0.8, position = "identity") +
      geom_vline(xintercept = longest_tx_atg,
                 color = "#e84a5f", linewidth = 1.5) +
      annotate("text", x = longest_tx_atg,
               y = max(gene_sum_agg$rpm_sum) * 1.1,
               label = paste0("ATG: ", longest_tx_atg),
               color = "#e84a5f", size = 4.5, fontface = "bold") +
      geom_vline(xintercept = longest_tx_cds_end,
                 color = "#2A9D8F", linewidth = 1, linetype = "dashed") +
      annotate("text", x = Inf, y = Inf, label = sum_label,
               hjust = 1.1, vjust = 1.1, size = 3.8, fontface = "bold") +
      labs(x = "PSite",
           y = "RPM Sum (per group total reads)",
           title = paste0(gene, " | Sum RPM"),
           fill = "Group") +
      scale_fill_brewer(palette = "Set1") +
      theme_bw()

    ggsave(file.path(gene_dir, paste0("gene_", gene, "_sum_rpm.pdf")),
           p_sum, width = 16, height = 7)
    total_sum_plots <- total_sum_plots + 1

    gene_mean_dt <- gene_all_dt[, .(raw_sum = .N),
                                by = .(group_name, psite, frame)]
    n_tx <- length(unique(gene_all_dt$transcript))
    gene_mean_dt[, raw_mean := raw_sum / n_tx]
    gene_mean_dt <- merge(gene_mean_dt, global_group_total, by = "group_name")
    gene_mean_dt[, rpm_mean := (raw_mean / global_total) * 1e6]

    gene_mean_agg <- gene_mean_dt[, .(raw_mean = sum(raw_mean),
                                      rpm_mean = sum(rpm_mean)),
                                 by = .(group_name, psite)]

    gene_mean_total <- gene_mean_dt[, .(gene_reads = sum(raw_sum)), by = group_name]
    gene_mean_total <- merge(gene_mean_total, global_group_total, by = "group_name")
    gene_mean_total[, rpm_total := (gene_reads / global_total) * 1e6]
    gene_mean_total[, label := paste0(group_name, ": ",
                                     format(gene_reads, big.mark = ","),
                                     " (", format(global_total, big.mark = ","),
                                     ") | RPM: ", round(rpm_total, 2))]
    mean_label <- paste(gene_mean_total$label, collapse = "\n")

    p_mean <- ggplot(gene_mean_agg,
                     aes(x = psite, y = rpm_mean, fill = group_name)) +
      geom_col(width = 0.8, alpha = 0.8, position = "identity") +
      geom_vline(xintercept = longest_tx_atg,
                 color = "#e84a5f", linewidth = 1.5) +
      annotate("text", x = longest_tx_atg,
               y = max(gene_mean_agg$rpm_mean) * 1.1,
               label = paste0("ATG: ", longest_tx_atg),
               color = "#e84a5f", size = 4.5, fontface = "bold") +
      geom_vline(xintercept = longest_tx_cds_end,
                 color = "#2A9D8F", linewidth = 1, linetype = "dashed") +
      annotate("text", x = Inf, y = Inf, label = mean_label,
               hjust = 1.1, vjust = 1.1, size = 3.8, fontface = "bold") +
      labs(x = "PSite",
           y = "RPM Mean (per group total reads)",
           title = paste0(gene, " | Mean RPM"),
           fill = "Group") +
      scale_fill_brewer(palette = "Set1") +
      theme_bw()

    ggsave(file.path(gene_dir, paste0("gene_", gene, "_mean_rpm.pdf")),
           p_mean, width = 16, height = 7)
    total_mean_plots <- total_mean_plots + 1

    fwrite(gene_sum[, .(group_name, psite, frame, raw_sum, rpm_sum)],
           file.path(gene_dir, paste0("gene_", gene, "_sum_rpm_table.txt")),
           sep = "\t")
    fwrite(gene_mean_dt[, .(group_name, psite, frame, raw_sum, raw_mean, rpm_mean)],
           file.path(gene_dir, paste0("gene_", gene, "_mean_rpm_table.txt")),
           sep = "\t")
  }

  cat("Gene-specific footprint plots finished\n")
  cat(paste0("Single transcript RPM plots: ", total_tx_plots, "\n"))
  cat(paste0("Gene SUM RPM plots: ", total_sum_plots, "\n"))
  cat(paste0("Gene MEAN RPM plots: ", total_mean_plots, "\n\n"))
}

plot_APE_codon_ratio_hist <- function(all_codon_occupancy_group, diff_result, output_dir) {
  cat("7. Generate A/P/E Codon Occupancy Ratio Histograms\n")
  cat("------------------------------------------------------------\n")

  group_names <- unique(all_codon_occupancy_group$group_name)
  codon_types <- unique(all_codon_occupancy_group$codon_type)

  if (length(group_names) != 2) {
    cat("Only two-group ratio analysis is supported. Skip this step.\n\n")
    return()
  }

  for (ctype in codon_types) {
    cat(paste0("   Plotting: ", ctype, "\n"))

    group_sub <- all_codon_occupancy_group[codon_type == ctype]
    ratio_dt <- dcast(group_sub, codon ~ group_name,
                      value.var = c("mean_occupancy", "se_occupancy"))
    colnames(ratio_dt) <- gsub("mean_occupancy_", "", colnames(ratio_dt))
    colnames(ratio_dt) <- gsub("se_occupancy_", "se_", colnames(ratio_dt))

    for (col in colnames(ratio_dt)[-1]) {
      ratio_dt[[col]] <- as.numeric(ratio_dt[[col]])
    }

    mean_col1 <- which(colnames(ratio_dt) == group_names[1])
    mean_col2 <- which(colnames(ratio_dt) == group_names[2])
    se_col1   <- which(colnames(ratio_dt) == paste0("se_", group_names[1]))
    se_col2   <- which(colnames(ratio_dt) == paste0("se_", group_names[2]))

    ratio_dt$codon_ratio <- ratio_dt[[mean_col2]] / ratio_dt[[mean_col1]]
    ratio_dt$se_ratio <- ratio_dt$codon_ratio *
      sqrt((ratio_dt[[se_col2]] / ratio_dt[[mean_col2]])^2 +
           (ratio_dt[[se_col1]] / ratio_dt[[mean_col1]])^2)
    ratio_dt$expr_status <- ifelse(ratio_dt$codon_ratio > 1, "up",
                                   ifelse(ratio_dt$codon_ratio < 1, "down", "non_diff"))

    ratio_dt <- merge(ratio_dt, codon2aa, by = "codon", all.x = TRUE)
    ratio_dt$aa <- ifelse(is.na(ratio_dt$aa), "*", ratio_dt$aa)
    ratio_dt$codon_aa <- paste0(ratio_dt$codon, " (", ratio_dt$aa, ")")

    if (!is.null(diff_result) && nrow(diff_result) > 0) {
      diff_sub <- diff_result[codon_type == ctype]
      ratio_dt <- merge(ratio_dt, diff_sub[, .(codon, diff_type)],
                        by = "codon", all.x = TRUE)
      ratio_dt$diff_type[is.na(ratio_dt$diff_type)] <-
        ratio_dt$expr_status[is.na(ratio_dt$diff_type)]
    } else {
      ratio_dt$diff_type <- ratio_dt$expr_status
    }

    ratio_dt <- ratio_dt[order(-codon_ratio)]
    ratio_dt$codon_aa <- factor(ratio_dt$codon_aa,
                                levels = ratio_dt$codon_aa)

    text_color <- ifelse(grepl("[CG]$", as.character(ratio_dt$codon)),
                         "red", "black")

    p <- ggplot(ratio_dt, aes(x = codon_aa, y = codon_ratio, fill = diff_type)) +
      geom_col(alpha = 0.8, width = 0.8) +
      geom_errorbar(aes(ymin = codon_ratio - se_ratio,
                        ymax = codon_ratio + se_ratio),
                    width = 0.3, color = "black") +
      geom_hline(yintercept = 1, linetype = "solid",
                 color = "black", linewidth = 1.2, alpha = 0.8) +
      scale_fill_manual(values = c("up" = "#E63946",
                                   "down" = "#457B9D",
                                   "non_diff" = "#95A5A6"),
                        labels = c("Up-regulated (ratio>1)",
                                   "Down-regulated (ratio<1)",
                                   "Not significant (ratio=1)"),
                        breaks = c("up", "down", "non_diff")) +
      labs(x = paste0(ctype, " - Codon (Amino acid) Sorted by Ratio"),
           y = paste0("Codon Occupancy Ratio (",
                      group_names[2], "/", group_names[1], ")"),
           fill = "Expression Change",
           title = paste0(ctype, " - Codon Occupancy Ratio Histogram")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                       size = 9, color = text_color))

    ggsave(file.path(output_dir, paste0("4.", ctype, "_codon_ratio_hist.pdf")),
           p, width = 18, height = 7, dpi = 300)

    syn_dt <- ratio_dt[order(aa, -codon_ratio)]
    syn_dt$codon <- factor(syn_dt$codon, levels = unique(syn_dt$codon))
    syn_dt$x_pos <- as.integer(syn_dt$codon)
    syn_dt$x_pos <- syn_dt$x_pos +
      c(0, cumsum(diff(as.integer(factor(syn_dt$aa))) != 0)) * 0.8

    syn_color <- ifelse(grepl("[CG]$", as.character(syn_dt$codon)),
                        "red", "black")

    p_syn <- ggplot(syn_dt, aes(x = x_pos, y = codon_ratio, fill = diff_type)) +
      geom_col(alpha = 0.8, width = 0.8) +
      geom_errorbar(aes(ymin = codon_ratio - se_ratio,
                        ymax = codon_ratio + se_ratio),
                    width = 0.3, color = "black") +
      geom_hline(yintercept = 1, linetype = "solid",
                 color = "black", linewidth = 1.2, alpha = 0.8) +
      scale_fill_manual(values = c("up" = "#E63946",
                                   "down" = "#457B9D",
                                   "non_diff" = "#95A5A6")) +
      scale_x_continuous(breaks = syn_dt$x_pos,
                         labels = syn_dt$codon_aa,
                         expand = c(0.01, 0.01)) +
      labs(x = paste0(ctype, " - Codon (Grouped by Synonym)"),
           y = paste0("Codon Occupancy Ratio (",
                      group_names[2], "/", group_names[1], ")"),
           fill = "Expression Change",
           title = paste0(ctype, " - Synonymous Codon Occupancy Ratio")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                       size = 9, color = syn_color),
            panel.grid = element_blank())

    ggsave(file.path(output_dir, paste0("4.", ctype, "_codon_ratio_synonymous.pdf")),
           p_syn, width = 20, height = 7, dpi = 300)
  }

  cat("A/P/E ratio histograms generated\n\n")
}

cat("9. Filter Codon Data (retain codons present in all samples)\n")
cat("------------------------------------------------------------\n")

codon_filter <- all_codon_occupancy_raw[, list(n_sample = length(unique(sample_name))),
                                        by = list(codon_type, codon)]
codon_keep <- codon_filter[n_sample == length(sample_names)]

all_codon_occupancy_raw_filt <-
  all_codon_occupancy_raw[codon_type %in% codon_keep$codon_type &
                            codon %in% codon_keep$codon]
all_codon_occupancy_group_filt <-
  all_codon_occupancy_group[codon_type %in% codon_keep$codon_type &
                              codon %in% codon_keep$codon]

cat(paste0("Codon filter completed: ", nrow(codon_keep),
           " codons retained\n\n"))

if (length(sample_names) >= 3) {
  plot_Psite_PCA(all_codon_occupancy_raw_filt, group_dt, output_root)
} else {
  cat("Sample count =", length(sample_names),
      "(1vs1 comparison) -> skip PCA\n\n")
}

plot_group_footprint(all_psite_raw_dt, output_root)
plot_gene_specific_footprint(all_psite_raw_dt, gene_list_path,
                             tx2gene, "./5.ribo_group_diff_output")
plot_APE_codon_ratio_hist(all_codon_occupancy_group_filt, NULL, output_root)

get_codon_aa_map <- function() {
  list(
    "TTT" = "F", "TTC" = "F",
    "TTA" = "L", "TTG" = "L", "CTT" = "L", "CTC" = "L",
    "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I",
    "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "AGT" = "S", "AGC" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y",
    "TAA" = "*", "TAG" = "*", "TGA" = "*",
    "CAT" = "H", "CAC" = "H",
    "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N",
    "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D",
    "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C",
    "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
}

get_codon_color <- function(codon_aa_map) {
  aa_main_color <- list(
    "F" = "#E41A1C", "L" = "#F28E8E", "I" = "#FAAA3A", "M" = "#F7E14F",
    "V" = "#A65628", "S" = "#C49A6C", "P" = "#FF7F00", "T" = "#FDB462",
    "A" = "#377EB8", "Y" = "#75B7E6", "*" = "#984EA3", "H" = "#C17BDB",
    "Q" = "#4DAF4A", "N" = "#92CD76", "K" = "#1B9E77", "D" = "#66C2A5",
    "E" = "#8DD3C7", "C" = "#FFFFB3", "W" = "#BEBADA", "R" = "#FB8072",
    "G" = "#80B1D3"
  )

  get_shade_color <- function(base_color, n, idx) {
    rgb_vals <- tryCatch({
      col2rgb(base_color) / 255
    }, error = function(e) {
      matrix(c(0.5, 0.5, 0.5), nrow = 3)
    })
    if (n <= 1) return(base_color)
    lightness <- 0.8 + (idx - 1) * (0.4 / max(n - 1, 1))
    new_rgb <- pmin(pmax(rgb_vals * lightness, 0), 1)
    new_rgb[is.na(new_rgb) | is.infinite(new_rgb)] <- 0.5
    return(rgb(new_rgb[1], new_rgb[2], new_rgb[3]))
  }

  codon_colors <- c()
  aa_codons <- split(names(codon_aa_map), unlist(codon_aa_map))

  for (aa in names(aa_codons)) {
    codons <- aa_codons[[aa]]
    n <- length(codons)
    base_color <- aa_main_color[[aa]]
    if (is.null(base_color)) base_color <- "#888888"
    for (i in seq_along(codons)) {
      codon_colors[codons[i]] <- get_shade_color(base_color, n, i)
    }
  }
  return(codon_colors)
}

get_aa_to_class <- function() {
  list(
    "D" = "acidic", "E" = "acidic",
    "K" = "basic", "R" = "basic", "H" = "basic",
    "N" = "polar", "Q" = "polar", "S" = "polar", "T" = "polar", "C" = "polar",
    "F" = "aromatic", "W" = "aromatic", "Y" = "aromatic",
    "G" = "aliphatic", "A" = "aliphatic", "V" = "aliphatic",
    "L" = "aliphatic", "I" = "aliphatic", "M" = "aliphatic", "P" = "aliphatic",
    "*" = "stop"
  )
}

get_class_color_map <- function() {
  list(
    "acidic"   = "#E41A1C",
    "aliphatic" = "#000000",
    "aromatic" = "#4DAF4A",
    "basic"    = "#377EB8",
    "polar"    = "#FF7F00",
    "stop"     = "#984EA3"
  )
}

codon_aa_map    <- get_codon_aa_map()
codon_colors    <- get_codon_color(codon_aa_map)
aa_to_class     <- get_aa_to_class()
class_color_map <- get_class_color_map()

cat("10. Analyze Codon Motifs\n")
cat("------------------------------------------------------------\n")

analyze_codon_motif <- function(sample_dir, transcriptome_fasta, output_dir,
                                gene_list_path, tx2gene,
                                motif_type = c("5codon", "20codon")) {
  library(Biostrings)
  library(data.table)

  target_genes <- character()
  if (file.exists(gene_list_path)) {
    target_genes <- readLines(gene_list_path, warn = FALSE)
    target_genes <- gsub("\\s+", "", target_genes)
    target_genes <- unique(target_genes[target_genes != ""])
    cat(paste0("Target genes loaded: ", length(target_genes), "\n"))
  } else {
    cat("Gene list file not found: ", gene_list_path, "\n")
  }

  target_tx_all <- character()
  if (length(target_genes) > 0) {
    for (gene in target_genes) {
      tx_for_gene <- tx2gene[gene_name == gene, unique(transcript)]
      target_tx_all <- c(target_tx_all, tx_for_gene)
    }
    target_tx_all <- unique(target_tx_all)
    cat(paste0("Target transcripts: ", length(target_tx_all), "\n"))
  }

  psite_path <- file.path(sample_dir, "reads_psite.rds")
  if (!file.exists(psite_path)) {
    psite_path <- file.path(dirname(sample_dir), "reads_psite.rds")
  }

  psite_data <- readRDS(psite_path)
  if (is(psite_data, "list")) {
    psite_dt <- NULL
    for (elem in psite_data) {
      if (is(elem, "data.frame") || is(elem, "data.table")) {
        psite_dt <- as.data.table(elem)
        break
      }
    }
    if (is.null(psite_dt)) stop("Valid data.table not found in RDS list")
  } else {
    psite_dt <- as.data.table(psite_data)
  }

  required_cols <- c("cds_start", "cds_stop", "psite", "P_codon", "transcript")
  missing_cols <- setdiff(required_cols, colnames(psite_dt))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing required columns: ",
                paste(missing_cols, collapse = ", ")))
  }

  rpf_tx <- unique(psite_dt$transcript)
  target_tx <- intersect(target_tx_all, rpf_tx)
  cat(paste0("Target transcripts present in RPF data: ",
             length(target_tx), "\n"))

  tx_seqs <- readDNAStringSet(transcriptome_fasta)
  tx_names_full <- names(tx_seqs)
  tx_names_clean <- sapply(strsplit(tx_names_full, "\\|"), `[`, 1)

  parse_fasta_header <- function(header) {
    cds_match <- regmatches(header, regexpr("CDS:\\d+-\\d+", header))
    if (length(cds_match) == 0) {
      return(list(cds_start = NA, cds_stop = NA))
    }
    cds_range <- gsub("CDS:", "", cds_match)
    parts <- strsplit(cds_range, "-")[[1]]
    return(list(cds_start = as.integer(parts[1]),
                cds_stop  = as.integer(parts[2])))
  }

  seq_meta <- lapply(tx_names_full, parse_fasta_header)
  seq_map <- data.table(
    transcript       = tx_names_clean,
    seq              = as.character(tx_seqs),
    seq_len          = width(tx_seqs),
    fasta_cds_start  = sapply(seq_meta, function(x) x$cds_start),
    fasta_cds_stop   = sapply(seq_meta, function(x) x$cds_stop)
  )
  setkey(seq_map, transcript)

  target_tx <- intersect(target_tx, seq_map$transcript)
  cat(paste0("Target transcripts present in FASTA: ",
             length(target_tx), "\n"))

  motif_config <- list()
  if ("5codon" %in% motif_type) {
    motif_config[["5codon"]] <- list(
      flank_up   = 6,
      flank_down = 8,
      total_len  = 15,
      nt_labels  = -6:8,
      codon_labels = -2:2,
      site_pos   = 7,
      suffix     = "5codons"
    )
  }
  if ("20codon" %in% motif_type) {
    motif_config[["20codon"]] <- list(
      flank_up   = 60,
      flank_down = 2,
      total_len  = 62,
      nt_labels  = -60:-1,
      codon_labels = -20:-1,
      site_pos   = 61,
      suffix     = "20codons_upstream"
    )
  }

  plot_aa_class_stacked <- function(sorted_freq_list, positions, config,
                                    codon_output_dir, seq_type, codon) {
    plot_base_name <- paste0(codon, "_", config$suffix, "_",
                             seq_type, "_AA_class_color")
    png_file  <- file.path(codon_output_dir, "plots",
                           paste0(plot_base_name, ".png"))
    pdf_file  <- file.path(codon_output_dir, "plots",
                           paste0(plot_base_name, ".pdf"))
    table_file <- file.path(codon_output_dir, "tables",
                            paste0(plot_base_name, "_freq_table.csv"))

    plot_width <- ifelse(config$total_len > 30, 22, 14)
    plot_height <- 8

    aa_freq_table <- data.table()
    for (pos in positions) {
      if (!pos %in% names(sorted_freq_list) ||
          length(sorted_freq_list[[pos]]) == 0) next

      aa_names <- sapply(names(sorted_freq_list[[pos]]), function(c) {
        aa <- codon_aa_map[[c]]
        if (is.null(aa)) "*" else aa
      })
      aa_classes <- sapply(aa_names, function(aa) {
        cls <- aa_to_class[[aa]]
        if (is.null(cls)) "unknown" else cls
      })

      temp_dt <- data.table(
        position  = as.numeric(pos),
        codon     = names(sorted_freq_list[[pos]]),
        aa        = aa_names,
        frequency = sorted_freq_list[[pos]],
        class     = aa_classes
      )
      aa_freq_table <- rbind(aa_freq_table, temp_dt)
    }

    if (nrow(aa_freq_table) == 0) return(NULL)
    aa_freq_table <- aa_freq_table[class != "unknown"]
    if (nrow(aa_freq_table) == 0) return(NULL)
    aa_freq_table <- aa_freq_table[order(position)]

    write.csv(aa_freq_table, table_file, row.names = FALSE)

    positions_unique <- sort(unique(aa_freq_table$position))
    y_max <- max(aa_freq_table[, sum(frequency), by = position]$V1) * 1.1

    png(png_file, width = plot_width, height = plot_height,
        units = "in", res = 300)
    par(mar = c(10, 6, 4, 10) + 0.1)

    plot(0, 0, type = "n",
         xlim = c(0.5, length(positions_unique) + 0.5),
         ylim = c(0, y_max),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ",
                       config$suffix, " | AA by Chemical Class"),
         xaxt = "n", yaxt = "s",
         cex.lab = 1.1, cex.axis = 1.0)

    x_coords <- seq_along(positions_unique)
    names(x_coords) <- positions_unique

    for (pos in positions_unique) {
      pos_data <- aa_freq_table[position == pos][order(-frequency)]
      y_bottom <- 0
      for (i in seq_len(nrow(pos_data))) {
        aa   <- pos_data$aa[i]
        freq <- pos_data$frequency[i]
        cls  <- pos_data$class[i]
        color <- class_color_map[[cls]]
        if (is.null(color)) color <- "#888888"

        rect(x_coords[as.character(pos)] - 0.15, y_bottom,
             x_coords[as.character(pos)] + 0.15, y_bottom + freq,
             col = color, border = "white", lwd = 0.5)

        text_size <- ifelse(freq < 5, 0.6, 0.8)
        text(x_coords[as.character(pos)], y_bottom + freq / 2,
             labels = aa, cex = text_size, font = 2, col = "white")
        y_bottom <- y_bottom + freq
      }
    }

    axis(1, at = x_coords, labels = positions_unique, las = 1,
         cex.axis = 0.9, mgp = c(3, 0.5, 0))
    mtext("Position relative to PSite", side = 1, line = 7, cex = 1.0)

    legend(x = length(positions_unique) + 0.8, y = y_max / 2,
           legend = names(class_color_map),
           fill = unlist(class_color_map),
           title = "Chemical Class",
           cex = 0.9, bty = "n", xpd = TRUE, title.cex = 1.0)
    dev.off()

    pdf(pdf_file, width = plot_width, height = plot_height)
    par(mar = c(10, 6, 4, 10) + 0.1)

    plot(0, 0, type = "n",
         xlim = c(0.5, length(positions_unique) + 0.5),
         ylim = c(0, y_max),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ",
                       config$suffix, " | AA by Chemical Class"),
         xaxt = "n", yaxt = "s",
         cex.lab = 1.1, cex.axis = 1.0)

    for (pos in positions_unique) {
      pos_data <- aa_freq_table[position == pos][order(-frequency)]
      y_bottom <- 0
      for (i in seq_len(nrow(pos_data))) {
        aa   <- pos_data$aa[i]
        freq <- pos_data$frequency[i]
        cls  <- pos_data$class[i]
        color <- class_color_map[[cls]]
        if (is.null(color)) color <- "#888888"

        rect(x_coords[as.character(pos)] - 0.15, y_bottom,
             x_coords[as.character(pos)] + 0.15, y_bottom + freq,
             col = color, border = "white", lwd = 0.5)

        text_size <- ifelse(freq < 5, 0.6, 0.8)
        text(x_coords[as.character(pos)], y_bottom + freq / 2,
             labels = aa, cex = text_size, font = 2, col = "white")
        y_bottom <- y_bottom + freq
      }
    }

    axis(1, at = x_coords, labels = positions_unique, las = 1,
         cex.axis = 0.9, mgp = c(3, 0.5, 0))
    mtext("Position relative to PSite", side = 1, line = 7, cex = 1.0)

    legend(x = length(positions_unique) + 0.8, y = y_max / 2,
           legend = names(class_color_map),
           fill = unlist(class_color_map),
           title = "Chemical Class",
           cex = 0.9, bty = "n", xpd = TRUE, title.cex = 1.0)
    dev.off()

    return(TRUE)
  }

  plot_stacked_manual <- function(freq_table, type = c("nt", "codon"),
                                   seq_type, codon, config,
                                   codon_output_dir) {
    if (is.null(freq_table) || nrow(freq_table) == 0) return(NULL)

    freq_table$position <- as.character(freq_table$position)
    freq_table <- freq_table[!is.na(position) &
                               position != "" & frequency > 0]
    if (nrow(freq_table) == 0) return(NULL)

    if (type == "nt") {
      col_name   <- "nt"
      color_map  <- c("A" = "#4DAF4A", "T" = "#E41A1C",
                      "C" = "#377EB8", "G" = "#FFFF33")
      text_size  <- 0.8
      plot_width <- ifelse(config$suffix == "5codons", 12,
                           ifelse(config$total_len > 30, 15, 10))
      plot_height <- 6
      rect_width <- 0.4
    } else {
      col_name   <- "codon"
      color_map  <- codon_colors
      text_size  <- 0.5
      plot_width <- ifelse(config$suffix == "5codons", 14,
                           ifelse(config$total_len > 30, 20, 10))
      plot_height <- 7
      rect_width <- 0.2
    }

    positions <- sort(unique(freq_table$position))
    elements  <- unique(freq_table[[col_name]])

    freq_matrix <- matrix(0, nrow = length(positions), ncol = length(elements))
    rownames(freq_matrix) <- positions
    colnames(freq_matrix) <- elements

    for (i in seq_len(nrow(freq_table))) {
      pos  <- as.character(freq_table$position[i])
      elem <- as.character(freq_table[[col_name]][i])
      if (pos %in% rownames(freq_matrix) &&
          elem %in% colnames(freq_matrix)) {
        freq_matrix[pos, elem] <- freq_table$frequency[i]
      }
    }

    sorted_freq_list <- list()
    sorted_color_list <- list()
    for (pos in positions) {
      pos_freq <- freq_matrix[pos, ]
      pos_freq <- pos_freq[pos_freq > 0]
      if (length(pos_freq) == 0) next
      pos_freq_sorted <- sort(pos_freq, decreasing = TRUE)
      sorted_freq_list[[pos]] <- pos_freq_sorted
      sorted_color_list[[pos]] <- color_map[names(pos_freq_sorted)]
    }

    if (length(sorted_freq_list) == 0) return(NULL)

    plot_base_name <- paste0(codon, "_", config$suffix, "_",
                             seq_type, "_", type, "_manual")
    png_file  <- file.path(codon_output_dir, "plots",
                           paste0(plot_base_name, ".png"))
    pdf_file  <- file.path(codon_output_dir, "plots",
                           paste0(plot_base_name, ".pdf"))
    table_file <- file.path(codon_output_dir, "tables",
                            paste0(plot_base_name, "_freq_table.csv"))

    write.csv(freq_table, table_file, row.names = FALSE)

    png(png_file, width = plot_width, height = plot_height,
        units = "in", res = 300)
    par(mar = c(8, 5, 4, 4) + 0.1)

    plot(0, 0, type = "n",
         xlim = c(0.5, length(positions) + 0.5),
         ylim = c(0, 100),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ",
                       config$suffix, " | Frequency descending"),
         xaxt = "n", yaxt = "s")

    x_pos <- seq_along(positions)
    names(x_pos) <- positions
    y_bottom <- rep(0, length(positions))
    names(y_bottom) <- positions

    for (pos in positions) {
      pos_idx <- which(names(x_pos) == pos)
      if (length(pos_idx) == 0 ||
          !pos %in% names(sorted_freq_list)) next

      pos_freq <- sorted_freq_list[[pos]]
      pos_colors <- sorted_color_list[[pos]]
      if (length(pos_freq) == 0) next

      for (i in seq_along(pos_freq)) {
        freq <- pos_freq[i]
        color <- pos_colors[i]
        elem_name <- names(pos_freq)[i]

        rect(x_pos[pos_idx] - rect_width, y_bottom[pos_idx],
             x_pos[pos_idx] + rect_width, y_bottom[pos_idx] + freq,
             col = color, border = "white")

        if (type == "codon") {
          aa <- ifelse(elem_name %in% names(codon_aa_map),
                       codon_aa_map[[elem_name]], "")
          display_text <- ifelse(aa != "",
                                 paste0(elem_name, "(", aa, ")"),
                                 elem_name)
        } else {
          display_text <- elem_name
        }

        text(x_pos[pos_idx], y_bottom[pos_idx] + freq / 2,
             labels = display_text,
             cex = text_size, font = 2, col = "white")
        y_bottom[pos_idx] <- y_bottom[pos_idx] + freq
      }
    }

    axis(1, at = x_pos, labels = positions, las = 2, cex.axis = 0.7)
    mtext(ifelse(type == "nt",
                 "Position relative to PSite (nt)",
                 "Position relative to PSite (codon)"),
          side = 1, line = 6, cex = 0.9)
    dev.off()

    pdf(pdf_file, width = plot_width, height = plot_height)
    par(mar = c(8, 5, 4, 4) + 0.1)

    plot(0, 0, type = "n",
         xlim = c(0.5, length(positions) + 0.5),
         ylim = c(0, 100),
         xlab = "", ylab = "Frequency (%)",
         main = paste0("Codon: ", codon, " (", seq_type, ") - ",
                       config$suffix, " | Frequency descending"),
         xaxt = "n", yaxt = "s")

    y_bottom <- rep(0, length(positions))
    names(y_bottom) <- positions
    for (pos in positions) {
      pos_idx <- which(names(x_pos) == pos)
      if (length(pos_idx) == 0 ||
          !pos %in% names(sorted_freq_list)) next

      pos_freq <- sorted_freq_list[[pos]]
      pos_colors <- sorted_color_list[[pos]]
      if (length(pos_freq) == 0) next

      for (i in seq_along(pos_freq)) {
        freq <- pos_freq[i]
        color <- pos_colors[i]
        elem_name <- names(pos_freq)[i]

        rect(x_pos[pos_idx] - rect_width, y_bottom[pos_idx],
             x_pos[pos_idx] + rect_width, y_bottom[pos_idx] + freq,
             col = color, border = "white")

        if (type == "codon") {
          aa <- ifelse(elem_name %in% names(codon_aa_map),
                       codon_aa_map[[elem_name]], "")
          display_text <- ifelse(aa != "",
                                 paste0(elem_name, "(", aa, ")"),
                                 elem_name)
        } else {
          display_text <- elem_name
        }

        text(x_pos[pos_idx], y_bottom[pos_idx] + freq / 2,
             labels = display_text,
             cex = text_size, font = 2, col = "white")
        y_bottom[pos_idx] <- y_bottom[pos_idx] + freq
      }
    }

    axis(1, at = x_pos, labels = positions, las = 2, cex.axis = 0.7)
    mtext(ifelse(type == "nt",
                 "Position relative to PSite (nt)",
                 "Position relative to PSite (codon)"),
          side = 1, line = 6, cex = 0.9)
    dev.off()

    if (type == "codon") {
      plot_aa_class_stacked(sorted_freq_list, positions, config,
                            codon_output_dir, seq_type, codon)
    }

    return(TRUE)
  }

  all_codons <- sort(unique(psite_dt[!is.na(P_codon), P_codon]))
  cat(paste0("Total unique P-codons: ", length(all_codons), "\n"))

  for (codon_idx in seq_along(all_codons)) {
    codon <- all_codons[codon_idx]
    cat(paste0("\nProcessing codon ", codon_idx, "/",
               length(all_codons), ": ", codon, "\n"))

    codon_output_dir <- file.path(output_dir, paste0("codon_", codon))
    codon_plot_dir   <- file.path(codon_output_dir, "plots")
    codon_table_dir  <- file.path(codon_output_dir, "tables")
    dir.create(codon_plot_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(codon_table_dir, showWarnings = FALSE, recursive = TRUE)

    codon_dt <- psite_dt[!is.na(P_codon) & P_codon == codon,
                         .(transcript, psite, cds_start, cds_stop)]
    if (nrow(codon_dt) == 0) next

    codon_dt <- merge(codon_dt, seq_map, by = "transcript", all.x = TRUE)
    codon_dt <- codon_dt[!is.na(seq) & !is.na(seq_len)]
    if (nrow(codon_dt) == 0) next

    for (motif_type_name in names(motif_config)) {
      config <- motif_config[[motif_type_name]]
      cat(paste0("  Motif type: ", motif_type_name, "\n"))

      codon_dt[, psite_in_transcript := psite >= 1 & psite <= seq_len]
      codon_dt[, psite_in_cds := psite >= cds_start & psite <= cds_stop]
      codon_dt[, psite_codon_start := psite - ((psite - cds_start) %% 3)]
      codon_dt[, motif_start := psite - config$flank_up]
      codon_dt[, motif_end   := psite + config$flank_down]
      codon_dt[, motif_in_seq := (motif_start >= 1) & (motif_end <= seq_len)]
      codon_dt[, motif_in_cds := (motif_start >= cds_start) & (motif_end <= cds_stop)]
      codon_dt[, motif_valid := psite_in_transcript &
                 psite_in_cds & motif_in_seq & motif_in_cds]

      valid_dt <- codon_dt[motif_valid == TRUE]
      if (nrow(valid_dt) == 0) next

      target_dt <- valid_dt[transcript %in% target_tx]
      non_target_dt <- valid_dt[!transcript %in% target_tx]
      set.seed(123)
      random_n <- min(1000, nrow(non_target_dt))
      random_dt <- non_target_dt[sample(.N, random_n)]

      extract_sequences <- function(data_dt, seq_type, codon, config) {
        if (nrow(data_dt) < 1) return(NULL)

        extracted_seqs <- character(nrow(data_dt))
        for (i in seq_len(nrow(data_dt))) {
          m_start <- data_dt$motif_start[i]
          m_end   <- data_dt$motif_end[i]
          full_seq <- data_dt$seq[i]
          motif_seq <- substr(full_seq, m_start, m_end)
          if (nchar(motif_seq) > config$total_len) {
            motif_seq <- substr(motif_seq, 1, config$total_len)
          }
          extracted_seqs[i] <- motif_seq
        }

        valid_nt_seqs <- extracted_seqs[
          nchar(extracted_seqs) == config$total_len]
        if (length(valid_nt_seqs) == 0) return(NULL)

        nt_freq_table <- data.table()
        for (pos in seq_len(config$total_len)) {
          pos_nts <- substr(valid_nt_seqs, pos, pos)
          pos_freq <- as.data.frame(table(pos_nts))
          colnames(pos_freq) <- c("nt", "count")
          pos_freq$frequency <-
            round((pos_freq$count / sum(pos_freq$count)) * 100, 2)
          pos_freq$position <- config$nt_labels[pos]
          pos_freq <- pos_freq[order(-pos_freq$frequency)]
          nt_freq_table <- rbind(nt_freq_table, pos_freq)
        }

        codon_freq_table <- NULL
        if (grepl("codon", motif_type_name)) {
          codon_seqs_list <- lapply(valid_nt_seqs, function(seq) {
            codon_starts <- seq(1, nchar(seq) - 2, 3)
            codons <- substring(seq, codon_starts, codon_starts + 2)
            if (length(codons) == length(config$codon_labels)) {
              names(codons) <- config$codon_labels
            }
            return(codons)
          })

          codon_freq_table <- data.table()
          codon_pos_matrix <- do.call(rbind, codon_seqs_list)
          for (pos in colnames(codon_pos_matrix)) {
            pos_codons <- codon_pos_matrix[, pos]
            pos_freq <- as.data.frame(table(pos_codons))
            colnames(pos_freq) <- c("codon", "count")
            pos_freq$frequency <-
              round((pos_freq$count / sum(pos_freq$count)) * 100, 2)
            pos_freq$position <- pos
            pos_freq <- pos_freq[order(-pos_freq$frequency)]
            codon_freq_table <- rbind(codon_freq_table, pos_freq)
          }
        }

        return(list(nt_freq_table = nt_freq_table,
                    codon_freq_table = codon_freq_table))
      }

      if (nrow(target_dt) > 0) {
        target_result <- extract_sequences(target_dt, "target_genes",
                                            codon, config)
        if (!is.null(target_result)) {
          if (!is.null(target_result$nt_freq_table)) {
            plot_stacked_manual(target_result$nt_freq_table, "nt",
                               "target_genes", codon, config,
                               codon_output_dir)
          }
          if (!is.null(target_result$codon_freq_table)) {
            plot_stacked_manual(target_result$codon_freq_table, "codon",
                               "target_genes", codon, config,
                               codon_output_dir)
          }
        }
      }

      if (nrow(random_dt) > 0) {
        random_result <- extract_sequences(random_dt, "random",
                                           codon, config)
        if (!is.null(random_result)) {
          if (!is.null(random_result$nt_freq_table)) {
            plot_stacked_manual(random_result$nt_freq_table, "nt",
                               "random", codon, config,
                               codon_output_dir)
          }
          if (!is.null(random_result$codon_freq_table)) {
            plot_stacked_manual(random_result$codon_freq_table, "codon",
                               "random", codon, config,
                               codon_output_dir)
          }
        }
      }
    }

    rm(codon_dt, valid_dt, target_dt, random_dt)
    gc()
  }

  cat("All codon motif analysis finished\n")
}

if (exists("all_codon_occupancy_group_filt") &&
    nrow(all_codon_occupancy_group_filt) > 0) {
  cat(paste0("Transcript-gene mappings loaded: ",
             nrow(tx2gene), "\n"))
  first_sample_dir <- dirname(rds_paths[1])

  analyze_codon_motif(
    sample_dir         = first_sample_dir,
    transcriptome_fasta = transcriptome_fasta,
    output_dir         = output_root,
    gene_list_path     = gene_list_path,
    tx2gene            = tx2gene,
    motif_type         = c("5codon", "20codon")
  )
} else {
  cat("No filtered codon data available. Skip motif analysis\n")
}

cat("\n============================================================\n")
cat("ALL ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(paste0("Results saved to: ", output_root, "\n"))
cat("============================================================\n")

sink()
