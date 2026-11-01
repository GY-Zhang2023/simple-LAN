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

Sys.setenv(LANGUAGE='en')
options(stringsAsFactors = FALSE)
rm(list = ls())

RUN_PARAMS_PATH <- file.path(getwd(), "4.riboWaltz_output", "riboWaltz_batch_run_report.txt")

selected_species <- NULL
rpf_threshold <- 1e6

if (file.exists(RUN_PARAMS_PATH)) {
  con <- file(RUN_PARAMS_PATH, "r")
  param_lines <- c()
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (grepl("^species\\s*=|^rpf_threshold\\s*=", line)) {
      param_lines <- c(param_lines, line)
    }
  }
  close(con)
  
  species_line <- grep("^species\\s*=", param_lines, value = TRUE)
  if (length(species_line) > 0) {
    selected_species <- gsub("^species\\s*=\\s*|['\"]", "", species_line)
    selected_species <- trimws(selected_species)
  }
  
  threshold_line <- grep("^rpf_threshold\\s*=", param_lines, value = TRUE)
  if (length(threshold_line) > 0) {
    rpf_threshold <- as.numeric(gsub("^rpf_threshold\\s*=\\s*", "", threshold_line))
    if (is.na(rpf_threshold)) rpf_threshold <- 1e6
  }
}

if (is.null(selected_species) || !selected_species %in% c("human", "mouse")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1 || !args[1] %in% c("human", "mouse")) {
    cat("❌ Error: valid species parameter not found!\n")
    cat("Reason: parameter log malformed or no valid command line argument\n")
    cat("Usage1: Rscript this_script.R human|mouse\n")
    cat("Usage2: ensure riboWaltz log contains 'species=human' or 'species=mouse'\n")
    quit(save = "no", status = 1)
  }
  selected_species <- args[1]
  cat(paste0("⚠️ Using command line argument: species = ", selected_species, "\n"))
}

PATH_HUMAN_FASTA <- "~/reference/human_transcriptome/gencode.v49.pc_transcripts.fa"
PATH_MOUSE_FASTA <- "~/reference/mouse_transcriptome/gencode.vM35.pc_transcripts.fa"
species_fasta_map <- list(
  human = PATH_HUMAN_FASTA,
  mouse = PATH_MOUSE_FASTA
)

if (!selected_species %in% names(species_fasta_map)) {
  cat(paste0("❌ Error: unsupported species '", selected_species, "'\n"))
  cat("Supported species: ", paste(names(species_fasta_map), collapse = ", "), "\n")
  quit(save = "no", status = 1)
}
TRANSCRIPTOME_FASTA <- species_fasta_map[[selected_species]]

cat(paste0("\n✅ Loaded config for 【", selected_species, "】 (inherited from riboWaltz):\n"))
cat(paste0("   Transcriptome FASTA: ", TRANSCRIPTOME_FASTA, "\n"))
cat(paste0("   RPF threshold: ", format(rpf_threshold, scientific = FALSE), "\n\n"))

N_CORES <- 8
N_THREADS <- 2
MEMORY_LIMIT <- "50G"

options(
  stringsAsFactors = FALSE,
  datatable.verbose = FALSE,
  datatable.nthreads = N_THREADS,
  mc.cores = N_CORES
)
options(future.globals.maxSize = 50 * 1024 * 1024 * 1024)
options(future.rng.onMisuse = "ignore")
WORK_DIR <- getwd()
OUTPUT_ROOT <- paste0(WORK_DIR,"/4.riboWaltz_output")
MOTIF_OUTPUT_ROOT <- file.path(OUTPUT_ROOT, "1.motif")
TARGET_SITE <- "P"
MOTIF_OUTPUT_ROOT_ASITE <- file.path(OUTPUT_ROOT, "2.motif_asite")
TARGET_SITE_ASITE <- "A"
META_PEAKS_CENTER_ASITE <- "A_codon"
META_PEAKS_RANGE_ASITE <- c(-30, 30)
MOTIF_RANGE <- 9
AA_MOTIF_CODON_NUM <- 4

DISABLE_TOP3 <- FALSE
BATCH_SIZE <- 3
CHUNK_SIZE <- 5e5

META_PEAKS_CENTER <- "P_codon"
META_PEAKS_RANGE <- c(-30, 30)

if (!require("dplyr")) {
  install.packages("dplyr", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(dplyr)
}
if (!require("future")) {
  install.packages("future", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(future)
}
if (!require("furrr")) {
  install.packages("furrr", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(furrr)
}
if (!require("purrr")) {
  install.packages("purrr", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(purrr)
}
if (!require("BiocParallel")) {
  if(!require("BiocManager")) install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  BiocManager::install("BiocParallel", update = F)
  library(BiocParallel)
}
library(data.table)
library(ggplot2)
library(Biostrings)
library(utils)
library(R.utils)
library(ggseqlogo)

system(paste0("ulimit -v ", as.numeric(gsub("G", "", MEMORY_LIMIT)) * 1024 * 1024 * 1024))
if (.Platform$OS.type == "unix") {
  plan(multicore, workers = N_CORES)
} else {
  plan(multisession, workers = N_CORES)
}

codon2aa <- c(
  "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
  "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
  "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*",
  "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",
  "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
  "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
  "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
  "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
  "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
  "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
  "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
  "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
  "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
  "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
  "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
  "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G",
  "---"="-"
)
codon2aa_hash <- new.env(hash = TRUE, size = 64L)
for (codon in names(codon2aa)) {
  codon2aa_hash[[codon]] <- codon2aa[[codon]]
}

dna2aa_seq_fast <- function(dna_seqs) {
  valid_dna <- dna_seqs[!is.na(dna_seqs) & nchar(dna_seqs) %% 3 == 0 & nchar(dna_seqs) > 0]
  if (length(valid_dna) == 0) return(character(0))
  
  codons_list <- lapply(valid_dna, function(seq) {
    substring(seq, seq(1, nchar(seq)-2, 3), seq(3, nchar(seq), 3))
  })
  
  aa_seqs <- lapply(codons_list, function(codons) {
    aa <- unlist(mget(codons, envir = codon2aa_hash, ifnotfound = "X"))
    paste(aa, collapse = "")
  }) %>% unlist()
  
  aa_seqs <- aa_seqs[nchar(aa_seqs) >= 3]
  return(aa_seqs)
}

assign("dna2aa_seq_fast", dna2aa_seq_fast, envir = .GlobalEnv)
assign("codon2aa_hash", codon2aa_hash, envir = .GlobalEnv)

read_ribo_rds_fast <- function(sample_dir) {
  rds_path <- file.path(sample_dir, "reads_psite.rds")
  if (!file.exists(rds_path)) return(NULL)
  
  file_size <- file.info(rds_path)$size / 1024 / 1024
  if (file_size > 5000) {
    warning(paste("Sample", basename(sample_dir), "too large (>5G), skip"))
    return(NULL)
  }
  
  psite_data <- readRDS(rds_path)
  sample_name <- names(psite_data)[1]
  psite_dt <- as.data.table(psite_data[[sample_name]])
  
  unique_cols <- unique(colnames(psite_dt))
  psite_dt <- psite_dt[, ..unique_cols]
  
  if (nrow(psite_dt) > CHUNK_SIZE) {
    psite_dt <- psite_dt[1:CHUNK_SIZE]
  }
  
  return(list(sample_name = sample_name, psite_dt = psite_dt))
}

load_transcriptome_fast <- function(fasta_path) {
  if (!file.exists(fasta_path)) stop(paste("FASTA not found:", fasta_path))
  transcriptome <- readDNAStringSet(fasta_path, use.names = TRUE, format = "fasta")
  fasta_enst_ids <- sub("\\|.*$", "", names(transcriptome))
  fasta_seqs <- as.character(transcriptome)
  
  fasta_map <- data.table(
    transcript = fasta_enst_ids,
    seq_str = fasta_seqs
  )
  fasta_map[, transcript := as.factor(transcript)]
  setkey(fasta_map, transcript)
  
  rm(transcriptome); invisible(gc())
  return(fasta_map)
}

calculate_codon_pos_fast <- function(psite_dt, target_site) {
  pos_col <- switch(target_site, P = "psite", A = "A_site", E = "E_site")
  
  psite_dt[, cds_pos := get(pos_col) - cds_start + 1]
  psite_dt[, cds_codon_pos := floor((cds_pos - 1)/3) + 1]
  psite_dt[, total_codon_num := floor((cds_stop - cds_start + 1)/3), by = transcript]
  psite_dt[, stop_codon_pos := total_codon_num - cds_codon_pos]
  
  psite_dt <- psite_dt[!is.na(cds_codon_pos) & cds_pos > 0]
  return(psite_dt)
}

extract_motif_seqs_fast <- function(psite_dt, fasta_map, target_site, motif_range, aa_codon_num) {
  site_config <- list(
    pos_col = switch(target_site, P = "psite", A = "A_site", E = "E_site"),
    valid_col = switch(target_site, P = "P_site_valid", A = "A_site_valid", E = "E_site_valid")
  )
  
  psite_dt[, `:=`(
    motif_start = get(site_config$pos_col) - motif_range,
    motif_end = get(site_config$pos_col) + motif_range
  )]
  
  psite_dt[, `:=`(
    raw_motif_seq = NA_character_,
    dna_motif_seq = NA_character_,
    aa_motif_seq = NA_character_
  )]
  
  valid_idx <- which(psite_dt[[site_config$valid_col]] &
                       !is.na(psite_dt$seq_str) &
                       psite_dt$motif_start >= 1 &
                       psite_dt$motif_end <= nchar(psite_dt$seq_str))
  
  if (length(valid_idx) > 0) {
    psite_dt[valid_idx, raw_motif_seq := lapply(
      X = valid_idx,
      FUN = function(i) {
        substr(psite_dt$seq_str[i], psite_dt$motif_start[i], psite_dt$motif_end[i])
      }
    ) %>% unlist()]
    
    dna_valid_idx <- which(psite_dt[[site_config$valid_col]] &
                             !is.na(psite_dt$seq_str) &
                             (psite_dt[[site_config$pos_col]] - 6) >= 1 &
                             (psite_dt[[site_config$pos_col]] + 6) <= nchar(psite_dt$seq_str))
    psite_dt[dna_valid_idx, dna_motif_seq := substr(
      seq_str,
      start = get(site_config$pos_col) - 6,
      stop = get(site_config$pos_col) + 6
    )]
    
    aa_valid_idx <- which(psite_dt[[site_config$valid_col]] &
                           !is.na(psite_dt$seq_str) &
                           !is.na(psite_dt$raw_motif_seq) &
                           nchar(psite_dt$raw_motif_seq) == 19)
    psite_dt[aa_valid_idx, aa_motif_seq := substr(raw_motif_seq, 7, 18)]
  }
  
  valid_dna_seqs <- unique(psite_dt[nchar(dna_motif_seq) == 13 & !is.na(dna_motif_seq), dna_motif_seq])
  valid_aa_seqs <- unique(psite_dt[nchar(aa_motif_seq) == 12 & !is.na(aa_motif_seq), aa_motif_seq])
  
  return(list(
    psite_dt = psite_dt,
    valid_dna_seqs = valid_dna_seqs,
    valid_aa_seqs = valid_aa_seqs
  ))
}

plot_motif_logo_fast <- function(valid_seqs, sample_name, target_site, is_amino = FALSE, is_stop = FALSE) {
  if (length(valid_seqs) < 10) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No valid sequences", size = 4) +
             ggtitle(paste0(sample_name, " ", target_site, "-site ", ifelse(is_amino, "AA ", "DNA "), "Motif")) +
             theme_void() +
             theme(plot.title = element_text(size = 10, hjust = 0.5)))
  }
  
  if (is_amino) {
    if (is_stop) {
      x_labels <- -2:1
      site_pos_in_motif <- 3
    } else {
      x_labels <- -1:2
      site_pos_in_motif <- 2
    }
    seq_length <- length(x_labels)
    
    plot <- ggseqlogo(valid_seqs, seq_type = "aa", col_scheme = "chemistry", font = "helvetica_bold") +
      labs(
        title = paste0(sample_name, " ", target_site, "-site AA Motif (", ifelse(is_stop, "Stop Region", "Normal Region"), ")"),
        x = "Position relative to P-site codon (AA)",
        y = "Frequency"
      ) +
      geom_vline(xintercept = site_pos_in_motif - 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = margin(5,5,5,5)
      ) +
      scale_x_continuous(
        breaks = 1:seq_length,
        labels = x_labels,
        expand = c(0.01, 0.01)
      )
  } else {
    x_labels <- -6:6
    site_pos_in_motif <- 7
    seq_length <- length(x_labels)
    
    plot <- ggseqlogo(valid_seqs, seq_type = "dna", col_scheme = "nucleotide", font = "helvetica_bold") +
      labs(
        title = paste0(sample_name, " ", target_site, "-site DNA Motif (", ifelse(is_stop, "Stop Region", "Normal Region"), ")"),
        x = "Position relative to P-site codon (nt)",
        y = "Frequency"
      ) +
      geom_vline(xintercept = site_pos_in_motif - 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = margin(5,5,5,5)
      ) +
      scale_x_continuous(
        breaks = 1:seq_length,
        labels = x_labels,
        expand = c(0.01, 0.01)
      )
  }
  return(plot)
}

extract_codon_metapeak_counts <- function(psite_dt, fasta_map) {
  valid_dt <- psite_dt[P_site_valid & !is.na(P_codon) & nchar(P_codon)==3 & !is.na(seq_str) & !is.na(cds_start) & !is.na(cds_stop)]
  if(nrow(valid_dt) == 0) {
    warning("No valid data for MetaPeak analysis")
    return(list())
  }
  
  all_codons <- unique(valid_dt$P_codon)
  codon_count_list <- list()
  
  for (target_codon in all_codons) {
    codon_dt <- valid_dt[P_codon == target_codon]
    if(nrow(codon_dt) < 5) next
    
    target_transcripts <- unique(codon_dt$transcript)
    tx_num <- length(target_transcripts)
    if(tx_num == 0) next
    
    codon_dt[, codon_pos := psite]
    
    meta_density_dt <- data.table()
    for (tx in target_transcripts) {
      tx_all_psites <- valid_dt[transcript == tx & P_site_valid & !is.na(psite)]
      if(nrow(tx_all_psites) < 3) next
      
      tx_codon_pos <- mean(codon_dt[transcript == tx]$codon_pos)
      tx_all_psites[, relative_pos := psite - tx_codon_pos]
      tx_all_psites <- tx_all_psites[relative_pos >= META_PEAKS_RANGE[1] & relative_pos <= META_PEAKS_RANGE[2]]
      meta_density_dt <- rbind(meta_density_dt, tx_all_psites[, .(relative_pos)])
    }
    
    if(nrow(meta_density_dt) == 0) next
    
    density_stats <- meta_density_dt[, .(raw_count = .N), by = .(relative_pos)]
    density_stats[, avg_density := raw_count / tx_num]
    
    full_pos <- data.table(relative_pos = META_PEAKS_RANGE[1]:META_PEAKS_RANGE[2])
    density_stats <- merge(full_pos, density_stats, by = "relative_pos", all.x = TRUE)
    density_stats[is.na(avg_density), avg_density := 0]
    density_stats[is.na(raw_count), raw_count := 0]
    
    attr(density_stats, "tx_num") <- tx_num
    codon_count_list[[target_codon]] <- density_stats
  }
  
  return(codon_count_list)
}

extract_codon_metapeak_counts_asite <- function(psite_dt, fasta_map) {
  valid_dt <- psite_dt[A_site_valid & !is.na(A_codon) & nchar(A_codon)==3 & !is.na(seq_str) & !is.na(cds_start) & !is.na(cds_stop)]
  if(nrow(valid_dt) == 0) {
    warning("No valid A-site data for MetaPeak analysis")
    return(list())
  }
  
  all_codons <- unique(valid_dt$A_codon)
  codon_count_list <- list()
  
  for (target_codon in all_codons) {
    codon_dt <- valid_dt[A_codon == target_codon]
    if(nrow(codon_dt) < 5) next
    
    target_transcripts <- unique(codon_dt$transcript)
    tx_num <- length(target_transcripts)
    if(tx_num == 0) next
    
    codon_dt[, codon_pos := A_site]
    
    meta_density_dt <- data.table()
    for (tx in target_transcripts) {
      tx_all_asites <- valid_dt[transcript == tx & A_site_valid & !is.na(A_site)]
      if(nrow(tx_all_asites) < 3) next
      
      tx_codon_pos <- mean(codon_dt[transcript == tx]$codon_pos)
      tx_all_asites[, relative_pos := A_site - tx_codon_pos]
      tx_all_asites <- tx_all_asites[relative_pos >= META_PEAKS_RANGE_ASITE[1] & relative_pos <= META_PEAKS_RANGE_ASITE[2]]
      meta_density_dt <- rbind(meta_density_dt, tx_all_asites[, .(relative_pos)])
    }
    
    if(nrow(meta_density_dt) == 0) next
    
    density_stats <- meta_density_dt[, .(raw_count = .N), by = .(relative_pos)]
    density_stats[, avg_density := raw_count / tx_num]
    
    full_pos <- data.table(relative_pos = META_PEAKS_RANGE_ASITE[1]:META_PEAKS_RANGE_ASITE[2])
    density_stats <- merge(full_pos, density_stats, by = "relative_pos", all.x = TRUE)
    density_stats[is.na(avg_density), avg_density := 0]
    density_stats[is.na(raw_count), raw_count := 0]
    
    attr(density_stats, "tx_num") <- tx_num
    codon_count_list[[target_codon]] <- density_stats
  }
  return(codon_count_list)
}

save_metapeak_counts <- function(count_list, sample_name, out_dir) {
  if(length(count_list) == 0) return(NULL)
  
  walk(names(count_list), function(codon) {
    count_dt <- count_list[[codon]]
    aa <- codon2aa[codon]
    tx_num <- attr(count_dt, "tx_num")
    total_psites <- sum(count_dt$raw_count, na.rm = TRUE)
    
    count_dt$sample <- sample_name
    count_dt$codon <- codon
    count_dt$aa <- aa
    count_dt$transcript_number <- tx_num
    count_dt$total_psites <- total_psites
    
    fwrite(
      count_dt,
      file.path(out_dir, paste0("7_MetaPeak_", codon, "_counts.txt")),
      sep = "\t",
      quote = FALSE
    )
  })
}

plot_codon_metapeak <- function(count_list, sample_name, out_dir) {
  if(length(count_list) == 0) return(NULL)
  
  save_metapeak_counts(count_list, sample_name, out_dir)
  
  walk(names(count_list), function(codon) {
    count_dt <- count_list[[codon]]
    aa <- codon2aa[codon]
    
    tx_num <- attr(count_dt, "tx_num")
    total_psites <- sum(count_dt$raw_count, na.rm = TRUE)
    
    plot_title <- paste0(sample_name, " | MetaPeak: ", codon, " (", aa, ")\n",
                         "Transcripts = ", tx_num, " | Total P-sites = ", total_psites)
    
    p <- ggplot(count_dt, aes(x = relative_pos, y = avg_density)) +
      geom_line(color = "#2c7fb8", linewidth = 1.1, alpha = 0.9) +
      geom_point(color = "#2c7fb8", size = 1.2, alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#d95f02", linewidth = 1.2, alpha = 0.9) +
      labs(
        title = plot_title,
        x = "Relative position to target codon (nt)",
        y = "Average P-site density (per transcript)"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10, face = "bold"),
        plot.margin = margin(8,8,8,8),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5)
      ) +
      scale_x_continuous(
        limits = META_PEAKS_RANGE,
        breaks = seq(-30, 30, 5),
        expand = c(0.01, 0.01)
      ) +
      scale_y_continuous(expand = c(0, 0))
    
    ggsave(
      filename = file.path(out_dir, paste0("7_MetaPeak_", codon, "_avg_density.pdf")),
      plot = p, width = 10, height = 4, dpi = 300, device = "pdf", bg = "white"
    )
  })
  invisible(gc())
}

plot_codon_metapeak_asite <- function(count_list, sample_name, out_dir) {
  if(length(count_list) == 0) return(NULL)
  save_metapeak_counts(count_list, sample_name, out_dir)
  
  walk(names(count_list), function(codon) {
    count_dt <- count_list[[codon]]
    aa <- codon2aa[codon]
    tx_num <- attr(count_dt, "tx_num")
    total_psites <- sum(count_dt$raw_count, na.rm = TRUE)
    
    plot_title <- paste0(sample_name, " | A-site MetaPeak: ", codon, " (", aa, ")\n",
                         "Transcripts = ", tx_num, " | Total A-sites = ", total_psites)
    
    p <- ggplot(count_dt, aes(x = relative_pos, y = avg_density)) +
      geom_line(color = "#d95f02", linewidth = 1.1, alpha = 0.9) +
      geom_point(color = "#d95f02", size = 1.2, alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#2c7fb8", linewidth = 1.2) +
      labs(x = "Relative position to A-codon (nt)", y = "Average A-site density") +
      theme_bw() +
      scale_x_continuous(limits = META_PEAKS_RANGE_ASITE, breaks = seq(-30,30,5))
    
    ggsave(file.path(out_dir, paste0("7_MetaPeak_A_", codon, ".pdf")),
           p, width=10, height=4, dpi=300, bg="white")
  })
  invisible(gc())
}

process_sample_fast <- function(sample_dir) {
  fasta_map <- get("fasta_map", envir = .GlobalEnv)
  sample_name <- basename(sample_dir)
  sample_motif_dir <- file.path(MOTIF_OUTPUT_ROOT, sample_name)
  dir.create(sample_motif_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    psite_result <- read_ribo_rds_fast(sample_dir)
    if (is.null(psite_result)) return(paste(sample_name, "read failed"))
    psite_dt <- psite_result$psite_dt
    
    psite_dt <- calculate_codon_pos_fast(psite_dt, TARGET_SITE)
    
    psite_dt <- merge(psite_dt, fasta_map, by = "transcript", all.x = TRUE)
    psite_dt <- psite_dt[!is.na(seq_str)]
    if (nrow(psite_dt) ==0) return(paste(sample_name, "no sequence"))
    
    general_motif <- extract_motif_seqs_fast(psite_dt, fasta_map, TARGET_SITE, MOTIF_RANGE, AA_MOTIF_CODON_NUM)
    ggsave(
      file.path(sample_motif_dir, "1_P_site_general_motif.pdf"),
      plot = plot_motif_logo_fast(general_motif$valid_dna_seqs, sample_name, TARGET_SITE, FALSE),
      width=8, height=3, dpi=200, device="pdf", bg="white"
    )
    aa_seqs <- dna2aa_seq_fast(general_motif$valid_aa_seqs)
    ggsave(
      file.path(sample_motif_dir, "1_P_site_general_aminomotif.pdf"),
      plot = plot_motif_logo_fast(aa_seqs, sample_name, TARGET_SITE, TRUE),
      width=8, height=3, dpi=200, device="pdf", bg="white"
    )
    
    start_psite_dt <- psite_dt[
      cds_codon_pos ==1 &
        !is.na(cds_codon_pos) &
        P_site_valid
    ]
    start_motif <- extract_motif_seqs_fast(start_psite_dt, fasta_map, TARGET_SITE, MOTIF_RANGE, AA_MOTIF_CODON_NUM)
    ggsave(
      file.path(sample_motif_dir, "2_P_site_start_region_motif.pdf"),
      plot = plot_motif_logo_fast(start_motif$valid_dna_seqs, sample_name, TARGET_SITE, FALSE),
      width=8, height=3, dpi=200, device="pdf", bg="white"
    )
    start_aa_seqs <- dna2aa_seq_fast(start_motif$valid_aa_seqs)
    ggsave(
      file.path(sample_motif_dir, "2_P_site_start_region_aminomotif.pdf"),
      plot = plot_motif_logo_fast(start_aa_seqs, sample_name, TARGET_SITE, TRUE),
      width=8, height=3, dpi=200, device="pdf", bg="white"
    )
    
    stop_psite_dt <- psite_dt[
      stop_codon_pos ==1 & !is.na(stop_codon_pos) & P_site_valid
    ]
    stop_motif <- extract_motif_seqs_fast(stop_psite_dt, fasta_map, TARGET_SITE, MOTIF_RANGE, AA_MOTIF_CODON_NUM)
    ggsave(
      file.path(sample_motif_dir, "3_P_site_stop_region_motif.pdf"),
      plot = plot_motif_logo_fast(stop_motif$valid_dna_seqs, sample_name, TARGET_SITE, FALSE, TRUE),
      width=8, height=3, dpi=200, device="pdf", bg="white"
    )
    stop_aa_seqs <- dna2aa_seq_fast(stop_motif$valid_aa_seqs)
    ggsave(
      file.path(sample_motif_dir, "3_P_site_stop_region_aminomotif.pdf"),
      plot = plot_motif_logo_fast(stop_aa_seqs, sample_name, TARGET_SITE, TRUE, TRUE),
      width=8, height=3, dpi=200, device="pdf", bg="white"
    )
    
    if (!DISABLE_TOP3) {
      if ("P_codon" %in% colnames(psite_dt)) {
        codon_count <- psite_dt[!is.na(P_codon) & nchar(P_codon)==3, .N, by=P_codon][order(-N)]
        top3_codon <- head(codon_count$P_codon,3)
        
        for (i in 1:3) {
          if (i <= length(top3_codon) && !is.na(top3_codon[i])) {
            top_dt <- psite_dt[P_codon == top3_codon[i]]
            top_motif <- extract_motif_seqs_fast(top_dt, fasta_map, TARGET_SITE, MOTIF_RANGE, AA_MOTIF_CODON_NUM)
            ggsave(
              file.path(sample_motif_dir, paste0(i+3, "_P_site_TOP",i,"_codon_motif.pdf")),
              plot = plot_motif_logo_fast(top_motif$valid_dna_seqs, sample_name, TARGET_SITE, FALSE),
              width=8, height=3, dpi=200, device="pdf", bg="white"
            )
          } else {
            ggsave(
              file.path(sample_motif_dir, paste0(i+3, "_P_site_TOP",i,"_codon_motif.pdf")),
              plot = ggplot() +
                annotate("text", x=0.5, y=0.5, label="No valid sequences", size=4) +
                ggtitle(paste0(sample_name, " TOP",i," Codon")) +
                theme_void() +
                theme(plot.title = element_text(size=10, hjust=0.5)),
              width=8, height=3, dpi=200, device="pdf", bg="white"
            )
          }
        }
      }
    }
    
    cat(paste0("  >>> ", sample_name, " - Running Codon Meta-peaks (-30~+30nt) analysis...\n"))
    codon_meta_counts <- extract_codon_metapeak_counts(psite_dt, fasta_map)
    plot_codon_metapeak(codon_meta_counts, sample_name, sample_motif_dir)
    sample_motif_dir_asite <- file.path(MOTIF_OUTPUT_ROOT_ASITE, sample_name)
    dir.create(sample_motif_dir_asite, recursive = TRUE, showWarnings = FALSE)
    
    cat(paste0("  >>> ", sample_name, " - Running A-site Codon Meta-peaks...\n"))
    codon_meta_counts_asite <- extract_codon_metapeak_counts_asite(psite_dt, fasta_map)
    plot_codon_metapeak_asite(codon_meta_counts_asite, sample_name, sample_motif_dir_asite)
    rm(psite_dt, general_motif, start_motif, stop_motif); invisible(gc())
    return(paste(sample_name, "success"))
  }, error = function(e) {
    return(paste(sample_name, "error:", substr(e$message,1,100)))
  })
}

cat("============================================================\n")
cat("Ribosome Motif Analysis - inherit riboWaltz params + Codon Meta-peaks (Reads Density)\n")
cat("Start Time:", Sys.time(), "\n")
cat("Cores:", N_CORES, " Threads:", N_THREADS, " Memory:", MEMORY_LIMIT, "\n")
cat("Future Global Limit:", format(object.size(getOption("future.globals.maxSize")), units = "GB"), "\n")
cat("Meta-peaks Type: P-site reads density line plot (-30~+30nt)\n")
cat("============================================================\n\n")

cat("Step 1: Loading transcriptome (memory mapped)...\n")
fasta_map <- load_transcriptome_fast(TRANSCRIPTOME_FASTA)
assign("fasta_map", fasta_map, envir = .GlobalEnv)
cat("Transcriptome loaded (", nrow(fasta_map), "transcripts, size:", format(object.size(fasta_map), units = "MB"), ")\n\n")

cat("Step 2: Getting sample list...\n")
sample_dirs <- list.dirs(OUTPUT_ROOT, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[!grepl(paste0(MOTIF_OUTPUT_ROOT, "|\\.txt$|\\.pdf$|\\.rds$"), sample_dirs)]
sample_dirs <- sample_dirs[file.exists(file.path(sample_dirs, "reads_psite.rds"))]
cat("Found", length(sample_dirs), "valid samples\n\n")

cat("Step 3: Processing samples (", N_CORES, " threads)...\n")
sample_batches <- split(sample_dirs, ceiling(seq_along(sample_dirs)/BATCH_SIZE))

global_vars <- c(
  "fasta_map", "TARGET_SITE", "MOTIF_RANGE", "AA_MOTIF_CODON_NUM",
  "DISABLE_TOP3", "META_PEAKS_RANGE", "codon2aa", "codon2aa_hash",
  "read_ribo_rds_fast", "calculate_codon_pos_fast", "extract_motif_seqs_fast",
  "dna2aa_seq_fast", "plot_motif_logo_fast", "extract_codon_metapeak_counts",
  "plot_codon_metapeak",
  "save_metapeak_counts",
  "extract_codon_metapeak_counts_asite",
  "plot_codon_metapeak_asite"
)

for (batch_idx in seq_along(sample_batches)) {
  batch <- sample_batches[[batch_idx]]
  cat(paste0("Processing batch ", batch_idx, "/", length(sample_batches), " (", length(batch), " samples)...\n"))
  
  batch_results <- future_map(
    .x = batch,
    .f = ~process_sample_fast(.x),
    .options = furrr_options(
      seed = TRUE,
      globals = global_vars
    )
  )
  
  cat("Batch", batch_idx, "results:\n")
  print(unlist(batch_results))
  cat("\n")
  
  invisible(gc())
}

cat("============================================================\n")
cat("Analysis Completed! Time:", Sys.time(), "\n")
cat("============================================================\n")
cat("Summary:\n")
cat("  Output Directory: ", MOTIF_OUTPUT_ROOT, "\n")
cat("  Total Samples:    ", length(sample_dirs), "\n")
cat("  Cores Used:       ", N_CORES, "\n")
cat("  Threads:          ", N_THREADS, "\n")
cat("  Memory Limit:     ", MEMORY_LIMIT, "\n")
cat("  TOP3 Analysis:    ", ifelse(DISABLE_TOP3, "Disabled", "Enabled"), "\n")
cat("  Codon Meta-peaks: P-site reads density line plot (-30nt ~ +30nt)\n")
cat("\nAll results saved to", MOTIF_OUTPUT_ROOT, "!\n")

plan(sequential)
invisible(gc())
