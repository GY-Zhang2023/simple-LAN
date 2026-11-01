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
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = FALSE)
rm(list = ls())

options(datatable.auto.index = TRUE)
options(datatable.optimize = 2)

human_gtf <- "~/reference/human_transcriptome/gencode.v49.pc_annotation.gtf"
human_fasta <- "~/reference/human_transcriptome/gencode.v49.pc_transcripts.fa"
mouse_gtf <- "~/reference/mouse_transcriptome/gencode.vM35.pc_annotation.gtf"
mouse_fasta <- "~/reference/mouse_transcriptome/gencode.vM35.pc_transcripts.fa"

species_config <- list(
  human = list(
    gtf = human_gtf,
    fasta = human_fasta
  ),
  mouse = list(
    gtf = mouse_gtf,
    fasta = mouse_fasta
  )
)

selected_species <- ""
min_len_default <- NA
max_len_default <- NA
min_len <- NA
max_len <- NA
selected_analysis <- "monosome"
trim_nt <- 0

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 7) {
  selected_species <- args[1]
  min_len_default <- as.integer(args[2])
  max_len_default <- as.integer(args[3])
  min_len <- as.integer(args[4])
  max_len <- as.integer(args[5])
  selected_analysis <- ifelse(length(args)>=6 & tolower(args[6])=="disome", "disome", "monosome")
  trim_nt <- as.integer(args[7])
  cat(paste0("\n✅ Load run parameters from bash:\n"))
  cat(paste0("   Species: ", selected_species, "\n"))
  cat(paste0("   Base RPF length range: ", min_len_default, "-", max_len_default, "nt\n"))
  cat(paste0("   Strict RPF length range: ", min_len, "-", max_len, "nt\n"))
  cat(paste0("   Analysis pipeline: ", ifelse(selected_analysis=="monosome", "Standard Monosome", paste0("Disome - 5' trim ",trim_nt," nt\n"))))

  cat("\n🔍 Start parameter validation...\n")
  if (!selected_species %in% c("human", "mouse")) {
    stop(paste0("❌ Parameter error: species must be human / mouse, received: ", selected_species))
  }
  if (any(is.na(c(min_len_default, max_len_default, min_len, max_len)))) {
    stop("❌ Parameter error: RPF length parameters must be positive integers!")
  }
  if (max_len_default < min_len_default | max_len < min_len | min_len < min_len_default | max_len > max_len_default) {
    stop("❌ Parameter error: unreasonable RPF length range! Check base / strict length bounds")
  }
  if (selected_analysis == "disome") {
    if (is.na(trim_nt) | trim_nt <= 0) {
      stop(paste0("❌ Disome pipeline error: trim length must be positive integer, received: ", trim_nt))
    }
    cat(paste0("✅ Disome validation passed: 5' trim length = ", trim_nt, " nt\n"))
  } else {
    trim_nt <- 0
    cat("✅ Monosome validation passed: no read trimming\n")
  }
  cat("✅ All parameters validated, start analysis...\n")

} else if (length(args) >= 5) {
  selected_species <- args[1]
  min_len_default <- as.integer(args[2])
  max_len_default <- as.integer(args[3])
  min_len <- as.integer(args[4])
  max_len <- as.integer(args[5])
  selected_analysis <- ifelse(length(args)>=6 & tolower(args[6])=="disome", "disome", "monosome")
  trim_nt <- 33
  cat(paste0("\n✅ Load run parameters from bash:\n"))
  cat(paste0("   Species: ", selected_species, "\n"))
  cat(paste0("   Base RPF length range: ", min_len_default, "-", max_len_default, "nt\n"))
  cat(paste0("   Strict RPF length range: ", min_len, "-", max_len, "nt\n"))
  cat(paste0("   Analysis pipeline: ", ifelse(selected_analysis=="monosome", "Standard Monosome", "Disome - 5' trim 33 nt\n")))

} else if (length(args) == 1 && args[1] %in% names(species_config)) {
  selected_species <- args[1]
  selected_analysis <- "monosome"
  trim_nt <- 0
  min_len_default <- 21
  max_len_default <- 40
  min_len <- 26
  max_len <- 34
  cat(paste0("\n✅ Only species supplied, load default config:\n"))
  cat(paste0("   Species: ", selected_species, "\n"))
  cat(paste0("   Base RPF length range: 21-40 nt\n"))
  cat(paste0("   Strict RPF length range: 26-34 nt\n"))
  cat(paste0("   Analysis pipeline: Standard Monosome - no trimming\n"))

} else if (interactive()) {
  cat("=====================================\n")
  cat("          Select species              \n")
  cat("=====================================\n")
  cat("1) human\n")
  cat("2) mouse\n")

  repeat {
    selected <- readline(prompt = "Input selection (1/2): ")
    if (selected == "1") {
      selected_species <- "human"
      break
    } else if (selected == "2") {
      selected_species <- "mouse"
      break
    } else {
      cat("❌ Invalid selection! Only 1(human)/2(mouse) allowed\n")
    }
  }

  cat("\n=====================================\n")
  cat("          Select pipeline             \n")
  cat("=====================================\n")
  cat("1) Monosome standard pipeline (no trimming)\n")
  cat("2) Disome pipeline (trim 5' end then infer Rear Codon)\n")
  repeat {
    analysis_choice <- readline(prompt = "Input pipeline selection (1/2): ")
    if (analysis_choice == "1") {
      selected_analysis <- "monosome"
      trim_nt <- 0
      break
    } else if (analysis_choice == "2") {
      selected_analysis <- "disome"
      cat("\n⚠️ Disome required: set 5' trim length ⚠️\n")
      repeat {
        trim_nt_input <- readline(prompt = "Input 5' trim length (nt, recommend 33): ")
        if (trim_nt_input == "") {
          trim_nt <- 33
          break
        } else if (grepl("^\\d+$", trim_nt_input) && as.integer(trim_nt_input) > 0) {
          trim_nt <- as.integer(trim_nt_input)
          break
        } else {
          cat("❌ Please input positive integer!\n")
        }
      }
      break
    } else {
      cat("❌ Invalid selection! Only 1(Monosome)/2(Disome) allowed\n")
    }
  }

  if(selected_analysis == "disome"){
    confirm_text = paste0("【", selected_species, "】 + Disome(5' trim ",trim_nt," nt) + Rear Codon inference")
  }else{
    confirm_text = paste0("【", selected_species, "】 + Standard Monosome")
  }
  cat(paste0("\n⚠️ Ready to run: ", confirm_text, "\n"))
  repeat {
    confirm <- readline(prompt = "Confirm continue? (y/n): ")
    if (tolower(confirm) %in% c("y", "n")) {
      break
    }
    cat("❌ Only y/n accepted!\n")
  }
  if (tolower(confirm) != "y") {
    cat("\n🔴 User cancelled, exit\n")
    quit(save = "no", status = 0)
  }

  cat("\n=====================================\n")
  cat("        Set RPF length filters        \n")
  cat("        Refer to 3.alignments/Report.txt        \n")
  cat("=====================================\n")
  repeat {
    min_len_input <- readline(prompt = "Base filter minimum RPF length (default 21): ")
    if (min_len_input == "") {
      min_len_default <- 21
      break
    }
    if (grepl("^\\d+$", min_len_input)) {
      min_len_default <- as.integer(min_len_input)
      break
    } else {
      cat("❌ Input positive integer!\n")
    }
  }
  repeat {
    max_len_input <- readline(prompt = "Base filter maximum RPF length (default 40): ")
    if (max_len_input == "") {
      max_len_default <- 40
      break
    }
    if (grepl("^\\d+$", max_len_input)) {
      max_len_default <- as.integer(max_len_input)
      if (max_len_default >= min_len_default) {
        break
      } else {
        cat(paste0("❌ Max length cannot be smaller than min (", min_len_default, ")\n"))
      }
    } else {
      cat("❌ Input positive integer!\n")
    }
  }
  repeat {
    min_len_strict_input <- readline(prompt = "Strict filter minimum length (>1M RPF, default 26): ")
    if (min_len_strict_input == "") {
      min_len <- 26
      break
    }
    if (grepl("^\\d+$", min_len_strict_input)) {
      min_len <- as.integer(min_len_strict_input)
      if (min_len >= min_len_default && min_len <= max_len_default) {
        break
      } else {
        cat(paste0("❌ Strict min must fall within [", min_len_default, ",", max_len_default, "]\n"))
      }
    } else {
      cat("❌ Input positive integer!\n")
    }
  }
  repeat {
    max_len_strict_input <- readline(prompt = "Strict filter maximum length (>1M RPF, default 34): ")
    if (max_len_strict_input == "") {
      max_len <- 34
      break
    }
    if (grepl("^\\d+$", max_len_strict_input)) {
      max_len <- as.integer(max_len_strict_input)
      if (max_len >= min_len && max_len <= max_len_default) {
        break
      } else {
        cat(paste0("❌ Strict max must fall within [", min_len, ",", max_len_default, "]\n"))
      }
    } else {
      cat("❌ Input positive integer!\n")
    }
  }
  cat(paste0("\n⚠️ Length parameters confirmed:\n"))
  cat(paste0("   Base range: ", min_len_default, "-", max_len_default, "nt\n"))
  cat(paste0("   High-RPF strict range: ", min_len, "-", max_len, "nt\n"))
  repeat {
    len_confirm <- readline(prompt = "Confirm these lengths? (y/n): ")
    if (tolower(len_confirm) %in% c("y", "n")) {
      break
    }
    cat("❌ Only y/n accepted!\n")
  }
  if (tolower(len_confirm) != "y") {
    cat("\n🔴 User cancelled, exit\n")
    quit(save = "no", status = 0)
  }

} else {
  cat("\n❌ Error: species argument required in non-interactive mode!\n")
  cat("Usage1(monosome default): Rscript this_script.R human|mouse\n")
  cat("Usage2(length args): Rscript this_script.R human|mouse 21 40 26 34\n")
  cat("Usage3(disome): Rscript this_script.R human|mouse 21 40 26 34 disome 33\n")
  quit(save = "no", status = 1)
}

gtf_path <- species_config[[selected_species]]$gtf
transcriptome_fasta <- species_config[[selected_species]]$fasta

if(selected_analysis == "disome"){
  flow_text = paste0("Disome pipeline - 5' trim ",trim_nt,"nt + 5' Rear Codon inference + 3nt periodicity correction")
}else{
  flow_text = "Standard Monosome pipeline"
}
cat(paste0("\n✅ Loaded config for 【", selected_species, "】:\n"))
cat(paste0("   GTF: ", gtf_path, "\n"))
cat(paste0("   Transcriptome FASTA: ", transcriptome_fasta, "\n"))
cat(paste0("   Selected pipeline: ", flow_text, "\n"))

project_dir <- getwd()
output_dir <- file.path(project_dir, "4.riboWaltz_output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("✅ Created output directory: ", output_dir, "\n")
}

WORK_DIR <- getwd()
main_dir <- paste0(WORK_DIR, "/3.alignments")
output_root <- paste0(WORK_DIR, "/4.riboWaltz_output")
report_path <- file.path(output_root, "riboWaltz_batch_run_report.txt")

sink(report_path, split = TRUE)

cat("============================================================\n")
cat(paste0("riboWaltz Batch Analysis Report (A/P/E site separation) - ", flow_text))
cat("\nReport generated at: ", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

library(riboWaltz)
library(data.table)
library(GenomicFeatures)
library(ggplot2)
library(Biostrings)
library(GenomicAlignments)
library(IRanges)
library(txdbmaker)
library(GenomicRanges)

threads <- 8
rpf_count_threshold <- 1e6
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

cat("Environment and parameter settings\n")
cat("------------------------------------------------------------\n")
cat("- System language: English\n")
cat("- Max DLL count: 999\n")
cat("- stringsAsFactors: disabled\n")
cat("- Loaded packages: riboWaltz, data.table, GenomicFeatures, ggplot2, Biostrings, GenomicAlignments, IRanges, txdbmaker, GenomicRanges\n")
cat("- Pipeline: ", flow_text, "\n")

cat("Core parameters\n")
cat("------------------------------------------------------------\n")
cat("- BAM directory: ", main_dir, "\n")
cat("- GTF: ", gtf_path, "\n")
cat("- Transcriptome FASTA: ", transcriptome_fasta, "\n")
cat("- Output root: ", output_root, "\n")
cat("- Base RPF length filter: ", min_len_default, "-", max_len_default, "nt\n")
cat("- RPF threshold: ", rpf_count_threshold/1e6, "M\n")
cat("- Strict high-RPF length filter: ", min_len, "-", max_len, "nt\n")
cat("- Threads: ", threads, "\n")
cat("- Plot rule: avg_occupancy centered at Median; Normalized_* exported only, not plotted\n\n")

trans_bam_paths <- list.files(
  path = main_dir,
  pattern = "Aligned.toTranscriptome.out.bam$",
  full.names = TRUE,
  recursive = FALSE
)

cat("BAM file discovery\n")
cat("------------------------------------------------------------\n")
cat("Found ", length(trans_bam_paths), " target BAM files:\n")
print(trans_bam_paths)
cat("\n")

if (length(trans_bam_paths) == 0) {
  sink()
  stop("Aligned.toTranscriptome.out.bam not found, check directory path!")
}

anno_path <- file.path(output_root, "ribo_annotation.txt")
CLASSIFY_ANNOT <- file.path(output_root, "tx_annotation_for_classify.tsv")
cat("Annotation generation / loading\n")
cat("------------------------------------------------------------\n")

if (!file.exists(CLASSIFY_ANNOT)) {
  cat("🚀 Start building standard riboWaltz annotation...\n")
  tmp_r_script <- tempfile(fileext = ".R")
  r_script_content <- sprintf(
    "
assignInNamespace(
  'makeTxDbFromGFF',
  function(...) {
    txdbmaker::makeTxDbFromGFF(...)
  },
  ns = 'GenomicFeatures'
)

library(S4Vectors)
library(txdbmaker)
library(riboWaltz)
library(data.table)

gtf_path <- '%s'
anno_path <- '%s'
CLASSIFY_ANNOT <- '%s'

if (!file.exists(anno_path)) {
  cat('Start generating riboWaltz annotation file...\\n')
  stopifnot(file.exists(gtf_path))
  anno_dt <- create_annotation(gtf = gtf_path)
  fwrite(anno_dt, anno_path, sep='\\t', quote=FALSE)
} else {
  anno_dt <- fread(anno_path)
}

fwrite(anno_dt[, .(transcript, l_tr, l_utr5, l_cds, l_utr3)],
       CLASSIFY_ANNOT, sep='\\t', quote=FALSE)
cat('Annotation file finished: ', CLASSIFY_ANNOT, '\\n')
",
    gtf_path, anno_path, CLASSIFY_ANNOT
  )
  writeLines(r_script_content, con = tmp_r_script)
  sys_ret <- system2(command = "Rscript", args = tmp_r_script, stdout = TRUE, stderr = TRUE)
  cat("R annotation script log:\n", paste(sys_ret, collapse = "\n"), "\n")
  if (!file.exists(anno_path)) {
    stop(paste0("❌ Annotation generation failed, subprocess log above\nanno_path: ", anno_path, "\ntmp script: ", tmp_r_script))
  }
  unlink(tmp_r_script)
  anno_dt <- fread(anno_path)
}
if (!file.exists(CLASSIFY_ANNOT)) {
  anno_dt[, .(transcript, l_tr, l_utr5, l_cds, l_utr3)] %>% fwrite(CLASSIFY_ANNOT, sep = "\t")
}

anno_dt <- anno_dt[, .(transcript, l_utr5, l_cds, l_utr3)]
anno_dt <- anno_dt[!is.na(l_utr5) & !is.na(l_cds) & l_cds > 0 & !is.na(l_utr3),]
anno_dt[, `:=`(
  cds_start = l_utr5 + 1,
  cds_stop = l_utr5 + l_cds
)]

c_transcripts <- anno_dt[
  l_utr5 > 0 & l_cds > 0 & l_cds %% 3 == 0 & l_utr3 > 0,
  unique(transcript)
]
cat(paste0("✅ Valid transcript set for frame_psite: ", length(c_transcripts), " transcripts\n"))
cat(paste0("✅ Filter rule: 5'UTR>0 + CDS>0 + CDS length divisible by 3 + 3'UTR>0\n"))

bam_to_ribo_reads <- function(bam_path, anno_dt) {
  gal <- GenomicAlignments::readGAlignments(
    file = bam_path,
    use.names = TRUE,
    param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
  )

  reads_dt <- data.table(
    transcript = as.character(GenomicRanges::seqnames(gal)),
    end5 = GenomicRanges::start(gal),
    end3 = GenomicRanges::end(gal),
    length = GenomicRanges::width(gal)
  )

  if (selected_analysis == "disome") {
    reads_dt <- reads_dt[length >= trim_nt,]
    reads_dt[, `:=`(
      ori_end5 = end5,
      ori_length = length,
      trim_nt_used = trim_nt
    )]
    reads_dt[, end5 := end5 + trim_nt]
    reads_dt[, length := length - trim_nt]
    cat(paste0("   ⚡ Disome trim ",trim_nt," nt done\n"))
  } else {
    reads_dt[, `:=`(ori_end5 = end5, ori_length = length, trim_nt_used = 0)]
  }

  data.table::setkey(reads_dt, transcript)
  reads_dt <- merge(
    reads_dt,
    anno_dt[, .(transcript, cds_start, cds_stop, l_utr3)],
    by.x = "transcript",
    by.y = "transcript",
    all.x = TRUE
  )[, .(
    transcript,
    end5,end3,length,ori_end5,ori_length,trim_nt_used,cds_start,cds_stop,l_utr3
  )]

  reads_dt <- reads_dt[!is.na(cds_start) & !is.na(cds_stop) & cds_stop > cds_start]
  return(reads_dt)
}

ribo_unified_theme <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10.5, face = "bold", margin = margin(b=5)),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, margin = margin(b=8)),
      axis.title.x = element_text(size = 9, margin = margin(t=5)),
      axis.title.y = element_text(size = 9, margin = margin(r=5)),
      axis.text = element_text(size = 8, color = "black"),
      legend.position = "none",
      panel.border = element_rect(color = "black", linewidth = 0.3),
      panel.grid = element_blank()
    )
}

plot_length_group_periodicity <- function(reads_psite_dt, sample_name, sample_output_dir, current_min, current_max,trim_nt) {
  cat(paste0("   ✅ Start generating Frame distribution plot | Trimmed ",current_min,"-",current_max,"nt] | RAW Length | \n"))

  plot_dt <- reads_psite_dt
  plot_dt[, frame_val := (psite - cds_start) %% 3]
  plot_dt[, frame := factor(frame_val, levels = c(0,1,2), labels = c("Frame 0", "Frame 1", "Frame 2"))]
  plot_dt <- plot_dt[frame_val %in% c(0,1,2)]

  valid_frame_reads <- nrow(plot_dt)
  total_filtered_reads <- nrow(reads_psite_dt[!is.na(length) & length >= current_min & length <= current_max])
  frame_valid_ratio <- round(valid_frame_reads / total_filtered_reads * 100, 2)
  cat(paste0("   ✅ Frame validity filter: only CDS-resident complete codon P-sites counted, valid rate = ", frame_valid_ratio, "%\n"))

  len_stats <- plot_dt[, .(
    read_count = .N,
    trimed_length = unique(length)
  ), by = .(raw_length = ori_length, frame)]
  len_stats[, total_reads := sum(read_count), by = raw_length]
  len_stats[, frame_percent := round(read_count / total_reads * 100, 2)]
  len_stats[, lab_text := paste0(read_count, "\n(", frame_percent, "%)")]

  total_psite_reads <- nrow(reads_psite_dt)
  valid_cds_psite_reads <- nrow(plot_dt)
  invalid_cds_psite_reads <- total_psite_reads - valid_cds_psite_reads
  cds_psite_efficiency <- round(valid_cds_psite_reads / total_psite_reads * 100, 2)

  cat(paste0("📊 CDS P-site efficiency (trim ", trim_nt, "nt):\n"))
  cat(paste0("   - Total P-site corrected reads: ", total_psite_reads, "\n"))
  cat(paste0("   - Valid CDS P-site reads: ", valid_cds_psite_reads, "\n"))
  cat(paste0("   - Invalid out-of-CDS P-site reads: ", invalid_cds_psite_reads, "\n"))
  cat(paste0("   - CDS P-site efficiency: ", cds_psite_efficiency, "%\n"))

  merged_plot <- ggplot(len_stats, aes(x = as.factor(raw_length), y = read_count, fill = frame)) +
    geom_col(width = 0.7, alpha = 0.9, color = "black", size = 0.3, position = position_dodge(width = 0.75)) +
    geom_text(aes(label = lab_text), position = position_dodge(width = 0.75),
              vjust = -0.4, hjust = 0.5, size = 3.0, color = "black", fontface = "bold") +
    scale_fill_manual(values = c("Frame 0" = "#2E86AB", "Frame 1" = "#E74C3C", "Frame 2" = "#FECA57")) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = scales::comma) +
    labs(
      title = paste0(sample_name, " | Frame Distribution (X-axis: Raw RPF Length | Trimmed: ", current_min, "-", current_max, "nt)"),
      subtitle = paste0("5' Trim = ", trim_nt, " nt | Total Valid CDS Reads: ", format(valid_cds_psite_reads, big.mark = ","), " | CDS Efficiency: ", cds_psite_efficiency, "%"),
      x = "Raw RPF Length (nt)",
      y = "RPF Read Counts",
      fill = "Frame Phase |"
    ) +
    ribo_unified_theme() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#666666"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.key.width = unit(1.2, "cm"),
      panel.grid.major.x = element_blank()
    )

  ggsave(file.path(sample_output_dir, paste0("FRAME_DISTRIBUTION_RAW_LENGTH_3Bars_Per_Length_Trim",trim_nt,"nt.pdf")),
         plot = merged_plot, width = 18, height = 7, dpi = 300)
  cat(paste0("   ✅ Length-stratified frame bar plot saved\n"))

  len_stats_out <- len_stats[, .(raw_length, trimed_length, frame, read_count, frame_percent, total_reads)]
  fwrite(len_stats_out,
         file.path(sample_output_dir, paste0("All_Trim",trim_nt,"nt_Frame_Stats_RAW_LENGTH.txt")),
         sep = "\t", row.names = F)
  cat(paste0("   ✅ Frame statistics table exported\n"))

  ribowaltz_overall_3nt <- frame_psite(
    data = reads_psite,
    annotation = anno_dt,
    sample = sample_name,
    multisamples = "independent",
    plot_style = "split",
    region = "cds",
    length_range = current_min:current_max,
    colour = "#333f50"
  )

  overall_stats <- ribowaltz_overall_3nt$count_dt
  fwrite(overall_stats, file.path(sample_output_dir, paste0(sample_name, "_3nt_OVERALL_stats.txt")), sep="\t")

  overall_plot <- ggplot(overall_stats, aes(x = factor(frame), y = count, fill = factor(frame))) +
    geom_col(width = 0.6, alpha = 0.8, color = "black", size = 0.2) +
    scale_fill_manual(values = c("0" = "#333f50", "1" = "#333f50", "2" = "#333f50")) +
    labs(
      title = paste0(sample_name, " - CDS"),
      x = "Frame",
      y = "Count",
      fill = "Frame"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", size = 0.2)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

  ggsave(file.path(sample_output_dir, paste0(sample_name, "_riboWaltz_OFFICIAL_3nt_OVERALL.png")),
         plot = overall_plot, width =7, height=4.5, dpi=300)

  cat(paste0("   ✅ riboWaltz native 3nt periodic plot saved | Frame0=",overall_stats[frame==0]$count,", Frame1=",overall_stats[frame==1]$count,", Frame2=",overall_stats[frame==2]$count,"\n"))

  return(len_stats_out)
}

plot_global_trinucleotide_periodicity <- function(reads_psite_dt, sample_name, sample_output_dir) {
  cat(paste0("   ✅ Start generating GLOBAL trinucleotide periodicity plot\n"))

  global_dt <- reads_psite_dt
  global_dt[, frame_val := (psite - cds_start) %% 3]
  global_dt[, frame := factor(frame_val, levels = c(0,1,2), labels = c("Frame 0", "Frame 1", "Frame 2"))]

  global_stats <- global_dt[, .(count = .N), by = frame]
  global_stats[, total_rpf := sum(count)]
  global_stats[, percent := round(count / total_rpf * 100, 2)]
  global_stats[, lab_text := paste0(count, "\n(", percent, "%)")]

  global_trim_nt <- unique(global_dt$trim_nt_used)
  global_trimed_len_range <- paste0(range(global_dt$length), collapse = "-")
  global_raw_len_range <- paste0(range(global_dt$ori_length), collapse = "-")

  global_plot <- ggplot(global_stats, aes(x = frame, y = count, fill = frame)) +
    geom_col(width = 0.65, alpha = 0.8, color = "black", size = 0.3) +
    scale_fill_manual(values = c("Frame 0" = "#2E86AB", "Frame 1" = "#E74C3C", "Frame 2" = "#FECA57")) +
    geom_text(aes(label = lab_text), vjust = -0.3, size = 3.6, color = "black", fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = paste0(sample_name, " | Trinucleotide Periodicity (OVERALL)"),
      subtitle = paste0("5' Trim = ", global_trim_nt, " nt | Trimmed Length: ", global_trimed_len_range, "nt | Real Raw Length: ", global_raw_len_range, "nt | Total CDS-mapped Reads: ", format(global_stats$total_rpf[1], big.mark = ","), " reads"),
      x = "Frame",
      y = "RPF Read Counts",
      fill = "Frame Phase"
    ) +
    ribo_unified_theme()

  ggsave(file.path(sample_output_dir, "trinucleotide_periodicity_OVERAL.pdf"),
         plot = global_plot, width = 7, height = 4.5, dpi = 300)

  fwrite(global_stats, file.path(sample_output_dir, "trinucleotide_periodicity_OVERALL_Stats.txt"), sep="\t")
  cat(paste0("   ✅ Global periodicity plot and statistics saved\n"))
  return(global_stats)
}

error_log <- data.table(sample_name = character(), error_message = character())
length_choice_log <- data.table(sample_name = character(), total_rpf = numeric(), length_range = character(), note = character())

cat("Sample batch processing log\n")
cat("------------------------------------------------------------\n")

for (bam_path in trans_bam_paths) {
  sample_name <- sub("_genome_Aligned.toTranscriptome.out.bam$", "", basename(bam_path))
  cat("\n=======================================\n")
  cat("Start processing sample: ", sample_name, "\n")
  cat("=======================================\n")

  sample_output_dir <- file.path(output_root, sample_name)
  dir.create(sample_output_dir, showWarnings = FALSE, recursive = TRUE)

  tryCatch({
    cat("1. Reading BAM and convert to riboWaltz format...\n")
    reads_dt <- bam_to_ribo_reads(bam_path, anno_dt)
    total_rpf <- nrow(reads_dt)
    cat("➡️ Total filtered RPF reads: >", trim_nt," nt ",round(total_rpf/1e6, 2), "M\n")

    if (total_rpf >= rpf_count_threshold) {
      current_min <- min_len
      current_max <- max_len
      length_note <- paste0("High quality (≥1M RPF), strict filter (", current_min, "-", current_max, "nt)")
      cat("   ✅ RPF ≥1M, apply strict length filter\n")
    } else {
      current_min <- min_len_default
      current_max <- max_len_default
      length_note <- paste0("Low RPF count (<1M), wide filter (", current_min, "-", current_max, "nt)")
      cat("   ⚠️ RPF <1M, apply wide length filter\n")
    }
    length_choice_log <- rbind(length_choice_log, data.table(
      sample_name = sample_name,total_rpf = total_rpf,length_range = paste0(current_min,"-",current_max,"nt"),note = length_note
    ))

    reads_list <- setNames(list(reads_dt), sample_name)

    cat("2. Filter RPF by length...\n")
    reads_filtered <- riboWaltz::length_filter(data = reads_list,length_filter_mode = "custom",length_range = current_min:current_max)
    reads_filtered_dt <- reads_filtered[[sample_name]]
    pre_psite_reads <- nrow(reads_filtered_dt)
    cat(paste0("➡️ After length filter, before P-site correction: ", pre_psite_reads, "\n"))

    cat("3. Calibrate P-site offset per RPF length...\n")
    if (selected_analysis == "disome") {
      psite_offset <- riboWaltz::psite(
        data = reads_filtered,
        flanking = 6,
        start = FALSE,
        extremity = "auto",
        plot = TRUE
      )
    } else {
      psite_offset <- riboWaltz::psite(
        data = reads_filtered,
        flanking = 6,
        extremity = "auto",
        plot = TRUE
      )
    }
    offset_match <- psite_offset[sample == sample_name, .(length, corrected_offset_from_5)]
    offset_match <- unique(na.omit(offset_match))
    colnames(offset_match) <- c("rpf_len", "offset")
    reads_psite_dt <- copy(reads_filtered_dt)
    reads_psite_dt[, offset := offset_match$offset[match(length, offset_match$rpf_len)]]
    reads_psite_dt <- reads_psite_dt[!is.na(offset)]
    psite_corrected_reads <- nrow(reads_psite_dt)
    cat(paste0("➡️ After P-site offset assignment: ", psite_corrected_reads, "\n"))
    fwrite(psite_offset, file.path(sample_output_dir, "psite_offset.txt"), sep = "\t")

    auto_plot_files <- list.files(path = getwd(),pattern = paste0("psite_offset_", sample_name, ".*\\.pdf$|psite_plot_.*\\.pdf$"),full.names = TRUE)
    if (length(auto_plot_files) > 0) {
      file.copy(auto_plot_files, sample_output_dir, overwrite = TRUE)
      file.remove(auto_plot_files)
    }
    cat("P-site offset calibration finished!\n")

    cat("4. Extract P-site coordinates with calibrated offset...\n")
    reads_psite <- riboWaltz::psite_info(data = reads_filtered,offset = psite_offset)
    reads_psite_dt <- reads_psite[[sample_name]]
    total_psite_reads_v2 <- nrow(reads_psite_dt)
    cat(paste0("➡️ After psite_info, before annotation merge: ", total_psite_reads_v2, "\n"))

    cols2del <- c("cds_start","cds_stop","end5")[c("cds_start","cds_stop","end5") %in% colnames(reads_psite_dt)]
    if(length(cols2del)>0) reads_psite_dt[, (cols2del) := NULL]

    reads_dt_filtered <- reads_dt[length %in% current_min:current_max,]
    reads_dt_filtered_unique <- reads_dt_filtered[, .SD[1], by = .(transcript, end5, end3)]

    reads_psite_dt <- merge(reads_psite_dt,
                            reads_dt_filtered_unique[, .(transcript, end5, end3, ori_end5, ori_length, cds_start, cds_stop, l_utr3)],
                            all.x=TRUE, sort=FALSE)

    cat("Rows before merge: ", nrow(reads_psite_dt), "\n")
    cat("Rows after merge: ", nrow(reads_psite_dt), "\n")

    reads_psite_dt <- reads_psite_dt[!is.na(cds_start) & !is.na(cds_stop),]
    total_psite_reads_v2<-nrow(reads_psite_dt)
    cat(paste0("➡️ After merging raw trim metadata: ", total_psite_reads_v2, "\n"))
    if(!"ori_end5" %in% colnames(reads_psite_dt)) reads_psite_dt[, ori_end5 := end5 - trim_nt]
    if(!"ori_length" %in% colnames(reads_psite_dt)) reads_psite_dt[, ori_length := length + trim_nt]
    reads_psite_dt[, psite_from_start := psite - cds_start]
    saveRDS(reads_psite_dt, file.path(sample_output_dir, "reads_psite_raw_version.rds"))

    cat("5. Generate standard QC plots + global periodicity + length-specific periodicity analysis...\n")

    cat("   ✅ Plot P-site distribution around start codon (ATG ± 30/60 nt)...\n")
    start_codon_data <- reads_psite_dt[
      psite_from_start %in% seq(-30, 60) & !is.na(psite_from_start)
    ][
      , .(mean_count = .N / length(unique(transcript))),
      by = .(x = psite_from_start)
    ][order(x)]
    start_codon_psite_lineplot <- ggplot(start_codon_data, aes(x = x, y = mean_count)) +
      geom_line(color = "#e84a5f", linewidth = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#333f50", linewidth=1) +
      labs(x = "relative to ATG position", y = "average P-site count",
           title = paste(sample_name, "RAW P-site distribution relative to ATG"),
           subtitle = paste0("Length filter: ", current_min, "-", current_max, "nt", ifelse(selected_analysis=="disome", paste0(" | 5' trim ",trim_nt,"nt"), ""))) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
            plot.subtitle = element_text(hjust=0.5, size=10),
            axis.title = element_text(size=11),
            axis.text = element_text(size=9))
    ggsave(
      file.path(sample_output_dir, "start_codon_psite_lineplot.pdf"),
      plot = start_codon_psite_lineplot,
      width = 12, height = 6,
      device = "pdf", dpi = 300
    )

    psite_cds_dt <- reads_psite_dt
    psite_cds_dt[, l_cds := cds_stop - cds_start + 1]
    psite_cds_dt[, cds_relative_pos := (psite - cds_start + 1) / l_cds * 100]
    bin_size <- 0.5
    cds_density_dt <- psite_cds_dt[
      !is.na(cds_relative_pos)
    ][
      , .(psite_count = .N),
      by = .(cds_bin = round(cds_relative_pos / bin_size) * bin_size)
    ][
      order(cds_bin)
    ][
      cds_bin >= -25 & cds_bin <= 125
    ]
    global_avg <- mean(cds_density_dt$psite_count)
    cds_density_dt[, relative_density := psite_count / global_avg]

    output_files <- file.path(sample_output_dir, "FootprintDensity.txt")
    fwrite(
      x = cds_density_dt,
      file = output_files,
      sep = "\t",
      na = "",
      quote = FALSE
    )

    final_plot <- ggplot(cds_density_dt, aes(x = cds_bin, y = relative_density)) +
      geom_line(color = "#2E8B01", linewidth = 0.8, alpha = 0.9) +
      geom_vline(xintercept = 0, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
      geom_vline(xintercept = 100, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
      annotate("text", x = 2, y = max(cds_density_dt$relative_density)*0.9,
               label = "Start Codon (ATG)", color = "#333f50", size = 3.5, hjust = 0) +
      annotate("text", x = 98, y = max(cds_density_dt$relative_density)*0.9,
               label = "Stop Codon", color = "#333f50", size = 3.5, hjust = 1) +
      labs(
        x = "Normalized CDS Length (%)",
        y = "Relative Footprint Density (a.u.)",
        title = paste(sample_name, "Ribosome Footprint Density (Raw Periodicity)"),
        subtitle = paste0("Binsize = ", bin_size, "% | Length filter: ", current_min, "-", current_max, "nt", ifelse(selected_analysis=="disome", paste0(" | 5' trim ",trim_nt,"nt (Disome)"), ""))
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#666666"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(alpha = 0.3)
      ) +
      ylim(0, max(cds_density_dt$relative_density) * 1.1)

    ggsave(
      file.path(sample_output_dir, "FootprintDensity.pdf"),
      plot = final_plot,
      width = 12, height = 6,
      device = "pdf", dpi = 300
    )
    cat("   ✅ Plot footprint density along normalized CDS\n")

    psite_cds_dt[, start_nt_pos := psite - cds_start]
    start_codon_dt <- psite_cds_dt[
      !is.na(start_nt_pos) & start_nt_pos >= -25 & start_nt_pos <= 125
    ][
      , .(psite_count = .N), by = start_nt_pos
    ][
      order(start_nt_pos)
    ][, relative_density := psite_count / global_avg]

    start_plot <- ggplot(start_codon_dt, aes(x = start_nt_pos, y = relative_density)) +
      geom_line(color = "#2E8B01", linewidth = 1, alpha = 0.9) +
      geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
      annotate("text", x = 5, y = max(start_codon_dt$relative_density)*0.9,
               label = "ATG (Start Codon)", color = "#e84a5f", size = 4, hjust = 0, fontface="bold") +
      labs(
        x = "Position relative to Start Codon (nt)",
        y = "Relative Footprint Density (a.u.)",
        title = paste(sample_name, "Footprint Density (Start Codon -25 nt to +125 nt)"),
        subtitle = paste0("No binning | Length filter: ", current_min, "-", current_max, "nt", ifelse(selected_analysis=="disome", paste0(" | 5' trim ",trim_nt,"nt"), ""))
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#666666"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(alpha = 0.3)
      ) +
      ylim(0, max(start_codon_dt$relative_density) * 1.1)

    ggsave(
      file.path(sample_output_dir, "FootprintDensity_StartCodon_Minus25_Plus125.pdf"),
      plot = start_plot, width = 10, height = 5, device = "pdf", dpi = 300
    )
    cat("   ✅ Plot start codon local density (-25 ~ +125 nt, no binning)\n")

    psite_cds_dt[, stop_nt_pos := psite - cds_stop]
    stop_codon_dt <- psite_cds_dt[
      !is.na(stop_nt_pos) & stop_nt_pos >= -125 & stop_nt_pos <= 25
    ][
      , .(psite_count = .N), by = stop_nt_pos
    ][
      order(stop_nt_pos)
    ][, relative_density := psite_count / global_avg]

    stop_plot <- ggplot(stop_codon_dt, aes(x = stop_nt_pos, y = relative_density)) +
      geom_line(color = "#2E8B01", linewidth = 1, alpha = 0.9) +
      geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
      annotate("text", x = -5, y = max(stop_codon_dt$relative_density)*0.9,
               label = "Stop Codon", color = "#e84a5f", size = 4, hjust = 1, fontface="bold") +
      labs(
        x = "Position relative to Stop Codon (nt)",
        y = "Relative Footprint Density (a.u.)",
        title = paste(sample_name, "Footprint Density (Stop Codon -125 nt to +25 nt)"),
        subtitle = paste0("No binning | Length filter: ", current_min, "-", current_max, "nt", ifelse(selected_analysis=="disome", paste0(" | 5' trim ",trim_nt,"nt"), ""))
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#666666"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(alpha = 0.3)
      ) +
      ylim(0, max(stop_codon_dt$relative_density) * 1.1)

    ggsave(
      file.path(sample_output_dir, "FootprintDensity_StopCodon_Minus125_Plus25.pdf"),
      plot = stop_plot, width = 10, height = 5, device = "pdf", dpi = 300
    )
    cat("   ✅ Plot stop codon local density (-125 ~ +25 nt, no binning)\n")

    cat("\n🔴 Apply global valid P-site filter (riboWaltz frame_psite rule)\n")
    region = "cds"
    length_range = current_min:current_max

    reads_psite_dt_before_filter = nrow(reads_psite_dt)
    reads_psite_dt <- reads_psite_dt[
      transcript %in% c_transcripts &
        length %in% length_range &
        !is.na(psite_region) & psite_region == region &
        !is.na(psite_from_start) &
        psite >= cds_start & psite <= (cds_stop - 3)
    ]
    reads_psite_dt_after_filter = nrow(reads_psite_dt)
    filter_remove_count = reads_psite_dt_before_filter - reads_psite_dt_after_filter
    filter_valid_ratio = round(reads_psite_dt_after_filter / reads_psite_dt_before_filter * 100
    cat(paste0("✅ Global valid P-site filter | Before: ", reads_psite_dt_before_filter, " | After: ", reads_psite_dt_after_filter, "\n"))
    cat(paste0("✅ Removed invalid P-sites: ", filter_remove_count, " | Valid rate: ", filter_valid_ratio, "%\n"))
    cat(paste0("✅ All downstream analyses will use this valid P-site set\n"))

    plot_global_trinucleotide_periodicity(reads_psite_dt, sample_name, sample_output_dir)
    plot_length_group_periodicity(reads_psite_dt, sample_name, sample_output_dir, current_min, current_max, trim_nt)

    cat("\n🔴 Apply Frame 0 strict filter | downstream codon / APE analyses use Frame 0 only\n")
    reads_psite_dt[, frame_val := (psite - cds_start) %% 3]

    reads_psite_TE <- setNames(list(reads_psite_dt), sample_name)
    saveRDS(reads_psite_TE, file.path(sample_output_dir, "reads_psite_TE_version.rds"))
    cat(paste0("✅ TE calculation RDS saved: reads_psite_TE_version.rds (all frames retained)\n"))

    frame0_before <- nrow(reads_psite_dt)
    reads_psite_dt <- reads_psite_dt[frame_val == 0]
    frame0_after <- nrow(reads_psite_dt)
    frame0_ratio <- round(frame0_after / frame0_before * 100, 2)
    cat(paste0("✅ Frame 0 filter | Before: ", frame0_before, " | After: ", frame0_after, "\n"))
    cat(paste0("✅ Frame 0 proportion: ", frame0_ratio, "% | Codon / APE / occupancy based on Frame 0 subset\n"))

    reads_psite_dt[, `:=`(
      rear_seq_start = NA_integer_, rear_seq_end = NA_integer_, rear_seq_length = NA_integer_,
      Rear_first_codon_pos = NA_integer_, Rear_first_cds_pos = NA_integer_, Rear_first_codon_num = NA_integer_,
      Rear_seq = NA_character_, Rear_Codon = NA_character_, Rear_Codon_next = NA_character_, Rear_Codon_prev = NA_character_,
      Rear_Codon_valid = FALSE, Rear_A_site = NA_integer_, Rear_P_site = NA_integer_, Rear_E_site = NA_integer_,
      Rear_P_cds_pos = NA_integer_, Rear_A_cds_pos = NA_integer_, Rear_E_cds_pos = NA_integer_,
      Rear_P_phase = NA_integer_, Rear_P_codon = NA_character_, Rear_A_codon = NA_character_, Rear_E_codon = NA_character_,
      Rear_APE_valid = FALSE
    )]
    setDT(reads_psite_dt)
    setkey(reads_psite_dt, transcript)

    cat("Calculating core CDS position for P/A/E sites...\n")
    reads_psite_dt[, P_cds_pos := ifelse(psite >= cds_start & psite <= (cds_stop - 3), psite - cds_start + 1, NA_integer_)]
    reads_psite_dt[, A_cds_pos := NA_integer_]
    reads_psite_dt[, E_cds_pos := NA_integer_]

    transcriptome <- Biostrings::readDNAStringSet(transcriptome_fasta)
    fasta_full_names <- names(transcriptome)
    fasta_enst_ids <- sub("\\|.*$", "", fasta_full_names)
    fasta_seqs <- as.character(transcriptome)

    cat("✅ FASTA ID processing completed!\n")
    cat("Original FASTA ID example: ", head(fasta_full_names, 1), "\n")
    cat("Extracted ENST ID example: ", head(fasta_enst_ids, 1), "\n")
    cat("BAM transcript ID example: ", head(reads_psite_dt$transcript, 1), "\n")

    cat("Calculating A/E site coordinates based on P-site...\n")
    reads_psite_dt[, `:=`(
      A_site = psite + 3,
      E_site = psite - 3
    )]
    reads_psite_dt[, `:=`(
      A_cds_pos = ifelse(A_site >= (cds_start + 3) & A_site <= cds_stop, P_cds_pos + 3, NA_integer_),
      E_cds_pos = ifelse(E_site >= (cds_start - 3) & E_site <= (cds_stop - 6), P_cds_pos - 3, NA_integer_)
    )]
    reads_psite_dt[, `:=`(
      P_codon_pos = ifelse(!is.na(P_cds_pos), floor((P_cds_pos - 1) / 3) + 1, NA_integer_),
      A_codon_pos = ifelse(!is.na(A_cds_pos), floor((A_cds_pos - 1) / 3) + 1, NA_integer_),
      E_codon_pos = ifelse(!is.na(E_cds_pos), floor((E_cds_pos - 1) / 3) + 1, NA_integer_)
    )]
    reads_psite_dt[, `:=`(A_site_valid = A_site, E_site_valid = E_site, P_site_valid = psite)]
    setDT(reads_psite_dt)

    cat("Extracting codon sequences for A/P/E sites...\n")
    reads_psite_dt[, fasta_idx := match(transcript, fasta_enst_ids)]
    reads_psite_dt <- reads_psite_dt[!is.na(fasta_idx)]
    reads_psite_dt[, seq_str := fasta_seqs[fasta_idx]]

    reads_psite_dt[, P_codon := ifelse(P_site_valid & !is.na(seq_str) & psite >= 1 & (psite + 2) <= nchar(seq_str), substr(seq_str, psite, psite + 2), NA_character_)]
    reads_psite_dt[, A_codon := ifelse(A_site_valid & !is.na(seq_str) & A_site >= 1 & (A_site + 2) <= nchar(seq_str), substr(seq_str, A_site, A_site + 2), NA_character_)]
    reads_psite_dt[, E_codon := ifelse(E_site_valid & !is.na(seq_str) & E_site >= 1 & (E_site + 2) <= nchar(seq_str), substr(seq_str, E_site, E_site + 2), NA_character_)]

    valid_p_ratio <- sum(!is.na(reads_psite_dt$P_codon)) / nrow(reads_psite_dt) * 100
    valid_a_ratio <- sum(!is.na(reads_psite_dt$A_codon)) / nrow(reads_psite_dt) * 100
    valid_e_ratio <- sum(!is.na(reads_psite_dt$E_codon)) / nrow(reads_psite_dt) * 100
    cat("✅ Codon sequence extraction completed!\n")
    cat("   Valid P-site codon ratio: ", sprintf("%.2f", valid_p_ratio), "%\n")
    cat("   Valid A-site codon ratio: ", sprintf("%.2f", valid_a_ratio), "%\n")
    cat("   Valid E-site codon ratio: ", sprintf("%.2f", valid_e_ratio), "%\n")

    cat("\n🔴 Start codon-specific density analysis for all valid Frame 0 codons\n")
    codon_plot_dir <- file.path(sample_output_dir, "codon_specific_analysis")
    dir.create(codon_plot_dir, showWarnings = FALSE, recursive = TRUE)

    if (!"l_cds" %in% colnames(reads_psite_dt)) reads_psite_dt[, l_cds := cds_stop - cds_start + 1]
    if (!"cds_relative_pos" %in% colnames(reads_psite_dt)) reads_psite_dt[, cds_relative_pos := (psite - cds_start + 1) / l_cds * 100]
    if (!"start_nt_pos" %in% colnames(reads_psite_dt)) reads_psite_dt[, start_nt_pos := psite - cds_start]
    if (!"stop_nt_pos" %in% colnames(reads_psite_dt)) reads_psite_dt[, stop_nt_pos := psite - cds_stop]

    standard_codons <- c(
      "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
      "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
      "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
      "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
      "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
      "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
      "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
      "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
    )

    library(RColorBrewer)
    color_palettes <- c(
      brewer.pal(12, "Paired"),
      brewer.pal(12, "Set3"),
      brewer.pal(11, "Spectral"),
      brewer.pal(8, "Dark2"),
      brewer.pal(8, "Accent"),
      brewer.pal(9, "Set1"),
      brewer.pal(4, "Pastel1")
    )
    fixed_codon_colors <- color_palettes[1:length(standard_codons)]
    names(fixed_codon_colors) <- standard_codons

    valid_codons <- unique(reads_psite_dt[!is.na(P_codon) & nchar(P_codon) == 3]$P_codon)
    valid_codons <- valid_codons[sapply(valid_codons, function(c) {
      nrow(reads_psite_dt[P_codon == c & !is.na(cds_relative_pos)]) > 0
    })]

    codon_color <- fixed_codon_colors[valid_codons]
    if (length(codon_color) < length(valid_codons)) {
      missing_codons <- setdiff(valid_codons, names(codon_color))
      missing_colors <- sample(colors(distinct = TRUE), length(missing_codons))
      codon_color <- c(codon_color, setNames(missing_colors, missing_codons))
    }
    cat(paste0("✅ Detected ", length(valid_codons), " valid Frame 0 codons with data\n"))

    global_cds_dt <- data.table(cds_bin = numeric(), psite_count = integer(), relative_density = numeric(), codon = character())
    global_start_dt <- data.table(start_nt_pos = numeric(), psite_count = integer(), relative_density = numeric(), codon = character())
    global_stop_dt <- data.table(stop_nt_pos = numeric(), psite_count = integer(), relative_density = numeric(), codon = character())

    for (target_codon in valid_codons) {
      codon_specific_dt <- reads_psite_dt[P_codon == target_codon & !is.na(cds_relative_pos)]
      if (nrow(codon_specific_dt) == 0) {
        cat(paste0("⚠️ Codon ", target_codon, " has no data, skip\n"))
        next
      }

      codon_cds_dt <- codon_specific_dt[, cds_bin := round(cds_relative_pos / bin_size) * bin_size][
        cds_bin >= -25 & cds_bin <= 125][
          , .(psite_count = .N), by = cds_bin][order(cds_bin)]

      if (nrow(codon_cds_dt) > 0) {
        codon_avg <- mean(codon_cds_dt$psite_count)
        codon_cds_dt[, relative_density := psite_count / codon_avg]
        codon_cds_dt[, codon := target_codon]
        if ("relative_density" %in% colnames(codon_cds_dt)) global_cds_dt <- rbind(global_cds_dt, codon_cds_dt, fill = TRUE)

        cds_plot <- ggplot(codon_cds_dt, aes(x = cds_bin, y = relative_density)) +
          geom_line(color = codon_color[target_codon], linewidth = 0.8, alpha = 0.9) +
          geom_vline(xintercept = 0, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
          geom_vline(xintercept = 100, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
          annotate("text", x = 2, y = max(codon_cds_dt$relative_density) * 0.9, label = "Start Codon (ATG)", color = "#333f50", size = 3.5, hjust = 0) +
          annotate("text", x = 98, y = max(codon_cds_dt$relative_density) * 0.9, label = "Stop Codon", color = "#333f50", size = 3.5, hjust = 1) +
          labs(x = "Normalized CDS Length (%)", y = "Relative Footprint Density (a.u.)",
               title = paste0(sample_name, " - Codon: ", target_codon, " (Full CDS | Frame 0)"),
               subtitle = paste0("Length filter: ", current_min, "-", current_max, "nt | ", flow_text)) +
          ribo_unified_theme() +
          theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 9, color = "#666666")) +
          ylim(0, max(codon_cds_dt$relative_density) * 1.1)

        ggsave(file.path(codon_plot_dir, paste0("Codon_", target_codon, "_FullCDS_density.pdf")),
               plot = cds_plot, width = 10, height = 5, dpi = 300)
      } else {
        codon_cds_dt <- data.table(cds_bin = numeric(), psite_count = integer(), relative_density = numeric())
      }

      codon_start_dt <- codon_specific_dt[!is.na(start_nt_pos) & start_nt_pos >= -25 & start_nt_pos <= 125][
        , .(psite_count = .N), by = start_nt_pos][order(start_nt_pos)]

      if (nrow(codon_start_dt) > 0) {
        start_avg <- mean(codon_start_dt$psite_count)
        codon_start_dt[, relative_density := psite_count / start_avg]
        codon_start_dt[, codon := target_codon]
        if ("relative_density" %in% colnames(codon_start_dt)) global_start_dt <- rbind(global_start_dt, codon_start_dt, fill = TRUE)

        start_plot <- ggplot(codon_start_dt, aes(x = start_nt_pos, y = relative_density)) +
          geom_line(color = codon_color[target_codon], linewidth = 1, alpha = 0.9) +
          geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
          annotate("text", x = 5, y = max(codon_start_dt$relative_density) * 0.9, label = "ATG (Start Codon)", color = "#e84a5f", size = 4, hjust = 0, fontface = "bold") +
          labs(x = "Position relative to Start Codon (nt)", y = "Relative Footprint Density (a.u.)",
               title = paste0(sample_name, " - Codon: ", target_codon, " (Start Codon ± nt | Frame 0)"),
               subtitle = paste0("Range: -25 nt to +125 nt | Length filter: ", current_min, "-", current_max, "nt")) +
          ribo_unified_theme() +
          theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 9, color = "#666666")) +
          ylim(0, max(codon_start_dt$relative_density) * 1.1)

        ggsave(file.path(codon_plot_dir, paste0("Codon_", target_codon, "_StartCodon_density.pdf")),
               plot = start_plot, width = 10, height = 5, dpi = 300)
      } else {
        codon_start_dt <- data.table(start_nt_pos = numeric(), psite_count = integer(), relative_density = numeric())
      }

      codon_stop_dt <- codon_specific_dt[!is.na(stop_nt_pos) & stop_nt_pos >= -125 & stop_nt_pos <= 25][
        , .(psite_count = .N), by = stop_nt_pos][order(stop_nt_pos)]

      if (nrow(codon_stop_dt) > 0) {
        stop_avg <- mean(codon_stop_dt$psite_count)
        codon_stop_dt[, relative_density := psite_count / stop_avg]
        codon_stop_dt[, codon := target_codon]
        if ("relative_density" %in% colnames(codon_stop_dt)) global_stop_dt <- rbind(global_stop_dt, codon_stop_dt, fill = TRUE)

        stop_plot <- ggplot(codon_stop_dt, aes(x = stop_nt_pos, y = relative_density)) +
          geom_line(color = codon_color[target_codon], linewidth = 1, alpha = 0.9) +
          geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
          annotate("text", x = -5, y = max(codon_stop_dt$relative_density) * 0.9, label = "Stop Codon", color = "#e84a5f", size = 4, hjust = 1, fontface = "bold") +
          labs(x = "Position relative to Stop Codon (nt)", y = "Relative Footprint Density (a.u.)",
               title = paste0(sample_name, " - Codon: ", target_codon, " (Stop Codon ± nt | Frame 0)"),
               subtitle = paste0("Range: -125 nt to +25 nt | Length filter: ", current_min, "-", current_max, "nt")) +
          ribo_unified_theme() +
          theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 9, color = "#666666")) +
          ylim(0, max(codon_stop_dt$relative_density) * 1.1)

        ggsave(file.path(codon_plot_dir, paste0("Codon_", target_codon, "_StopCodon_density.pdf")),
               plot = stop_plot, width = 10, height = 5, dpi = 300)
      } else {
        codon_stop_dt <- data.table(stop_nt_pos = numeric(), psite_count = integer(), relative_density = numeric())
      }

      codon_all_stats <- rbind(
        codon_cds_dt[, .(position_type = "Normalized_CDS", position = cds_bin, psite_count, relative_density)],
        codon_start_dt[, .(position_type = "Start_Codon_nt", position = start_nt_pos, psite_count, relative_density)],
        codon_stop_dt[, .(position_type = "Stop_Codon_nt", position = stop_nt_pos, psite_count, relative_density)]
      )
      if (nrow(codon_all_stats) > 0) {
        fwrite(codon_all_stats, file.path(codon_plot_dir, paste0("Codon_", target_codon, "_all_stats.txt")), sep = "\t", na = "")
      }
      cat(paste0("✅ Codon ", target_codon, " done (", nrow(codon_specific_dt), " reads)\n"))
    }

    if (nrow(global_cds_dt) > 0 | nrow(global_start_dt) > 0 | nrow(global_stop_dt) > 0) {
      cat("\n🔴 Generate merge plots for all codons\n")

      get_peak_pos <- function(dt, x_col, y_col, codon_col) {
        if (nrow(dt) == 0 || !all(c(x_col, y_col, codon_col) %in% colnames(dt))) {
          return(data.table(codon = character(), x_peak = numeric(), y_peak = numeric()))
        }
        peak_dt <- dt[, .SD[which.max(get(y_col))], by = get(codon_col)]
        colnames(peak_dt)[colnames(peak_dt) == paste0("get(", codon_col, ")")] <- codon_col
        return(peak_dt[, .(codon = get(codon_col), x_peak = get(x_col), y_peak = get(y_col))])
      }

      if (nrow(global_cds_dt) > 0 && "relative_density" %in% colnames(global_cds_dt)) {
        cds_peak_dt <- get_peak_pos(global_cds_dt, "cds_bin", "relative_density", "codon")
        merge_cds_plot <- ggplot(global_cds_dt, aes(x = cds_bin, y = relative_density, color = codon, group = codon)) +
          geom_line(linewidth = 0.6, alpha = 0.8) +
          {if (nrow(cds_peak_dt) > 0) geom_text(data = cds_peak_dt, aes(x = x_peak, y = y_peak * 1.05, label = codon, color = codon), size = 3.5, fontface = "bold", alpha = 0.9, hjust = 0.5)} +
          geom_vline(xintercept = 0, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
          geom_vline(xintercept = 100, color = "#333f50", linetype = "dashed", linewidth = 1, alpha = 0.7) +
          annotate("text", x = 2, y = max(global_cds_dt$relative_density) * 0.9, label = "Start Codon (ATG)", color = "#333f50", size = 3.5, hjust = 0) +
          annotate("text", x = 98, y = max(global_cds_dt$relative_density) * 0.9, label = "Stop Codon", color = "#333f50", size = 3.5, hjust = 1) +
          scale_color_manual(values = fixed_codon_colors) +
          labs(x = "Normalized CDS Length (%)", y = "Relative Footprint Density (a.u.)",
               title = paste0(sample_name, " - All Codons Merge (Full CDS | Frame 0)"),
               subtitle = paste0("Total valid codons: ", length(valid_codons), " | Length filter: ", current_min, "-", current_max, "nt"),
               color = "Codon") +
          ribo_unified_theme() +
          theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#666666"),
                legend.position = "right", legend.key.width = unit(1, "cm"), legend.text = element_text(size = 8)) +
          ylim(0, max(global_cds_dt$relative_density) * 1.2)

        ggsave(file.path(codon_plot_dir, paste0("Merge_AllCodons_FullCDS_density.pdf")),
               plot = merge_cds_plot, width = 14, height = 7, dpi = 300)
      }

      if (nrow(global_start_dt) > 0 && "relative_density" %in% colnames(global_start_dt)) {
        start_peak_dt <- get_peak_pos(global_start_dt, "start_nt_pos", "relative_density", "codon")
        merge_start_plot <- ggplot(global_start_dt, aes(x = start_nt_pos, y = relative_density, color = codon, group = codon)) +
          geom_line(linewidth = 0.7, alpha = 0.8) +
          {if (nrow(start_peak_dt) > 0) geom_text(data = start_peak_dt, aes(x = x_peak, y = y_peak * 1.05, label = codon, color = codon), size = 3.5, fontface = "bold", alpha = 0.9, hjust = 0.5)} +
          geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
          annotate("text", x = 5, y = max(global_start_dt$relative_density) * 0.9, label = "ATG (Start Codon)", color = "#e84a5f", size = 4, hjust = 0, fontface = "bold") +
          scale_color_manual(values = fixed_codon_colors) +
          labs(x = "Position relative to Start Codon (nt)", y = "Relative Footprint Density (a.u.)",
               title = paste0(sample_name, " - All Codons Merge (Start Codon ± nt | Frame 0)"),
               subtitle = paste0("Total valid codons: ", length(valid_codons), " | Range: -25 nt to +125 nt"),
               color = "Codon") +
          ribo_unified_theme() +
          theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#666666"),
                legend.position = "right", legend.key.width = unit(1, "cm"), legend.text = element_text(size = 8)) +
          ylim(0, max(global_start_dt$relative_density) * 1.2)

        ggsave(file.path(codon_plot_dir, paste0("Merge_AllCodons_StartCodon_density.pdf")),
               plot = merge_start_plot, width = 14, height = 7, dpi = 300)
      }

      if (nrow(global_stop_dt) > 0 && "relative_density" %in% colnames(global_stop_dt)) {
        stop_peak_dt <- get_peak_pos(global_stop_dt, "stop_nt_pos", "relative_density", "codon")
        merge_stop_plot <- ggplot(global_stop_dt, aes(x = stop_nt_pos, y = relative_density, color = codon, group = codon)) +
          geom_line(linewidth = 0.7, alpha = 0.8) +
          {if (nrow(stop_peak_dt) > 0) geom_text(data = stop_peak_dt, aes(x = x_peak, y = y_peak * 1.05, label = codon, color = codon), size = 3.5, fontface = "bold", alpha = 0.9, hjust = 0.5)} +
          geom_vline(xintercept = 0, color = "#e84a5f", linetype = "solid", linewidth = 1.2, alpha = 0.9) +
          annotate("text", x = -5, y = max(global_stop_dt$relative_density) * 0.9, label = "Stop Codon", color = "#e84a5f", size = 4, hjust = 1, fontface = "bold") +
          scale_color_manual(values = fixed_codon_colors) +
          labs(x = "Position relative to Stop Codon (nt)", y = "Relative Footprint Density (a.u.)",
               title = paste0(sample_name, " - All Codons Merge (Stop Codon ± nt | Frame 0)"),
               subtitle = paste0("Total valid codons: ", length(valid_codons), " | Range: -125 nt to +25 nt"),
               color = "Codon") +
          ribo_unified_theme() +
          theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#666666"),
                legend.position = "right", legend.key.width = unit(1, "cm"), legend.text = element_text(size = 8)) +
          ylim(0, max(global_stop_dt$relative_density) * 1.2)

        ggsave(file.path(codon_plot_dir, paste0("Merge_AllCodons_StopCodon_density.pdf")),
               plot = merge_stop_plot, width = 14, height = 7, dpi = 300)
      }

      global_all_stats <- rbind(
        if (nrow(global_cds_dt) > 0) global_cds_dt[, .(position_type = "Normalized_CDS", position = cds_bin, psite_count, relative_density, codon)],
        if (nrow(global_start_dt) > 0) global_start_dt[, .(position_type = "Start_Codon_nt", position = start_nt_pos, psite_count, relative_density, codon)],
        if (nrow(global_stop_dt) > 0) global_stop_dt[, .(position_type = "Stop_Codon_nt", position = stop_nt_pos, psite_count, relative_density, codon)]
      )
      if (nrow(global_all_stats) > 0) {
        fwrite(global_all_stats, file.path(codon_plot_dir, "Merge_AllCodons_all_stats.txt"), sep = "\t", na = "")
      }

      codon_color_map <- data.table(codon = names(codon_color), color_hex = as.character(codon_color))
      fwrite(codon_color_map, file.path(codon_plot_dir, "Codon_Color_Map.txt"), sep = "\t", row.names = FALSE)
      cat(paste0("✅ All codon merge plots saved\n"))
    }

    if (length(valid_codons) > 0) {
      all_codon_stats <- reads_psite_dt[P_codon %in% valid_codons][
        , .(total_reads = .N, transcript_count = uniqueN(transcript), avg_occupancy = .N / uniqueN(transcript)), by = P_codon][order(-total_reads)]
      fwrite(all_codon_stats, file.path(codon_plot_dir, "all_codon_summary_stats.txt"), sep = "\t")
    }

    cat(paste0("\n✅ Codon-specific analysis finished\n"))
    cat(paste0("📁 Output dir: ", codon_plot_dir, "\n"))
    cat(paste0("📊 PDF plots: ", length(list.files(codon_plot_dir, pattern = "*.pdf")), "\n"))
    cat(paste0("📈 TXT tables: ", length(list.files(codon_plot_dir, pattern = "*.txt")), "\n"))
    cat(paste0("🔢 Valid Frame 0 codons: ", length(valid_codons), "\n"))

    if (tolower(selected_analysis) == "disome") {
      cat(paste0("\n✅ Disome Rear inference: trim ", trim_nt, " nt, infer Rear Codon and Rear APE sites\n"))
      reads_psite_dt[, trim_3x := (trim_nt %/% 3) * 3]
      cat(paste0("   ✅ Trim rule: trim_nt=", trim_nt, " → 3x multiple=", unique(reads_psite_dt$trim_3x), "\n"))

      reads_psite_dt[, `:=`(
        rear_seq_start = ori_end5,
        rear_seq_end = end5 - 1,
        Rear_first_codon_pos = psite - trim_3x,
        Rear_first_cds_pos = ifelse(!is.na(P_cds_pos), P_codon_pos - trim_3x, NA_integer_),
        Rear_first_codon_num = ifelse(!is.na(P_codon_pos), P_codon_pos - floor(trim_3x / 3) - 1, NA_integer_)
      )]
      reads_psite_dt[, rear_seq_length := rear_seq_end - rear_seq_start + 1]
      reads_psite_dt[, seq_str := fasta_seqs[fasta_idx]]
      reads_psite_dt[, Rear_seq := ifelse(!is.na(seq_str) & rear_seq_start >= 1 & rear_seq_end <= nchar(seq_str), substr(seq_str, rear_seq_start, rear_seq_end), NA_character_)]
      reads_psite_dt[, Rear_Codon := ifelse(!is.na(Rear_seq) & nchar(Rear_seq) >= 3, substr(Rear_seq, nchar(Rear_seq) - 2, nchar(Rear_seq)), NA_character_)]
      reads_psite_dt[, Rear_Codon_next := ifelse(!is.na(Rear_seq) & nchar(Rear_seq) >= 6, substr(Rear_seq, nchar(Rear_seq) - 5, nchar(Rear_seq) - 3), NA_character_)]
      reads_psite_dt[, Rear_Codon_prev := ifelse(!is.na(seq_str) & (rear_seq_start + 3) <= rear_seq_end, substr(seq_str, rear_seq_start, rear_seq_start + 2), NA_character_)]
      reads_psite_dt[, Rear_Codon_valid := !is.na(Rear_Codon) & nchar(Rear_Codon) == 3]

      reads_psite_dt[, `:=`(
        Rear_P_site = psite - trim_3x,
        Rear_A_site = Rear_P_site + 3,
        Rear_E_site = Rear_P_site - 3
      )]
      reads_psite_dt[, `:=`(
        Rear_P_cds_pos = Rear_P_site - cds_start + 1,
        Rear_A_cds_pos = Rear_A_site - cds_start + 1,
        Rear_E_cds_pos = Rear_E_site - cds_start + 1
      )]
      reads_psite_dt[, real_phase := ifelse(Rear_P_site >= cds_start & Rear_P_site <= (cds_stop - 3), (Rear_P_site - cds_start) %% 3, NA_integer_)]
      reads_psite_dt[, corrected_phase := fifelse(real_phase < 0, real_phase + 3, real_phase)]
      reads_psite_dt[, `:=`(
        Rear_A_site = Rear_P_site + 3,
        Rear_E_site = Rear_P_site - 3,
        Rear_A_cds_pos = Rear_A_site - cds_start + 1,
        Rear_E_cds_pos = Rear_E_site - cds_start + 1
      )]
      reads_psite_dt[, `:=`(
        Rear_P_codon = ifelse(!is.na(corrected_phase) & !is.na(seq_str) & Rear_P_site >= cds_start & Rear_P_site <= (cds_stop - 3) & (Rear_P_site + 2) <= nchar(seq_str), substr(seq_str, Rear_P_site, Rear_P_site + 2), NA_character_),
        Rear_A_codon = ifelse(!is.na(Rear_P_codon) & !is.na(seq_str) & Rear_A_site >= cds_start & Rear_A_site <= (cds_stop - 3) & (Rear_A_site + 2) <= nchar(seq_str), substr(seq_str, Rear_A_site, Rear_A_site + 2), NA_character_),
        Rear_E_codon = ifelse(!is.na(Rear_P_codon) & !is.na(seq_str) & Rear_E_site >= cds_start & Rear_E_site <= (cds_stop - 3) & (Rear_E_site + 2) <= nchar(seq_str), substr(seq_str, Rear_E_site, Rear_E_site + 2), NA_character_)
      )]
      reads_psite_dt[, Rear_APE_valid := !is.na(Rear_P_codon) & nchar(Rear_P_codon) == 3 & !is.na(Rear_A_codon) & nchar(Rear_A_codon) == 3 & !is.na(Rear_E_codon) & nchar(Rear_E_codon) == 3]

      reads_psite_dt[, c("fasta_idx", "seq_str", "trim_3x") := NULL]
      valid_rear_ratio <- sum(reads_psite_dt$Rear_Codon_valid, na.rm = TRUE) / sum(!is.na(reads_psite_dt$Rear_Codon)) * 100
      valid_rear_ape_ratio <- sum(reads_psite_dt$Rear_APE_valid, na.rm = TRUE) / nrow(reads_psite_dt) * 100
      cat(paste0("✅ Rear Codon inference done | Valid Rear Codon ratio: ", sprintf("%.2f", valid_rear_ratio), "%\n"))
      cat(paste0("✅ Rear APE derivation done | Valid Rear APE ratio: ", sprintf("%.2f", valid_rear_ape_ratio), "%\n"))
    } else {
      reads_psite_dt[, c("fasta_idx", "seq_str") := NULL]
    }

    reads_psite <- setNames(list(reads_psite_dt), sample_name)
    saveRDS(reads_psite, file.path(sample_output_dir, "reads_psite.rds"))

    cat("Calculating codon occupancy index...\n")
    annotation <- anno_dt
    psite_data <- reads_psite

    p_codon_occupancy <- reads_psite_dt[!is.na(P_codon), .(total_reads = .N, avg_occupancy = .N / uniqueN(transcript), transcript_count = uniqueN(transcript)), by = P_codon]
    p_avg <- mean(p_codon_occupancy$avg_occupancy)
    p_total <- sum(p_codon_occupancy$avg_occupancy)
    p_codon_occupancy[, `:=`(Normalized_avg = avg_occupancy / p_avg, Normalized_avg_by_total = avg_occupancy / p_total)]
    setnames(p_codon_occupancy, "P_codon", "codon")
    p_codon_occupancy[, site := "P"]

    a_codon_occupancy <- reads_psite_dt[!is.na(A_codon), .(total_reads = .N, avg_occupancy = .N / uniqueN(transcript), transcript_count = uniqueN(transcript)), by = A_codon]
    a_avg <- mean(a_codon_occupancy$avg_occupancy)
    a_total <- sum(a_codon_occupancy$avg_occupancy)
    a_codon_occupancy[, `:=`(Normalized_avg = avg_occupancy / a_avg, Normalized_avg_by_total = avg_occupancy / a_total)]
    setnames(a_codon_occupancy, "A_codon", "codon")
    a_codon_occupancy[, site := "A"]

    e_codon_occupancy <- reads_psite_dt[!is.na(E_codon), .(total_reads = .N, avg_occupancy = .N / uniqueN(transcript), transcript_count = uniqueN(transcript)), by = E_codon]
    e_avg <- mean(e_codon_occupancy$avg_occupancy)
    e_total <- sum(e_codon_occupancy$avg_occupancy)
    e_codon_occupancy[, `:=`(Normalized_avg = avg_occupancy / e_avg, Normalized_avg_by_total = avg_occupancy / e_total)]
    setnames(e_codon_occupancy, "E_codon", "codon")
    e_codon_occupancy[, site := "E"]

    if (selected_analysis == "disome") {
      cat("Calculating Rear codon occupancy...\n")
      rear_codon_occupancy <- reads_psite_dt[!is.na(Rear_Codon) & Rear_Codon_valid == TRUE, .(total_reads = .N, avg_occupancy = .N / uniqueN(transcript), transcript_count = uniqueN(transcript)), by = Rear_Codon]
      l_avg <- mean(rear_codon_occupancy$avg_occupancy)
      l_total <- sum(rear_codon_occupancy$avg_occupancy)
      rear_codon_occupancy[, `:=`(Normalized_avg = avg_occupancy / l_avg, Normalized_avg_by_total = avg_occupancy / l_total)]
      setnames(rear_codon_occupancy, "Rear_Codon", "codon")
      rear_codon_occupancy[, site := "Rear_Original"]

      rear_p_occupancy <- reads_psite_dt[!is.na(Rear_P_codon) & Rear_APE_valid == TRUE, .(total_reads = .N, avg_occupancy = .N / uniqueN(transcript), transcript_count = uniqueN(transcript)), by = Rear_P_codon]
      lp_avg <- mean(rear_p_occupancy$avg_occupancy)
      lp_total <- sum(rear_p_occupancy$avg_occupancy)
      rear_p_occupancy[, `:=`(Normalized_avg = avg_occupancy / lp_avg, Normalized_avg_by_total = avg_occupancy / lp_total)]
      setnames(rear_p_occupancy, "Rear_P_codon", "codon")
      rear_p_occupancy[, site := "Rear_P"]

      rear_a_occupancy <- reads_psite_dt[!is.na(Rear_A_codon) & Rear_APE_valid == TRUE, .(total_reads = .N, avg_occupancy = .N / uniqueN(transcript), transcript_count = uniqueN(transcript)), by = Rear_A_codon]
      la_avg <- mean(rear_a_occupancy$avg_occupancy)
      la_total <- sum(rear_a_occupancy$avg_occupancy)
      rear_a_occupancy[, `:=`(Normalized_avg = avg_occupancy / la_avg, Normalized_avg_by_total = avg_occupancy / la_total)]
      setnames(rear_a_occupancy, "Rear_A_codon", "codon")
      rear_a_occupancy[, site := "Rear_A"]

      rear_e_occupancy <- reads_psite_dt[!is.na(Rear_E_codon) & Rear_APE_valid == TRUE, .(total_reads = .N, avg_occupancy = .N / uniqueN(transcript), transcript_count = uniqueN(transcript)), by = Rear_E_codon]
      le_avg <- mean(rear_e_occupancy$avg_occupancy)
      le_total <- sum(rear_e_occupancy$avg_occupancy)
      rear_e_occupancy[, `:=`(Normalized_avg = avg_occupancy / le_avg, Normalized_avg_by_total = avg_occupancy / le_total)]
      setnames(rear_e_occupancy, "Rear_E_codon", "codon")
      rear_e_occupancy[, site := "Rear_E"]

      all_site_occupancy <- rbind(p_codon_occupancy, a_codon_occupancy, e_codon_occupancy, rear_codon_occupancy, rear_p_occupancy, rear_a_occupancy, rear_e_occupancy)
      fwrite(rear_codon_occupancy, file.path(sample_output_dir, "Rear_Codon_occupancy.txt"), sep = "\t")
      fwrite(rear_p_occupancy, file.path(sample_output_dir, "Rear_P_Codon_occupancy_corrected.txt"), sep = "\t")
      fwrite(rear_a_occupancy, file.path(sample_output_dir, "Rear_A_Codon_occupancy_corrected.txt"), sep = "\t")
      fwrite(rear_e_occupancy, file.path(sample_output_dir, "Rear_E_Codon_occupancy_corrected.txt"), sep = "\t")
    } else {
      all_site_occupancy <- rbind(p_codon_occupancy, a_codon_occupancy, e_codon_occupancy)
    }

    fwrite(all_site_occupancy, file.path(sample_output_dir, "all_site_codon_occupancy_final.txt"), sep = "\t")
    fwrite(p_codon_occupancy, file.path(sample_output_dir, "P_site_codon_occupancy.txt"), sep = "\t")
    fwrite(a_codon_occupancy, file.path(sample_output_dir, "A_site_codon_occupancy.txt"), sep = "\t")
    fwrite(e_codon_occupancy, file.path(sample_output_dir, "E_site_codon_occupancy.txt"), sep = "\t")
    cat("✅ Codon occupancy calculation completed\n")

    cat("Generating codon occupancy plots with Normalized_avg centered at 1.00...\n")
    p_valid <- p_codon_occupancy[, .(codon, norm_val = Normalized_avg)][, site := "P"]
    a_valid <- a_codon_occupancy[, .(codon, norm_val = Normalized_avg)][, site := "A"]
    e_valid <- e_codon_occupancy[, .(codon, norm_val = Normalized_avg)][, site := "E"]

    if (selected_analysis == "disome") {
      l_valid <- rear_codon_occupancy[, .(codon, norm_val = Normalized_avg)][, site := "Rear_Original"]
      lp_valid <- rear_p_occupancy[, .(codon, norm_val = Normalized_avg)][, site := "Rear_P"]
      la_valid <- rear_a_occupancy[, .(codon, norm_val = Normalized_avg)][, site := "Rear_A"]
      le_valid <- rear_e_occupancy[, .(codon, norm_val = Normalized_avg)][, site := "Rear_E"]
      all_valid_norm <- rbind(p_valid, a_valid, e_valid, l_valid, lp_valid, la_valid, le_valid)
      fill_color <- c("P" = "#333f50", "A" = "#39827c", "E" = "#e84a5f",
                      "Rear_Original" = "#feca57", "Rear_P" = "#ff9f43",
                      "Rear_A" = "#00d2d3", "Rear_E" = "#a29bfe")
    } else {
      all_valid_norm <- rbind(p_valid, a_valid, e_valid)
      fill_color <- c("P" = "#333f50", "A" = "#39827c", "E" = "#e84a5f")
    }
    all_valid_norm_sorted <- all_valid_norm[order(codon, site)]

    codon2aa_map <- c(
      "TTT"="F","TTC"="F","TTA"="L","TTG"="L","TCT"="S","TCC"="S","TCA"="S","TCG"="S","TAT"="Y","TAC"="Y","TAA"="*","TAG"="*","TGT"="C","TGC"="C","TGA"="*","TGG"="W",
      "CTT"="L","CTC"="L","CTA"="L","CTG"="L","CCT"="P","CCC"="P","CCA"="P","CCG"="P","CAT"="H","CAC"="H","CAA"="Q","CAG"="Q","CGT"="R","CGC"="R","CGA"="R","CGG"="R",
      "ATT"="I","ATC"="I","ATA"="I","ATG"="M","ACT"="T","ACC"="T","ACA"="T","ACG"="T","AAT"="N","AAC"="N","AAA"="K","AAG"="K","AGT"="S","AGC"="S","AGA"="R","AGG"="R",
      "GTT"="V","GTC"="V","GTA"="V","GTG"="V","GCT"="A","GCC"="A","GCA"="A","GCG"="A","GAT"="D","GAC"="D","GAA"="E","GAG"="E","GGT"="G","GGC"="G","GGA"="G","GGG"="G"
    )

    codon_label <- function(cds) paste0(cds, "(", codon2aa_map[as.character(cds)], ")")
    codon_color_fun <- function(cds) ifelse(grepl("[CG]$", as.character(cds)), "red", "black")

    plot_all_sites <- ggplot(all_valid_norm_sorted, aes(x = reorder(codon, -norm_val), y = norm_val, fill = site)) +
      geom_col(position = "dodge", alpha = 0.8, width = 0.8) +
      geom_hline(yintercept = 1.00, color = "black", linewidth = 0.8, linetype = "dashed") +
      scale_fill_manual(values = fill_color) +
      scale_x_discrete(labels = codon_label) +
      labs(x = "Codon", y = "Normalized Codon Occupancy (Centered at 1.00)",
           title = paste(sample_name, "A/P/E Codon Occupancy"),
           subtitle = flow_text, fill = "Ribosome Site") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(all_valid_norm_sorted$codon, -all_valid_norm_sorted$norm_val)))))

    ggsave(file.path(sample_output_dir, "A_P_E_codon_occupancy_Centered_at_Median.pdf"),
           plot = plot_all_sites, width = 16, height = 7, dpi = 300)

    syn_all <- all_valid_norm_sorted
    syn_all$aa <- codon2aa_map[as.character(syn_all$codon)]
    syn_all <- syn_all[order(syn_all$aa, -syn_all$norm_val)]
    syn_all$codon <- factor(syn_all$codon, levels = unique(syn_all$codon))
    syn_all$x_pos <- as.integer(syn_all$codon)
    syn_all$x_pos <- syn_all$x_pos + c(0, cumsum(diff(as.integer(factor(syn_all$aa))) != 0)) * 0.8

    p_syn_all <- ggplot(syn_all, aes(x = x_pos, y = norm_val, fill = site)) +
      geom_col(position = "dodge", alpha = 0.8, width = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
      scale_fill_manual(values = fill_color) +
      scale_x_continuous(breaks = syn_all$x_pos, labels = codon_label(syn_all$codon), expand = c(0.01, 0.01)) +
      labs(x = "Codon (Grouped by Synonym)", y = "Normalized Codon Occupancy (Centered at 1.00)",
           title = paste(sample_name, "Synonymous A/P/E/Rear"), fill = "Ribosome Site") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_all$codon)))

    ggsave(file.path(sample_output_dir, "A_P_E_codon_occupancy_SYNONYMOUS.pdf"),
           plot = p_syn_all, width = 18, height = 7, dpi = 300)

    p_plot <- ggplot(p_valid[order(-norm_val)], aes(x = reorder(codon, -norm_val), y = norm_val, fill = ifelse(norm_val > 1.00, "#333f50", "#e84a5f"))) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 1.00, linetype = "dashed", color = "black", linewidth = 0.8) +
      scale_x_discrete(labels = codon_label) +
      labs(x = "Codon", y = "Normalized Occupancy (Centered at 1.00)", title = "P-site Codon Occupancy") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(p_valid$codon, -p_valid$norm_val)))))

    ggsave(file.path(sample_output_dir, "P_site_codon_occupancy_Centered_at_Median.pdf"),
           plot = p_plot, width = 12, height = 6, dpi = 300)

    syn_p <- p_valid
    syn_p$aa <- codon2aa_map[as.character(syn_p$codon)]
    syn_p <- syn_p[order(syn_p$aa, -syn_p$norm_val)]
    syn_p$codon <- factor(syn_p$codon, levels = unique(syn_p$codon))
    syn_p$x_pos <- as.integer(syn_p$codon)
    syn_p$x_pos <- syn_p$x_pos + c(0, cumsum(diff(as.integer(factor(syn_p$aa))) != 0)) * 0.8

    p_syn_p <- ggplot(syn_p, aes(x = x_pos, y = norm_val, fill = ifelse(norm_val > 1, "#333f50", "#e84a5f"))) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
      scale_x_continuous(breaks = syn_p$x_pos, labels = codon_label(syn_p$codon), expand = c(0.01, 0.01)) +
      labs(x = "Codon (Grouped by Synonym)", y = "Normalized Occupancy", title = "P-site Synonymous Codon Occupancy") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_p$codon)))

    ggsave(file.path(sample_output_dir, "P_site_codon_occupancy_SYNONYMOUS.pdf"),
           plot = p_syn_p, width = 14, height = 6, dpi = 300)

    a_plot <- ggplot(a_valid[order(-norm_val)], aes(x = reorder(codon, -norm_val), y = norm_val, fill = ifelse(norm_val > 1.00, "#39827c", "#e84a5f"))) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 1.00, linetype = "dashed", color = "black", linewidth = 0.8) +
      scale_x_discrete(labels = codon_label) +
      labs(x = "Codon", y = "Normalized Occupancy (Centered at 1.00)", title = "A-site Codon Occupancy") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(a_valid$codon, -a_valid$norm_val)))))

    ggsave(file.path(sample_output_dir, "A_site_codon_occupancy_Centered_at_Median.pdf"),
           plot = a_plot, width = 12, height = 6, dpi = 300)

    syn_a <- a_valid
    syn_a$aa <- codon2aa_map[as.character(syn_a$codon)]
    syn_a <- syn_a[order(syn_a$aa, -syn_a$norm_val)]
    syn_a$codon <- factor(syn_a$codon, levels = unique(syn_a$codon))
    syn_a$x_pos <- as.integer(syn_a$codon)
    syn_a$x_pos <- syn_a$x_pos + c(0, cumsum(diff(as.integer(factor(syn_a$aa))) != 0)) * 0.8

    p_syn_a <- ggplot(syn_a, aes(x = x_pos, y = norm_val, fill = ifelse(norm_val > 1, "#39827c", "#e84a5f"))) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
      scale_x_continuous(breaks = syn_a$x_pos, labels = codon_label(syn_a$codon), expand = c(0.01, 0.01)) +
      labs(x = "Codon (Grouped by Synonym)", y = "Normalized Occupancy", title = "A-site Synonymous Codon Occupancy") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_a$codon)))

    ggsave(file.path(sample_output_dir, "A_site_codon_occupancy_SYNONYMOUS.pdf"),
           plot = p_syn_a, width = 14, height = 6, dpi = 300)

    e_plot <- ggplot(e_valid[order(-norm_val)], aes(x = reorder(codon, -norm_val), y = norm_val, fill = ifelse(norm_val > 1.00, "#333f50", "#e84a5f"))) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 1.00, linetype = "dashed", color = "black", linewidth = 0.8) +
      scale_x_discrete(labels = codon_label) +
      labs(x = "Codon", y = "Normalized Occupancy (Centered at 1.00)", title = "E-site Codon Occupancy") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(e_valid$codon, -e_valid$norm_val)))))

    ggsave(file.path(sample_output_dir, "E_site_codon_occupancy_Centered_at_Median.pdf"),
           plot = e_plot, width = 12, height = 6, dpi = 300)

    syn_e <- e_valid
    syn_e$aa <- codon2aa_map[as.character(syn_e$codon)]
    syn_e <- syn_e[order(syn_e$aa, -syn_e$norm_val)]
    syn_e$codon <- factor(syn_e$codon, levels = unique(syn_e$codon))
    syn_e$x_pos <- as.integer(syn_e$codon)
    syn_e$x_pos <- syn_e$x_pos + c(0, cumsum(diff(as.integer(factor(syn_e$aa))) != 0)) * 0.8

    p_syn_e <- ggplot(syn_e, aes(x = x_pos, y = norm_val, fill = ifelse(norm_val > 1, "#333f50", "#e84a5f"))) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
      scale_x_continuous(breaks = syn_e$x_pos, labels = codon_label(syn_e$codon), expand = c(0.01, 0.01)) +
      labs(x = "Codon (Grouped by Synonym)", y = "Normalized Occupancy", title = "E-site Synonymous Codon Occupancy") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_e$codon)))

    ggsave(file.path(sample_output_dir, "E_site_codon_occupancy_SYNONYMOUS.pdf"),
           plot = p_syn_e, width = 14, height = 6, dpi = 300)

    if (selected_analysis == "disome") {
      l_plot <- ggplot(l_valid[order(-norm_val)], aes(x = reorder(codon, -norm_val), y = norm_val, fill = ifelse(norm_val > 1.00, "#feca57", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1.00, linetype = "dashed", color = "black", linewidth = 0.8) +
        scale_x_discrete(labels = codon_label) +
        labs(x = "Codon", y = "Normalized Occupancy (Centered at 1.00)", title = "Rear Codon Occupancy") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(l_valid$codon, -l_valid$norm_val)))))

      ggsave(file.path(sample_output_dir, "Rear_Codon_occupancy_Centered_at_Median.pdf"),
             plot = l_plot, width = 12, height = 6, dpi = 300)

      syn_l <- l_valid
      syn_l$aa <- codon2aa_map[as.character(syn_l$codon)]
      syn_l <- syn_l[order(syn_l$aa, -syn_l$norm_val)]
      syn_l$codon <- factor(syn_l$codon, levels = unique(syn_l$codon))
      syn_l$x_pos <- as.integer(syn_l$codon)
      syn_l$x_pos <- syn_l$x_pos + c(0, cumsum(diff(as.integer(factor(syn_l$aa))) != 0)) * 0.8

      p_syn_l <- ggplot(syn_l, aes(x = x_pos, y = norm_val, fill = ifelse(norm_val > 1, "#feca57", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
        scale_x_continuous(breaks = syn_l$x_pos, labels = codon_label(syn_l$codon), expand = c(0.01, 0.01)) +
        labs(x = "Codon (Grouped by Synonym)", title = "Rear Synonymous Codon Occupancy") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_l$codon)))

      ggsave(file.path(sample_output_dir, "Rear_Codon_occupancy_SYNONYMOUS.pdf"),
             plot = p_syn_l, width = 14, height = 6, dpi = 300)

      lp_plot <- ggplot(lp_valid[order(-norm_val)], aes(x = reorder(codon, -norm_val), y = norm_val, fill = ifelse(norm_val > 1.00, "#ff9f43", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1.00, linetype = "dashed", color = "black", linewidth = 0.8) +
        scale_x_discrete(labels = codon_label) +
        labs(x = "Codon", y = "Normalized Occupancy (Centered at 1.00)", title = "Rear P-Codon Occupancy (3nt Corrected)") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(lp_valid$codon, -lp_valid$norm_val)))))

      ggsave(file.path(sample_output_dir, "Rear_P_Codon_occupancy_Centered_at_Median.pdf"),
             plot = lp_plot, width = 12, height = 6, dpi = 300)

      syn_lp <- lp_valid
      syn_lp$aa <- codon2aa_map[as.character(syn_lp$codon)]
      syn_lp <- syn_lp[order(syn_lp$aa, -syn_lp$norm_val)]
      syn_lp$codon <- factor(syn_lp$codon, levels = unique(syn_lp$codon))
      syn_lp$x_pos <- as.integer(syn_lp$codon)
      syn_lp$x_pos <- syn_lp$x_pos + c(0, cumsum(diff(as.integer(factor(syn_lp$aa))) != 0)) * 0.8

      p_syn_lp <- ggplot(syn_lp, aes(x = x_pos, y = norm_val, fill = ifelse(norm_val > 1, "#ff9f43", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
        scale_x_continuous(breaks = syn_lp$x_pos, labels = codon_label(syn_lp$codon), expand = c(0.01, 0.01)) +
        labs(x = "Codon (Grouped by Synonym)", title = "Rear P Synonymous Codon Occupancy") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_lp$codon)))

      ggsave(file.path(sample_output_dir, "Rear_P_Codon_occupancy_SYNONYMOUS.pdf"),
             plot = p_syn_lp, width = 14, height = 6, dpi = 300)

      la_plot <- ggplot(la_valid[order(-norm_val)], aes(x = reorder(codon, -norm_val), y = norm_val, fill = ifelse(norm_val > 1.00, "#00d2d3", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1.00, linetype = "dashed", color = "black", linewidth = 0.8) +
        scale_x_discrete(labels = codon_label) +
        labs(x = "Codon", y = "Normalized Occupancy (Centered at 1.00)", title = "Rear A-Codon Occupancy (3nt Corrected)") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(la_valid$codon, -la_valid$norm_val)))))

      ggsave(file.path(sample_output_dir, "Rear_A_Codon_occupancy_Centered_at_Median.pdf"),
             plot = la_plot, width = 12, height = 6, dpi = 300)

      syn_la <- la_valid
      syn_la$aa <- codon2aa_map[as.character(syn_la$codon)]
      syn_la <- syn_la[order(syn_la$aa, -syn_la$norm_val)]
      syn_la$codon <- factor(syn_la$codon, levels = unique(syn_la$codon))
      syn_la$x_pos <- as.integer(syn_la$codon)
      syn_la$x_pos <- syn_la$x_pos + c(0, cumsum(diff(as.integer(factor(syn_la$aa))) != 0)) * 0.8

      p_syn_la <- ggplot(syn_la, aes(x = x_pos, y = norm_val, fill = ifelse(norm_val > 1, "#00d2d3", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
        scale_x_continuous(breaks = syn_la$x_pos, labels = codon_label(syn_la$codon), expand = c(0.01, 0.01)) +
        labs(x = "Codon (Grouped by Synonym)", title = "Rear A Synonymous Codon Occupancy") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_la$codon)))

      ggsave(file.path(sample_output_dir, "Rear_A_Codon_occupancy_SYNONYMOUS.pdf"),
             plot = p_syn_la, width = 14, height = 6, dpi = 300)

      le_plot <- ggplot(le_valid[order(-norm_val)], aes(x = reorder(codon, -norm_val), y = norm_val, fill = ifelse(norm_val > 1.00, "#a29bfe", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1.00, linetype = "dashed", color = "black", linewidth = 0.8) +
        scale_x_discrete(labels = codon_label) +
        labs(x = "Codon", y = "Normalized Occupancy (Centered at 1.00)", title = "Rear E-Codon Occupancy (3nt Corrected)") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(levels(reorder(le_valid$codon, -le_valid$norm_val)))))

      ggsave(file.path(sample_output_dir, "Rear_E_Codon_occupancy_Centered_at_Median.pdf"),
             plot = le_plot, width = 12, height = 6, dpi = 300)

      syn_le <- le_valid
      syn_le$aa <- codon2aa_map[as.character(syn_le$codon)]
      syn_le <- syn_le[order(syn_le$aa, -syn_le$norm_val)]
      syn_le$codon <- factor(syn_le$codon, levels = unique(syn_le$codon))
      syn_le$x_pos <- as.integer(syn_le$codon)
      syn_le$x_pos <- syn_le$x_pos + c(0, cumsum(diff(as.integer(factor(syn_le$aa))) != 0)) * 0.8

      p_syn_le <- ggplot(syn_le, aes(x = x_pos, y = norm_val, fill = ifelse(norm_val > 1, "#a29bfe", "#e84a5f"))) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
        scale_x_continuous(breaks = syn_le$x_pos, labels = codon_label(syn_le$codon), expand = c(0.01, 0.01)) +
        labs(x = "Codon (Grouped by Synonym)", title = "Rear E Synonymous Codon Occupancy") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, color = codon_color_fun(syn_le$codon)))

      ggsave(file.path(sample_output_dir, "Rear_E_Codon_occupancy_SYNONYMOUS.pdf"),
             plot = p_syn_le, width = 14, height = 6, dpi = 300)
    }

    cat("✅ All codon occupancy plots generated\n")

  }, error = function(e) {
    error_log <- rbind(error_log, data.table(sample_name = sample_name, error_message = e$message))
    cat(paste0("❌ Error processing sample ", sample_name, ": ", e$message, "\n"))
  })
}

cat("\n============================================================\n")
cat("                    Batch Analysis Summary\n")
cat("============================================================\n")
cat("✅ Successfully processed samples: ", length(trans_bam_paths) - nrow(error_log), "\n")
cat("❌ Failed samples: ", nrow(error_log), "\n")
if (nrow(error_log) > 0) {
  cat("Failed sample details:\n")
  print(error_log)
}
cat("\n✅ All results saved to: ", output_root, "\n")

fwrite(length_choice_log, file.path(output_root, "sample_length_choice_log.txt"), sep = "\t")
if (nrow(error_log) > 0) fwrite(error_log, file.path(output_root, "sample_error_log.txt"), sep = "\t")

cat("\n🎉 riboWaltz batch analysis pipeline completed successfully!\n")

sink()
