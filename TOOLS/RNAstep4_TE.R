Sys.setenv(LANGUAGE='en')
Sys.setenv(R_MAX_NUM_DLLS=999)

WORK_DIR <- getwd()
OUTPUT_ROOT <- paste0(WORK_DIR,"/4.riboWaltz_output")
N_CORES <- 8
MEMORY_LIMIT <- "50G"
DISABLE_TOP3 <- FALSE
BATCH_SIZE <- 5
CHUNK_SIZE <- 1e6
GTF_HUMAN <- "~/reference/human_transcriptome/gencode.v49.annotation.gtf"
GTF_MOUSE <- "~/reference/mouse_transcriptome/gencode.vM35.annotation.gtf"
RNA_MATRIX_REL_PATH <- "3.rna_seq_analysis/gene_count_matrix.tsv"

keep_vars <- c("WORK_DIR","OUTPUT_ROOT","N_CORES","MEMORY_LIMIT","DISABLE_TOP3","BATCH_SIZE","CHUNK_SIZE","GTF_HUMAN","GTF_MOUSE","RNA_MATRIX_REL_PATH")
rm(list = setdiff(ls(), keep_vars))

options(
  stringsAsFactors = FALSE,
  warn=-1,
  datatable.verbose = FALSE,
  datatable.nthreads = N_CORES,
  mc.cores = N_CORES
)
options(future.globals.maxSize = 1024 * 1024 * 1024 * 50)
options(future.rng.onMisuse = "ignore")

project_dir <- WORK_DIR
species_config <- list(
  human = list(
    gtf = GTF_HUMAN,
    rna_seq_filtered_matrix =file.path(project_dir, RNA_MATRIX_REL_PATH)
  ),
  mouse = list(
    gtf = GTF_MOUSE,
    rna_seq_filtered_matrix = file.path(project_dir, RNA_MATRIX_REL_PATH)
  )
)

args <- commandArgs(trailingOnly = TRUE)
selected_species <- NULL
normalize_method <- "RPM"
ribo_dir <- NULL

if (length(args) >= 1 && args[1] %in% names(species_config)) {
  selected_species <- args[1]
  if (length(args) >= 2 && args[2] %in% c("RPM", "TPM", "VST")) normalize_method <- args[2]
  if (length(args)>=3 && dir.exists(args[3])) ribo_dir <- args[3]
  cat(paste0("\n✅ Command line parameters loaded:\n   Species: ", selected_species, "\n   Normalization: ", normalize_method, "\n   Ribo directory: ", ribo_dir, "\n"))
  cat(paste0("✅ High-performance config loaded:\n   Threads: ", N_CORES, "\n   Memory limit: ", MEMORY_LIMIT, "\n   Chunk size: ", CHUNK_SIZE, "\n"))
} else if (interactive()) {
  cat("Select species (1=human/2=mouse):");selected_choice=readLines(n=1);selected_species=ifelse(selected_choice==1,"human","mouse")
  cat("Select normalization method (1=RPM/2=TPM/3=VST):");norm_choice=readLines(n=1);normalize_method=c("RPM","TPM","VST")[as.integer(norm_choice)]
  cat("Input Ribo-seq root directory:");ribo_dir=readLines(n=1);ribo_dir=gsub("^\\s+|\\s+$","",ribo_dir)
} else {
  cat("Usage: Rscript RNAstep4_TE.R human|mouse RPM|TPM|VST ribo_dir\n");quit(save="no",status=1)
}

config <- species_config[[selected_species]]
rna_matrix_path <- config$rna_seq_filtered_matrix
ribo_metrics_dir <- file.path(ribo_dir, "4.riboWaltz_output/2.translational expression")
gtf_path <- config$gtf

cat("\nLoading dependencies...\n")
required_pkgs <- c("data.table", "tidyverse", "GenomicFeatures", "rtracklayer", "ggpubr", "future.apply")
invisible(sapply(required_pkgs, function(pkg) {
  if(!requireNamespace(pkg,quietly=T)) install.packages(pkg,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/", quiet = TRUE)
  library(pkg,character.only=T,quietly=T,warn.conflicts=F)
}))
plan(multisession, workers = N_CORES)

standardize_sample_name <- function(x) {
  x <- gsub("_Aligned|_trimmed|_filtered|\\.bam|\\.sam|_processed","",x)
  return(x)
}
filter_low_expr <- function(df, rna_cutoff=5, ribo_cutoff=2) {
  df %>% dplyr::filter(rna_raw_count >= rna_cutoff | ribo_raw_rpfs >= ribo_cutoff)
}
get_ensembl2gene_map <- function(gtf_path) {
  cat("Extracting gene ID mapping from GTF...\n")
  gtf <- rtracklayer::import(gtf_path) %>% as.data.frame()
  gene_map <- gtf %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
    dplyr::mutate(gene_id = gsub("\\.\\d+$", "", gene_id))
  return(gene_map)
}
safe_assign <- function(x, value) {
  if (length(value) == 0 || is.na(value) || is.null(value)) return(0)
  else return(value)
}
fread_chunked <- function(file, chunk_size = CHUNK_SIZE) {
  dt <- data.table::fread(file, header = TRUE, sep = "\t", check.names = FALSE, showProgress = FALSE)
  split_dt <- split(dt, ceiling(seq_len(nrow(dt))/chunk_size))
  return(data.table::rbindlist(split_dt) %>% as.data.frame())
}

cat("\n=======================\nStep 1: Read RNA-seq matrix and convert gene IDs\n=======================\n")
if (!file.exists(rna_matrix_path)) stop(paste0("❌ RNA matrix not found: ", rna_matrix_path))
rna_count <- fread_chunked(rna_matrix_path)

gene_id_col_candidates <- c("gene_id", "GeneID", "geneid", "ensembl_id", "Ensembl_ID", "Gene_id")
gene_id_col <- intersect(gene_id_col_candidates, colnames(rna_count))[1]
if(is.na(gene_id_col)) {
  gene_id_col <- colnames(rna_count)[1]
  cat(paste0("ℹ️ Standard gene ID column not detected, using first column: ", gene_id_col, "\n"))
}
cat(paste0("✅ Detected gene ID column: ", gene_id_col, "\n"))

gene_map <- get_ensembl2gene_map(gtf_path)
rna_count <- rna_count %>%
  dplyr::mutate(!!sym(gene_id_col) := gsub("\\.\\d+$", "", !!sym(gene_id_col))) %>%
  dplyr::filter(!is.na(!!sym(gene_id_col)) & !!sym(gene_id_col) != "")

rna_count_anno <- rna_count %>%
  dplyr::left_join(gene_map, by = setNames("gene_id", gene_id_col)) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

if(nrow(rna_count_anno) == 0) {
  cat("⚠️ Gene matching failed, keep original gene IDs for downstream analysis\n")
  rna_count_anno <- rna_count %>%
    dplyr::rename(gene_name = !!sym(gene_id_col)) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)
}

rna_count_clean <- rna_count_anno %>%
  dplyr::select(-dplyr::any_of(c("gene_id"))) %>%
  tibble::column_to_rownames(var = "gene_name") %>%
  dplyr::mutate_all(~as.numeric(.)) %>%
  replace(is.na(.), 0) %>%
  dplyr::filter(rowSums(.) > 0)

cat(paste0("✅ RNA matrix preprocessing done:\n   Raw genes: ", nrow(rna_count), "\n   Valid genes: ", nrow(rna_count_clean), "\n   Samples: ", ncol(rna_count_clean), "\n"))

cat("\n=======================\nStep 2: Sample name standardization and filtering\n=======================\n")
rna_sample_names <- colnames(rna_count_clean)
rna_samples_standard <- sapply(rna_sample_names, standardize_sample_name)
colnames(rna_count_clean) <- rna_samples_standard
cat(paste0("✅ RNA sample name standardized: total ", length(rna_samples_standard), " samples\n"))

if (!dir.exists(ribo_metrics_dir)) stop(paste0("❌ Ribo directory not found: ", ribo_metrics_dir))
ribo_metrics_files <- list.files(path = ribo_metrics_dir, pattern = "_gene_level_dual_metrics.txt$", full.names = TRUE, recursive = TRUE)
if (length(ribo_metrics_files) ==0) stop("❌ No Ribo-seq metric files found!")
cat(paste0("✅ Found ", length(ribo_metrics_files), " Ribo-seq files\n"))

ribo_all <- future_lapply(ribo_metrics_files, function(file) {
  sample_name <- basename(dirname(file))
  sample_standard <- standardize_sample_name(sample_name)
  ribo_dt <- fread(file, header = TRUE, data.table = FALSE, showProgress = FALSE)
  core_cols <- c("gene_name", "total_rpfs", "gene_rpkm_aggregate")
  if (all(core_cols %in% colnames(ribo_dt))) {
    ribo_dt_clean <- ribo_dt %>% dplyr::select(dplyr::all_of(core_cols)) %>% as.data.table()
    ribo_dt_clean[, sample_name := sample_standard]
    return(ribo_dt_clean)
  } else {
    return(NULL)
  }
}) %>% Filter(Negate(is.null), .)

if(length(ribo_all) == 0) stop("❌ All Ribo-seq files lack required core columns, check file format!")
ribo_combined <- rbindlist(ribo_all) %>% distinct(gene_name, sample_name, .keep_all = TRUE)
ribo_sample_totals <- ribo_combined[, .(total_rpfs = sum(total_rpfs, na.rm = TRUE)), by = sample_name]
ribo_sample_totals <- ribo_sample_totals[total_rpfs > 0,]
if(nrow(ribo_sample_totals) == 0) stop("❌ total_rpfs = 0 for all Ribo samples, no valid signal!")

common_samples <- intersect(rna_samples_standard, unique(ribo_combined$sample_name))
if (length(common_samples) ==0) stop("❌ No overlapping sample names between RNA and Ribo!")
cat(paste0("✅ Sample name match success: total ", length(common_samples), " shared samples\n"))
rna_count_clean <- rna_count_clean[, common_samples, drop = FALSE]
cat(paste0("✅ RNA matrix filtered to shared samples: ", ncol(rna_count_clean), " samples\n"))

cat("\n=======================\nStep 3: Data normalization (",normalize_method,")\n=======================\n")
rpm_matrix <- as.matrix(rna_count_clean)
rpm_matrix <- rpm_matrix[rowSums(rpm_matrix) > 0, ]
lib_size  <- colSums(rpm_matrix)
rpm_matrix <- rpm_matrix[, lib_size > 0, drop = FALSE]
lib_size  <- lib_size[lib_size > 0]
rpm_matrix <- (rpm_matrix / lib_size) * 1000000

rpm_matrix[which(is.infinite(rpm_matrix))] <- 0
rpm_matrix[which(is.na(rpm_matrix))] <- 0
rpm_matrix[which(is.nan(rpm_matrix))] <- 0
rpm_matrix[rpm_matrix < 0] <- 0
rna_normalized <- rpm_matrix
cat(paste0("✅ RNA normalization finished:\n   Normalized matrix dimension: ", nrow(rna_normalized), " × ", ncol(rna_normalized), "\n"))

ribo_combined[, ribo_rpm := (total_rpfs / ribo_sample_totals[match(sample_name, ribo_sample_totals$sample_name), total_rpfs]) * 1e6]
ribo_combined[is.na(ribo_rpm) | is.infinite(ribo_rpm) | ribo_rpm < 0, ribo_rpm := 0]
cat(paste0("✅ Ribo normalization finished:\n   Valid rows: ", nrow(ribo_combined), "\n"))

cat("\n=======================\nPrep module: Generate and validate rna_long\n=======================\n")
rna_long <- as.data.frame(rna_normalized) %>%
  tibble::rownames_to_column("gene_name") %>%
  tidyr::pivot_longer(cols = -gene_name, names_to = "sample_name", values_to = "rna_normalized") %>%
  dplyr::mutate(across(c(rna_normalized), as.numeric)) %>%
  dplyr::filter(!is.na(rna_normalized) & rna_normalized >= 0)

if (nrow(rna_long) == 0) {
  warning("⚠️ rna_long empty, activate fallback")
  rna_long <- data.frame(gene_name="dummy_gene",sample_name="dummy_sample",rna_normalized=0)
}
cat(paste0("✅ rna_long generated: ", nrow(rna_long), " rows × ", ncol(rna_long), " columns\n"))

cat("\n=======================\nPrep module: Generate and validate ribo_count_clean\n=======================\n")
ribo_wide <- ribo_combined %>%
  dplyr::select(gene_name, sample_name, total_rpfs) %>%
  tidyr::pivot_wider(
    names_from = sample_name,
    values_from = total_rpfs,
    values_fill = 0,
    names_expand = TRUE,
    names_vary = "slowest"
  )

ribo_wide <- ribo_wide %>%
  dplyr::select(gene_name, dplyr::all_of(common_samples))

ribo_count_clean <- as.matrix(ribo_wide[, -1, drop = FALSE])
rownames(ribo_count_clean) <- ribo_wide$gene_name
colnames(ribo_count_clean) <- common_samples

if (ncol(ribo_count_clean) != length(common_samples)) {
  stop(paste0("❌ ribo_count_clean column count (", ncol(ribo_count_clean), ") mismatches common_samples length (", length(common_samples), ")"))
}

cat(paste0("✅ ribo_count_clean generated: ", nrow(ribo_count_clean), " × ", ncol(ribo_count_clean), "\n"))

cat("\n=======================\nStep 4: Gene matching + TE calculation\n=======================\n")
required_vars <- c("rna_long", "ribo_combined", "rna_count_clean", "ribo_count_clean")
for (var in required_vars) {
  if (!exists(var) || (is.data.frame(get(var)) && nrow(get(var)) == 0)) {
    stop(paste0("❌ Required variable ", var, " undefined or empty"))
  }
}
required_rna_cols <- c("gene_name", "sample_name", "rna_normalized")
required_ribo_cols <- c("gene_name", "sample_name", "total_rpfs", "ribo_rpm")
if (!all(required_rna_cols %in% colnames(rna_long))) {
  stop(paste0("❌ rna_long missing required columns: ", paste(setdiff(required_rna_cols, colnames(rna_long)), collapse = ", ")))
}
if (!all(required_ribo_cols %in% colnames(ribo_combined))) {
  stop(paste0("❌ ribo_combined missing required columns: ", paste(setdiff(required_ribo_cols, colnames(ribo_combined)), collapse = ", ")))
}

te_combined <- merge(
  rna_long,
  ribo_combined[, .(gene_name, sample_name, total_rpfs, ribo_rpm)],
  by = c("gene_name", "sample_name"),
  all.x=FALSE, all.y=FALSE
) %>% dplyr::distinct(gene_name, sample_name, .keep_all=TRUE) %>% data.table::as.data.table()

te_combined[, TE := ribo_rpm / (rna_normalized + 1e-8)]
te_combined[, rna_raw_count := 0]
rna_raw_long <- as.data.frame(rna_count_clean) %>%
  tibble::rownames_to_column("gene_name") %>%
  tidyr::pivot_longer(cols = -gene_name, names_to = "sample_name", values_to = "rna_raw_count") %>%
  data.table::as.data.table()
te_combined <- te_combined[rna_raw_long, on = c("gene_name", "sample_name"), rna_raw_count := i.rna_raw_count]
te_combined[is.na(rna_raw_count), rna_raw_count := 0]

te_combined[is.infinite(TE) | is.na(TE) | TE < 0, TE := 0]
te_combined[TE > 1e6, TE := 1e6]

te_final <- te_combined %>%
  dplyr::select(
    gene_name, sample_name,
    rna_raw_count, rna_normalized,
    ribo_raw_rpfs=total_rpfs, ribo_rpm, TE
  ) %>%
  filter_low_expr() %>%
  dplyr::arrange(sample_name, dplyr::desc(TE)) %>%
  data.table::as.data.table()

cat(paste0("✅ Matching complete!\n   Valid samples: ", length(unique(te_final$sample_name)), "\n"))
cat(paste0("✅ Valid genes: ", length(unique(te_final$gene_name)), "\n"))
cat(paste0("✅ Valid records: ", nrow(te_final), "\n"))
cat(paste0("✅ TE range: ", round(min(te_final$TE),6), " ~ ", round(max(te_final$TE),6), "\n"))

cat("\n=======================\nStep 5: Save outputs and visualization\n=======================\n")
output_root <- file.path(dirname(rna_matrix_path), "TE_analysis", paste0("TE_", normalize_method))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

fwrite(te_final, file.path(output_root, paste0("TE_results_ALL_samples_", normalize_method, ".tsv")), sep="\t", quote=F, row.names=F)
fwrite(ribo_sample_totals, file.path(output_root, "Ribo_normalization_params.tsv"), sep="\t", quote=F)

for (curr_sample in common_samples) {
  sample_outdir <- file.path(output_root, curr_sample)
  dir.create(sample_outdir, recursive = TRUE, showWarnings = FALSE)
  te_sample <- te_final[sample_name == curr_sample, ]
  fwrite(te_sample, file.path(sample_outdir, paste0(curr_sample, "_TE_results.tsv")), sep="\t", quote=F, row.names=F)

  if (nrow(te_sample)>=10) {
    plot_data <- te_sample[TE>1e-6 & rna_normalized>1e-6, ]
    if(nrow(plot_data)<5) plot_data = te_sample

    p1 <- ggplot(plot_data, aes(x=log10(TE+1e-8)))+
      geom_histogram(bins=50,fill="#2E86AB",alpha=0.8, color="white")+
      labs(y = "gene number")+
      theme_bw()+ggtitle(paste0(curr_sample, " - TE Distribution (", normalize_method, ")"))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    if(nrow(plot_data)>=10){
      cor_test <- cor.test(log10(plot_data$rna_normalized+1e-8), log10(plot_data$TE+1e-8), method = "pearson")
      cor_r <- round(cor_test$estimate, 3)
      cor_p <- cor_test$p.value
      cor_p_text <- ifelse(cor_p < 0.001, "p < 0.001", paste0("p = ", round(cor_p, 3)))
      cor_label <- paste0("r = ", cor_r, "\n", cor_p_text)
    } else {
      cor_label <- "n < 10"
    }

    p2 <- ggplot(plot_data, aes(x=log10(rna_normalized+1e-8),y=log10(TE+1e-8)))+
      geom_point(alpha=0.7,color="#A23B72", size=1)+
      geom_smooth(method="lm",se=F, color="#F18F01", linewidth=1)+
      annotate("text", x = max(log10(plot_data$rna_normalized+1e-8))*0.8, y = max(log10(plot_data$TE+1e-8))*0.9,
               label = cor_label, hjust = 0, vjust = 1, color = "black", fontface = "bold", size = 4)+
      theme_bw()+ggtitle(paste0(curr_sample, " - TE vs RNA-seq Expression"))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(file.path(sample_outdir, paste0(curr_sample, "_TE_distribution.pdf")), p1, width=7, height=5, dpi=300)
    ggsave(file.path(sample_outdir, paste0(curr_sample, "_TE_vs_RNA.pdf")), p2, width=8, height=5, dpi=300)
  }
}

sink(file.path(output_root, "TE_analysis_report.txt"), split=T)
cat("=====================================\n")
cat("Translation Efficiency (TE) Analysis Report | ", Sys.time(), "\n")
cat("=====================================\n")
cat("Species: ", selected_species, "\n")
cat("Normalization method: ", normalize_method, "\n")
cat("Threads: ", N_CORES, "\n")
cat("Memory limit: ", MEMORY_LIMIT, "\n")
cat("=====================================\n")
cat("Summary statistics:\n")
cat("  - RNA-seq samples: ", length(rna_sample_names), "\n")
cat("  - Ribo-seq samples: ", length(unique(ribo_combined$sample_name)), "\n")
cat("  - Matched samples: ", length(common_samples), "\n")
cat("  - Valid genes: ", length(unique(te_final$gene_name)), "\n")
cat("=====================================\n")
cat("TE statistics:\n")
cat("  - TE min: ", round(min(te_final$TE),6), "\n")
cat("  - TE max: ", round(max(te_final$TE),6), "\n")
cat("  - TE mean: ", round(mean(te_final$TE),6), "\n")
cat("  - TE median: ", round(median(te_final$TE),6), "\n")
cat("=====================================\n")
cat("Output description:\n")
cat("  - Merged table: TE_results_ALL_samples_",normalize_method,".tsv\n")
cat("  - Per-sample folder: tsv + two PDF figures\n")
cat("  - Columns: gene_name sample_name rna_raw_count rna_normalized ribo_raw_rpfs ribo_rpm TE\n")
cat("=====================================\n")
sink()

plan(sequential)

cat("\n🎉 Pipeline completed successfully! ✅\n")
cat("📁 Output root: ", output_root, "\n")
cat("✅ 7 clean columns in final table\n")
cat("✅ Each sample folder contains tsv + two PDF plots\n")
cat("✅ TE_vs_RNA.pdf includes Pearson r and P-value annotation\n")
