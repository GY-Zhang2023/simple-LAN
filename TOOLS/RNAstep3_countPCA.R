#!/usr/bin/env Rscript
Sys.setenv(LANGUAGE='en')
Sys.setenv(R_MAX_NUM_DLLS=999)
options(
  stringsAsFactors = FALSE,
  warn=-1,
  datatable.verbose = FALSE,
  datatable.nthreads = 8,
  mc.cores = 8,
  scipen = 999
)
options(future.globals.maxSize = 1024 * 1024 * 1024 * 50)
options(future.rng.onMisuse = "ignore")

WORK_DIR <- getwd()
N_CORES <- 8
MEMORY_LIMIT <- "50G"
CHUNK_SIZE <- 1e30
padj_cutoff_def <- 0.01
fc_cutoff_def <- 1

rm(list = setdiff(ls(), c("WORK_DIR","N_CORES","MEMORY_LIMIT","CHUNK_SIZE","padj_cutoff_def","fc_cutoff_def")))
library(scales)
options(repos = c(CRAN = 'http://10.10.11.66:8182'))
options(BioC_mirror = 'http://10.10.11.66:8183')

color_palette <- list(
  group_colors = c("#E63946", "#457B9D"),
  volcano_colors = c("#457B9D", "#AAAAAA", "#CCCCCC", "#E63946"),
  heatmap_colors = colorRampPalette(c("#219EBC", "#FFFFFF", "#FB8500"))(100)
)

library(BiocManager)

if (BiocManager::version() != "3.20") {
  BiocManager::install(version = "3.20", update = FALSE, ask = FALSE)
}

install_load_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressMessages(library(pkg, character.only = TRUE))
    message(paste("✅", pkg, "loaded (installed)"))
    return(TRUE)
  }
  tryCatch({
    install.packages(pkg)
    suppressMessages(library(pkg, character.only = TRUE))
    message(paste("✅", pkg, "downloaded and loaded (CRAN)"))
  }, error = function(e) {
    BiocManager::install(pkg, version = "3.20", update = FALSE, quiet = TRUE, ask = FALSE)
    suppressMessages(library(pkg, character.only = TRUE))
    message(paste("✅", pkg, "downloaded and loaded (Bioconductor 3.20)"))
  })
}

required_packages <- c("dplyr", "tibble", "tidyselect", "tidyverse", "data.table", "future.apply",
                       "DESeq2", "limma", "edgeR", "pheatmap", "ggplot2", "ggrepel",
                       "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db",
                       "pathview", "enrichplot", "ggpubr", "AnnotationDbi","codetools","future")

loaded_success <- c()
for (pkg in required_packages) {
  tryCatch({
    suppressMessages(library(pkg, character.only = TRUE))
    loaded_success <- c(loaded_success, pkg)
    message(paste("✅", pkg, "explicitly loaded"))
  }, error = function(e) {
    message(paste("❌", pkg, "loading failed: ", e$message))
  })
}

message("\n📋 Successfully loaded packages:")
print(loaded_success)

library(future)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)

plan(multisession, workers = N_CORES)

species_config <- list(
  human = list(kegg_org = "hsa", go_db = org.Hs.eg.db, gene_type = "ENSEMBL"),
  mouse = list(kegg_org = "mmu", go_db = org.Mm.eg.db, gene_type = "ENSEMBL")
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("\n❌ Invalid arguments: Rscript script.R human|mouse group_file.csv\n");quit(save = "no", status = 1)
}
selected_species <- args[1]
group_csv_path <- args[2]

if (!selected_species %in% names(species_config)) {
  cat(paste0("\n❌ Invalid species: only support ", paste(names(species_config), collapse = "/"), "\n"));quit(save = "no", status = 1)
}
if (!file.exists(group_csv_path)) {
  cat(paste0("❌ File not found: ", group_csv_path, "\n"));quit(save = "no", status = 1)
}

kegg_org <- species_config[[selected_species]]$kegg_org
go_db <- species_config[[selected_species]]$go_db
gene_id_type <- species_config[[selected_species]]$gene_type

group_content <- readLines(group_csv_path)
group_content <- group_content[group_content != "" & !is.na(group_content) & !grepl("^#", group_content)]
sample_info <- as.data.frame(do.call(rbind, strsplit(group_content, ",")), stringsAsFactors = FALSE)
colnames(sample_info) <- c("sample_name", "group_name")
sample_info$sample_name <- trimws(sample_info$sample_name)
sample_info$group_name <- trimws(sample_info$group_name)
sample_info <- sample_info %>% dplyr::distinct(sample_name, .keep_all = TRUE) %>% dplyr::filter(!is.na(sample_name))
sample_info$group_name <- factor(sample_info$group_name, levels = c("control","test"))
rownames(sample_info) <- sample_info$sample_name

cat(paste0("\n✅ Run config: ", selected_species, " | ", nrow(sample_info)," samples | ", paste(levels(sample_info$group_name), collapse = " VS "),"\n"))
cat(paste0("✅ Gene ID type: ", gene_id_type, " | KEGG organism: ", kegg_org, "\n"))
cat(paste0("✅ Compute config: ", N_CORES, " threads | ", MEMORY_LIMIT, " memory limit\n\n"))

project_dir <- getwd()
count_matrix_path <- file.path(project_dir, "3.rna_seq_analysis/gene_count_matrix.tsv")
tpm_matrix_path <- file.path(project_dir, "3.rna_seq_analysis/gene_tpm_matrix.tsv")
output_dir <- file.path(project_dir, "4.rna_seq_qc_diff")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
report_path <- file.path(output_dir, "RNAseq_Analysis_Report.txt")
sink(report_path, split = TRUE)

read_expr_matrix <- function(file) {
  dt <- fread(file, header = TRUE, sep = "\t", check.names = FALSE, showProgress = FALSE, stringsAsFactors = FALSE) %>% as.data.frame()
  gene_id_col <- colnames(dt)[1]
  expr_cols <- colnames(dt)[-1]
  
  has_na <- apply(dt[, expr_cols, drop=FALSE], 1, function(x) any(is.na(x)))
  na_count <- sum(has_na)
  dt <- dt[!has_na, , drop=FALSE]
  
  dt[, expr_cols] <- lapply(dt[, expr_cols], function(x) {
    x[x == ""] <- NA
    return(x)
  })
  has_empty <- apply(dt[, expr_cols, drop=FALSE], 1, function(x) any(is.na(x)))
  empty_count <- sum(has_empty)
  dt <- dt[!has_empty, , drop=FALSE]
  
  dt[, expr_cols] <- lapply(dt[, expr_cols], as.numeric)
  has_negative <- apply(dt[, expr_cols, drop=FALSE], 1, function(x) any(x < 0, na.rm = TRUE))
  negative_count <- sum(has_negative)
  dt <- dt[!has_negative, , drop=FALSE]
  
  cat(paste0("🔍 Filter stats: NA genes=", na_count, " | empty genes=", empty_count, " | negative genes=", negative_count, "\n"))
  
  dt[, expr_cols] <- apply(dt[, expr_cols], 2, function(x) {
    x[is.na(x)] <- 0
    round(x)
  })
  
  dt <- dt %>% dplyr::group_by(!!sym(gene_id_col)) %>% dplyr::summarise_all(sum) %>% as.data.frame()
  
  rownames(dt) <- dt[[gene_id_col]]
  mat <- as.matrix(dt[, expr_cols, drop=FALSE])
  rownames(mat) <- dt[[gene_id_col]]
  return(mat)
}

group_normalize_matrix <- function(mat, sample_info) {
  mat <- mat[, rownames(sample_info), drop = FALSE]
  normalized_mat <- mat
  all_genes <- rownames(mat)
  groups <- unique(sample_info$group_name)
  
  cat("📌 Start group-wise per-gene normalization...\n")
  for (g in groups) {
    cat(paste0("\n===== Process group: ", g, " =====\n"))
    group_samples <- rownames(sample_info[sample_info$group_name == g, , drop = FALSE])
    if (length(group_samples) < 2) {
      cat(paste0("⚠️  Group ", g, " sample size <2, skip normalization\n"))
      next
    }
    cat(paste0("📌 Group ", g, " sample count: ", length(group_samples), " | genes to normalize: ", length(all_genes), "\n"))
    
    group_mat <- mat[, group_samples, drop = FALSE]
    
    for (gene in all_genes) {
      gene_vals <- group_mat[gene, ]
      gene_total <- sum(gene_vals, na.rm = TRUE)
      gene_mean_total <- gene_total / length(group_samples)
      gene_coef <- ifelse(gene_vals == 0, 1, gene_mean_total / gene_total)
      normalized_mat[gene, group_samples] <- gene_vals * gene_coef
    }
    
    cat(paste0("✅ Group ", g, " per-gene normalization done!\n"))
    test_gene <- sample(all_genes, 1)
    test_before <- round(sum(group_mat[test_gene, ], na.rm = TRUE)/length(group_samples), 4)
    test_after <- round(sum(normalized_mat[test_gene, group_samples], na.rm = TRUE)/length(group_samples), 4)
    cat(paste0("🔍 Validation gene[", test_gene, "]: group mean total (before/after) = ", test_before, " / ", test_after, "\n"))
  }
  
  q99 <- quantile(normalized_mat, 0.99, na.rm = TRUE)
  q01 <- quantile(normalized_mat, 0.01, na.rm = TRUE)
  normalized_mat[normalized_mat > q99] <- q99
  normalized_mat[normalized_mat < q01] <- q01
  cat(paste0("\n📌 Global outlier clipping at 1%/99% quantiles done\n"))
  
  return(normalized_mat)
}

tpm_matrix <- read_expr_matrix(tpm_matrix_path)
count_matrix <- read_expr_matrix(count_matrix_path)
count_matrix <- round(count_matrix)
count_matrix[count_matrix < 0] <- 0

common_samples <- intersect(colnames(count_matrix), sample_info$sample_name) %>% unique()
count_matrix <- count_matrix[, common_samples, drop=FALSE]
tpm_matrix <- tpm_matrix[, common_samples, drop=FALSE]
sample_info <- sample_info[common_samples, , drop=FALSE]

cat(paste0("   TPM matrix: ", nrow(tpm_matrix), " genes × ", ncol(tpm_matrix), " samples\n"))
cat(paste0("   Count matrix: ", nrow(count_matrix), " genes × ", ncol(count_matrix), " samples\n\n"))

cat("Part 2. PCA\n------------------------------------------------------------\n")
tpm_norm <- log2(tpm_matrix + 1)
tpm_norm[is.infinite(tpm_norm) | is.na(tpm_norm)] <- 0
pca <- prcomp(t(tpm_norm), scale. = FALSE)
pca_var <- pca$sdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 2)

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  tibble::rownames_to_column("sample_name") %>%
  dplyr::left_join(sample_info, by = "sample_name") %>%
  na.omit()

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group_name, label = sample_name)) +
  geom_point(size = 5, alpha = 0.8) + geom_text_repel(size = 4, max.overlaps = 20) +
  labs(x = paste0("PC1 (", pca_var_percent[1], "%)"), y = paste0("PC2 (", pca_var_percent[2], "%)"), title = "PCA of RNA-seq Samples", color = "Group") +
  scale_color_manual(values = color_palette$group_colors) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"))
ggsave(file.path(output_dir, "PCA_Plot.pdf"), pca_plot, width = 8, height = 6, dpi = 300)
write.table(pca_df, file.path(output_dir, "PCA_Result.tsv"), sep = "\t", row.names = F, quote = F)
cat("✅ PCA finished, results saved!\n\n")

cat("Part 3. Sample correlation heatmap\n------------------------------------------------------------\n")
sample_cor <- cor(tpm_norm, method = "pearson")
sample_cor[is.na(sample_cor)] <- 0
heatmap_plot <- pheatmap(
  sample_cor, main = "Sample Correlation Heatmap (Pearson)",
  treeheight_row=15, treeheight_col=15, show_rownames=T, show_colnames=T,
  annotation_col = sample_info[, "group_name", drop=FALSE],
  annotation_colors = list(group_name = c(control=color_palette$group_colors[2], test=color_palette$group_colors[1])),
  col = color_palette$heatmap_colors, border_color = "white", fontsize=8
)
ggsave(file.path(output_dir, "Sample_Correlation_Heatmap.pdf"), heatmap_plot, width =8, height=6, dpi=300)
cat("✅ Correlation heatmap finished, saved!\n\n")

cat("Part 4. DESeq2 DE analysis\n------------------------------------------------------------\n")
fc_cutoff <- fc_cutoff_def
padj_cutoff <- padj_cutoff_def
cat(paste0("📌 DE threshold: padj < ", padj_cutoff, " & |log2FC| ≥ ", fc_cutoff, " (", round(2^fc_cutoff,1)," fold)\n"))

keep_genes <- rowSums(count_matrix >= 100) >= (1 * ncol(count_matrix))
count_matrix_filtered <- count_matrix[keep_genes, , drop=FALSE]
cat(paste0("✅ Low-expression filter done: retain ", nrow(count_matrix_filtered), " / ", nrow(count_matrix), " genes\n"))

dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered, colData = sample_info, design = ~ group_name)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
de_result <- results(dds, contrast = c("group_name", "test", "control"))

de_result_df <- as.data.frame(de_result) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(gene_id = as.character(gene_id)) %>%
  dplyr::mutate(across(c(log2FoldChange, pvalue, padj), ~replace_na(.x, 0))) %>%
  dplyr::select(gene_id, log2FoldChange, pvalue, padj, tidyselect::everything()) %>%
  dplyr::arrange(padj)

de_result_df <- de_result_df %>%
  dplyr::mutate(gene_id_no_version = gsub("\\..*$", "", gene_id))
all_ensembl_use <- if(selected_species == "mouse") {
  de_result_df$gene_id_no_version[grepl("^ENSMUSG", de_result_df$gene_id_no_version)]
} else {
  de_result_df$gene_id_no_version[grepl("^ENSG", de_result_df$gene_id_no_version)]
}
gene_symbol_df <- tryCatch({
  bitr(all_ensembl_use, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = go_db) %>%
    dplyr::rename(gene_id_no_version = ENSEMBL, gene_symbol = SYMBOL) %>%
    dplyr::distinct(gene_id_no_version, .keep_all = TRUE)
}, error = function(e) {
  data.frame(gene_id_no_version = all_ensembl_use, gene_symbol = "NA", stringsAsFactors = F)
})
de_result_df <- dplyr::left_join(de_result_df, gene_symbol_df, by = "gene_id_no_version")
de_result_df$gene_symbol[is.na(de_result_df$gene_symbol)] <- "NA"

de_result_df <- de_result_df %>%
  dplyr::mutate(gene_type = dplyr::case_when(
    pvalue == 0 ~ "Filtered_by_DESeq2",
    pvalue < padj_cutoff & log2FoldChange >= fc_cutoff ~ "Up-regulated",
    pvalue < padj_cutoff & log2FoldChange <= -fc_cutoff ~ "Down-regulated",
    TRUE ~ "No significant"
  )) %>%
  dplyr::mutate(FoldChange = 2^log2FoldChange) %>%
  dplyr::select(gene_id, gene_symbol, log2FoldChange, FoldChange, pvalue, padj, everything())

up_gene <- nrow(dplyr::filter(de_result_df, gene_type == "Up-regulated"))
down_gene <- nrow(dplyr::filter(de_result_df, gene_type == "Down-regulated"))
cat(paste0("✅ DESeq2 DE done: up = ", up_gene, " | down = ", down_gene, " | total tested = ", nrow(de_result_df),"\n\n"))

write.table(de_result_df, file.path(output_dir, "DEG_All_Result.tsv"), sep = "\t", row.names = F, quote = F)
sig_deg <- dplyr::filter(de_result_df, gene_type %in% c("Up-regulated", "Down-regulated"))
if(nrow(sig_deg) > 0){
  write.table(sig_deg, file.path(output_dir, "DEG_Sig_List.tsv"), sep = "\t", row.names = F, quote = F)
  cat(paste0("✅ Significant DEG saved as DEG_Sig_List.tsv (", nrow(sig_deg), " genes)\n"))
}

cat("Part 5. DESeq2 volcano plot\n------------------------------------------------------------\n")
de_result_df_filtered <- de_result_df[!is.na(de_result_df$pvalue) & de_result_df$pvalue != 0, ]

volcano_plot <- ggplot(de_result_df_filtered, aes(x = log2FoldChange, y = -log10(pvalue), color = gene_type)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = color_palette$volcano_colors,
                     breaks = c("Down-regulated", "Filtered_by_DESeq2", "No significant", "Up-regulated"),
                     labels = c("Down-regulated", "Filtered", "No significant", "Up-regulated")) +
  labs(x = "log2(Fold Change)", y = "-log10(pvalue)", title = "Volcano Plot of Differentially Expressed Genes", color = "Gene Type") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"))

ggsave(file.path(output_dir, "Volcano_Plotpvalue.pdf"), volcano_plot, width =8, height=6, dpi=300)
cat("✅ DESeq2 volcano plot saved!\n\n")

cat("Part 6. edgeR DE analysis\n------------------------------------------------------------\n")
fc_cutoff <- fc_cutoff_def
padj_cutoff <- padj_cutoff_def
cat(paste0("📌 DE threshold: FDR < ", padj_cutoff, " & |log2FC| ≥ ", fc_cutoff, " (", round(2^fc_cutoff,1)," fold)\n"))

keep_genes <- rowSums(count_matrix >= 10) >= (0.75 * ncol(count_matrix))
count_matrix_filtered <- count_matrix[keep_genes, , drop=FALSE]
cat(paste0("✅ Low-expression filter done: retain ", nrow(count_matrix_filtered), " / ", nrow(count_matrix), " genes\n"))

y <- DGEList(counts = count_matrix_filtered)
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0 + group_name, data = sample_info)
colnames(design) <- levels(sample_info$group_name)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmFit(y, design)
contrasts <- makeContrasts(test - control, levels = design)
lrt <- glmLRT(fit, contrast = contrasts)
de_result <- topTags(lrt, n = Inf, adjust.method = "fdr")
de_result <- de_result$table

de_result_df <- as.data.frame(de_result) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(gene_id = as.character(gene_id)) %>%
  dplyr::rename(
    log2FoldChange = logFC,
    pvalue = PValue,
    padj = FDR
  ) %>%
  dplyr::mutate(across(c(log2FoldChange, pvalue, padj), ~replace_na(.x, 0))) %>%
  dplyr::select(gene_id, log2FoldChange, pvalue, padj, tidyselect::everything()) %>%
  dplyr::arrange(padj)

de_result_df <- de_result_df %>%
  dplyr::mutate(gene_id_no_version = gsub("\\..*$", "", gene_id))
all_ensembl_use <- if(selected_species == "mouse") {
  de_result_df$gene_id_no_version[grepl("^ENSMUSG", de_result_df$gene_id_no_version)]
} else {
  de_result_df$gene_id_no_version[grepl("^ENSG", de_result_df$gene_id_no_version)]
}
gene_symbol_df <- tryCatch({
  bitr(all_ensembl_use, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = go_db) %>%
    dplyr::rename(gene_id_no_version = ENSEMBL, gene_symbol = SYMBOL) %>%
    dplyr::distinct(gene_id_no_version, .keep_all = TRUE)
}, error = function(e) {
  data.frame(gene_id_no_version = all_ensembl_use, gene_symbol = "NA", stringsAsFactors = F)
})
de_result_df <- dplyr::left_join(de_result_df, gene_symbol_df, by = "gene_id_no_version")
de_result_df$gene_symbol[is.na(de_result_df$gene_symbol)] <- "NA"

de_result_df <- de_result_df %>%
  dplyr::mutate(gene_type = dplyr::case_when(
    pvalue < padj_cutoff & log2FoldChange >= fc_cutoff ~ "Up-regulated",
    pvalue < padj_cutoff & log2FoldChange <= -fc_cutoff ~ "Down-regulated",
    TRUE ~ "No significant"
  )) %>%
  dplyr::mutate(FoldChange = 2^log2FoldChange) %>%
  dplyr::select(gene_id, gene_symbol, log2FoldChange, FoldChange, pvalue, padj, everything())

up_gene <- nrow(dplyr::filter(de_result_df, gene_type == "Up-regulated"))
down_gene <- nrow(dplyr::filter(de_result_df, gene_type == "Down-regulated"))
cat(paste0("✅ edgeR DE done: up = ", up_gene, " | down = ", down_gene, " | total tested = ", nrow(de_result_df),"\n\n"))

write.table(de_result_df, file.path(output_dir, "edgeRDEG_All_Result.tsv"), sep = "\t", row.names = F, quote = F)
sig_deg <- dplyr::filter(de_result_df, gene_type %in% c("Up-regulated", "Down-regulated"))
if(nrow(sig_deg) > 0){
  write.table(sig_deg, file.path(output_dir, "edgeRDEG_Sig_List.tsv"), sep = "\t", row.names = F, quote = F)
  cat(paste0("✅ Significant DEG saved as DEG_Sig_List.tsv (", nrow(sig_deg), " genes)\n"))
}

cat("Part 7. edgeR volcano plot\n------------------------------------------------------------\n")
de_result_df_filtered <- de_result_df[!is.na(de_result_df$pvalue) & de_result_df$pvalue != 0, ]

volcano_plot <- ggplot(de_result_df_filtered, aes(x = log2FoldChange, y = -log10(pvalue), color = gene_type)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = color_palette$volcano_colors[ c(1,3,4) ],
                     breaks = c("Down-regulated", "No significant", "Up-regulated"),
                     labels = c("Down-regulated", "No significant", "Up-regulated")) +
  labs(x = "log2(Fold Change)", y = "-log10(pvalue)", title = "Volcano Plot of Differentially Expressed Genes (edgeR)", color = "Gene Type") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"))

ggsave(file.path(output_dir, "edgeRVolcano_Plotpvalue.pdf"), volcano_plot, width =8, height=6, dpi=300)
cat("✅ edgeR volcano plot saved!\n\n")

cat("Part 8. GO / KEGG enrichment\n------------------------------------------------------------\n")
if(nrow(sig_deg) > 0){
  sig_gene_ids <- sig_deg$gene_id %>% unique() %>% na.omit()
  cat(paste0("✅ Input DEG count: ", length(sig_gene_ids), "\n"))
  
  final_gene_ids <- sig_gene_ids %>% na.omit()
  final_keytype <- gene_id_type
  final_gene_ids <- gsub("\\..*$", "", final_gene_ids)
  if(selected_species == "mouse"){
    final_gene_ids <- final_gene_ids[grepl("^ENSMUSG", final_gene_ids)]
  }else{
    final_gene_ids <- final_gene_ids[grepl("^ENSG", final_gene_ids)]
  }
  if(selected_species == "human" && gene_id_type == "ENSEMBL" && all(grepl("^[0-9]+$", final_gene_ids))){
    final_gene_ids <- paste0("ENSG", final_gene_ids)
    cat(paste0("✅ Prepend ENSG prefix for human ENSEMBL IDs\n"))
  }
  
  gene_symbol_df_sig <- tryCatch({
    bitr(final_gene_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = go_db) %>% na.omit() %>% dplyr::distinct()
  }, error = function(e) {NULL})
  
  gene_entrez_df <- NULL
  if(!is.null(gene_symbol_df_sig) && nrow(gene_symbol_df_sig) > 0){
    gene_entrez_df <- tryCatch({
      bitr(gene_symbol_df_sig$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = go_db) %>%
        na.omit() %>% dplyr::distinct()
    }, error = function(e) {NULL})
  }
  
  go_gene_ids <- NULL
  kegg_gene_ids <- NULL
  
  if(!is.null(gene_symbol_df_sig) && nrow(gene_symbol_df_sig) > 0){
    go_gene_ids <- gene_symbol_df_sig$SYMBOL %>% unique()
    final_keytype <- "SYMBOL"
    cat(paste0("✅ ID conversion: ENSEMBL → SYMBOL (", length(go_gene_ids), " genes)\n"))
  }else{
    go_gene_ids <- final_gene_ids
    final_keytype <- "ENSEMBL"
    cat(paste0("⚠️ SYMBOL conversion failed, GO uses raw ENSEMBL (", length(go_gene_ids), " genes)\n"))
  }
  
  if(!is.null(gene_entrez_df) && nrow(gene_entrez_df) > 0){
    kegg_gene_ids <- gene_entrez_df$ENTREZID %>% unique()
    cat(paste0("✅ ID conversion: SYMBOL → EntrezID (", length(kegg_gene_ids), " genes for KEGG)\n"))
  }else{
    cat("⚠️ EntrezID conversion failed, skip KEGG\n")
  }
  
  go_enrich <- tryCatch({
    enrichGO(gene = go_gene_ids, OrgDb = go_db, keyType = final_keytype, ont = "ALL",
             pvalueCutoff = 0.2, qvalueCutoff = 0.25, readable = TRUE) %>% na.omit()
  }, error = function(e) {NULL})
  
  kegg_enrich <- NULL
  if(!is.null(kegg_gene_ids) && length(kegg_gene_ids) > 0){
    kegg_enrich <- tryCatch({
      enrichKEGG(gene = kegg_gene_ids,
                 organism = kegg_org,
                 pvalueCutoff = 0.2,
                 qvalueCutoff = 0.25,
                 use_internal_data = TRUE) %>% na.omit()
    }, error = function(e) {
      cat(paste0("❌ KEGG failed: ", e$message, "\n"))
      NULL
    })
  }
  
  if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
    go_df <- as.data.frame(go_enrich)
    go_df$p.adjust <- format(go_df$p.adjust, scientific = TRUE, digits = 2)
    go_df$qvalue <- format(go_df$qvalue, scientific = TRUE, digits = 2)
    write.table(go_df, file.path(output_dir, "GO_Enrich_Result.tsv"), sep = "\t", row.names = F, quote = F)
    cat(paste0("✅ GO enrichment done: ", nrow(go_enrich), " terms (BP+CC+MF)\n"))
  } else {
    cat("ℹ️ GO: no significant terms\n")
  }
  
  if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
    kegg_df <- as.data.frame(kegg_enrich)
    kegg_df$p.adjust <- format(kegg_df$p.adjust, scientific = TRUE, digits = 2)
    kegg_df$qvalue <- format(kegg_df$qvalue, scientific = TRUE, digits = 2)
    write.table(kegg_df, file.path(output_dir, "KEGG_Enrich_Result.tsv"), sep = "\t", row.names = F, quote = F)
    cat(paste0("✅ KEGG enrichment done: ", nrow(kegg_enrich), " pathways\n"))
  } else {
    cat("ℹ️ KEGG: no significant pathways\n")
  }
  
  cat("Part 9. GO / KEGG enrichment plots\n------------------------------------------------------------\n")
  if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
    go_bubble <- dotplot(go_enrich, showCategory = 10, split="ONTOLOGY") +
      facet_grid(ONTOLOGY~., scale='free') +
      scale_color_gradient(low = color_palette$group_colors[2],
                           high = color_palette$group_colors[1],
                           trans = "reverse",
                           labels = scientific_format(digits=2)) +
      labs(title = "GO Enrichment Bubble Plot (TOP10 per ONTOLOGY)", color = "p.adjust") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=12, face="bold"), axis.text.y = element_text(size=8))
    ggsave(file.path(output_dir, "GO_Enrich_Bubble.pdf"), go_bubble, width = 10, height = 8, dpi = 300)
    cat("✅ GO bubble plot saved\n")
  }
  if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
    kegg_bubble <- dotplot(kegg_enrich, showCategory = 15) +
      scale_color_gradient(low = color_palette$group_colors[2],
                           high = color_palette$group_colors[1],
                           trans = "reverse",
                           labels = scientific_format(digits=2)) +
      labs(title = "KEGG Pathway Enrichment Bubble Plot (TOP15)", color = "p.adjust") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=12, face="bold"), axis.text.y = element_text(size=8))
    ggsave(file.path(output_dir, "KEGG_Enrich_Bubble.pdf"), kegg_bubble, width = 9, height = 7, dpi = 300)
    kegg_bar <- barplot(kegg_enrich, showCategory = 15, fill = "p.adjust") +
      scale_fill_gradient(low = color_palette$group_colors[2],
                         high = color_palette$group_colors[1],
                         trans = "reverse",
                         labels = scientific_format(digits=2)) +
      labs(title = "KEGG Pathway Enrichment Bar Plot (TOP15)", fill = "p.adjust") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=12, face="bold"), axis.text.y = element_text(size=8))
    ggsave(file.path(output_dir, "KEGG_Enrich_Bar.pdf"), kegg_bar, width = 9, height = 7, dpi = 300)
    cat("✅ KEGG bubble + bar plot saved\n")
  }
  
} else {
  cat("⚠️ No significant DEG, skip enrichment\n")
}

sink()
plan(sequential)
options(warn=0)
cat("\n============================================================\n")
cat("🎉 RNA-seq / ribo-seq pipeline finished!\n")
cat(paste0("📁 Output directory: ", output_dir, "\n"))
cat(paste0("📊 Summary: up-regulated ", up_gene, " | down-regulated ", down_gene, "\n"))
cat("============================================================\n")
