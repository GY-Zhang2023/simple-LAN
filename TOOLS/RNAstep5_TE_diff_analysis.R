#!/usr/bin/env Rscript
Sys.setenv(LANGUAGE='en')
options(stringsAsFactors = FALSE, warn=-1, datatable.verbose = FALSE)
options(digits = 15)

args <- commandArgs(trailingOnly = TRUE)
te_matrix_path <- NULL
group_csv_path <- NULL

if (length(args) >= 2 && file.exists(args[1]) && file.exists(args[2])) {
  te_matrix_path <- args[1]
  group_csv_path <- args[2]
} else if (interactive()) {
  cat("TE matrix TSV path: ");te_matrix_path=readLines(n=1);te_matrix_path=gsub("^\\s+|\\s+$","",te_matrix_path)
  cat("Group CSV path: ");group_csv_path=readLines(n=1);group_csv_path=gsub("^\\s+|\\s+$","",group_csv_path)
} else {
  cat("Usage: Rscript TE_diff_final_correct.R TE_data.tsv group_info.csv\n");quit(save="no",status=1)
}

OUTPUT_ROOT <- file.path(dirname(te_matrix_path), "TE_diff_result_correct")
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c("data.table", "tidyverse", "ggpubr", "pheatmap", "ggrepel", "qvalue", "car", "limma")
invisible(sapply(required_pkgs, function(pkg) {
  if(!requireNamespace(pkg,quietly=T)) install.packages(pkg,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(pkg,character.only=T,quietly=T,warn.conflicts=F)
}))

process_te <- function(df, low_expr_threshold = 0.1) {
  df <- df %>%
    dplyr::mutate(TE = ifelse(rna_raw_count == 0 & TE == 1e6, NA, TE)) %>%
    dplyr::mutate(TE = ifelse(TE < 1e-4 | is.infinite(TE) | TE > 1e8, NA, TE)) %>%
    dplyr::mutate(TE = ifelse(is.na(TE), NA, as.numeric(TE))) %>%
    group_by(gene_name) %>%
    mutate(mean_te = mean(TE, na.rm = TRUE)) %>%
    filter(mean_te > low_expr_threshold) %>%
    ungroup() %>%
    select(-mean_te)
  return(df)
}

limma_te_analysis <- function(te_matrix, ctrl_samples, test_samples) {
  te_quantile <- quantile(te_matrix, c(0.01, 0.99), na.rm = TRUE)
  te_matrix_clip <- te_matrix
  te_matrix_clip[te_matrix_clip < te_quantile[1]] <- te_quantile[1]
  te_matrix_clip[te_matrix_clip > te_quantile[2]] <- te_quantile[2]
  cat(paste0("Trimming thresholds (for statistics only): ", round(te_quantile[1], 4), " ~ ", round(te_quantile[2], 4), "\n"))

  te_matrix_log <- log2(te_matrix_clip + 1e-30)

  gene_cv <- apply(te_matrix_log, 1, function(x) sd(x)/mean(x))
  te_matrix_log_filter <- te_matrix_log[!is.infinite(gene_cv) & gene_cv <= 5, ]
  cat(paste0("CV filter: retain ", nrow(te_matrix_log_filter), "/", nrow(te_matrix_log), " genes\n"))

  group <- factor(c(rep("Ctrl", length(ctrl_samples)), rep("Test", length(test_samples))))
  design <- model.matrix(~0 + group)
  colnames(design) <- c("Ctrl", "Test")
  rownames(design) <- c(ctrl_samples, test_samples)

  contrast.matrix <- makeContrasts(Test - Ctrl, levels = design)

  fit <- lmFit(te_matrix_log_filter, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)

  limma_res <- topTable(fit2, coef=1, adjust="fdr", number=Inf) %>%
    rownames_to_column("gene_name") %>%
    rename(
      log2FC_limma = logFC,
      pvalue = P.Value,
      padj_fdr_bh = adj.P.Val,
      t_stat = t,
      B_stat = B
    ) %>%
    rowwise() %>%
    mutate(
      mean_ctrl_raw = round(mean(te_matrix[gene_name, ctrl_samples]), 4),
      mean_test_raw = round(mean(te_matrix[gene_name, test_samples]), 4),
      fc_raw = round(mean_test_raw / mean_ctrl_raw, 4),
      log2fc_raw = round(log2(fc_raw + 1e-30), 4),
      mean_ctrl_clip = round(2^mean(te_matrix_log_filter[gene_name, ctrl_samples]), 4),
      mean_test_clip = round(2^mean(te_matrix_log_filter[gene_name, test_samples]), 4)
    ) %>%
    ungroup() %>%
    select(gene_name,
           mean_ctrl_raw, mean_test_raw, fc_raw, log2fc_raw,
           mean_ctrl_clip, mean_test_clip, log2FC_limma,
           t_stat, pvalue, padj_fdr_bh, B_stat)

  limma_res <- limma_res %>%
    left_join(
      adjust_pvalues_simple(limma_res$pvalue, limma_res$gene_name) %>% select(-pvalue),
      by = "gene_name"
    )

  return(limma_res)
}

adjust_pvalues_simple <- function(p_values, gene_names) {
  p_values[is.na(p_values)] <- 1
  p_values[p_values < 1e-20] <- 1e-20

  res_df <- tibble(
    gene_name = gene_names,
    pvalue = p_values,
    padj_qvalue = NA_real_
  )

  res_df$padj_qvalue <- tryCatch({
    q_obj <- qvalue::qvalue(p = p_values)
    q_obj$qvalues
  }, error = function(e) {
    message("qvalue failed, fallback to BH: ", e$message)
    p.adjust(p_values, method = "fdr")
  })

  res_df[is.na(res_df)] <- 1
  return(res_df)
}

format_numeric <- function(df, cols) {
  for (col in cols) {
    df[[col]] <- as.numeric(df[[col]])
    df[[col]] <- round(df[[col]], 15)
  }
  return(df)
}

cat("\n=======================\nStep 1: Load data + preprocessing\n=======================\n")
group_info <- fread(group_csv_path, header = T, sep = ",", check.names = F) %>%
  distinct(sample_name, .keep_all = T) %>%
  filter(!is.na(group_name))

group_names <- unique(group_info$group_name)
if (length(group_names) != 2) {
  stop(paste0("Group file must contain exactly 2 groups! Current: ", paste(group_names, collapse = ","), " (n=", length(group_names), ")"))
}
ctrl_group <- group_names[1]
test_group <- group_names[2]
ctrl_samples <- group_info$sample_name[group_info$group_name == ctrl_group]
test_samples <- group_info$sample_name[group_info$group_name == test_group]
target_samples <- group_info$sample_name

cat(paste0("Group info (NA filtered):\n"))
cat(paste0("  -> Control(", ctrl_group, "): ", paste(ctrl_samples, collapse = ","),"\n"))
cat(paste0("  -> Test(", test_group, "): ", paste(test_samples, collapse = ","),"\n"))

te_all <- fread(te_matrix_path, header = T, sep = "\t", check.names = F) %>% as.data.frame()
te_all <- te_all %>% mutate(gene_name = as.character(gene_name))
te_filtered <- process_te(te_all, low_expr_threshold = 0.1)

cat(paste0("Preprocessing done:\n"))
cat(paste0("  -> Raw genes: ", length(unique(te_all$gene_name)),"\n"))
cat(paste0("  -> Filtered genes: ", length(unique(te_filtered$gene_name)),"\n"))
cat(paste0("  -> NA TE entries: ", sum(is.na(te_filtered$TE)),"\n"))

cat("\n=======================\nStep 2: Build wide matrix\n=======================\n")
if (!all(target_samples %in% unique(te_filtered$sample_name))) {
  missing <- setdiff(target_samples, unique(te_filtered$sample_name))
  cat(paste0("Sample mismatch, missing: ", paste(missing, collapse = ","), ", auto filtered\n"))
  target_samples <- intersect(target_samples, unique(te_filtered$sample_name))
}

te_wide <- te_filtered %>%
  select(gene_name, sample_name, TE) %>%
  pivot_wider(names_from = sample_name, values_from = TE) %>%
  select(gene_name, all_of(target_samples))

te_wide_filter <- te_wide %>%
  rowwise() %>%
  mutate(na_count = sum(is.na(c_across(all_of(target_samples))))) %>%
  filter(na_count == 0) %>%
  select(-na_count) %>%
  ungroup()

te_matrix <- te_wide_filter %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
mode(te_matrix) <- "numeric"

valid_genes <- rownames(te_matrix)
cat(paste0("Wide matrix ready:\n"))
cat(paste0("  -> Dimension: ", nrow(te_matrix), " × ", ncol(te_matrix), "\n"))
cat(paste0("  -> Valid genes: ", length(valid_genes),"\n"))
if (length(valid_genes) == 0) stop("No valid genes left!")

write.table(te_matrix, file.path(OUTPUT_ROOT, "TE_wide_matrix.tsv"), sep="\t", quote=F, row.names=T)

cat("\n=======================\nStep 3: Limma differential analysis\n=======================\n")
diff_result_raw <- limma_te_analysis(te_matrix, ctrl_samples, test_samples)

diff_result <- diff_result_raw %>%
  format_numeric(cols = c("pvalue", "padj_fdr_bh", "padj_qvalue",
                          "fc_raw", "log2fc_raw", "log2FC_limma",
                          "mean_ctrl_raw", "mean_test_raw", "mean_ctrl_clip", "mean_test_clip")) %>%
  mutate(
    diff_type = case_when(
      padj_fdr_bh < 0.1 & log2fc_raw > 1  ~ "Up-regulated",
      padj_fdr_bh < 0.1 & log2fc_raw < -1 ~ "Down-regulated",
      TRUE ~ "No significant"
    ),
    diff_type_qvalue = case_when(
      padj_qvalue < 0.1 & log2fc_raw > 1  ~ "Up-regulated",
      padj_qvalue < 0.1 & log2fc_raw < -1 ~ "Down-regulated",
      TRUE ~ "No significant"
    )
  ) %>%
  arrange(padj_fdr_bh, abs(log2fc_raw), desc(log2fc_raw))

up_bh <- nrow(filter(diff_result, diff_type == "Up-regulated"))
down_bh <- nrow(filter(diff_result, diff_type == "Down-regulated"))
up_q <- nrow(filter(diff_result, diff_type_qvalue == "Up-regulated"))
down_q <- nrow(filter(diff_result, diff_type_qvalue == "Down-regulated"))

cat(paste0("Diff summary (threshold: padj<0.1, |raw TE log2FC|>1):\n"))
cat(paste0("  -> FDR(BH): Up ", up_bh, " | Down ", down_bh, "\n"))
cat(paste0("  -> qvalue: Up ", up_q, " | Down ", down_q, "\n"))

mean_ctrl_all_raw <- round(mean(apply(te_matrix[, ctrl_samples], 2, mean)), 4)
mean_test_all_raw <- round(mean(apply(te_matrix[, test_samples], 2, mean)), 4)
cat(paste0("Global raw TE summary:\n"))
cat(paste0("  -> Control mean TE: ", mean_ctrl_all_raw, "\n"))
cat(paste0("  -> Test mean TE: ", mean_test_all_raw, "\n"))
cat(paste0("  -> Test TE is ", ifelse(mean_test_all_raw < mean_ctrl_all_raw, "lower", "higher"), " than Control\n"))

if (down_bh > 0) {
  top_down <- diff_result %>% filter(diff_type=="Down-regulated") %>% arrange(padj_fdr_bh) %>% head(3)
  cat(paste0("  -> Top down genes (raw TE):\n"))
  for (i in 1:nrow(top_down)) {
    cat(paste0("     - ", top_down$gene_name[i], ": FC=", top_down$fc_raw[i], " | log2FC=", top_down$log2fc_raw[i], " | FDR(BH)=", top_down$padj_fdr_bh[i], "\n"))
  }
} else if (up_bh > 0) {
  top_up <- diff_result %>% filter(diff_type=="Up-regulated") %>% arrange(padj_fdr_bh) %>% head(3)
  cat(paste0("  -> Top up genes (raw TE):\n"))
  for (i in 1:nrow(top_up)) {
    cat(paste0("     - ", top_up$gene_name[i], ": FC=", top_up$fc_raw[i], " | log2FC=", top_up$log2fc_raw[i], " | FDR(BH)=", top_up$padj_fdr_bh[i], "\n"))
  }
} else {
  cat("  -> No significant genes found. Try looser threshold: padj<0.6 + |raw TE log2FC|>0.2\n")
}

cat("\n=======================\nStep 4: Save outputs\n=======================\n")
anno_info <- te_filtered %>%
  select(gene_name, sample_name, rna_raw_count, rna_normalized, ribo_raw_rpfs, ribo_rpm) %>%
  group_by(gene_name, sample_name) %>% slice(1) %>% ungroup()

diff_result_full <- diff_result %>%
  left_join(anno_info, by = "gene_name") %>%
  select(gene_name,
         mean_ctrl_raw, mean_test_raw, fc_raw, log2fc_raw,
         pvalue, padj_fdr_bh, padj_qvalue,
         diff_type, diff_type_qvalue,
         sample_name, rna_raw_count, rna_normalized, ribo_raw_rpfs, ribo_rpm) %>%
  arrange(padj_fdr_bh, abs(log2fc_raw), desc(log2fc_raw))

write.table(diff_result_full,
            file.path(OUTPUT_ROOT, "diff_TE_raw_FC_full_anno.tsv"),
            sep="\t", quote=F, row.names=F, na="", dec=".")

write.table(filter(diff_result, diff_type!="No significant") %>%
              select(gene_name, mean_ctrl_raw, mean_test_raw, fc_raw, log2fc_raw, pvalue, padj_fdr_bh, diff_type),
            file.path(OUTPUT_ROOT, "diff_TE_raw_FC_fdr_bh.tsv"),
            sep="\t", quote=F, row.names=F, na="", dec=".")

write.table(filter(diff_result, diff_type_qvalue!="No significant") %>%
              select(gene_name, mean_ctrl_raw, mean_test_raw, fc_raw, log2fc_raw, pvalue, padj_qvalue, diff_type_qvalue),
            file.path(OUTPUT_ROOT, "diff_TE_raw_FC_qvalue.tsv"),
            sep="\t", quote=F, row.names=F, na="", dec=".")

summary_table <- data.frame(
  run_time = Sys.time(),
  raw_gene_count = length(unique(te_all$gene_name)),
  filtered_gene_count = length(unique(te_filtered$gene_name)),
  valid_gene_count = length(valid_genes),
  ctrl_sample_count = length(ctrl_samples),
  test_sample_count = length(test_samples),
  fdr_bh_up = up_bh,
  fdr_bh_down = down_bh,
  qvalue_up = up_q,
  qvalue_down = down_q,
  ctrl_mean_te_raw = mean_ctrl_all_raw,
  test_mean_te_raw = mean_test_all_raw,
  global_trend = ifelse(mean_test_all_raw < mean_ctrl_all_raw, "Test<Control", "Test>Control"),
  threshold = "padj < 0.1 + |raw TE log2FC| > 1",
  test_method = "limma (trend=TRUE+robust=TRUE)",
  adjust_method = "limma FDR(BH) + qvalue",
  preprocessing = "quantile trimming + CV filter (statistics only)",
  metric = "raw TE fold change + adjusted p-value"
)
write.table(summary_table, file.path(OUTPUT_ROOT, "diff_TE_raw_FC_summary.tsv"), sep="\t", quote=F, row.names=F, na="", dec=".")

cat("\n=======================\nStep 5: Visualization\n=======================\n")
plot_theme <- theme_bw() + theme(
  plot.title = element_text(hjust=0.5, size=14, face="bold"),
  axis.title = element_text(size=12, face="bold"),
  panel.grid = element_blank(),
  legend.position = "top"
)

plot_data <- te_filtered %>%
  left_join(group_info, by = "sample_name") %>%
  filter(!is.na(TE)) %>%
  filter(!is.na(group_name)) %>%
  filter(group_name %in% c(ctrl_group, test_group))

cat(paste0("Boxplot data stats (NA removed):\n"))
cat(paste0("  -> Total rows: ", nrow(plot_data),"\n"))
cat(paste0("  -> Group distribution:\n"))
for (g in group_names) {
  cat(paste0("     - ", g, ": ", nrow(filter(plot_data, group_name == g)), " rows\n"))
}

p1 <- ggboxplot(plot_data, x="group_name", y="TE", fill="group_name",
                palette = c("#2E86AB", "#E74C3C"), alpha=0.7, add="jitter", add.params=list(size=0.3)) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x=1.5) +
  labs(x="Group", y="Raw TE Value", title="Raw TE Expression Distribution (No NA Groups)") +
  plot_theme + theme(legend.position="none")
ggsave(file.path(OUTPUT_ROOT, "TE_raw_boxplot_no_na.pdf"), p1, width=6, height=5, dpi=300)

volcano_data <- diff_result %>%
  mutate(
    color = case_when(
      diff_type == "Up-regulated" ~ "Up",
      diff_type == "Down-regulated" ~ "Down",
      TRUE ~ "NS"
    ),
    y = -log10(padj_fdr_bh)
  ) %>%
  filter(!is.na(log2fc_raw) & !is.na(y))

top20_genes <- volcano_data %>%
  filter(diff_type != "No significant") %>%
  arrange(padj_fdr_bh, desc(abs(log2fc_raw))) %>%
  head(20) %>%
  mutate(label = gene_name)

x_max_abs <- ceiling(max(abs(volcano_data$log2fc_raw), na.rm = TRUE) + 0.5)
x_min <- -x_max_abs
x_max <- x_max_abs

p2 <- ggplot(volcano_data, aes(x=log2fc_raw, y=y, color=color)) +
  geom_point(alpha=0.8, stroke=0.1, size=1) +
  geom_text_repel(
    data = top20_genes,
    aes(label = label),
    size = 3.5,
    color = "black",
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "gray50",
    segment.size = 0.5,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="gray50", linewidth=0.8) +
  geom_hline(yintercept = -log10(0.1), linetype="dashed", color="gray50", linewidth=0.8) +
  scale_color_manual(values = c("#27AE60", "#95A5A6", "#E74C3C"),
                     labels = c("Down-regulated", "No significant", "Up-regulated")) +
  xlim(c(x_min, x_max)) +
  labs(x="log2(Fold Change) (Raw TE)", y="-log10(FDR-BH adjusted p-value)",
       title="TE Volcano Plot (Raw FC + Adjusted p-value) - Top 20 Genes Labeled") +
  plot_theme +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=11),
    axis.title = element_text(size=13)
  )

ggsave(file.path(OUTPUT_ROOT, "TE_raw_FC_volcano_top20.pdf"), p2, width=12, height=7, dpi=300)

total_diff <- up_bh + down_bh
if (total_diff >=5) {
  te_matrix_log_raw <- log2(te_matrix + 1e-30)

  top_genes <- diff_result %>% filter(diff_type!="No significant") %>% arrange(padj_fdr_bh) %>% head(min(50, total_diff)) %>% pull(gene_name)

  plot_mat <- te_matrix_log_raw[top_genes, ]

  row_normalize <- function(mat) {
    row_means <- apply(mat, 1, mean, na.rm = TRUE)
    row_sds <- apply(mat, 1, sd, na.rm = TRUE)
    row_sds[row_sds == 0] <- 1
    norm_mat <- (mat - row_means) / row_sds
    return(norm_mat)
  }

  heatmap_mat <- row_normalize(plot_mat)

  heatmap_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

  anno_col <- data.frame(Group = group_info$group_name)
  rownames(anno_col) <- group_info$sample_name
  anno_colors <- list(Group = setNames(c("#2E86AB", "#E74C3C"), group_names))

  pheatmap(
    mat = heatmap_mat,
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    na_col = "white",
    color = heatmap_colors,
    border_color = NA,
    main = paste0("Top ", min(50, total_diff), " Differential TE Genes (Row Normalized Z-score)"),
    filename = file.path(OUTPUT_ROOT, "TE_raw_FC_heatmap_row_normalized.pdf"),
    width = 8,
    height = 10,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    scale = "none"
  )

  cat(paste0("Heatmap row Z-score ready: ", nrow(heatmap_mat), " genes, ", ncol(heatmap_mat), " samples\n"))
  cat(paste0("Z-score range: ", round(min(heatmap_mat), 2), " ~ ", round(max(heatmap_mat), 2), "\n"))
}

cat("\nPipeline finished!\n")
cat("Output directory: ", OUTPUT_ROOT, "\n")
cat("Main outputs:\n")
cat("  - diff_TE_raw_FC_full_anno.tsv\n")
cat("  - TE_raw_boxplot_no_na.pdf\n")
cat("  - TE_raw_FC_volcano_top20.pdf\n")
cat("  - TE_raw_FC_heatmap_row_normalized.pdf\n")
