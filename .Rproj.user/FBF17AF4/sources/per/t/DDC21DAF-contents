# packages and functions ----
rm(list = ls())
source("./code/utils.R")
library(sva)
library(preprocessCore)
## analysis init settings ----
NumberPSM_threshold <- 2
missing_threshold <- 0.5
fill_na_value <- 0.8
output_path <- "./output/step1_QC_of_data/"
# decide whether to use quantiles to normalize the data
use_quantiles <- FALSE
# decide whether to perform batch correction
perform_batch_correction <- F
bad_samples <- c("b11_13", "b18_13", "b17_12", "b20_11", "b3_14", "b19_3")
# bad_samples <- c()


# read data and preprocess ----
raw_abundance <- read.csv("./data/samples/Protein_exp/abundance/abundance_protein_None.tsv", sep = "\t") %>%
  .[.$NumberPSM >= NumberPSM_threshold, ]
# generate index to gene symbol
if (!file.exists("./data/index_to_gene_symbol.rds")) {
  index_to_gene_symbol <- raw_abundance[, c(1, 3)]
  saveRDS(index_to_gene_symbol, "./data/index_to_gene_symbol.rds")
  index_symbol_mapping <- readRDS("./data/index_to_gene_symbol.rds")
}
## separate the pool samples and normal(batch) samples ----
separate_masks <- separater_pool_batch(colnames(raw_abundance))
pool_abundance <- cbind(raw_abundance[, c(1:5)], raw_abundance[, separate_masks$pool])
batch_abundance <- cbind(raw_abundance[, c(1:5)], raw_abundance[, separate_masks$batch])
# if there are bad samples, remove them
if (length(bad_samples) > 0) {
  pool_abundance <- pool_abundance[, !colnames(pool_abundance) %in% bad_samples]
  batch_abundance <- batch_abundance[, !colnames(batch_abundance) %in% bad_samples]
}
## remove and fill the missing values ----
remaining_cols <- c(1:5)
pool_abundance_filled <- remove_fill_na(remaining_cols, pool_abundance, missing_threshold, fill_na_value)
batch_abundance_filled <- remove_fill_na(remaining_cols, batch_abundance, missing_threshold, fill_na_value)
# here we get the pool_abundance_filled and batch_abundance_filled which are the data we need for the following analysis if we want to
# use quantiles to normalize the data
if (use_quantiles) {
  pool_abundance_filled[, -remaining_cols] <- normalize.quantiles(as.matrix(pool_abundance_filled[, -remaining_cols]))
  batch_abundance_filled[, -remaining_cols] <- normalize.quantiles(as.matrix(batch_abundance_filled[, -remaining_cols]))
}


# QC ---- Next we perform the QC of the data QC: CV calculate ---- CV among pool samples(CV_pool) should be smaller than it in batch
# samples(CV_bx) Don't use filled data for CV calculation
CV_pool <- calculate_CV(pool_abundance[, -remaining_cols], verbose = T)
CV_bx <- calculate_CV(batch_abundance[, -remaining_cols], verbose = T)
### CV plot and save ----
df <- data.frame(CV = c(CV_pool, CV_bx), type = c(rep("pool", length(CV_pool)), rep("batch", length(CV_bx))))
df$type <- factor(df$type)
stat.test.CV <- df %>%
  wilcox_test(CV ~ type, p.adjust.method = "BH", detailed = T)
# final plot
colors <- c("#eb4b3a", "#48bad0")
p <- ggplot(df) +
  geom_boxplot(aes(type, CV, color = type), outlier.shape = 21) +
  scale_color_manual(values = c("#eb4b3a", "#48bad0")) +
  scale_x_discrete(labels = c("Batch samples", "Pool samples")) +
  xlab("") +
  ylab("CV(%)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.8, face = "bold", color = colors, size = 11, hjust = 0.5), axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 21)
  ) +
  scale_y_continuous(breaks = c(
    0, 10, 15,
    20, 25, 30, 40
  ), limits = c(0, 45))
x_value <- c(2)
y_value_data <- unlist(lapply(levels(df$type), function(x) {
  return(max(df[df$type == x, 1], na.rm = T))
}))[2]
names(y_value_data) <- levels(df$group)[2]
y_value <- rep(y_value_data, 1) + 0.01
y_value <- y_value + c(2.5)
color_value <- c(colors[1])
for (i in 1:nrow(stat.test.CV)) {
  if (stat.test.CV$p[i] != "ns") {
    y_tmp <- y_value[i]
    p <- p + annotate(geom = "text", label = stat.test.CV$p[i], x = x_value[i], y = y_tmp, color = color_value[i], size = 4.5)
  }
}
save_folder <- paste(output_path, "/", sep = "")
if (!dir.exists(save_folder)) {
  dir.create(save_folder)
}
ggsave(paste(output_path, "/", "CV_boxplot.tiff", sep = ""), height = 10, width = 10, dpi = 600, units = "cm")

## QC: corr calculate ---- Correlation between pool samples Don't use unfilled data for corr calculation
cor_res <- rcorr(as.matrix(pool_abundance_filled[, -remaining_cols][, paste0("pool", 1:22)]), type = "spearman")
cor_mat <- cor_res$r
p_mat <- cor_res$P
### corr dot plot and save ----
cor_mat_up <- cor_mat
# cor_mat_up[lower.tri(cor_mat_up, diag = T)] <- NA
cor_mat_up_long <- cor_mat_up %>%
  as.data.frame() %>%
  mutate(x = factor(rownames(cor_mat_up), levels = rev(paste0("pool", 1:22)))) %>%
  tidyr::pivot_longer(cols = !x, names_to = "y", values_to = "cor") %>%
  mutate(y = factor(y, levels = paste0("pool", 1:22)))
ggplot() +
  geom_point(data = cor_mat_up_long, aes(x, y, fill = cor), shape = 21, color = "#c4bcba", size = 5) +
  scale_fill_gradient2(
    low = "#4178a7",
    mid = "#ffffff", high = "#e2201c", midpoint = 0.92
  ) +
  scale_x_discrete(expand = c(0.025, 0.025), position = "top") +
  scale_y_discrete(expand = c(
    0.025,
    0.025
  ), position = "right") +
  geom_vline(aes(xintercept = seq(0.5, 22.5, 1)), color = "#bbbbbb") +
  geom_hline(aes(yintercept = seq(
    0.5,
    22.5, 1
  )), color = "#bbbbbb") +
  theme_classic() +
  theme(
    panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 15), axis.text.y = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 15)
  )
save_folder <- paste(output_path, "/", sep = "")
if (!dir.exists(save_folder)) {
  dir.create(save_folder)
}
ggsave(paste(output_path, "/", "corr_dot_plot.tiff", sep = ""), width = 22, height = 20, units = "cm", dpi = 600)

## QC: batch effect ----
raw_ratio <- read.csv("./data/samples/Protein_exp/ratio/ratio_protein_MD.tsv", sep = "\t") %>%
  .[.$NumberPSM >= NumberPSM_threshold, ]
# if there are bad samples, remove them
if (length(bad_samples) > 0) {
  raw_ratio <- raw_ratio[, !colnames(raw_ratio) %in% bad_samples]
}
# remove and fill the missing values
remaining_cols <- c(1:5)
ratio_filled <- remove_fill_na(remaining_cols, raw_ratio, missing_threshold, fill_na_value) %>%
  .[, -c(2, 4, 5)]
rownames(ratio_filled) <- ratio_filled$Index
ratio_filled <- ratio_filled[, -c(1, 2)]
# read group information and remove samples with missing group information
group_info <- readRDS("./data/samples/info/group_info_of_samples.rds")
ratio_filled <- preprocess_ratio_data(ratio_filled, group_info)
ratio_filled$batch <- factor(ratio_filled$batch, levels = paste0("b", 1:22))
# if we use quantiles to normalize the data
if (use_quantiles) {
  ratio_filled[, 1:(ncol(ratio_filled) - 2)] <- normalize.quantiles(t(as.matrix(t(ratio_filled[, 1:(ncol(ratio_filled) - 2)]))))
}
### use barplot to show the batch effect ----
colors <- get_color(seed = 1226, n = 22)
names(colors) <- levels(ratio_filled$batch)
for (type in c("ov", "cys", "ctrl")) {
  msk <- unlist(lapply(ratio_filled$group, function(x) {
    return(grepl(type, x))
  }))
  tmp_data <- as.data.frame(t(ratio_filled[msk, 1:(ncol(ratio_filled) - 2)]))
  tmp_data <- pivot_longer(data = tmp_data, cols = colnames(tmp_data), names_to = "sample", values_to = "feature_value")
  tmp_data$batch <- unlist(lapply(tmp_data$sample, function(x) {
    return(str_split_1(x, "_")[1])
  }))
  tmp_data$batch <- factor(tmp_data$batch, levels = paste0("b", 1:22))
  tmp_data$color <- unlist(lapply(tmp_data$batch, function(x) {
    return(colors[[x]])
  }))
  p <- ggplot(data = tmp_data, aes(x = sample, y = feature_value, color = batch)) +
    geom_boxplot() +
    scale_color_manual(values = colors) +
    xlab("") +
    ylab("Feature") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(), axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 26),
      plot.title = element_text(hjust = 0.5, size = 31)
    ) +
    ggtitle(type)

  save_folder <- paste(output_path, "/", sep = "")
  if (!dir.exists(save_folder)) {
    dir.create(save_folder)
  }
  ggsave(paste(output_path, "/", "batch_box_plot_", type, ".tiff", sep = ""), width = 40, height = 25, units = "cm", dpi = 600)
}
### use heatmap to show the batch effect ---- TO BE IMPLEMENTED batch effect removal (based combat algorithm) ----
mat <- t(as.matrix(ratio_filled[, 1:(dim(ratio_filled)[2] - 2)]))
batch <- as.factor(ratio_filled$batch)
if (perform_batch_correction) {
  mat <- ComBat(mat, batch = batch)
} else {
  message("Skip batch correction")
}
ratio_corrected <- ratio_filled
ratio_corrected[, 1:(dim(ratio_corrected)[2] - 2)] <- t(mat)
### use barplot to show the batch effect after correction ----
if (perform_batch_correction) {
  for (type in c("ov", "cys", "ctrl")) {
    msk <- unlist(lapply(ratio_corrected$group, function(x) {
      return(grepl(type, x))
    }))
    tmp_data <- as.data.frame(t(ratio_corrected[msk, 1:(ncol(ratio_corrected) - 2)]))
    tmp_data <- pivot_longer(data = tmp_data, cols = colnames(tmp_data), names_to = "sample", values_to = "feature_value")
    tmp_data$batch <- unlist(lapply(tmp_data$sample, function(x) {
      return(str_split_1(x, "_")[1])
    }))
    tmp_data$batch <- factor(tmp_data$batch, levels = paste0("b", 1:22))
    tmp_data$color <- unlist(lapply(tmp_data$batch, function(x) {
      return(colors[[x]])
    }))
    p <- ggplot(data = tmp_data, aes(x = sample, y = feature_value, color = batch)) +
      geom_boxplot() +
      scale_color_manual(values = colors) +
      xlab("") +
      ylab("Feature") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(), axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 26),
        plot.title = element_text(hjust = 0.5, size = 31)
      ) +
      ggtitle(type)

    save_folder <- paste(output_path, "/", sep = "")
    if (!dir.exists(save_folder)) {
      dir.create(save_folder)
    }
    ggsave(paste(output_path, "/", "corrected_box_plot_", type, ".tiff", sep = ""), width = 40, height = 25, units = "cm", dpi = 600)
  }
}
### save the corrected data ----
saveRDS(ratio_corrected, "./data/samples/Protein_exp/ratio/ratio_after_QC_corrected.rds")


# checkpoint ----
save.image("./data/analysis_check_point/step1_QC_of_data.RData")
rm(list = ls())
