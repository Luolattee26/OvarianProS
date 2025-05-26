# packages and functions ----
rm(list = ls())
library(readxl)
library(RColorBrewer)
source("./code/utils.R")
## analysis init settings ----
fill_na_value <- 0.8
adjust_method <- "BH"
pvalue_threshold <- 0.05
logfc_threshold <- 0
output_path <- "./output/step4_PRM_eval/"
DEPs_path <- "./output/step2_DEPs/"
input_exp_file <- "./data/samples/Protein_exp/ratio/ratio_after_QC_corrected.rds"
PRM_data <- "./data/samples/Protein_exp/PRM_res_24.12.12/wosp24127_prot_matrix.csv"
PRM_sample_info <- "./data/samples/info/sample_info_evalCohort.xlsx"
remove_list <- c()
remove_list <- c("CA3", "CA22", "CA23", "CA1", "CA47", "CA53", "CA56", "CA66")
if (!dir.exists(output_path)) {
  dir.create(output_path)
}
index_symbol_mapping <- readRDS("./data/index_to_gene_symbol.rds")


# read data and preprocess ----
PRM_data <- read.csv(PRM_data, header = T, row.names = 1)
PRM_sample_info <- read_excel(PRM_sample_info)
s1_list <- PRM_sample_info$SampleID...2[grepl("T", PRM_sample_info$S1)]
# use a for loop to find the columns with 'pool' in the name
batch_samples <- vector(mode = "list")
pool_cols <- vector(mode = "list")
tmp_cols <- c()
for (colname in colnames(PRM_data)) {
  if (grepl("pool", colname)) {
    pool_cols <- c(pool_cols, colname)
    batch_samples <- c(batch_samples, list(tmp_cols))
    tmp_cols <- c()
  } else {
    tmp_cols <- c(tmp_cols, colname)
  }
}
# get ratio (in log2 scale)
for (i in 1:length(pool_cols)) {
  PRM_data[, batch_samples[[i]]] <- PRM_data[, batch_samples[[i]]] / PRM_data[, pool_cols[[i]]]
}
PRM_data <- log2(subset(PRM_data, select = unlist(batch_samples)))
# remove samples
PRM_data <- PRM_data %>% select(-remove_list)

# perform DEP test on PRM data ----
data_input <- as.data.frame(apply(PRM_data, 1, function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE) * fill_na_value
  return(x)
}))
data_input$group <- unlist(lapply(rownames(data_input), function(x) {
  if (x %in% s1_list) {
    return("ov_s1")
  } else if (grepl("CA", x)) {
    return("ov")
  } else if (grepl("CYS", x)) {
    return("cys")
  } else if (grepl("NC", x)) {
    return("ctrl")
  }
}))
colnames(data_input) <- lapply(colnames(data_input), function(x) {
  return(str_split_1(x, "_")[1])
})
# DEPs
DEPs_ov_s1_cys <- compare_groups(data = data_input, group1 = "ov_s1", group2 = "cys", p_adjust_method = adjust_method)
DEPs_ov_cys <- compare_groups(data = data_input, group1 = "ov", group2 = "cys", p_adjust_method = adjust_method)
DEPs_ov_s1_ctrl <- compare_groups(data = data_input, group1 = "ov_s1", group2 = "ctrl", p_adjust_method = adjust_method)
DEPs_cys_ctrl <- compare_groups(data = data_input, group1 = "cys", group2 = "ctrl", p_adjust_method = adjust_method)
DEPs_ov_ctrl <- compare_groups(data = data_input, group1 = "ov", group2 = "ctrl", p_adjust_method = adjust_method)
# select DEPs
DEPs_ov_s1_cys <- DEPs_ov_s1_cys[abs(DEPs_ov_s1_cys$logFC) > logfc_threshold & DEPs_ov_s1_cys$adj.P.Val < pvalue_threshold, ]
DEPs_ov_cys <- DEPs_ov_cys[abs(DEPs_ov_cys$logFC) > logfc_threshold & DEPs_ov_cys$adj.P.Val < pvalue_threshold, ]
DEPs_ov_s1_ctrl <- DEPs_ov_s1_ctrl[abs(DEPs_ov_s1_ctrl$logFC) > logfc_threshold & DEPs_ov_s1_ctrl$adj.P.Val < pvalue_threshold, ]
DEPs_cys_ctrl <- DEPs_cys_ctrl[abs(DEPs_cys_ctrl$logFC) > logfc_threshold & DEPs_cys_ctrl$adj.P.Val < pvalue_threshold, ]
DEPs_ov_ctrl <- DEPs_ov_ctrl[abs(DEPs_ov_ctrl$logFC) > logfc_threshold & DEPs_ov_ctrl$adj.P.Val < pvalue_threshold, ]

# compare DEP results with training cohort data ----
compare_list <- c("DEPs_ov_cys", "DEPs_cys_ctrl", "DEPs_ov_ctrl", "DEPs_ov_s1_cys", "DEPs_ov_s1_ctrl")
compare_res <- c()
for (item in compare_list) {
  train_res <- readRDS(paste0(DEPs_path, item, ".rds"))
  train_res_raw <- readRDS(paste0(DEPs_path, item, "_raw.rds"))
  feas_num <- nrow(train_res)
  prm_res <- get(item)
  common_proteins <- intersect(rownames(train_res), rownames(prm_res))
  train_res <- train_res[common_proteins, ]
  prm_res <- prm_res[common_proteins, ]

  train_res_raw$Feature <- lapply(train_res_raw$Feature, function(x) {
    return(index_symbol_mapping[index_symbol_mapping$Index == x, 2])
  })
  common_proteins_symbol <- index_symbol_mapping[index_symbol_mapping$Index %in% common_proteins, 2]

  train_res_raw <- train_res_raw %>% arrange(desc(abs(logFC)))

  # generate volcano plot
  if (item == "DEPs_ov_s1_cys") {
    plot_title <- "stage1-OV / CYS"
    custom_label_num <- 20
  } else if (item == "DEPs_ov_cys") {
    plot_title <- "later-stage OV / CYS"
    custom_label_num <- 20
  } else if (item == "DEPs_ov_s1_ctrl") {
    plot_title <- "stage1-OV / CTRL"
    custom_label_num <- 20
  } else if (item == "DEPs_cys_ctrl") {
    plot_title <- "CYS / CTRL"
    custom_label_num <- 20
  } else if (item == "DEPs_ov_ctrl") {
    plot_title <- "OV / CTRL"
    custom_label_num <- 20
  }

  if (custom_label_num > length(common_proteins_symbol)) {
    custom_label_num <- length(common_proteins_symbol)
  }
  plot_title <- paste0("Top ", custom_label_num, " ", plot_title, " Features")

  p1 <- gradual_volcano(train_res_raw,
    x = "logFC", y = "adj.P.Val",
    fills = brewer.pal(5, "RdYlBu"),
    colors = "black",
    label = "Feature",
    output = FALSE, custom_label = common_proteins_symbol[1:custom_label_num],
    log2FC_cut = 0.25, x_lab = paste0("log2FC"),
    label_size = 3, label_max_overlaps = 35,
    pointSizeRange = c(0.5, 4)
  )
  p1 <- p1 + theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    legend.position = "none",
    title = element_text(size = 15)
  ) + ggtitle(plot_title)
  ggsave(paste0(output_path, item, "_volcano.tiff"),
    units = "cm", height = 15, width = 15, dpi = 600
  )
  # compare
  colnames(prm_res) <- paste0(colnames(prm_res), ".prm")
  final_res <- cbind(train_res, prm_res)
  final_res$same <- final_res$logFC * final_res$logFC.prm > 0
  compare_res <- c(compare_res, list(final_res))
  saveRDS(final_res, paste0(output_path, item, "_compare.rds"))
  message(
    item, " have ", sum(final_res$same), " same direction DEPs, ", "total ",
    nrow(final_res), " proteins get compared."
  )
}
saveRDS(data_input, paste0(output_path, "PRM_data_filled.rds"))


# checkpoint ----
save.image("./data/analysis_check_point/step4_PRM_eval.RData")
rm(list = ls())
