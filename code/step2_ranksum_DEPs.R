# packages and functions ----
rm(list = ls())
source("./code/utils.R")
## analysis init settings ----
output_path <- "./output/step2_DEPs/"
input_exp_file <- "./data/samples/Protein_exp/ratio/ratio_after_QC_corrected.rds"
groups <- c("ov_s1", "ov", "cys", "ctrl")
logfc_threshold <- 0.25
pvalue_threshold <- 0.05
adjust_method <- "BH"
if (!dir.exists(output_path)) {
  dir.create(output_path)
}


# read data and preprocess ----
exp_data <- readRDS(input_exp_file)
exp_data$group <- factor(exp_data$group, levels = groups)
data_to_input <- exp_data[, -ncol(exp_data)]


# DEPs analysis ----
DEPs_ov_s1_cys <- compare_groups(data = data_to_input, group1 = "ov_s1", group2 = "cys", p_adjust_method = adjust_method)
DEPs_ov_cys <- compare_groups(data = data_to_input, group1 = "ov", group2 = "cys", p_adjust_method = adjust_method)
DEPs_ov_s1_ctrl <- compare_groups(data = data_to_input, group1 = "ov_s1", group2 = "ctrl", p_adjust_method = adjust_method)
DEPs_cys_ctrl <- compare_groups(data = data_to_input, group1 = "cys", group2 = "ctrl", p_adjust_method = adjust_method)
DEPs_ov_ctrl <- compare_groups(data = data_to_input, group1 = "ov", group2 = "ctrl", p_adjust_method = adjust_method)
DEPs_ov_s1_ov <- compare_groups(data = data_to_input, group1 = "ov_s1", group2 = "ov", p_adjust_method = adjust_method)
# save raw DEPs
saveRDS(DEPs_ov_s1_cys, paste0(output_path, "DEPs_ov_s1_cys_raw.rds"))
saveRDS(DEPs_ov_cys, paste0(output_path, "DEPs_ov_cys_raw.rds"))
saveRDS(DEPs_ov_s1_ctrl, paste0(output_path, "DEPs_ov_s1_ctrl_raw.rds"))
saveRDS(DEPs_cys_ctrl, paste0(output_path, "DEPs_cys_ctrl_raw.rds"))
saveRDS(DEPs_ov_ctrl, paste0(output_path, "DEPs_ov_ctrl_raw.rds"))
## filter the DEPs by logFC and pvalue ----
DEPs_ov_s1_cys <- DEPs_ov_s1_cys[abs(DEPs_ov_s1_cys$logFC) > logfc_threshold & DEPs_ov_s1_cys$adj.P.Val < pvalue_threshold, ]
DEPs_ov_cys <- DEPs_ov_cys[abs(DEPs_ov_cys$logFC) > logfc_threshold & DEPs_ov_cys$adj.P.Val < pvalue_threshold, ]
DEPs_ov_s1_ctrl <- DEPs_ov_s1_ctrl[abs(DEPs_ov_s1_ctrl$logFC) > logfc_threshold & DEPs_ov_s1_ctrl$adj.P.Val < pvalue_threshold, ]
DEPs_cys_ctrl <- DEPs_cys_ctrl[abs(DEPs_cys_ctrl$logFC) > logfc_threshold & DEPs_cys_ctrl$adj.P.Val < pvalue_threshold, ]
DEPs_ov_ctrl <- DEPs_ov_ctrl[abs(DEPs_ov_ctrl$logFC) > logfc_threshold & DEPs_ov_ctrl$adj.P.Val < pvalue_threshold, ]
## save the DEPs ----
saveRDS(DEPs_ov_s1_cys, paste0(output_path, "DEPs_ov_s1_cys.rds"))
saveRDS(DEPs_ov_cys, paste0(output_path, "DEPs_ov_cys.rds"))
saveRDS(DEPs_ov_s1_ctrl, paste0(output_path, "DEPs_ov_s1_ctrl.rds"))
saveRDS(DEPs_cys_ctrl, paste0(output_path, "DEPs_cys_ctrl.rds"))
saveRDS(DEPs_ov_ctrl, paste0(output_path, "DEPs_ov_ctrl.rds"))
# show the number of DEPs
message("DEPs_ov_s1_cys: ", dim(DEPs_ov_s1_cys)[1])
message("DEPs_ov_cys: ", dim(DEPs_ov_cys)[1])
message("DEPs_ov_s1_ctrl: ", dim(DEPs_ov_s1_ctrl)[1])
message("DEPs_cys_ctrl: ", dim(DEPs_cys_ctrl)[1])
message("DEPs_ov_ctrl: ", dim(DEPs_ov_ctrl)[1])
# save features
write.csv(DEPs_ov_s1_cys$Feature, file = paste0(output_path, "DEPs_ov_s1_cys.csv"), row.names = FALSE)
write.csv(DEPs_ov_cys$Feature, file = paste0(output_path, "DEPs_ov_cys.csv"), row.names = FALSE)
write.csv(DEPs_ov_s1_ctrl$Feature, file = paste0(output_path, "DEPs_ov_s1_ctrl.csv"), row.names = FALSE)
write.csv(DEPs_cys_ctrl$Feature, file = paste0(output_path, "DEPs_cys_ctrl.csv"), row.names = FALSE)
write.csv(DEPs_ov_ctrl$Feature, file = paste0(output_path, "DEPs_ov_ctrl.csv"), row.names = FALSE)


# checkpoint ----
save.image("./data/analysis_check_point/step2_DEPs.RData")
rm(list = ls())
