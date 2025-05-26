# packages and functions ----
rm(list = ls())
source("./code/utils.R")
library(limma)
## analysis init settings ----
output_path <- "./output/step2_DEPs/"
input_exp_file <- "./data/samples/Protein_exp/ratio/ratio_after_QC_corrected.rds"
groups <- c("ov_s1", "ov", "cys", "ctrl")
logfc_threshold <- 0.25
pvalue_threshold <- 0.05
adjust_method <- "BH"


# read data and preprocess ---- First to make sure the group is in factor type
exp_data <- readRDS(input_exp_file)
exp_data$group <- factor(exp_data$group, levels = groups)
exp <- t(exp_data[, 1:(dim(exp_data)[2] - 2)])
# Then generate design matrix and contrast matrix
design <- model.matrix(~ 0 + exp_data$group)
colnames(design) <- groups
rownames(design) <- rownames(exp_data)
contrast <- makeContrasts(ov_s1 - cys, ov - cys, ov_s1 - ctrl, cys - ctrl, ov - ctrl, ov_s1 - ov, levels = design)


# DEPs analysis ---- Next we perform the DEPs analysis
fit <- lmFit(exp, design)
contrast_fit <- contrasts.fit(fit, contrast)
fit_final <- eBayes(contrast_fit)
## get the DEPs of each comparison ----
DEPs_ov_s1_cys <- topTable(fit_final, coef = 1, number = Inf, adjust.method = adjust_method)
DEPs_ov_cys <- topTable(fit_final, coef = 2, number = Inf, adjust.method = adjust_method)
DEPs_ov_s1_ctrl <- topTable(fit_final, coef = 3, number = Inf, adjust.method = adjust_method)
DEPs_cys_ctrl <- topTable(fit_final, coef = 4, number = Inf, adjust.method = adjust_method)
DEPs_ov_ctrl <- topTable(fit_final, coef = 5, number = Inf, adjust.method = adjust_method)
DEPs_ov_s1_ov <- topTable(fit_final, coef = 6, number = Inf, adjust.method = adjust_method)
## filter the DEPs by logFC and pvalue ----
DEPs_ov_s1_cys <- DEPs_ov_s1_cys[abs(DEPs_ov_s1_cys$logFC) > logfc_threshold & DEPs_ov_s1_cys$adj.P.Val < pvalue_threshold, ]
DEPs_ov_cys <- DEPs_ov_cys[abs(DEPs_ov_cys$logFC) > logfc_threshold & DEPs_ov_cys$adj.P.Val < pvalue_threshold, ]
DEPs_ov_s1_ctrl <- DEPs_ov_s1_ctrl[abs(DEPs_ov_s1_ctrl$logFC) > logfc_threshold & DEPs_ov_s1_ctrl$adj.P.Val < pvalue_threshold, ]
DEPs_cys_ctrl <- DEPs_cys_ctrl[abs(DEPs_cys_ctrl$logFC) > logfc_threshold & DEPs_cys_ctrl$adj.P.Val < pvalue_threshold, ]
DEPs_ov_ctrl <- DEPs_ov_ctrl[abs(DEPs_ov_ctrl$logFC) > logfc_threshold & DEPs_ov_ctrl$adj.P.Val < pvalue_threshold, ]
## save the DEPs ----
if (!dir.exists(output_path)) {
  dir.create(output_path)
}
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


# checkpoint ----
save.image("./data/analysis_check_point/step2_DEPs.RData")
rm(list = ls())
