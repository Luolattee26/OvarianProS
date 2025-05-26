# packages and functions ----
rm(list = ls())
source("./code/utils.R")
## analysis init settings ----
output_path <- "./output/step3_ML_data/"
DEPs_path <- "./output/step2_DEPs/"
input_exp_file <- "./data/samples/Protein_exp/ratio/ratio_after_QC_corrected.rds"
if (!dir.exists(output_path)) {
  dir.create(output_path)
}


# read data and preprocess ---- First to read the DEPs data
DEPs_ov_s1_cys <- readRDS(paste0(DEPs_path, "DEPs_ov_s1_cys.rds"))
DEPs_ov_cys <- readRDS(paste0(DEPs_path, "DEPs_ov_cys.rds"))
DEPs_ov_s1_ctrl <- readRDS(paste0(DEPs_path, "DEPs_ov_s1_ctrl.rds"))
DEPs_cys_ctrl <- readRDS(paste0(DEPs_path, "DEPs_cys_ctrl.rds"))
DEPs_ov_ctrl <- readRDS(paste0(DEPs_path, "DEPs_ov_ctrl.rds"))
# read exp data
exp_data <- readRDS(input_exp_file)

# get the DEPs data for ML ----

# all 4 models: ov_s1_cys, ov_cys, ov_s1_ctrl, cys_ctrl ov_s1_cys
overlap <- intersect(rownames(DEPs_ov_s1_cys), rownames(DEPs_ov_s1_ctrl))
message("ov_s1_cys: ", length(overlap))
ml_data <- cbind(exp_data[, overlap], exp_data$group)
colnames(ml_data)[ncol(ml_data)] <- "label"
# save to csv file
write.csv(ml_data, paste0(output_path, "ML_data_ov_s1_cys.csv"), row.names = T)
## ov_cys ----
overlap <- intersect(rownames(DEPs_ov_cys), rownames(DEPs_ov_ctrl))
message("ov_cys: ", length(overlap))
ml_data <- cbind(exp_data[, overlap], exp_data$group)
colnames(ml_data)[ncol(ml_data)] <- "label"
# save to csv file
write.csv(ml_data, paste0(output_path, "ML_data_ov_cys.csv"), row.names = T)
## ov_s1_ctrl ----
overlap <- rownames(DEPs_ov_s1_ctrl)
message("ov_s1_ctrl: ", length(overlap))
ml_data <- cbind(exp_data[, overlap], exp_data$group)
colnames(ml_data)[ncol(ml_data)] <- "label"
# save to csv file
write.csv(ml_data, paste0(output_path, "ML_data_ov_s1_ctrl.csv"), row.names = T)
## cys_ctrl ----
overlap <- intersect(rownames(DEPs_cys_ctrl), rownames(DEPs_ov_ctrl))
message("cys_ctrl: ", length(overlap))
ml_data <- cbind(exp_data[, overlap], exp_data$group)
colnames(ml_data)[ncol(ml_data)] <- "label"
# save to csv file
write.csv(ml_data, paste0(output_path, "ML_data_cys_ctrl.csv"), row.names = T)


# checkpoint ----
save.image("./data/analysis_check_point/step3_get_ML_data.RData")
rm(list = ls())
