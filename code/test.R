ml_data <- readRDS('/Users/calluo/Desktop/ml_data_all.rds')
data_s1 <- data.frame(ml_data$DEPs_ov_s1_cys)
data_all <- data.frame(ml_data$DEPs_ov_cys)
length(intersect(colnames(data_s1), colnames(data_all)))

data_s1 <- data_s1[, intersect(colnames(data_s1), colnames(data_all))]
data_all <- data_all[, intersect(colnames(data_s1), colnames(data_all))]


write.csv(data_all, '/Users/calluo/Desktop/ml_data_cys_ov.csv')
