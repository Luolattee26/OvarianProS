# custom functions for data analysis---- packages used ----
if (!requireNamespace("rstatix", quietly = TRUE)) {
  install.packages("rstatix")
}
if (!requireNamespace("Hmisc", quietly = TRUE)) {
  install.packages("Hmisc")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel")
}
if (!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}
if (!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach")
}
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)
library(rstatix)
library(Hmisc)
library(logger)
library(ggrepel)
library(latex2exp)


# read index and symbol mapping automatically ----
if (file.exists("./data/index_to_gene_symbol.rds")) {
  index_symbol_mapping <- readRDS("./data/index_to_gene_symbol.rds")
}

# function to separate the pool samples and normal(batch) samples ----
separater_pool_batch <- function(colnames_to_separate) {
  pool_mask <- unlist(lapply(colnames_to_separate, function(x) {
    if (grepl("pool", x)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))

  batch_mask <- unlist(lapply(colnames_to_separate, function(x) {
    if (grepl("^b\\d+_\\d+$", x)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))

  res <- list(pool = pool_mask, batch = batch_mask)

  return(res)
}

# function to remove and fill the missing values ----
remove_fill_na <- function(remaining_cols, data, missing_threshold, fill_na_value) {
  remain_data <- data[, remaining_cols]
  data <- data[, -remaining_cols]

  delete_mask <- apply(data, 1, function(x) sum(is.na(x)) > missing_threshold * length(x))
  data <- data[!delete_mask, ]
  remain_data <- remain_data[!delete_mask, ]

  data_filled <- t(apply(data, 1, function(x) {
    x[is.na(x)] <- min(x, na.rm = TRUE) * fill_na_value
    return(x)
  }))

  return(cbind(remain_data, data_filled))
}

# function to calculate the CV ----
calculate_CV <- function(data, verbose = T, percent = T) {
  CV <- unlist(apply(data, 1, function(x) {
    if (length(unique(x)) == 1) {
      return(0)
    } else {
      return(sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
    }
  }))

  if (percent) {
    CV <- CV * 100
  }

  if (verbose) {
    print(paste("CV mean:", mean(CV, na.rm = T)))
    print(paste("CV median:", median(CV, na.rm = T)))
  }

  return(CV)
}

# function to perform pre-processing of ratio data ----
preprocess_ratio_data <- function(data, group_data, transpose = T) {
  if (transpose) {
    data <- as.data.frame(t(data))
  }

  data$group <- unlist(lapply(rownames(data), function(x) {
    if (x %in% group_data$X) {
      return(group_data$group[group_data$X == x])
    } else {
      return("BAD")
    }
  }))
  data <- data[data$group != "BAD", ]

  data$batch <- unlist(lapply(rownames(data), function(x) {
    return(strsplit(x, "_")[[1]][1])
  }))

  return(data)
}

# function to get color for the plot ----
get_color <- function(seed = 2028, n = 5) {
  allcolour <- c(
    "#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#DEB887", "#808000",
    "#FF00FF", "#FA8072", "#7B68EE", "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493",
    "#0000CD", "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", "#FF6347",
    "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", "#FFFF00"
  )

  set.seed(seed)

  if (n > length(allcolour)) {
    message("The number of colors required is larger than the number of colors available")
    message("So the colors will be randomly selected and will be repeated")
    return(sample(allcolour, n, replace = T))
  } else {
    return(sample(allcolour, n))
  }
}

# function to perform wilcox test ----
compare_groups <- function(data, group1, group2, p_adjust_method = "BH") {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data frame or matrix.")
  }

  data[, ncol(data)] <- as.factor(data[, ncol(data)])

  groups <- unique(data[, ncol(data)])
  if (!(group1 %in% groups) || !(group2 %in% groups)) {
    stop("Specified groups do not exist in the dataset.")
  }

  p_values <- numeric(ncol(data) - 1)
  log_fold_changes <- numeric(ncol(data) - 1)
  names(p_values) <- colnames(data)[1:(ncol(data) - 1)]
  names(log_fold_changes) <- colnames(data)[1:(ncol(data) - 1)]

  for (i in 1:(ncol(data) - 1)) {
    feature_data <- data[, i]

    group1_data <- feature_data[data[, ncol(data)] == group1]
    group2_data <- feature_data[data[, ncol(data)] == group2]

    test_result <- wilcox.test(group1_data, group2_data, alternative = "two.sided")
    p_values[i] <- test_result$p.value

    mean_group1 <- mean(2^group1_data, na.rm = TRUE)
    mean_group2 <- mean(2^group2_data, na.rm = TRUE)
    log_fold_changes[i] <- log2(mean_group1 / mean_group2)
  }

  adjusted_p_values <- p.adjust(p_values, method = p_adjust_method)

  results <- data.frame(Feature = names(p_values), P_Value = p_values, adj.P.Val = adjusted_p_values, logFC = log_fold_changes)

  return(results)
}



gradual_volcano <- function(data,
                            x = "log2FoldChange", y = "padj",
                            pointSizeRange = c(0.5, 4),
                            fills = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25"),
                            colors = c("#17194e", "#68bfe7", "#f9ed36", "#a22f27", "#211f1f"),
                            x_lab = NULL, y_lab = NULL,
                            legend_title = NULL, legend_position = NULL,
                            add_line = TRUE, log2FC_cut = 1, FDR_cut = 0.05,
                            add_label = TRUE, label = "row",
                            label_number = 10, custom_label = NULL,
                            output = TRUE, filename = "volcano_plot",
                            label_size = 4.5, label_max_overlaps = 100) {
  colnames(data)[colnames(data) == x] <- "x"
  colnames(data)[colnames(data) == y] <- "y"
  colnames(data)[colnames(data) == label] <- "geneName"

  # 默认排序标签 or 自定义标签
  if (is.null(custom_label)) {
    if (label_number != 0) {
      data$label <- rep("", nrow(data))
      data$label[order(data$y)[1:label_number]] <- data$geneName[order(data$y)[1:label_number]]
    } else {
      data$label <- rep("", nrow(data))
    }
  } else {
    data$label <- rep("", nrow(data))
    data$label[match(custom_label, data$geneName)] <- custom_label
  }

  p <- ggplot(data, aes(x, -log10(y)))

  p <- p +
    # 散点图:
    geom_point(aes(size = -log10(y), fill = -log10(y), color = -log10(y)), shape = 21) +
    # 指定颜色渐变模式：
    scale_fill_gradientn(
      values = seq(0, 1, 0.2),
      colors = fills
    ) +
    scale_color_gradientn(
      values = seq(0, 1, 0.2),
      colors = colors
    ) +
    # 指定散点大小渐变模式：
    scale_size_continuous(range = pointSizeRange) +
    # 主题调整：
    theme_bw() +
    # 调整主题和图例位置：
    theme(
      panel.grid = element_blank(),
      legend.background = element_blank()
    ) +
    # 设置部分图例不显示：
    guides(
      fill = guide_colourbar(title = legend_title %||% "-Log10_q-value"),
      color = "none",
      size = "none"
    ) +
    # 修改坐标轴：
    labs(
      x = x_lab %||% TeX("$Log_2 \\textit{FC}$"),
      y = y_lab %||% TeX("$-Log_{10} \\textit{FDR} $")
    )

  # 添加虚线：
  if (add_line == TRUE) {
    p <- p +
      geom_vline(xintercept = c(-log2FC_cut, log2FC_cut), linetype = "dashed") +
      geom_hline(yintercept = -log10(FDR_cut), linetype = "dashed")
  }

  # 添加标签：
  if (add_label == TRUE) {
    p <- p +
      geom_text_repel(aes(label = label, color = -log10(y)),
        size = label_size, max.overlaps = label_max_overlaps, key_glyph = draw_key_point, color = "red",
        fontface = "bold"
      )
  }

  # 图例位置：
  if (is.null(legend_position)) {
    p <- p +
      theme(
        legend.position = c(0.01, 0.7),
        legend.justification = c(0, 1)
      )
  } else if (legend_position == "UL") {
    p <- p +
      theme(
        legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
      )
  } else if (legend_position == "UR") {
    p <- p +
      theme(
        legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
      )
  } else if (legend_position == "DL") {
    p <- p +
      theme(
        legend.position = c(0.01, 0.01),
        legend.justification = c(0, 0),
      )
  } else {
    p <- p +
      theme(
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
      )
  }

  # 保存图片：
  if (output == TRUE) {
    ggsave(paste0(filename, ".pdf"), plot = p, height = 5, width = 6)
  }
  return(p)
}
