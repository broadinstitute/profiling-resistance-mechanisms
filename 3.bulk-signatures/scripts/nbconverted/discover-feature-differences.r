suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))

viz_file <- file.path("..", "scripts", "visualization_utils.R")
source(viz_file)

util_file <- file.path("..", "scripts", "processing_utils.R")
source(util_file)

doses <- c(0.7, 7)

axis_title_size <- 10
axis_text_size <- 9
strip_text_size <- 8
ggrepel_label_size <- 1.9
title_text_size <- 10
ymax <- 9

file <- file.path("..", "data", "merged_intersected_variable_selected.csv")
data_df <- load_data(file)

print(dim(data_df))
head(data_df, 3)

tstat_dose <- c()
pval_dose <- c()

tstat_cell <- c()
pval_cell <- c()

all_features <- c()
for (feature in colnames(data_df)) {
    if (!grepl("Metadata_", feature)) {
        # Perform Dose Experiment
        doseA <- data_df %>%
            dplyr::filter(Metadata_Dosage == doses[1],
                          Metadata_CellLine != "WT") %>%
            dplyr::pull(!!feature)
        
        doseB <- data_df %>%
            dplyr::filter(Metadata_Dosage == doses[2],
                          Metadata_CellLine != "WT") %>%
            dplyr::pull(!!feature)
        
        result <- t.test(doseA, doseB, var.equal = FALSE)
        
        tstat_dose <- c(tstat_dose, as.numeric(paste(result$statistic)))
        pval_dose <- c(pval_dose, result$p.value)
        
        # Perform Cell Line Experiment at 0.7uM Dose
        resistant_clones <- data_df %>%
            dplyr::filter(Metadata_Dosage == 0.7,
                          Metadata_CellLine != "WT") %>%
            dplyr::pull(!!feature)
        
        wt_cells <- data_df %>%
            dplyr::filter(Metadata_Dosage == 0.7,
                          Metadata_CellLine == "WT") %>%
            dplyr::pull(!!feature)
        
        result <- t.test(resistant_clones, wt_cells, var.equal = FALSE)
        
        tstat_cell <- c(tstat_cell, as.numeric(paste(result$statistic)))
        pval_cell <- c(pval_cell, result$p.value)
        
        # Track which feature is being tested
        all_features <- c(all_features, feature)
    }
}

# Obtain Results for Dose Differences
result_dose_df <- as.data.frame(cbind(tstat_dose, pval_dose))
result_dose_df$neglog10p <- -log10(result_dose_df$pval_dose)
result_dose_df$feature <- all_features
result_dose_df <- result_dose_df %>% dplyr::arrange(desc(neglog10p))

dim(result_dose_df)
head(result_dose_df, 3)

# Obtain Results for Cell Line Differences
result_cell_df <- as.data.frame(cbind(tstat_cell, pval_cell))
result_cell_df$neglog10p <- -log10(result_cell_df$pval_cell)
result_cell_df$feature <- all_features
result_cell_df <- result_cell_df %>% dplyr::arrange(desc(neglog10p))

dim(result_cell_df)
head(result_cell_df, 3)

alpha_correction <- -log10(0.05 / dim(result_dose_df)[1])
repel_logic <- result_dose_df$neglog10p > alpha_correction * 1.5

ttest_dose_gg <- ttest_volcano(
    df = result_dose_df,
    x_string = "tstat_dose",
    title = "0.7nM vs. 7nM Dose (Resistant Clones)",
    title_text_size = title_text_size,
    yintercept = alpha_correction,
    repel_logic = repel_logic,
    ggrepel_label_size = ggrepel_label_size,
    axis_text_size = axis_text_size,
    axis_title_size = axis_title_size,
    ymax = ymax
)

alpha_correction <- -log10(0.05 / dim(result_cell_df)[1])
repel_logic <- result_cell_df$neglog10p > alpha_correction * 1.1

ttest_cell_gg <- ttest_volcano(
    df = result_cell_df,
    x_string = "tstat_cell",
    title = "Resistant vs. Wildtype Cells (0.7nM)",
    title_text_size = title_text_size,
    yintercept = alpha_correction,
    repel_logic = repel_logic,
    ggrepel_label_size = ggrepel_label_size,
    axis_text_size = axis_text_size,
    axis_title_size = axis_title_size,
    ymax = ymax
)

# The top feature is something to do with nuclear area
top_feature <- paste(result_dose_df$feature[1])

append_batch <- function(string) paste("Batch:", string)
append_dose <- function(string) paste0("Dose: ", string, "nM")

distrib_gg <- ggplot(data_df, aes_string(x = `top_feature`)) +
    geom_density(aes(fill = Metadata_CellLine),
                 alpha = 0.6) +
    geom_rug(aes(color = Metadata_CellLine),
             alpha = 0.8,
             size = 0.5) +
    facet_grid(Metadata_Dosage ~ Metadata_Batch_Number,
               scales = "fixed",
               labeller = labeller(Metadata_Batch_Number = as_labeller(append_batch),
                                   Metadata_Dosage = as_labeller(append_dose))) +
    scale_color_manual(name = "Cell Line",
                      labels = c("CloneA" = "Clone A",
                                 "CloneE" = "Clone E",
                                 "WT" = "Wild type"),
                      values = c("CloneA" = "#1b9e77",
                                 "CloneE" = "#d95f02",
                                 "WT" = "#7570b3")) +
    scale_fill_manual(name = "Cell Line",
                      labels = c("CloneA" = "Clone A",
                                 "CloneE" = "Clone E",
                                 "WT" = "Wild type"),
                      values = c("CloneA" = "#1b9e77",
                                 "CloneE" = "#d95f02",
                                 "WT" = "#7570b3")) +
    theme_bw() +
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          strip.text = element_text(size = strip_text_size),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

table(data_df$Metadata_CellLine,
      data_df$Metadata_Dosage,
      data_df$Metadata_Batch_Number)

feature_plot <- (
    cowplot::plot_grid(
        ttest_cell_gg,
        ttest_dose_gg,
        labels = c("a", "b"),
        ncol = 1,
        nrow = 2
    )
)

main_plot <- (
    cowplot::plot_grid(
        feature_plot,
        distrib_gg,
        labels = c("", "c"),
        ncol = 2,
        nrow = 1,
        rel_widths = c(0.7, 1)
    )
)

main_plot

file_base <- file.path("figures", "dosage_feature_figure")
save_figure(main_plot, file_base, height = 6, width = 8)

batch <- "2019_06_25_Batch3"
file <- file.path("..", "data", paste0(batch, "_merged_intersected_variable_selected.csv"))

data_df <- load_data(file)

print(dim(data_df))
head(data_df, 3)

tstat_cell <- c()
pval_cell <- c()

all_features <- c()
for (feature in colnames(data_df)) {
    if (!grepl("Metadata_", feature)) {
        
        # Perform Cell Line Experiment at 0.7uM Dose
        resistant_clones <- data_df %>%
            dplyr::filter(Metadata_Plate == "MutClones") %>%
            dplyr::pull(!!feature)
        
        wt_cells <- data_df %>%
            dplyr::filter(Metadata_Plate == "WTClones") %>%
            dplyr::pull(!!feature)
        
        result <- t.test(resistant_clones, wt_cells, var.equal = FALSE)
        
        tstat_cell <- c(tstat_cell, as.numeric(paste(result$statistic)))
        pval_cell <- c(pval_cell, result$p.value)
        
        # Track which feature is being tested
        all_features <- c(all_features, feature)
    }
}

# Obtain Results for Cell Line Differences
result_cell_df <- as.data.frame(cbind(tstat_cell, pval_cell))
result_cell_df$neglog10p <- -log10(result_cell_df$pval_cell)
result_cell_df$feature <- all_features
result_cell_df <- result_cell_df %>% dplyr::arrange(desc(neglog10p))

dim(result_cell_df)
head(result_cell_df, 3)

alpha_correction <- -log10(0.05 / dim(result_cell_df)[1])
repel_logic <- (
    result_cell_df$neglog10p > alpha_correction | 
    result_cell_df$tstat_cell > 2.65 |
    result_cell_df$tstat_cell < -2.4
    )

ttest_cell_gg <- ttest_volcano(
    df = result_cell_df,
    x_string = "tstat_cell",
    title = "Resistant vs. Wildtype Cells",
    title_text_size = title_text_size,
    yintercept = alpha_correction,
    repel_logic = repel_logic,
    ggrepel_label_size = ggrepel_label_size,
    axis_text_size = axis_text_size,
    axis_title_size = axis_title_size,
    ymax = ymax
)

# The top feature is something to do with nuclear area
top_feature <- paste(result_cell_df$feature[1])

distrib_gg <- ggplot(data_df,
                     aes_string(x = `top_feature`)) +
    geom_density(aes(fill = Metadata_Plate_Map_Name),
                 alpha = 0.6) +
    geom_rug(aes(color = Metadata_Plate_Map_Name),
             alpha = 0.8,
             size = 0.5) +
    scale_color_manual(name = "Cell Line",
                      labels = c("MutClones" = "Mutant",
                                 "WTClones" = "Wild-type"),
                      values = c("MutClones" = "#1b9e77",
                                 "WTClones" = "#d95f02")) +
    scale_fill_manual(name = "Cell Line",
                      labels = c("MutClones" = "Mutant",
                                 "WTClones" = "Wild-type"),
                      values = c("MutClones" = "#1b9e77",
                                 "WTClones" = "#d95f02")) +
    theme_bw() +
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          strip.text = element_text(size = strip_text_size),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

top_neg <- result_cell_df %>%
    dplyr::arrange(tstat_cell) %>%
    select(feature)

top_feature = top_neg$feature[1]

distrib_neg_gg <- ggplot(data_df,
                         aes_string(x = `top_feature`)) +
    geom_density(aes(fill = Metadata_Plate_Map_Name),
                 alpha = 0.6) +
    geom_rug(aes(color = Metadata_Plate_Map_Name),
             alpha = 0.8,
             size = 0.5) +
    scale_color_manual(name = "Cell Line",
                      labels = c("MutClones" = "Mutant",
                                 "WTClones" = "Wild-type"),
                      values = c("MutClones" = "#1b9e77",
                                 "WTClones" = "#d95f02")) +
    scale_fill_manual(name = "Cell Line",
                      labels = c("MutClones" = "Mutant",
                                 "WTClones" = "Wild-type"),
                      values = c("MutClones" = "#1b9e77",
                                 "WTClones" = "#d95f02")) +
    theme_bw() +
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          strip.text = element_text(size = strip_text_size),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

full_distrib_gg <- (
    cowplot::plot_grid(
        distrib_gg,
        distrib_neg_gg,
        labels = c("b", "c"),
        nrow = 2
    )
)

main_plot <- (
    cowplot::plot_grid(
        ttest_cell_gg,
        full_distrib_gg,
        labels = c("a", ""),
        ncol = 2,
        nrow = 1,
        rel_widths = c(1, 1)
    )
)

main_plot

file_base <- file.path("figures", paste0(batch, "_feature_figure"))
save_figure(main_plot, file_base, height = 6, width = 8.5)
