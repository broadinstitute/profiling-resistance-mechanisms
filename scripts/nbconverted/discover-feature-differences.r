
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))

# Load Plotting Function for ttest Volcano Plots
ttest_volcano <- function(df, x_string, title, yintercept,
                          repel_logic, ggrepel_label_size,
                          title_text_size, axis_text_size,
                          axis_title_size) {
    # Plot the results of a t-test on various samples for common features
    #
    # Arguments:
    # df - the dataframe storing the t-test results
    # x_string - string indicating what variable to plot on x axis
    # title - string indicating the title of the plot
    # yintercept - an alpha corrected value to plot a red dotted line
    # repel_logic - which features to highlight and name
    # ggrepel_label_size - int of the size of the feature names
    # title_text_size - int of the size of the title text
    # axis_text_size - int of the size of the text on the ggplot axes
    # axis_title_size - int of the size of the titles on the ggplot axes
    #
    # Output:
    # The ggplot2 object for downstream saving
    ttest_gg <- ggplot(df,
                       aes_string(x = x_string,
                                  y = "neglog10p")) +
        geom_point(alpha = 0.5,
                   size = 0.8,
                   color = ifelse(repel_logic, "red", "grey50")) +
        geom_hline(yintercept = yintercept,
                   color = "red",
                   linetype = "dashed") +
        xlab("t Statistic") +
        ylab("-log10 P") +
        geom_text_repel(data = subset(df, repel_logic),
                        arrow = arrow(length = unit(0.01, "npc")),
                        size = ggrepel_label_size,
                        segment.size = 0.1,
                        segment.alpha = 0.8,
                        force = 20,
                        aes_string(x = x_string,
                                   y = "neglog10p",
                                   label = "feature")) +
        ggtitle(title) +
        theme_bw() +
        theme(plot.title = element_text(size = title_text_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))

    return(ttest_gg)
}

doses <- c(0.7, 7)

axis_title_size <- 10
axis_text_size <- 9
strip_text_size <- 8
ggrepel_label_size <- 1.9
title_text_size <- 10

# Set column types for reading in data
batch_cols = readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_Assay_Plate_Barcode = readr::col_character(),
    Metadata_Plate_Map_Name = readr::col_character(),
    Metadata_Batch_Number = readr::col_integer(),
    Metadata_well_position = readr::col_character(),
    Metadata_CellLine = readr::col_character()
)

file <- file.path("data", "merged_intersected_variable_selected.csv")
data_df <- readr::read_csv(file, col_types = batch_cols) 

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
    axis_title_size = axis_title_size
)

ttest_dose_gg

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
    axis_title_size = axis_title_size
)

ttest_cell_gg

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
               scales = "free_y",
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

distrib_gg

table(data_df$Metadata_CellLine, data_df$Metadata_Dosage, data_df$Metadata_Batch_Number)

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
for (extension in c('.png', '.pdf', '.svg')) {
    ggsave(main_plot,
           filename = paste0(file_base, extension),
           height = 6,
           width = 8)
}
