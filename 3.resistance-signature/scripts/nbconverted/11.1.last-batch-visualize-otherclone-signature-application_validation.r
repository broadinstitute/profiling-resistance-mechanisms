suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gridExtra))

source(file.path("utils", "viz.R"))

seed <- 1234
set.seed(seed)

# Set constants
dataset <- "otherclones"
sig_dir <- file.path("results", "singscore")
perf_dir <- file.path("results", "performance")

sig_results_file <- file.path(sig_dir, paste0("singscore_results_LAST_BATCH_VALIDATIONotherclones.tsv.gz"))
roc_auc_file <- file.path(perf_dir, paste0(dataset, "_bortezomibsignature_roc_auc_LAST_BATCH_VALIDATION.tsv"))
roc_curve_file <- file.path(perf_dir, paste0(dataset, "_bortezomibsignature_roc_curve_LAST_BATCH_VALIDATION.tsv"))

output_figure_dir <- file.path("figures", "singscore")

# Load signature scores
# Load raw signature results
sig_cols <- readr::cols(
    .default = readr::col_character(),
    Metadata_clone_type_indicator = readr::col_integer(),
    Metadata_celltype_shorthand_from_plate_graph = readr::col_integer(),
    Metadata_date = readr::col_integer(),
    Metadata_treatment_shorthand_from_plate_graph = readr::col_integer(),
    TotalScore = readr::col_double(),
    TotalDispersion = readr::col_double(),
    UpScore = readr::col_double(),
    UpDispersion = readr::col_double(),
    DownScore = readr::col_double(),
    DownDispersion = readr::col_double(),
    Metadata_permuted_p_value = readr::col_double(),
    min_permuted_value = readr::col_double(),
    max_permuted_value = readr::col_double()
)

#inference_df <- readr::read_tsv(sig_results_file, col_types = sig_cols) %>%
    #dplyr::filter(Metadata_model_split == "inference")

signature_df <- readr::read_tsv(sig_results_file, col_types = sig_cols) %>%
    #dplyr::filter(Metadata_model_split != "inference")

print(dim(signature_df))
head(signature_df, 4)

# Load performance metrics info
metric_levels <- c("dataset")
shuffle_levels <- c("", "shuffle")

performance_data <- list()
for (metric_level in metric_levels) {
    performance_data[[metric_level]] <- list()
    for (shuffle_level in shuffle_levels) {
        if (shuffle_level == "shuffle") {
            shuffle_id <- "_shuffle"
        } else {
            shuffle_id <- ""
        }
        
        # Note I set file name here
        input_file <- paste0(metric_level, "_otherclones", shuffle_id, "_metric_performance.tsv")
        input_file <- file.path(perf_dir, input_file)
        
        performance_data[[metric_level]][[shuffle_level]] <-
            readr::read_tsv(input_file, col_types = readr::cols())
    }
}

# Load ROC info
roc_curve_cols <- readr::cols(
    fpr = readr::col_double(),
    tpr = readr::col_double(),
    threshold = readr::col_double(),
    model_split = readr::col_character(),
    shuffled = readr::col_character()
)

roc_curve_df <- readr::read_tsv(roc_curve_file, col_types=roc_curve_cols)

roc_cols <- readr::cols(
    roc_auc = readr::col_double(),
    model_split = readr::col_character(),
    shuffled = readr::col_character()
)

roc_df <- readr::read_tsv(roc_auc_file, col_types=roc_cols) %>%
    dplyr::arrange(shuffled) %>%
    dplyr::mutate(roc_auc = round(roc_auc, 3))

roc_df <- roc_df[, c("shuffled", "roc_auc")]
colnames(roc_df) <- c("AUROC\nShuffled:", "roc_auc")

#roc_df <- roc_df %>%
    #tidyr::spread(`AUROC\nShuffled:`, roc_auc, sep="") %>%
    #dplyr::arrange(desc(Split))

#roc_df$Split <- dplyr::recode(roc_df$Split, !!!legend_labels)
#roc_df

table_theme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.3)),
    colhead = list(fg_params=list(cex = 0.3)),
    padding = unit(c(1.5, 1.5), "mm")
)

table_gg <- gridExtra::tableGrob(roc_df,
                                 theme = table_theme,
                                 rows = NULL)

threshold_points_df <- roc_curve_df %>%
    dplyr::filter(shuffled == "False") %>%
    #dplyr::group_by(model_split) %>%
    dplyr::filter(abs(threshold) == min(abs(threshold))) %>%
    dplyr::ungroup()

#legend_colors <- c("#D41159")
#legend_labels <- c("WT" = "WT", "BZ" = "BZ")

roc_gg <- (
    ggplot(roc_curve_df, aes(x = fpr, y = tpr))
    + geom_line(aes(linetype = shuffled), size=0.2)
    + geom_point(data = threshold_points_df, size = 1)
    + geom_abline(intercept=0, slope=1, size=0.2, linetype="dashed", color="black")
    + scale_color_manual(name="Dataset", labels=legend_labels, values=legend_colors)
    + scale_linetype_manual(name="Shuffled", labels=linetype_labels, values=linetype_default)
    + coord_fixed()
    + xlab("False positive rate")
    + ylab("True positive rate")
    + annotation_custom(table_gg, ymin = -0.82, xmin = 0.42)
    + custom_theme
    + theme(legend.position = "bottom")
    + guides(col = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE))
)

min_val <- unique(signature_df$min_permuted_value)
max_val <- unique(signature_df$max_permuted_value)

box_plot_gg <- (
    ggplot(signature_df, aes(x = Metadata_clone_number, y = TotalScore))
    + geom_boxplot(aes(fill = Metadata_clone_type), lwd = 0.1, alpha = 0.9, outlier.alpha = 0)
    + geom_jitter(
        position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75),
        aes(fill = Metadata_clone_type),
        stroke = 0.2,
        alpha = 0.6,
        shape = 21,
        lwd = 0.5
    )
    #+ facet_wrap(
       # "~Metadata_dataset",
       # ncol = 1,
       # scales = "free_y")
    
    + xlab("")
    + ylab("Signature score\n(singscore)")
    + scale_fill_manual("Resistance status", values = resistance_status_colors, labels = resistance_status_labels)
    + scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), minor_breaks = seq(-10, 10, 0.1))
    + geom_hline(yintercept = 0, linetype = "dashed", size = 0.2, color = "black")
    + coord_flip()
    + annotate(
        "rect",
         ymin = min_val,
         ymax = max_val,
         xmin = 0.3,
         xmax = 11,
         alpha = 0.2,
         color = "red",
         linetype = "dashed",
         fill = "grey",
         size = 0.1
    )
    + custom_theme
    + guides(fill = guide_legend(reverse = TRUE))
)

output_fig_file <- file.path("figures", "newbatch_validation_roc.png")
ggsave(output_fig_file, roc_gg, dpi = 500, height = 6, width = 7)

roc_gg

output_fig_file <- file.path("figures", "newbatch_validation_boxplot.png")
ggsave(output_fig_file, box_plot_gg, dpi = 500, height = 6, width = 7)

box_plot_gg


