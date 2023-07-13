suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork))

# Output file
output_file <- file.path("figures", "singscore_misclassified_samples.png")

# Load singscore summary
accuracy_summary_file <- file.path("results", "singscore_accuracy_summary.tsv")

summary_df <- readr::read_tsv(
    accuracy_summary_file, show_col_types = FALSE
) %>%
    dplyr::arrange(prop_high_confidence, desc(prop_inaccurate))

clone_number_order <- unique(summary_df$Metadata_clone_number)

# Pivot the input dataframe for plotting
summary_df <- summary_df %>%
    dplyr::select(!c(
        total_samples,
        completely_incorrect,
        high_confidence,
        accurate,
        incorrect
    )) %>%
    tidyr::pivot_longer(
        cols = c(
            "prop_completely_incorrect",
            "prop_high_confidence",
            "prop_accurate",
            "prop_inaccurate"
        ),
        names_to = "category",
        values_to = "metric"
    ) 

summary_df$Metadata_clone_number <- factor(
    summary_df$Metadata_clone_number,
    levels = clone_number_order
    )

summary_df$category <- dplyr::recode(
    summary_df$category,
    prop_completely_incorrect = "Completely wrong",
    prop_high_confidence = "High confidence",
    prop_accurate = "Accurate",
    prop_inaccurate = "Inaccurate"
)

summary_df$category <- factor(
    summary_df$category,
    levels = rev(c("High confidence", "Accurate", "Inaccurate", "Completely wrong"))
)

print(dim(summary_df))
head(summary_df, 3)

accuracy_summary_gg <- (
    ggplot(summary_df, aes(x = Metadata_clone_number, y = category))
    + geom_point(aes(fill = metric), shape = 22, size = 5)
    + theme_bw()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
    + labs(x = "Clone (no training set samples)", y = "Categories")
    + scale_fill_continuous(name = "Proportion", type = "viridis")
    + guides(fill = guide_colorbar(barwidth = 0.8, barheight = 3.5))
)

accuracy_summary_gg

ks_test_file <- file.path("results", "ks_test_misclassified_differences.tsv")
ks_test_df <- readr::read_tsv(ks_test_file, show_col_types = FALSE) %>%
    tidyr::separate(
        feature,
        into = c(
            "compartment",
            "feature_group",
            "measurement",
            "channel", 
            "parameter1", 
            "parameter2"
        ),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::mutate(channel_cleaned = channel)

ks_test_df$channel_cleaned <- dplyr::recode(ks_test_df$channel_cleaned,
    "DNA" = "Nucleus",
    "ER" = "ER",
    "AGP" = "AGP",
    "Mito" = "Mito",
    "RNA" = "RNA",
    .default = "other",
    .missing = "other"
)

print(dim(ks_test_df))
head(ks_test_df)

ks_test_group_mean_df <- ks_test_df %>%
    dplyr::group_by(compartment, feature_group, channel_cleaned, clone_type) %>%
    dplyr::summarise(mean_ks_test = mean(ks_stat))

ks_test_group_mean_df$feature_group <- factor(
    ks_test_group_mean_df$feature_group,
    levels = rev(c("AreaShape", "Texture", "RadialDistribution", "Intensity", "Granularity", "Correlation"))
)

ks_test_group_mean_df$channel_cleaned <- factor(
    ks_test_group_mean_df$channel_cleaned,
    levels = rev(c("AGP", "Nucleus", "ER", "Mito", "RNA", "other"))
)

ks_test_group_mean_df$clone_type <- dplyr::recode(
    ks_test_group_mean_df$clone_type,
    resistant = "Misclassified resistant",
    wildtype = "Misclassified wildtype",
    prop_accurate = "Accurate",
    prop_inaccurate = "Inaccurate"
)

print(dim(ks_test_group_mean_df))
head(ks_test_group_mean_df)

misclassified_summary_gg <- (
    ggplot(ks_test_group_mean_df, aes(x = channel_cleaned, y = feature_group))
    + geom_point(aes(fill = mean_ks_test), shape = 22, size = 5)
    + facet_grid("clone_type~compartment")
    + theme_bw()
    + theme(
        axis.text.x = element_text(angle = 90),
        strip.background = element_rect(
            colour = "black",
            fill = "#fdfff4"
        ),
        legend.spacing.y = unit(0.1, "cm"),
        legend.box.spacing = unit(0.5, "cm"),
        legend.key.size = unit(1, "lines"),
        legend.key.width = unit(1, "lines")
    )
    + labs(x = "Channel", y = "Feature group")
    + scale_fill_viridis_c(name = "KS\nstatistic\n(mean)", option = "magma")
)

misclassified_summary_gg

ks_test_comparison_df <- ks_test_df %>%
    dplyr::select(feature, clone_type, channel_cleaned, ks_stat) %>%
    tidyr::pivot_wider(names_from = clone_type, values_from = ks_stat)

ks_test_gg <- (
    ggplot(ks_test_comparison_df, aes(x = wildtype, y = resistant))
    + geom_point(
        aes(fill = channel_cleaned),
        shape = 21,
        size = 3,
        alpha = 0.8,
        color = "black"
    )
    + theme_bw()
    + labs(x = "Wildtype KS test statistic", y = "Resistant KS test statistic")
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue")
    + scale_fill_manual(
        name = "Channel",
        labels = c(
            "AGP" =  "AGP",
            "ER" = "ER",
            "Mito" = "Mito",
            "Nucleus" = "Nucleus",
            "other" = "other",
            "RNA" = "RNA"
        ),
        values = c(
            "AGP" =  "#332288",
            "ER" = "#117733",
            "Mito" = "#44AA99",
            "Nucleus" = "#DDCC77",
            "other" = "black",
            "RNA" = "#CC6677"
        )
    )
)

ks_test_gg

patchwork_plot <- (
    accuracy_summary_gg /
    (
        (
            misclassified_summary_gg
            + ks_test_gg
        ) + plot_layout(widths = c(1, 0.5))
    )
    )

patchwork_plot <- (
    patchwork_plot 
    + plot_annotation(tag_levels = "A")
    + plot_layout(heights = c(0.28, 1))
)

ggsave(output_file, patchwork_plot, height = 5.9, width = 10, dpi = 500)

patchwork_plot
