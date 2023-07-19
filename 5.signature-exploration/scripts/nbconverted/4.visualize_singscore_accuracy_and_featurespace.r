suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork))

# Output file
output_file <- file.path("figures", "singscore_misclassified_samples.png")

# Clones for manuscript
# We selected these clones as our initial training, test, holdout and validation sets
# for the manuscript. We selected these after a period of data collection where we
# classified proliferation under Bortezomib and other drugs to determine
# bortezomib and multi-drug resistance. We do not have proliferation values for
# clones not included in this set
select_clones <- c(
    "WT_parental",
    "CloneA",
    "CloneE",
    "WT001",
    "WT002",
    "WT003",
    "WT004",
    "WT005",
    "WT006",
    "WT007",
    "WT010",
    "WT012",
    "WT013",
    "WT014",
    "WT015",
    "BZ001",
    "BZ002",
    "BZ003",
    "BZ004",
    "BZ005",
    "BZ006",
    "BZ007",
    "BZ008",
    "BZ009",
    "BZ010"
)

# Load singscore summary
accuracy_summary_file <- file.path("results", "singscore_accuracy_summary.tsv")

summary_df <- readr::read_tsv(
    accuracy_summary_file, show_col_types = FALSE
) %>%
    dplyr::arrange(prop_high_confidence, desc(prop_inaccurate)) %>%
    dplyr::filter(Metadata_clone_number %in% select_clones)

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
    prop_completely_incorrect = "High incorrect",
    prop_high_confidence = "High correct",
    prop_accurate = "Low correct",
    prop_inaccurate = "Low incorrect"
)

summary_df$category <- factor(
    summary_df$category,
    levels = rev(c("High correct", "Low correct", "Low incorrect", "High incorrect"))
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
    levels = rev(c(
        "AreaShape",
        "Texture",
        "RadialDistribution",
        "Intensity",
        "Granularity",
        "Correlation"
    ))
)

ks_test_group_mean_df$channel_cleaned <- factor(
    ks_test_group_mean_df$channel_cleaned,
    levels = rev(c("AGP", "Nucleus", "ER", "Mito", "RNA", "other"))
)

ks_test_group_mean_df$clone_type <- dplyr::recode(
    ks_test_group_mean_df$clone_type,
    resistant = "Resistant",
    wildtype = "Wildtype",
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
    + scale_fill_viridis_c(name = "mean KS\nstatistic\n(wrong vs.\nhigh\nconfident\nsamples)", option = "magma")
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
    + labs(
        x = "Misclassified wildtype\nKS test statistic",
        y = "Misclassified resistant\nKS test statistic"
    )
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
    (
        (
            accuracy_summary_gg | plot_spacer()
        ) + plot_layout(widths = c(1, 0.25))
    )
    /
    (
        (
            misclassified_summary_gg
            + ks_test_gg
        ) + plot_layout(widths = c(1, 0.55))
    )
    )

patchwork_plot <- (
    patchwork_plot 
    + plot_annotation(tag_levels = "A")
    + plot_layout(heights = c(0.36, 1))
)

ggsave(output_file, patchwork_plot, height = 5.8, width = 10, dpi = 500)

patchwork_plot
