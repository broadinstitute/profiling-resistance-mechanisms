suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(heatmap3))

source(file.path("utils", "viz.R"))

seed <- 1234
set.seed(seed)

dataset <- "bortezomib"

signature_dir <- file.path("results", "signatures")
data_dir <- "data"

anova_file <- file.path(signature_dir, paste0("anova_results_", dataset, "_signature.tsv.gz"))
tukey_file <- file.path(signature_dir, paste0("tukey_results_", dataset, "_signature.tsv.gz"))
cell_count_file <- file.path(signature_dir, paste0("lm_cell_count_results_", dataset, "_signature.tsv.gz"))
summary_file <- file.path(signature_dir, paste0("signature_summary_", dataset, "_signature.tsv.gz"))

data_file <- file.path(data_dir, "bortezomib_signature_analytical_set.tsv.gz")

output_fig_dir <- file.path("figures", "signature_features")

# Load data
anova_cols <- readr::cols(
    term = readr::col_character(),
    df = readr::col_integer(),
    sumsq = readr::col_double(),
    meansq = readr::col_double(),
    statistic = readr::col_double(),
    p.value = readr::col_double(),
    feature = readr::col_character(),
    neg_log_p = readr::col_double(),
    dataset = readr::col_character()
)

anova_df <- readr::read_tsv(anova_file, col_types = anova_cols)

tukey_cols <- readr::cols(
    term = readr::col_character(),
    comparison = readr::col_character(),
    estimate = readr::col_double(),
    conf.low = readr::col_double(),
    conf.high = readr::col_double(),
    adj.p.value = readr::col_double(),
    feature = readr::col_character(),
    neg_log_adj_p = readr::col_double(),
    dataset = readr::col_character()
)

tukey_df <- readr::read_tsv(tukey_file, col_types = tukey_cols)

cell_count_cols <- readr::cols(
    term = readr::col_character(),
    estimate = readr::col_double(),
    std.error = readr::col_double(),
    statistic = readr::col_double(),
    p.value = readr::col_double(),
    feature = readr::col_character(),
    rsquared = readr::col_double(),
    neg_log_p = readr::col_double()
)

cell_count_df <- readr::read_tsv(cell_count_file, col_types = cell_count_cols)

summary_cols <- readr::cols(
    features = readr::col_character(),
    non_status_significant_exclude = readr::col_character(),
    cell_count_exclude = readr::col_character(),
    batch_exclude = readr::col_character(),
    non_specific_exclude = readr::col_character(),
    final_signature = readr::col_character(),
    dataset = readr::col_character()
)

summary_df <- readr::read_tsv(summary_file, col_types = summary_cols)

# Recode the factors in plotting datasets for improved viz
recode_terms <- c(
    "Metadata_batch" = "Batch",
    "Metadata_clone_number" = "Within same clone type",
    "Metadata_clone_type_indicator" = "Resistance status",
    "scale(Metadata_cell_count)" = "Cell count",
    "Metadata_treatment_time" = "Treatment time"
)

anova_df$term <- dplyr::recode(anova_df$term, !!!recode_terms)
cell_count_df$term <- dplyr::recode(cell_count_df$term, !!!recode_terms)
tukey_df$term <- dplyr::recode(tukey_df$term, !!!recode_terms)

# Set colors and labels
term_labels = c(
    "Batch" = "Batch",
    "Within same clone type" = "Within same clone type",
    "Resistance status" = "Resistance status",
    "Cell count" = "Cell count",
    "Plate" = "Plate",
    "Treatment time" = "Treatment time"
)
term_values = c(
    "Batch" = "#75BE3A",
    "Within same clone type" = "#EF7333",
    "Resistance status" = "#146CB1",
    "Cell count" = "black",
    "Plate" = "#301B65",
    "Treatment time" = "purple"
)

# Identify significance line
num_cp_features <- length(unique(anova_df$feature))
signif_line <- -log10(0.05 / num_cp_features)
signif_line

# Get the final resistance status signature features
final_sig <- summary_df %>%
    dplyr::filter(
        dataset == !!dataset,
        final_signature == "TRUE"
    ) %>%
    dplyr::pull(features)

# Split features into categories
tukey_subset_df <- tukey_df %>%
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
    dplyr::mutate(compartment_feature_group = paste(compartment, feature_group, channel, sep=" - "))

tukey_subset_df <- tukey_subset_df %>%
    dplyr::mutate(in_signature = feature %in% final_sig)

# Call infinite values the max (for plotting only)
max_val <- max(tukey_subset_df$neg_log_adj_p[!is.infinite(tukey_subset_df$neg_log_adj_p)])
tukey_subset_df[is.infinite(tukey_subset_df$neg_log_adj_p), "neg_log_adj_p"] <- max_val

# Split off the biological factor we want to isolate
status_subset_df <- tukey_subset_df %>% dplyr::filter(term == "Resistance status")
no_status_subset_df <- tukey_subset_df %>% dplyr::filter(term != "Resistance status")

# Determine if the clone number comparison is between like-clones
wt_clone_count <- stringr::str_count(
    no_status_subset_df %>%
    dplyr::filter(term == "Within same clone type") %>%
    dplyr::pull("comparison"), "WT"
)

clone_id_feature_drop <- no_status_subset_df %>%
    dplyr::filter(term == "Within same clone type") %>%
    dplyr::mutate(wt_clone_count = wt_clone_count) %>%
    dplyr::filter(neg_log_adj_p > !!signif_line, wt_clone_count != 1) %>%
    dplyr::count(feature) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::filter(n > 1) %>%
    dplyr::pull(feature)

clone_id_unique_df <- no_status_subset_df %>%
    dplyr::filter(term == "Within same clone type") %>%
    dplyr::mutate(wt_clone_count = wt_clone_count) %>%
    dplyr::filter(!(feature %in% clone_id_feature_drop))

no_status_subset_df <- no_status_subset_df %>%
    dplyr::filter(term != "Within same clone type")

head(status_subset_df)

cell_count_subset_df <- cell_count_df %>%
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
    dplyr::mutate(compartment_feature_group = paste(compartment, feature_group, channel, sep=" - ")) %>%
    dplyr::filter(term == "Cell count") %>%
    dplyr::mutate(in_signature = feature %in% final_sig)

head(cell_count_subset_df)

# Define which points to highlight
repel_logic <- (
    status_subset_df$dataset == dataset &
    status_subset_df$feature %in% final_sig
)

point_size <- 0.2
point_alpha <- 0.5
point_shape <- 1

main_effect_gg <- (
    ggplot(status_subset_df, aes(x = estimate, y = neg_log_adj_p))
    + geom_point(
        color = ifelse(repel_logic, "red", "darkgrey"),
        size = point_size,
        alpha = point_alpha,
        shape = point_shape
    )
    + geom_hline(yintercept = signif_line, linetype = "dashed", color = "red", lwd = 0.1)
    + ylim(c(0, max_val + 0.5))
    + custom_theme
    + ylab("-log10 p value")
    + xlab("Fold change")
    + facet_wrap("~term")
)

point_size = 0.05
point_alpha = 0.3
point_shape <- 1

covariate_gg <- (
    ggplot(no_status_subset_df, aes(x = estimate, y = neg_log_adj_p))
    + geom_point(
        aes(color = in_signature),
        size = point_size,
        alpha = 0,
        shape = point_shape
    )
    + geom_point(
        data = no_status_subset_df %>% dplyr::filter(!(feature %in% final_sig)),
        size = point_size,
        alpha = point_alpha,
        shape = point_shape
    )
    + geom_point(
        data = no_status_subset_df %>% dplyr::filter(feature %in% final_sig),
        color = "red",
        size = point_size,
        alpha = point_alpha,
        shape = point_shape
    )
    
    + geom_point(
        data = clone_id_unique_df %>% dplyr::filter(!(feature %in% final_sig)),
        size = point_size,
        alpha = point_alpha,
        shape = point_shape
    )
    + geom_point(
        data = clone_id_unique_df %>% dplyr::filter(feature %in% final_sig),
        color = "red",
        size = point_size,
        alpha = point_alpha,
        shape = point_shape
    )
    + geom_point(
        data = cell_count_subset_df %>% dplyr::filter(!(feature %in% final_sig)),
        aes(x = estimate, y = neg_log_p),
        size = point_size,
        alpha = point_alpha,
        shape = point_shape
    )
    + geom_point(
        data = cell_count_subset_df %>% dplyr::filter(feature %in% final_sig),
        aes(x = estimate, y = neg_log_p),
        color = "red",
        size = point_size,
        alpha = point_alpha,
        shape = point_shape
    )
    + scale_color_manual(
        name = "Bortezomib\nsensitivity\nsignature\nfeature",
        values = c("TRUE" = "red", "FALSE" = "black"),
        labels = c("TRUE" = "TRUE", "FALSE" = "FALSE")
    )
    + facet_wrap("~term", nrow = 2)
    + geom_hline(yintercept = signif_line, linetype = "dashed", color = "red", lwd = 0.1)
    + ylim(c(0, max_val + 0.5))
    + custom_theme
    + ylab("-log10 p value")
    + xlab("Fold change")
    + guides(color = guide_legend(override.aes = list(size = 1, alpha = point_alpha) ) )
)

volcano_plot_gg <- cowplot::plot_grid(
    main_effect_gg,
    covariate_gg,
    ncol = 2,
    rel_widths = c(1, 1),
    align = "hv",
    axis = "l"
)

output_fig_file <- file.path(output_fig_dir, "bortezomib_signature_feature_volcano.png")
ggsave(output_fig_file, volcano_plot_gg, dpi = 500, height = 2.5, width = 5.25)

volcano_plot_gg

final_sig_df <- status_subset_df %>%
    dplyr::filter(feature %in% final_sig) %>%
    dplyr::select(
        feature,
        compartment,
        feature_group,
        measurement,
        channel,
        estimate,
        neg_log_adj_p,
        compartment_feature_group
    )

head(final_sig_df)

sort(final_sig_df$feature)

profile_cols <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_batch = readr::col_character(),
    Metadata_cell_count = readr::col_integer(),
    Metadata_cell_density = readr::col_character(),
    Metadata_clone_number = readr::col_character(),
    Metadata_plate_map_name = readr::col_character(),
    Metadata_time_to_adhere = readr::col_character(),
    Metadata_treatment = readr::col_character(),
    Metadata_treatment_time = readr::col_character(),
    Metadata_dataset = readr::col_character(),
    Metadata_clone_type = readr::col_character(),
    Metadata_clone_type_indicator = readr::col_integer(),
    Metadata_model_split = readr::col_character(),
    Metadata_unique_sample_name = readr::col_character()
)

profiles_df <- readr::read_tsv(data_file, col_types = profile_cols) %>%
    dplyr::filter(Metadata_model_split == "validation")

feature_cor_df <- profiles_df %>%
    dplyr::select(!!final_sig_df$feature) %>%
    as.matrix() %>%
    cor()

head(feature_cor_df)

feature_metadata_df <- as.data.frame(rownames(feature_cor_df)) %>%
    dplyr::left_join(final_sig_df, keep = TRUE, by = c("rownames(feature_cor_df)" = "feature"))

full_feature_cor_df <- profiles_df %>%
    dplyr::select(!dplyr::starts_with("Metadata_")) %>%
    as.matrix() %>%
    cor()

full_feature_cor_df[is.na(full_feature_cor_df)] <- 0

output_file <- file.path(output_fig_dir, "bortezomib_signature_feature_correlation_heatmap.pdf")
pdf(output_file)
heatmap3::heatmap3(
    feature_cor_df,
    labRow = rownames(feature_cor_df),
    labCol = feature_metadata_df$compartment_feature_group,
    sym = TRUE,
    margins = c(8, 11),
    cexRow = 0.5,
    cexCol = 0.5
)
dev.off()

output_file <- file.path(output_fig_dir, "full_feature_heatmap_validation.pdf")
pdf(output_file)
heatmap3::heatmap3(
    full_feature_cor_df,
    sym = TRUE,
    labRow = NA,
    labCol = NA,
    margins = c(2, 2)
)
dev.off()
