suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(ggrepel))

seed <- 1234
set.seed(seed)

datasets <- c(
    "cloneAE" = "Bortezomib",
    "ixazomib" = "Ixazomib",
    "cb5083" = "CB-5083"
)

signature_dir <- file.path("results", "signatures")

anova_file <- file.path(signature_dir, "anova_results_full_bulk_signature.tsv.gz")
tukey_file <- file.path(signature_dir, "tukey_results_full_bulk_signature.tsv.gz")
cell_count_file <- file.path(signature_dir, "lm_cell_count_results_full_bulk_signature.tsv.gz")
summary_file <- file.path(signature_dir, "signature_summary_full_bulk_signature.tsv")

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
    plate_exclude = readr::col_character(),
    batch_exclude = readr::col_character(),
    non_specific_exclude = readr::col_character(),
    final_signature = readr::col_character(),
    dataset = readr::col_character()
)

summary_df <- readr::read_tsv(summary_file, col_types = summary_cols)

# Identify significance line
num_cp_features <- length(unique(anova_df$feature))
signif_line <- -log10(0.05 / num_cp_features)
signif_line

# Split features into categories
anova_df <- anova_df %>%
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
    dplyr::mutate(compartment_feature_group = paste(compartment, feature_group, sep=" - "))

head(anova_df, 2)

cell_count_df <- cell_count_df %>%
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
    dplyr::mutate(compartment_feature_group = paste(compartment, feature_group, sep=" - "))

head(cell_count_df)

# Recode the factors in plotting datasets for improved viz
recode_terms <- c(
    "Metadata_batch" = "Batch",
    "Metadata_clone_number" = "Clone ID",
    "Metadata_clone_type_indicator" = "Resistance Status",
    "scale(Metadata_cell_count)" = "Cell Count",
    "Metadata_Plate" = "Plate"
)

anova_df$term <- dplyr::recode(anova_df$term, !!!recode_terms)
cell_count_df$term <- dplyr::recode(cell_count_df$term, !!!recode_terms)
tukey_df$term <- dplyr::recode(tukey_df$term, !!!recode_terms)

# Set colors and labels
term_labels = c(
    "Batch" = "Batch",
    "Clone ID" = "Clone ID",
    "Resistance Status" = "Resistance Status",
    "Cell Count" = "Cell Count",
    "Plate" = "Plate"
)
term_values = c(
    "Batch" = "#75BE3A",
    "Clone ID" = "#EF7333",
    "Resistance Status" = "#146CB1",
    "Cell Count" = "black",
    "Plate" = "#301B65"
)

for (dataset in names(datasets)) {
    perturbation <- datasets[[dataset]]
    
    anova_subset <- anova_df %>%
        dplyr::filter(dataset == !!dataset, neg_log_p < 200)
    cell_count_subset <- cell_count_df %>%
        dplyr::filter(dataset == !!dataset, term == "Cell Count", neg_log_p < 200)
    
    feature_significance_gg <- ggplot(anova_subset, aes(x = feature, y = neg_log_p)) +
        geom_point(aes(color = term), size = 0.5, alpha = 0.5) +
        geom_point(data = cell_count_subset,
                   aes(color = term), size = 0.7, alpha = 0.5) +
        geom_hline(yintercept = signif_line, color = "red", linetype = "dashed") +
        theme_bw() +
        facet_wrap("compartment~feature_group", scales = "free_x", ncol = 7) +
        xlab("CellProfiler Feature") +
        ylab("-log10 p value") +
        ggtitle(paste("Perturbation:", perturbation)) +
        scale_color_manual("Model term", labels = term_labels, values = term_values) +
        theme(
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 5),
            legend.key.size = unit(0.4, "cm"),
            axis.title = element_text(size = 8),
            axis.text.y = element_text(size = 6),
            strip.text = element_text(size = 5),
            strip.background = element_rect(colour="black", fill="#fdfff4")
        )
    
    
    output_fig_file <- file.path(output_fig_dir, paste0("signature_linear_model_results_", dataset, ".png"))
    ggsave(output_fig_file, feature_significance_gg, dpi = 500, height = 4, width = 7)
    
    print(feature_significance_gg)
}

exclude_criteria <- colnames(summary_df %>% dplyr::select(ends_with("_exclude")))

for (dataset in names(datasets)) {
    exclude_list <- list()
    subset_summary_df <- summary_df %>% dplyr::filter(dataset == !!dataset)
    for (exclude in exclude_criteria) {
        exclude_name <- stringr::str_remove(exclude, "_exclude")

        exclude_list[[exclude_name]] <- subset_summary_df %>%
            dplyr::select(features, !!exclude) %>%
            dplyr::filter(!!as.symbol(exclude) == "TRUE") %>%
            dplyr::pull(features)
    }
    
    filename <- file.path(output_fig_dir, paste0("signature_upset_dataset_", dataset, ".pdf"))
    
    pdf(file = filename, onefile = FALSE)
    sig_size <- length(
        subset_summary_df %>% dplyr::filter(final_signature == "TRUE") %>% dplyr::pull(features)
    )
    title <- paste0("Dataset: ", dataset, "\nSignature size: ", sig_size)
    print(upset(fromList(exclude_list), order.by="freq"))
    grid::grid.text(title, x = 0.65, y = 0.95, gp=grid::gpar(fontsize = 12))
    dev.off()
}

signature_terms <- list()
tukey_point_size = 0.25
for (dataset in names(datasets)) {
    perturbation <- datasets[[dataset]]

    # Select the tukey results per dataset
    tukey_subset_df <- tukey_df %>%
        dplyr::filter(
            feature != "Cytoplasm_Correlation_Costes_RNA_Mito",
            dataset == !!dataset
        )

    # Pull the results of the cell count linear model
    cell_count_subset <- cell_count_df %>%
        dplyr::filter(
            dataset == !!dataset,
            term == "Cell Count",
            feature != "Cytoplasm_Correlation_Costes_RNA_Mito"
        )

    # Get the resistance status signature
    final_sig <- summary_df %>%
        dplyr::filter(
            dataset == !!dataset,
            final_signature == "TRUE"
        ) %>%
        dplyr::pull(features)
    
    signature_terms[[dataset]] <- final_sig

    # Remove the resistance status term - will be plotted with repel logic
    status_subset_df <- tukey_subset_df %>% dplyr::filter(term == "Resistance Status")
    tukey_no_status_subset_df <- tukey_subset_df %>% dplyr::filter(term != "Resistance Status")

    repel_logic <- (
        status_subset_df$dataset == dataset &
        status_subset_df$neg_log_adj_p > signif_line &
        status_subset_df$feature %in% final_sig
    )

    tukey_gg <- (
        ggplot(tukey_no_status_subset_df,
               aes(x = estimate, y = neg_log_adj_p)) +
        geom_point(aes(color = term), size = tukey_point_size) +
        geom_point(
            data = cell_count_subset,
            aes(x = statistic, y = neg_log_p, color = term),
            size = tukey_point_size
        ) +
        geom_point(
            data = status_subset_df,
            aes(x = estimate, y = neg_log_adj_p),
            color = ifelse(repel_logic, "red", paste(term_values["Resistance Status"])),
            size = tukey_point_size
        ) +
        geom_hline(yintercept = signif_line, linetype = "dashed", color = "red", lwd = 0.5) +
        xlab("Fold change") +
        ylab("-log 10 p value") +
        ggtitle(paste("Perturbation:", perturbation)) +
        theme_bw() +
        scale_color_manual("Model term", labels = term_labels, values = term_values) +
        facet_wrap("~term", ncol = length(unique(tukey_subset_df$term)) + 1) +
        geom_text_repel(
            data = subset(status_subset_df, repel_logic),
            arrow = arrow(length = unit(0.01, "npc")),
            size = 1.5,
            segment.size = 0.1,
            segment.alpha = 0.8,
            force = 20,
            aes(
                x = estimate,
                y = neg_log_adj_p,
                label = feature,
            )
        ) +
        theme(
            legend.position = "none",
            strip.text = element_text(size = 5),
            strip.background = element_rect(colour="black", fill="#fdfff4")
        )
    )
    
        
    output_fig_file <- file.path(output_fig_dir, paste0("tukey_foldchange_results_", dataset, ".png"))
    ggsave(output_fig_file, tukey_gg, dpi = 500, height = 3, width = 8)
    
    print(tukey_gg)
}

signature_overlap <- upset(fromList(signature_terms), order.by="freq")
title <- "Signature overlap"

filename <- file.path(output_fig_dir, paste0("signature_overlap_upset.pdf"))
    
pdf(file = filename, onefile = FALSE)
print(signature_overlap)
grid::grid.text(title, x = 0.65, y = 0.95, gp = grid::gpar(fontsize = 12))
dev.off()

print(signature_overlap)

all_signature_summary_df <- summary_df %>%
    dplyr::filter(final_signature == "TRUE") %>%
    tidyr::separate(
        features,
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
    dplyr::mutate(compartment_feature_group = paste(compartment, feature_group, sep=" - ")) %>%
    dplyr::left_join(
        tukey_df %>% dplyr::filter(term == "Resistance Status"),
        by = c("features" = "feature", "dataset" = "dataset")
    )

channels <- c("ER", "Mito", "DNA", "AGP", "RNA")

all_signature_summary_df$channel <- ifelse(
    all_signature_summary_df$channel %in% channels,
    all_signature_summary_df$channel,
    "N/A"
)

all_signature_summary_df$channel <- factor(
    all_signature_summary_df$channel, levels = c(channels, "N/A")
)

head(all_signature_summary_df)

final_signature_gg <- (
    ggplot(all_signature_summary_df, aes(x = estimate, y = neg_log_adj_p)) +
    geom_point(aes(color = dataset, shape = channel)) +
    facet_wrap("~compartment_feature_group") +
    theme_bw() +
    xlab("Fold change") +
    ylab("-log 10 p value") +
    ggtitle("Morphology signatures") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    theme(
        strip.text = element_text(size = 5),
        strip.background = element_rect(colour="black", fill="#fdfff4"),
        axis.text = element_text(size = 6)
    ) +
    geom_text_repel(
        data = all_signature_summary_df,
        arrow = arrow(length = unit(0.01, "npc")),
        size = 1.5,
        segment.size = 0.1,
        segment.alpha = 0.8,
        force = 20,
        aes(
            x = estimate,
            y = neg_log_adj_p,
            label = measurement,
        )
    )
)

output_fig_file <- file.path(output_fig_dir, paste0("signature_summary.png"))
ggsave(output_fig_file, final_signature_gg, dpi = 500, height = 6, width = 7)

final_signature_gg
