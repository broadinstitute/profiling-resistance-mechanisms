suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

source(file.path("scripts", "signature_utils.R"))

col_types <- readr::cols(
    .default = readr::col_double(),
    Metadata_CellLine = readr::col_character(),
    Metadata_Dosage = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_batch = readr::col_character(),
    Metadata_plate_map_name = readr::col_character(),
    Metadata_clone_type = readr::col_character()
)

profile_dir <- file.path("..", "2.describe-data", "data", "merged")
profile_file <- file.path(profile_dir, "combined_cloneAcloneE_dataset_feature_select.csv")

# Load Data and subset to only untreated samples
cloneAE_data_df <- readr::read_csv(profile_file, col_types = col_types) %>%
    dplyr::filter(Metadata_Dosage == 0)

print(dim(cloneAE_data_df))
head(cloneAE_data_df, 3)

table(
    cloneAE_data_df$Metadata_batch,
    cloneAE_data_df$Metadata_CellLine,
    cloneAE_data_df$Metadata_Dosage
    )

cloneAE_formula_terms <- paste(
    "~",
    "Metadata_clone_type", "+",
    "Metadata_batch", "+",
    "Metadata_Plate", "+",
    "Metadata_CellLine"
)

cloneAE_full_results <- perform_anova(cloneAE_data_df, cloneAE_formula_terms)

cloneAE_full_results_df <- cloneAE_full_results[["full_results_df"]]

cloneAE_full_results_df$term <- factor(
    cloneAE_full_results_df$term, levels = c(
        "Metadata_batch",
        "Metadata_Plate",
        "Metadata_CellLine",
        "Metadata_clone_type"
    )
)

cloneAE_full_results_df <- cloneAE_full_results_df %>%
    dplyr::arrange(desc(neg_log_p))

print(dim(cloneAE_full_results_df))
head(cloneAE_full_results_df, 6)

cloneAE_aov <- cloneAE_full_results[["aovs"]]
full_tukey_results_df <- process_tukey(cloneAE_aov, unique(cloneAE_full_results_df$feature))

print(dim(full_tukey_results_df))
head(full_tukey_results_df)

num_cp_features <- length(unique(cloneAE_full_results_df$feature))

difference_contribution_gg <- ggplot(cloneAE_full_results_df,
                                     aes(x = neg_log_p)) +
    geom_density(aes(fill = term), alpha = 0.4) +
    theme_bw() +
    xlab("-log10 p value") +
    ylab("") +
    ggtitle(
        paste("Comparing Clone A, E, and Wildtype\nANOVA effects for all", 
               num_cp_features,
              "CP features")) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          title = element_text(size = 8))

out_file <- file.path("figures", "cloneAE_anova_effect_term_distributions.png")
ggsave(out_file, dpi = 500, height = 5, width = 4)

difference_contribution_gg

difference_contribution_gg <- ggplot(cloneAE_full_results_df %>% dplyr::filter(neg_log_p < 200),
                                     aes(x = neg_log_p)) +
    geom_density(aes(fill = term), alpha = 0.4) +
    theme_bw() +
    xlab("-log10 p value") +
    ylab("") +
    ggtitle(
        paste("Comparing Clone A, E, and Wildtype\nANOVA effects for all", 
               num_cp_features,
              "CP features")) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          title = element_text(size = 8))

out_file <- file.path("figures", "cloneAE_anova_effect_term_distributions_cutoff.png")
ggsave(out_file, dpi = 500, height = 5, width = 4)

difference_contribution_gg

signif_line <- -log10(0.05 / nrow(full_tukey_results_df))

tukey_volcano_gg <- ggplot(full_tukey_results_df, aes(x = estimate, y = neg_log_p)) +
    geom_point(aes(color = comparison)) +
    geom_hline(yintercept = signif_line, color="red", linetype="dashed") +
    facet_wrap("~term", nrow=length(unique(full_tukey_results_df$term))) +
    theme_bw() +
    xlab("Statistic") +
    ylab("-log10 P Value") +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

out_file <- file.path("figures", "cloneAE_tukey_volcano.png")
ggsave(out_file, dpi = 500, height = 5, width = 6)

tukey_volcano_gg

feature_exclude_batch <- full_tukey_results_df %>%
    dplyr::filter(term == "Metadata_batch", neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

feature_exclude_cellline <- full_tukey_results_df %>%
    dplyr::filter(term == "Metadata_CellLine", neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

signature_features <- full_tukey_results_df %>%
    dplyr::filter(term == "Metadata_clone_type", neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

signature_features <- setdiff(signature_features, feature_exclude_cellline)
signature_features <- setdiff(signature_features, feature_exclude_batch)

signature_features <- sort(signature_features)

cloneAE_signature_df <- full_tukey_results_df %>%
    dplyr::filter(feature %in% signature_features,
                  term == "Metadata_clone_type")

output_file <- file.path("results", "cloneAE_signature_tukey.tsv")
cloneAE_signature_df %>% readr::write_tsv(output_file)

cloneAE_signature_df

col_types <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_plate_map_name = readr::col_character(),
    Metadata_clone_number = readr::col_character(),
    Metadata_clone_type = readr::col_character(),
    Metadata_plate_ID = readr::col_character(),
    Metadata_plate_filename = readr::col_character(),
    Metadata_treatment = readr::col_character(),
    Metadata_batch = readr::col_character()
)

profile_dir <- file.path("..", "2.describe-data", "data", "merged")
profile_file <- file.path(profile_dir, "combined_four_clone_dataset_feature_select.csv")

# Select only DMSO treatment
fourclone_data_df <- readr::read_csv(profile_file, col_types = col_types) %>%
    dplyr::filter(Metadata_treatment == "DMSO")

print(dim(fourclone_data_df))
head(fourclone_data_df)

table(
    fourclone_data_df$Metadata_clone_number,
    fourclone_data_df$Metadata_clone_type,
    fourclone_data_df$Metadata_treatment,
    fourclone_data_df$Metadata_batch
)

four_clone_formula_terms <- paste(
    "~",
    "Metadata_clone_type", "+",
    "Metadata_batch", "+",
    "Metadata_clone_number"
)

full_results <- perform_anova(fourclone_data_df, four_clone_formula_terms)

full_results_df <- full_results[["full_results_df"]]

full_results_df$term <- factor(
    full_results_df$term, levels = c(
        "Metadata_batch",
        "Metadata_clone_number",
        "Metadata_clone_type"
    )
)

full_results_df <- full_results_df %>%
    dplyr::arrange(desc(neg_log_p))

print(dim(full_results_df))
head(full_results_df, 6)

fourclone_aov <- full_results[["aovs"]]
full_tukey_results_df <- process_tukey(fourclone_aov, unique(full_results_df$feature))

print(dim(full_tukey_results_df))
head(full_tukey_results_df)

num_cp_features <- length(unique(full_results_df$feature))

difference_contribution_gg <- ggplot(full_results_df,
                                     aes(x = neg_log_p)) +
    geom_density(aes(fill = term), alpha = 0.4) +
    theme_bw() +
    xlab("-log10 p value") +
    ylab("") +
    ggtitle(
        paste("Comparing Four Wildtype and Resistant Clones\nANOVA effects for all", 
               num_cp_features,
              "CP features")) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          title = element_text(size = 8))

out_file <- file.path("figures", "fourclone_anova_effect_term_distributions.png")
ggsave(out_file, dpi = 500, height = 3.5, width = 4)

difference_contribution_gg

signif_line <- -log10(0.05 / nrow(full_tukey_results_df))

tukey_volcano_gg <- ggplot(full_tukey_results_df, aes(x = estimate, y = neg_log_p)) +
    geom_point(aes(color = comparison)) +
    geom_hline(yintercept = signif_line, color="red", linetype="dashed") +
    facet_wrap("~term", nrow=length(unique(full_tukey_results_df$term))) +
    theme_bw() +
    xlab("Statistic") +
    ylab("-log10 P Value") +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4")) +
    theme(legend.position = "none")

out_file <- file.path("figures", "fourclone_tukey_volcano.png")
ggsave(out_file, dpi = 500, height = 5, width = 4)

tukey_volcano_gg

# We want to exclude same clone-type comparisons, but not different clone type comparisons
clone_comparisons <- unique(
    full_tukey_results_df %>%
        dplyr::filter(term == "Metadata_clone_number") %>%
        dplyr::pull(comparison)
    )

same_clonetype_compare <- c(
    'BZ008-BZ001', 'BZ017-BZ001', 'BZ018-BZ001', 'BZ017-BZ008',
    'BZ018-BZ008', 'BZ018-BZ017', 'WT002-WT_parental', 'WT008-WT_parental',
    'WT009-WT_parental', 'WT011-WT_parental', 'WT008-WT002', 'WT009-WT002',
    'WT011-WT002', 'WT009-WT008', 'WT011-WT008', 'WT011-WT009'
)

setdiff(clone_comparisons, same_clonetype_compare)

feature_exclude_batch <- full_tukey_results_df %>%
    dplyr::filter(term == "Metadata_batch", neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

feature_exclude_cellline <- full_tukey_results_df %>%
    dplyr::filter(term == "Metadata_clone_number",
                  comparison %in% same_clonetype_compare,
                  neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

signature_features <- full_tukey_results_df %>%
    dplyr::filter(term == "Metadata_clone_type", neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

signature_features <- setdiff(signature_features, feature_exclude_cellline)
signature_features <- setdiff(signature_features, feature_exclude_batch)

signature_features <- sort(signature_features)

signature_df <- full_tukey_results_df %>%
    dplyr::filter(feature %in% signature_features,
                  term == "Metadata_clone_type")

output_file <- file.path("results", "fourclone_signature_tukey.tsv")
signature_df %>% readr::write_tsv(output_file)

signature_df

feature_focus_df <- full_results_df %>%
    reshape2::dcast(feature ~ term, value.var = "neg_log_p")

head(feature_focus_df)

label_logic <- (
    (
        feature_focus_df$Metadata_clone_type > 20
    ) | (
        feature_focus_df$Metadata_clone_number > 20
    )
    )

feature_interpret_gg <- ggplot(feature_focus_df,
       aes(x = Metadata_clone_type,
           y = Metadata_clone_number,
           size = Metadata_batch)) +
    geom_point(alpha = 0.7) + 
    geom_text_repel(data = subset(feature_focus_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.7,
                    point.padding = 0.2,
                    segment.size = 0.5,
                    segment.alpha = 0.9,
                    size = 1.45,
                    fontface = "italic",
                    aes(label = feature,
                        x = Metadata_clone_type,
                        y = Metadata_clone_number)) +
    scale_size_continuous(name = "Batch effect", range = c(0.1, 3)) +
    theme_bw() +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7))

out_file <- file.path("figures", "fourclone_anova_effect_term_features_clone_vs_treatment.png")
ggsave(out_file, dpi = 400, height = 5, width = 5)

feature_interpret_gg

combined_gg <- process_signature_features(
    signature_df, plot_title = "Four Clone Resistance Signature", visualize_metric="Max"
)

output_file <- file.path("figures", "fourclone_signature_feature_interpret.png")
cowplot::save_plot(output_file, combined_gg, base_height = 5, base_width = 6, dpi = 500)

combined_gg

combined_gg <- process_signature_features(
    cloneAE_signature_df, plot_title = "Clone A and E Resistance Signature", visualize_metric="Max"
)

output_file <- file.path("figures", "cloneAE_signature_feature_interpret.png")
cowplot::save_plot(output_file, combined_gg, base_height = 5, base_width = 6, dpi = 500)

combined_gg
