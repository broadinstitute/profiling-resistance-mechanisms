suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

source(file.path("../3.bulk-signatures/scripts", "signature_utils.R"))

set.seed(123)

dataset <- "bortezomib"
input_data_dir <- "data"
data_file <- file.path(input_data_dir, "bortezomib_signature_analytical_set.tsv.gz")
feat_file <- file.path(input_data_dir, "dataset_features_selected.tsv")

output_fig_dir = file.path("figures", "anova")
output_results_dir = file.path("results", "signatures")

alpha = 0.05

# Load profiles
bulk_col_types <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_cell_count = readr::col_integer(),
    Metadata_batch = readr::col_character(),
    Metadata_clone_number = readr::col_character(),
    Metadata_plate_map_name = readr::col_character(),
    Metadata_treatment = readr::col_character(),
    Metadata_dataset = readr::col_character(),
    Metadata_clone_type = readr::col_character(),
    Metadata_clone_type_indicator = readr::col_character(),
    Metadata_model_split = readr::col_character(),
    Metadata_cell_density = readr::col_character(),
    Metadata_treatment_time = readr::col_character(),
    Metadata_unique_sample_name = readr::col_character(),
    Metadata_time_to_adhere = readr::col_character()
)

data_df <- readr::read_tsv(data_file, col_types = bulk_col_types)

print(dim(data_df))
head(data_df, 4)

# Load feature selected features
all_selected_features_df <- readr::read_tsv(feat_file, col_types = readr::cols())
head(all_selected_features_df, 3)

# Subset dataset
bulk_subset_df <- data_df %>%
    dplyr::filter(
        Metadata_dataset == !!dataset,
        Metadata_model_split == "training",
    )

# Apply feature selection performed in 0.compile-bulk-datasets
selected_features <- all_selected_features_df %>%
    dplyr::filter(dataset == !!dataset) %>%
    dplyr::pull(features)

bulk_subset_df <- bulk_subset_df %>%
    dplyr::select(starts_with("Metadata"), all_of(selected_features))

# Populate the list for signature building
bulk_subset_df$Metadata_clone_type_indicator <- factor(
    bulk_subset_df$Metadata_clone_type_indicator, levels = c("0", "1")
)

# Print dataset description
print(paste("Training dataset:", dataset))
print(table(
    bulk_subset_df$Metadata_clone_number,
    bulk_subset_df$Metadata_batch
))

formula_terms <- paste(
    "~",
    "Metadata_clone_type_indicator", "+",
    "Metadata_batch", "+",
    "Metadata_treatment_time", "+",
    "Metadata_clone_number"
)

cell_count_formula <- paste(
    "~",
    "Metadata_clone_type_indicator", "+",
    "scale(Metadata_cell_count)"
)

# Fit ANOVA to determine sources of variation and process results
lm_results <- perform_anova(bulk_subset_df, formula_terms)

# Order the full results data frame by significance and extract feature names
full_results_df <- lm_results[["full_results_df"]] %>%
    dplyr::arrange(desc(neg_log_p))

features <- unique(full_results_df$feature)

# Perform TukeyHSD posthoc test
tukey_results <- process_tukey(
    aov_list = lm_results[["aovs"]],
    features = features
)

# Fit a linear model on cell counts
cell_count_results <- perform_linear_model(bulk_subset_df, cell_count_formula) %>%
    dplyr::arrange(desc(neg_log_p)) %>%
    dplyr::mutate(dataset = dataset)

# Isolate features that represent significant differences between resistant and senstive clones
anova_results_df <- lm_results[["full_results_df"]] %>%
    dplyr::mutate(dataset = dataset)

tukey_results_df <- tukey_results %>%
    dplyr::mutate(dataset = dataset)

features <- unique(anova_results_df$feature)

num_cp_features <- length(features)
signif_line <- -log10(alpha / num_cp_features)

# Derive signature by systematically removing features influenced by technical artifacts
signature_features <- tukey_results_df %>%
    dplyr::filter(term == "Metadata_clone_type_indicator", neg_log_adj_p > !!signif_line) %>%
    dplyr::pull(feature)

# Determine if the clone number comparison is between like-clones
wt_clone_count <- stringr::str_count(
    tukey_results_df %>%
    dplyr::filter(term == "Metadata_clone_number") %>%
    dplyr::pull("comparison"), "WT"
)

# Exclude features with very high within sensitivity-type clones
feature_exclude_nonspecific_variation <- unique(
    tukey_results_df %>%
    dplyr::filter(term == "Metadata_clone_number") %>%
    dplyr::mutate(wt_clone_count = wt_clone_count) %>%
    dplyr::filter(neg_log_adj_p > !!signif_line, wt_clone_count != 1) %>%
    dplyr::pull(feature)
    )

# Exclude features that are significantly different as explained by batch
feature_exclude_batch <- tukey_results_df %>%
    dplyr::filter(term == "Metadata_batch", neg_log_adj_p > !!signif_line) %>%
    dplyr::pull(feature)

# Exclude features that are significantly impacted by cell count
feature_exclude_cell_count <- cell_count_results %>%
    dplyr::filter(term == "scale(Metadata_cell_count)", neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

# Exclude features that are significantly impacted by cell count
feature_exclude_time <- cell_count_results %>%
    dplyr::filter(term == "Metadata_treatment_time", neg_log_p > !!signif_line) %>%
    dplyr::pull(feature)

# Restrict signature
final_signature_features <- setdiff(
    signature_features, unique(feature_exclude_cell_count)
)
final_signature_features <- setdiff(
    final_signature_features, unique(feature_exclude_nonspecific_variation)
)
final_signature_features <- setdiff(
    final_signature_features, unique(feature_exclude_batch)
)
final_signature_features <- setdiff(
    final_signature_features, unique(feature_exclude_cell_count)
)

# Create a summary of the signatures
signature_summary_df <- tibble(features)

signature_summary_df <- signature_summary_df %>%
    dplyr::mutate(
        non_status_significant_exclude = !(signature_summary_df$features %in% signature_features),
        batch_exclude = signature_summary_df$features %in% feature_exclude_batch,
        cell_count_exclude = signature_summary_df$features %in% feature_exclude_cell_count,
        non_specific_exclude = signature_summary_df$features %in% feature_exclude_nonspecific_variation,
        treatment_time_exclude = signature_summary_df$features %in% feature_exclude_time,
        final_signature = signature_summary_df$features %in% final_signature_features,
        dataset = dataset
    )

print(paste("For the dataset:", dataset))
print(paste("the number of features in the core signature:", sum(signature_summary_df$final_signature)))

# Determine feature direction
final_signature <- signature_summary_df %>% dplyr::filter(final_signature)

tukey_subset_results_df <- tukey_results_df %>%
    dplyr::filter(
        dataset == !!dataset,
        term == "Metadata_clone_type_indicator",
        feature %in% final_signature$features
    )

up_features <- tukey_subset_results_df %>% dplyr::filter(estimate > 0) %>% dplyr::pull(feature)
down_features <- tukey_subset_results_df %>% dplyr::filter(estimate < 0) %>% dplyr::pull(feature)

# Store signature for downstream analyses
signature_features <- list("up" = up_features, "down" = down_features)

signature_features

anova_output_file <- file.path(output_results_dir, paste0("anova_results_", dataset, "_signature.tsv.gz"))
tukey_output_file <- file.path(output_results_dir, paste0("tukey_results_", dataset, "_signature.tsv.gz"))
cell_count_output_file <- file.path(output_results_dir, paste0("lm_cell_count_results_", dataset, "_signature.tsv.gz"))
signature_output_file <- file.path(output_results_dir, paste0("signature_summary_", dataset, "_signature.tsv.gz"))

anova_results_df %>% readr::write_tsv(anova_output_file)
tukey_results_df %>% readr::write_tsv(tukey_output_file)
cell_count_results %>% readr::write_tsv(cell_count_output_file)
signature_summary_df %>% readr::write_tsv(signature_output_file)
