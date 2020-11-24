suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

source(file.path("scripts", "signature_utils.R"))

set.seed(123)

datasets <- c(
    "cloneAE",
    "ixazomib",
    "cb5083"
)

input_data_dir <- "data"
data_file <- file.path(input_data_dir, "bulk_profiles_analytical_set.csv.gz")

output_fig_dir = file.path("figures", "anova")
output_results_dir = file.path("results", "signatures")

# Load feature selected features
feat_file <- file.path(input_data_dir, "dataset_features_selected.tsv")
all_selected_features_df <- readr::read_tsv(feat_file, col_types = readr::cols())
head(all_selected_features_df, 3)

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
    Metadata_plate_filename = readr::col_character(),
    Metadata_treatment_time = readr::col_character(),
    Metadata_unique_sample_name = readr::col_character(),
    Metadata_time_to_adhere = readr::col_character()
)

bulk_df <- readr::read_csv(data_file, col_types = bulk_col_types)

print(dim(bulk_df))
head(bulk_df, 4)

training_data <- list()
for (dataset in datasets) {
    # Subset dataset
    bulk_subset_df <- bulk_df %>%
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
    training_data[[dataset]] <- bulk_subset_df
    
    # Print dataset description
    print(paste("Training dataset:", dataset))
    print(table(
        bulk_subset_df$Metadata_clone_number,
        bulk_subset_df$Metadata_batch
    ))
}

formula_terms <- paste(
    "~",
    "Metadata_clone_type_indicator", "+",
    "Metadata_batch", "+",
    "Metadata_Plate", "+",
    "Metadata_clone_number"
)

cell_count_formula <- paste(
    "~",
    "Metadata_clone_type_indicator", "+",
    "scale(Metadata_cell_count)"
)

# Fit two models:
# 1) ANOVA using categorical covariates followed by TukeyHSD posthoc test
# 2) Linear model using cell count as a continuous variable
lm_results <- list()
cell_count_results <- list()
tukey_results <- list()
for (dataset in datasets) {
    print(paste("Now processing...", dataset))
    
    # Extract the dataset used to train
    analytical_df <- training_data[[dataset]]
    
    # Fit linear model to determine sources of variation and process results
    lm_results[[dataset]] <- perform_anova(analytical_df, formula_terms)
    
    # Order the full results data frame by significance and extract feature names
    full_results_df <- lm_results[[dataset]][["full_results_df"]] %>%
        dplyr::arrange(desc(neg_log_p))
    
    features <- unique(full_results_df$feature)
    
    # Perform TukeyHSD posthoc test
    tukey_results[[dataset]] <- process_tukey(
        aov_list = lm_results[[dataset]][["aovs"]],
        features = features
    )
    
    # Fit a linear model on cell counts
    cell_count_results[[dataset]] <- perform_linear_model(analytical_df, cell_count_formula) %>%
        dplyr::arrange(desc(neg_log_p)) %>%
        dplyr::mutate(dataset = dataset)
}

all_anova_results <- list()
all_tukey_results <- list()
all_signature_results <- list()
for (dataset in datasets) {
    # Process ANOVA results
    anova_results_df <- lm_results[[dataset]][["full_results_df"]] %>%
        dplyr::mutate(dataset = dataset)

    all_anova_results[[dataset]] <- anova_results_df
    
    # Process tukey results
    tukey_results_df <- tukey_results[[dataset]] %>%
        dplyr::mutate(dataset = dataset)
    
    all_tukey_results[[dataset]] <- tukey_results_df
    
    # Build signature
    features <- unique(anova_results_df$feature)

    # Note that TukeyHSD() p value is already adjusted for multiple within comparisons,
    # but not across multiple features
    num_cp_features <- length(features)
    signif_line <- -log10(0.05 / num_cp_features)

    # Derive signature by systematically removing features influenced by technical artifacts
    signature_features <- tukey_results_df %>%
        dplyr::filter(term == "Metadata_clone_type_indicator", neg_log_adj_p > !!signif_line) %>%
        dplyr::pull(feature)

    feature_exclude_plate <- tukey_results_df %>%
        dplyr::filter(term == "Metadata_Plate", neg_log_adj_p > !!signif_line) %>%
        dplyr::pull(feature)
    
    feature_exclude_batch <- tukey_results_df %>%
        dplyr::filter(term == "Metadata_batch", neg_log_adj_p > !!signif_line) %>%
        dplyr::pull(feature)

    # Determine if the clone number comparison is between like-clones
    wt_clone_count <- stringr::str_count(
        tukey_results_df %>%
        dplyr::filter(term == "Metadata_clone_number") %>%
        dplyr::pull("comparison"), "WT"
    )

    # Exclude features with very high within sensitivity-type clones
    feature_exclude_nonspecific_variation <- tukey_results_df %>%
        dplyr::filter(term == "Metadata_clone_number") %>%
        dplyr::mutate(wt_clone_count = wt_clone_count) %>%
        dplyr::filter(neg_log_adj_p > !!signif_line, wt_clone_count != 1) %>%
        dplyr::pull(feature)
    
    # Exclude features that are significantly impacted by cell count
    feature_exclude_cell_count <- cell_count_results[[dataset]] %>%
        dplyr::filter(term == "scale(Metadata_cell_count)", neg_log_p > !!signif_line) %>%
        dplyr::pull(feature)

    final_signature_features <- setdiff(
        signature_features, unique(feature_exclude_cell_count)
    )
    final_signature_features <- setdiff(
        final_signature_features, unique(feature_exclude_plate)
    )
    final_signature_features <- setdiff(
        final_signature_features, unique(feature_exclude_batch)
    )
    final_signature_features <- setdiff(
        final_signature_features, unique(feature_exclude_nonspecific_variation)
    )
    
    # Create a summary of the signatures
    signature_summary_df <- tibble(features)

    signature_summary_df <- signature_summary_df %>%
        dplyr::mutate(
            non_status_significant_exclude = !(signature_summary_df$features %in% signature_features),
            cell_count_exclude = signature_summary_df$features %in% feature_exclude_cell_count,
            plate_exclude = signature_summary_df$features %in% feature_exclude_plate,
            batch_exclude = signature_summary_df$features %in% feature_exclude_batch,
            non_specific_exclude = signature_summary_df$features %in% feature_exclude_nonspecific_variation,
            final_signature = signature_summary_df$features %in% final_signature_features,
            dataset = dataset
        )
    
    print(paste("For the dataset:", dataset))
    print(paste("the number of features in the core signature:", sum(signature_summary_df$final_signature)))
    
    all_signature_results[[dataset]] <- signature_summary_df
}

# Output files
anova_output_file <- file.path(output_results_dir, "anova_results_full_bulk_signature.tsv.gz")
tukey_output_file <- file.path(output_results_dir, "tukey_results_full_bulk_signature.tsv.gz")
cell_count_output_file <- file.path(output_results_dir, "lm_cell_count_results_full_bulk_signature.tsv.gz")
signature_output_file <- file.path(output_results_dir, "signature_summary_full_bulk_signature.tsv")

dplyr::bind_rows(all_anova_results) %>% readr::write_tsv(anova_output_file)
dplyr::bind_rows(all_tukey_results) %>% readr::write_tsv(tukey_output_file)
dplyr::bind_rows(cell_count_results) %>% readr::write_tsv(cell_count_output_file)
dplyr::bind_rows(all_signature_results) %>% readr::write_tsv(signature_output_file)
