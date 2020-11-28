suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(singscore))

source(file.path("utils", "singscore_utils.R"))

seed <- 1234
num_permutations <- 1000
datasets <- c(
    "cloneAE",
    "ixazomib",
    "cb5083"
)

data_dir <- "data"
data_file <- file.path(data_dir, "bulk_profiles_analytical_set.csv.gz")

input_results_dir <- file.path("results", "signatures")
signature_file <- file.path(input_results_dir, "signature_summary_full_bulk_signature.tsv")
tukey_file <- file.path(input_results_dir, "tukey_results_full_bulk_signature.tsv.gz")

output_dir <- file.path("results", "singscore")
output_results_file <- file.path(output_dir, "full_bulk_signature_singscore_results.tsv.gz")

set.seed(seed)

# Load profiles
bulk_col_types <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
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

# Load signatures
sig_col_types <- readr::cols(
    signature_features = readr::col_character(),
    plate_exclude = readr::col_logical(),
    batch_exclude = readr::col_logical(),
    non_specific_exclude = readr::col_logical(),
    final_signature = readr::col_logical(),
    dataset = readr::col_character()
)

signature_df <- readr::read_tsv(signature_file, col_types = sig_col_types)

print(dim(signature_df))
head(signature_df, 4)

# Load Tukey results (to determine if feature is "up" or "down")
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

print(dim(tukey_df))
head(tukey_df, 4)

signature_features <- list()
for (dataset in datasets) {
    # Subset data to process dataset-specific signature
    signature_subset_df <- signature_df %>%
        dplyr::filter(dataset == !!dataset, final_signature)
    
    tukey_subset_df <- tukey_df %>%
        dplyr::filter(
            dataset == !!dataset,
            term == "Metadata_clone_type_indicator",
            feature %in% signature_subset_df$features
        )
    
    # Ensure that the comparison is always resistant vs. senstive
    # and never the other way around!
    stopifnot(length(table(tukey_subset_df$comparison)) == 1)
    
    # Determine feature direction
    up_features <- tukey_subset_df %>% dplyr::filter(estimate > 0) %>% dplyr::pull(feature)
    down_features <- tukey_subset_df %>% dplyr::filter(estimate < 0) %>% dplyr::pull(feature)
    
    # Store signature for downstream analyses
    signature_features[[dataset]] <- list("up" = up_features, "down" = down_features)
}

# Print signature size info
for (dataset in datasets) {
    print(dataset)
    print(length(signature_features[[dataset]][["up"]]))
    print(length(signature_features[[dataset]][["down"]]))
}

sing_score_results <- list()
for (dataset in datasets) {
    
    bulk_subset_df <- bulk_df %>% dplyr::filter(Metadata_dataset == !!dataset)
    sing_score_results[[dataset]] <- list()
    
    for (signature in datasets) {
        signature_info <- signature_features[[signature]]

        singscore_output = singscorePipeline(
            df = bulk_subset_df,
            sig_feature_list = signature_info,
            num_permutations = num_permutations
        )
        
        full_results_df <- singscore_output[["results"]]
        permuted <- singscore_output[["permuted"]]

        # Get max and minimum values of permutation results
        min_val <- quantile(as.vector(as.matrix(permuted)), 0.05)
        max_val <- quantile(as.vector(as.matrix(permuted)), 0.95)
        
        # Annotate some key metadata and store to list
        full_results_df <- full_results_df %>%
            dplyr::mutate(
                dataset = dataset,
                signature = signature,
                min_permuted_value = min_val,
                max_permuted_value = max_val
            )
        
        # Store results
        sing_score_results[[paste0(dataset, "-", signature)]] <- full_results_df
    }
}

all_singscore_results_df <- dplyr::bind_rows(sing_score_results)

table(
    all_singscore_results_df$dataset,
    all_singscore_results_df$signature
)

all_singscore_results_df %>% readr::write_tsv(output_results_file)

print(dim(all_singscore_results_df))
head(all_singscore_results_df)
