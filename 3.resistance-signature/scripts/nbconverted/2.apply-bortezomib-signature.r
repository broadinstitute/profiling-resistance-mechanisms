suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(singscore))

source(file.path("../3.bulk-signatures/utils", "singscore_utils.R"))

seed <- 1234
num_permutations <- 1000
dataset <- "bortezomib"

data_dir <- "data"
input_results_dir <- file.path("results", "signatures")
output_dir <- file.path("results", "singscore")

data_file <- file.path(data_dir, paste0(dataset, "_signature_analytical_set.tsv.gz"))
feat_file <- file.path(data_dir, "dataset_features_selected.tsv")
signature_file <- file.path(input_results_dir, paste0("signature_summary_", dataset, "_signature.tsv.gz"))
tukey_file <- file.path(input_results_dir, paste0("tukey_results_", dataset, "_signature.tsv.gz"))
output_results_file <- file.path(output_dir, paste0("singscore_results", dataset, ".tsv.gz"))

set.seed(seed)

# Load feature selected features
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
    Metadata_treatment_time = readr::col_character(),
    Metadata_unique_sample_name = readr::col_character(),
    Metadata_time_to_adhere = readr::col_character()
)

data_df <- readr::read_tsv(data_file, col_types = bulk_col_types)

# Apply feature selection performed in 0.compile-bulk-datasets
selected_features <- all_selected_features_df %>%
    dplyr::filter(dataset == !!dataset) %>%
    dplyr::pull(features)

data_df <- data_df %>%
    dplyr::select(starts_with("Metadata"), all_of(selected_features))

print(dim(data_df))
head(data_df, 4)

# Load signatures
sig_col_types <- readr::cols(
    features = readr::col_character(),
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

# Subset data to process dataset-specific signature
signature_subset_df <- signature_df %>%
    dplyr::filter(dataset == !!dataset, final_signature)

tukey_subset_df <- tukey_df %>%
    dplyr::filter(
        dataset == !!dataset,
        term == "Metadata_clone_type_indicator",
        feature %in% signature_subset_df$features
    ) %>%
    dplyr::arrange(desc(estimate))

# Ensure that the comparison is always resistant vs. senstive
# and never the other way around!
stopifnot(length(table(tukey_subset_df$comparison)) == 1)

# Determine feature direction
up_features <- tukey_subset_df %>% dplyr::filter(estimate > 0) %>% dplyr::pull(feature)
down_features <- tukey_subset_df %>% dplyr::filter(estimate < 0) %>% dplyr::pull(feature)

# Store signature for downstream analyses
signature_features <- list("up" = up_features, "down" = down_features)

signature_features

singscore_output = singscorePipeline(
    df = data_df,
    sig_feature_list = signature_features,
    num_permutations = num_permutations
)

full_results_df <- singscore_output[["results"]]
permuted <- singscore_output[["permuted"]]

# Get max and minimum values of permutation results
min_val <- quantile(as.vector(as.matrix(permuted)), 0.05)
max_val <- quantile(as.vector(as.matrix(permuted)), 0.95)

# Annotate some key metadata and store to list
sing_score_results_df <- full_results_df %>%
    dplyr::mutate(
        dataset = dataset,
        min_permuted_value = min_val,
        max_permuted_value = max_val
    )

sing_score_results_df %>% readr::write_tsv(output_results_file)

print(dim(sing_score_results_df))
head(sing_score_results_df)
