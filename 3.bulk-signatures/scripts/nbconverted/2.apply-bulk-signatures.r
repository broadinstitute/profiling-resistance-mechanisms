suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbeeswarm))

suppressPackageStartupMessages(library(singscore))

source(file.path("utils", "singscore_utils.R"))

seed <- 1234
num_permutations <- 1000
datasets <- c("cloneAE", "four_clone")
data_dir <- "data"
input_results_dir = file.path("results", "signatures")
output_dir <- file.path("results", "singscore")
figure_dir <- file.path("figures", "singscore")

set.seed(seed)

# Load train test status
status_file <- file.path(input_results_dir, "train_test_status.csv")
status_df <- readr::read_csv(status_file, col_types = readr::cols()) %>%
    dplyr::select(Metadata_sample_index, Metadata_signature_train_test, Metadata_dataset)

head(status_df, 3)

data_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_Plate = readr::col_character(),
  Metadata_Well = readr::col_character(),
  Metadata_batch = readr::col_character(),
  Metadata_clone_number = readr::col_character(),
  Metadata_plate_ID = readr::col_integer(),
  Metadata_plate_map_name = readr::col_character(),
  Metadata_treatment = readr::col_character(),
  Metadata_clone_type = readr::col_character(),
  Metadata_sample_index = readr::col_character()
)

dataset_dfs <- list()
for (dataset in datasets) {
    data_file <- file.path(data_dir, paste0("bulk_profiles_", dataset, ".csv.gz"))
    data_df <- readr::read_csv(data_file, col_types=data_cols)
    
    # Generate unique sample names (for downstream merging of results)
    sample_names <- paste(
        data_df$Metadata_clone_number,
        data_df$Metadata_Plate,
        data_df$Metadata_Well,
        data_df$Metadata_batch,
        sep = "_"
    )

    data_df <- data_df %>%
        dplyr::mutate(Metadata_unique_sample_name = sample_names)

    # Merge with status identifiers
    dataset_status_df <- status_df %>% dplyr::filter(Metadata_dataset == !!dataset)
    
    # Also note that I fill missing values here
    data_update_df <- data_df %>%
        dplyr::left_join(dataset_status_df, by = "Metadata_sample_index") %>%
        tidyr::replace_na(list(Metadata_signature_train_test = "test", Metadata_dataset = dataset))
    
    dataset_dfs[[dataset]] <- data_update_df
}

sig_cols <- readr::cols(
  feature = readr::col_character(),
  estimate_tukey = readr::col_double(),
  adj.p.value_tukey = readr::col_double()
)

signature_dfs <- list()
signature_features <- list()
for (dataset in datasets) {
    # Load signature
    sig_file <- file.path("results", "signatures", paste0("bulk_signature_", dataset, ".tsv"))
    signature_df <- readr::read_tsv(sig_file, col_types=sig_cols)
    
    # Extract features that are up and down in the signature
    up_features <- signature_df %>% dplyr::filter(estimate_tukey > 0) %>% dplyr::pull(feature)
    down_features <- signature_df %>% dplyr::filter(estimate_tukey < 0) %>% dplyr::pull(feature)
    
    signature_features[[dataset]] <- list("up" = up_features, "down" = down_features)
    signature_dfs[[dataset]] <- signature_df
}

sing_score_results <- list()
for (dataset in datasets) {
    data_df = dataset_dfs[[dataset]]
    sing_score_results[[dataset]] <- list()
    for (signature in datasets) {
        signature_info <- signature_features[[signature]]

        # Rank the features per sample
        rank_df <- getRankData(df = data_df)
        
        # Get the scores
        simple_score_df <- applySimpleScore(
            df = data_df,
            rank_df = rank_df,
            sig_feature_list = signature_info
        )
    
        # Permute the data to generate null scores
        permuted_output <- getPermutedRanks(
            df = data_df,
            rank_df = rank_df,
            simple_score_df = simple_score_df,
            sig_feature_list = signature_info,
            num_permutations = num_permutations
        )
        
        full_results_df <- permuted_output[["results"]]
        permuted <- permuted_output[["permuted"]]

        # Annotate some key metadata and store to list
        full_results_df <- full_results_df %>%
            dplyr::mutate(dataset = dataset, signature = signature)
        
        # Output file
        output_file = file.path(
            output_dir, paste0("data_", dataset, "_signature_", signature, ".tsv.gz")
        )
        full_results_df %>% readr::write_tsv(output_file)
        
        # Store results in a list for downstream plotting
        sing_score_results[[dataset]][[signature]] <- list(
            "results" = full_results_df, "permuted" = permuted
        )
    }
}

for (dataset in datasets) {
    for (signature in datasets) {
        output_file <- file.path(
            figure_dir, paste0("data_", dataset, "_signature_", signature, "_apply_singscore.png")
            )
        result <- sing_score_results[[dataset]][[signature]][["results"]]
        permute_result <- sing_score_results[[dataset]][[signature]][["permuted"]]

        min_val <- quantile(as.vector(as.matrix(permute_result)), 0.05)
        max_val <- quantile(as.vector(as.matrix(permute_result)), 0.95)
        
        results_gg <- ggplot(result,
           aes(y = TotalScore,
               x = Metadata_clone_number,
               group = paste(Metadata_clone_number, Metadata_treatment))) +
            facet_wrap("~Metadata_signature_train_test")
        
        if (dataset == "cloneAE") {
            results_gg <- results_gg +
                geom_boxplot(aes(fill = Metadata_treatment), outlier.alpha = 0) +
                geom_quasirandom(
                    aes(shape = Metadata_Plate),
                    dodge.width = 0.75,
                    size = 0.5) +
            theme_bw()
        } else {
            results_gg <- results_gg +
                geom_boxplot(aes(color = Metadata_treatment), outlier.alpha = 0) +
                geom_quasirandom(
                    aes(fill = Metadata_batch),
                    dodge.width = 0.75,
                    size = 0.5,
                    shape = 21) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90))
        }
        
        results_gg <- results_gg +
            annotate("rect",
                     ymin = min_val,
                     ymax = max_val,
                     xmin = 0,
                     xmax = length(unique(result$Metadata_clone_number)) + 1,
                     alpha = 0.2,
                     color = "red",
                     linetype = "dashed",
                     fill = "grey") +
            xlab("") +
            ylab("TotalScore (singscore)") +
            ggtitle(paste("Signature:", signature)) +
            theme(strip.text = element_text(size = 8, color = "black"),
                  strip.background = element_rect(colour = "black", fill = "#fdfff4"))
        
        print(results_gg)
        ggsave(output_file, height = 4, width = 8, dpi = 400)
    }
}
