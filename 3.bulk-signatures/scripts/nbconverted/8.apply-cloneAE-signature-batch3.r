suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(singscore))

source(file.path("utils", "singscore_utils.R"))

set.seed(863)

# Load batch 3 normalized profiles
batch <- "2019_06_25_Batch3"
batch3_file <- file.path("data", paste0(batch, "_combined_normalized.csv.gz"))

profile_cols <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_batch = readr::col_character(),
    Metadata_cell_count = readr::col_integer(),
    Metadata_clone_number = readr::col_character(),
    Metadata_unique_sample_name = readr::col_character(),
    Metadata_plate_map_name = readr::col_character(),
    Metadata_treatment = readr::col_character()
)

df <- readr::read_csv(batch3_file, col_types = profile_cols)

print(dim(df))
head(df)

# Load signature
sig_dir <- file.path("results", "signatures")
signature_file <- file.path(sig_dir, "signature_summary_full_bulk_signature.tsv")

sig_col_types <- readr::cols(
    features = readr::col_character(),
    plate_exclude = readr::col_logical(),
    batch_exclude = readr::col_logical(),
    non_specific_exclude = readr::col_logical(),
    final_signature = readr::col_logical(),
    dataset = readr::col_character()
)

signature_df <- readr::read_tsv(signature_file, col_types = sig_col_types) %>%
    dplyr::filter(dataset == "cloneAE", final_signature)

print(dim(signature_df))
head(signature_df, 4)

# Load Tukey results (to determine if feature is "up" or "down")
tukey_file <- file.path(sig_dir, "tukey_results_full_bulk_signature.tsv.gz")

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

tukey_df <- readr::read_tsv(tukey_file, col_types = tukey_cols) %>%
    dplyr::filter(dataset == "cloneAE")

print(dim(tukey_df))
head(tukey_df, 4)

tukey_subset_df <- tukey_df %>%
    dplyr::filter(
        term == "Metadata_clone_type_indicator",
        feature %in% signature_df$features
)

# Determine feature direction
up_features <- tukey_subset_df %>% dplyr::filter(estimate > 0) %>% dplyr::pull(feature)
down_features <- tukey_subset_df %>% dplyr::filter(estimate < 0) %>% dplyr::pull(feature)
    
# Store signature for downstream analyses
cloneAE_signature <- list("up" = up_features, "down" = down_features)

cloneAE_signature

seed <- 1234

singscore_output = singscorePipeline(
    df = df,
    sig_feature_list = cloneAE_signature,
    num_permutations = 1000
)

full_results_df <- singscore_output[["results"]]

print(dim(full_results_df))
head(full_results_df)

min_val <- quantile(as.vector(as.matrix(singscore_output[["permuted"]])), 0.05)
max_val <- quantile(as.vector(as.matrix(singscore_output[["permuted"]])), 0.95)

append_plate <- function(string) paste0("Plate: ", string)

batch3_gg <- (
    ggplot(full_results_df, aes(y = TotalScore, x = Metadata_clone_number)) +
    geom_boxplot(lwd = 0.3) +
    geom_jitter(size = 0.2) +
    facet_wrap(
        "~Metadata_Plate",
        scales = "free_x",
        labeller = labeller(Metadata_Plate = as_labeller(append_plate))
    ) +
    xlab("") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, size = 6),
        strip.text = element_text(size = 7),
        strip.background = element_rect(colour="black", fill="#fdfff4"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.4, "cm"),
        axis.title = element_text(size = 8),
    ) +
    annotate(
            "rect",
             ymin = min_val,
             ymax = max_val,
             xmin = 0,
             xmax = 20,
             alpha = 0.2,
             color = "red",
             linetype = "dashed",
             fill = "grey",
             lwd = 0.25
        )
)

output_fig_file <- file.path("figures", "batch3", paste0("singscore_validation_", batch, ".png"))
ggsave(output_fig_file, batch3_gg, dpi = 500, height = 2.5, width = 4)

batch3_gg
