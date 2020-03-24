suppressPackageStartupMessages(library(singscore))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

seed <- 1234
num_permutations <- 1000

sig_cols <- readr::cols(
  feature = readr::col_character(),
  estimate = readr::col_double(),
  adj.p.value = readr::col_double()
)

sig_file <- file.path("results", "cloneAE_signature_tukey.tsv")
psmb_signature_scores <- readr::read_tsv(sig_file, col_types=sig_cols)

head(psmb_signature_scores, 2)

# Extract features that are up and down in the signature
up_features <- psmb_signature_scores %>% dplyr::filter(estimate > 0) %>% dplyr::pull(feature)
down_features <- psmb_signature_scores %>% dplyr::filter(estimate < 0) %>% dplyr::pull(feature)

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

# Do not load the feature selected data
profile_dir <- file.path("..", "2.describe-data", "data", "merged")
profile_file <- file.path(profile_dir, "combined_four_clone_dataset.csv")

fourclone_data_df <- readr::read_csv(profile_file, col_types = col_types)

print(dim(fourclone_data_df))
head(fourclone_data_df, 2)

# Generate unique sample names (for downstream merging of results)
sample_names <- paste(
    fourclone_data_df$Metadata_clone_number,
    fourclone_data_df$Metadata_Plate,
    fourclone_data_df$Metadata_Well,
    fourclone_data_df$Metadata_batch,
    sep = "_"
)

fourclone_data_df <- fourclone_data_df %>%
    dplyr::mutate(Metadata_unique_sample_name = sample_names)

# Convert the four clone dataset into a feature x sample matrix without metadata
features_only_df <- t(fourclone_data_df %>% dplyr::select(!starts_with("Metadata_")))

# Apply the `rankGenes()` method to get feature rankings per feature for each sample
rankData <- rankGenes(features_only_df)
colnames(rankData) <- fourclone_data_df$Metadata_unique_sample_name

print(dim(rankData))
head(rankData, 3)

# Using the rank dataframe, up, and down features, get the sample scores
scoredf <- simpleScore(rankData, upSet = up_features, downSet = down_features)

# Merge scores with metadata features
full_result_df <- dplyr::bind_cols(
    fourclone_data_df %>% dplyr::select(starts_with("Metadata_")),
    scoredf
    )

print(dim(full_result_df))
head(full_result_df, 2)

# Generate a null distribution of scores by randomly shuffling ranks
permuteResult <- generateNull(
    upSet = up_features,
    downSet = down_features, 
    rankData = rankData,
    centerScore = TRUE,
    knownDirection = TRUE,
    B = num_permutations,
    seed = seed,
    useBPPARAM = NULL
)

# Calculate p values and add to list
pvals <- getPvals(permuteResult, scoredf)
pval_tidy <- broom::tidy(pvals)
colnames(pval_tidy) <- c("names", "Metadata_permuted_p_value")

full_result_df <- full_result_df %>%
    dplyr::left_join(
        pval_tidy,
        by = c("Metadata_unique_sample_name" = "names")
    )

# Are there differences in quantiles across batch?
batch_info <- gsub("^.*_", "", rownames(t(permuteResult)))
batch_permute <- t(permuteResult) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(batch = batch_info)

permute_bounds <- list()
for (batch_id in unique(batch_permute$batch)) {
    subset_permute <- batch_permute %>% dplyr::filter(batch == !!batch_id) %>% dplyr::select(!batch)
    min_val <- quantile(as.vector(as.matrix(subset_permute)), 0.005)
    max_val <- quantile(as.vector(as.matrix(subset_permute)), 0.995)
    permute_bounds[[batch_id]] <- c(batch_id, min_val, max_val)
}

do.call(rbind, permute_bounds)

min_val <- quantile(as.vector(as.matrix(permuteResult)), 0.05)
max_val <- quantile(as.vector(as.matrix(permuteResult)), 0.95)

apply_psmb_signature_gg <- ggplot(full_result_df,
       aes(y = TotalScore,
           x = Metadata_clone_number)) +
    geom_boxplot(aes(fill = Metadata_treatment), outlier.alpha = 0) +
    geom_point(
        aes(fill = Metadata_treatment, group = Metadata_treatment),
        position = position_dodge(width=0.75),
        size = 0.9,
        alpha = 0.7,
        shape = 21) +
    scale_fill_manual(name = "Treatment",
                      labels = c("bortezomib" = "Bortezomib", "DMSO" = "DMSO"),
                      values = c("bortezomib" = "#9e0ba3", "DMSO" = "#fcba03")) +
    theme_bw() +
    annotate("rect", ymin = min_val,
              ymax = max_val,
              xmin = 0,
              xmax = length(unique(full_result_df$Metadata_clone_number)) + 1,
              alpha = 0.2,
              color = "red",
              linetype = "dashed",
              fill = "grey") +
    xlab("") +
    ylab("PSMB5 Signature Score") +
    theme(axis.text.x = element_text(angle=90)) +
    facet_wrap("~Metadata_batch", nrow=3) +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

output_fig <- file.path("figures", "signature", "psmb5_signature_apply_fourclone.png")
ggsave(output_fig, dpi = 500, height = 5, width = 7)
apply_psmb_signature_gg

sig_file <- file.path("results", "fourclone_signature_tukey.tsv")
resistance_signature_scores <- readr::read_tsv(sig_file, col_types=sig_cols)

head(resistance_signature_scores, 2)

# Extract features that are up and down in the signature
up_resistance_features <- resistance_signature_scores %>%
    dplyr::filter(estimate > 0) %>%
    dplyr::pull(feature)

down_resistance_features <- resistance_signature_scores %>%
    dplyr::filter(estimate < 0) %>%
    dplyr::pull(feature)

# Do not load the feature selected data
profile_file <- file.path(profile_dir, "combined_cloneAcloneE_dataset.csv")

cloneae_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_CellLine = readr::col_character(),
  Metadata_Plate = readr::col_character(),
  Metadata_Well = readr::col_character(),
  Metadata_batch = readr::col_character(),
  Metadata_plate_map_name = readr::col_character(),
  Metadata_clone_type = readr::col_character()
)


cloneAE_data_df <- readr::read_csv(profile_file, col_types = cloneae_cols)

print(dim(cloneAE_data_df))
head(cloneAE_data_df, 2)

# Generate unique sample names (for downstream merging of results)
cloneae_sample_names <- paste(
    cloneAE_data_df$Metadata_CellLine,
    cloneAE_data_df$Metadata_Plate,
    cloneAE_data_df$Metadata_Well,
    cloneAE_data_df$Metadata_batch,
    sep = "_"
)

cloneAE_data_df <- cloneAE_data_df %>%
    dplyr::mutate(Metadata_unique_sample_name = cloneae_sample_names)

# Convert the four clone dataset into a feature x sample matrix without metadata
features_only_res_df <- t(cloneAE_data_df %>% dplyr::select(!starts_with("Metadata_")))

# Apply the `rankGenes()` method to get feature rankings per feature for each sample
rankData_res <- rankGenes(features_only_res_df)
colnames(rankData_res) <- cloneAE_data_df$Metadata_unique_sample_name

print(dim(rankData_res))
head(rankData_res, 3)

# Using the rank dataframe, up, and down features, get the sample scores
scoredf_res <- simpleScore(rankData_res,
                           upSet = up_resistance_features,
                           downSet = down_resistance_features)

# Merge scores with metadata features
full_res_result_df <- dplyr::bind_cols(
    cloneAE_data_df %>% dplyr::select(starts_with("Metadata_")),
    scoredf_res
    )

print(dim(full_res_result_df))
head(full_res_result_df, 2)

# Generate a null distribution of scores by randomly shuffling ranks
permuteResult_res <- generateNull(
    upSet = up_resistance_features,
    downSet = down_resistance_features, 
    rankData = rankData_res,
    centerScore = TRUE,
    knownDirection = TRUE,
    B = num_permutations,
    seed = seed,
    useBPPARAM = NULL
)

# Calculate p values and add to list
pvals_res <- getPvals(permuteResult_res, scoredf_res)
pval_res_tidy <- broom::tidy(pvals)
colnames(pval_res_tidy) <- c("names", "Metadata_permuted_p_value")

full_res_result_df <- full_res_result_df %>%
    dplyr::left_join(
        pval_res_tidy,
        by = c("Metadata_unique_sample_name" = "names")
    )

min_val <- quantile(as.vector(as.matrix(permuteResult_res)), 0.05)
max_val <- quantile(as.vector(as.matrix(permuteResult_res)), 0.95)

append_dose <- function(string) paste0("Dose: ", string, "nM")

apply_res_signature_gg <- ggplot(full_res_result_df,
       aes(y = TotalScore,
           x = Metadata_CellLine)) +
    geom_boxplot(aes(fill = Metadata_clone_type), outlier.alpha = 0) +
    geom_point(
        aes(fill = Metadata_clone_type, group = Metadata_clone_type),
        position = position_dodge(width=0.75),
        size = 0.9,
        alpha = 0.7,
        shape = 21) +
    scale_fill_manual(name = "Clone Type",
                      labels = c("resistant" = "Resistant", "wildtype" = "WildType"),
                      values = c("resistant" = "#f5b222", "wildtype" = "#4287f5")) +
    theme_bw() +
    annotate("rect", ymin = min_val,
                  ymax = max_val,
                  xmin = 0,
                  xmax = length(unique(full_res_result_df$Metadata_CellLine)) + 1,
                  alpha = 0.2,
                  color = "red",
                  linetype = "dashed",
                  fill = "grey") +
    xlab("") +
    ylab("Generic Resistance Signature Score") +
    theme(axis.text.x = element_text(angle=90)) +
    facet_grid("Metadata_Dosage~Metadata_batch",
               labeller = labeller(Metadata_Dosage = as_labeller(append_dose))) +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))


output_fig <- file.path("figures", "signature", "generic_resistance_signature_apply_cloneAE.png")
ggsave(output_fig, dpi = 500, height = 5, width = 5)
apply_res_signature_gg

full_res_result_df$Metadata_Dosage <- factor(
    full_res_result_df$Metadata_Dosage, levels = unique(sort(full_res_result_df$Metadata_Dosage))
)
full_res_result_df <- full_res_result_df %>%
    dplyr::mutate(Metadata_group = paste0(Metadata_batch, Metadata_CellLine))

ggplot(full_res_result_df, aes(x = Metadata_Dosage, y = TotalScore, color = Metadata_CellLine, group = Metadata_group)) +
    geom_point(size = 1) +
    geom_smooth(aes(fill = Metadata_clone_type), method = "loess", lwd = 0.5) +
    facet_wrap("~Metadata_batch", nrow = 2) +
    theme_bw() +
    scale_fill_manual(name = "Clone Type",
                          labels = c("resistant" = "Resistant", "wildtype" = "WildType"),
                          values = c("resistant" = "#f5b222", "wildtype" = "#4287f5")) +
    ylab("Generic Resistance Signature Score") +
    annotate("rect", ymin = min_val,
                      ymax = max_val,
                      xmin = 0,
                      xmax = length(unique(full_res_result_df$Metadata_CellLine)) + 2,
                      alpha = 0.2,
                      color = "red",
                      linetype = "dashed",
                      fill = "grey") +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

output_fig <- file.path("figures", "signature", "generic_resistance_signature_apply_cloneAE_xaxis_dosage.png")
ggsave(output_fig, dpi = 500, height = 5, width = 5)
