suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(singscore))

singscorePipeline <- function(
  df, sig_feature_list, num_permutations, permute = TRUE
) {
  # Rank the features per sample
  rank_df <- getRankData(df = df)

  # Get the scores
  simple_score_df <- applySimpleScore(
      df = df,
      rank_df = rank_df,
      sig_feature_list = sig_feature_list
  )

  # Permute the data to generate null scores
  full_output <- getPermutedRanks(
      df = df,
      rank_df = rank_df,
      simple_score_df = simple_score_df,
      sig_feature_list = sig_feature_list,
      num_permutations = num_permutations
  )
  return(full_output)
}


getRankData <- function(df) {
    # Convert the four clone dataset into a feature x sample matrix without metadata
    features_only_df <- t(df %>% dplyr::select(!starts_with("Metadata_")))

    # Apply the `rankGenes()` method to get feature rankings per feature for each sample
    rankData <- rankGenes(features_only_df)
    colnames(rankData) <- df$Metadata_unique_sample_name

    return(rankData)
}

applySimpleScore <- function(df, rank_df, sig_feature_list) {
    scoredf <- simpleScore(
        rank_df,
        upSet = sig_feature_list[["up"]],
        downSet = sig_feature_list[["down"]]
    )

    return(scoredf)
}

getPermutedRanks <- function(
    df, rank_df, simple_score_df, sig_feature_list, num_permutations
) {
    permuteResult <- generateNull(
        upSet = sig_feature_list[["up"]],
        downSet = sig_feature_list[["down"]],
        rankData = rank_df,
        centerScore = TRUE,
        knownDirection = TRUE,
        B = num_permutations,
        seed = seed,
        useBPPARAM = NULL
    )

    # Merge scores with metadata features
    annotated_scores_df <- dplyr::bind_cols(
        df %>% dplyr::select(starts_with("Metadata_")),
        simple_score_df
    )

    # Calculate p values and add to list
    pvals <- getPvals(permuteResult, simple_score_df)
    pval_tidy <- broom::tidy(pvals)
    colnames(pval_tidy) <- c("names", "Metadata_permuted_p_value")

    full_result_df <- annotated_scores_df %>%
        dplyr::left_join(
            pval_tidy,
            by = c("Metadata_unique_sample_name" = "names")
        )

    return(list("results" = full_result_df, "permuted" = permuteResult))
}
