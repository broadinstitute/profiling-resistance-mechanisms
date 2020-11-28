suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(broom))


perform_anova <- function(df, formula_terms) {
    cp_features <- colnames(
        df %>% dplyr::select(-starts_with("Metadata_"))
        )

    all_results <- list()
    aovout_list <- list()
    for (feature in cp_features) {
        # Build formula call
        formula_call = paste(
            feature, formula_terms
        )

        aov.out <- aov(
            formula = as.formula(formula_call),
            data = df
        )
        aovout_list[[feature]] <- aov.out
        results <- broom::tidy(aov.out) %>%
            dplyr::mutate(feature = feature)

        all_results[[feature]] <- results
    }

    full_results_df <- do.call(rbind, all_results)

    full_results_df <- full_results_df %>%
        dplyr::mutate(neg_log_p = -log10(p.value)) %>%
        tidyr::drop_na()

    return_list <- list(
      "full_results_df" = full_results_df,
      "aovs" = aovout_list
  )
    return(return_list)
}


perform_linear_model <- function(df, formula_terms) {
  cp_features <- colnames(
      df %>% dplyr::select(-starts_with("Metadata_"))
      )

  lmout_results <- list()
  for (feature in cp_features) {
      # Build formula call
      formula_call = paste(
          feature, formula_terms
      )

      lm.out <- lm(
          formula = as.formula(formula_call),
          data = df
      )

      lm_summary <- summary(lm.out)
      rsquared <- lm_summary$r.squared

      results <- broom::tidy(lm_summary) %>%
          dplyr::mutate(feature = feature, rsquared = rsquared)

      lmout_results[[feature]] <- results
  }

  full_results_df <- do.call(rbind, lmout_results)

  full_results_df <- full_results_df %>%
      dplyr::mutate(neg_log_p = -log10(p.value)) %>%
      tidyr::drop_na()

  return(full_results_df)
}


process_tukey <- function(aov_list, features) {
  all_tukey_results <- list()
  for (feature in features) {
      aov.out <- aov_list[[feature]]
      tukey_results <- TukeyHSD(aov.out)
      tukey_df <- broom::tidy(tukey_results) %>%
          dplyr::mutate(feature=feature)
      all_tukey_results[[feature]] <- tukey_df
  }

  full_tukey_results_df <- do.call(rbind, all_tukey_results) %>%
    dplyr::mutate(neg_log_adj_p = -log10(adj.p.value))

  return(full_tukey_results_df)
}


process_signature_features <- function(signature_df, plot_title, visualize_metric="sum") {
  split_feature_df <- signature_df %>%
        tidyr::separate(feature,
                        into=c("compartment",
                               "feature_group",
                               "measurement",
                               "channel",
                               "parameter1",
                               "parameter2"), sep="_", remove=FALSE)

  area_split_df <- split_feature_df %>%
      dplyr::filter(feature_group  == "AreaShape") %>%
      dplyr::group_by(term, compartment, feature_group)

  not_area_split_df <- split_feature_df %>%
      dplyr::filter(
          feature_group %in% c("Texture", "Intensity", "RadialDistribution", "Correlation", "Granularity")) %>%
      dplyr::group_by(term, compartment, feature_group, channel)

  if (visualize_metric == "sum") {
    area_split_df <- area_split_df %>%
      dplyr::mutate(metric_fill = max(abs(estimate)))

    not_area_split_df <- not_area_split_df %>%
      dplyr::mutate(metric_fill = max(abs(estimate)))
  } else {
    area_split_df <- area_split_df %>%
      dplyr::mutate(metric_fill = max(abs(estimate)))

    not_area_split_df <- not_area_split_df %>%
      dplyr::mutate(metric_fill = max(abs(estimate)))
  }

  # Visualize
  area_gg <- ggplot(area_split_df,
             aes(x = compartment, y = feature_group)) +
      geom_point(aes(fill = metric_fill), size = 5, pch = 21) +
      ggtitle(plot_title) +
      ylab("") +
      theme_bw() +
      scale_fill_continuous(name = visualize_metric)

  other_feature_gg <- ggplot(not_area_split_df,
         aes(x = channel, y = feature_group)) +
      geom_point(aes(fill = metric_fill), size = 5, pch = 21) +
      facet_grid(rows = vars(compartment)) +
      theme_bw() +
      scale_fill_continuous(name = visualize_metric) +
      ylab("") +
      xlab("") +
      theme(strip.text = element_text(size = 8, color = "black"),
            strip.background = element_rect(colour = "black", fill = "#fdfff4"))

  combined_gg <- cowplot::plot_grid(
      area_gg,
      other_feature_gg,
      nrow = 2,
      rel_heights = c(0.5, 1)
  )

  return(combined_gg)

}
