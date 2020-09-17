suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

col_types <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_Assay_Plate_Barcode = readr::col_character(),
    Metadata_Plate_Map_Name = readr::col_character(),
    Metadata_well_position = readr::col_character(),
    Metadata_clone_number = readr::col_character(),
    Metadata_clone_type = readr::col_character(),
    Metadata_plate_ID = readr::col_character(),
    Metadata_plate_filename = readr::col_character(),
    Metadata_treatment = readr::col_character(),
    Metadata_batch = readr::col_character()
)

profile_file <- file.path("data", "core_batch_profiles.tsv")
data_df <- readr::read_tsv(profile_file, col_types = col_types)

print(dim(data_df))
head(data_df)

cp_features <- colnames(
    data_df %>% dplyr::select(-starts_with("Metadata_"))
    )

num_cp_features <- length(cp_features)
print(num_cp_features)

all_results <- list()
for (feature in cp_features) {
    # Build formula call
    formula_call = paste(
        feature, "~",
        "Metadata_clone_type", "+",
        "Metadata_treatment", "+",
        "Metadata_batch", "+",
        "Metadata_plate_ID", "+",
        "Metadata_clone_number", "+",
        "Metadata_Well"
    )
    
    aov.out <- aov(
        formula = as.formula(formula_call),
        data = data_df
    )
    
    results <- broom::tidy(aov.out) %>%
        dplyr::mutate(feature = feature)
    
    all_results[[feature]] <- results
}

full_results_df <- do.call(rbind, all_results)

full_results_df <- full_results_df %>%
    dplyr::mutate(neg_log_p = -log10(p.value)) %>%
    tidyr::drop_na()

full_results_df$term <- factor(
    full_results_df$term, levels = c(
        "Metadata_batch",
        "Metadata_plate_ID",
        "Metadata_Well",
        "Metadata_treatment",
        "Metadata_clone_number",
        "Metadata_clone_type"
    )
)

print(dim(full_results_df))
head(full_results_df, 6)

difference_contribution_gg <- ggplot(full_results_df, aes(x = neg_log_p, y = term)) +
    ggridges::geom_density_ridges(aes(fill = term),
                                  alpha = 0.4,
                                  scale = 4) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("-log10 p value") +
    ylab("") +
    ggtitle(paste("ANOVA effects for all", num_cp_features, "CP features")) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          title = element_text(size = 8))

out_file <- file.path("figures", "anova_effect_term_distributions.png")
ggsave(out_file, dpi = 400, height = 5, width = 4)

difference_contribution_gg

feature_focus_df <- full_results_df %>%
    reshape2::dcast(feature ~ term, value.var = "neg_log_p")

head(feature_focus_df)

label_logic <- (
    (
        feature_focus_df$Metadata_treatment > 50
    ) | (
        feature_focus_df$Metadata_clone_number > 45
    )
    )

feature_interpret_gg <- ggplot(feature_focus_df,
       aes(x = Metadata_treatment,
           y = Metadata_clone_number,
           size = Metadata_Well)) +
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
                        x = Metadata_treatment,
                        y = Metadata_clone_number)) +
    scale_size_continuous(name = "Well Effect", range = c(0.1, 3)) +
    theme_bw() +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7))

out_file <- file.path("figures", "anova_effect_term_features_clone_vs_treatment.png")
ggsave(out_file, dpi = 400, height = 5, width = 5)

feature_interpret_gg

label_logic <- (
    (
        feature_focus_df$Metadata_treatment > 50
    ) | (
        feature_focus_df$Metadata_clone_type > 35
    )
    )

feature_interpret_gg <- ggplot(feature_focus_df,
       aes(x = Metadata_treatment,
           y = Metadata_clone_type,
           size = Metadata_clone_number)) +
    geom_point(alpha = 0.7) + 
    geom_text_repel(data = subset(feature_focus_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.7,
                    point.padding = 0.5,
                    segment.size = 0.5,
                    segment.alpha = 0.9,
                    size = 1.45,
                    fontface = "italic",
                    aes(label = feature,
                        x = Metadata_treatment,
                        y = Metadata_clone_type)) +
    theme_bw() +
    scale_size_continuous(name = "Sample Effect", range = c(0.1, 3)) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7))

out_file <- file.path("figures", "anova_effect_term_features_clone_type_vs_treatment.png")
ggsave(out_file, dpi = 400, height = 5, width = 5)

feature_interpret_gg

split_feature_df <- full_results_df %>%
        tidyr::separate(feature,
                        into=c("compartment",
                               "feature_group",
                               "measurement",
                               "channel", 
                               "parameter1", 
                               "parameter2"), sep="_", remove=FALSE)

head(split_feature_df)

area_split_df <- split_feature_df %>%
    dplyr::filter(feature_group  == "AreaShape") %>%
    dplyr::group_by(term, compartment, feature_group) %>%
    dplyr::top_n(n = 1, wt = neg_log_p)

group <- "Metadata_treatment"

ggplot(area_split_df %>% dplyr::filter(term == !!group),
           aes(x = compartment, y = feature_group)) +
    geom_point(aes(fill = neg_log_p), size = 10, pch = 21) +
    ggtitle(group) +
    theme_bw() +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

out_file <- file.path("figures", "anova_effect_term_features_treatment_interpret_area.png")
ggsave(out_file, dpi = 400, height = 3, width = 5)

group <- "Metadata_clone_type"

ggplot(area_split_df %>% dplyr::filter(term == !!group),
           aes(x = compartment, y = feature_group)) +
    geom_point(aes(fill = neg_log_p), size = 10, pch = 21) +
    ggtitle(group) +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

out_file <- file.path("figures", "anova_effect_term_features_clone_type_interpret_area.png")
ggsave(out_file, dpi = 400, height = 2, width = 5)

not_area_split_df <- split_feature_df %>%
    dplyr::filter(
        feature_group %in% c("Texture", "Intensity", "RadialDistribution", "Correlation", "Granularity")) %>%
    dplyr::group_by(term, compartment, feature_group, channel) %>%
    dplyr::top_n(n = 1, wt = neg_log_p)

group <- "Metadata_clone_type"

ggplot(not_area_split_df %>% dplyr::filter(term == !!group),
           aes(x = channel, y = feature_group)) +
    geom_point(aes(fill = neg_log_p), size = 10, pch = 21) +
    facet_grid(rows = vars(compartment)) +
    ggtitle(group) +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

out_file <- file.path("figures", "anova_effect_term_features_clone_type_interpret.png")
ggsave(out_file, dpi = 400, height = 6, width = 5)

group <- "Metadata_Well"

ggplot(not_area_split_df %>% dplyr::filter(term == !!group),
           aes(x = channel, y = feature_group)) +
    geom_point(aes(fill = neg_log_p), size = 10, pch = 21) +
    facet_grid(rows = vars(compartment)) +
    ggtitle(group) +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

out_file <- file.path("figures", "anova_effect_term_features_well_interpret.png")
ggsave(out_file, dpi = 400, height = 6, width = 5)

group <- "Metadata_clone_type"

ggplot(not_area_split_df %>% dplyr::filter(term == !!group),
           aes(x = channel, y = feature_group)) +
    geom_point(aes(fill = neg_log_p), size = 10, pch = 21) +
    facet_grid(rows = vars(compartment)) +
    ggtitle(group) +
    theme_bw() +
    theme(strip.text = element_text(size = 8, color = "black"),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"))

out_file <- file.path("figures", "anova_effect_term_features_clone_type_interpret.png")
ggsave(out_file, dpi = 400, height = 6, width = 5)
