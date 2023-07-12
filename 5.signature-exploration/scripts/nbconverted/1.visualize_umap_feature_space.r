suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))

# Output file
output_file <- file.path("figures", "umap_feature_comparison.png")

# Load singscore results
singscore_file <- file.path("../3.resistance-signature/results/singscore/singscore_resultsbortezomib.tsv.gz")

singscore_df <- readr::read_tsv(
    singscore_file,
    col_types = readr::cols(
        Metadata_Plate = "c",
        Metadata_Well = "c",
        Metadata_batch = "c",
        TotalScore = "d"
    )
) %>%
    dplyr::select(Metadata_Plate, Metadata_Well, Metadata_batch, Metadata_dataset, TotalScore)

print(dim(singscore_df))
head(singscore_df, 3)

# Load umap summary
umap_file <- file.path("results", "umap_feature_summary.tsv.gz")
umap_df <- readr::read_tsv(
    umap_file,
    col_types = readr::cols(
        Metadata_Plate = "c",
        Metadata_Well = "c",
        Metadata_batch = "c"
    )
    )%>%
    dplyr::inner_join(
        singscore_df,
        by = c("Metadata_Plate", "Metadata_Well", "Metadata_batch"),
        relationship = "many-to-many"
    )

umap_df$Metadata_umap_category <- dplyr::recode(
    umap_df$Metadata_umap_category,
    all_features = "All features",
    feature_selected = "Feature selected",
    all_except_bortezomib_signature_features = "All except BZ",
    bortezomib_signature_features = "BZ features"
)

umap_df$Metadata_umap_category <- factor(
    umap_df$Metadata_umap_category, levels = c("All features", "Feature selected", "All except BZ", "BZ features")
)

print(dim(umap_df))
head(umap_df, 3)

colnames(umap_df)

panel_status_gg <- (
    ggplot(
        umap_df,
        aes(x = umap_0, y = umap_1)
    )
    + geom_point(aes(color = Metadata_clone_type), alpha = 0.7)
    + facet_wrap("~Metadata_umap_category", scales="free")
    + theme_bw()
    + labs(x = "UMAP 0", y = "UMAP 1")
    + theme(
        strip.background = element_rect(
            colour = "black",
            fill = "#fdfff4"
        )
    )
    + scale_color_manual(
        name = "Status",
        labels = c("resistant" = "Resistant", "wildtype" = "Wild-type"),
        values = c("resistant" = "blue", "wildtype" = "orange")
    )
    # Decrease spacing in legend
    + theme(
        legend.spacing.y = unit(0.1, "cm"),
        legend.box.spacing = unit(0.5, "cm"),
        legend.key.size = unit(1, "lines"),
        legend.key.width = unit(1, "lines")
    )
)

panel_status_gg

panel_score_gg <- (
    ggplot(
        umap_df,
        aes(x = umap_0, y = umap_1)
    )
    + geom_point(aes(color = TotalScore))
    + facet_wrap("~Metadata_umap_category", scales="free")
    + theme_bw()
    + labs(x = "UMAP 0", y = "UMAP 1")
    + theme(
        strip.background = element_rect(
            colour = "black",
            fill = "#fdfff4"
        )
    )
    + scale_colour_gradientn(
        name = "BZ Score",
        colors = terrain.colors(10)
    )
    # Decrease spacing in legend
    + theme(
        legend.spacing.y = unit(0.1, "cm"),
        legend.box.spacing = unit(0.5, "cm"),
        legend.key.size = unit(1, "lines"),
        legend.key.width = unit(1, "lines")
    )
)

panel_score_gg

supplementary_umap_gg <- cowplot::plot_grid(
    panel_status_gg,
    panel_score_gg,
    ncol = 2,
    labels = c("A", "B")
)

ggsave(output_file, supplementary_umap_gg, height = 4, width = 10, dpi = 500)

supplementary_umap_gg
