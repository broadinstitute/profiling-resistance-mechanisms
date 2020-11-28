suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(devtools))

# The correct version of ComplexHeatmap was not available on anaconda
# https://anaconda.org/bioconda/bioconductor-complexheatmap
# Install it here
if (!("ComplexHeatmap" %in% rownames(installed.packages()))) {
    install_github("jokergoo/ComplexHeatmap@a387b860186be1d09249128be1ff46d13101e45d")
}

suppressPackageStartupMessages(library(ComplexHeatmap))

# Set constants
datasets <- c("cloneAE", "ixazomib", "cb5083")
data_dir <- "data"

output_dir <- file.path("figures", "heatmaps")

lgd_title_fontsize = 9
lgd_label_fontsize = 6.5
anno_name_height = 0.45
legend_scale_cols = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Load profiles
dataset_file <- file.path(data_dir, "bulk_profiles_analytical_set.csv.gz")

data_cols <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_batch = readr::col_character(),
    Metadata_cell_count = readr::col_integer(),
    Metadata_clone_number = readr::col_character(),
    Metadata_plate_map_name = readr::col_character(),
    Metadata_treatment = readr::col_character(),
    Metadata_dataset = readr::col_character(),
    Metadata_clone_type = readr::col_character(),
    Metadata_clone_type_indicator = readr::col_integer(),
    Metadata_model_split = readr::col_character(),
    Metadata_cell_density = readr::col_character(),
    Metadata_plate_filename = readr::col_character(),
    Metadata_treatment_time = readr::col_character(),
    Metadata_unique_sample_name = readr::col_character(),
    Metadata_time_to_adhere = readr::col_character()
)

dataset_df <- readr::read_csv(dataset_file, col_types=data_cols)

print(dim(dataset_df))
head(dataset_df, 3)

# Load feature selection results
feat_file <- file.path(data_dir, "dataset_features_selected.tsv")

feat_cols <- readr::cols(
  features = readr::col_character(),
  dataset = readr::col_character()
)

features_df <- readr::read_tsv(feat_file, col_types = feat_cols)

print(dim(features_df))
head(features_df, 3)

# Load signature scores
score_file <- file.path("results", "singscore", "full_bulk_signature_singscore_results.tsv.gz")

score_cols <- readr::cols(
    .default = readr::col_character(),
    Metadata_cell_count = readr::col_integer(),
    Metadata_clone_type_indicator = readr::col_integer(),
    Metadata_plate_ID = readr::col_integer(),
    Metadata_celltype_shorthand_from_plate_graph = readr::col_integer(),
    Metadata_date = readr::col_integer(),
    Metadata_treatment_shorthand_from_plate_graph = readr::col_integer(),
    TotalScore = readr::col_double(),
    TotalDispersion = readr::col_double(),
    UpScore = readr::col_double(),
    UpDispersion = readr::col_double(),
    DownScore = readr::col_double(),
    DownDispersion = readr::col_double(),
    Metadata_permuted_p_value = readr::col_double(),
    min_permuted_value = readr::col_double(),
    max_permuted_value = readr::col_double()
)

score_df <- readr::read_tsv(score_file, col_types = score_cols)

print(dim(score_df))
head(score_df, 3)

for (dataset in datasets) {
    # Determine features that were selected for the dataset
    selected_features <- features_df %>% dplyr::filter(dataset == !!dataset) %>% dplyr::pull(features)

    # Subset the signature results for later merge
    subset_score_df <- score_df %>%
        dplyr::filter(signature == !!dataset) %>%
        dplyr::select(Metadata_unique_sample_name, TotalScore)

    # Subset the data to only these features
    subset_data_df <- dataset_df %>%
        dplyr::filter(Metadata_dataset == !!dataset) %>%
        dplyr::select(starts_with("Metadata_"), !!!selected_features) %>%
        dplyr::filter(Metadata_treatment == "0.1% DMSO")  %>%
        dplyr::left_join(subset_score_df, by = "Metadata_unique_sample_name")

    # Obtain correlation matrix
    correlation_matrix_df <- t(subset_data_df %>% dplyr::select(!!!selected_features)) %>% cor()

    # Create heatmap
    ht = Heatmap(
        correlation_matrix_df,
        name = "Correlation",

        top_annotation = HeatmapAnnotation(
            Clone = subset_data_df$Metadata_clone_number,
            Sensitivity = subset_data_df$Metadata_clone_type,
            Plate = subset_data_df$Metadata_Plate,
            Batch = subset_data_df$Metadata_batch,
            CellCount = anno_barplot(
                subset_data_df$Metadata_cell_count,
                height = unit(anno_name_height * 1.5, "cm")
            ),
            SignatureScore = subset_data_df$TotalScore,

            annotation_legend_param = list(
                Clone = list(
                    title_gp = gpar(fontsize = lgd_title_fontsize),
                    labels_gp = gpar(fontsize = lgd_label_fontsize)
                ),
                Sensitivity = list(
                    title_gp = gpar(fontsize = lgd_title_fontsize),
                    labels_gp = gpar(fontsize = lgd_label_fontsize)
                ),
                Plate = list(
                    title_gp = gpar(fontsize = lgd_title_fontsize),
                    labels_gp = gpar(fontsize = lgd_label_fontsize)
                ),
                Batch = list(
                    title_gp = gpar(fontsize = lgd_title_fontsize),
                    labels_gp = gpar(fontsize = lgd_label_fontsize)
                ),
                SignatureScore = list(
                    title_gp = gpar(fontsize = lgd_title_fontsize),
                    labels_gp = gpar(fontsize = lgd_label_fontsize),
                    color_bar = "continuous"
                )
            ),
            simple_anno_size = unit(anno_name_height, "cm"),
            annotation_name_gp = gpar(fontsize = lgd_label_fontsize)
        ),
        
        heatmap_legend_param = list(
            title = "Pearson\ncorrelation",
            color_bar = "continuous",
            col_fun = legend_scale_cols,
            title_gp = gpar(fontsize = lgd_title_fontsize),
            title_position = "topleft",
            labels_gp = gpar(fontsize = lgd_label_fontsize),
            at = c(-1, 0, 1)
        )
    )

    # Save heatmap to file
    fig_file <- file.path(output_dir, paste0("heatmap_", dataset, ".pdf"))
    pdf(fig_file)
    draw(ht)
    dev.off()
    
    draw(ht)
}
