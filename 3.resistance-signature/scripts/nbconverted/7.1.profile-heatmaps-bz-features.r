suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(devtools))
source(file.path("utils/viz.R"))

# The correct version of ComplexHeatmap was not available on anaconda
# https://anaconda.org/bioconda/bioconductor-complexheatmap
# Install it here
if (!("ComplexHeatmap" %in% rownames(installed.packages()))) {
    install_github("jokergoo/ComplexHeatmap")
}

suppressPackageStartupMessages(library(ComplexHeatmap))

# Set constants
dataset <- "bortezomib"
data_dir <- "data"

output_dir <- file.path("figures", "heatmaps")

lgd_title_fontsize = 9
lgd_label_fontsize = 6.5
anno_name_height = 0.45

# Load profiles
dataset_file <- file.path(data_dir, paste0(dataset, "_signature_analytical_full_set.csv"))

data_cols <- readr::cols(
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

dataset_df <- readr::read_csv(dataset_file, col_types=data_cols)

print(dim(dataset_df))
head(dataset_df, 3)

# Load signatures
sig_dir <- file.path("results", "signatures")
signature_file <- file.path(sig_dir, paste0("signature_summary_", dataset, "_signature.tsv.gz"))

sig_col_types <- readr::cols(
    features = readr::col_character(),
    non_specific_exclude = readr::col_logical(),
    final_signature = readr::col_logical(),
    dataset = readr::col_character()
)

signature_df <- readr::read_tsv(signature_file, col_types = sig_col_types)

print(dim(signature_df))
head(signature_df, 4)

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
score_dir <- file.path("results", "singscore")
score_file <- file.path(score_dir, paste0("singscore_results_", dataset, "_full_set.csv"))

score_cols <- readr::cols(
    .default = readr::col_character(),
    Metadata_clone_type_indicator = readr::col_integer(),
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

score_df <- readr::read_csv(score_file, col_types = score_cols)
print(dim(score_df))
head(score_df, 3)

# Set colors
legend_scale_cols = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

legend_colors = c(
        "otherclone"="#66C2A5",
        "test"="#FC8D62",
        "training"="#8DA0CB",
        "validation"="#E78AC3",
        "holdout"="#A6D854"
        )


plate_col = c(
    "219901" = "#E1DAAE",
    "219907" = "#FF934F",
    "219956" = "#CC2D35",
    "219973" = "#058ED9",
    "220039" = "#848FA2",
    "220040" = "#2D3142",
    "220055" = "#FFC857",
    "221093" = "#E495A5",
    "221094" = "#ABB065",
    "221058" = "#39BEB1",
    "221057" = "#ACA4E2"
    
)

clonetype_col = c("Clone" = "#785EF0", "Parental" = "#DC267F")

# From colorbrewer2.org
spectral_palette <- c(
    '#9e0142',
    '#d53e4f',
    '#f46d43',
    '#fdae61',
    '#fee08b',
    '#ffffbf',
    '#e6f598',
    '#abdda4',
    '#66c2a5',
    '#3288bd',
    '#5e4fa2'
)

spectral_limits <- c(-1, 1) * max(abs(score_df$TotalScore))
spectral_breaks <- seq(spectral_limits[1] * 100, spectral_limits[2] * 100, length(spectral_palette) + 1) / 100

signature_col = circlize::colorRamp2(rev(spectral_breaks), spectral_palette)
sensitivity_col = c("resistant" = "#332a2a", "sensitive" = "#c3c7c9")

#for (feature_select_type in c("selected", "signature")) {
    
    feature_select_type <- "signature"
    # Determine features that were selected for the dataset
    if (feature_select_type == "signature") {
        # Signature only features
        selected_features <- signature_df %>%
            dplyr::filter(final_signature) %>%
            dplyr::pull(features)
    } else {
        # Features that survived the standard feature selection steps
        selected_features <- features_df %>%
            dplyr::filter(dataset == !!dataset) %>%
            dplyr::pull(features)
    }

    # Subset the signature results for later merge
    subset_score_df <- score_df %>%
        dplyr::filter(dataset == !!dataset) %>%
        dplyr::select(Metadata_unique_sample_name, TotalScore)

    # Subset the data to only these features
    subset_data_df <- dataset_df %>%
        dplyr::filter(Metadata_dataset == !!dataset) %>%
        dplyr::select(starts_with("Metadata_"), !!!selected_features) %>%
        dplyr::filter(Metadata_treatment == "0.1% DMSO")  %>%
        dplyr::left_join(subset_score_df, by = "Metadata_unique_sample_name") %>%
        dplyr::mutate(Metadata_parental = "Clone")


total_score<-subset_data_df$TotalScore
subset_data_df = subset(subset_data_df, select = -c(TotalScore) )
subset_data_df[order(as.character(subset_data_df$Metadata_clone_type)),]

ht=
Heatmap(
   subset_data_df%>% select(-matches("Metadata")) %>% select(1:45) %>% as.matrix() %>% t(),
       top_annotation = HeatmapAnnotation(
       resistance_type = subset_data_df %>% pull(Metadata_clone_type),
       model_split= subset_data_df %>% pull(Metadata_model_split),
       total_score=total_score
      ),
      
       heatmap_legend_param = list(
            title = "Normalized\nfeature\nvalue",
            color_bar = "continuous",
            title_gp = gpar(fontsize = lgd_title_fontsize),
            title_position = "topleft",
            col_fun = legend_scale_cols,
            labels_gp = gpar(fontsize = lgd_label_fontsize)
            
      ),
    row_names_gp = gpar(fontsize = 5),
    clustering_method_columns = "average",
    clustering_method_rows = "average",
    column_split = subset_data_df %>% pull(Metadata_clone_type)
    )
    

draw(ht,annotation_legend_side = "left", heatmap_legend_side = "left")


