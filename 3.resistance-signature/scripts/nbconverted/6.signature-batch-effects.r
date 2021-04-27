library(dplyr)
library(platetools)
library(ggplot2)
library(optparse)
library(yaml)

audit_dir <- file.path("..", "1.profiling-audit")
source(file.path(audit_dir, "scripts", "plate_utils.R"))

plot_total_score_platemap <- function(sig_df, output_file, plate, batch, spectral_limits) {
    
    title_base <- paste0("Batch: ", batch, "\nPlate: ", plate)
    
    platemap_sig_gg <- (
        platetools::raw_map(
            data = platemap_signature_df$TotalScore,
            well = platemap_signature_df$Metadata_Well,
            plate = 96,
            size = 4
        )
        + ggtitle(title_base)
        + scale_fill_distiller(name = "Bortezomib\nsignature", palette = "Spectral", limits = spectral_limits)
        + scale_shape_discrete(name = "Clone type")
        + geom_point(
            aes(shape = platemap_signature_df$Metadata_clone_type),
            size = 1,
            position = position_nudge(x = 0.1, y = -0.1)
        )
        + platemap_theme
        + theme(
            legend.title = element_text(size = 5),
            legend.position = "right",
            legend.key.size = unit(0.5, "line"),
            legend.margin = margin(unit="mm"),
            legend.title.align = 0
        )
        + guides(
          shape = guide_legend(order = 1, override.aes = list(size = 1))
        )
    )
    
    ggsave(
      plot = platemap_sig_gg,
      filename = output_file,
      height = 1.75,
      width = 3,
      dpi = 500
    )
}

# Load signature data
dataset <- "bortezomib"
sig_dir <- file.path("results", "singscore")
profile_dir <- "../0.generate-profiles"
metadata_dir <- file.path(profile_dir, "metadata")
output_dir <- file.path("figures", "plate_effects")

sig_results_file <- file.path(sig_dir, paste0("singscore_results", dataset, ".tsv.gz"))

bortezomib_signature_data <- list(
    "2021_03_03_Batch12" = c("219907"),
    "2021_03_03_Batch13" = c("219973"),
    "2021_03_03_Batch14" = c("219901"),
    "2021_03_03_Batch15" = c("219956"),
    "2021_03_05_Batch16" = c("220039"),
    "2021_03_05_Batch17" = c("220040"),
    "2021_03_12_Batch18" = c("220055")
)

clone_label_recode <- c(
  "Clone A" = "CloneA",
  "Clone E" = "CloneE",
  "WT parental" = "WT_parental"
)

sig_cols <- readr::cols(
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

signature_df <- readr::read_tsv(sig_results_file, col_types = sig_cols) %>%
    dplyr::filter(Metadata_model_split != "inference")

spectral_limits <- c(-1, 1) * max(abs(signature_df$TotalScore))

print(dim(signature_df))
head(signature_df, 4)

for (batch in names(bortezomib_signature_data)) {
    meta_batch_dir <- file.path(metadata_dir, batch)
    platemap_barcode_file <- file.path(meta_batch_dir, "barcode_platemap.csv")

    platemap_barcode_df <- readr::read_csv(
      platemap_barcode_file, col_types = readr::cols()
    )
    
    for (plate in bortezomib_signature_data[[batch]]) {
        plate_name <- platemap_barcode_df %>%
            dplyr::filter(Assay_Plate_Barcode == !!plate) %>%
            dplyr::pull(Plate_Map_Name)

        platemap_file <- file.path(
          meta_batch_dir, "platemap", paste0(plate_name, ".txt")
        )
        platemap_df <- readr::read_tsv(
            platemap_file, col_types = readr::cols(.default = readr::col_character())
        )
        colnames(platemap_df) <- paste0("Metadata_", colnames(platemap_df))
        
        platemap_df$Metadata_clone_number <-
            dplyr::recode(platemap_df$Metadata_clone_number, !!!clone_label_recode)
        
        platemap_signature_df <- signature_df %>%
            dplyr::select(
                Metadata_Plate,
                Metadata_Well,
                Metadata_clone_number,
                Metadata_clone_type,
                Metadata_treatment,
                TotalScore
            ) %>%
            dplyr::right_join(
                platemap_df,
                by = c(
                    "Metadata_Plate" = "Metadata_plate_map_name",
                    "Metadata_Well" = "Metadata_well_position",
                    "Metadata_clone_number" = "Metadata_clone_number",
                    "Metadata_treatment" = "Metadata_treatment"
                )
            )
        
        output_file <- file.path(output_dir, paste0(batch, "_", plate, "_signature_plate_effect.png"))
        
        plot_total_score_platemap(
            sig_df = platemap_signature_df,
            output_file = output_file,
            batch = batch,
            plate = plate,
            spectral_limits = spectral_limits
        )
    }
}
