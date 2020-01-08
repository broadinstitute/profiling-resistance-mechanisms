# Inspect plate layout to identify potential plate artifacts

library(dplyr)
library(platetools)
library(ggplot2)
library(optparse)
library(yaml)

source(file.path("scripts", "plate_utils.R"))

option_list <- list(
  optparse::make_option(c("-p", "--directory"), help = "the location of the project backend directory"),
  optparse::make_option(c("-c", "--config"), help = "the location of the config yaml")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
sql_directory <- opt$directory
config <- opt$config
metadata_directory <- file.path("..", "1.process-profiles", "metadata")

yaml_list <- load_config_yaml(yaml_file=config)

#sql_directory <- "/home/ubuntu/bucket/projects/2018_05_30_ResistanceMechanisms_Kapoor/workspace/backend"

for (batch_idx in seq(1, length(yaml_list))) {
    batch_config <- yaml_list[[batch_idx]]

    batch <- batch_config$batch
    plates <- batch_config$plates
    audit_cols <- batch_config$auditcols

    meta_batch_dir <- file.path(metadata_directory, batch)
    platemap_barcode_file <- file.path(meta_batch_dir, "barcode_platemap.csv")
    platemap_barcode_df <- readr::read_csv(platemap_barcode_file, col_types = readr::cols())

    for (plate in plates) {

        cell_count_df <- count_cells(
            plate_name = plate,
            batch_id = batch,
            project_directory = sql_directory
       )

        plate_name <- platemap_barcode_df %>%
            dplyr::filter(Assay_Plate_Barcode == !!plate) %>%
            dplyr::pull(Plate_Map_Name)

        platemap_file <- file.path(meta_batch_dir, "platemap", paste0(plate_name, ".txt"))
        platemap_df <- readr::read_tsv(platemap_file, col_types = readr::cols())
        colnames(platemap_df) <- paste0("Metadata_", colnames(platemap_df))

        platemap_info_df <- cell_count_df %>%
            dplyr::right_join(platemap_df, by = c("Metadata_Well" = "Metadata_well_position"))

        if (length(audit_cols) == 2) {
            replicate_info <- paste(
                platemap_info_df %>% dplyr::pull(audit_cols[1]),
                platemap_info_df %>% dplyr::pull(audit_cols[2])
            )
        } else {
            replicate_info <- platemap_df %>% dplyr::pull(audit_cols[1])
        }

        platemap_info_df <- platemap_info_df %>%
            dplyr::mutate(Metadata_plate_replicate = replicate_info)

        # Visualize platemaps
        visualize_platemaps(platemap_info_df = platemap_info_df, batch = batch, plate = plate)
    }
}
