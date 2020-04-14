# Inspect plate layout to identify potential plate artifacts

library(dplyr)
library(platetools)
library(ggplot2)
library(optparse)
library(yaml)

source(file.path("scripts", "plate_utils.R"))

option_list <- list(
  optparse::make_option(
    c("-p", "--profile_dir"), help = "the location of the profile directory"
  ),
  optparse::make_option(
    c("-c", "--config"), help = "the location of the config yaml")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

config <- opt$config
profile_dir <- opt$profile_dir
metadata_directory <- file.path(profile_dir, "metadata")
cell_count_directory <- file.path(profile_dir, "cell_counts")

yaml_list <- load_config_yaml(yaml_file=config)

for (batch_idx in seq(1, length(yaml_list))) {

    batch_config <- yaml_list[[batch_idx]]

    batch <- batch_config$batch
    plates <- batch_config$plates
    audit_cols <- batch_config$auditcols

    meta_batch_dir <- file.path(metadata_directory, batch)
    platemap_barcode_file <- file.path(meta_batch_dir, "barcode_platemap.csv")

    platemap_barcode_df <- readr::read_csv(
      platemap_barcode_file, col_types = readr::cols()
    )

    for (plate in plates) {
        plate <- as.character(plate)

        output_figure_dir <- file.path("figures", batch, plate)

        cell_count_file <- file.path(
          cell_count_directory,
          paste(batch, plate, "cell_count.tsv", sep="_")
        )
        cell_count_df <- readr::read_tsv(cell_count_file, col_types = readr::cols())

        plate_name <- platemap_barcode_df %>%
            dplyr::filter(Assay_Plate_Barcode == !!plate) %>%
            dplyr::pull(Plate_Map_Name)

        platemap_file <- file.path(
          meta_batch_dir, "platemap", paste0(plate_name, ".txt")
        )
        platemap_df <- readr::read_tsv(platemap_file, col_types = readr::cols())
        colnames(platemap_df) <- paste0("Metadata_", colnames(platemap_df))

        # Confirm that the wells are not duplicated (a common mistake)
        duplicated_well <- sum(duplicated(platemap_df$Metadata_well_position))
        if (duplicated_well > 0) {
            print(
              paste(
                "WARNING! Duplicate well information detected in: Batch:",
                batch,
                "Plate:",
                plate
              )
            )
        }

       if ("Metadata_Site" %in% colnames(cell_count_df)) {
           cell_merge_count_df <- cell_count_df %>%
               dplyr::select(Metadata_Well, Metadata_Site, cell_count)
       } else {
           cell_merge_count_df <- cell_count_df %>%
               dplyr::select(Metadata_Well, !!plate) %>%
               dplyr::rename(cell_count = !!plate)
       }

        platemap_info_df <- cell_merge_count_df %>%
            dplyr::right_join(
              platemap_df, by = c("Metadata_Well" = "Metadata_well_position")
            )

        if (length(audit_cols) == 2) {
            replicate_info <- paste(
                platemap_info_df %>% dplyr::pull(audit_cols[1]),
                platemap_info_df %>% dplyr::pull(audit_cols[2])
            )
        } else {
            replicate_info <- platemap_info_df %>% dplyr::pull(audit_cols[1])
        }

        platemap_info_df <- platemap_info_df %>%
            dplyr::mutate(Metadata_plate_replicate = replicate_info)

        if ("Metadata_Site" %in% colnames(cell_count_df)) {

          visualize_site_counts(
            platemap_info_df = platemap_info_df,
            batch = batch,
            plate = plate,
            output_dir = output_figure_dir
          )

          platemap_info_df <- platemap_info_df %>%
            dplyr::group_by(Metadata_Well) %>%
            dplyr::mutate(cell_count = sum(cell_count)) %>%
            dplyr::select(!Metadata_Site) %>%
            dplyr::distinct()
        }

        visualize_platemaps(
            platemap_info_df = platemap_info_df,
            batch = batch,
            plate = plate,
            output_dir = output_figure_dir
        )

    }
}
