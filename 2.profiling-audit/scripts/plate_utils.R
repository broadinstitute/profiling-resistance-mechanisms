# Utility functions to facilitate platemap layout analysis and visualization

count_cells <- function(plate_name, batch_id, project_directory) {

    require(RSQLite)
    require(dplyr)

    # Determine SQL file name
    sqlite_file <- file.path(project_directory, batch_id, plate_name, paste0(plate_name, ".sqlite"))

    # Setup a connection
    con <- RSQLite::dbConnect(drv = RSQLite::SQLite(),
                              dbname = sqlite_file)

    # Connect to image table
    image_df <- dplyr::tbl(
        con,
        dbplyr::sql("SELECT TableNumber, ImageNumber, Metadata_Plate, Metadata_Well from image")
    )

    # Connect to cells table
    cells_df <- dplyr::tbl(
        con,
        dbplyr::sql("SELECT TableNumber, ImageNumber from cells")
    )

    # Merge them together and count cells per well
    count_df <- cells_df %>%
        dplyr::inner_join(image_df,
                          by = c("TableNumber", "ImageNumber")) %>%
        dplyr::group_by(Metadata_Plate, Metadata_Well) %>%
        dplyr::count() %>%
        dplyr::as_tibble()

    RSQLite::dbDisconnect(con)
    return(count_df)
}


load_config_yaml <- function(yaml_file) {
    # Load the file for the audit to determine plate information

    require(yaml)

    # Read in yaml config file
    con <- file(config_file, "r")
    yaml_result <- readLines(con, n = -1)
    unlink(con)

    # Process yaml compound file
    yaml_list <- list()
    yaml_idx <- 0
    for (a in yaml_result) {
        if (a == "---") {
            if (yaml_idx > 0) {
                yaml_list[[yaml_idx]] <- yaml.load(yaml_string)
            }
            yaml_string <- ""
            a <- ""
            yaml_idx <- yaml_idx + 1
        } else {
            yaml_string <- paste0(yaml_string, "\n", a)
        }
    }
    yaml_list[[yaml_idx]] <- yaml.load(yaml_string)

    return(yaml_list)
}


visualize_platemaps <- function(platemap_info_df, batch, plate) {

    require(platetools)
    require(ggplot2)

    title_base <- paste0(batch, ": ", plate)
    figure_directory <- file.path("figures", "plate_effects")

    plate_replicate_gg <-
        platetools::raw_map(data = platemap_info_df$Metadata_plate_replicate,
                            well = platemap_info_df$Metadata_Well,
                            plate = 96) +
          ggplot2::ggtitle(paste0(title_base, "\nReplicate Info")) +
          ggplot2::theme_dark() +
          ggplot2::scale_fill_discrete() +
          ggplot2::theme(legend.position = "none")

    replicate_file <- file.path(figure_directory, paste0(batch, plate, "_plate_effects_replicates.png"))
    ggplot2::ggsave(plot = plate_replicate_gg, filename = replicate_file, height = 4, width = 6)

    cell_count_gg <-
        platetools::raw_map(data = platemap_info_df$n,
                            well = platemap_info_df$Metadata_Well,
                            plate = 96) +
          ggplot2::ggtitle(paste0(title_base, "\nCell Count")) +
          ggplot2::theme_dark() +
          ggplot2::scale_fill_continuous(name = "Cells")

    count_file <- file.path(figure_directory, paste0(batch, plate, "_plate_effects_cell_count.png"))
    ggplot2::ggsave(plot = cell_count_gg, filename = count_file, height = 4, width = 6)
}
