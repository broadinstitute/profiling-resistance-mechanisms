# Utility functions to facilitate platemap layout analysis and visualization

platemap_theme <- ggplot2::theme(
    title = element_text(size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "line")
)


count_cells <- function(plate_name, batch_id, project_directory) {

    require(RSQLite)
    require(dplyr)

    # Determine SQL file name
    sqlite_file <- file.path(
      project_directory,
      batch_id,
      plate_name,
      paste0(plate_name, ".sqlite")
    )

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
    con <- file(yaml_file, "r")
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


visualize_platemaps <- function(platemap_info_df, batch, plate, output_dir) {

    require(platetools)
    require(ggplot2)

    title_base <- paste0(batch, "\nPlate: ", plate)

    plate_replicate_gg <-
        platetools::raw_map(
            data = platemap_info_df$Metadata_plate_replicate,
            well = platemap_info_df$Metadata_Well,
            plate = 96,
            size = 4
        ) +
        ggplot2::ggtitle(paste0(title_base, "\nReplicate Info")) +
        ggplot2::theme_dark() +
        ggplot2::scale_fill_discrete() +
        platemap_theme +
        ggplot2::guides(fill = guide_legend(ncol = 5))

    replicate_file <- file.path(
      output_dir,
      paste0(batch, plate, "_plate_effects_replicates.png")
    )
    ggplot2::ggsave(
      plot = plate_replicate_gg,
      filename = replicate_file,
      height = 4,
      width = 3
    )

    cell_count_gg <-
        platetools::raw_map(
            data = platemap_info_df$cell_count,
            well = platemap_info_df$Metadata_Well,
            plate = 96,
            size = 4
        ) +
        ggplot2::ggtitle(paste0(title_base, "\nCell Count")) +
        ggplot2::theme_dark() +
        ggplot2::scale_fill_continuous(name = "Cells") +
        platemap_theme +
        ggplot2::guides(fill = guide_legend(ncol = 5))

    count_file <- file.path(
      output_dir,
      paste0(batch, plate, "_plate_effects_cell_count.png")
    )
    ggplot2::ggsave(
      plot = cell_count_gg,
      filename = count_file,
      height = 4,
      width = 3
    )
}


visualize_site_counts <- function(platemap_info_df, batch, plate, output_dir) {
  require(platetools)
  require(ggplot2)

  plate_title <- paste0(batch, "\nPlate: ", plate)

  site_count_gg <- ggplot2::ggplot(platemap_info_df,
                                   aes(x = factor(Metadata_Site), y = cell_count)) +
      ggplot2::geom_bar(stat = "identity", aes_string(fill = audit_cols[1])) +
      ggplot2::facet_wrap("~Metadata_Well") +
      ggplot2::scale_fill_discrete(name = audit_cols[1]) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.text = element_text(size = 5, color = "black"),
                     strip.background = element_rect(colour = "black", fill = "#fdfff4")) +
      ggplot2::ylab("Cell Count") +
      ggplot2::xlab("Site") +
      ggplot2::ggtitle(plate_title) +
      ggplot2::theme(
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        title = element_text(size = 10)
      )

  site_count_file <- file.path(
        output_dir,
        paste0(batch, plate, "_plate_effects_cell_count_by_site.png")
      )

  ggplot2::ggsave(
    plot = site_count_gg,
    filename = site_count_file,
    height = 4.5,
    width = 6
  )
}
