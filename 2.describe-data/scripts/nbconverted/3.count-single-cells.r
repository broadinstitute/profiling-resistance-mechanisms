suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

cell_count_dir <- file.path("..", "0.generate-profiles", "cell_counts")
total_cell_count <- 0
for (cell_count_file in list.files(cell_count_dir, full.names = TRUE)) {
    cell_count_df <- readr::read_tsv(cell_count_file, col_types=readr::cols())
    
    if (
        any(
            grepl(cell_count_file,
                  c("2019_02_15_Batch1_20X", "2019_02_15_Batch1_40X", "2019_03_20_Batch2"),
                  fixed=TRUE)
        )
    ) {
        use_cols <- c("Metadata_CellLine", "Metadata_Dosage", "Metadata_Plate")
    } else if (
        any(
            grepl(cell_count_file,
                  c("2019_06_25_Batch3"),
                  fixed=TRUE)
        )
    ) {
        use_cols <- c("Metadata_clone_number", "Metadata_Plate")
    } else {
        use_cols <- c("Metadata_clone_number", "Metadata_treatment", "Metadata_Plate")
    }
    
    total_cell_count <- total_cell_count + sum(cell_count_df$cell_count)
}

cell_count_file <- file.path("results", "total_cell_count.txt")
names(total_cell_count) <- "total_cell_count"

write.table(total_cell_count, sep=",", cell_count_file, col.names=FALSE)

total_cell_count
