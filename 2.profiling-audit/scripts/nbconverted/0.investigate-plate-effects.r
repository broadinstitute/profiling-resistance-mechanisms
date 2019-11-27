suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(platetools))
suppressPackageStartupMessages(library(ggplot2))

count_cells <- function(plate_name) {
    # Determine SQL file name
    project <- "2018_05_30_ResistanceMechanisms_Kapoor"
    bucket_dir <- file.path("~", "bucket", "projects", project, "workspace", "backend")
    sqlite_file <- file.path(bucket_dir, batch, plate_name, paste0(plate_name, ".sqlite"))
    
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
        dplyr::count()
    
    return(count_df)
}

batch <- "2019_06_25_Batch3"
backend_dir <- file.path("..", "..", "..", "backend", batch)

backend_folders = list.files(backend_dir, full.names = FALSE)

backend_files = list()
for (backend_file in backend_folders) {
    backend_files[[backend_file]] <-
        file.path(backend_dir, backend_file, paste0(backend_file, "_normalized_variable_selected.csv"))
}

backend_files

mut_df <- readr::read_csv(backend_files[["MutClones"]],
                          col_types=readr::cols())

print(dim(mut_df))
head(mut_df, 2)

wt_df <- readr::read_csv(backend_files[["WTClones"]],
                         col_types=readr::cols())

print(dim(wt_df))
head(wt_df, 2)

cor_file <- file.path("results", "WT_median_replicate_correlation.tsv")
wt_median_cor_df <- readr::read_tsv(cor_file, col_types = readr::cols())

cor_file <- file.path("results", "MUT_median_replicate_correlation.tsv")
mut_median_cor_df  <- readr::read_tsv(cor_file, col_types = readr::cols())

wt_median_cor_df <- wt_df %>%
    dplyr::select(Metadata_Well, Metadata_Plate_Map_Name, Metadata_clone_number) %>%
    dplyr::inner_join(wt_median_cor_df, by = "Metadata_clone_number")

mut_median_cor_df <- mut_df %>%
    dplyr::select(Metadata_Well, Metadata_Plate_Map_Name, Metadata_clone_number) %>%
    dplyr::inner_join(mut_median_cor_df, by = "Metadata_clone_number")

# Mutant Clone Plate
mut_count_file <- file.path("results", "MutClones_cell_count_by_well.tsv")
if (!file.exists(mut_count_file)) {
    mut_count <- count_cells("MutClones") %>%
        dplyr::collect()
    mut_count %>% readr::write_tsv(mut_count_file)
} else {
    mut_count <- readr::read_tsv(mut_count_file, col_types=readr::cols())
}

head(mut_count)

# Wildtype Clone Plate
wt_count_file <- file.path("results", "WTClones_cell_count_by_well.tsv")
if (!file.exists(wt_count_file)) {
    wt_count <- count_cells("WTClones") %>%
        dplyr::collect()
    wt_count %>% readr::write_tsv(wt_count_file)
} else {
    wt_count <- readr::read_tsv(wt_count_file, col_types=readr::cols())
}

head(wt_count)

plate_mut_gg <- 
    platetools::raw_map(data = mut_df$Metadata_clone_number,
                        well = mut_df$Metadata_Well,
                        plate = 96) +
    ggtitle("Mutant Clones - Replicates") +
    theme_dark() +
    scale_fill_discrete() +
    theme(legend.position = "none")

fig_file <- file.path("figures", "plate_effects", paste0(batch, "_mut_replicates.png"))
ggsave(fig_file, height = 4, width = 6)

plate_mut_gg

plate_wt_gg <- 
    platetools::raw_map(data = wt_df$Metadata_clone_number,
                        well = wt_df$Metadata_Well,
                        plate = 96) +
    ggtitle("Wild-type Clones - Replicates") +
    theme_dark() +
    scale_fill_discrete() +
    theme(legend.position = "none")

fig_file <- file.path("figures", "plate_effects", paste0(batch, "_wt_replicates.png"))
ggsave(fig_file, height = 4, width = 6)

plate_wt_gg

cell_count_mut_gg <- 
    platetools::raw_map(data = mut_count$n,
                        well = mut_count$Metadata_Well,
                        plate = 96) +
    ggtitle("Mutant Clones - Cell Count") +
    theme_dark() +
    scale_fill_continuous(name = "Cells")

fig_file <- file.path("figures", "plate_effects", paste0(batch, "_mut_cell_counts.png"))
ggsave(fig_file, height = 4, width = 6)

cell_count_mut_gg

cell_count_wt_gg <- 
    platetools::raw_map(data = wt_count$n,
                        well = wt_count$Metadata_Well,
                        plate = 96) +
    ggtitle("Wild-type Clones - Cell Count") +
    theme_dark() +
    scale_fill_continuous(name = "Cells")

fig_file <- file.path("figures", "plate_effects", paste0(batch, "_wt_cell_counts.png"))
ggsave(fig_file, height = 4, width = 6)

cell_count_wt_gg

mut_replicate <- mut_median_cor_df %>% dplyr::filter(replicate == "TRUE")

cor_mut_gg <- 
    platetools::raw_map(data = mut_replicate$correlation,
                        well = mut_replicate$Metadata_Well,
                        plate = 96) +
    ggtitle("Mutant Clones - Replicate Correlation") +
    theme_dark() +
    scale_fill_continuous(type = "viridis", name = "Pearson\nCorrelation")

fig_file <- file.path("figures", "plate_effects", paste0(batch, "_mut_replicate_correlation.png"))
ggsave(fig_file, height = 4, width = 6)

cor_mut_gg

wt_replicate <- wt_median_cor_df %>% dplyr::filter(replicate == "TRUE")

cor_wt_gg <- 
    platetools::raw_map(data = wt_replicate$correlation,
                        well = wt_replicate$Metadata_Well,
                        plate = 96) +
    ggtitle("Wild-type Clones - Replicate Correlation") +
    theme_dark() +
    scale_fill_continuous(type = "viridis", name = "Pearson\nCorrelation")

fig_file <- file.path("figures", "plate_effects", paste0(batch, "_wt_replicate_correlation.png"))
ggsave(fig_file, height = 4, width = 6)

cor_wt_gg
