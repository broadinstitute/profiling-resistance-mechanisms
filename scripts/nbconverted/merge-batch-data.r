suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(curl))

con <- curl("https://raw.githubusercontent.com/broadinstitute/cytominer_scripts/master/write_gct.R")
source(con)
close(con)

# Set column types for reading in data
batch_cols = readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_Assay_Plate_Barcode = readr::col_character(),
    Metadata_Plate_Map_Name = readr::col_character(),
    Metadata_Batch_Number = readr::col_integer(),
    Metadata_well_position = readr::col_character(),
    Metadata_CellLine = readr::col_character()
)

file <- file.path("data", "2019_02_15_Batch1_20X", "HCT116bortezomib_normalized_variable_selected.csv")
batch1_df <- readr::read_csv(file, col_types = batch_cols)

print(dim(batch1_df))
head(batch1_df)

file <- file.path("data", "2019_03_20_Batch2", "207106_exposure320_normalized_variable_selected.csv")
batch2_df <- readr::read_csv(file, col_types = batch_cols)

print(dim(batch2_df))
head(batch2_df)

# What are the common columns
common_cols <- intersect(colnames(batch1_df), colnames(batch2_df))
length(common_cols)

# Combine dataframes together by common columns
batch1_commoncols_df <- batch1_df %>%
    dplyr::select(common_cols)

batch2_commoncols_df <- batch2_df %>%
    dplyr::select(common_cols)

all_commoncols_df <- dplyr::bind_rows(batch1_commoncols_df,
                                      batch2_commoncols_df)

dim(all_commoncols_df)

file <- file.path("data", "merged_intersected_variable_selected.csv")
readr::write_csv(all_commoncols_df, file)

# Output GCT file with combined columns
file <- paste0(tools::file_path_sans_ext(file), ".gct")
write_gct(all_commoncols_df, file)
