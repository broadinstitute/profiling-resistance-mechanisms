suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))

util_file <- file.path("..", "scripts", "processing_utils.R")
source(util_file)

backend_dir <- file.path("..", "..", "..", "backend")

batch <- "2019_02_15_Batch1_20X"
plate <- "HCT116bortezomib"

file <- file.path(backend_dir, batch, plate, paste0(plate, "_normalized_variable_selected.csv"))
batch1_df <- load_data(data_file = file)

print(dim(batch1_df))
head(batch1_df, 2)

batch <- "2019_03_20_Batch2"
plate <- "207106_exposure320"

file <- file.path(backend_dir, batch, plate, paste0(plate, "_normalized_variable_selected.csv"))
batch2_df <- load_data(data_file = file)

print(dim(batch2_df))
head(batch2_df, 2)

merge_file <- file.path("data", "merged_intersected_variable_selected.csv")
full_df <- merge_data(batch1_df, batch2_df, output_file = merge_file, output_gct = FALSE)

dim(full_df)

batch <- "2019_06_25_Batch3"
plate <- "MutClones"

file <- file.path(backend_dir, batch, plate, paste0(plate, "_normalized_variable_selected.csv"))
batch3_mut_df <- load_data(file)

print(dim(batch3_mut_df))
head(batch3_mut_df, 2)

plate <- "WTClones"

file <- file.path(backend_dir, batch, plate, paste0(plate, "_normalized_variable_selected.csv"))
batch3_wt_df <- load_data(file)

print(dim(batch3_wt_df))
head(batch3_wt_df, 2)

merge_file <- file.path("data", paste0(batch, "_merged_intersected_variable_selected.csv"))
full_df <- merge_data(batch3_mut_df, batch3_wt_df, output_file = merge_file, output_gct = FALSE)

dim(full_df)
