# Functions to facilitate data processing
#
# Usage:
# Import Only

require(dplyr)

load_data <- function(data_file, col_types = "infer") {
  # Load profiling data
  #
  # Arguments:
  # data_file - the relative or absolute file path to load profiling data
  # col_types - a readr::cols() object specifying the type of each column
  #             [default: "infer"] - if infer, use readr internals to infer column type
  #
  # Output:
  # The data stored in the file location
  
  if (col_types == "infer") {
    col_types = readr::cols()
  }

  return(readr::read_csv(data_file, col_types = col_types))
}


merge_data <- function(left_df, right_df, output_file = "none", output_gct = FALSE) {
  # Merge two profiling datasets based on common columns. Currently only supports two
  # dataframes.
  #
  # Arguments:
  # left_df - cell painting profiles and metadata
  # right_df - cell painting profiles and metadata
  # output_file - file name to output merged data
  # output_gct - boolean if gct file should also be generated and saved
  #
  # Output:
  # A single merged dataframe with common columns

  common_cols <- intersect(colnames(left_df), colnames(right_df))
  
  left_df <- left_df %>%
    dplyr::select(common_cols)
  
  right_df <- right_df %>%
    dplyr::select(common_cols)
  
  full_df <- dplyr::bind_rows(left_df, right_df)
  
  if (output_file != "none") {
    readr::write_csv(full_df, output_file)
  }
  
  if (output_gct) {
    require(curl)
    
    con <- curl("https://raw.githubusercontent.com/broadinstitute/cytominer_scripts/master/write_gct.R")
    source(con)
    close(con)

    gct_file <- paste0(tools::file_path_sans_ext(output_file), ".gct")
    write_gct(full_df, gct_file)
  }
  
  return(full_df)
}
