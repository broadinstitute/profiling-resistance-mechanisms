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

opt <- optparse::parse_arge(optparse::OptionParser(option_list = option_list))
directory <- opt$directory
config <- opt$config

yaml_info <- load_config_yaml(yaml_file=config)


