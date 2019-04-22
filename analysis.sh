#!/bin/bash
#
# Gregory Way 2019
# Cellular Morphology Resistance Mechanisms
#
# Analysis Pipeline

# Step 0 - Convert all notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# Step 1 - Merge two batches of data together into one dataset
Rscript --vanilla scripts/nbconverted/merge-batch-data.r
