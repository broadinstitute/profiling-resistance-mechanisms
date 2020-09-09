#!/bin/bash
#
# Gregory Way 2020
# Identifying morphology markers of drug resistance
# 2.describe-data/run_pipeline.sh
#
# Perform the full pipeline describing data

# Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 0 - Generate summary figures about the existing data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.describe-data.ipynb

# Step 1 - Merge data and output gct files
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.merge-datasets-gct.ipynb

# Step 2 - Generate UMAP figures
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.umap-aggregate-profiles.ipynb

# Step 3 - Generate UMAP figures
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.count-single-cells.ipynb
