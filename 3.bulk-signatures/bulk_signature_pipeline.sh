#!/bin/bash
#
# Gregory Way 2020
# The full analysis pipeline to reproduce the bulk signature analysis in which
# we identify a morphology signature that differentiates sensitive and resistant
# cell line clones before treatment.
#
# Note that all results, figures, and datasets are already included in the repository.
#
# 3.bulk-signatures/bulk_signature_pipeline.sh

# Exit after error
set -e

# Setup conda commands
eval "$(conda shell.bash hook)"

# Make sure the correct environment is activated; we will switch in this script
conda activate resistance-mechanisms

# Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 0 - Compile profiles across batches into single datasets
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.compile-bulk-datasets.ipynb

# Step 1 - Identify features that are most descriptive of resistant vs. sensitive clones
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.derive-bulk-signatures.ipynb

# Switch to conda environment for singscore signature application
conda deactivate
conda activate singscore

# Step 2 - Apply signature features to identify enrichment
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.apply-bulk-signatures.ipynb

# Switch back to original conda environment
conda deactivate
conda activate resistance-mechanisms

# Step 3 - Assess signature performance
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.performance-metrics.ipynb

# Step 4 - Visualize performance metrics
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 4.visualize-performance.ipynb

# Step 5 - Visualize signature features
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 5.visualize-signature-features.ipynb

# Step 6 - Investigate if resistance score is associated with cell count
# Note: we removed features that were signifantly associated with cell count in building
# the signature - it is still possible though that cell count plays a role!
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 6.score-by-cell-count.ipynb

# Step 7 - Compile batch 3 data into a separate dataset for additional testing
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 7.combine-batch3.ipynb

# Switch to conda environment for another singscore signature application
conda deactivate
conda activate singscore

# Step 8 - Apply the cloneAE (bortezomib) signature to batch 3 data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 8.apply-cloneAE-signature-batch3.ipynb

# Switch back to original conda environment
conda deactivate
conda activate resistance-mechanisms

# Step 9 - Compile profile heatmaps
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 9.profile-heatmaps.ipynb
