To generate profiles, two CellProfiler pipelines were run for each batch.

Images are first loaded into the `illumination_correction` pipeline, which is used to generate per-plate illumination correction functions.
This pipeline will export .npy files for use in the subsequent analysis pipeline.

Following this, the images and the illumination correction .npy files are loaded into the `analysis` pipeline. This pipeline will perform
the measurements and export them into .csv output files.

Note: The final module in each pipeline is `CreateBatchFiles`. If this is enabled, running the pipeline will generate a CellProfiler batch file 
which can be used to perform the analysis on a computational cluster.
If this module is disabled the analysis will be performed on the local machine.
