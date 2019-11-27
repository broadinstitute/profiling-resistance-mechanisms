# Morpheus Figures

[Morpheus](https://software.broadinstitute.org/morpheus/) was developed at the Broad Institute and is a web based application.
The application can be used to visualize large datasets.
Here, we use Morpheus to visualize hierarchically clustered similarity heatmaps of our cell painting data.

## Datasets Used

We input the following `.gct` files into Morpheus:

* Merged Batch - `data/merged_intersected_variable_selected.gct`
* Batch 1 - `data/2019_02_15_Batch1_20X/HCT116bortezomib_normalized_variable_selected.gct`
* Batch 2 - `data/2019_03_20_Batch2/2017106_exposure320_normalized_variable_selected.gct`

## Morpheus Parameters

We performed the following procedure to generate the heatmaps.

1. Compute similarity matrix with `Tools` > `Similarity Matrix` with options `"Metric"="Pearson correlation"` and `"Compute matrix for"="Columns"`
2. Configure `Options` > `Annotations` > `Row annotations` to display `"CellLine"`, `"Dosage"`, and `"Batch Number"`
3. Configure `Options` > `Annotations` > `Column annotations` to display `"CellLine"`, `"Dosage"`, and `"Batch Number"`
4. Perform `Hierarchical Clustering` with `Tools` > `Hierarchical Clustering` with options `"Metric": "Matrix values (for a precomputed similarity matrix)"`, `"Linkage method": "Average"`, `"Cluster": "Rows and Columns"`
5. Configured the `Batch Number` to display colors and text; `Dosage` to display a bar chart and text; `CellLine` to display color and text.

## Output

We saved the heatmaps as `.png` and `.pdf` files (shown below)

### Merged Batch

![merged batch](https://raw.githubusercontent.com/broadinstitute/2018_05_30_ResistanceMechanisms_Kapoor/master/figures/morpheus/merged_batch_morpheus_heatmap.png)

### Batch 1

![batch 1](https://raw.githubusercontent.com/broadinstitute/2018_05_30_ResistanceMechanisms_Kapoor/master/figures/morpheus/batch1_morpheus_heatmap.png)

### Batch 2

![batch 2](https://raw.githubusercontent.com/broadinstitute/2018_05_30_ResistanceMechanisms_Kapoor/master/figures/morpheus/batch2_morpheus_heatmap.png)
