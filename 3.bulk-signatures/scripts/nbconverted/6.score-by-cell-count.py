#!/usr/bin/env python
# coding: utf-8

# ## Determine if resistance score is associated with cell count
# 
# Gregory Way, 2020
# 
# For every cell line clone and treatment combination, we know two things:
# 
# 1. Profile resistance score (via singscore)
# 2. Cell count before and after treatment (must assume uniform seeding density)
# 
# In this notebook, I perform a series of systematic visualizations to detect if there is any association between the resistance score and cell count after treatment.
# The resistance score was calculated in untreated clones only.
# An association between the resistance score and cell count after treatment would provide evidence that we are infact detecting a resistance signature.
# 
# ### Score interpretation
# 
# A high resistance score indicates that the profile has high concordance with samples annotated as resistant to treatment.
# A low score indicates that the profile has high concordance with samples annotated as sensitive.
# This score is bounded between -1 and 1. 
# Therefore, a sample with a high resistance score is hypothesized to be very resistant to the treatment, whereas a sample with a low resistance score is hypothesized to be very sensitive to the treatment.

# In[1]:


import pathlib
import numpy as np
import pandas as pd

import plotnine as gg


# In[2]:


# Set file names and directories
results_file = pathlib.Path("results/singscore/full_bulk_signature_singscore_results.tsv.gz")
cell_count_dir = pathlib.Path("../0.generate-profiles/cell_counts/")

output_fig_dir = pathlib.Path("figures/cell_count_associations")
output_fig_dir.mkdir(exist_ok=True)


# In[3]:


# Load singscore results
scores_df = pd.read_csv(results_file, sep="\t")

# This results file includes model predictions for all three signatures
# Only subset the targeted signature subset
scores_df = scores_df.loc[scores_df.dataset == scores_df.signature, :]
scores_df.Metadata_Plate = scores_df.Metadata_Plate.astype(str)

print(scores_df.shape)
scores_df.head()


# In[4]:


# This is the number of wells in each plate
pd.crosstab(scores_df.Metadata_Plate, scores_df.Metadata_batch)


# In[5]:


# Setup a dictionary to load cell count files
plate_info = dict(zip(scores_df.Metadata_Plate, scores_df.Metadata_batch))
plate_info


# In[6]:


# Load cell count data
cell_count_data = []
for plate in plate_info:
    batch = plate_info[plate]
    cell_count_file = pathlib.Path(f"{cell_count_dir}/{batch}_{plate}_cell_count.tsv")
    cell_count_df = pd.read_csv(cell_count_file, sep="\t").assign(Metadata_batch=batch)
    cell_count_data.append(cell_count_df)
    
cell_count_data = pd.concat(cell_count_data, sort=True).reset_index(drop=True)
cell_count_data.Metadata_Plate = cell_count_data.Metadata_Plate.astype(str)

print(cell_count_data.shape)
cell_count_data.head()


# In[7]:


# Merge scores and counts
scores_with_counts_df = scores_df.merge(
    cell_count_data,
    on=["Metadata_Well", "Metadata_batch", "Metadata_Plate", "Metadata_cell_density"],
    suffixes=["", "_count"]
)

# Set factors for plotting downstream
scores_with_counts_df.Metadata_model_split = pd.Categorical(
    scores_with_counts_df.Metadata_model_split,
    categories=["training", "test", "validation", "perturbation"]
)

print(scores_with_counts_df.shape)
scores_with_counts_df.head()


# In[8]:


for dataset in scores_with_counts_df.dataset.unique():
    plot_df = (
        scores_with_counts_df
        .query("Metadata_model_split != 'perturbation'")
        .query("dataset == @dataset")
    )
    
    fig_file = pathlib.Path(f"{output_fig_dir}/{dataset}_untreated_score_by_count.png")
    count_score_gg = (
        gg.ggplot(plot_df, gg.aes(y = "cell_count", x = "TotalScore"))
        + gg.geom_point(gg.aes(fill = "Metadata_clone_number", shape="Metadata_Plate"), size=2, alpha=0.8)
        + gg.facet_grid("Metadata_clone_type~Metadata_model_split")
        + gg.theme_bw()
        + gg.scale_fill_discrete(name="Clone")
        + gg.scale_shape_discrete(name="Plate")
        + gg.ylab("Untreated cell count (median replicate)\n(per well)")
        + gg.xlab("Untreated resistance score (median replicate)\n(singscore)")
        + gg.ggtitle(f"Dataset: {dataset}")
        + gg.theme(
            strip_text=gg.element_text(size=6, color="black"),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            axis_text_x=gg.element_text(rotation=90)
        )
    )
    
    count_score_gg.save(fig_file, width=6, height=4, dpi=500)
    print(count_score_gg)


# In[9]:


# Plot several diagnostic plot of cell count by well position
scores_with_counts_df = scores_with_counts_df.assign(
    row=scores_with_counts_df.Metadata_Well.str[0],
    column=scores_with_counts_df.Metadata_Well.str[1:3].astype(int)
)

scores_with_counts_df.row = pd.Categorical(
    scores_with_counts_df.row, categories=(sorted(scores_with_counts_df.row.unique(), reverse=True))
)

scores_with_counts_df.cell_count = scores_with_counts_df.cell_count.astype(float)

plate_layout_plots = []
for dataset in scores_with_counts_df.dataset.unique():
    plot_df = (
        scores_with_counts_df
        .query("Metadata_model_split != 'perturbation'")
        .query("dataset == @dataset")
    )
    for plate in plot_df.Metadata_Plate.unique():

        count_score_gg = (
            gg.ggplot(plot_df.query("Metadata_Plate == @plate"), gg.aes(y="row", x="column"))
                + gg.geom_point(gg.aes(color = "cell_count"), size = 8)
                + gg.facet_grid("Metadata_clone_type~Metadata_model_split")
                + gg.coord_fixed()
                + gg.theme_bw()
                + gg.scale_fill_discrete(name="Clone")
                + gg.ylab("Untreated cell count\n(per well)")
                + gg.xlab("Untreated resistance score\n(singscore)")
                + gg.ggtitle(f"Dataset: {dataset}; Plate: {plate}")
                + gg.theme(
                    strip_text=gg.element_text(size=6, color="black"),
                    strip_background=gg.element_rect(colour="black", fill="#fdfff4")
                )
        )
        
        plate_layout_plots.append(count_score_gg)
        
fig_file = pathlib.Path(f"{output_fig_dir}/plate_layout_summary_full.pdf")
gg.save_as_pdf_pages(plate_layout_plots, filename=fig_file)


# ## Compare untreated score to perturbation cell count drop
# 
# For each clone and treatment, we know how many cells are present before and after perturbation.
# Does the resistance score correlate with this differential?
# 
# We need to average the scores and counts across replicates, since different cells were treated.
# 
# ## Relative cell count
# 
# To quantify cell count changes upon perturbation, we simply take the following ratio:
# 
# > Number of Cells after Perturbation / Number of Cells with 0.1% DMSO treatment
# 
# A value above 1 indicates that the clone was resistant to treatment.
# A value significantly below 1 indicates that the clone was sensitive.
# 
# This metric controls for plate seeding density differences, and, since it is an median value, will minimize the within plate across well differences.
# 
# ## We perform this analysis differently depending on the dataset
# 
# The CloneAE data contained treated wells in the same plate as untreated wells.
# The Ixazomib and CB-5083 treatments were collected in different plates.
# 
# ### Step 1: CloneAE Dataset

# In[10]:


# Get the median score and cell count per clone and model splits
scores_with_counts_df.Metadata_model_split = scores_with_counts_df.Metadata_model_split.astype(str)

cloneae_summary_df = (
    scores_with_counts_df.query("dataset == 'cloneAE'")
    .groupby(
        [
            "Metadata_batch",
            "Metadata_Plate",
            "Metadata_clone_number",
            "Metadata_treatment",
            "Metadata_model_split",
        ]
    )
    .agg(
        {
            "TotalScore": "median",
            "cell_count": "median"
        }
    )
    .reset_index()
)

# Set factors for plotting downstream
cloneae_summary_df.Metadata_model_split = pd.Categorical(
    cloneae_summary_df.Metadata_model_split,
    categories=["training", "test", "validation", "perturbation"]
)

print(cloneae_summary_df.shape)
cloneae_summary_df.head()


# In[11]:


# Merge together treated and untreated cell counts
clonae_perturbation_summary_df = (
    cloneae_summary_df
    .query("Metadata_treatment != '0.1% DMSO'")
    .reset_index(drop=True)
)
clonae_untreated_summary_df = (
    cloneae_summary_df
    .query("Metadata_treatment == '0.1% DMSO'")
    .reset_index(drop=True)
)

cloneae_result = (
    clonae_perturbation_summary_df
    .merge(
        clonae_untreated_summary_df,
        on=["Metadata_Plate", "Metadata_batch", "Metadata_clone_number"],
        how="left",
        suffixes=["_perturbed", "_untreated"]
    )
)

# Generate a simple relative cell count metric
cloneae_result = cloneae_result.assign(
    relative_cell_count=cloneae_result.cell_count_perturbed / cloneae_result.cell_count_untreated
)

print(cloneae_result.shape)
cloneae_result.head()


# In[12]:


cloneAE_score_count_gg = (
    gg.ggplot(cloneae_result, gg.aes(x="TotalScore_untreated", y="relative_cell_count"))
    + gg.geom_point(gg.aes(fill="Metadata_clone_number", shape="Metadata_Plate"), size=2, alpha=0.8)
    + gg.xlab("Untreated resistance signature (median replicate)\n(singscore)")
    + gg.ylab("Relative cell count (median replicate)\n (Compared to untreated)")
    + gg.ggtitle("Dataset: CloneAE; Bortezomib treatment")
    + gg.geom_hline(yintercept=1, linetype="dashed", color="red")
    + gg.scale_fill_discrete(name="Clone")
    + gg.scale_shape_discrete(name="Plate")
    + gg.facet_grid("Metadata_model_split_untreated~Metadata_treatment_perturbed")
    + gg.theme_bw()
    + gg.theme(
        strip_text=gg.element_text(size=6, color="black"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4")
    )
)

fig_file = pathlib.Path(f"{output_fig_dir}/cloneAE_treatment_score_by_count.png")
cloneAE_score_count_gg.save(fig_file, width=5, height=4, dpi=500)

cloneAE_score_count_gg


# ### Step 2: Ixazomib and CB-5083
# 
# Both of these cell lines were acquired with a different plate layout than the cloneAE dataset.
# Specifically, each plate were subject to indepedent treatments: DMSO and each drug treatment were collected on separate plates.
# 
# To compare cell counts before and after treatment, we must assume similar seeding density (this is a relatively safe assumption, since this info was captured) and compare treated plates to DMSO plates.
# To protect against seeding density artifacts, we compare the median cell count with the median resistance signature.

# In[13]:


# Load treatment plates
treatment_plate_info = {
    '218775': "2020_08_24_Batch9",
    '218699': '2020_08_24_Batch9',
    '218697': '2020_08_24_Batch9',
    '218853': '2020_09_08_Batch10',
    '218855': '2020_09_08_Batch10',
    '218857': '2020_09_08_Batch10',
    '218859': '2020_09_08_Batch10',
}

treatment_cell_count_data = []
for plate in treatment_plate_info:
    batch = treatment_plate_info[plate]
    cell_count_file = pathlib.Path(f"{cell_count_dir}/{batch}_{plate}_cell_count.tsv")
    cell_count_df = pd.read_csv(cell_count_file, sep="\t").assign(Metadata_batch=batch)
    treatment_cell_count_data.append(cell_count_df)
    
treatment_cell_count_data = pd.concat(treatment_cell_count_data, sort=True).reset_index(drop=True)
treatment_cell_count_data.Metadata_Plate = treatment_cell_count_data.Metadata_Plate.astype(str)

treatment_cell_count_summary_df = (
    treatment_cell_count_data.groupby(
        ["Metadata_batch", "Metadata_treatment", "Metadata_clone_number", "Metadata_cell_density"]
    ).agg(
        {
            "cell_count": "median"
        }
    )
    .reset_index()
    .assign(dataset="cb5083")
)

treatment_cell_count_summary_df.loc[
        treatment_cell_count_summary_df.Metadata_treatment.str.contains("Ixazomib"), "dataset"
] = "ixazomib"


print(treatment_cell_count_summary_df.shape)
treatment_cell_count_summary_df.head()


# In[14]:


# Calculate median resistance scores and merge with perturbation cell counts
other_treatment_summary_df = (
    scores_with_counts_df.query("dataset != 'cloneAE'")
    .groupby(
        [
            "Metadata_batch",
            "Metadata_Plate",
            "Metadata_clone_number",
            "Metadata_clone_type",
            "Metadata_treatment",
            "Metadata_model_split",
            "Metadata_cell_density",
            "dataset"
        ]
    )
    .agg(
        {
            "TotalScore": "median",
            "cell_count": "median"
        }
    )
    .reset_index()
)


other_treatment_summary_df = other_treatment_summary_df.merge(
    treatment_cell_count_summary_df,
    on=["Metadata_batch", "Metadata_clone_number", "Metadata_cell_density", "dataset"],
    suffixes=["_untreated", "_perturbed"]
)

other_treatment_summary_df = other_treatment_summary_df.assign(
    relative_cell_count=(
        other_treatment_summary_df.cell_count_perturbed / other_treatment_summary_df.cell_count_untreated
    )
)

other_treatment_summary_df.Metadata_model_split = pd.Categorical(
    other_treatment_summary_df.Metadata_model_split,
    categories=["training", "test", "validation"]
)

print(other_treatment_summary_df.shape)
other_treatment_summary_df.head()


# In[15]:


for dataset in other_treatment_summary_df.dataset.unique():

    score_count_gg = (
        gg.ggplot(other_treatment_summary_df.query("dataset == @dataset"),
                  gg.aes(x="TotalScore", y="relative_cell_count"))
        + gg.geom_point(gg.aes(fill="Metadata_clone_type", shape="Metadata_Plate"), size=2, alpha=0.8)
        + gg.coord_fixed()
        + gg.xlab("Untreated resistance signature (median replicate)\n(singscore)")
        + gg.ylab("Relative cell count (median replicate)\n (Compared to untreated)")
        + gg.ggtitle(f"Dataset: {dataset} treatment")
        + gg.geom_hline(yintercept=1, linetype="dashed", color="red")
        + gg.facet_grid("Metadata_treatment_perturbed~Metadata_model_split")
        + gg.theme_bw()
        + gg.theme(
            strip_text=gg.element_text(size=6, color="black"),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4")
        )
    )
    
    fig_file = pathlib.Path(f"{output_fig_dir}/{dataset}_treatment_score_by_count.png")
    score_count_gg.save(fig_file, width=5, height=4, dpi=500)

    print(score_count_gg)

