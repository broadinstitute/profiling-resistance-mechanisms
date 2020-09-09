#!/usr/bin/env python
# coding: utf-8

# # Apply and Visualize UMAP
# 
# **Gregory Way, 2019**
# 
# We are interested in visualizing the relationship among samples according to several variables.
# These variables include `batch`, `dosage`, and `cell line`.

# In[1]:


import os
import numpy as np
import pandas as pd
import umap

import plotnine as gg

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


np.random.seed(123)


# In[4]:


def process_umap(data_df):    
    # Prepare UMAP input by removing metadata columns
    metadata_cols = infer_cp_features(data_df, metadata=True)

    metadata_df = data_df.loc[:, metadata_cols]
    umap_data_df = data_df.drop(metadata_cols, axis="columns")
    
    # Apply UMAP
    reducer = umap.UMAP(random_state=123)
    embedding = reducer.fit_transform(umap_data_df)
    
    # Setup plotting logic
    embedding_df = pd.DataFrame(embedding, columns=['x', 'y'])
    embedding_df = embedding_df.merge(metadata_df, left_index=True, right_index=True)
    
    return embedding_df


# In[5]:


save_file_extensions = ['.png']
data_dir = os.path.join("data", "merged")


# ## Load Data

# In[6]:


# Load and process data
file = os.path.join(data_dir, "all_merged_profiles.csv.gz")
data_df = pd.read_csv(file)

embedding_df = process_umap(data_df)
embedding_df.head()


# ## Visualize a Series of UMAP Representations

# In[13]:


umap_resistant_type_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
    + gg.geom_point(
        gg.aes(shape="Metadata_clone_type", fill="Metadata_treatment"),
        color='black', alpha=0.6)
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.ggtitle("All Batches")
    + gg.scale_shape_manual(name="Clone Type", values=[".", "+"])
)

file = os.path.join("figures", "umap", "umap_resistant_type")
for extension in save_file_extensions:
    umap_resistant_type_gg.save(filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400)

umap_resistant_type_gg


# In[14]:


umap_cell_count_gg = (
    gg.ggplot(
        embedding_df.rename({"Metadata_cell_count": "Cell Count"}, axis="columns")
    )
    + gg.geom_point(
        gg.aes(x="x", y="y", fill="Cell Count"),
        color='black', alpha=0.6)
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.ggtitle("All Batches")
)

file = os.path.join("figures", "umap", "umap_cell_count")
for extension in save_file_extensions:
    umap_cell_count_gg.save(filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400)
    
umap_cell_count_gg


# In[16]:


umap_batch_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
    + gg.geom_point(
        gg.aes(fill="Metadata_batch"),
        color='black',
        alpha=0.6
    )
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.ggtitle("All Batches")
)

file = os.path.join("figures", "umap", "umap_batch")
for extension in save_file_extensions:
    umap_batch_gg.save(filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400)

umap_batch_gg


# In[20]:


umap_batch_facet_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
     + gg.geom_point(
        gg.aes(fill="Metadata_treatment", shape="Metadata_clone_type"),
        color='black',
        alpha=0.8,
        size=0.6
    )
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.ggtitle("All Batches")
    + gg.facet_wrap("~Metadata_Plate")
    + gg.scale_shape_manual(name="Clone Type", values=[".", "+"])
    + gg.theme(
        strip_text=gg.element_text(size=6, color="black"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
)
    
file = os.path.join("figures", "umap", "umap_plate_facet")
for extension in save_file_extensions:
    umap_batch_facet_gg.save(filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400)

umap_batch_facet_gg


# In[21]:


# Visualize UMAP results
clone_facet_gg = (
    gg.ggplot(embedding_df, gg.aes('x', 'y'))
    + gg.geom_point(
        gg.aes(fill="Metadata_treatment", shape="Metadata_clone_type"),
        alpha=0.6
    )
    + gg.theme_bw()
    + gg.xlab("UMAP X")
    + gg.ylab("UMAP Y")
    + gg.scale_shape_manual(name="Clone Type", values=[".", "+"])
    + gg.scale_fill_discrete(name="Treatment")
    + gg.facet_wrap("~Metadata_clone_number")
    + gg.ggtitle("All Batches")
    + gg.theme(
        legend_key=gg.element_rect(color="black", fill = "white"),
        strip_text=gg.element_text(size=6, color="black"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4")
    )
)
    
file = os.path.join("figures", "umap", "umap_facet_clone")
for extension in save_file_extensions:
    clone_facet_gg.save(filename='{}{}'.format(file, extension), height=4, width=4.5, dpi=400)

clone_facet_gg


# In[24]:


umap_well_embedding_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
    + gg.geom_point(
        gg.aes(fill="Metadata_batch", shape="Metadata_clone_type"),
        color='black',
        alpha=0.6
    )
    + gg.theme_bw()
    + gg.scale_shape_manual(name="Clone Type", values=[".", "+"])
    + gg.scale_fill_discrete(name="Clone")
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.ggtitle("All Batches")
)

file = os.path.join("figures", "umap", "four_clone_umap_clone_sample")
for extension in save_file_extensions:
    umap_well_embedding_gg.save(filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400)

umap_well_embedding_gg


# ## For Clone A and E Data

# In[13]:


# Load and process data
file = os.path.join("data", "merged", "combined_cloneAcloneE_dataset_feature_select.csv.gz")
cloneAE_data_df = pd.read_csv(file)

embedding_cloneAE_df = process_umap(cloneAE_data_df)
embedding_cloneAE_df.head()


# In[14]:


# Visualize UMAP results
clone_ae_umap_gg = (
    gg.ggplot(embedding_cloneAE_df)
    + gg.geom_point(
        gg.aes('x', 'y',
               shape="Metadata_Plate", 
               size='factor(Metadata_Dosage)',
               color="Metadata_CellLine"),
        alpha=0.8
    )
    + gg.theme_bw()
    + gg.scale_shape_manual(name="Plate", values=[".", "+"])
    + gg.scale_color_discrete(name="Clone")
    + gg.scale_size_manual(name="Dosage", values=[1, 3, 5, 7])
    + gg.xlab("UMAP X")
    + gg.ylab("UMAP Y")
    + gg.ggtitle("Clone A and E - Merged")
    + gg.theme(
        legend_key=gg.element_rect(color="black", fill = "white"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4")
    )
)

file = os.path.join("figures", "umap", "cloneAE_umap")
for extension in save_file_extensions:
    clone_ae_umap_gg.save(
        filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400
    )

clone_ae_umap_gg


# In[15]:


clone_ae_umap_cell_count_gg = (
    gg.ggplot(
        embedding_cloneAE_df.rename({"Metadata_cell_count": "Cell Count"}, axis="columns")
    )
    + gg.geom_point(
        gg.aes(x="x", y="y", fill="Cell Count"),
        color='black', alpha=0.6)
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.ggtitle("Clone A and E - Merged")
)

file = os.path.join("figures", "umap", "cloneAE_umap_cell_count")
for extension in save_file_extensions:
    clone_ae_umap_cell_count_gg.save(
        filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400
    )
    
clone_ae_umap_cell_count_gg


# ## Merged Data

# In[16]:


drop_cols = ["Metadata_plate_ID", "Metadata_plate_filename"]

fourclone_data_recode_df = (
    fourclone_data_df
    .drop(drop_cols, axis="columns")
    .rename(
        {
            "Metadata_clone_number": "Metadata_CellLine"
        }, axis="columns"
    )
    .assign(Metadata_Dosage=0.7)
    .assign(Metadata_Dataset="FourClone")
)

cloneAE_data_recode_df = (
    cloneAE_data_df.assign(Metadata_treatment="bortezomib")
    .assign(Metadata_Dataset="CloneAE")
)

cloneAE_data_recode_df.loc[cloneAE_data_recode_df.Metadata_Dosage == 0, "Metadata_treatment"] = "DMSO"


# In[17]:


combined_df = pd.concat([fourclone_data_recode_df, cloneAE_data_recode_df], sort=True).reset_index(drop=True)
combined_df = feature_select(combined_df, operation="drop_na_columns")

print(combined_df.shape)
combined_df.head()


# In[18]:


embedding_combined_df = process_umap(combined_df)
embedding_combined_df.head()


# In[19]:


# Visualize UMAP results
merged_umap_gg = (
    gg.ggplot(embedding_combined_df)
    + gg.geom_point(
        gg.aes(
            'x', 'y',
            shape="Metadata_treatment",
            fill='Metadata_Dataset'
        ),
        size=1,
        alpha=0.4
    )
    + gg.theme_bw()
    + gg.xlab("x")
    + gg.ylab("y")
    + gg.scale_shape_manual(name="Treatment", values=["^", "+"])
    + gg.scale_fill_manual(name="Dataset", values=["#F0C13E", "#BD4135"])
    + gg.ggtitle("Cell Painting Merged (UMAP)")
    + gg.theme(legend_key=gg.element_rect(color="black", fill = "white"))
)

file = os.path.join("figures", "umap", "clone_compare_batch_effect")
for extension in save_file_extensions:
    merged_umap_gg.save(filename='{}{}'.format(file, extension), height=3, width=3.5, dpi=400)

merged_umap_gg

