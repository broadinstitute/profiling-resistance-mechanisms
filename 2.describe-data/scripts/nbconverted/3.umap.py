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


# ## For Combined Batches of Four WT + Resistant Clones

# In[5]:


# Load and process data
file = os.path.join("data", "merged", "combined_four_clone_dataset.csv")
data_df = pd.read_csv(file)

embedding_df = process_umap(data_df)
embedding_df.head()


# ## Visualize a Series of UMAP Representations

# In[6]:


umap_resistant_type_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
    + gg.geom_point(
        gg.aes(fill="Metadata_clone_type", shape="Metadata_treatment"),
        color='black', alpha=0.8)
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.scale_shape_manual(name="Treatment", values=[".", "+"])
)

file = os.path.join("figures", "umap", "four_clone_umap_resistant_type")
for extension in ['.png', '.pdf', '.svg']:
    umap_resistant_type_gg.save(filename='{}{}'.format(file, extension), height=3, width=4, dpi=400)

umap_resistant_type_gg


# In[7]:


umap_batch_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
    + gg.geom_point(gg.aes(fill="Metadata_batch"), color='black', alpha=0.8)
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
)

file = os.path.join("figures", "umap", "four_clone_umap_batch")
for extension in ['.png', '.pdf', '.svg']:
    umap_batch_gg.save(filename='{}{}'.format(file, extension), height=3, width=4, dpi=400)

umap_batch_gg


# In[8]:


umap_batch_facet_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
    + gg.geom_point(gg.aes(fill="Metadata_batch"), color='black', alpha=0.8, size=0.8)
    + gg.theme_bw()
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
    + gg.facet_wrap("~Metadata_plate_ID")
    + gg.theme(
        strip_text=gg.element_text(size=6, color="black"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
)
    
file = os.path.join("figures", "umap", "four_clone_umap_plate_facet")
for extension in ['.png', '.pdf', '.svg']:
    umap_batch_facet_gg.save(filename='{}{}'.format(file, extension), height=3, width=4, dpi=400)

umap_batch_facet_gg


# In[9]:


# Visualize UMAP results
clone_facet_gg = (
    gg.ggplot(embedding_df,
              gg.aes('x', 'y', color='factor(Metadata_Plate)', shape="Metadata_treatment"))
    + gg.geom_point()
    + gg.theme_bw()
    + gg.xlab("UMAP X")
    + gg.ylab("UMAP Y")
    + gg.scale_shape_manual(name="Treatment", values=[".", "+"])
    + gg.scale_color_discrete(name="Plate")
    + gg.facet_wrap("~Metadata_clone_number")
    + gg.ggtitle("Cell Painting Four Clones")
    + gg.theme(
        legend_key=gg.element_rect(color="black", fill = "white"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4")
    )
)
    
file = os.path.join("figures", "umap", "four_clone_umap_facet_clone_sample")
for extension in ['.png', '.pdf', '.svg']:
    clone_facet_gg.save(filename='{}{}'.format(file, extension), height=3, width=4, dpi=400)

clone_facet_gg


# In[10]:


umap_well_embedding_gg = (
    gg.ggplot(embedding_df, gg.aes(x="x", y="y"))
    + gg.geom_point(
        gg.aes(fill="Metadata_clone_number", shape="Metadata_treatment"),
        color='black', alpha=0.8
    )
    + gg.theme_bw()
    + gg.scale_shape_manual(name="Treatment", values=[".", "+"])
    + gg.xlab("UMAP (X)")
    + gg.ylab("UMAP (Y)")
)

file = os.path.join("figures", "umap", "four_clone_umap_clone_sample")
for extension in ['.png', '.pdf', '.svg']:
    umap_well_embedding_gg.save(filename='{}{}'.format(file, extension), height=3, width=4, dpi=400)

umap_well_embedding_gg


# ## For Clone A and E Data

# In[11]:


# Load and process data
file = os.path.join("data", "merged", "combined_cloneAcloneE_dataset.csv")
cloneAE_data_df = pd.read_csv(file)

embedding_df = process_umap(cloneAE_data_df)
embedding_df.head()


# In[12]:


# Visualize UMAP results
clone_ae_umap_gg = (
    gg.ggplot(embedding_df,
              gg.aes('x', 'y',
                     shape="Metadata_Plate", 
                     size='factor(Metadata_Dosage)',
                     color="Metadata_CellLine"))
    + gg.geom_point()
    + gg.theme_bw()
    + gg.scale_shape_manual(name="Plate", values=[".", "+"])
    + gg.scale_color_discrete(name="Clone")
    + gg.scale_size_discrete(name="Dosage")
    + gg.xlab("UMAP X")
    + gg.ylab("UMAP Y")
    + gg.ggtitle("Cell Painting Clone A and E")
    + gg.theme(
        legend_key=gg.element_rect(color="black", fill = "white"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4")
    )
)

file = os.path.join("figures", "umap", "cloneAE_umap")
for extension in ['.png', '.pdf', '.svg']:
    clone_ae_umap_gg.save(filename='{}{}'.format(file, extension), height=3, width=4, dpi=400)

clone_ae_umap_gg


# ## Merged Data

# In[13]:


drop_cols = ["Metadata_plate_ID", "Metadata_plate_filename"]

data_recode_df = (
    data_df
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


# In[14]:


combined_df = pd.concat([data_recode_df, cloneAE_data_recode_df], sort=True).reset_index(drop=True)
combined_df = feature_select(combined_df, operation="drop_na_columns")

print(combined_df.shape)
combined_df.head()


# In[15]:


embedding_df = process_umap(combined_df)
embedding_df.head()


# In[16]:


# Visualize UMAP results
merged_umap_gg = (
    gg.ggplot(
        embedding_df,
        gg.aes(
            'x', 'y',
            color="Metadata_treatment",
            shape='Metadata_Dataset'
        )
    )
    + gg.geom_point()
    + gg.theme_bw()
    + gg.xlab("x")
    + gg.ylab("y")
    + gg.scale_shape_manual(name="Plate", values=[".", "+"])
    + gg.ggtitle("Cell Painting Merged (UMAP)")
    + gg.theme(legend_key=gg.element_rect(color="black", fill = "white"))
)

file = os.path.join("figures", "umap", "clone_compare_batch_effect")
for extension in ['.png', '.pdf', '.svg']:
    merged_umap_gg.save(filename='{}{}'.format(file, extension), height=3, width=4, dpi=400)

merged_umap_gg

