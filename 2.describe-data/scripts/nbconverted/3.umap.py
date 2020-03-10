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


# In[6]:


# Visualize UMAP results
p = (
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
    
p


# In[7]:


file = os.path.join("figures", "umap", "four_clone_umap")
for extension in ['.png', '.pdf', '.svg']:
    gg.ggsave(p, filename='{}{}'.format(file, extension), height=6, width=7)


# ## For Clone A and E Data

# In[8]:


# Load and process data
file = os.path.join("data", "merged", "combined_cloneAcloneE_dataset.csv")
cloneAE_data_df = pd.read_csv(file)

embedding_df = process_umap(cloneAE_data_df)
embedding_df.head()


# In[9]:


# Visualize UMAP results
p = (
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

p


# In[10]:


file = os.path.join("figures", "umap", "merged_umap_clone_ae")
for extension in ['.png', '.pdf', '.svg']:
    gg.ggsave(p, filename='{}{}'.format(file, extension), height=6, width=7)


# ## Merged Data

# In[11]:


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


# In[12]:


combined_df = pd.concat([data_recode_df, cloneAE_data_recode_df], sort=True).reset_index(drop=True)
combined_df = feature_select(combined_df, operation="drop_na_columns")

print(combined_df.shape)
combined_df.head()


# In[13]:


embedding_df = process_umap(combined_df)
embedding_df.head()


# In[14]:


# Visualize UMAP results
p = (
    gg.ggplot(embedding_df,
              gg.aes('x', 'y', color="Metadata_treatment", size='Metadata_Dataset'))
    + gg.geom_point()
    + gg.theme_bw()
    + gg.xlab("x")
    + gg.ylab("y")
    + gg.ggtitle("Cell Painting Merged (UMAP)")
    + gg.theme(legend_key=gg.element_rect(color="black", fill = "white"))
)

p


# In[15]:


file = os.path.join("figures", "umap", "clone_compare_batch_effect")
for extension in ['.png', '.pdf', '.svg']:
    gg.ggsave(p, filename='{}{}'.format(file, extension), height=6, width=7)

