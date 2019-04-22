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


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


np.random.seed(123)


# In[4]:


# Load data
file = os.path.join("data", "merged_intersected_variable_selected.csv")
data_df = pd.read_csv(file)

# Remove metadata columns
metadata_columns = data_df.columns.str.contains("Metadata_")

metadata_df = data_df.loc[:, metadata_columns]
umap_data_df = data_df.loc[:, ~metadata_columns]

print(umap_data_df.shape)
umap_data_df.head()


# In[5]:


# Apply UMAP
reducer = umap.UMAP(random_state=123)
embedding = reducer.fit_transform(umap_data_df)


# In[6]:


# Setup plotting logic
embedding_df = pd.DataFrame(embedding, columns=['x', 'y'])
embedding_df = embedding_df.merge(metadata_df, left_index=True, right_index=True)

embedding_df.head(3)


# In[7]:


# Visualize UMAP results
p = gg.ggplot(embedding_df, gg.aes('x',
                                    'y',
                                    shape='factor(Metadata_Batch_Number)',
                                    color='Metadata_CellLine',
                                    size='factor(Metadata_Dosage)')) + \
    gg.geom_point() + \
    gg.theme_bw() + \
    gg.xlab("x") + \
    gg.ylab("y") + \
    gg.ggtitle("Cell Painting Merged (UMAP)") + \
    gg.scale_color_manual(name="CellLine", values=["#1b9e77", "#d95f02", "#7570b3"]) + \
    gg.scale_shape_manual(name="Batch", labels=[1, 2], values=['o', '+']) + \
    gg.scale_size_manual(name="Dosage", values=[2, 3, 4, 5]) + \
    gg.theme(legend_key=gg.element_rect(color="black", fill = "white"))

p


# In[8]:


file = os.path.join("figures", "merged_umap")
for extension in ['.png', '.pdf']:
    gg.ggsave(p, filename='{}{}'.format(file, extension), height=6, width=7)

