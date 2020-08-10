#!/usr/bin/env python
# coding: utf-8

# In[1]:


import umap
import pathlib
import warnings
import pandas as pd
import plotnine as gg

from numba.core.errors import NumbaWarning

from utils.data_utils import load_data


# In[2]:


def apply_umap(x_df, meta_df):
    reducer = umap.UMAP(random_state=123)
    embedding_df = reducer.fit_transform(x_df)

    # Setup plotting logic
    embedding_df = pd.DataFrame(embedding_df, columns=['x', 'y'])
    embedding_df = embedding_df.merge(meta_df, left_index=True, right_index=True)
    
    return embedding_df


def plot_umap_cell_line(embedding_df, fig_file, cell_line_column, color_labels, color_values):
    cell_line_gg = (
        gg.ggplot(embedding_df, gg.aes(x="x", y="y")) +
        gg.geom_point(gg.aes(color=cell_line_column), size = 0.2, shape = ".", alpha = 0.2) +
        gg.theme_bw() +
        gg.scale_color_manual(name="Cell Line", labels=color_labels, values=color_values)
        )

    cell_line_gg.save(filename=fig_file, height=4, width=5, dpi=500)
    return cell_line_gg

    
def plot_umap_well(embedding_df, fig_file, well_column):
    well_gg = (
        gg.ggplot(embedding_df, gg.aes(x="x", y="y")) +
        gg.geom_point(gg.aes(color=well_column), size = 0.2, shape = ".", alpha = 0.2) +
        gg.theme_bw() 
    )

    well_gg.save(filename=fig_file, height=4, width=5, dpi=500)
    return well_gg


# ## Load Data

# In[3]:


data_dict = load_data(
    return_meta=True,
    shuffle_row_order=True,
    holdout=True,
    othertreatment=True
)

print(data_dict["train"]["x"].shape)
print(data_dict["test"]["x"].shape)
print(data_dict["holdout"]["x"].shape)
print(data_dict["othertreatment"]["x"].shape)

data_dict["test"]["x"].head(3)


# In[4]:


warnings.filterwarnings("ignore", category=NumbaWarning)


# ## Apply and visualize UMAP

# In[5]:


embedding_dict = {}

embedding_dict["train"] = apply_umap(data_dict["train"]["x"], meta_train_df)
embedding_dict["test"] = apply_umap(data_dict["test"]["x"], meta_test_df)
embedding_dict["holdout"] = apply_umap(data_dict["holdout"]["x"], meta_holdout_df)
embedding_dict["othertreatment"] = apply_umap(data_dict["othertreatment"]["x"], meta_other_df)


# In[6]:


cell_line_column = "Metadata_clone_number"
well_column = "Metadata_Well"

cell_line_labels = {"Clone A": "Clone A", "Clone E": "Clone E", "WT parental": "WT parental"}
cell_line_colors = {"Clone A": "#1b9e77", "Clone E": "#d95f02", "WT parental": "#7570b3"}


# In[7]:


for data_fit, embedding_df in embedding_dict.items():
    fig_file = pathlib.Path("figures", "umap", f"single_cell_umap_{data_fit}.png")
    cell_gg = plot_umap_cell_line(
        embedding_df, fig_file, cell_line_column, cell_line_labels, cell_line_colors
    )
    print(cell_gg)
    
    fig_file = pathlib.Path("figures", "umap", f"single_cell_umap_well_{data_fit}.png")
    well_gg = plot_umap_well(embedding_df, fig_file, well_column)
    print(well_gg)


# In[8]:


treatment_gg = (
    gg.ggplot(embedding_dict["othertreatment"], gg.aes(x="x", y="y")) +
    gg.geom_point(gg.aes(color="Metadata_treatment"), size = 0.2, shape = ".", alpha = 0.2) +
    gg.theme_bw() 
)

fig_file = pathlib.Path("figures", "umap", "single_cell_othertreatment.png")
treatment_gg.save(filename=fig_file, height=4, width=5, dpi=500)

treatment_gg

