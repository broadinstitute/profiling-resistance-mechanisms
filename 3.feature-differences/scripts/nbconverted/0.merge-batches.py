#!/usr/bin/env python
# coding: utf-8

# # Merge Select Batches Together
# 
# Merge batches 5, 6, and 7 together.
# These batches contain 5 total plates that all have roughly the same platemap.
# 
# There are 8 unique samples measured on each plate in 3 replicates each, plus one sample measured on each plate in 9 replicates.
# The 8 samples are: `BZ001`, `BZ008`, `BZ017`, `BZ018`, `WT002`, `WT008`, `WT009`, `WT011` and `WT_parental`.
# All samples have been treated with both DMSO (control) and 7nm Bortezomib.

# In[1]:


import os
import numpy as np
import pandas as pd
import umap
import plotnine as gg

from pycytominer import feature_select
from pycytominer import write_gct
from pycytominer.cyto_utils import infer_cp_features
from pycytominer import audit
from pycytominer import aggregate


# In[2]:


np.random.seed(123)


# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[4]:


def load_data(batch,
              profile_dir="profiles",
              file_select="normalized_variable_selected"):
    batch_dir = os.path.join(profile_dir, batch)

    backend_folders = os.listdir(batch_dir)
    plate_files = [
        os.path.join(batch_dir, x, "{}_{}.csv".format(x, file_select))
        for x in backend_folders if ".DS_Store" not in x
    ]
    
    plate_dfs = []
    for plate_file in plate_files:
        plate_dfs.append(pd.read_csv(plate_file))
        
    plate_df = pd.concat(plate_dfs, axis="rows")
    plate_df = plate_df.assign(Metadata_batch=batch)
    
    # Reorder columns
    meta_col_bool = plate_df.columns.str.contains("Metadata_")
    meta_cols = plate_df.columns[meta_col_bool].tolist()
    cp_cols =  plate_df.columns[~meta_col_bool].tolist()
    col_order = meta_cols + cp_cols
    
    return plate_df.loc[:, col_order]


# In[5]:


profile_dir = os.path.join("..", "1.process-profiles", "profiles")


# In[6]:


# Only use certain batches
# These batches have the same platemap
use_batches = ["2019_11_19_Batch5", "2019_11_20_Batch6", "2019_11_22_Batch7"]


# In[7]:


# Merge batches
batch_df = list()
for batch in use_batches:
    data = load_data(batch, profile_dir=profile_dir, file_select="normalized")
    batch_df.append(data)
    
batch_df = pd.concat(batch_df, axis="rows")

# Add column of clonal resistance
batch_df = batch_df.assign(Metadata_clone_type="resistant")
batch_df.loc[batch_df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type"] = "wildtype"

meta_col_bool = batch_df.columns.str.contains("Metadata_")
meta_cols = batch_df.columns[meta_col_bool].tolist()
cp_cols =  batch_df.columns[~meta_col_bool].tolist()
col_order = meta_cols + cp_cols
batch_df = batch_df.loc[:, col_order]

print(batch_df.shape)
batch_df.head()


# In[8]:


pd.crosstab(batch_df.Metadata_clone_number, batch_df.Metadata_treatment)


# In[9]:


# Apply feature selection using all features
batch_feature_select_df = feature_select(
    profiles=batch_df,
    operation="drop_na_columns"
)

ops = [
    "variance_threshold",
    "correlation_threshold",
    "blacklist",
    "drop_outliers"
]

batch_feature_select_df = feature_select(
    profiles=batch_feature_select_df,
    operation=ops
)


# In[10]:


print(batch_feature_select_df.shape)
batch_feature_select_df.head()


# In[11]:


output_real_file = os.path.join("data", "core_batch_profiles.tsv")
batch_feature_select_df.to_csv(output_real_file, sep="\t", index=False)


# # Generate consensus profiles by well

# In[12]:


consensus_well_df = (
    aggregate(batch_feature_select_df,
              strata="Metadata_well_position",
              operation="median")
)

consensus_well_df = (
    batch_feature_select_df.loc[:,
                                [
                                    "Metadata_well_position",
                                    "Metadata_clone_number",
                                    "Metadata_clone_type",
                                    "Metadata_treatment"
                                ]
                               ]
    .merge(
        consensus_well_df,
        on="Metadata_well_position"
    )
    .drop_duplicates()
    .reset_index(drop=True)
)

consensus_well_df.head()


# ## Two methods to visualize if batch effect correction is required

# In[13]:


# Step 1: Output GCT file
output_file = os.path.join("data", "core_batch.gct")
write_gct(
    profiles=batch_feature_select_df,
    output_file=output_file
)

# Step 1: Output GCT file
output_file = os.path.join("data", "core_batch_well_collapsed.gct")
write_gct(
    profiles=consensus_well_df,
    output_file=output_file
)


# In[14]:


# Step 2: Apply UMAP
cp_features = infer_cp_features(batch_feature_select_df)

reducer = umap.UMAP(random_state=1234, n_components=2)

metadata_df = batch_feature_select_df.drop(cp_features, axis="columns").reset_index(drop=True)

real_embedding_df = pd.DataFrame(
    reducer.fit_transform(batch_feature_select_df.loc[:, cp_features]),
    columns=["umap_x", "umap_y"]
).reset_index(drop=True)

real_embedding_df = (
    metadata_df
    .merge(real_embedding_df,
           left_index=True,
           right_index=True)
)

print(real_embedding_df.shape)
real_embedding_df.head()


# In[15]:


output_real_file = os.path.join("data", "core_batch_umap.tsv")
real_embedding_df.to_csv(output_real_file, sep="\t", index=False)


# In[16]:


reducer = umap.UMAP(random_state=1234, n_components=2)

metadata_df = consensus_well_df.drop(cp_features, axis="columns").reset_index(drop=True)

real_well_embedding_df = pd.DataFrame(
    reducer.fit_transform(consensus_well_df.loc[:, cp_features]),
    columns=["umap_x", "umap_y"]
).reset_index(drop=True)

real_well_embedding_df = (
    metadata_df
    .merge(real_well_embedding_df,
           left_index=True,
           right_index=True)
)

print(real_well_embedding_df.shape)
real_well_embedding_df.head()


# In[17]:


output_real_file = os.path.join("data", "core_batch_well_collapsed_umap.tsv")
real_well_embedding_df.to_csv(output_real_file, sep="\t", index=False)


# ## Visualize UMAP

# In[18]:


umap_batch_gg = (
    gg.ggplot(real_embedding_df, gg.aes(x="umap_x", y="umap_y")) +
    gg.geom_point(gg.aes(fill="Metadata_batch"), color='black', alpha=0.8) +
    gg.theme_bw() +
    gg.xlab("UMAP (X)") +
    gg.ylab("UMAP (Y)")
)

out_file = os.path.join("figures", "core_batch_umap_batch.png")
umap_batch_gg.save(out_file, width = 4, height = 3, dpi = 400)

umap_batch_gg


# In[19]:


umap_batch_facet_gg = (
    gg.ggplot(real_embedding_df, gg.aes(x="umap_x", y="umap_y")) +
    gg.geom_point(gg.aes(fill="Metadata_batch"), color='black', alpha=0.8, size=0.8) +
    gg.theme_bw() +
    gg.xlab("UMAP (X)") +
    gg.ylab("UMAP (Y)") +
    gg.facet_wrap("~Metadata_plate_ID") +
    gg.theme(
            strip_text=gg.element_text(size=6, color="black"),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
)

out_file = os.path.join("figures", "core_batch_umap_batch_facet.png")
umap_batch_facet_gg.save(out_file, width=4, height=3, dpi=400, verbose=False)

umap_batch_facet_gg


# In[20]:


umap_resistant_type_gg = (
    gg.ggplot(real_embedding_df, gg.aes(x="umap_x", y="umap_y")) +
    gg.geom_point(gg.aes(fill="Metadata_clone_type", shape="Metadata_treatment"),
                  color='black', alpha=0.8) +
    gg.theme_bw() +
    gg.xlab("UMAP (X)") +
    gg.ylab("UMAP (Y)")
)

out_file = os.path.join("figures", "core_batch_umap_resistant_type.png")
umap_resistant_type_gg.save(out_file, width=4, height=3, dpi=400, verbose=False)

umap_resistant_type_gg


# ## Visualize Consensus Treatments across Wells

# In[21]:


umap_well_embedding_gg = (
    gg.ggplot(real_well_embedding_df, gg.aes(x="umap_x", y="umap_y")) +
    gg.geom_point(gg.aes(fill="Metadata_clone_number", shape="Metadata_treatment"),
                  color='black', alpha=0.8) +
    gg.theme_bw() +
    gg.xlab("UMAP (X)") +
    gg.ylab("UMAP (Y)")
)

out_file = os.path.join("figures", "core_batch_umap_well_sample.png")
umap_well_embedding_gg.save(out_file, width=4, height=3, dpi=400, verbose=False)

umap_well_embedding_gg

