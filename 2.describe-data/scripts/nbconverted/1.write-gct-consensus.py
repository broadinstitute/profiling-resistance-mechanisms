#!/usr/bin/env python
# coding: utf-8

# ## Write all Profiles to GCT for Heatmap Visualization
# 
# **Gregory Way, 2020**
# 
# I also build consensus signatures for all unique treatments and output associated files.

# In[1]:


import os
import pandas as pd

from pycytominer import (
    feature_select,
    write_gct
)

from pycytominer.consensus import modz
from pycytominer.cyto_utils import infer_cp_features

from scripts.processing_utils import load_data


# In[2]:


# Set constants
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blacklist",
    "drop_outliers"
]
gct_dir = os.path.join("data", "gct_files")
profile_dir = os.path.join("..", "0.generate-profiles", "profiles")
cell_count_dir = os.path.join("..", "0.generate-profiles", "cell_counts")
output_dir = os.path.join("data", "merged")

suffix = "normalized.csv.gz"

batches = [x for x in os.listdir(profile_dir) if x != ".DS_Store"]
batches


# In[3]:


profile_batches = {}
for batch in batches:
    # Build output information
    output_gct_dir = os.path.join(gct_dir, batch)
    os.makedirs(output_gct_dir, exist_ok=True)
    output_gct_file = os.path.join(
        output_gct_dir, "{}_feature_select.gct".format(batch)
    )
    
    # Load the profile data and add cell counts
    df = load_data(
        batch=batch,
        suffix=suffix,
        profile_dir=profile_dir,
        combine_dfs=True,
        add_cell_count=True,
        cell_count_dir=cell_count_dir
    )
    
    # Save normalized and non-feature selected data
    profile_batches[batch] = df
    
    # Apply feature selection again - this is particularly important for batches
    # with multiple plates
    df = feature_select(df, operation=feature_select_ops)
    
    # Write the dataframe as a gct file for input into Morpheus
    write_gct(profiles=df, output_file=output_gct_file)


# ## Merge Profiles Together and Output

# In[4]:


all_profiles_df = pd.concat(profile_batches.values(), sort=True).reset_index(drop=True)

meta_features = infer_cp_features(all_profiles_df, metadata=True)
cp_cols = infer_cp_features(all_profiles_df, metadata=False)

all_profiles_df = all_profiles_df.reindex(meta_features + cp_cols, axis="columns")

print(all_profiles_df.shape)
all_profiles_df.head()


# In[5]:


all_profiles_df = feature_select(all_profiles_df, operation=feature_select_ops)

print(all_profiles_df.shape)
all_profiles_df.head()


# In[6]:


output_file = os.path.join(output_dir, "all_merged_profiles.csv.gz")
all_profiles_df.to_csv(output_file, index=False, compression="gzip")


# ## Generate Consensus Signatures

# In[7]:


consensus_data = {}
for batch in profile_batches:
    meta_features = infer_cp_features(profile_batches[batch], metadata=True)
    meta_features = [x for x in meta_features if "well" not in x.lower()]
    meta_features = [x for x in meta_features if "site" not in x.lower()]
    
    consensus_df = (
        profile_batches[batch]
        .groupby(meta_features)
        .median()
        .drop("Metadata_Site", axis="columns")
        .reset_index(drop=False)
    )
    
    consensus_data[batch] = consensus_df.reset_index()


# In[8]:


full_consensus_df = (
    pd.concat(consensus_data.values(), sort=True)
    .reset_index(drop=True)
)

meta_features = infer_cp_features(full_consensus_df, metadata=True)
cp_cols = infer_cp_features(full_consensus_df, metadata=False)

full_consensus_df = (
    full_consensus_df
    .reindex(meta_features + cp_cols, axis="columns")
    .drop("Metadata_cell_count", axis="columns")
)

print(full_consensus_df.shape)
full_consensus_df.head()


# In[9]:


consensus_df = feature_select(full_consensus_df, operation=feature_select_ops)

print(consensus_df.shape)
consensus_df.head()


# In[10]:


output_gct_file = os.path.join(gct_dir, "consensus_feature_select.gct")
write_gct(profiles=consensus_df, output_file=output_gct_file)

