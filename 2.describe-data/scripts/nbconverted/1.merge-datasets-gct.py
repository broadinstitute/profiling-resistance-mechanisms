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

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features, write_gct

from scripts.processing_utils import load_data


# In[2]:


# Set constants
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers"
]

gct_dir = os.path.join("data", "gct_files")
profile_dir = os.path.join("..", "0.generate-profiles", "profiles")
cell_count_dir = os.path.join("..", "0.generate-profiles", "cell_counts")
output_dir = os.path.join("data", "merged")

suffix = "normalized.csv.gz"

# Ignore 40X batch
batches = [x for x in os.listdir(profile_dir) if x not in [".DS_Store", "2019_02_15_Batch1_40X"]]
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
        harmonize_cols=True,
        cell_count_dir=cell_count_dir
    )

    # Save normalized and non-feature selected data
    profile_batches[batch] = df
    
    # Apply feature selection
    feature_select_df = feature_select(df, operation=feature_select_ops)
        
    # Write the dataframe as a gct file for input into Morpheus
    write_gct(profiles=feature_select_df, output_file=output_gct_file)


# ## Merge Profiles Together and Output

# In[4]:


all_profiles_df = pd.concat(profile_batches.values(), sort=True).reset_index(drop=True)

all_profiles_df = all_profiles_df.assign(Metadata_clone_type="resistant")
all_profiles_df.loc[all_profiles_df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type"] = "wildtype"

meta_features = infer_cp_features(all_profiles_df, metadata=True)
cp_cols = infer_cp_features(all_profiles_df, metadata=False)

all_profiles_df = all_profiles_df.reindex(meta_features + cp_cols, axis="columns")

output_file = os.path.join(output_dir, "all_merged_profiles_before_feature_selection.csv.gz")
all_profiles_df.to_csv(output_file, index=False, compression="gzip")

print(all_profiles_df.shape)
all_profiles_df.head()


# In[5]:


all_profiles_df = feature_select(all_profiles_df, operation=feature_select_ops)

all_profiles_df = all_profiles_df.drop(["Metadata_plate_ID", "Metadata_plate_filename"], axis="columns")

print(all_profiles_df.shape)
all_profiles_df.head()


# In[6]:


output_file = os.path.join(output_dir, "all_merged_profiles.csv.gz")
all_profiles_df.to_csv(output_file, index=False, compression="gzip")

output_gct_file = os.path.join(gct_dir, "all_merged_profiles.gct")
write_gct(profiles=all_profiles_df, output_file=output_gct_file)


# ## Collapse replicates into consensus profiles

# In[7]:


median_consensus_df = (
    all_profiles_df.groupby(["Metadata_clone_number", "Metadata_treatment"])
    .median()
    .reset_index()
)

print(median_consensus_df.shape)
median_consensus_df.head()


# In[8]:


output_file = os.path.join(output_dir, "consensus_profiles.csv.gz")
median_consensus_df.to_csv(output_file, index=False, compression="gzip")

output_gct_file = os.path.join(gct_dir, "consensus_profiles.gct")
write_gct(profiles=median_consensus_df, output_file=output_gct_file)

