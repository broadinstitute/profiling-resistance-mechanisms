#!/usr/bin/env python
# coding: utf-8

# # Combine Specific Batches for Downstream Experiments
# 
# Here, I combine batches 1 (20X) and 2 because they measure clones A and E.
# I also combine batches 5, 6, and 7. These platemaps are all the same, and each measure 8 different wild-type and mutant clones.

# In[1]:


import os
import pandas as pd

from pycytominer import feature_select, write_gct

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
profile_dir = os.path.join("..", "0.generate-profiles", "profiles")
gct_dir = os.path.join("data", "gct_files")
output_dir = os.path.join("data", "merged")

suffix = "normalized.csv.gz"

batches = [x for x in os.listdir(profile_dir) if x != ".DS_Store"]
batches


# In[3]:


dataset_a = ["2019_02_15_Batch1_20X", "2019_03_20_Batch2"]
dataset_b = ["2019_11_19_Batch5", "2019_11_20_Batch6", "2019_11_22_Batch7"]


# In[4]:


dataset_a_dict = {}
dataset_b_dict = {}
for batch in batches:    

    if batch in dataset_a:
        dataset_a_dict[batch] = load_data(
            batch=batch, suffix=suffix, profile_dir=profile_dir, combine_dfs=True
        )
        
    if batch in dataset_b:
        dataset_b_dict[batch] = load_data(
            batch=batch, suffix=suffix, profile_dir=profile_dir, combine_dfs=True
        )


# ## Process and Output Dataset A

# In[5]:


dataset_a_df = pd.concat(dataset_a_dict.values()).reset_index(drop=True)

print(dataset_a_df.shape)
dataset_a_df.head()


# In[6]:


pd.crosstab(
    dataset_a_df.Metadata_CellLine,
    dataset_a_df.Metadata_Dosage
)


# In[7]:


dataset_a_name = "combined_cloneAcloneE_dataset"


# In[8]:


dataset_a_featureselect_df = feature_select(dataset_a_df, operation=feature_select_ops)

output_file = os.path.join(output_dir, "{}.csv".format(dataset_a_name))
dataset_a_featureselect_df.to_csv(output_file, index=False)

output_gct_file = os.path.join(gct_dir, "{}.gct".format(dataset_a_name))
write_gct(profiles=dataset_a_featureselect_df, output_file=output_gct_file)

print(dataset_a_featureselect_df.shape)
dataset_a_featureselect_df.head()


# ## Process and Output Dataset B

# In[9]:


dataset_b_df = pd.concat(dataset_b_dict.values()).reset_index(drop=True)

print(dataset_b_df.shape)
dataset_b_df.head()


# In[10]:


dataset_b_df.Metadata_clone_number.value_counts()


# In[11]:


dataset_b_name = "combined_four_clone_dataset"


# In[12]:


dataset_b_featureselect_df = feature_select(dataset_b_df, operation=feature_select_ops)

output_file = os.path.join(output_dir, "{}.csv".format(dataset_b_name))
dataset_b_featureselect_df.to_csv(output_file, index=False)

output_gct_file = os.path.join(gct_dir, "{}.gct".format(dataset_b_name))
write_gct(profiles=dataset_b_featureselect_df, output_file=output_gct_file)

print(dataset_b_featureselect_df.shape)
dataset_b_featureselect_df.head()

