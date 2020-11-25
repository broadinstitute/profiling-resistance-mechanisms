#!/usr/bin/env python
# coding: utf-8

# ## Compile dataset "Batch 3"
# 
# **Gregory Way, 2020**
# 
# We acquired two plates in batch 3 - in each we collected either WT or Mutant clones
# 
# Here, we attempt to combine their raw measurements, normalize and output to apply the cloneAE signature in a later notebook.

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd

from pycytominer import normalize
from pycytominer.cyto_utils import infer_cp_features, output

sys.path.insert(0, "../2.describe-data/scripts")
from processing_utils import load_data


# In[2]:


batch = "2019_06_25_Batch3"
plates = ["MutClones", "WTClones"]
suffix = "augmented.csv.gz"

data_dir = pathlib.Path("../0.generate-profiles/profiles")
cell_count_dir = pathlib.Path("../0.generate-profiles/cell_counts/")

output_file = pathlib.Path(f"data/{batch}_combined_normalized.csv.gz")


# In[3]:


# Load and harmonize data for the given plates
df = load_data(
    batch=batch,
    plates=plates,
    profile_dir=data_dir,
    suffix=suffix,
    combine_dfs=True,
    harmonize_cols=True,
    add_cell_count=True,
    cell_count_dir=cell_count_dir
)

df = df.assign(
    Metadata_unique_sample_name=[f"profile_{x}_{batch}" for x in range(0, df.shape[0])]
)

print(df.shape)
df.head()


# In[4]:


normalized_df = df.groupby("Metadata_Plate").apply(
    lambda x: normalize(
        profiles=x,
        features="infer",
        samples="Metadata_clone_number == 'WT_parental'",
        method="standardize"
    )
)


output(
    df=normalized_df,
    output_filename=output_file,
    compression="gzip"
)

print(normalized_df.shape)
normalized_df.head()

