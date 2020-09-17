#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import pathlib
import pandas as pd

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features, write_gct

sys.path.insert(0, "../2.describe-data/scripts")
from processing_utils import load_data


# In[2]:


data_dir = pathlib.Path("../0.generate-profiles/profiles")
cell_count_dir = pathlib.Path("../0.generate-profiles/cell_counts/")

output_dir = pathlib.Path("data")
batches = [
    "2019_11_11_Batch4",
    "2019_11_19_Batch5",
    "2019_11_20_Batch6",
    "2019_11_22_Batch7",
    "2020_07_02_Batch8",
]

profile_suffix = "normalized.csv.gz"
drop_cols = ["Metadata_plate_filename"]

feature_select_opts = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers",
]
corr_threshold = 0.95
na_cutoff = 0


# In[3]:


dfs = {"four_clone": [], "cloneAE": []}
for batch in batches:
    # Load and harmonize data
    df = load_data(
        batch=batch,
        profile_dir=data_dir,
        suffix=profile_suffix,
        combine_dfs=True,
        harmonize_cols=True,
        cell_count_dir=cell_count_dir
    )
    
    # Add important metadata features
    df = df.assign(Metadata_batch=batch, Metadata_clone_type="resistant")
    df.loc[df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type"] = "sensitive"

    # Store in dictionary
    if batch == "2020_07_02_Batch8":
        df_index = "cloneAE"
    else:
        df_index = "four_clone"

    dfs[df_index].append(df)


# In[4]:


for dataset in dfs:
    bulk_df = pd.concat(dfs[dataset], sort=False).reset_index(drop=True)
    
    # Reorder features
    feat = infer_cp_features(bulk_df)
    meta = infer_cp_features(bulk_df, metadata=True)
    bulk_df = bulk_df.reindex(meta + feat, axis="columns").drop(drop_cols, axis="columns")

    dfs[dataset] = bulk_df
    print(dataset)
    print(bulk_df.shape)


# In[5]:


# Apply feature selection in only the four_clone dataset and reindex features in clone AE
dfs["four_clone"] = feature_select(
    dfs["four_clone"],
    operation=feature_select_opts,
    na_cutoff=na_cutoff,
)

dfs["cloneAE"] = dfs["cloneAE"].reindex(dfs["four_clone"].columns, axis="columns")
dfs["cloneAE"].head()


# In[6]:


print(dfs["four_clone"].shape)
print(dfs["cloneAE"].shape)


# In[7]:


for dataset in dfs:
    output_file = pathlib.Path(f"{output_dir}/bulk_profiles_{dataset}.csv.gz")
    output_gct_file = pathlib.Path(f"{output_dir}/bulk_profiles_{dataset}.gct")
    
    dfs[dataset].to_csv(output_file, sep=",", compression="gzip", index=False)
    write_gct(profiles=dfs[dataset], output_file=output_gct_file)

