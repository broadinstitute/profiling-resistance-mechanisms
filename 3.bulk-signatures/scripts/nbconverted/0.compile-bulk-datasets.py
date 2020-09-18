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


datasets = ["four_clone", "cloneAE"]
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


dfs = {x: [] for x in datasets}
common_meta = []
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
    df = df.assign(
        Metadata_batch=batch,
        Metadata_clone_type="resistant",
        Metadata_clone_type_indicator=1
    )
    df.loc[df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type"] = "sensitive"
    df.loc[df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type_indicator"] = 0

    # Get metadata features    
    meta = infer_cp_features(df, metadata=True)
    common_meta.append(set(meta))

    # Store in dictionary
    if batch == "2020_07_02_Batch8":
        df_index = "cloneAE"
    else:
        df_index = "four_clone"

    dfs[df_index].append(df)

common_metadata = list(set.intersection(*common_meta)) + ["Metadata_sample_index"]


# In[4]:


for dataset in dfs:
    bulk_df = pd.concat(dfs[dataset], sort=False).reset_index(drop=True)
    
    # Set sample index metadata feature
    bulk_df = bulk_df.assign(
        Metadata_sample_index=[f"sample_index_{x}" for x in range(0, bulk_df.shape[0])]
    )
    
    # Reorder features
    feat = infer_cp_features(bulk_df)

    bulk_df = (
        bulk_df
        .reindex(common_metadata + feat, axis="columns")
    )

    dfs[dataset] = bulk_df
    
    print(dataset)
    print(bulk_df.shape)

common_metadata = set.intersection(*common_meta)


# In[5]:


# Apply feature selection and track which features are selected for the analysis
feature_select_dfs = {}
selected_features = []
for dataset in datasets:
    # Apply feature selection
    feature_select_dfs[dataset] = feature_select(
        dfs[dataset],
        operation=feature_select_opts,
        na_cutoff=na_cutoff,
    )

    dataset_features = infer_cp_features(feature_select_dfs[dataset])

    selected_features.append(
        pd.DataFrame(dataset_features, columns=["features"])
        .assign(dataset=dataset)
    )


# In[6]:


# Output results of feature selection
all_selected_features = pd.concat(selected_features).reset_index(drop=True)

output_file = pathlib.Path(f"{output_dir}/dataset_features_selected.tsv")
all_selected_features.to_csv(output_file, sep="\t", index=False)

all_selected_features.head()


# In[7]:


print(feature_select_dfs["four_clone"].shape)
print(feature_select_dfs["cloneAE"].shape)


# ## Output Datasets
# 
# We output the feature selected datasets as .gct files, but the full feature set as compressed csvs.

# In[8]:


for dataset in dfs:
    output_file = pathlib.Path(f"{output_dir}/bulk_profiles_{dataset}.csv.gz")
    output_gct_file = pathlib.Path(f"{output_dir}/bulk_profiles_feature_select_{dataset}.gct")
    
    dfs[dataset].to_csv(output_file, sep=",", compression="gzip", index=False)
    write_gct(profiles=feature_select_dfs[dataset], output_file=output_gct_file)

