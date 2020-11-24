#!/usr/bin/env python
# coding: utf-8

# ## Combine specific plates of profiles together to form analytical datasets
# 
# In image-based profiling, we collect data on a per-plate basis.
# This means that to increase sample size of cell lines and treatments, we often need to collect data from multiple plates.
# 
# We collected many different plates in this experiment, and many of these plates contained different cell line clones and different treatments.
# In this notebook, we combine the following plates and batches to form three complete datasets to perform downstream analysis.
# 
# The datasets are:
# 
# 1. cloneAE (bortezomib)
# 2. ixazomib (no treatment and treatment)
# 3. CB-5083 (no treatment and treatment)
# 
# ### Datasets by plate
# 
# | Plate | Number of Profiles | Batch | Dataset | 
# | :---- | :----------------- | :---- | :------ | 
# | 218363 | 60 | 2020_07_02_Batch8 | cloneAE |
# | 218362 | 60 | 2020_07_02_Batch8 | cloneAE |
# | 218361 | 60 | 2020_07_02_Batch8 | cloneAE |
# | 218360 | 60 | 2020_07_02_Batch8 | cloneAE |
# | 218361 | 60 | 2020_07_02_Batch8 | cloneAE |
# | 218360 | 60 | 2020_07_02_Batch8 | cloneAE |
# | HCT116bortezomib | 36 | 2019_02_15_Batch1_20X | cloneAE_validation |
# | 207106_exposure320 | 36 | 2019_03_20_Batch2 | cloneAE_validation |
# | 218858 | 60 | 2020_09_08_Batch10 | ixazomib_no_trt |
# | 218854 | 60 | 2020_09_08_Batch10 | ixazomib_no_trt |
# | 218698 | 60 | 2020_08_24_Batch9 | ixazomib_no_trt_validation |
# | 218859 | 60 | 2020_09_08_Batch10 | ixazomib_trt |
# | 218855 | 60 | 2020_09_08_Batch10 | ixazomib_trt |
# | 218699 | 60 | 2020_08_24_Batch9 | ixazomib_trt_validation |
# | 218856 | 60 | 2020_09_08_Batch10 | cb5083_no_trt |
# | 218852 | 60 | 2020_09_08_Batch10 | cb5083_no_trt |
# | 218774 | 60 | 2020_08_24_Batch9 | cb5083_no_trt_validation |
# | 218696 | 60 | 2020_08_24_Batch9 | cb5083_no_trt_validation |
# | 218857 | 60 | 2020_09_08_Batch10 | cb5083_trt |
# | 218853 | 60 | 2020_09_08_Batch10 | cb5083_trt |
# | 218775 | 60 | 2020_08_24_Batch9 | cb5083_trt_validation |
# | 218697 | 60 | 2020_08_24_Batch9 | cb5083_trt_validation |
# 
# Note that for the cb5083 and ixazomib datasets, we only compile untreated plates.
# See https://github.com/broadinstitute/profiling-resistance-mechanisms/issues/89 for more details.
# 
# We specifically selected validation sets that were collected on an entirely different plate.
# No info from this plate was included in training or testing.
# 
# | Plate | Dataset | Batch |
# | :--- | :----- | :--------------: |
# | 218360 | cloneAE | 2020_07_02_Batch8 |
# | 218858 | ixazomib | 2020_09_08_Batch10 |
# | 218774 | cb5083 | 2020_08_24_Batch9 |
# 
# 
# ### Other enhancements
# 
# In this notebook, I also perform several other operations to prepare the data for downstream analyses:
# 
# * Add columns `Metadata_unique_sample_name`, `Metadata_clone_type`, and `Metadata_clone_type_indicator` to all profiles.
#   * These columns determine if the clone that was passaged was resistant or senstive to the given drug and will help with identification later.
# * Perform feature selection within each dataset.
#   * I save the outcome of feature selection in a separate table ("dataset_features_selected.tsv")
# * I output the combined and processed datasets to two files:
#   * Dataset-specific compressed CSV files (no feature selection)
#   * Dataset-specific .gct files (with feature selection)

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features, write_gct

sys.path.insert(0, "../2.describe-data/scripts")
from processing_utils import load_data


# In[2]:


np.random.seed(1233)


# In[3]:


data_dir = pathlib.Path("../0.generate-profiles/profiles")
cell_count_dir = pathlib.Path("../0.generate-profiles/cell_counts/")

output_dir = pathlib.Path("data")

profile_suffix = "normalized.csv.gz"

feature_select_opts = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers",
]

test_set_size = 0.25
corr_threshold = 0.90
na_cutoff = 0


# In[4]:


datasets = {
    "cloneAE": {
        "2019_02_15_Batch1_20X": ["HCT116bortezomib"],
        "2019_03_20_Batch2": ["207106_exposure320"],
        "2020_07_02_Batch8": ["218360", "218361"],
    },
    "ixazomib": {
        "2020_08_24_Batch9": ["218698"],
        "2020_09_08_Batch10": ["218854", "218858"]
    },
    "cb5083": {
        "2020_08_24_Batch9": ["218774", "218696"],
        "2020_09_08_Batch10": ["218852", "218856"]
    }
}

# The rest of the plates belong to part of the training set
validation_plates = {
    "cloneAE": "218360",
    "ixazomib": "218858",
    "cb5083": "218774"
}


# In[5]:


# Load and concatenate profiles into a single file
all_datasets_df = []
for dataset in datasets:
    dataset_df = []
    validation_plate = validation_plates[dataset]
    for batch in datasets[dataset]:
        plates = datasets[dataset][batch]

        # Load and harmonize data for the given plates
        df = load_data(
            batch=batch,
            plates=plates,
            profile_dir=data_dir,
            suffix=profile_suffix,
            combine_dfs=True,
            harmonize_cols=True,
            add_cell_count=True,
            cell_count_dir=cell_count_dir
        )
        
        # Add important metadata features
        df = df.assign(
            Metadata_dataset=dataset,
            Metadata_batch=batch,
            Metadata_clone_type="resistant",
            Metadata_clone_type_indicator=1,
            Metadata_model_split="training"
        )

        df.loc[df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type"] = "sensitive"
        df.loc[df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type_indicator"] = 0
        dataset_df.append(df)
        
    # Merge plates of the same dataset together
    dataset_df = pd.concat(dataset_df, axis="rows", sort=False).reset_index(drop=True)
    
    # Generate a unique sample ID
    # (This will be used in singscore calculation)
    dataset_df = dataset_df.assign(
        Metadata_unique_sample_name=[f"profile_{x}_{dataset}" for x in range(0, dataset_df.shape[0])]
    )
    
    # Determine training/test/validation splits
    if dataset != "cloneAE":
        dataset_df.loc[
            dataset_df.Metadata_clone_number.astype(str) == "WT_parental", "Metadata_model_split"
        ] = "test"
        
    dataset_df.loc[
        dataset_df.Metadata_Plate.astype(str) == validation_plate, "Metadata_model_split"
    ] = "validation"
    
    training_df = (
        dataset_df
        .query("Metadata_model_split == 'training'")
        .query("Metadata_treatment == '0.1% DMSO'")
    )
    
    if dataset != "cloneAE":
        training_df = training_df.query("Metadata_clone_number != 'WT_parental'")

    train_samples, test_samples = train_test_split(
        training_df.Metadata_unique_sample_name,
        random_state=9876,
        test_size=test_set_size,
        stratify=training_df.Metadata_clone_number.astype(str)
    )

    dataset_df.loc[
        dataset_df.Metadata_unique_sample_name.isin(test_samples), "Metadata_model_split"
    ] = "test"
    dataset_df.loc[
        dataset_df.Metadata_treatment != '0.1% DMSO', "Metadata_model_split"
    ] = "perturbation"

    all_datasets_df.append(dataset_df)

all_datasets_df = pd.concat(all_datasets_df, axis="rows", sort=False).reset_index(drop=True)


# In[6]:


# Remember, we are not training with WT_parental lines
pd.crosstab(all_datasets_df.Metadata_dataset, all_datasets_df.Metadata_model_split)


# In[7]:


pd.crosstab(all_datasets_df.Metadata_clone_number, all_datasets_df.Metadata_model_split)


# In[8]:


pd.crosstab(all_datasets_df.Metadata_clone_number, all_datasets_df.Metadata_dataset)


# In[9]:


pd.crosstab(all_datasets_df.Metadata_treatment, all_datasets_df.Metadata_model_split)


# In[10]:


# Only the cloneAE (bortezomib) dataset should contain WT_parental lines for training
parental_splits = all_datasets_df.query("Metadata_clone_number == 'WT_parental'")
pd.crosstab(parental_splits.Metadata_dataset, parental_splits.Metadata_model_split)


# In[11]:


# We see a very large difference in cell count across profiles
# Remember that profiles were generated from averaging feature values for all single cells
all_datasets_df.Metadata_cell_count.hist()


# In[12]:


# Reorder features
common_metadata = infer_cp_features(all_datasets_df, metadata=True)
morph_features = infer_cp_features(all_datasets_df)

all_datasets_df = all_datasets_df.reindex(common_metadata + morph_features, axis="columns")

print(all_datasets_df.shape)
all_datasets_df.head()


# ## Apply feature selection
# 
# We apply feature selection per dataset using only the training set.
# We track which features are selected per dataset and subset.

# In[13]:


selected_features = []
for dataset in datasets:
        
    # Apply feature selection
    feature_select_df = feature_select(
        all_datasets_df.query("Metadata_dataset == @dataset").query("Metadata_model_split == 'training'"),
        operation=feature_select_opts,
        na_cutoff=na_cutoff,
    )

    dataset_features = infer_cp_features(feature_select_df)

    selected_features.append(
        pd.DataFrame(dataset_features, columns=["features"])
        .assign(dataset=dataset)
    )
    
# Output results of feature selection
all_selected_features = pd.concat(selected_features).reset_index(drop=True)

output_file = pathlib.Path(f"{output_dir}/dataset_features_selected.tsv")
all_selected_features.to_csv(output_file, sep="\t", index=False)

all_selected_features.head()


# In[14]:


# How many features are selected?
all_selected_features.groupby("dataset")["features"].count()


# ## Output Datasets
# 
# We output the feature selected datasets as .gct files, but the full feature set as compressed csvs.

# In[15]:


output_file = pathlib.Path(f"{output_dir}/bulk_profiles_analytical_set.csv.gz")
output_gct_file = pathlib.Path(f"{output_dir}/bulk_profiles_feature_select_analytical_set.gct")
    
all_datasets_df.to_csv(output_file, sep=",", compression="gzip", index=False)

# Write the gct using selected features only
feature_selected_df = all_datasets_df.loc[:, common_metadata + all_selected_features.features.unique().tolist()]
write_gct(profiles=feature_selected_df, output_file=output_gct_file)


# In[16]:


print(feature_selected_df.shape)
feature_selected_df.head()

