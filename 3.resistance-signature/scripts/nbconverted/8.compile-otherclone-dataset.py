#!/usr/bin/env python
# coding: utf-8

# ## Compile dataset of clones resistant to other drugs
# 
# **Gregory Way, 2021**
# 
# These clones are resistant to Ixazomib and CB-5083.
# I create a dataset of DMSO-treated clones to later apply the bortezomib resistance signature.

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features

sys.path.insert(0, "../2.describe-data/scripts")
from processing_utils import load_data


# In[2]:


np.random.seed(1234)


# In[3]:


data_dir = pathlib.Path("../0.generate-profiles/profiles")
cell_count_dir = pathlib.Path("../0.generate-profiles/cell_counts/")

output_dir = pathlib.Path("data")

profile_suffix = "normalized.csv.gz"


# In[4]:


datasets = {
    "ixazomib": {
        "2020_08_24_Batch9": ["218698"],
        "2020_09_08_Batch10": ["218854", "218858"],
    },
    "cb5083": {
        "2020_08_24_Batch9": ["218696", "218774"],
        "2020_09_08_Batch10": ["218852", "218856"],
    }
}


# In[5]:


full_df = []
for dataset in datasets:
    dataset_df = []
    for batch in datasets[dataset]:
        plates = datasets[dataset][batch]
        
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
            Metadata_model_split="otherclone"
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
    
    full_df.append(dataset_df)

full_df = pd.concat(full_df, axis="rows", sort=False).reset_index(drop=True)


# In[6]:


full_df


# In[7]:


# Reorder features
common_metadata = infer_cp_features(full_df, metadata=True)
morph_features = infer_cp_features(full_df)

full_df = full_df.reindex(common_metadata + morph_features, axis="columns")

print(full_df.shape)
full_df.head()


# In[8]:


pd.crosstab(full_df.Metadata_dataset, full_df.Metadata_model_split)


# In[9]:


pd.crosstab(full_df.Metadata_clone_number, full_df.Metadata_model_split)


# In[10]:


output_file = pathlib.Path(f"{output_dir}/otherclones_normalized_profiles.tsv.gz")
full_df.to_csv(output_file, sep="\t", index=False)


# In[ ]:





# In[ ]:




