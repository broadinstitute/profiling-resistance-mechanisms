#!/usr/bin/env python
# coding: utf-8

# ## Compile dataset of clones resistant to other drugs
# 
# **Gregory Way, 2021**
# 
# **Yu Han, 2021**
# 
# This script is modified from Greg Way's original scripts of 8.compile-otherclone-dataset.
# 
# This dataset includes new batches of 24~27 including WT (10, 12-15) and Bz-resistant (6-10) clones.
# They are DMSO-treated clones to later apply the bortezomib resistance signature.

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


datasets_val = {
        "2021_08_02_Batch24": ["221057"],
        "2021_08_02_Batch25": ["221058"],
        "2021_08_03_Batch26": ["221093"],
        "2021_08_03_Batch27": ["221094"],
}


# In[5]:


#added 'val' to each original variable names from 8.0

full_df_val = []
for dataset_val in datasets_val:
    dataset_df_val = []
    for batch in datasets_val:
        plates = datasets_val[batch]
        
        df_val = load_data(
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
        df_val = df_val.assign(
            Metadata_dataset=dataset_val,
            Metadata_batch=batch,
            Metadata_clone_type="resistant",
            Metadata_clone_type_indicator=1,
            Metadata_model_split="otherclone"
        )

        df_val.loc[df_val.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type"] = "sensitive"
        df_val.loc[df_val.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type_indicator"] = 0
        dataset_df_val.append(df_val)

    # Merge plates of the same dataset together
    dataset_df_val = pd.concat(dataset_df_val, axis="rows", sort=False).reset_index(drop=True)
    
    # Generate a unique sample ID
    # (This will be used in singscore calculation)
    dataset_df_val = dataset_df_val.assign(
        Metadata_unique_sample_name=[f"profile_{x}_{dataset_val}" for x in range(0, dataset_df_val.shape[0])]
    )
    
    full_df_val.append(dataset_df_val)

#remove other clone types from the df
full_df_val = pd.concat(full_df_val, axis="rows", sort=False).reset_index(drop=True)
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'WT_parental']
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'CloneA']
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'CloneE']
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'TX clone 2']
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'TX clone 3']
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'TX clone 4']
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'TX clone 5']
full_df_val = full_df_val[full_df_val.Metadata_clone_number != 'TX clone 6']


# In[6]:


full_df_val.head()


# In[7]:


# Reorder features
common_metadata = infer_cp_features(full_df_val, metadata=True)
morph_features = infer_cp_features(full_df_val)

full_df_val = full_df_val.reindex(common_metadata + morph_features, axis="columns")

print(full_df_val.shape)


# In[8]:


pd.crosstab(full_df_val.Metadata_clone_type_indicator, full_df_val.Metadata_model_split)


# In[9]:


pd.crosstab(full_df_val.Metadata_clone_number, full_df_val.Metadata_model_split)


# In[10]:


#saved output file as 'otherclones_normalized_profiles_LAST_BATCH_VALIDATION.tsv.gz'
output_file = pathlib.Path(f"{output_dir}/otherclones_normalized_profiles_LAST_BATCH_VALIDATION.tsv.gz")
full_df_val.to_csv(output_file, sep="\t", index=False)


# In[ ]:





# In[ ]:




