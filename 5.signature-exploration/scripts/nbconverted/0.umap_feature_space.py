#!/usr/bin/env python
# coding: utf-8

# In[1]:


import umap
import pathlib
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


def process_umap(data_df, umap_category, features="infer"):    
    # Prepare UMAP input by removing metadata columns
    metadata_cols = infer_cp_features(data_df, metadata=True)
    metadata_df = data_df.loc[:, metadata_cols]

    if features == "infer":
        umap_data_df = data_df.drop(metadata_cols, axis="columns")
    else:
        umap_data_df = data_df.loc[:, features]

    # Apply UMAP
    reducer = umap.UMAP(random_state=123)
    embedding = reducer.fit_transform(umap_data_df)
    
    # Setup plotting logic
    embedding_df = pd.DataFrame(embedding, columns=['umap_0', 'umap_1'])
    embedding_df = embedding_df.merge(metadata_df, left_index=True, right_index=True)
    embedding_df = embedding_df.assign(Metadata_umap_category=umap_category)
    
    return embedding_df


# In[3]:


data_dir = pathlib.Path("..", "2.describe-data", "data", "merged")
signature_dir = pathlib.Path("..", "3.resistance-signature")

profile_file = pathlib.Path(f"{data_dir}/all_merged_profiles_before_feature_selection.csv.gz")
profile_feature_selection_file = pathlib.Path(f"{data_dir}/all_merged_profiles.csv.gz")

bz_signature_file = pathlib.Path(f"{signature_dir}/results/signatures/signature_summary_bortezomib_signature.tsv.gz")

output_umap_file = pathlib.Path("results", "umap_feature_summary.tsv.gz")


# In[4]:


profile_df = pd.read_csv(profile_file, low_memory=False)

print(profile_df.shape)
profile_df.head()


# In[5]:


# Perform feature selection on the profile dataframe
profile_feature_select_df = pd.read_csv(profile_feature_selection_file)

print(profile_feature_select_df.shape)
profile_feature_select_df.head()


# In[6]:


bz_sig_df = pd.read_csv(bz_signature_file, sep="\t")


bz_sig_features = bz_sig_df.query("final_signature").features.to_list()

print(bz_sig_df.shape)
print(len(bz_sig_features))
bz_sig_df.head()


# In[7]:


all_features_except_bz_sig_features = infer_cp_features(profile_df)
all_features_except_bz_sig_features = [x for x in all_features_except_bz_sig_features if x not in bz_sig_features]
len(all_features_except_bz_sig_features)


# ## Calculate UMAP coordinates using four different feature spaces:
# 
# 1. All features
# 2. Feature selected features by traditional methods
# 3. Bortezomib signature features
# 4. All features except bortezomib signature features

# In[8]:


# 1) All feature umap
umap_all_feature_df = process_umap(
    profile_df,
    umap_category="all_features"
)

# 2) Feature selected umap
umap_fs_feature_df = process_umap(
    profile_feature_select_df,
    umap_category="feature_selected"
)

# 3) Bortezomib signature features
umap_bz_sig_df = process_umap(
    profile_df,
    features=bz_sig_features,
    umap_category="bortezomib_signature_features"
)

# 4) All features except bortezomib signature
umap_non_bz_sig_df = process_umap(
    profile_df,
    features=all_features_except_bz_sig_features,
    umap_category="all_except_bortezomib_signature_features"
)


# In[9]:


# Output umap summary
umap_df = pd.concat(
    [
        umap_all_feature_df,
        umap_fs_feature_df,
        umap_bz_sig_df,
        umap_non_bz_sig_df
    ]
).reset_index(drop=True)

umap_df.to_csv(output_umap_file, sep="\t", index=False)

print(umap_df.shape)
umap_df.head()

