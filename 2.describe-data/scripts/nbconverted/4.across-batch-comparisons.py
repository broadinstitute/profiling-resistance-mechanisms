#!/usr/bin/env python
# coding: utf-8

# ## Determine across batch correlations
# 
# **Gregory Way, 2021**
# 
# We've collected Cell Painting readouts for the same cell line clones and perturbations across many different batches. Determine across batch correlations here.

# In[1]:


import pathlib
import pandas as pd
import plotnine as gg

from cytominer_eval import evaluate
from cytominer_eval.transform import metric_melt
from cytominer_eval.operations.util import assign_replicates
from pycytominer.cyto_utils import infer_cp_features

from scripts.processing_utils import load_data


# In[2]:


data_dir = pathlib.Path("../0.generate-profiles/profiles")

cloneAE_data = {
    "2019_03_20_Batch2": ["207106_exposure320"],
    "2020_07_02_Batch8": ["218360", "218361"],
    "2021_02_08_Batch11": ["219814"]
}

profile_suffix = "normalized.csv.gz"


# In[3]:


full_datasets_df = []
for batch in cloneAE_data:
    plates = cloneAE_data[batch]
    
    # Load and harmonize data for the given plates
    df = load_data(
        batch=batch,
        plates=plates,
        profile_dir=data_dir,
        suffix=profile_suffix,
        combine_dfs=True,
        harmonize_cols=True,
        add_cell_count=False
    )
    
    full_datasets_df.append(df)


# In[4]:


full_datasets_df = pd.concat(full_datasets_df, axis="rows", sort=False).reset_index(drop=True)

features = infer_cp_features(full_datasets_df)
meta_features = infer_cp_features(full_datasets_df, metadata=True)

print(full_datasets_df.shape)
full_datasets_df.head()


# In[5]:


# Melt a pairwise similarity matrix
similarity_melted_df = metric_melt(
    df=full_datasets_df,
    features=features,
    metadata_features=meta_features,
    similarity_metric="pearson",
    eval_metric="grit",
)

# Determine which samples are replicates
replicate_groups = ["Metadata_batch", "Metadata_clone_number", "Metadata_treatment"]
similarity_melted_df = assign_replicates(
    similarity_melted_df=similarity_melted_df,
    replicate_groups=replicate_groups
)


# In[6]:


same_treatment_across_batch_replicate = (
    similarity_melted_df.Metadata_clone_number_replicate &
    similarity_melted_df.Metadata_treatment_replicate &
    ~(similarity_melted_df.Metadata_batch_replicate)
)

different_treatment_within_batch = (
    similarity_melted_df.Metadata_batch_replicate &
    ~(similarity_melted_df.group_replicate) 
)


# In[7]:


similarity_melted_df = similarity_melted_df.assign(comparison_category="across_batch_nonreplicate")
similarity_melted_df.loc[similarity_melted_df.group_replicate, "comparison_category"] = "same_batch_replicate"
similarity_melted_df.loc[same_treatment_across_batch_replicate, "comparison_category"] = "across_batch_replicate"
similarity_melted_df.loc[different_treatment_within_batch, "comparison_category"] = "same_batch_nonreplicate"


# In[8]:


similarity_melted_df.comparison_category.value_counts()


# In[15]:


similarity_melted_df.Metadata_clone_number_pair_a.value_counts()


# In[9]:


(
    gg.ggplot(similarity_melted_df, gg.aes(x="similarity_metric"))
    + gg.geom_density(gg.aes(fill="comparison_category"), alpha=0.5)
    + gg.theme_bw()
)


# In[10]:


(
    gg.ggplot(similarity_melted_df, gg.aes(x="similarity_metric"))
    + gg.geom_density(gg.aes(fill="comparison_category"), alpha=0.5)
    + gg.theme_bw()
    + gg.facet_wrap("~Metadata_batch_pair_a")
)


# In[17]:


(
    gg.ggplot(
        similarity_melted_df.query("Metadata_clone_number_pair_a in ['WT_parental', 'CloneA', 'CloneE']"),
        gg.aes(x="similarity_metric"))
    + gg.geom_density(gg.aes(fill="comparison_category"), alpha=0.5)
    + gg.theme_bw()
    + gg.facet_wrap("~Metadata_batch_pair_a")
)


# In[18]:


(
    gg.ggplot(
        similarity_melted_df.query("Metadata_clone_number_pair_a not in ['WT_parental', 'CloneA', 'CloneE']"),
        gg.aes(x="similarity_metric"))
    + gg.geom_density(gg.aes(fill="comparison_category"), alpha=0.5)
    + gg.theme_bw()
    + gg.facet_wrap("~Metadata_batch_pair_a")
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# Load and process data
data_dir = pathlib.Path("data/merged")
file = pathlib.Path(data_dir, "all_merged_profiles.csv.gz")

data_df = pd.read_csv(file)

# Create columns for grit calculation
data_df = data_df.assign(
    Metadata_treatment_group_id=data_df.Metadata_clone_number + "_" + data_df.Metadata_treatment
)
data_df = data_df.assign(
    Metadata_treatment_profile_id=data_df.Metadata_treatment_group_id + "_" + data_df.Metadata_batch
)

features = infer_cp_features(data_df)
meta_features = infer_cp_features(data_df, metadata=True)

print(data_df.shape)
data_df.head()


# In[3]:


# Select only batches that have "WT_parental_0.1% DMSO" as a treatment group
grit_control_pert = "WT_parental_0.1% DMSO"

select_batches = sorted(
    data_df.query("Metadata_treatment_group_id == @grit_control_pert").Metadata_batch.unique().tolist()
)

select_batches


# In[4]:


all_grit_control_perts = data_df.query("Metadata_batch in @select_batches").query("Metadata_treatment_group_id == @grit_control_pert").Metadata_treatment_profile_id.unique()
all_grit_control_perts


# In[5]:


meta_features


# In[6]:


data_df.Metadata_treatment_profile_id.value_counts()


# In[7]:


# Get replicate correlation
percent_strong, corr_df = evaluate(
    profiles=data_df,
    features=features,
    meta_features=meta_features,
    replicate_groups=["Metadata_clone_number", "Metadata_treatment"],
    operation="replicate_reproducibility",
    replicate_reproducibility_return_median_cor=True
)


# In[8]:


percent_strong


# In[9]:


corr_df.head()


# In[11]:


# Get technical grit for batch
replicate_groups = {
    "profile_col": "Metadata_treatment_profile_id",
    "replicate_group_col": "Metadata_treatment_group_id"
}

grit_df = evaluate(
    profiles=data_df.query("Metadata_batch in @select_batches"),
    features=features,
    meta_features=meta_features,
    replicate_groups=replicate_groups,
    operation="grit",
    grit_control_perts=all_grit_control_perts
)


# In[12]:


grit_df.sort_values(by="grit", ascending=False)


# In[ ]:




