#!/usr/bin/env python
# coding: utf-8

# # Describing Data by Batch

# In[1]:


import os
import numpy as np
import pandas as pd
import plotnine as gg

from pycytominer.cyto_utils import infer_cp_features

from scripts.processing_utils import load_data


# In[2]:


def get_count_per_batch(df, batch_name):
    result = (
        df
        .Metadata_Plate
        .value_counts()
        .reset_index()
        .rename({
            "index": "Metadata_Plate",
            "Metadata_Plate": "profile_count"
        }, axis="columns")
        .assign(batch=batch_name)
    )
    return result

def count_treatments_per_plate(df, batch_name):
    
    if batch_name in ["2019_02_15_Batch1_20X", "2019_02_15_Batch1_40X", "2019_03_20_Batch2"]:
        group_cols = ["Metadata_CellLine", "Metadata_Dosage", "Metadata_Plate"]
    elif batch_name in ["2019_06_25_Batch3"]:
        group_cols = ["Metadata_clone_number", "Metadata_Plate"]
    else:
        group_cols = ["Metadata_clone_number", "Metadata_treatment", "Metadata_Plate"]

    result = (
        df
        .groupby(group_cols)
        ["Metadata_Well"]
        .count()
        .reset_index()
        .rename({
            "Metadata_Well": "profile_count",
            group_cols[0]: "Metadata_clone"
        }, axis="columns")
        .assign(batch=batch_name)
    )
    
    if batch_name not in ["2019_06_25_Batch3"]:
        result = (
            result.rename({
                group_cols[1]: "Metadata_treatment"
            }, axis="columns")
        )
    return result

def process_counts(batch_name, profile_dir="profiles"):
    df = load_data(batch_name, profile_dir, combine_dfs=True)
    batch_count = get_count_per_batch(df, batch_name)
    treatment_count = count_treatments_per_plate(df, batch_name)
    return df, batch_count, treatment_count


# In[3]:


profile_dir = os.path.join("..", "0.generate-profiles", "profiles")
batches = [x for x in os.listdir(profile_dir) if x != ".DS_Store"]

batches


# In[4]:


batch_data = {}
all_clones = list()
profile_count_list = list()
for batch in batches:
    print("Now processing... {}".format(batch))
    df, batch_count, treatment_count = process_counts(batch, profile_dir=profile_dir)
    
    batch_data[batch] = {
            "dataframe": df,
            "metafeatures": infer_cp_features(df, metadata=True),
            "batch_count": batch_count,
            "treatment_count": treatment_count
        }
    
    all_clones += treatment_count.Metadata_clone.unique().tolist()
    profile_count_list.append(
        treatment_count.loc[:, ["Metadata_clone", "Metadata_treatment", "profile_count"]]
    )


# In[5]:


sample_count_df = (
    pd.DataFrame(
        pd.concat(profile_count_list, axis="rows")
        .fillna("DMSO")
        .reset_index(drop=True)
        .groupby(["Metadata_clone", "Metadata_treatment"])
        ["profile_count"]
        .sum()
    )
    .sort_values("profile_count", ascending=False)
    .reset_index()
)

sample_count_df


# In[6]:


sample_treatment_count_df = (
    sample_count_df
    .pivot_table(values="profile_count", index="Metadata_clone", columns="Metadata_treatment")
    .fillna(0)
    .astype(int)
)

sample_treatment_count_df


# In[7]:


len(set(all_clones))


# In[8]:


all_profile_counts = []
for key, value in batch_data.items():
    all_profile_counts.append(batch_data[key]["batch_count"])

profile_counts_df = pd.concat(all_profile_counts, axis="rows")
profile_counts_df


# In[9]:


all_treatment_counts = []
for key, value in batch_data.items():
    all_treatment_counts.append(batch_data[key]["treatment_count"])

treatment_counts_df = pd.concat(all_treatment_counts, axis="rows", sort=True)
treatment_counts_df.head()


# In[10]:


clone_counts_df = (
    treatment_counts_df
    .groupby(["Metadata_clone", "Metadata_treatment"])
    ["profile_count"]
    .sum()
    .reset_index()
)

output_file = os.path.join("tables", "clone_counts_bortezomib.csv")
clone_counts_df.to_csv(output_file, sep=',', index=False)

clone_counts_df


# ## Visualize Counts

# In[11]:


total_count = profile_counts_df.profile_count.sum()
total_label = "Total Profile Count: {}".format(total_count)

batch_count_gg = (
    gg.ggplot(profile_counts_df, gg.aes(y="profile_count", x="batch")) +
    gg.geom_bar(gg.aes(fill="Metadata_Plate"), stat="identity") +
    gg.coord_flip() +
    gg.theme_bw() +
    gg.ylab("Profile Count") +
    gg.xlab("Batch") +
    gg.ggtitle(total_label)
)

output_figure = os.path.join("figures", "batch_count.png")
batch_count_gg.save(output_figure, height=4, width=5.5, dpi=400, verbose=False)

batch_count_gg


# ## Output Metadata Counts for Each Batch
# 
# For quick description

# In[12]:


batch1_40x_df = treatment_counts_df.query("batch == '2019_02_15_Batch1_40X'").dropna(axis="columns")
batch1_40x_df


# In[13]:


batch1_20x_df = treatment_counts_df.query("batch == '2019_02_15_Batch1_20X'").dropna(axis="columns")
batch1_20x_df


# In[14]:


batch2_df = treatment_counts_df.query("batch == '2019_03_20_Batch2'").dropna(axis="columns")
batch2_df


# In[15]:


batch3_df = treatment_counts_df.query("batch == '2019_06_25_Batch3'").dropna(axis="columns")
batch3_df


# In[16]:


batch4_df = treatment_counts_df.query("batch == '2019_11_11_Batch4'").dropna(axis="columns")
batch4_df


# In[17]:


batch5_df = treatment_counts_df.query("batch == '2019_11_19_Batch5'").dropna(axis="columns")
batch5_df


# In[18]:


batch6_df = treatment_counts_df.query("batch == '2019_11_20_Batch6'").dropna(axis="columns")
batch6_df


# In[19]:


batch7_df = treatment_counts_df.query("batch == '2019_11_22_Batch7'").dropna(axis="columns")
batch7_df

