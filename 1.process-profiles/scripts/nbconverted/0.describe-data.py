#!/usr/bin/env python
# coding: utf-8

# # Describing Data by Batch

# In[1]:


import os
import numpy as np
import pandas as pd
import plotnine as gg

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


batches = [x for x in os.listdir("profiles") if x != ".DS_Store"]
batches


# In[3]:


def load_data(batch, profile_dir="profiles"):
    batch_dir = os.path.join(profile_dir, batch)

    backend_folders = os.listdir(batch_dir)
    plate_files = [
        os.path.join(batch_dir, x, "{}_normalized_variable_selected.csv".format(x))
        for x in backend_folders if ".DS_Store" not in x
    ]
    
    plate_dfs = []
    for plate_file in plate_files:
        plate_dfs.append(pd.read_csv(plate_file))

    return pd.concat(plate_dfs, axis="rows")

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
        }, axis="columns")
        .assign(batch=batch_name)
    )
    return result

def process_counts(batch_name, profile_dir="profiles"):
    df = load_data(batch_name, profile_dir)
    batch_count = get_count_per_batch(df, batch_name)
    treatment_count = count_treatments_per_plate(df, batch_name)
    return df, batch_count, treatment_count


# In[4]:


#load_data("2019_06_25_Batch3")


# In[5]:


batch_data = {}
for batch in batches:
    print(batch)
    df, batch_count, treatment_count = process_counts(batch)
    
    batch_data[batch] = {
            "dataframe": df,
            "metafeatures": infer_cp_features(df, metadata=True),
            "batch_count": batch_count,
            "treatment_count": treatment_count
        }


# In[6]:


all_profile_counts = []
for key, value in batch_data.items():
    all_profile_counts.append(batch_data[key]["batch_count"])

profile_counts_df = pd.concat(all_profile_counts, axis="rows")
profile_counts_df


# In[7]:


all_treatment_counts = []
for key, value in batch_data.items():
    all_treatment_counts.append(batch_data[key]["treatment_count"])

treatment_counts_df = pd.concat(all_treatment_counts, axis="rows", sort=True)
treatment_counts_df.head()


# ## Visualize Counts

# In[8]:


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


# ## Describe Metadata Counts for Each Batch

# In[9]:


batch1_40x_df = treatment_counts_df.query("batch == '2019_02_15_Batch1_40X'").dropna(axis="columns")
batch1_40x_df


# In[10]:


batch1_20x_df = treatment_counts_df.query("batch == '2019_02_15_Batch1_20X'").dropna(axis="columns")
batch1_20x_df


# In[11]:


batch2_df = treatment_counts_df.query("batch == '2019_03_20_Batch2'").dropna(axis="columns")
batch2_df


# In[12]:


batch3_df = treatment_counts_df.query("batch == '2019_06_25_Batch3'").dropna(axis="columns")
batch3_df


# In[13]:


batch4_df = treatment_counts_df.query("batch == '2019_11_11_Batch4'").dropna(axis="columns")
batch4_df


# In[14]:


batch5_df = treatment_counts_df.query("batch == '2019_11_19_Batch5'").dropna(axis="columns")
batch5_df


# In[15]:


batch6_df = treatment_counts_df.query("batch == '2019_11_20_Batch6'").dropna(axis="columns")
batch6_df


# In[16]:


batch7_df = treatment_counts_df.query("batch == '2019_11_22_Batch7'").dropna(axis="columns")
batch7_df

