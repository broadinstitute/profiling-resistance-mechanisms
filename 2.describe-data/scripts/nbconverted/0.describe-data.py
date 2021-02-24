#!/usr/bin/env python
# coding: utf-8

# # Describing Data by Batch

# In[1]:


import os
import pathlib
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
    
    return result

def process_counts(batch_name, profile_dir="profiles"):
    df = load_data(
        batch=batch_name,
        plates="all",
        profile_dir=profile_dir,
        suffix="normalized_feature_selected.csv.gz",
        combine_dfs=True,
        harmonize_cols=True,
        add_cell_count=False,
    )
    
    batch_count = get_count_per_batch(df, batch_name)
    treatment_count = count_treatments_per_plate(df, batch_name)
    return df, batch_count, treatment_count


# In[3]:


profile_dir = pathlib.Path("../0.generate-profiles/profiles")
batches = sorted([x for x in os.listdir(profile_dir) if x != ".DS_Store"])

batches


# In[4]:


profile_dir


# In[5]:


batch_data = {}
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
    
    profile_count_list.append(
        treatment_count.loc[:, ["Metadata_clone", "Metadata_treatment", "profile_count"]]
    )


# In[6]:


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

sample_treatment_count_df = (
    sample_count_df
    .pivot_table(
        values="profile_count",
        index="Metadata_clone",
        columns="Metadata_treatment",
        aggfunc=lambda x: x.sum()
    )
    .fillna(0)
    .astype(int)
)

sample_treatment_count_df.to_csv(
    pathlib.Path("results/sample_summary_profile_counts.tsv"), sep="\t", index=True
)

sample_treatment_count_df


# In[7]:


plot_ready_df = (
    sample_treatment_count_df
    .reset_index()
    .melt(
        id_vars="Metadata_clone",
        value_vars=sample_count_df.Metadata_treatment.unique(),
        value_name="profile_count"
    )
)

clone_order = (
    plot_ready_df
    .groupby("Metadata_clone")
    .sum()
    .reset_index()
    .sort_values(by="profile_count")
    .Metadata_clone
)

plot_ready_df.Metadata_clone = pd.Categorical(
    plot_ready_df.Metadata_clone,
    categories=clone_order
)
plot_ready_df.head()


# In[8]:


total_count = plot_ready_df.profile_count.sum()
total_label = "Total Profile Count: {}".format(total_count)

treatment_count_gg = (
    gg.ggplot(plot_ready_df, gg.aes(y="profile_count", x="Metadata_clone")) +
    gg.geom_bar(gg.aes(fill="Metadata_treatment"), position="stack", stat="identity") +
    gg.coord_flip() +
    gg.theme_bw() +
    gg.theme(axis_text_y=gg.element_text(size=5)) +
    gg.ylab("Profile Count") +
    gg.xlab("Clone") +
    gg.ggtitle(total_label)
)

output_figure = pathlib.Path("figures", "treatment_count.png")
treatment_count_gg.save(output_figure, height=4, width=5.5, dpi=400, verbose=False)

treatment_count_gg


# In[9]:


# How many unique clones
len(sample_treatment_count_df.index.unique())


# In[10]:


all_profile_counts = []
for key, value in batch_data.items():
    all_profile_counts.append(batch_data[key]["batch_count"])

profile_counts_df = pd.concat(all_profile_counts, axis="rows").reset_index(drop=True)
profile_counts_df


# In[11]:


all_treatment_counts = []
for key, value in batch_data.items():
    all_treatment_counts.append(batch_data[key]["treatment_count"])

treatment_counts_df = pd.concat(all_treatment_counts, axis="rows", sort=True).reset_index(drop=True)

treatment_counts_df.head()


# In[12]:


treatment_counts_df.query("batch == '2019_06_25_Batch3'")


# In[13]:


clone_counts_df = (
    treatment_counts_df
    .groupby(["Metadata_clone", "Metadata_treatment"])
    ["profile_count"]
    .sum()
    .reset_index()
    .sort_values(by=["Metadata_clone", "Metadata_treatment"])
)

output_file = pathlib.Path("tables/clone_counts_bortezomib.csv")
clone_counts_df.to_csv(output_file, sep=',', index=False)

clone_counts_df


# ## Visualize Counts

# In[14]:


total_count = profile_counts_df.profile_count.sum()
total_label = "Total Profile Count: {}".format(total_count)

profile_counts_df.Metadata_Plate = profile_counts_df.Metadata_Plate.astype(str)

batch_count_gg = (
    gg.ggplot(profile_counts_df, gg.aes(y="profile_count", x="batch")) +
    gg.geom_bar(gg.aes(fill="Metadata_Plate"), stat="identity") +
    gg.coord_flip() +
    gg.theme_bw() +
    gg.ylab("Profile Count") +
    gg.xlab("Batch") +
    gg.ggtitle(total_label)
)

output_figure = pathlib.Path("figures/batch_count.png")
batch_count_gg.save(output_figure, height=4, width=5.5, dpi=400, verbose=False)

batch_count_gg


# ## Output Metadata Counts for Each Batch
# 
# For quick description

# In[15]:


suspect_batches = [
    "2019_06_25_Batch3", # Too confluent, not even DMSO control
    "2019_11_11_Batch4", # Too confluent
    "2019_11_19_Batch5", # Too confluent
]


# In[16]:


non_suspect_counts = treatment_counts_df.loc[~treatment_counts_df.batch.isin(suspect_batches), :]
treatment_counts_df.Metadata_clone = pd.Categorical(
    treatment_counts_df.Metadata_clone,
    categories=clone_order
)

total_count = non_suspect_counts.profile_count.sum()
total_label = "Total Usable Profile Count: {}".format(total_count)

treatment_count_by_batch_gg = (
    gg.ggplot(treatment_counts_df, gg.aes(y="profile_count", x="Metadata_clone")) +
    gg.geom_bar(gg.aes(fill="Metadata_treatment"), position="stack", stat="identity") +
    gg.coord_flip() +
    gg.facet_wrap("~batch") +
    gg.theme_bw() +
    gg.theme(
        axis_text_y=gg.element_text(size=3.5),
        axis_text_x=gg.element_text(size=6),
        strip_text=gg.element_text(size=6, color="black"),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4")
    ) +
    gg.ylab("Profile Count") +
    gg.xlab("Clones") +
    gg.ggtitle(total_label)
)

output_figure = pathlib.Path("figures/treatment_count_by_batch.png")
treatment_count_by_batch_gg.save(output_figure, height=6.5, width=5.5, dpi=400, verbose=False)

treatment_count_by_batch_gg


# In[17]:


batch1_40x_df = treatment_counts_df.query("batch == '2019_02_15_Batch1_40X'").dropna(axis="columns")
batch1_40x_df


# In[18]:


batch1_20x_df = treatment_counts_df.query("batch == '2019_02_15_Batch1_20X'").dropna(axis="columns")
batch1_20x_df


# In[19]:


batch2_df = treatment_counts_df.query("batch == '2019_03_20_Batch2'").dropna(axis="columns")
batch2_df


# In[20]:


batch3_df = treatment_counts_df.query("batch == '2019_06_25_Batch3'").dropna(axis="columns")
batch3_df


# In[21]:


batch4_df = treatment_counts_df.query("batch == '2019_11_11_Batch4'").dropna(axis="columns")
batch4_df


# In[22]:


batch5_df = treatment_counts_df.query("batch == '2019_11_19_Batch5'").dropna(axis="columns")
batch5_df


# In[23]:


batch6_df = treatment_counts_df.query("batch == '2019_11_20_Batch6'").dropna(axis="columns")
batch6_df


# In[24]:


batch7_df = treatment_counts_df.query("batch == '2019_11_22_Batch7'").dropna(axis="columns")
batch7_df


# In[25]:


batch8_df = treatment_counts_df.query("batch == '2020_07_02_Batch8'").dropna(axis="columns")
batch8_df


# In[26]:


batch9_df = treatment_counts_df.query("batch == '2020_08_24_Batch9'").dropna(axis="columns")
batch9_df


# In[27]:


batch10_df = treatment_counts_df.query("batch == '2020_09_08_Batch10'").dropna(axis="columns")
batch10_df


# In[28]:


batch11_df = treatment_counts_df.query("batch == '2021_02_08_Batch11'").dropna(axis="columns")
batch11_df

