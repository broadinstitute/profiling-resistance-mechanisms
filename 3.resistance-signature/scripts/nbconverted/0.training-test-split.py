#!/usr/bin/env python
# coding: utf-8

# ## Determine Data Splits
# 
# We use batch 11 data to demonstrate proof of concept that we can identify features that distinguish resistant from sensitive clones.
# 
# In this notebook, I create a single dataset, using the following procedure:
# 
# 1. Load batch 11 normalized (level 4a) data
# 2. Split the five wildtype and and five sensitive clones into training/testing sets
# 3. Keep the wildtype parental and clone A/E held out
# 4. Perform feature selection using the training data only
# 5. Also load the batch 3 data and add to the analytical set as an inference set

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split

from pycytominer import normalize, feature_select
from pycytominer.cyto_utils import infer_cp_features, write_gct

sys.path.insert(0, "../2.describe-data/scripts")
from processing_utils import load_data


# In[2]:


np.random.seed(1233)


# In[3]:


data_dir = pathlib.Path("../0.generate-profiles/profiles")
cell_count_dir = pathlib.Path("../0.generate-profiles/cell_counts/")

output_dir = pathlib.Path("data")

profile_suffix = "augmented.csv.gz"

feature_select_opts = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers",
]

corr_threshold = 0.90
na_cutoff = 0

test_set_size = 0.25


# In[4]:


dataset = "bortezomib"

batch = "2021_02_08_Batch11"
plate = "219814"


# In[5]:


clones = [
    "BZ001",
    "BZ002",
    "BZ003",
    "BZ004",
    "BZ005",
    "WT clone 01",
    "WT clone 02",
    "WT clone 03",
    "WT clone 04",
    "WT clone 05"
]


# In[6]:


# Load and harmonize data for the given plates
df = load_data(
    batch=batch,
    plates=plate,
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
df = df.assign(
    Metadata_unique_sample_name=[f"profile_{x}_{dataset}" for x in range(0, df.shape[0])]
)

df.head()


# In[7]:


# Normalize with respect to WT controls
df = normalize(
    df,
    features="infer",
    meta_features="infer",
    samples="Metadata_clone_number == 'WT_parental'",
    method="standardize",
    output_file="none"
)


# In[8]:


# Select only the uncharacterized clones
training_df = df.query("Metadata_clone_number in @clones").reset_index(drop=True)


# In[9]:


# Split data
train_samples, test_samples = train_test_split(
    training_df.Metadata_unique_sample_name,
    random_state=9876,
    test_size=test_set_size,
    stratify=training_df.Metadata_clone_number.astype(str)
)

print(len(train_samples))
print(len(test_samples))


# In[10]:


# Apply feature selection using only the training samples
feature_select_training_df = feature_select(
    training_df.query("Metadata_unique_sample_name in @train_samples"),
    operation=feature_select_opts,
    na_cutoff=na_cutoff,
)


# In[11]:


# Identify testing data
testing_df = (
    training_df
    .query("Metadata_unique_sample_name in @test_samples")
    .reindex(feature_select_training_df.columns, axis="columns")
)

testing_df.loc[:, "Metadata_model_split"] = "testing"


# In[12]:


# Load inference data (a different hold out)
inference_batch = "2019_06_25_Batch3"
inference_file = pathlib.Path(f"../3.bulk-signatures/data/{inference_batch}_combined_normalized.csv.gz")
inference_df = pd.read_csv(inference_file)

inference_df = inference_df.assign(
    Metadata_dataset="untreated_mystery_clones",
    Metadata_batch=inference_batch,
    Metadata_clone_type="resistant",
    Metadata_clone_type_indicator=1,
    Metadata_model_split="inference"
)

inference_df.loc[inference_df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type"] = "sensitive"
inference_df.loc[inference_df.Metadata_clone_number.str.contains("WT"), "Metadata_clone_type_indicator"] = 0
inference_df = inference_df.assign(
    Metadata_unique_sample_name=[f"profile_{x}_inference" for x in range(0, inference_df.shape[0])]
)

inference_df =  inference_df.reindex(feature_select_training_df.columns, axis="columns")

inference_df.Metadata_clone_number.value_counts()


# In[13]:


# Combine profiles into a single dataset and output
heldout_df = (
    df.query("Metadata_clone_number not in @clones")
    .reindex(feature_select_training_df.columns, axis="columns")
)

heldout_df.loc[:, "Metadata_model_split"] = "holdout"

bortezomib_df = pd.concat(
    [
        feature_select_training_df,
        testing_df,
        heldout_df,
        inference_df
    ],
    axis="rows"
).reset_index(drop=True)

output_file = pathlib.Path(f"{output_dir}/bortezomib_signature_analytical_set.tsv.gz")
bortezomib_df.to_csv(output_file, sep="\t", index=False)


# In[14]:


print(bortezomib_df.shape)
bortezomib_df.head()


# In[15]:


assert len(bortezomib_df.Metadata_unique_sample_name.unique()) == bortezomib_df.shape[0]

