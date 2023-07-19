#!/usr/bin/env python
# coding: utf-8

# # Evaluate misclassified samples in their different feature spaces
# 
# We compare the samples with the highest incorrect predictions against those with the highest confident accurate predictions.
# 
# We compare Wildtype and Resistant clones separately, and then compare the feature spaces together.

# In[1]:


import pathlib
import pandas as pd
import numpy as np
from scipy import stats

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


# Output file
output_ks_test_file = pathlib.Path("results", "ks_test_misclassified_differences.tsv")


# In[3]:


# Define paths
data_dir = pathlib.Path("..", "2.describe-data", "data", "merged")
signature_dir = pathlib.Path("..", "3.resistance-signature")

profile_file = pathlib.Path(f"{data_dir}/all_merged_profiles_before_feature_selection.csv.gz")
bz_signature_file = pathlib.Path(f"{signature_dir}/results/signatures/signature_summary_bortezomib_signature.tsv.gz")
accuracy_summary_file = pathlib.Path("results", "singscore_accuracy_summary.tsv")


# In[4]:


# Load profile data
profile_df = pd.read_csv(profile_file, low_memory=False)

print(profile_df.shape)
profile_df.head(3)


# In[5]:


# Load bortezomib signature features
bz_sig_df = pd.read_csv(bz_signature_file, sep="\t")

bz_sig_features = bz_sig_df.query("final_signature").features.to_list()

print(bz_sig_df.shape)
print(len(bz_sig_features))
bz_sig_df.head()


# In[6]:


# Load singscore summary
summary_df = pd.read_csv(accuracy_summary_file, sep="\t")

print(summary_df.shape)
summary_df.head()


# In[7]:


# Select samples with higher than 75 percent completely incorrect
incorrect_samples = summary_df.head(3).Metadata_clone_number.tolist()
incorrect_samples


# In[8]:


# Select samples with higher than 70 percent high confidence
correct_samples = (
    summary_df
    .sort_values(by="prop_high_confidence", ascending=False)
    .head(6)
    .Metadata_clone_number
    .tolist()
)

correct_samples


# In[9]:


# Manually define these samples in specific dictionaries
sample_comparison_dict = {
    "wildtype": {
        "correct": ["WT002", "WT012", "WT013", "WT014"],
        "incorrect": ["WT015", "WT010"]
    },
    "resistant": {
        "correct": ["BZ003", "BZ007"],
        "incorrect": ["BZ006"]
    }
}


# In[10]:


# Perform KS test for each feature for these mischaracterized columns
all_ks_results = []
for sig_feature in bz_sig_features:

    for clone_type in sample_comparison_dict.keys():
        correct_samples = sample_comparison_dict[clone_type]["correct"]
        incorrect_samples = sample_comparison_dict[clone_type]["incorrect"]

        # Subset the profile dataframe
        correct_feature_values = (
            profile_df
            .query("Metadata_clone_number in @correct_samples")
            .loc[:, sig_feature]
            .tolist()
        )

        incorrect_feature_values = (
            profile_df
            .query("Metadata_clone_number in @incorrect_samples")
            .loc[:, sig_feature]
            .tolist()
        )

        ks_stat, p_value = stats.ks_2samp(correct_feature_values, incorrect_feature_values)
        all_ks_results.append([sig_feature, clone_type, ks_stat, p_value])

# Save results to file for downstream visualization
all_ks_results = pd.DataFrame(all_ks_results)
all_ks_results.columns = ["feature", "clone_type", "ks_stat", "ks_pval"]

all_ks_results.to_csv(output_ks_test_file, sep="\t", index=False)

print(all_ks_results.shape)
all_ks_results.head()

