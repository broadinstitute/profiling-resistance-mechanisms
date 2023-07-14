#!/usr/bin/env python
# coding: utf-8

# # Summarize singscore accuracy by sample
# 
# While our signature was able to successfully predict many samples well, even those it had never seen before, many samples were miscategorized.
# 
# Here, we will summarize predictions for all wildtype and bortezomib resistant profiles.
# 
# There are four categories, with some overlap.
# 
# 1. High confidence predictions (greater than or less than the permuted Singscore, depending on status)
# 2. Accurate predictions (greater than or equal to zero, depending on status)
# 3. Inaccurate predictions (greater than or equal to zero in the opposite direction, depending on status)
# 4. Completely wrong predictions (greater or less than the permuted Singscore in the opposite direction, depending on status)
# 
# For example, a resistant profile with score 1.7 would be categorized as "high confident" if the permuted singscore is <1.7.
# A resistant profile with score -1.7 would be categorized as "completely incorrect" if the permuted singscore is >-1.7.
# 
# Completely wrong predictions are also inaccurate, and high confidence predictions are also accurate.

# In[1]:


import pathlib
import pandas as pd


# In[2]:


# Define paths
singscore_dir = pathlib.Path("..", "3.resistance-signature", "results", "singscore")

singscore_files = {
    "singscore_resultsbortezomib.tsv.gz": "initial",
    "singscore_results_valotherclones.tsv.gz": "validation",
    "singscore_results_otherclones.tsv.gz": "other",
    "singscore_results_LAST_BATCH_VALIDATIONotherclones.tsv.gz": "last_batch"
}

# Output files
output_dir = "results"
output_singscore_file = pathlib.Path(output_dir, "singscore_accuracy_per_sample.tsv")
output_singscore_summary_file = pathlib.Path(output_dir, "singscore_accuracy_summary.tsv")


# In[3]:


# Load and process singscore results
full_singscore_df = []
for singscore_file, singscore_category in singscore_files.items():
    singscore_file = pathlib.Path(f"{singscore_dir}/{singscore_file}")
    singscore_df = pd.read_csv(singscore_file, sep="\t")
    full_singscore_df.append(singscore_df)

full_singscore_df = (
    pd.concat(full_singscore_df)
    .reset_index(drop=True)
    .drop(columns="Metadata_unique_sample_name") 
    .drop_duplicates(subset=["TotalScore", "Metadata_Plate", "Metadata_Well", "Metadata_batch"])
    .reset_index()
    .rename(
        columns={"index": "profile_id"}
    )
)

# Remove ixazomib and cb5083
full_singscore_df = full_singscore_df.query(
    "Metadata_dataset not in ['cb5083', 'ixazomib']"
)

# Drop samples from the training set
full_singscore_df = full_singscore_df.query(
    "Metadata_model_split != 'training'"
)

# Create a unique id per profile
full_singscore_df.profile_id = [f"profile_id_{x}" for x in full_singscore_df.profile_id]

full_singscore_df = full_singscore_df.loc[
    :, [
        "profile_id",
        "TotalScore",
        "min_permuted_value",
        "max_permuted_value",
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_batch",
        "Metadata_model_split",
        "Metadata_clone_number",
        "Metadata_dataset",
        "Metadata_clone_type",
    ]
]

all_pass_sign = []
all_pass_permuted = []
all_completely_incorrect = []  # Must pass oppositive permuted value
for profile_idx, profile in full_singscore_df.iterrows():
    if profile["Metadata_clone_type"] == "sensitive":
        if profile["TotalScore"] < profile["min_permuted_value"]:
            pass_sign = 1
            pass_permuted = 1
        else:
            pass_permuted = 0
            if profile["TotalScore"] < 0:
                pass_sign = 1
            else:
                pass_sign = 0
            
        completely_incorrect = 0
        if profile["TotalScore"] >= profile["max_permuted_value"]:
            completely_incorrect = 1

    if profile["Metadata_clone_type"] == "resistant":
        if profile["TotalScore"] > profile["max_permuted_value"]:
            pass_sign = 1
            pass_permuted = 1
        else:
            pass_permuted = 0
            if profile["TotalScore"] > 0:
                pass_sign = 1
            else:
                pass_sign = 0

        completely_incorrect = 0
        if profile["TotalScore"] <= profile["min_permuted_value"]:
            completely_incorrect = 1

    all_pass_sign.append(pass_sign)
    all_pass_permuted.append(pass_permuted)
    all_completely_incorrect.append(completely_incorrect)
    
full_singscore_df = full_singscore_df.assign(
    pass_sign = all_pass_sign,
    pass_permuted = all_pass_permuted,
    completely_incorrect = all_completely_incorrect
)

# Create a column indicating incorrect predictions
full_singscore_df = full_singscore_df.assign(incorrect = full_singscore_df.pass_sign + full_singscore_df.pass_permuted)
full_singscore_df.incorrect = ~(full_singscore_df.incorrect > 0)
full_singscore_df.incorrect = full_singscore_df.incorrect.astype(int)

# Output to file
full_singscore_df.to_csv(output_singscore_file, index=False, sep="\t")

print(full_singscore_df.shape)
full_singscore_df.head()


# In[4]:


# Summary of completely incorrect samples
all_sample_groups_df = (
    full_singscore_df
    .groupby(["Metadata_model_split", "Metadata_clone_number"])
    ["completely_incorrect"]
    .count()
    .reset_index()
    .rename(
        columns={"completely_incorrect": "total_samples"}
    )
)

confidently_incorrect_counts_df = (
    full_singscore_df
    .query("completely_incorrect == 1")
    .groupby(["Metadata_model_split", "Metadata_clone_number"])
    ["completely_incorrect"]
    .sum()
    .reset_index()
)

all_sample_groups_df = (
    all_sample_groups_df.merge(
        confidently_incorrect_counts_df,
        on=["Metadata_model_split", "Metadata_clone_number"],
        how="left"
    )
    .fillna(0)
)

all_sample_groups_df = (
    all_sample_groups_df.assign(
        proportion_completely_incorrect = all_sample_groups_df.completely_incorrect / all_sample_groups_df.total_samples
    )
    .sort_values("proportion_completely_incorrect", ascending=False)
    .reset_index(drop=True)
)

all_sample_groups_df.head(10)


# In[5]:


full_sample_summary_df = (
    all_sample_groups_df
    .groupby("Metadata_clone_number")
    [["total_samples", "completely_incorrect"]]
    .sum()
    .reset_index()
)

full_sample_summary_df.head(5)


# In[6]:


high_confidence_summary_df = (
    full_singscore_df
    .groupby(["Metadata_clone_number"])
    ["pass_permuted"]
    .sum()
    .reset_index()
    .rename(
        columns={"pass_permuted": "high_confidence"}
    )
)

high_confidence_summary_df.head()


# In[7]:


correct_summary_df = (
    full_singscore_df
    .groupby(["Metadata_clone_number"])
    ["pass_sign"]
    .sum()
    .reset_index()
    .rename(
        columns={"pass_sign": "accurate"}
    )
)

correct_summary_df.head()


# In[8]:


incorrect_summary_df = (
    full_singscore_df
    .groupby(["Metadata_clone_number"])
    ["incorrect"]
    .sum()
    .reset_index()
)

incorrect_summary_df.head()


# In[9]:


overall_summary_df = full_sample_summary_df.merge(
    high_confidence_summary_df,
    on="Metadata_clone_number"
).merge(
    correct_summary_df,
    on="Metadata_clone_number"
).merge(
    incorrect_summary_df,
    on="Metadata_clone_number"
)

overall_summary_df.head()


# In[10]:


# Calculate proportions
overall_summary_df = (
    overall_summary_df.assign(
        prop_completely_incorrect = overall_summary_df.completely_incorrect / overall_summary_df.total_samples,
        prop_high_confidence = overall_summary_df.high_confidence / overall_summary_df.total_samples,
        prop_accurate = overall_summary_df.accurate / overall_summary_df.total_samples,
        prop_inaccurate = overall_summary_df.incorrect / overall_summary_df.total_samples
    )
    .sort_values("prop_completely_incorrect", ascending=False)
    .reset_index(drop=True)
)

overall_summary_df.to_csv(output_singscore_summary_file, index=False, sep="\t")

overall_summary_df.head()

