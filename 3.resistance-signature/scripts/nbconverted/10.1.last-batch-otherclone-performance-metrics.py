#!/usr/bin/env python
# coding: utf-8

# ## Output ROC scores for otherclone datasets applied bortezomib signature
# 
# **Gregory Way, 2021**
# 
# How well is a bortezomib specific classifier able to separate clones resistant to cb5083 and ixazomib?
# 
# **Yu Han, 2021**
# 
# How well is a bortezomib specific classifier able to separate WT and BZ clone types?
# 
# Copied from Greg's original script of 10.0; Changed variable names to load in the new batch data results from scripts 8.1 and 9.1. Commented out/removed all lines associated with 'split'. 

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, average_precision_score

import plotnine as gg

from utils.metrics import get_metrics, get_metric_pipeline


# In[2]:


np.random.seed(56789)


# In[3]:


# Set constants
dataset = "otherclones"

sig_dir = pathlib.Path("results", "singscore")
results_file = pathlib.Path(sig_dir, f"singscore_results_LAST_BATCH_VALIDATION{dataset}.tsv.gz")

output_dir = pathlib.Path("results", "performance")

num_permutations = 100
threshold = 0

metric_comparisons = {
    "dataset": ["Metadata_clone_type"],
}

#roc_model_split_focus = ["resistant", "sensitive"]


# In[4]:


# Load data and double check if there are two classes (e.g., resistant and sensitve)
results_df = pd.read_csv(results_file, sep="\t")

print(results_df.shape)


# In[5]:


# Get performance metrics using real predictions
real_metric_results = get_metric_pipeline(
    results_df,
    metric_comparisons,
    [dataset],
    shuffle=False,
    signature=False,
    threshold=threshold
)


# In[6]:


# Get performance metrics using shuffled predictions
all_shuffle_results = {compare: [] for compare in metric_comparisons}
for i in range(0, num_permutations):
    np.random.seed(i)
    shuffle_metric_results = get_metric_pipeline(
        results_df,
        metric_comparisons,
        datasets=[dataset],
        shuffle=True,
        signature=False,
        threshold=threshold
    )
    for compare in metric_comparisons:
        metric_df = shuffle_metric_results[compare].assign(permutation=i)
        all_shuffle_results[compare].append(metric_df)


# In[7]:


#double check the data counts are correct
results_df.Metadata_clone_type_indicator.value_counts()


# In[8]:


results_df.groupby(['Metadata_clone_type', 'Metadata_clone_type_indicator','Metadata_clone_number' ]).size().reset_index(name='counts')


# In[9]:


# Get ROC curve information for model sets
roc_scores = []
roc_curve_data = []
#for split in roc_model_split_focus:
    #results_subset_df = results_df.query("Metadata_clone_type == @split")
for shuffle in [True, False]:
    roc_auc_val, roc_df = get_metrics(df=results_df, return_roc_curve=True, shuffle=shuffle)

    roc_scores.append(pd.Series([roc_auc_val, shuffle]))
    roc_curve_data.append(roc_df.assign(shuffled=shuffle))

roc_scores_df = pd.DataFrame(roc_scores)
roc_scores_df.columns = ["roc_auc", "shuffled"]
roc_curve_data_df = pd.concat(roc_curve_data).reset_index(drop=True)


# In[10]:


roc_curve_data_df.head()


# In[11]:


roc_scores_df


# In[12]:


# Output performance results
for compare in metric_comparisons:
    full_results_df = real_metric_results[compare]
    shuffle_results_df = pd.concat(all_shuffle_results[compare]).reset_index(drop=True)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_metric_performance_LAST_BATCH_VALIDATION.tsv")
    full_results_df.to_csv(output_file, sep="\t", index=False)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_shuffle_metric_performance_LAST_BATCH_VALIDATION.tsv")
    shuffle_results_df.to_csv(output_file, sep="\t", index=False)
    
# Output ROC results
output_file = pathlib.Path(f"{output_dir}/{dataset}_bortezomibsignature_roc_auc_LAST_BATCH_VALIDATION.tsv")
roc_scores_df.to_csv(output_file, sep="\t", index=False)

output_file = pathlib.Path(f"{output_dir}/{dataset}_bortezomibsignature_roc_curve_LAST_BATCH_VALIDATION.tsv")
roc_curve_data_df.to_csv(output_file, sep="\t", index=False)

