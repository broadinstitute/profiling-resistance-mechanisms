#!/usr/bin/env python
# coding: utf-8

# ## Output ROC scores for otherclone datasets applied bortezomib signature
# 
# **Gregory Way, 2021**
# 
# How well is a bortezomib specific classifier able to separate clones resistant to cb5083 and ixazomib?

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
results_file = pathlib.Path(sig_dir, f"singscore_results_{dataset}.tsv.gz")

output_dir = pathlib.Path("results", "performance")

num_permutations = 100
threshold = 0

metric_comparisons = {
    "dataset": ["Metadata_dataset"],
}

roc_model_split_focus = ["ixazomib", "cb5083"]


# In[4]:


# Load data
results_df = pd.read_csv(results_file, sep="\t")

print(results_df.shape)
results_df.head()


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


# Get ROC curve information for model sets
roc_scores = []
roc_curve_data = []
for split in roc_model_split_focus:
    results_subset_df = results_df.query("Metadata_dataset == @split")
    for shuffle in [True, False]:
        roc_auc_val, roc_df = get_metrics(df=results_subset_df, return_roc_curve=True, shuffle=shuffle)

        roc_scores.append(pd.Series([roc_auc_val, split, shuffle]))
        roc_curve_data.append(roc_df.assign(model_split=split, shuffled=shuffle))

roc_scores_df = pd.DataFrame(roc_scores)
roc_scores_df.columns = ["roc_auc", "model_split", "shuffled"]
roc_curve_data_df = pd.concat(roc_curve_data).reset_index(drop=True)


# In[8]:


# Output performance results
for compare in metric_comparisons:
    full_results_df = real_metric_results[compare]
    shuffle_results_df = pd.concat(all_shuffle_results[compare]).reset_index(drop=True)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_metric_performance.tsv")
    full_results_df.to_csv(output_file, sep="\t", index=False)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_shuffle_metric_performance.tsv")
    shuffle_results_df.to_csv(output_file, sep="\t", index=False)
    
# Output ROC results
output_file = pathlib.Path(f"{output_dir}/{dataset}_bortezomibsignature_roc_auc.tsv")
roc_scores_df.to_csv(output_file, sep="\t", index=False)

output_file = pathlib.Path(f"{output_dir}/{dataset}_bortezomibsignature_roc_curve.tsv")
roc_curve_data_df.to_csv(output_file, sep="\t", index=False)

