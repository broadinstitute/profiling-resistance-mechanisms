#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd

from sklearn.metrics import accuracy_score, average_precision_score

import plotnine as gg

sys.path.insert(0, "../3.bulk-signatures/")
from utils.metrics import get_metrics, get_metric_pipeline


# In[2]:


np.random.seed(5678)


# In[3]:


# Set constants
dataset = "bortezomib"

sig_dir = pathlib.Path("results", "singscore")
results_file = pathlib.Path(sig_dir, f"singscore_results{dataset}.tsv.gz")

output_dir = pathlib.Path("results", "performance")

num_permutations = 25
threshold = 0


metric_comparisons = {
    "total": ["Metadata_model_split"],
    "plate": ["Metadata_model_split", "Metadata_Plate"],
    "sample": ["Metadata_model_split", "Metadata_clone_number"]
}


# In[4]:


# Load data
results_df = pd.read_csv(results_file, sep="\t")

print(results_df.shape)
results_df.head()


# In[5]:


# Using real predictions
real_metric_results = get_metric_pipeline(
    results_df,
    metric_comparisons,
    [dataset],
    shuffle=False,
    signature=False,
    threshold=threshold
)


# In[6]:


# Using shuffled predictions
all_shuffle_results = {compare: [] for compare in metric_comparisons}
for i in range(0, num_permutations):
    np.random.seed(i)
    shuffle_metric_results = get_metric_pipeline(
        results_df,
        metric_comparisons,
        [dataset],
        shuffle=True,
        signature=False,
        threshold=threshold
    )
    for compare in metric_comparisons:
        metric_df = shuffle_metric_results[compare].assign(permutation=i)
        all_shuffle_results[compare].append(metric_df)


# In[7]:


for compare in metric_comparisons:
    full_results_df = real_metric_results[compare]
    shuffle_results_df = pd.concat(all_shuffle_results[compare]).reset_index(drop=True)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_metric_performance.tsv")
    full_results_df.to_csv(output_file, sep="\t", index=False)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_shuffle_metric_performance.tsv")
    shuffle_results_df.to_csv(output_file, sep="\t", index=False)

