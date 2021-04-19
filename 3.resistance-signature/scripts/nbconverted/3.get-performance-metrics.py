#!/usr/bin/env python
# coding: utf-8

# # Calculate performance of signature
# 
# Gregory Way, 2021
# 
# I previously identified a series of morphology features that were significantly different between sensitive and resistant clones.
# I also applied this signature to all profiles from training, testing, validation, and holdout sets.
# Here, I evaluate the performance of this signature.
# 
# ## Evaluation
# 
# * Accuracy
#   - The resistant and sensitive clones were balanced, so accuracy is an appropriate measure
# * Average precision
#   - How well are we able to classify the resistant samples (number correctly identified as resistant / total resistant)
#   
# ## Shuffled results
# 
# I also randomly permute the signature score 100 times and perform the full evaluation.
# I record performance in this shuffled set as a negative control.
# 
# ## Metadata stratification
# 
# Lastly, I calculate performance in a variety of different metadata subsets. I calculate performance separately for:
# 
# 1. Across model splits (training, test, validation, holdout)
# 2. Across model splits and plates (to identify plate-specific performance)
# 3. Across model splits and clone ID (to identify if certain clones are consistently predicted differentially)

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, average_precision_score

import plotnine as gg

from utils.metrics import get_metrics, get_metric_pipeline


# In[2]:


np.random.seed(5678)


# In[3]:


# Set constants
dataset = "bortezomib"

sig_dir = pathlib.Path("results", "singscore")
results_file = pathlib.Path(sig_dir, f"singscore_results{dataset}.tsv.gz")

output_dir = pathlib.Path("results", "performance")

num_permutations = 100
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
        datasets=[dataset],
        shuffle=True,
        signature=False,
        threshold=threshold
    )
    for compare in metric_comparisons:
        metric_df = shuffle_metric_results[compare].assign(permutation=i)
        all_shuffle_results[compare].append(metric_df)


# In[7]:


# Output performance results
for compare in metric_comparisons:
    full_results_df = real_metric_results[compare]
    shuffle_results_df = pd.concat(all_shuffle_results[compare]).reset_index(drop=True)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_metric_performance.tsv")
    full_results_df.to_csv(output_file, sep="\t", index=False)
    
    output_file = pathlib.Path(f"{output_dir}/{compare}_{dataset}_shuffle_metric_performance.tsv")
    shuffle_results_df.to_csv(output_file, sep="\t", index=False)

